"""Bundled SQLite-backed lookup and fuzzy search for microorganism profiles."""

from __future__ import annotations

import copy
import json
import re
import sqlite3
import unicodedata
from contextlib import ExitStack
from difflib import SequenceMatcher
from importlib import resources
from pathlib import Path
from typing import Any

SEARCH_MODES = {"auto", "name"}

DATASET_LABELS = {
    "ja": {
        "bacteria": "細菌",
        "bacteria_fish": "魚病細菌",
        "fungi": "真菌",
    },
    "en": {
        "bacteria": "Bacteria",
        "bacteria_fish": "Fish-pathogenic bacteria",
        "fungi": "Fungi",
    },
}

ANNOTATION_LABELS = {
    "ja": {
        "infection_law": "感染症法（特定病原体等）",
        "bsl_niid": "国立感染症研究所BSL",
        "bsl_bsj": "日本細菌学会BSL",
        "trba": "TRBAリスクグループ",
        "domestic_animal_law": "家畜伝染病予防法",
        "plant_protection": "植物防疫法",
        "cartagena": "カルタヘナ法（実験分類）",
        "foreign_exchange": "外為法",
        "fish_pathogen": "魚介類病原菌",
        "plant_pathogen": "植物病原菌",
        "housing_fungi": "住環境菌",
        "fish_pathogen_profile": "魚病情報",
    },
    "en": {
        "infection_law": "Infectious Diseases Control Law",
        "bsl_niid": "NIID BSL",
        "bsl_bsj": "Japanese Society for Bacteriology BSL",
        "trba": "TRBA risk group",
        "domestic_animal_law": "Domestic Animal Infectious Diseases Control Law",
        "plant_protection": "Plant Protection Act",
        "cartagena": "Cartagena Act experimental classification",
        "foreign_exchange": "Foreign Exchange and Foreign Trade Act",
        "fish_pathogen": "Fish pathogen",
        "plant_pathogen": "Plant pathogen",
        "housing_fungi": "Housing environment fungi",
        "fish_pathogen_profile": "Fish pathogen profile",
    },
}


def normalize_name(value: str) -> str:
    """Normalize organism names for deterministic exact and fuzzy matching."""
    normalized = unicodedata.normalize("NFKC", value or "")
    normalized = normalized.replace(" ", "").replace("\u3000", "")
    normalized = normalized.casefold()
    normalized = re.sub(r"\([^)]*\)", "", normalized)
    normalized = re.sub(r"[\-\u2010\u2011\u2012\u2013\u2014\u2015\u30FC\u2212\uFF0D]", "", normalized)
    normalized = re.sub(r"[.,;:'\"`´・･/\\\[\]{}<>]", "", normalized)
    return normalized.strip()


def _name_similarity(query_raw: str, query_normalized: str, text_raw: str, text_normalized: str) -> tuple[float, str]:
    """Return a bounded similarity score and match strategy for name search."""
    query_raw_norm = (query_raw or "").strip().casefold()
    text_raw_norm = (text_raw or "").strip().casefold()
    query_key = (query_normalized or "").strip()
    text_key = (text_normalized or "").strip()

    if query_raw_norm and query_raw_norm == text_raw_norm:
        return 1.0, "exact_raw"
    if query_key and query_key == text_key:
        return 1.0, "exact_normalized"
    if query_raw_norm and query_raw_norm in text_raw_norm:
        coverage = len(query_raw_norm) / max(len(text_raw_norm), 1)
        return max(0.8, min(0.98, 0.72 + coverage * 0.26)), "contains_raw"
    if query_key and query_key in text_key:
        coverage = len(query_key) / max(len(text_key), 1)
        return max(0.78, min(0.96, 0.7 + coverage * 0.24)), "contains_normalized"
    if text_key and text_key in query_key:
        coverage = len(text_key) / max(len(query_key), 1)
        return max(0.72, min(0.9, 0.66 + coverage * 0.2)), "contains_inverse"

    ratio_left = query_key or query_raw_norm
    ratio_right = text_key or text_raw_norm
    if not ratio_left or not ratio_right:
        return 0.0, "none"
    return SequenceMatcher(None, ratio_left, ratio_right).ratio(), "sequence_ratio"


def _safe_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


class BioDatabase:
    """In-memory runtime database loaded from the published `ra_bio` artifacts."""

    _instance: BioDatabase | None = None

    def __init__(self, db_path: str | Path | None = None):
        self._resource_stack = ExitStack()
        self._instance_key = self._normalize_instance_key(db_path)
        self.bio_root: Path
        self.sqlite_path: Path
        self._configure_paths(db_path)

        self._loaded = False
        self._profiles_by_cluster: dict[str, dict[str, Any]] = {}
        self._cluster_summary: dict[str, dict[str, Any]] = {}
        self._alias_rows: list[dict[str, str]] = []
        self._snapshots: dict[str, dict[str, Any]] = {}

    @classmethod
    def get_instance(cls, db_path: str | Path | None = None) -> BioDatabase:
        """Get a singleton instance bound to a specific DB path or the bundled DB."""
        path_key = cls._normalize_instance_key(db_path)
        if cls._instance is None or cls._instance._instance_key != path_key:
            if cls._instance is not None:
                cls._instance.close()
            cls._instance = cls(db_path)
            cls._instance._load_data()
        return cls._instance

    @classmethod
    def reset_instance(cls) -> None:
        """Reset the singleton instance, mainly for tests."""
        if cls._instance is not None:
            cls._instance.close()
        cls._instance = None

    @staticmethod
    def _normalize_instance_key(db_path: str | Path | None) -> str:
        if db_path is None:
            return "__bundled__"
        text = str(db_path).strip()
        if not text:
            return "__bundled__"
        return str(Path(text))

    def close(self) -> None:
        """Release extracted package-resource state."""
        self._resource_stack.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def _configure_paths(self, db_path: str | Path | None) -> None:
        if db_path is None or not str(db_path).strip():
            bundled_path = self._resolve_bundled_sqlite_path()
            self.bio_root = bundled_path.parent
            self.sqlite_path = bundled_path
            return

        candidate = Path(db_path)
        if candidate.suffix == ".sqlite3" or candidate.is_file():
            if not candidate.exists():
                raise FileNotFoundError(f"Bio DB SQLite not found: {candidate}")
            self.bio_root = candidate.parent
            self.sqlite_path = candidate
            return

        sqlite_path = candidate / "bio.sqlite3"
        if not sqlite_path.exists():
            raise FileNotFoundError(f"Bio DB SQLite not found under repo layout: {sqlite_path}")
        self.bio_root = candidate
        self.sqlite_path = sqlite_path

    def _resolve_bundled_sqlite_path(self) -> Path:
        resource = resources.files("ra_bio").joinpath("data").joinpath("bio.sqlite3")
        return Path(self._resource_stack.enter_context(resources.as_file(resource)))

    def _load_data(self) -> None:
        self._profiles_by_cluster.clear()
        self._cluster_summary.clear()
        self._alias_rows = []
        self._snapshots.clear()

        connection = sqlite3.connect(self.sqlite_path)
        connection.row_factory = sqlite3.Row
        try:
            for row in connection.execute("SELECT * FROM source_snapshots ORDER BY dataset_id"):
                record = {key: row[key] for key in row.keys()}
                self._snapshots[(record.get("dataset_id") or "").strip()] = record

            for row in connection.execute("SELECT * FROM organisms ORDER BY canonical_name, cluster_id"):
                record = {key: row[key] for key in row.keys()}
                cluster_id = (record.get("cluster_id") or "").strip()
                if not cluster_id:
                    continue

                datasets = json.loads(record["datasets_json"] or "[]")
                scientific_names = json.loads(record["scientific_names_json"] or "[]")
                aliases = json.loads(record["aliases_json"] or "[]")
                profile = json.loads(record["profile_json"] or "{}")

                self._profiles_by_cluster[cluster_id] = profile
                self._cluster_summary[cluster_id] = {
                    "cluster_id": cluster_id,
                    "canonical_name": (record.get("canonical_name") or "").strip(),
                    "preferred_scientific_name": (record.get("preferred_scientific_name") or "").strip(),
                    "datasets": datasets,
                    "datasets_set": set(datasets),
                    "scientific_names": scientific_names,
                    "aliases": aliases,
                }

            for row in connection.execute(
                "SELECT cluster_id, alias_name, normalized_alias, source_type, source_ref "
                "FROM organism_aliases ORDER BY cluster_id, alias_name"
            ):
                record = {key: (row[key] or "") for key in row.keys()}
                cluster_id = record["cluster_id"].strip()
                if cluster_id and cluster_id in self._cluster_summary:
                    self._alias_rows.append(record)
        finally:
            connection.close()

        self._loaded = True

    def _ensure_loaded(self) -> None:
        if not self._loaded:
            self._load_data()

    def get_source_snapshots(self) -> list[dict[str, Any]]:
        """Return source snapshot metadata for the bundled dataset."""
        self._ensure_loaded()
        return [
            copy.deepcopy(self._snapshots[key])
            for key in sorted(self._snapshots)
        ]

    def get_runtime_status(self) -> dict[str, Any]:
        """Return lightweight runtime metadata suitable for health checks."""
        self._ensure_loaded()
        return {
            "status": "ok",
            "sqlite_path": str(self.sqlite_path),
            "using_bundled_database": self._instance_key == "__bundled__",
            "organism_count": len(self._profiles_by_cluster),
            "alias_count": len(self._alias_rows),
            "source_count": len(self._snapshots),
            "datasets": sorted(self._snapshots),
        }

    def _add_hit(
        self,
        hits_by_cluster: dict[str, dict[str, Any]],
        *,
        cluster_id: str,
        score: float,
        match_type: str,
        matched_value: str,
        source: str,
    ) -> None:
        summary = self._cluster_summary[cluster_id]
        existing = hits_by_cluster.get(cluster_id)
        if existing is None:
            hits_by_cluster[cluster_id] = {
                "cluster_id": cluster_id,
                "canonical_name": summary["canonical_name"],
                "preferred_scientific_name": summary["preferred_scientific_name"],
                "datasets": summary["datasets"],
                "score": score,
                "match_type": match_type,
                "matched_value": matched_value,
                "match_sources": {source},
                "matched_terms": [matched_value],
            }
            return

        existing["match_sources"].add(source)
        if matched_value not in existing["matched_terms"]:
            existing["matched_terms"].append(matched_value)
        if score > existing["score"]:
            existing["score"] = score
            existing["match_type"] = match_type
            existing["matched_value"] = matched_value

    def _rank_hits(self, hits_by_cluster: dict[str, dict[str, Any]], limit: int) -> list[dict[str, Any]]:
        ranked = sorted(
            hits_by_cluster.values(),
            key=lambda item: (-item["score"], item["canonical_name"], item["preferred_scientific_name"]),
        )
        payload: list[dict[str, Any]] = []
        for item in ranked[:limit]:
            payload.append(
                {
                    "cluster_id": item["cluster_id"],
                    "canonical_name": item["canonical_name"],
                    "preferred_scientific_name": item["preferred_scientific_name"],
                    "datasets": item["datasets"],
                    "score": round(item["score"], 4),
                    "match_type": item["match_type"],
                    "matched_value": item["matched_value"],
                    "match_sources": sorted(item["match_sources"]),
                    "matched_terms": item["matched_terms"],
                }
            )
        return payload

    def search(
        self,
        query: str,
        mode: str = "auto",
        *,
        dataset: str | None = None,
        limit: int = 20,
        min_score: float = 0.6,
    ) -> dict[str, Any]:
        """Search organism profiles by scientific name, canonical name, or synonym."""
        self._ensure_loaded()

        query_value = (query or "").strip()
        requested_mode = (mode or "auto").strip().lower()
        if requested_mode not in SEARCH_MODES:
            requested_mode = "auto"
        effective_mode = "name"
        bounded_limit = min(max(_safe_int(limit, 20), 1), 100)
        bounded_min_score = min(max(_safe_float(min_score, 0.6), 0.0), 1.0)
        dataset_filter = (dataset or "").strip() or None

        if not query_value:
            return {
                "query": {
                    "value": "",
                    "requested_mode": requested_mode,
                    "effective_mode": effective_mode,
                    "normalized_name": "",
                    "dataset_filter": dataset_filter,
                    "limit": bounded_limit,
                    "min_score": bounded_min_score,
                },
                "hits": [],
                "total_hits": 0,
            }

        query_normalized = normalize_name(query_value)
        hits_by_cluster: dict[str, dict[str, Any]] = {}

        for alias_row in self._alias_rows:
            cluster_id = alias_row["cluster_id"].strip()
            summary = self._cluster_summary.get(cluster_id)
            if summary is None:
                continue
            if dataset_filter and dataset_filter not in summary["datasets_set"]:
                continue

            alias_name = alias_row["alias_name"].strip()
            normalized_alias = alias_row["normalized_alias"].strip()
            score, match_type = _name_similarity(query_value, query_normalized, alias_name, normalized_alias)
            if score < bounded_min_score:
                continue
            self._add_hit(
                hits_by_cluster,
                cluster_id=cluster_id,
                score=score,
                match_type=match_type,
                matched_value=alias_name,
                source=alias_row["source_type"].strip() or "alias",
            )

        hits = self._rank_hits(hits_by_cluster, bounded_limit)
        return {
            "query": {
                "value": query_value,
                "requested_mode": requested_mode,
                "effective_mode": effective_mode,
                "normalized_name": query_normalized,
                "dataset_filter": dataset_filter,
                "limit": bounded_limit,
                "min_score": bounded_min_score,
            },
            "hits": hits,
            "total_hits": len(hits),
        }

    def lookup(
        self,
        query: str | None = None,
        *,
        scientific_name: str | None = None,
        language: str = "ja",
    ) -> dict[str, Any]:
        """Resolve one best-matched organism profile by name or synonym."""
        self._ensure_loaded()

        query_value = (scientific_name or query or "").strip()
        normalized_language = "en" if language == "en" else "ja"
        if not query_value:
            return {
                "query": {"value": "", "language": normalized_language},
                "matched": False,
                "match": None,
                "profile": None,
            }

        search_result = self.search(query=query_value, mode="name", limit=1, min_score=0.6)
        if not search_result["hits"]:
            return {
                "query": {"value": query_value, "language": normalized_language},
                "matched": False,
                "match": None,
                "profile": None,
            }

        best_hit = search_result["hits"][0]
        profile = copy.deepcopy(self._profiles_by_cluster[best_hit["cluster_id"]])
        profile["dataset_labels"] = {
            dataset_id: DATASET_LABELS[normalized_language].get(dataset_id, dataset_id)
            for dataset_id in profile.get("datasets", [])
        }
        profile["annotation_labels"] = {
            key: ANNOTATION_LABELS[normalized_language].get(key, key)
            for key in profile.get("risk_annotations", {})
        }
        return {
            "query": {"value": query_value, "language": normalized_language},
            "matched": True,
            "match": {
                "cluster_id": best_hit["cluster_id"],
                "score": best_hit["score"],
                "match_type": best_hit["match_type"],
                "matched_value": best_hit["matched_value"],
                "match_sources": best_hit["match_sources"],
            },
            "profile": profile,
        }


def get_bio_database(db_path: str | None = None) -> BioDatabase:
    """Get the singleton bio database using either a custom path or the bundled DB."""
    return BioDatabase.get_instance(db_path)


def search_organisms(
    query: str,
    mode: str = "auto",
    *,
    dataset: str | None = None,
    limit: int = 20,
    min_score: float = 0.6,
    db_path: str | None = None,
) -> dict[str, Any]:
    """Search organisms through the canonical runtime library API."""
    db = get_bio_database(db_path)
    return db.search(
        query=query,
        mode=mode,
        dataset=dataset,
        limit=limit,
        min_score=min_score,
    )


def lookup_bio_profile(
    query: str | None = None,
    *,
    scientific_name: str | None = None,
    language: str = "ja",
    db_path: str | None = None,
) -> dict[str, Any]:
    """Resolve one canonical bio profile through the public runtime API."""
    db = get_bio_database(db_path)
    return db.lookup(
        query=query,
        scientific_name=scientific_name,
        language=language,
    )


def get_bio_source_snapshots(db_path: str | None = None) -> list[dict[str, Any]]:
    """Return source snapshot metadata through a top-level helper."""
    db = get_bio_database(db_path)
    return db.get_source_snapshots()


def get_bio_runtime_status(db_path: str | None = None) -> dict[str, Any]:
    """Return lightweight runtime metadata through a top-level helper."""
    db = get_bio_database(db_path)
    return db.get_runtime_status()
