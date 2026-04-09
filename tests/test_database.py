"""Tests for the public `ra_bio` runtime package."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

from ra_bio import (
    BioDatabase,
    get_bio_database,
    get_bio_runtime_status,
    get_bio_source_snapshots,
    lookup_bio_profile,
    search_organisms,
)


def _write_fixture_db(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path)
    try:
        connection.executescript(
            """
            CREATE TABLE source_snapshots (
                dataset_id TEXT PRIMARY KEY,
                source_url TEXT NOT NULL,
                source_filename TEXT NOT NULL,
                source_version TEXT NOT NULL,
                fetched_at TEXT NOT NULL,
                content_hash TEXT NOT NULL,
                row_count INTEGER NOT NULL
            );

            CREATE TABLE organisms (
                cluster_id TEXT PRIMARY KEY,
                canonical_name TEXT NOT NULL,
                canonical_name_normalized TEXT NOT NULL,
                preferred_scientific_name TEXT NOT NULL,
                preferred_scientific_name_normalized TEXT NOT NULL,
                datasets_json TEXT NOT NULL,
                scientific_names_json TEXT NOT NULL,
                aliases_json TEXT NOT NULL,
                profile_json TEXT NOT NULL
            );

            CREATE TABLE organism_aliases (
                cluster_id TEXT NOT NULL,
                alias_name TEXT NOT NULL,
                normalized_alias TEXT NOT NULL,
                source_type TEXT NOT NULL,
                source_ref TEXT NOT NULL
            );
            """
        )

        connection.execute(
            """
            INSERT INTO source_snapshots
            (dataset_id, source_url, source_filename, source_version, fetched_at, content_hash, row_count)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "fungi",
                "https://example.invalid/fungi",
                "risk_fungi_fixture.xlsx",
                "fixture",
                "2026-04-10T00:00:00Z",
                "fixture-hash",
                1,
            ),
        )

        profile = {
            "cluster_id": "ORG-000001",
            "canonical_name": "Actinomortierella wolfii",
            "preferred_scientific_name": "Actinomortierella wolfii",
            "scientific_names": ["Actinomortierella wolfii"],
            "aliases": ["Mortierella wolfii"],
            "datasets": ["fungi"],
            "dataset_versions": {"fungi": "fixture"},
            "canonical_statuses": ["canonical"],
            "hosts": [],
            "diseases": [],
            "risk_annotations": {
                "trba": [
                    {
                        "dataset": "fungi",
                        "value": "‡",
                        "source_label": "11)TRBA 460 [2023-12]",
                        "source_record_id": "fungi:7",
                    }
                ]
            },
            "source_records": [
                {
                    "entry_id": "fungi:7",
                    "dataset": "fungi",
                    "scientific_name": "Actinomortierella wolfii",
                    "canonical_name": "Actinomortierella wolfii",
                    "canonical_status": "canonical",
                    "aliases": ["Mortierella wolfii"],
                }
            ],
        }

        connection.execute(
            """
            INSERT INTO organisms
            (cluster_id, canonical_name, canonical_name_normalized, preferred_scientific_name,
             preferred_scientific_name_normalized, datasets_json, scientific_names_json,
             aliases_json, profile_json)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "ORG-000001",
                "Actinomortierella wolfii",
                "actinomortierellawolfii",
                "Actinomortierella wolfii",
                "actinomortierellawolfii",
                json.dumps(["fungi"], ensure_ascii=False),
                json.dumps(["Actinomortierella wolfii"], ensure_ascii=False),
                json.dumps(["Mortierella wolfii"], ensure_ascii=False),
                json.dumps(profile, ensure_ascii=False),
            ),
        )

        alias_rows = [
            ("ORG-000001", "Actinomortierella wolfii", "actinomortierellawolfii", "scientific_name", "fixture"),
            ("ORG-000001", "Actinomortierella wolfii", "actinomortierellawolfii", "canonical_name", "fixture"),
            ("ORG-000001", "Mortierella wolfii", "mortierellawolfii", "alias", "fixture"),
        ]
        connection.executemany(
            """
            INSERT INTO organism_aliases
            (cluster_id, alias_name, normalized_alias, source_type, source_ref)
            VALUES (?, ?, ?, ?, ?)
            """,
            alias_rows,
        )

        connection.commit()
    finally:
        connection.close()


def test_search_finds_synonym(tmp_path):
    """Fuzzy search should resolve a historical synonym."""
    db_path = tmp_path / "bio.sqlite3"
    _write_fixture_db(db_path)

    BioDatabase.reset_instance()
    db = get_bio_database(str(db_path))
    result = db.search(query="Mortierella wolfi", limit=5, min_score=0.6)

    assert result["hits"]
    assert result["hits"][0]["canonical_name"] == "Actinomortierella wolfii"
    assert "alias" in result["hits"][0]["match_sources"]


def test_lookup_returns_profile_and_labels(tmp_path):
    """Lookup should return the aggregated profile plus localized labels."""
    db_path = tmp_path / "bio.sqlite3"
    _write_fixture_db(db_path)

    BioDatabase.reset_instance()
    db = get_bio_database(str(db_path))
    result = db.lookup(query="Mortierella wolfii", language="ja")

    assert result["matched"] is True
    assert result["profile"]["canonical_name"] == "Actinomortierella wolfii"
    assert result["profile"]["dataset_labels"]["fungi"] == "真菌"
    assert result["profile"]["annotation_labels"]["trba"] == "TRBAリスクグループ"


def test_bundled_database_is_loadable():
    """The generated bundled SQLite database should be readable without a custom path."""
    BioDatabase.reset_instance()
    db = get_bio_database()
    result = db.search(query="Mortierella wolfii", limit=5, min_score=0.6)

    assert result["hits"], "expected generated bundled dataset to be present"


def test_source_snapshots_are_exposed(tmp_path):
    """Source snapshot metadata should be available through the runtime API."""
    db_path = tmp_path / "bio.sqlite3"
    _write_fixture_db(db_path)

    BioDatabase.reset_instance()
    db = get_bio_database(str(db_path))
    snapshots = db.get_source_snapshots()

    assert snapshots[0]["dataset_id"] == "fungi"
    assert snapshots[0]["source_version"] == "fixture"


def test_top_level_helpers_delegate_to_database(tmp_path):
    """Top-level helper functions should expose the same canonical library behavior."""
    db_path = tmp_path / "bio.sqlite3"
    _write_fixture_db(db_path)

    BioDatabase.reset_instance()
    search_result = search_organisms(query="Mortierella wolfii", db_path=str(db_path))
    lookup_result = lookup_bio_profile(query="Mortierella wolfii", db_path=str(db_path))
    snapshots = get_bio_source_snapshots(str(db_path))
    status = get_bio_runtime_status(str(db_path))

    assert search_result["hits"][0]["canonical_name"] == "Actinomortierella wolfii"
    assert lookup_result["profile"]["canonical_name"] == "Actinomortierella wolfii"
    assert snapshots[0]["dataset_id"] == "fungi"
    assert status["status"] == "ok"
    assert status["organism_count"] == 1
