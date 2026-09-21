"""Acceptance checks for the staged ATLAS2025 v1.1 database release."""

import hashlib
import json
import os
import sqlite3

import app as flask_app


EXPECTED_SHA256 = (
    "3d2356b6f82016d3d5a927a0d141a5803cb1ec45266c1d26e2607d075aecc9ed"
)


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def test_database_path_points_to_current_database():
    assert os.path.isfile(flask_app.ATLAS_DB_PATH)


def test_v1_1_database_identity_and_integrity():
    assert _sha256(flask_app.ATLAS_DB_PATH) == EXPECTED_SHA256
    connection = sqlite3.connect(
        f"file:{flask_app.ATLAS_DB_PATH}?mode=ro", uri=True
    )
    try:
        assert connection.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
        assert connection.execute("SELECT COUNT(*) FROM data").fetchone()[0] == 455264
        assert connection.execute("SELECT COUNT(*) FROM PK").fetchone()[0] == 7491
        data_pdb = {
            row[0].upper()
            for row in connection.execute(
                "SELECT DISTINCT pdbid FROM data WHERE pdbid IS NOT NULL"
            )
        }
        pk_pdb = {
            row[0].upper()
            for row in connection.execute(
                "SELECT DISTINCT pdbid FROM PK WHERE pdbid IS NOT NULL"
            )
        }
        assert len(data_pdb | pk_pdb) == 5846
    finally:
        connection.close()


def test_download_route_resolves_without_reading_multi_gigabyte_body():
    flask_app.app.config["TESTING"] = True
    with flask_app.app.test_client() as client:
        current = client.get("/download-database", buffered=False)
        try:
            assert current.status_code == 200
            assert "RNAdex.db" in current.headers["Content-Disposition"]
        finally:
            current.close()


def test_home_and_user_guide_expose_release_values():
    flask_app.app.config["TESTING"] = True
    with flask_app.app.test_client() as client:
        home = client.get("/")
        guide = client.get("/user-guide")
    assert home.status_code == 200
    assert guide.status_code == 200
    for response in (home, guide):
        assert b"462,755" in response.data
        assert b"5,846" in response.data
        assert b"RNAdex" in response.data
    assert b"MIT License" in home.data
    assert b"CC BY 4.0" in home.data


def test_v1_1_standard_and_pk_queries_return_results():
    hairpins, hairpin_metadata = flask_app.search_motif(
        "hairpin", 6, max_results=5
    )
    pseudoknots, pk_metadata = flask_app.search_motif(
        "pseudoknot_H", max_results=5
    )
    assert hairpins
    assert hairpin_metadata["result_count"] == len(hairpins)
    assert pseudoknots
    assert pk_metadata["result_count"] == len(pseudoknots)


def test_v1_1_kturn_evidence_is_consistent():
    workflow = flask_app.load_kturn_workflow_data()
    counts = {item["label"]: item["count"] for item in workflow["stage_counts"]}
    assert counts == {
        "Search pool": 290864,
        "Stage 1 graph-topology candidates": 70701,
        "Candidates with G·A pairs": 11240,
        "Representatives shown (of many candidates)": 4,
    }
    verification_path = os.path.join(
        flask_app.KTURN_WORKFLOW_DATA_DIR, "caseC_release_verification.json"
    )
    with open(verification_path, encoding="utf-8") as handle:
        verification = json.load(handle)
    assert verification["database_sha256"] == EXPECTED_SHA256
