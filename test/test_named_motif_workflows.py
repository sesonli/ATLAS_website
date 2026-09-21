import csv
import json
import os
import sys
import zipfile
from pathlib import Path


os.environ.setdefault("ATLAS_SKIP_STARTUP_CLEANUP", "1")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import app as flask_app


BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "static" / "data" / "named_motifs"
FR3D_COMMIT = "ed850c00df01616e58c643b0f84bdc69662b5d55"


def load_summary(motif):
    return json.loads(
        (DATA_DIR / motif / "summary.json").read_text(encoding="utf-8")
    )


def csv_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_sarcin_ricin_release_has_reproducible_acceptance_evidence():
    summary = load_summary("sarcin_ricin")
    counts = {item["label"]: item["count"] for item in summary["stage_counts"]}
    assert counts == {
        "WC-bounded RNAdex internal-loop search pool": 224796,
        "Seven-node projected-graph candidates": 3399,
        "Candidates in the formal core sequence space": 2750,
        "High-confidence canonical SRM representatives": 3,
    }
    assert len(csv_rows(DATA_DIR / "sarcin_ricin" / "candidates.csv")) == 2750
    assert {
        item["pdb"] for item in summary["confirmed_representatives"]
    } == {"1F84", "1JBR", "1K73"}
    assert all(
        item["high_confidence_canonical_srm"]
        and item["fr3d_all_required_interactions"]
        and item["fr3d_discrepancy"] <= 0.5
        and all(item["fr3d_interaction_checks"].values())
        for item in summary["confirmed_representatives"]
    )
    validation = summary["validation"]
    assert validation["status"] == "passed"
    assert validation["fr3d_python_commit"] == FR3D_COMMIT
    assert validation["representatives_passing"] == 3
    assert validation["candidate_audit"]["sample_size"] == 24
    assert validation["candidate_audit"]["passing_all_criteria"] == 24


def test_gnra_release_rejects_graph_candidates_that_fail_full_criteria():
    summary = load_summary("gnra")
    counts = {item["label"]: item["count"] for item in summary["stage_counts"]}
    assert counts == {
        "Four-nucleotide RNAdex hairpins": 35177,
        "GNRA sequence candidates": 16948,
        "Candidates with first-to-fourth non-WC edge": 15331,
        "High-confidence canonical GNRA representatives": 3,
    }
    assert len(csv_rows(DATA_DIR / "gnra" / "candidates.csv")) == 15331
    assert {
        item["pdb"] for item in summary["confirmed_representatives"]
    } == {"2JYJ", "5X8R", "8FLD"}
    assert all(
        item["high_confidence_canonical_gnra"]
        and item["fr3d_all_required_interactions"]
        and item["fr3d_discrepancy"] <= 0.8
        and all(item["fr3d_interaction_checks"].values())
        for item in summary["confirmed_representatives"]
    )
    validation = summary["validation"]
    assert validation["status"] == "passed"
    assert validation["fr3d_python_commit"] == FR3D_COMMIT
    assert validation["representatives_passing"] == 3
    assert validation["candidate_audit"]["sample_size"] == 24
    assert validation["candidate_audit"]["passing_all_criteria"] == 3


def test_named_motif_download_archives_match_displayed_representatives():
    expected = {
        "sarcin_ricin": {
            "1F84_sarcin_ricin.pdb",
            "1JBR_sarcin_ricin.pdb",
            "1K73_sarcin_ricin.pdb",
        },
        "gnra": {"2JYJ_gnra.pdb", "5X8R_gnra.pdb", "8FLD_gnra.pdb"},
    }
    for motif, filenames in expected.items():
        motif_dir = DATA_DIR / motif
        assert {path.name for path in (motif_dir / "motifs").glob("*.pdb")} == filenames
        with zipfile.ZipFile(motif_dir / "confirmed_structures.zip") as archive:
            assert set(archive.namelist()) == filenames
        assert (motif_dir / "fr3d_validation.json").is_file()
        assert len(csv_rows(motif_dir / "candidate_audit.csv")) == 24


def test_worked_examples_index_renders_all_three_workflows():
    flask_app.app.config["TESTING"] = True
    with flask_app.app.test_client() as client:
        response = client.get("/examples")

    assert response.status_code == 200
    assert b"Kink-turn" in response.data
    assert b"Sarcin-ricin motif" in response.data
    assert b"GNRA tetraloop" in response.data
    assert b"3 FR3D-validated representatives" in response.data
    assert b"Candidate counts are" in response.data


def test_sarcin_ricin_page_exposes_definition_validation_and_files():
    flask_app.app.config["TESTING"] = True
    with flask_app.app.test_client() as client:
        response = client.get("/examples/sarcin-ricin")

    assert response.status_code == 200
    assert b"Sarcin-ricin motif" in response.data
    assert b"Candidate boundary" in response.data
    assert b"All four formal Leontis-Westhof base pairs" in response.data
    assert b"24/24" in response.data
    assert b"discrepancy 0.437" in response.data
    assert b"1F84_sarcin_ricin.pdb" in response.data
    assert b"fr3d_validation.json" in response.data
    assert b"10.1093/nar/gkm069" in response.data


def test_gnra_page_exposes_interaction_signature_and_new_representatives():
    flask_app.app.config["TESTING"] = True
    with flask_app.app.test_client() as client:
        response = client.get("/examples/gnra")

    assert response.status_code == 200
    assert b"GNRA tetraloop" in response.data
    assert b"first-to-fourth G-A tSH pair" in response.data
    assert b"stricter than the published broad FR3D GNRA search" in response.data
    assert b"reproduce the published FR3D GNRA search" not in response.data
    assert b"3/24" in response.data
    for pdb in (b"2JYJ", b"5X8R", b"8FLD"):
        assert pdb in response.data
    for rejected_old_example in (b"1MZP", b"1ZIH", b"2GIS"):
        assert rejected_old_example not in response.data
    assert b"candidate set" in response.data
    assert b"10.1007/s00285-007-0110-x" in response.data


def test_unknown_named_motif_page_returns_404():
    flask_app.app.config["TESTING"] = True
    with flask_app.app.test_client() as client:
        response = client.get("/examples/not-a-motif")
    assert response.status_code == 404
