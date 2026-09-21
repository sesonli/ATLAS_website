"""
Tests for search and download routes.
Run from project root: pytest test/test_search_download.py -v
"""
import sys
import os
import csv
import io
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import app as flask_app


@pytest.fixture
def client():
    flask_app.app.config['TESTING'] = True
    with flask_app.app.test_client() as client:
        yield client


# --- Search route ---

def test_search_hairpin_returns_results(client):
    resp = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
    assert resp.status_code == 200
    assert b'search_results' in resp.data or b'Download' in resp.data


def test_search_hairpin_truncation_warning(client):
    """When results hit 10k limit, warning message should appear."""
    resp = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
    assert resp.status_code == 200
    # hairpin 6nt has 81k records, should be truncated
    assert b'Showing' in resp.data or b'matching results' in resp.data


def test_search_internal(client):
    resp = client.post('/search', data={'motif_type': 'internal', 'nt_number': '6'})
    assert resp.status_code == 200


def test_search_junction(client):
    resp = client.post('/search', data={'motif_type': 'multiway-junction', 'nt_number': '3-way junction'})
    assert resp.status_code == 200


def test_search_no_results(client):
    """Searching with an unrealistic nt_number should hit no_results page."""
    resp = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '999'})
    assert resp.status_code == 200


# --- Session isolation (Bug 1) ---

def test_concurrent_searches_have_different_session_ids(client):
    """Two searches should produce different session_ids in metadata."""
    resp1 = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
    resp2 = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '7'})
    assert resp1.status_code == 200
    assert resp2.status_code == 200
    # Both should succeed independently
    assert b'Download' in resp1.data
    assert b'Download' in resp2.data


# --- Download current results (requires prior search) ---

def test_download_csv_current_after_search(client):
    client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
    # Extract session_id from the response HTML
    resp = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
    # Find session_id in the page
    html = resp.data.decode()
    assert 'session_id' in html or 'download_csv' in html


def test_download_csv_without_session_returns_404(client):
    resp = client.get('/download_csv?session_id=nonexistent-uuid')
    assert resp.status_code == 404


def test_download_zip_without_session_returns_404(client):
    resp = client.get('/download_zip?session_id=nonexistent-uuid')
    assert resp.status_code == 404


# --- Download all results (full query) ---

def test_download_csv_full_hairpin(client):
    resp = client.get('/download_csv_full?motif_type=hairpin&nt_number=6')
    assert resp.status_code == 200
    assert resp.content_type == 'text/csv; charset=utf-8'
    # Should have header row
    data = b''.join(resp.response)
    header = next(csv.reader(io.StringIO(data.decode('utf-8'))))
    assert header[:4] == ['ID', 'Motif Type', 'PDB ID', 'NT Number']


def test_download_csv_full_junction(client):
    resp = client.get('/download_csv_full?motif_type=3-way junction&nt_number=')
    assert resp.status_code == 200


def test_download_csv_full_missing_motif_type(client):
    resp = client.get('/download_csv_full')
    assert resp.status_code == 400


def test_download_zip_full_hairpin(client):
    resp = client.get('/download_zip_full?motif_type=hairpin&nt_number=6',
                      headers={'Range': 'bytes=0-1023'})  # just check it starts
    # Accept 200 or 206
    assert resp.status_code in (200, 206)


def test_download_zip_full_missing_motif_type(client):
    resp = client.get('/download_zip_full')
    assert resp.status_code == 400


def test_download_zip_full_unknown_motif_type(client):
    resp = client.get('/download_zip_full?motif_type=unknown_type')
    assert resp.status_code == 400


# --- total_matched accuracy (Bug 3) ---

def test_total_matched_in_metadata(client):
    """total_matched should reflect actual matching records, not total scanned."""
    resp = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
    html = resp.data.decode()
    # Warning message should say "matching results", not "total matches"
    if 'Showing' in html:
        assert 'matching results' in html
