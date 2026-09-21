"""
Comprehensive project health checks.
Run: pytest test/test_project_health.py -v -s
"""
import sys, os, json, sqlite3
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import app as flask_app

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def client():
    flask_app.app.config['TESTING'] = True
    with flask_app.app.test_client() as c:
        yield c


# ── Standard search ──────────────────────────────────────────────────────────

class TestStandardSearch:
    def test_hairpin_returns_200(self, client):
        r = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
        assert r.status_code == 200

    def test_hairpin_result_tuple_has_8_elements(self, client):
        """result tuple: id, motif_type, pdbid, nt_number, filecontent, sequence, quality_ok, warning"""
        results, metadata = flask_app.search_motif('hairpin', 6, max_results=5)
        assert len(results) > 0
        for row in results:
            assert len(row) == 8, f"Expected 8 elements, got {len(row)}"

    def test_sequence_column_is_string_of_bases(self, client):
        results, _ = flask_app.search_motif('hairpin', 6, max_results=20)
        for row in results:
            seq = row[5]
            assert isinstance(seq, str)
            assert len(seq) > 0
            # Should be A/G/C/U or N/A
            valid = set('AGCU?[]NA/')
            bad = [c for c in seq if c not in valid]
            assert not bad, f"Unexpected chars {bad} in sequence '{seq}'"

    def test_truncation_warning_shown(self, client):
        r = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
        assert b'matching results' in r.data

    def test_session_id_in_response(self, client):
        r = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
        assert b'session_id' in r.data or b'download_csv' in r.data

    def test_sequence_column_in_results_page(self, client):
        r = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
        assert b'Sequence' in r.data

    def test_internal_search(self, client):
        r = client.post('/search', data={'motif_type': 'internal', 'nt_number': '6'})
        assert r.status_code == 200

    def test_junction_search(self, client):
        r = client.post('/search', data={'motif_type': 'multiway-junction', 'nt_number': '3-way junction'})
        assert r.status_code == 200

    def test_no_results_page(self, client):
        r = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '999'})
        assert r.status_code == 200
        assert b'No Results' in r.data


# ── Download routes ───────────────────────────────────────────────────────────

class TestDownloads:
    def _do_search(self, client):
        r = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
        html = r.data.decode()
        # Extract session_id from href
        import re
        m = re.search(r'session_id=([\w-]+)', html)
        return m.group(1) if m else None

    def test_csv_current_requires_valid_session(self, client):
        r = client.get('/download_csv?session_id=fake-id')
        assert r.status_code == 404

    def test_zip_current_requires_valid_session(self, client):
        r = client.get('/download_zip?session_id=fake-id')
        assert r.status_code == 404

    def test_csv_current_after_search(self, client):
        session_id = self._do_search(client)
        assert session_id is not None, "No session_id found in response"
        r = client.get(f'/download_csv?session_id={session_id}')
        assert r.status_code == 200
        assert b'Sequence' in r.data  # sequence column present

    def test_zip_current_after_search(self, client):
        session_id = self._do_search(client)
        assert session_id is not None
        r = client.get(f'/download_zip?session_id={session_id}')
        assert r.status_code == 200
        assert r.content_type == 'application/zip'

    def test_csv_full_has_sequence_column(self, client):
        r = client.get('/download_csv_full?motif_type=hairpin&nt_number=6')
        assert r.status_code == 200
        data = b''.join(r.response)
        assert b'Sequence' in data

    def test_csv_full_missing_param(self, client):
        r = client.get('/download_csv_full')
        assert r.status_code == 400

    def test_zip_full_missing_param(self, client):
        r = client.get('/download_zip_full')
        assert r.status_code == 400

    def test_zip_full_unknown_type(self, client):
        r = client.get('/download_zip_full?motif_type=unknown')
        assert r.status_code == 400


# ── Custom search templates ───────────────────────────────────────────────────

class TestCustomSearchTemplates:
    SIMPLE_GRAPH = json.dumps({
        "nodes": [
            {"id": 1, "x": 200, "y": 200, "name": "N1"},
            {"id": 2, "x": 300, "y": 200, "name": "N2"}
        ],
        "edges": [
            {"from": 1, "to": 2, "type": "covalent", "attribute": [0, 0, 1]}
        ]
    })

    def test_custom_results_shows_warning(self, client):
        r = client.post('/process-custom-motif', data={'graph_data': self.SIMPLE_GRAPH})
        assert r.status_code == 200
        assert b'subset of the RNA structure database' in r.data

    def test_custom_results_or_no_results_shown(self, client):
        r = client.post('/process-custom-motif', data={'graph_data': self.SIMPLE_GRAPH})
        assert r.status_code == 200
        # Either found results or no-results page
        assert b'Custom Motif Search Results' in r.data or b'No Results Found' in r.data

    def test_custom_no_results_has_warning(self, client):
        """custom_no_results.html must contain the subset warning."""
        tmpl_path = os.path.join(BASE_DIR, 'templates', 'custom_no_results.html')
        assert os.path.exists(tmpl_path), "custom_no_results.html missing"
        content = open(tmpl_path).read()
        assert 'subset of the RNA structure database' in content

    def test_no_results_html_has_no_custom_warning(self):
        """Generic no_results.html must NOT contain the custom search warning."""
        tmpl_path = os.path.join(BASE_DIR, 'templates', 'no_results.html')
        content = open(tmpl_path).read()
        assert 'subset of the RNA structure database' not in content, \
            "Generic no_results.html should not contain custom search warning"

    def test_custom_results_template_has_canvas(self):
        tmpl_path = os.path.join(BASE_DIR, 'templates', 'custom_results.html')
        content = open(tmpl_path).read()
        assert 'motif-canvas' in content

    def test_custom_no_results_template_has_canvas(self):
        tmpl_path = os.path.join(BASE_DIR, 'templates', 'custom_no_results.html')
        content = open(tmpl_path).read()
        assert 'motif-canvas' in content

    def test_custom_results_column_header_not_misleading(self):
        """Column headers must match actual content: nt_number -> 'Full Nucleotide List'."""
        tmpl_path = os.path.join(BASE_DIR, 'templates', 'custom_results.html')
        content = open(tmpl_path).read()
        assert 'Full Nucleotide Sequence' not in content, \
            "Column header 'Full Nucleotide Sequence' is misleading — content is nt_number not sequence"
        assert 'Full Nucleotide List' in content


# ── Database integrity ────────────────────────────────────────────────────────

class TestDatabaseIntegrity:
    def test_atlas_db_exists(self):
        assert os.path.exists(os.path.join(BASE_DIR, 'ATLAS.db'))

    def test_atlas_db_has_data_table(self):
        conn = sqlite3.connect(os.path.join(BASE_DIR, 'ATLAS.db'))
        c = conn.cursor()
        c.execute("SELECT COUNT(*) FROM data")
        count = c.fetchone()[0]
        conn.close()
        assert count > 100000, f"Expected >100k rows, got {count}"

    def test_atlas_db_has_pk_table(self):
        conn = sqlite3.connect(os.path.join(BASE_DIR, 'ATLAS.db'))
        c = conn.cursor()
        c.execute("SELECT COUNT(*) FROM PK")
        count = c.fetchone()[0]
        conn.close()
        assert count > 0

    def test_batch_pickle_exists(self):
        assert os.path.exists(os.path.join(BASE_DIR, 'batch_0000_graphs.pickle'))


# ── total_matched accuracy ────────────────────────────────────────────────────

class TestMetadata:
    def test_total_matched_less_than_total_scanned(self):
        _, metadata = flask_app.search_motif('hairpin', 6, max_results=100)
        assert metadata['total_matched'] <= metadata['total_scanned']

    def test_total_matched_equals_result_count_when_not_truncated(self):
        # Use a rare nt_number that won't hit 10k limit
        results, metadata = flask_app.search_motif('hairpin', 19, max_results=10000)
        if not metadata.get('truncated'):
            assert metadata['result_count'] == len(results)

    def test_warning_message_uses_total_matched(self, client):
        r = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
        html = r.data.decode()
        if 'Showing' in html:
            assert 'matching results' in html

    def test_total_matched_does_not_depend_on_the_display_cap(self):
        """The truncation notice quotes total_matched as the size of the whole
        result set, so it must describe the database, not however many rows
        happened to be seen before collection stopped at the cap."""
        _, small = flask_app.search_motif('bulge', 7, max_results=50)
        _, large = flask_app.search_motif('bulge', 7, max_results=500)
        assert small['truncated'] and large['truncated']
        assert small['total_matched'] == large['total_matched']

    def test_total_matched_equals_the_full_download_row_count(self, client):
        """'Download CSV/ZIP (All Results) for complete data' promises exactly
        total_matched rows; the two counts have to agree."""
        _, metadata = flask_app.search_motif('bulge', 7, max_results=50)
        assert metadata['truncated']
        response = client.get('/download_csv_full?motif_type=bulge&nt_number=7')
        rows = response.get_data(as_text=True).strip().splitlines()[1:]
        assert metadata['total_matched'] == len(rows)
