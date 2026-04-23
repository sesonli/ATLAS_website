"""
Extended coverage tests for previously untested paths.
Run: pytest test/test_extended_coverage.py -v -s
"""
import sys, os, sqlite3, zipfile, io
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import app as flask_app

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def client():
    flask_app.app.config['TESTING'] = True
    with flask_app.app.test_client() as c:
        yield c


# ── Bulge search ──────────────────────────────────────────────────────────────

class TestBulgeSearch:
    def test_bulge_returns_results(self, client):
        # UI sends display_value + 4, so bulge "1" -> nt_number=5
        r = client.post('/search', data={'motif_type': 'bulge', 'nt_number': '5'})
        assert r.status_code == 200
        assert b'Download' in r.data

    def test_bulge_result_tuple_has_8_elements(self):
        results, metadata = flask_app.search_motif('bulge', 5, max_results=5)
        assert len(results) > 0
        for row in results:
            assert len(row) == 8

    def test_bulge_sequence_column_present(self, client):
        r = client.post('/search', data={'motif_type': 'bulge', 'nt_number': '5'})
        assert b'Sequence' in r.data

    def test_bulge_csv_full(self, client):
        r = client.get('/download_csv_full?motif_type=bulge&nt_number=5')
        assert r.status_code == 200
        data = b''.join(r.response)
        assert b'Sequence' in data
        assert len(data) > 100


# ── Pseudoknot search ─────────────────────────────────────────────────────────

class TestPseudoknotSearch:
    def test_pseudoknot_H_returns_results(self, client):
        # pseudoknot type H
        r = client.post('/search', data={'motif_type': 'pseudoknot', 'nt_number': 'H'})
        assert r.status_code == 200

    def test_pseudoknot_result_tuple_has_8_elements(self):
        results, metadata = flask_app.search_motif('pseudoknot_H', max_results=5)
        if len(results) == 0:
            pytest.skip("No pseudoknot_H results in DB")
        for row in results:
            assert len(row) == 8

    def test_pseudoknot_csv_full(self, client):
        r = client.get('/download_csv_full?motif_type=pseudoknot_H&nt_number=')
        assert r.status_code == 200

    def test_pseudoknot_zip_full(self, client):
        r = client.get('/download_zip_full?motif_type=pseudoknot_H&nt_number=')
        assert r.status_code == 200


# ── Single PDB download ───────────────────────────────────────────────────────

class TestSinglePDBDownload:
    def _get_first_id(self, motif_type='hairpin'):
        conn = sqlite3.connect(os.path.join(BASE_DIR, 'ATLAS.db'))
        c = conn.cursor()
        c.execute("SELECT id FROM data WHERE motif_type=? LIMIT 1", (motif_type,))
        row = c.fetchone()
        conn.close()
        return row[0] if row else None

    def _get_first_pk_rowid(self):
        conn = sqlite3.connect(os.path.join(BASE_DIR, 'ATLAS.db'))
        c = conn.cursor()
        c.execute("SELECT rowid FROM PK LIMIT 1")
        row = c.fetchone()
        conn.close()
        return row[0] if row else None

    def test_download_single_standard_pdb(self, client):
        motif_id = self._get_first_id('hairpin')
        assert motif_id is not None
        r = client.get(f'/download_single_pdb/{motif_id}?type=standard&pdb_id=test')
        assert r.status_code == 200
        assert b'ATOM' in r.data

    def test_download_single_pseudoknot_pdb(self, client):
        rowid = self._get_first_pk_rowid()
        if rowid is None:
            pytest.skip("No PK records")
        r = client.get(f'/download_single_pdb/{rowid}?type=pseudoknot&pdb_id=test')
        assert r.status_code == 200
        assert b'ATOM' in r.data

    def test_download_single_pdb_not_found(self, client):
        r = client.get('/download_single_pdb/999999999?type=standard')
        assert r.status_code == 404

    def test_download_single_pdb_filename_format(self, client):
        motif_id = self._get_first_id('hairpin')
        r = client.get(f'/download_single_pdb/{motif_id}?type=standard&pdb_id=test')
        assert r.status_code == 200
        cd = r.headers.get('Content-Disposition', '')
        assert '.pdb' in cd


# ── ZIP content validation ────────────────────────────────────────────────────

class TestZIPContent:
    def _get_session_id(self, client):
        import re
        r = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
        m = re.search(r'session_id=([\w-]+)', r.data.decode())
        return m.group(1) if m else None

    def test_zip_current_contains_pdb_files(self, client):
        session_id = self._get_session_id(client)
        assert session_id
        r = client.get(f'/download_zip?session_id={session_id}')
        assert r.status_code == 200
        zf = zipfile.ZipFile(io.BytesIO(r.data))
        names = zf.namelist()
        assert len(names) > 0
        assert all(n.endswith('.pdb') for n in names)

    def test_zip_current_pdb_files_have_atom_records(self, client):
        session_id = self._get_session_id(client)
        r = client.get(f'/download_zip?session_id={session_id}')
        zf = zipfile.ZipFile(io.BytesIO(r.data))
        # Check first file has ATOM records
        first = zf.namelist()[0]
        content = zf.read(first).decode('utf-8', errors='ignore')
        assert 'ATOM' in content


# ── Session cleanup ───────────────────────────────────────────────────────────

class TestSessionCleanup:
    def test_search_sessions_dir_created(self, client):
        client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
        sessions_dir = os.path.join(BASE_DIR, 'search_sessions')
        assert os.path.exists(sessions_dir)

    def test_session_dir_contains_csv_and_pdb(self, client):
        import re
        r = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
        m = re.search(r'session_id=([\w-]+)', r.data.decode())
        assert m
        session_id = m.group(1)
        session_dir = os.path.join(BASE_DIR, 'search_sessions', session_id)
        assert os.path.exists(session_dir)
        files = os.listdir(session_dir)
        assert any(f.endswith('.csv') for f in files)
        assert any(f.endswith('.pdb') for f in files)


# ── CSV content validation ────────────────────────────────────────────────────

class TestCSVContent:
    def _get_session_id(self, client):
        import re
        r = client.post('/search', data={'motif_type': 'hairpin', 'nt_number': '6'})
        m = re.search(r'session_id=([\w-]+)', r.data.decode())
        return m.group(1) if m else None

    def test_csv_current_has_correct_columns(self, client):
        session_id = self._get_session_id(client)
        r = client.get(f'/download_csv?session_id={session_id}')
        assert r.status_code == 200
        lines = r.data.decode().split('\n')
        header = lines[0]
        assert 'ID' in header
        assert 'Sequence' in header
        assert 'Quality_OK' in header
        assert 'NT Number' in header

    def test_csv_current_has_data_rows(self, client):
        session_id = self._get_session_id(client)
        r = client.get(f'/download_csv?session_id={session_id}')
        lines = [l for l in r.data.decode().split('\n') if l.strip()]
        assert len(lines) > 1  # header + at least one data row

    def test_csv_full_sequence_values_are_bases(self, client):
        import csv as csvmod
        r = client.get('/download_csv_full?motif_type=hairpin&nt_number=6')
        data = b''.join(r.response).decode()
        reader = csvmod.reader(data.splitlines())
        header = next(reader)
        seq_idx = header.index('Sequence')
        for i, row in enumerate(reader):
            if i >= 5:
                break
            if len(row) > seq_idx:
                seq = row[seq_idx]
                valid = set('AGCU?[]NA/')
                bad = [c for c in seq if c not in valid]
                assert not bad, f"Invalid chars {bad} in sequence '{seq}'"


# ── Route existence ───────────────────────────────────────────────────────────

class TestRoutes:
    def test_homepage(self, client):
        r = client.get('/')
        assert r.status_code == 200

    def test_user_guide(self, client):
        r = client.get('/user-guide')
        assert r.status_code == 200

    def test_custom_motif_search_page(self, client):
        r = client.get('/custom-motif-search')
        assert r.status_code == 200

    def test_download_database_route_exists(self, client):
        # Just check route exists (don't download 8GB)
        r = client.get('/download-database')
        # 200 if file exists, 404 if not — either is acceptable, just not 500
        assert r.status_code in (200, 404)
