import os
import sys
import json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import app as flask_app


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_kturn_workflow_data_counts():
    workflow = flask_app.load_kturn_workflow_data()

    counts = {item['label']: item['count'] for item in workflow['stage_counts']}
    assert counts['Search pool'] == 290864
    assert counts['Stage 1 graph-topology candidates'] == 70701
    assert counts['Candidates with G·A pairs'] == 11240
    assert counts['Representatives shown (of many candidates)'] == 4

    pdbs = {item['pdb'] for item in workflow['confirmed_representatives']}
    assert pdbs == {'2GIS', '1E7K', '4BW0', '2OZB'}


def test_kturn_workflow_summary_controls_are_not_overclaimed():
    data_dir = os.path.join(BASE_DIR, 'static', 'data', 'kturn_workflow')
    with open(os.path.join(data_dir, 'caseC_kturn_summary.json'), 'r', encoding='utf-8') as handle:
        summary = json.load(handle)

    descriptions = {
        item['pdb']: item['description']
        for item in summary['representatives']
    }
    assert descriptions['4C4W'] == 'rare non-standard Kt-23 k-turn outside the standard tandem G.A confirmation set'
    assert descriptions['4LCK'] == 'T-box riboswitch/tRNA G.A-topology candidate rejected by the standard geometry rule'


def test_kturn_workflow_static_files_exist():
    data_dir = os.path.join(BASE_DIR, 'static', 'data', 'kturn_workflow')
    expected = {
        'caseC_kturn_summary.json',
        'caseC_geometry_confirmed.json',
        'caseC_kturn_hits.csv',
        'caseC_kturn_figure.png',
        'caseC_kturn_figure.pdf',
        'confirmed_kink_turn_structures.zip',
    }
    assert expected <= set(os.listdir(data_dir))
    motif_dir = os.path.join(data_dir, 'motifs')
    for pdb in ('2GIS', '1E7K', '4BW0', '2OZB'):
        assert os.path.exists(os.path.join(motif_dir, f'{pdb}_kink_turn.pdb'))


def test_kturn_workflow_page_renders():
    flask_app.app.config['TESTING'] = True
    with flask_app.app.test_client() as client:
        response = client.get('/examples/kink-turn')

    assert response.status_code == 200
    assert b'Kink-turn Worked Example' in response.data
    assert b'Stage 1 graph-topology candidates' in response.data
    assert b'Geometry-confirmed Representatives' in response.data
    assert b'2GIS' in response.data
    assert b'1E7K' in response.data
    assert b'4BW0' in response.data
    assert b'2OZB' in response.data
    assert b'isolated kink-turn RNA' not in response.data
    assert b'box C/D RNA kink-turn' not in response.data
    assert b'<td>riboswitch kink-turn</td>' not in response.data
    assert b'HMKt-7 kink-turn bound by L7Ae protein' in response.data
    assert b'U4 snRNA kink-turn within the human Prp31-15.5K-U4 snRNA complex' in response.data
    # found motif structures are downloadable
    assert b'Download found motif structures' in response.data
    assert b'confirmed_kink_turn_structures.zip' in response.data
    assert b'2GIS_kink_turn.pdb' in response.data


def test_custom_search_page_has_no_kturn_template():
    # The live drawing search cannot reproduce the curated kink-turn funnel
    # (it matches a sample of full-structure graphs, not the internal/bulge
    # motif library), so the kink-turn template was removed from this page.
    flask_app.app.config['TESTING'] = True
    with flask_app.app.test_client() as client:
        response = client.get('/custom-motif-search')

    assert response.status_code == 200
    assert b'Load Kink-turn Candidate Template' not in response.data
    assert b'Worked Example: Kink-turn Retrieval Workflow' not in response.data


def test_homepage_links_worked_examples():
    flask_app.app.config['TESTING'] = True
    with flask_app.app.test_client() as client:
        response = client.get('/')

    assert response.status_code == 200
    assert b'VIEW WORKED EXAMPLES' in response.data
    assert b'/examples' in response.data
