"""Regression checks for database-backed custom motif query examples."""

import pickle
import sqlite3

import networkx as nx
import pytest
from networkx.algorithms import isomorphism

import app as flask_app
import find_hairpin


EXPECTED_EXAMPLES = {
    "bulge-1x0": {
        "title": "Bulge: 1 × 0 loop",
        "source_pdb": "1AQ3",
        "nodes": 5,
        "edges": 5,
    },
    "internal-loop-1x1": {
        "title": "Internal loop: 1 × 1",
        "source_pdb": "17RA",
        "nodes": 6,
        "edges": 7,
    },
    "internal-loop-2x2": {
        "title": "Internal loop: 2 × 2",
        "source_pdb": "1A4D",
        "nodes": 8,
        "edges": 10,
    },
}


def as_networkx_graph(example):
    graph = nx.Graph()
    for node in example["nodes"]:
        graph.add_node(node["name"])
    names = {node["id"]: node["name"] for node in example["nodes"]}
    for edge in example["edges"]:
        graph.add_edge(
            names[edge["from"]],
            names[edge["to"]],
            attribute=edge["attribute"],
        )
    return graph


def test_custom_search_page_lists_database_backed_examples():
    flask_app.app.config["TESTING"] = True
    with flask_app.app.test_client() as client:
        response = client.get("/custom-motif-search")

    assert response.status_code == 200
    html = response.get_data(as_text=True)
    assert html.count('class="example-card"') == 5
    assert "Adjacent non-WC pairs" in html
    assert "Hairpin with non-WC pair" in html
    assert "Source topology:" not in html
    for key, expected in EXPECTED_EXAMPLES.items():
        assert f'id="{key}-button"'.encode() in response.data
        assert expected["title"] in html
        assert f'{expected["nodes"]} nt · {expected["edges"]} edges' in html


def as_posted_payload(example):
    """Re-create what custom_motif_draw.html actually POSTs.

    Its submit handler serialises `id: node.name`, so the wire format uses the
    display names ('N1', 'N2', ...) even though CUSTOM_MOTIF_EXAMPLES carries
    integer ids for the canvas. Tests must use this shape, not the raw example.
    """
    names = {node["id"]: node["name"] for node in example["nodes"]}
    return {
        "nodes": [
            {"id": node["name"], "x": node["x"], "y": node["y"]}
            for node in example["nodes"]
        ],
        "edges": [
            {
                "from": names[edge["from"]],
                "to": names[edge["to"]],
                "type": edge["type"],
                "attribute": edge["attribute"],
            }
            for edge in example["edges"]
        ],
    }


def test_posted_payload_keeps_node_names_as_ids():
    """Pins the wire format the drawing page produces."""
    graph = flask_app.create_networkx_graph(
        as_posted_payload(flask_app.CUSTOM_MOTIF_EXAMPLES[0])
    )
    assert sorted(graph.nodes) == ["N1", "N2", "N3", "N4", "N5"]


@pytest.mark.parametrize("node", [1, 7, 42])
def test_sort_key_tolerates_non_string_node_ids(node):
    """Hardening, not a live bug: the browser posts names, so sort_key only ever
    sees strings today. It still has to stay total, because CUSTOM_MOTIF_EXAMPLES
    carries integer ids and anything driving the endpoint from that data (a
    script, a future caller) would otherwise abort the search with AttributeError.
    """
    assert find_hairpin.sort_key(node) == ("", 0)


@pytest.mark.parametrize(
    "node,expected",
    [("A1", ("A", 1)), ("A-1", ("A", 1)), ("'0'1", ("0", 1)), ("'0'-1", ("0", 1))],
)
def test_sort_key_still_orders_library_node_ids(node, expected):
    """All four chain-id spellings must keep sorting the way they always have."""
    assert find_hairpin.sort_key(node) == expected


def test_search_engine_runs_against_a_browser_drawn_target():
    """End-to-end over the real engine with a real library graph.

    Guards the path a submitted drawing actually takes: posted JSON ->
    create_networkx_graph -> find_matching_subgraphs. Runs on the single
    structure the example is drawn from, so it stays fast.
    """
    example = next(
        item for item in flask_app.CUSTOM_MOTIF_EXAMPLES if item["key"] == "bulge-1x0"
    )
    target = flask_app.create_networkx_graph(as_posted_payload(example))

    with open(flask_app.os.path.join(flask_app.BASE_DIR, "batch_0000_graphs.pickle"), "rb") as handle:
        graphs = pickle.load(handle)
    graph_id = next(
        name for name in graphs if name[:4].upper() == example["source_pdb"]
    )

    conn = sqlite3.connect(":memory:")
    conn.execute(
        "CREATE TABLE files (id INTEGER PRIMARY KEY AUTOINCREMENT, motif_type TEXT,"
        " pdbid TEXT, paired_nt_number TEXT, nt_number TEXT, filecontent TEXT,"
        " created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP)"
    )
    try:
        find_hairpin.find_matching_subgraphs(
            {graph_id: graphs[graph_id]}, {"custom_motif_1": target}, conn
        )
        rows = conn.execute(
            "SELECT pdbid, nt_number, filecontent FROM files"
        ).fetchall()
    finally:
        conn.close()

    assert rows, "browser-drawn target found no match in its own source structure"
    assert all(pdbid == example["source_pdb"] for pdbid, _, _ in rows)
    assert all(len(nt.split(",")) == len(example["nodes"]) for _, nt, _ in rows)
    assert all("ATOM" in content for _, _, content in rows)


@pytest.mark.parametrize("example", flask_app.CUSTOM_MOTIF_EXAMPLES)
def test_custom_example_matches_its_source_structure(example):
    expected = EXPECTED_EXAMPLES[example["key"]]
    query = as_networkx_graph(example)
    assert query.number_of_nodes() == expected["nodes"]
    assert query.number_of_edges() == expected["edges"]

    with open(flask_app.os.path.join(flask_app.BASE_DIR, "batch_0000_graphs.pickle"), "rb") as handle:
        graphs = pickle.load(handle)
    source_graph = next(
        graph for name, graph in graphs.items()
        if name[:4].upper() == expected["source_pdb"]
    )
    matcher = isomorphism.GraphMatcher(
        source_graph,
        query,
        edge_match=lambda source, target: source["attribute"] == target["attribute"],
    )
    assert matcher.subgraph_is_isomorphic()
