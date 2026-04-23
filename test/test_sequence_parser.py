"""
Test PDBParser-based sequence extraction on real ATLAS.db data.
Run: pytest test/test_sequence_parser.py -v -s
"""
import sys, os, io, sqlite3
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
from Bio.PDB import PDBParser
from Bio import BiopythonWarning
import warnings

DB_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'ATLAS.db')

RESNAME_MAP = {
    'A': 'A', 'ADE': 'A',
    'G': 'G', 'GUA': 'G',
    'C': 'C', 'CYT': 'C',
    'U': 'U', 'URI': 'U', 'URA': 'U',
}

def extract_sequence_biopython(filecontent, nt_list):
    """Extract sequence using Biopython PDBParser."""
    parser = PDBParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        structure = parser.get_structure('motif', io.StringIO(filecontent))

    # Build lookup: chain_id+resseq -> resname
    # Also build with stripped quotes to handle '0'390 -> 0390 style keys
    residues = {}
    for model in structure:
        for chain in model:
            for res in chain:
                _, resseq, _ = res.get_id()
                key = f'{chain.id}{resseq}'
                residues[key] = res.resname.strip()

    seq = []
    for nt in nt_list:
        # Strip quotes from node IDs like '0'390 -> 0390
        nt_clean = nt.replace("'", "")
        resname = residues.get(nt_clean)
        if resname is None:
            seq.append('?')
        else:
            seq.append(RESNAME_MAP.get(resname, f'[{resname}]'))
    return ''.join(seq)


def get_samples(motif_type, nt_count=None, n=1000):
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT nt_number, filecontent FROM data WHERE motif_type=? LIMIT 5000", (motif_type,))
    results = []
    for row in c:
        nts = row[0].split(',')
        if nt_count is None or len(nts) == nt_count:
            results.append((row[0], row[1]))
        if len(results) >= n:
            break
    conn.close()
    return results


@pytest.mark.parametrize("motif_type,nt_count", [
    ("hairpin", 6),
    ("hairpin", 10),
    ("internal", 6),
    ("bulge", None),
])
def test_parse_success_rate(motif_type, nt_count):
    """Parse 1000 samples, report failure and unknown base rates."""
    samples = get_samples(motif_type, nt_count, n=1000)
    assert len(samples) > 0, f"No samples found for {motif_type} {nt_count}nt"

    failures = []
    unknown = []
    for nt_number, filecontent in samples:
        nt_list = nt_number.split(',')
        try:
            seq = extract_sequence_biopython(filecontent, nt_list)
            if '?' in seq or '[' in seq:
                unknown.append((nt_number, seq))
        except Exception as e:
            failures.append((nt_number, str(e)))

    total = len(samples)
    print(f"\n{motif_type} {nt_count}nt — {total} samples")
    print(f"  Parse failures:              {len(failures)} ({100*len(failures)/total:.1f}%)")
    print(f"  Unknown/non-standard bases:  {len(unknown)} ({100*len(unknown)/total:.1f}%)")
    for nt, seq in unknown[:3]:
        print(f"    {nt} -> {seq}")

    assert len(failures) == 0, f"Parse failures:\n" + "\n".join(f"  {n}: {e}" for n, e in failures[:5])


def test_sequence_values_are_valid_bases():
    """Sequences should only contain A/G/C/U (or non-standard markers)."""
    samples = get_samples("hairpin", 6, n=50)
    for nt_number, filecontent in samples:
        nt_list = nt_number.split(',')
        seq = extract_sequence_biopython(filecontent, nt_list)
        valid = set('AGCU?[]')
        invalid = [c for c in seq if c not in valid]
        assert not invalid, f"Unexpected chars {invalid} in seq '{seq}' for {nt_number}"


def test_sequence_length_matches_nt_count():
    """Sequence length should equal number of nucleotides."""
    samples = get_samples("hairpin", 6, n=20)
    for nt_number, filecontent in samples:
        nt_list = nt_number.split(',')
        seq = extract_sequence_biopython(filecontent, nt_list)
        assert len(seq) == len(nt_list), f"Length mismatch: {len(seq)} vs {len(nt_list)} for {nt_number}"


def test_negative_residue_numbers():
    """Check if negative residue numbers are handled."""
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("SELECT nt_number, filecontent FROM data WHERE nt_number LIKE '%--%' OR nt_number LIKE '%,-' LIMIT 10")
    rows = c.fetchall()
    conn.close()
    if not rows:
        pytest.skip("No negative residue number samples found")
    for nt_number, filecontent in rows:
        nt_list = nt_number.split(',')
        seq = extract_sequence_biopython(filecontent, nt_list)
        print(f"  {nt_number} -> {seq}")
        assert len(seq) == len(nt_list)


# ── Bug 2 fix: direct ATOM-line parser ───────────────────────────────────────

def extract_sequence_fixed(filecontent, nt_list):
    """
    Parse ATOM lines directly to avoid Biopython crash on overlong atom names.

    Root cause: atom names like H5'' (5 chars) overflow PDB col 12-15,
    shifting chain ID and resseq fields right by 1. Biopython reads
    line[22:26] as resseq and gets 'A215' instead of '2150', raising ValueError.

    This parser reads chain (col 21) and resseq (cols 22-26) directly,
    with a fallback shift of +1 when resseq is not a valid integer.
    """
    residues = {}
    for line in filecontent.splitlines():
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        if len(line) < 26:
            continue
        chain_id = line[21]
        resname = line[17:20].strip()
        resseq_raw = line[22:26].strip()
        try:
            resseq = int(resseq_raw)
        except ValueError:
            # Atom name overflow: shift +1
            chain_id = line[22]
            resname = line[18:21].strip()
            resseq_raw = line[23:27].strip()
            try:
                resseq = int(resseq_raw)
            except ValueError:
                continue
        key = f'{chain_id}{resseq}'
        if key not in residues:
            residues[key] = resname

    seq = []
    for nt in nt_list:
        nt_clean = nt.strip().replace("'", "")
        resname = residues.get(nt_clean)
        seq.append(RESNAME_MAP.get(resname, f'[{resname}]') if resname else '?')
    return ''.join(seq)


# ── Tests for Bug 2 ───────────────────────────────────────────────────────────

OVERFLOW_IDS = [292444, 293091, 293156, 293224, 293290]


@pytest.mark.parametrize("record_id", OVERFLOW_IDS)
def test_buggy_parser_fails_on_overflow_records(record_id):
    """Confirm the old Biopython parser returns ??? on these records."""
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute('SELECT filecontent, nt_number FROM data WHERE id = ?', (record_id,))
    row = c.fetchone()
    conn.close()
    filecontent, nt_number = row
    nt_list = nt_number.split(',')
    try:
        seq = extract_sequence_biopython(filecontent, nt_list)
        # If it didn't raise, it should have returned all '?'
        assert all(c == '?' for c in seq), f"id={record_id}: expected all '?', got {seq!r}"
    except Exception:
        pass  # Exception is also acceptable — confirms the bug


@pytest.mark.parametrize("record_id,expected", [
    (292444, 'CUG'),
    (293091, 'CUG'),
    (293156, 'CUG'),
    (293224, 'CUG'),
    (293290, 'CUG'),
])
def test_fixed_parser_resolves_overflow_records(record_id, expected):
    """Fixed parser must return real nucleotide sequence, no '?'."""
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute('SELECT filecontent, nt_number FROM data WHERE id = ?', (record_id,))
    row = c.fetchone()
    conn.close()
    filecontent, nt_number = row
    nt_list = nt_number.split(',')
    result = extract_sequence_fixed(filecontent, nt_list)
    assert '?' not in result, f"id={record_id}: still got '?' in {result!r}"
    assert result == expected, f"id={record_id}: expected {expected!r}, got {result!r}"


def test_fixed_parser_matches_biopython_on_normal_records():
    """Fixed parser must agree with Biopython on records without overflow atoms."""
    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()
    c.execute("""
        SELECT id, filecontent, nt_number FROM data
        WHERE motif_type = 'hairpin'
          AND filecontent NOT LIKE '%H5''''%'
          AND filecontent NOT LIKE '%HO2''%'
        LIMIT 200
    """)
    rows = c.fetchall()
    conn.close()
    mismatches = []
    for record_id, filecontent, nt_number in rows:
        nt_list = nt_number.split(',')
        try:
            buggy = extract_sequence_biopython(filecontent, nt_list)
        except Exception:
            continue  # skip records biopython can't handle
        fixed = extract_sequence_fixed(filecontent, nt_list)
        if buggy != fixed:
            mismatches.append((record_id, buggy, fixed))
    assert not mismatches, (
        f"{len(mismatches)} mismatches:\n" +
        "\n".join(f"  id={r}: biopython={b!r} fixed={f!r}" for r, b, f in mismatches[:5])
    )
