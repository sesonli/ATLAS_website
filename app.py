# -*- coding: utf-8 -*-
"""
Created on Nov 25 11:40:07 2024
@author: JingyiLi
"""

from flask import Flask, render_template, request, send_file, jsonify, Response, abort
import sqlite3
import csv
import warnings
from Bio.PDB import PDBParser
from Bio import BiopythonWarning
import pandas as pd
import os
import shutil
import fcntl
import io
import zipfile
from io import BytesIO
from flask import send_from_directory
import json
import pickle
import networkx as nx
import subprocess
import time
from threading import Lock
import uuid
import tempfile
import logging
from logging.handlers import RotatingFileHandler
from datetime import datetime, timedelta
from collections import Counter

# Configure logging with rotation (10MB max, 3 backups)
log_handler = RotatingFileHandler(
    'app.log',
    maxBytes=10*1024*1024,  # 10MB
    backupCount=3
)
log_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
log_handler.setLevel(logging.INFO)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(log_handler)

app = Flask(__name__, static_folder='static', static_url_path='/static')

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
KTURN_WORKFLOW_DATA_DIR = os.path.join(BASE_DIR, 'static', 'data', 'kturn_workflow')
NAMED_MOTIF_WORKFLOW_DIRS = {
    'sarcin-ricin': os.path.join(
        BASE_DIR, 'static', 'data', 'named_motifs', 'sarcin_ricin'
    ),
    'gnra': os.path.join(
        BASE_DIR, 'static', 'data', 'named_motifs', 'gnra'
    ),
}
ATLAS_DB_PATH = os.environ.get('ATLAS_DB_PATH', os.path.join(BASE_DIR, 'ATLAS.db'))

# Database-backed examples for the custom graph drawing interface. Each graph
# is an attributed topology taken from a motif record whose source structure is
# present in the 100-structure custom-search sample (batch_0000_graphs.pickle).
CUSTOM_MOTIF_EXAMPLES = [
    {
        'key': 'bulge-1x0',
        'title': 'Bulge: 1 × 0 loop',
        'source_pdb': '1AQ3',
        'nodes': [
            {'id': 1, 'x': 250, 'y': 160, 'name': 'N1'},
            {'id': 2, 'x': 400, 'y': 160, 'name': 'N2'},
            {'id': 3, 'x': 550, 'y': 160, 'name': 'N3'},
            {'id': 4, 'x': 550, 'y': 340, 'name': 'N4'},
            {'id': 5, 'x': 250, 'y': 340, 'name': 'N5'},
        ],
        'edges': [
            {'from': 1, 'to': 2, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 2, 'to': 3, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 4, 'to': 5, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 1, 'to': 5, 'type': 'wc', 'attribute': [1, 0, 0]},
            {'from': 3, 'to': 4, 'type': 'wc', 'attribute': [1, 0, 0]},
        ],
    },
    {
        'key': 'internal-loop-1x1',
        'title': 'Internal loop: 1 × 1',
        'source_pdb': '17RA',
        'nodes': [
            {'id': 1, 'x': 250, 'y': 150, 'name': 'N1'},
            {'id': 2, 'x': 400, 'y': 150, 'name': 'N2'},
            {'id': 3, 'x': 550, 'y': 150, 'name': 'N3'},
            {'id': 4, 'x': 550, 'y': 350, 'name': 'N4'},
            {'id': 5, 'x': 400, 'y': 350, 'name': 'N5'},
            {'id': 6, 'x': 250, 'y': 350, 'name': 'N6'},
        ],
        'edges': [
            {'from': 1, 'to': 2, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 2, 'to': 3, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 4, 'to': 5, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 5, 'to': 6, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 1, 'to': 6, 'type': 'wc', 'attribute': [1, 0, 0]},
            {'from': 3, 'to': 4, 'type': 'wc', 'attribute': [1, 0, 0]},
            {'from': 2, 'to': 5, 'type': 'non-wc', 'attribute': [0, 1, 0]},
        ],
    },
    {
        'key': 'internal-loop-2x2',
        'title': 'Internal loop: 2 × 2',
        'source_pdb': '1A4D',
        'nodes': [
            {'id': 1, 'x': 200, 'y': 155, 'name': 'N1'},
            {'id': 2, 'x': 335, 'y': 155, 'name': 'N2'},
            {'id': 3, 'x': 465, 'y': 155, 'name': 'N3'},
            {'id': 4, 'x': 600, 'y': 155, 'name': 'N4'},
            {'id': 5, 'x': 200, 'y': 345, 'name': 'N5'},
            {'id': 6, 'x': 335, 'y': 345, 'name': 'N6'},
            {'id': 7, 'x': 465, 'y': 345, 'name': 'N7'},
            {'id': 8, 'x': 600, 'y': 345, 'name': 'N8'},
        ],
        'edges': [
            {'from': 1, 'to': 2, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 2, 'to': 3, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 3, 'to': 4, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 5, 'to': 6, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 6, 'to': 7, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 7, 'to': 8, 'type': 'covalent', 'attribute': [0, 0, 1]},
            {'from': 1, 'to': 5, 'type': 'wc', 'attribute': [1, 0, 0]},
            {'from': 4, 'to': 8, 'type': 'wc', 'attribute': [1, 0, 0]},
            {'from': 2, 'to': 6, 'type': 'non-wc', 'attribute': [0, 1, 0]},
            {'from': 3, 'to': 7, 'type': 'non-wc', 'attribute': [0, 1, 0]},
        ],
    },
]

# Create lock for custom search to prevent concurrent executions (avoid memory issues on 4GB systems)
custom_search_lock = Lock()

def cleanup_on_startup():
    """Clean up temporary files and old database records on application startup"""
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. Clean temporary PDB folders
    temp_dirs = ['pdb_files', 'custom_pdb_files', 'search_sessions']
    for dir_name in temp_dirs:
        dir_path = os.path.join(base_dir, dir_name)
        if os.path.exists(dir_path):
            try:
                shutil.rmtree(dir_path)
                logger.info(f"Cleaned up {dir_name}/ directory")
            except Exception as e:
                logger.warning(f"Failed to delete {dir_name}/: {e}")
        os.makedirs(dir_path, exist_ok=True)

    # 2. Clean temporary CSV files
    temp_files = ['search_results.csv', 'custom_search_results.csv', 'temp_target_graphs.pickle']
    for file_name in temp_files:
        file_path = os.path.join(base_dir, file_name)
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
                logger.info(f"Deleted temporary file: {file_name}")
            except Exception as e:
                logger.warning(f"Failed to delete {file_name}: {e}")

    # 3. Truncate flask_output.log if it exceeds 10 MB
    flask_log = os.path.join(base_dir, 'flask_output.log')
    if os.path.exists(flask_log) and os.path.getsize(flask_log) > 10 * 1024 * 1024:
        try:
            with open(flask_log, 'w'):
                pass
            logger.info("Truncated flask_output.log (exceeded 10 MB)")
        except Exception as e:
            logger.warning(f"Failed to truncate flask_output.log: {e}")

    # 4. Clean hairpin.db records older than 7 days
    hairpin_db_path = os.path.join(base_dir, 'hairpin.db')
    if os.path.exists(hairpin_db_path):
        conn = None
        try:
            conn = sqlite3.connect(hairpin_db_path)
            cursor = conn.cursor()

            # Calculate date 7 days ago
            cutoff_date = (datetime.now() - timedelta(days=7)).strftime('%Y-%m-%d %H:%M:%S')

            # Check if created_at column exists
            cursor.execute("PRAGMA table_info(files)")
            columns = [col[1] for col in cursor.fetchall()]

            if 'created_at' in columns:
                # Delete old records
                cursor.execute("DELETE FROM files WHERE created_at < ?", (cutoff_date,))
                deleted_count = cursor.rowcount
                conn.commit()

                # Vacuum to reclaim space
                cursor.execute("VACUUM")
                logger.info(f"Cleaned up hairpin.db: deleted {deleted_count} records older than 7 days")
            else:
                logger.info("hairpin.db: created_at column not found, skipping cleanup")
        except Exception as e:
            logger.warning(f"Failed to clean hairpin.db: {e}")
        finally:
            if conn:
                conn.close()

    logger.info("Startup cleanup completed")

# Run cleanup on startup unless an isolated staging/test process requests a
# side-effect-free import.
if os.environ.get('ATLAS_SKIP_STARTUP_CLEANUP') != '1':
    cleanup_on_startup()

# Connect to SQLite database
def get_db_connection():
    conn = sqlite3.connect(ATLAS_DB_PATH)
    return conn

# Extract nucleotide sequence from PDB filecontent using Biopython PDBParser
_RESNAME_MAP = {
    'A': 'A', 'ADE': 'A',
    'G': 'G', 'GUA': 'G',
    'C': 'C', 'CYT': 'C',
    'U': 'U', 'URI': 'U', 'URA': 'U',
}

def extract_sequence(filecontent, nt_list):
    # Parse ATOM lines directly to avoid Biopython crash on overlong atom names.
    # Atom names like H5'' (5 chars) overflow PDB col 12-15, shifting chain ID
    # and resseq fields right by 1. Biopython reads line[22:26] as resseq and
    # gets 'A215' instead of '2150', raising ValueError -> returning '???'.
    try:
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
                # Atom name overflow: shift all fields +1
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
            seq.append(_RESNAME_MAP.get(resname, f'[{resname}]') if resname else '?')
        return ''.join(seq)
    except Exception:
        return '?' * len(nt_list)

# Check if residue data is complete (has sufficient atoms for each nucleotide)
def is_residue_data_complete(filecontent, nt_list, min_atoms_per_residue=10):
    """
    Verify that each nucleotide in the structure has complete atom data.

    Args:
        filecontent: PDB file content string
        nt_list: List of nucleotide identifiers (e.g., ['A726', 'A727', 'A728'])
        min_atoms_per_residue: Minimum number of atoms required (default 10)

    Returns:
        True if all residues have sufficient atoms, False otherwise
    """
    # Tally every residue in one pass. Scanning the file once per nucleotide
    # instead costs O(residues x lines), which dominates the search on the long
    # pseudoknots: a 252-nucleotide LR motif re-reads ~5,600 ATOM lines 252 times.
    atom_counts = Counter()
    for line in filecontent.strip().split('\n'):
        if line.startswith('ATOM') and len(line) >= 27:
            atom_counts[f'{line[21]}{line[22:27].strip()}'] += 1

    # Standard RNA nucleotides should have 15-23 atoms
    # Use conservative threshold of 10 to avoid false positives
    return all(atom_counts[nt] >= min_atoms_per_residue for nt in nt_list)

# Check C4' atom distances for coordinate quality
def check_c4_distances(filecontent, nt_list, min_dist=4.0, max_dist=8.0):
    """
    Check if C4' atom distances between adjacent nucleotides are within normal range.
    Only checks distances within continuous segments (residue numbers differ by ≤2).

    Args:
        filecontent: PDB file content string
        nt_list: List of nucleotide identifiers (e.g., ['A10', 'A11', 'B5'])
        min_dist: Minimum acceptable distance in Angstroms (default 4.0)
        max_dist: Maximum acceptable distance in Angstroms (default 8.0)

    Returns:
        tuple: (is_normal: bool, reason: str)
            - (True, "Normal") if all distances are within range
            - (False, "Poor geometry") if any distance < min_dist or > max_dist
    """
    import numpy as np
    import re

    try:
        lines = filecontent.strip().split('\n')
        c4_atoms = {}

        # Extract C4' atoms and store by nucleotide ID
        for line in lines:
            if not line.startswith('ATOM'):
                continue

            try:
                atom_name = line[12:16].strip()
                if atom_name != "C4'":
                    continue

                chain = line[21]
                res_num_str = line[22:27].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                full_id = f"{chain}{res_num_str}"
                c4_atoms[full_id] = {
                    'chain': chain,
                    'res_num': int(res_num_str),
                    'coords': np.array([x, y, z])
                }
            except (ValueError, IndexError):
                # Skip malformed lines
                continue

        # Parse nt_list and group by chain
        parsed_nts = []
        for nt in nt_list:
            # Match patterns like 'A10', "'0'10", 'A-5'
            match = re.match(r"('[^']+')?([\w]+?)(-?\d+)", nt.strip())
            if match:
                chain_prefix = match.group(1) if match.group(1) else ''
                chain = chain_prefix + match.group(2)
                res_num = int(match.group(3))
                parsed_nts.append({
                    'original': nt.strip(),
                    'chain': chain,
                    'res_num': res_num
                })

        if len(parsed_nts) < 2:
            return (True, "Normal")

        # Group by chain
        chain_groups = {}
        for pnt in parsed_nts:
            chain = pnt['chain']
            if chain not in chain_groups:
                chain_groups[chain] = []
            chain_groups[chain].append(pnt)

        # Process each chain separately
        for chain, nts in chain_groups.items():
            # Sort by residue number
            nts_sorted = sorted(nts, key=lambda x: x['res_num'])

            # Identify continuous segments (gap ≤ 2)
            segments = []
            current_segment = [nts_sorted[0]]

            for i in range(1, len(nts_sorted)):
                gap = nts_sorted[i]['res_num'] - nts_sorted[i-1]['res_num']
                if gap <= 2:
                    # Continuous or small insertion
                    current_segment.append(nts_sorted[i])
                else:
                    # Large gap, start new segment
                    segments.append(current_segment)
                    current_segment = [nts_sorted[i]]
            segments.append(current_segment)

            # Check C4' distances only within continuous segments
            for segment in segments:
                if len(segment) < 2:
                    continue

                for i in range(len(segment) - 1):
                    nt1 = segment[i]['original']
                    nt2 = segment[i+1]['original']

                    # Both C4' atoms must exist
                    if nt1 not in c4_atoms or nt2 not in c4_atoms:
                        continue

                    coords1 = c4_atoms[nt1]['coords']
                    coords2 = c4_atoms[nt2]['coords']
                    dist = np.linalg.norm(coords1 - coords2)

                    if dist < min_dist or dist > max_dist:
                        return (False, "Poor geometry")

        return (True, "Normal")

    except Exception as e:
        # On any error, assume normal to avoid false positives
        return (True, "Normal")

def has_usable_coordinates(filecontent):
    """Whether a row carries coordinates the app can serve.

    This is the same predicate the "All Results" CSV/ZIP downloads apply, so
    search counts and download counts stay in agreement.
    """
    return bool(filecontent) and 'No matching nodes found' not in filecontent


def count_remaining_matches(cursor, nt_number=None):
    """Finish counting matches on a cursor that result collection stopped early.

    Collection breaks at max_results, so anything past that point was never
    counted. Draining the rest here keeps the truncation notice honest; it costs
    well under a second because no coordinates are parsed.

    Returns:
        tuple: (matches_found, rows_scanned)
    """
    matched = 0
    scanned = 0
    for row in cursor:
        scanned += 1
        if not row[3]:
            continue
        if nt_number is not None and len(row[3].split(',')) != nt_number:
            continue
        if has_usable_coordinates(row[4]):
            matched += 1
    return matched, scanned


# Search RNA motifs based on motif_type and nt_number or subtype
def search_motif(motif_type, nt_number=None, max_results=10000):
    """
    Search for RNA motifs with memory-efficient cursor iteration.

    Args:
        motif_type: Type of motif to search for
        nt_number: Number of nucleotides (for hairpin/internal/bulge)
        max_results: Maximum number of rows returned for display (the full set
            stays available through the "All Results" downloads)

    Returns:
        tuple: (filtered_rows, metadata)
        metadata includes 'truncated' flag and 'total_scanned'. 'total_matched'
        counts every matching row in the database, including those past the cap.
    """
    conn = None
    try:
        conn = get_db_connection()
        cursor = conn.cursor()

        filtered_rows = []
        was_truncated = False
        total_scanned = 0
        total_matched = 0

        if motif_type in ["hairpin", "internal", "bulge"]:
            if nt_number is None:
                conn.close()
                return [], {'error': 'nt_number is required for this motif type', 'total_scanned': 0, 'total_matched': 0, 'result_count': 0, 'truncated': False}
            # Query for hairpin, internal, and bulge from "data" table
            # Use cursor iteration instead of fetchall() to avoid loading all rows into memory
            query = "SELECT id, motif_type, pdbid, nt_number, filecontent FROM data WHERE motif_type = ?"
            cursor.execute(query, (motif_type,))

            for row in cursor:
                total_scanned += 1
                if not row[3]:
                    continue
                if len(row[3].split(',')) == int(nt_number):
                    if has_usable_coordinates(row[4]):
                        total_matched += 1
                        is_atom_complete = is_residue_data_complete(row[4], row[3].split(','))
                        is_coord_normal, coord_reason = check_c4_distances(row[4], row[3].split(','))
                        is_quality_ok = is_atom_complete and is_coord_normal
                        warning_reason = ""
                        if not is_atom_complete:
                            warning_reason = "Incomplete atom data"
                        if not is_coord_normal:
                            warning_reason += (" + " if warning_reason else "") + coord_reason
                        if is_quality_ok:
                            warning_reason = ""
                        sequence = extract_sequence(row[4], row[3].split(','))
                        row_with_flag = list(row) + [sequence, is_quality_ok, warning_reason]
                        filtered_rows.append(tuple(row_with_flag))
                        if len(filtered_rows) >= max_results:
                            was_truncated = True
                            break

            if was_truncated:
                extra_matched, extra_scanned = count_remaining_matches(cursor, int(nt_number))
                total_matched += extra_matched
                total_scanned += extra_scanned

        elif 'junction' in motif_type:
            # Query for multiway junctions
            motif_type_db = motif_type.replace("-way", "").replace(" ", "_")
            query = "SELECT id, motif_type, pdbid, nt_number, filecontent FROM data WHERE motif_type = ?"
            cursor.execute(query, (motif_type_db,))

            # Use cursor iteration with limit
            for row in cursor:
                total_scanned += 1
                if not row[3]:
                    continue
                if has_usable_coordinates(row[4]):
                    total_matched += 1
                    is_atom_complete = is_residue_data_complete(row[4], row[3].split(','))
                    is_coord_normal, coord_reason = check_c4_distances(row[4], row[3].split(','))
                    is_quality_ok = is_atom_complete and is_coord_normal
                    warning_reason = ""
                    if not is_atom_complete:
                        warning_reason = "Incomplete atom data"
                    if not is_coord_normal:
                        warning_reason += (" + " if warning_reason else "") + coord_reason
                    if is_quality_ok:
                        warning_reason = ""
                    sequence = extract_sequence(row[4], row[3].split(','))
                    row_with_flag = list(row) + [sequence, is_quality_ok, warning_reason]
                    filtered_rows.append(tuple(row_with_flag))
                    if len(filtered_rows) >= max_results:
                        was_truncated = True
                        break

            if was_truncated:
                extra_matched, extra_scanned = count_remaining_matches(cursor)
                total_matched += extra_matched
                total_scanned += extra_scanned

        elif 'pseudoknot' in motif_type:
            # Query for pseudoknot subtypes from "PK" table
            # Note: PK table's id column is NULL, so we use rowid instead
            motif_type_db = motif_type
            query = "SELECT rowid, motif_type, pdbid, nt_number, file_content FROM PK WHERE motif_type = ?"
            cursor.execute(query, (motif_type_db.strip(),))

            # Use cursor iteration with limit
            for row in cursor:
                total_scanned += 1
                if not row[3]:
                    continue
                if has_usable_coordinates(row[4]):
                    total_matched += 1
                    is_atom_complete = is_residue_data_complete(row[4], row[3].split(','))
                    is_coord_normal, coord_reason = check_c4_distances(row[4], row[3].split(','))
                    is_quality_ok = is_atom_complete and is_coord_normal
                    warning_reason = ""
                    if not is_atom_complete:
                        warning_reason = "Incomplete atom data"
                    if not is_coord_normal:
                        warning_reason += (" + " if warning_reason else "") + coord_reason
                    if is_quality_ok:
                        warning_reason = ""
                    sequence = extract_sequence(row[4], row[3].split(','))
                    row_with_flag = list(row) + [sequence, is_quality_ok, warning_reason]
                    filtered_rows.append(tuple(row_with_flag))
                    if len(filtered_rows) >= max_results:
                        was_truncated = True
                        break

            if was_truncated:
                extra_matched, extra_scanned = count_remaining_matches(cursor)
                total_matched += extra_matched
                total_scanned += extra_scanned

        else:
            # If motif_type is unknown
            filtered_rows = []

        conn.close()

        # If results are found, return with metadata
        if filtered_rows:
            base_dir = os.path.dirname(os.path.abspath(__file__))

            # Use UUID session dir to avoid concurrent search conflicts
            session_id = str(uuid.uuid4())
            sessions_dir = os.path.join(base_dir, 'search_sessions')
            os.makedirs(sessions_dir, exist_ok=True)
            session_dir = os.path.join(sessions_dir, session_id)
            os.makedirs(session_dir)

            # Save PDB files for quick download
            for row in filtered_rows:
                with open(os.path.join(session_dir, f'{row[0]}.pdb'), 'w') as f:
                    f.write(row[4])

            # Save current-page CSV for quick download
            csv_path = os.path.join(session_dir, 'search_results_current.csv')
            df = pd.DataFrame(filtered_rows, columns=['ID', 'Motif Type', 'PDB ID', 'NT Number', 'Filecontent', 'Sequence', 'Quality_OK', 'Warning_Reason'])
            df[['ID', 'Motif Type', 'PDB ID', 'NT Number', 'Sequence', 'Quality_OK', 'Warning_Reason']].to_csv(csv_path, index=False, quoting=1)  # QUOTE_ALL prevents Excel from misreading pdbids like 8E35 as scientific notation

            metadata = {
                'truncated': was_truncated,
                'total_scanned': total_scanned,
                'total_matched': total_matched,
                'result_count': len(filtered_rows),
                'max_results': max_results,
                'session_id': session_id
            }
            return filtered_rows, metadata

        # No results found
        return [], {'truncated': False, 'total_scanned': total_scanned, 'total_matched': total_matched, 'result_count': 0}

    except sqlite3.Error as e:
        print(f"Error during search: {e}")
        return [], {'error': str(e)}
    finally:
        if conn:
            conn.close()

# Route to the home page
@app.route('/')
def home():
    return render_template('index.html')

# Browsers request /favicon.ico on every page; without this each view logs a 404
@app.route('/favicon.ico')
def favicon():
    return send_from_directory(
        os.path.join(BASE_DIR, 'static', 'images'),
        'favicon.ico',
        mimetype='image/vnd.microsoft.icon',
    )

# User Guide page
@app.route('/user-guide', methods=['GET'])
def user_guide():
    return render_template('user_guide.html')

# Route to handle motif search requests
@app.route('/search', methods=['POST'])
def search():
    motif_type = request.form.get('motif_type')
    nt_number = request.form.get('nt_number')
    if not motif_type:
        return "Missing motif_type parameter.", 400
    if motif_type == 'multiway-junction':
        # Here nt_number will be a string like '3-way junction'
        motif_type = nt_number
        if not motif_type:
            return "Missing junction type parameter.", 400
        nt_number_count = None
    elif motif_type == 'pseudoknot':
        if not nt_number:
            return "Missing pseudoknot subtype parameter.", 400
        motif_type = f'pseudoknot_{nt_number}'
        nt_number_count = None
    else:
        # For other types, it's a number, so convert to integer
        try:
            nt_number_count = int(nt_number)
        except (ValueError, TypeError):
            nt_number_count = None

    search_results, metadata = search_motif(motif_type, nt_number_count)

    if metadata.get('error'):
        return metadata['error'], 400

    if search_results:
        warning_message = None
        if metadata.get('truncated'):
            warning_message = f"Showing {metadata['result_count']:,} of {metadata['total_matched']:,} matching results. Download CSV/ZIP (All Results) for complete data."

        return render_template('results.html',
                               search_results=search_results,
                               motif_type=motif_type,
                               nt_number=nt_number_count or '',
                               session_id=metadata.get('session_id', ''),
                               warning=warning_message,
                               metadata=metadata)
    else:
        return render_template('no_results.html')

# Route to download CSV of current results (quick, from pre-saved file)
@app.route('/download_csv', methods=['GET'])
def download_csv():
    session_id = request.args.get('session_id', '')
    base_dir = os.path.dirname(os.path.abspath(__file__))
    sessions_dir = os.path.realpath(os.path.join(base_dir, 'search_sessions'))
    csv_path = os.path.realpath(os.path.join(sessions_dir, session_id, 'search_results_current.csv'))
    if not session_id or not csv_path.startswith(sessions_dir + os.sep) or not os.path.exists(csv_path):
        return "CSV file not found. Please run a search first.", 404
    return send_file(csv_path, mimetype='text/csv', as_attachment=True, download_name='search_results.csv')

# Route to download CSV of all results (full query, slow)
@app.route('/download_csv_full', methods=['GET'])
def download_csv_full():
    motif_type = request.args.get('motif_type')
    nt_number = request.args.get('nt_number')

    if not motif_type:
        return "Missing motif_type parameter.", 400

    try:
        nt_number_int = int(nt_number) if nt_number and nt_number != 'None' else None
    except ValueError:
        nt_number_int = None

    def generate_csv():
        conn = get_db_connection()
        try:
            cursor = conn.cursor()
            output = io.StringIO()
            writer = csv.writer(output, quoting=csv.QUOTE_ALL)
            writer.writerow(['ID', 'Motif Type', 'PDB ID', 'NT Number', 'Sequence', 'Quality_OK', 'Warning_Reason'])
            yield output.getvalue()
            output.truncate(0)
            output.seek(0)

            if motif_type in ["hairpin", "internal", "bulge"]:
                cursor.execute("SELECT id, motif_type, pdbid, nt_number, filecontent FROM data WHERE motif_type = ?", (motif_type,))
                for row in cursor:
                    if not row[3]:
                        continue
                    if nt_number_int and len(row[3].split(',')) != nt_number_int:
                        continue
                    if not row[4] or 'No matching nodes found' in row[4]:
                        continue
                    is_atom_complete = is_residue_data_complete(row[4], row[3].split(','))
                    is_coord_normal, coord_reason = check_c4_distances(row[4], row[3].split(','))
                    is_quality_ok = is_atom_complete and is_coord_normal
                    warning_reason = ""
                    if not is_atom_complete:
                        warning_reason = "Incomplete atom data"
                    if not is_coord_normal:
                        warning_reason += (" + " if warning_reason else "") + coord_reason
                    writer.writerow([row[0], row[1], row[2], row[3], extract_sequence(row[4], row[3].split(',')), is_quality_ok, warning_reason])
                    yield output.getvalue()
                    output.truncate(0)
                    output.seek(0)
            elif 'junction' in motif_type:
                motif_type_db = motif_type.replace("-way", "").replace(" ", "_")
                cursor.execute("SELECT id, motif_type, pdbid, nt_number, filecontent FROM data WHERE motif_type = ?", (motif_type_db,))
                for row in cursor:
                    if not row[3]:
                        continue
                    if not row[4] or 'No matching nodes found' in row[4]:
                        continue
                    is_atom_complete = is_residue_data_complete(row[4], row[3].split(','))
                    is_coord_normal, coord_reason = check_c4_distances(row[4], row[3].split(','))
                    is_quality_ok = is_atom_complete and is_coord_normal
                    warning_reason = ""
                    if not is_atom_complete:
                        warning_reason = "Incomplete atom data"
                    if not is_coord_normal:
                        warning_reason += (" + " if warning_reason else "") + coord_reason
                    writer.writerow([row[0], row[1], row[2], row[3], extract_sequence(row[4], row[3].split(',')), is_quality_ok, warning_reason])
                    yield output.getvalue()
                    output.truncate(0)
                    output.seek(0)
            elif 'pseudoknot' in motif_type:
                cursor.execute("SELECT rowid, motif_type, pdbid, nt_number, file_content FROM PK WHERE motif_type = ?", (motif_type.strip(),))
                for row in cursor:
                    if not row[3]:
                        continue
                    if not row[4] or 'No matching nodes found' in row[4]:
                        continue
                    is_atom_complete = is_residue_data_complete(row[4], row[3].split(','))
                    is_coord_normal, coord_reason = check_c4_distances(row[4], row[3].split(','))
                    is_quality_ok = is_atom_complete and is_coord_normal
                    warning_reason = ""
                    if not is_atom_complete:
                        warning_reason = "Incomplete atom data"
                    if not is_coord_normal:
                        warning_reason += (" + " if warning_reason else "") + coord_reason
                    writer.writerow([row[0], row[1], row[2], row[3], extract_sequence(row[4], row[3].split(',')), is_quality_ok, warning_reason])
                    yield output.getvalue()
                    output.truncate(0)
                    output.seek(0)
        finally:
            conn.close()

    return Response(
        generate_csv(),
        mimetype='text/csv',
        headers={'Content-Disposition': 'attachment; filename=search_results.csv'}
    )

# Route to download ZIP of current results (quick, from pre-saved files)
@app.route('/download_zip', methods=['GET'])
def download_zip():
    session_id = request.args.get('session_id', '')
    base_dir = os.path.dirname(os.path.abspath(__file__))
    sessions_dir = os.path.realpath(os.path.join(base_dir, 'search_sessions'))
    pdb_dir = os.path.realpath(os.path.join(sessions_dir, session_id))

    if not session_id or not pdb_dir.startswith(sessions_dir + os.sep) or not os.path.exists(pdb_dir) or not any(f.endswith('.pdb') for f in os.listdir(pdb_dir)):
        return "PDB files not found. Please run a search first.", 404

    temp_zip = tempfile.NamedTemporaryFile(delete=False, suffix='.zip')
    try:
        with zipfile.ZipFile(temp_zip.name, 'w', zipfile.ZIP_DEFLATED) as zf:
            for file in os.listdir(pdb_dir):
                if file.endswith('.pdb'):
                    zf.write(os.path.join(pdb_dir, file), arcname=file)
        response = send_file(temp_zip.name, mimetype='application/zip', as_attachment=True, download_name='pdb_files.zip')
        return response
    except Exception as e:
        return f"An error occurred during ZIP creation: {str(e)}", 500
    finally:
        if os.path.exists(temp_zip.name):
            os.unlink(temp_zip.name)

# Route to download ZIP of all results (full query, slow)
@app.route('/download_zip_full', methods=['GET'])
def download_zip_full():
    motif_type = request.args.get('motif_type')
    nt_number = request.args.get('nt_number')

    if not motif_type:
        return "Missing motif_type parameter.", 400

    try:
        nt_number_int = int(nt_number) if nt_number and nt_number != 'None' else None
    except ValueError:
        nt_number_int = None

    conn = None
    try:
        conn = get_db_connection()
        cursor = conn.cursor()

        if motif_type in ["hairpin", "internal", "bulge"]:
            cursor.execute("SELECT id, nt_number, filecontent FROM data WHERE motif_type = ?", (motif_type,))
        elif 'junction' in motif_type:
            motif_type_db = motif_type.replace("-way", "").replace(" ", "_")
            cursor.execute("SELECT id, nt_number, filecontent FROM data WHERE motif_type = ?", (motif_type_db,))
        elif 'pseudoknot' in motif_type:
            cursor.execute("SELECT rowid, nt_number, file_content FROM PK WHERE motif_type = ?", (motif_type.strip(),))
        else:
            return "Unknown motif type.", 400

        temp_zip = tempfile.NamedTemporaryFile(delete=False, suffix='.zip')
        try:
            with zipfile.ZipFile(temp_zip.name, 'w', zipfile.ZIP_DEFLATED) as zf:
                for row in cursor:
                    motif_id, nt_number_val, filecontent = row[0], row[1], row[2]
                    if not nt_number_val:
                        continue
                    if nt_number_int and len(nt_number_val.split(',')) != nt_number_int:
                        continue
                    if not filecontent or 'No matching nodes found' in filecontent:
                        continue
                    zf.writestr(f'{motif_id}.pdb', filecontent)
            response = send_file(temp_zip.name, mimetype='application/zip', as_attachment=True, download_name='pdb_files_all.zip')
            return response
        finally:
            if os.path.exists(temp_zip.name):
                os.unlink(temp_zip.name)
    except Exception as e:
        return f"An error occurred during ZIP creation: {str(e)}", 500
    finally:
        if conn:
            conn.close()


# Route to download the complete database
@app.route('/download-database', methods=['GET'])
def download_database():
    """Download the entire RNAdex database file."""
    try:
        if not os.path.exists(ATLAS_DB_PATH):
            return "Database file not found.", 404
        
        return send_file(
            ATLAS_DB_PATH,
            as_attachment=True,
            download_name='RNAdex.db',
            mimetype='application/x-sqlite3'
        )
    except Exception as e:
        return f"An error occurred during database download: {str(e)}", 500


# Route to custom motif search drawing page
@app.route('/custom-motif-search', methods=['GET'])
def custom_motif_search():
    """Render the custom motif drawing interface"""
    return render_template(
        'custom_motif_draw.html',
        custom_motif_examples=CUSTOM_MOTIF_EXAMPLES,
    )

def load_kturn_workflow_data():
    """Load the precomputed kink-turn worked example results."""
    summary_path = os.path.join(KTURN_WORKFLOW_DATA_DIR, 'caseC_kturn_summary.json')
    confirmed_path = os.path.join(KTURN_WORKFLOW_DATA_DIR, 'caseC_geometry_confirmed.json')

    with open(summary_path, 'r', encoding='utf-8') as f:
        summary = json.load(f)
    with open(confirmed_path, 'r', encoding='utf-8') as f:
        confirmed = json.load(f)

    representative_notes = {
        item.get('pdb'): item
        for item in summary.get('representatives', [])
    }
    confirmed_representatives = []
    for detail in confirmed.get('details', []):
        if not detail.get('confirmed_kink_turn'):
            continue
        note = representative_notes.get(detail.get('pdb'), {})
        confirmed_representatives.append({
            'pdb': detail.get('pdb'),
            'description': note.get('description', 'Confirmed kink-turn representative'),
            'sheared_pairs': detail.get('sheared_GA_pairs_SaengerXI', []),
            'n_sheared_ga': detail.get('n_sheared_GA'),
            'loop_strand_sizes': detail.get('loop_strand_sizes', []),
            'loop': note.get('loop', ''),
            'structure': f"motifs/{detail.get('pdb')}_kink_turn.pdb"
        })

    funnel = summary.get('funnel', {})
    stage_counts = [
        {
            'label': 'Search pool',
            'count': funnel.get('search_pool_internal_bulge_loops'),
            'description': 'All RNAdex internal and bulge loops considered for this example.'
        },
        {
            'label': 'Stage 1 graph-topology candidates',
            'count': funnel.get('stage1_graph_tandem_nonWC_topology'),
            'description': 'Loops carrying the tandem non-WC cross topology. These are candidates, not confirmed kink-turns.'
        },
        {
            'label': 'Candidates with G·A pairs',
            'count': funnel.get('with_GA_base_composition'),
            'description': 'Stage 1 candidates in which both non-WC pairs are G·A. Still candidates, not confirmed kink-turns.'
        },
        {
            'label': 'Representatives shown (of many candidates)',
            'count': confirmed.get('representatives_confirmed'),
            'description': (
                'Illustrative kink-turns confirmed by geometry on selected '
                'representatives; Stage 2 was not run exhaustively over all '
                f"{funnel.get('with_GA_base_composition', 0):,} composition "
                'candidates.'
            )
        }
    ]

    downloads = [
        {
            'label': 'Confirmed kink-turn structures (ZIP)',
            'filename': 'confirmed_kink_turn_structures.zip',
            'description': 'PDB coordinate files of the four geometry-confirmed kink-turn motifs.'
        },
        {
            'label': 'Stage 1 summary JSON',
            'filename': 'caseC_kturn_summary.json',
            'description': 'Query definition, funnel counts, representatives before geometry confirmation.'
        },
        {
            'label': 'Stage 2 confirmation JSON',
            'filename': 'caseC_geometry_confirmed.json',
            'description': 'MC-Annotate geometry checks for the representative structures.'
        },
        {
            'label': 'G.A-filtered candidate CSV',
            'filename': 'caseC_kturn_hits.csv',
            'description': 'Candidate rows after the G.A composition filter.'
        },
        {
            'label': 'Workflow figure PNG',
            'filename': 'caseC_kturn_figure.png',
            'description': 'Funnel figure used to explain the two-stage workflow.'
        },
        {
            'label': 'Workflow figure PDF',
            'filename': 'caseC_kturn_figure.pdf',
            'description': 'Publication-ready version of the workflow figure.'
        }
    ]

    return {
        'claim': summary.get('claim', ''),
        'query': summary.get('query', ''),
        'anchor': summary.get('anchor', ''),
        'boundary_note': summary.get('boundary_note', ''),
        'method': confirmed.get('method', ''),
        'stage_counts': stage_counts,
        'confirmed_representatives': confirmed_representatives,
        'downloads': downloads
    }

@app.route('/examples/kink-turn', methods=['GET'])
def kink_turn_workflow():
    """Render a worked example showing the two-stage kink-turn workflow."""
    try:
        workflow = load_kturn_workflow_data()
    except Exception as e:
        logger.error(f"Failed to load kink-turn workflow data: {e}")
        return render_template(
            'kink_turn_workflow.html',
            workflow=None,
            error_message='Kink-turn workflow data could not be loaded.'
        ), 500
    return render_template('kink_turn_workflow.html', workflow=workflow)


def load_named_motif_workflow(slug):
    """Load a precomputed named-motif worked example."""
    data_dir = NAMED_MOTIF_WORKFLOW_DIRS.get(slug)
    if not data_dir:
        raise KeyError(f"Unknown named-motif workflow: {slug}")

    summary_path = os.path.join(data_dir, 'summary.json')
    with open(summary_path, 'r', encoding='utf-8') as handle:
        workflow = json.load(handle)

    static_prefix = (
        f"data/named_motifs/{os.path.basename(data_dir)}"
    )
    for representative in workflow.get('confirmed_representatives', []):
        representative['structure_path'] = (
            f"{static_prefix}/motifs/{representative['structure']}"
        )
    for download in workflow.get('downloads', []):
        download['path'] = f"{static_prefix}/{download['filename']}"
    return workflow


@app.route('/examples', methods=['GET'])
def worked_examples():
    """Render the index of reviewer-requested named-motif examples."""
    examples = [
        {
            'slug': 'kink-turn',
            'title': 'Kink-turn',
            'eyebrow': 'Internal loop',
            'summary': (
                'Tandem non-WC G.A topology proposes candidates; stored '
                'geometry confirms selected kink-turns.'
            ),
            'metric': '4 confirmed representatives',
        }
    ]
    for slug in ('sarcin-ricin', 'gnra'):
        try:
            workflow = load_named_motif_workflow(slug)
        except Exception as exc:
            logger.error(f"Failed to load {slug} workflow data: {exc}")
            return render_template(
                'worked_examples.html',
                examples=[],
                error_message='Worked-example data could not be loaded.'
            ), 500
        final_stage = workflow['stage_counts'][-1]
        examples.append({
            'slug': slug,
            'title': workflow['title'],
            'eyebrow': (
                'Internal-loop graph'
                if slug == 'sarcin-ricin'
                else 'Four-nucleotide hairpin'
            ),
            'summary': workflow['subtitle'],
            'metric': (
                f"{final_stage['count']} FR3D-validated representatives"
            ),
        })
    return render_template(
        'worked_examples.html', examples=examples, error_message=None
    )


@app.route('/examples/<slug>', methods=['GET'])
def named_motif_workflow(slug):
    """Render a sarcin-ricin or GNRA worked example."""
    if slug not in NAMED_MOTIF_WORKFLOW_DIRS:
        abort(404)
    try:
        workflow = load_named_motif_workflow(slug)
    except Exception as exc:
        logger.error(f"Failed to load {slug} workflow data: {exc}")
        return render_template(
            'named_motif_workflow.html',
            workflow=None,
            error_message='Worked-example data could not be loaded.'
        ), 500
    return render_template(
        'named_motif_workflow.html',
        workflow=workflow,
        error_message=None,
    )


# Route to process custom motif search
@app.route('/process-custom-motif', methods=['POST'])
def process_custom_motif():
    """Process the user-drawn motif and search for similar patterns"""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    lock_file_path = os.path.join(base_dir, 'custom_search.lock')

    # Acquire cross-process file lock (protects against multiple gunicorn workers)
    lock_fd = None
    try:
        lock_fd = open(lock_file_path, 'w')
        fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except (IOError, OSError):
        if lock_fd is not None:
            lock_fd.close()
        return render_template('custom_results.html',
                               error_message="Another custom motif search is currently in progress. Please wait and try again in a few minutes.",
                               graph_info={'nodes': 0, 'edges': 0})

    # Also acquire threading lock (protects against threads within the same worker)
    if not custom_search_lock.acquire(blocking=False):
        fcntl.flock(lock_fd, fcntl.LOCK_UN)
        lock_fd.close()
        return render_template('custom_results.html',
                               error_message="Another custom motif search is currently in progress. Please wait and try again in a few minutes.",
                               graph_info={'nodes': 0, 'edges': 0})

    try:
        # Get the graph data from the form
        graph_data_str = request.form.get('graph_data')
        if not graph_data_str:
            return render_template('custom_results.html',
                                   error_message="No graph data received. Please draw a motif structure.",
                                   graph_info={'nodes': 0, 'edges': 0})

        # Parse the graph data
        graph_data = json.loads(graph_data_str)
        
        # Convert to NetworkX graph format
        target_graph = create_networkx_graph(graph_data)
        
        # Save target graph as temporary pickle file for find_hairpin.py
        base_dir = os.path.dirname(os.path.abspath(__file__))
        temp_target_file = os.path.join(base_dir, 'temp_target_graphs.pickle')
        
        # Create target graphs dictionary (find_hairpin.py expects this format)
        target_graphs = {
            'custom_motif_1': target_graph
        }
        
        with open(temp_target_file, 'wb') as f:
            pickle.dump(target_graphs, f)
        
        # Record start time
        start_time = time.time()
        
        # Call find_hairpin.py using subprocess
        result = run_find_hairpin()
        
        # Calculate processing time
        processing_time = round(time.time() - start_time, 2)
        
        if result['success']:
            # Read results from the database created by find_hairpin.py
            search_results = read_custom_search_results()
            
            # Calculate statistics
            unique_pdbs = list(set([result[2] for result in search_results])) if search_results else []
            graph_info = {
                'nodes': len(graph_data['nodes']),
                'edges': len(graph_data['edges'])
            }
            
            # Store results for download
            store_custom_results(search_results)

            if search_results:
                return render_template('custom_results.html',
                                       search_results=search_results,
                                       graph_data_str=graph_data_str,
                                       graph_info=graph_info,
                                       processing_time=processing_time,
                                       unique_pdbs=unique_pdbs,
                                       total_structures=result.get('total_structures', 'Unknown'))
            else:
                return render_template('custom_no_results.html',
                                       graph_data_str=graph_data_str,
                                       graph_info=graph_info,
                                       processing_time=processing_time)
        else:
            # Provide default graph_info even when there's an error
            graph_info = {
                'nodes': len(graph_data['nodes']),
                'edges': len(graph_data['edges'])
            }
            return render_template('custom_results.html',
                                   error_message=result['error'],
                                   graph_info=graph_info,
                                   processing_time=processing_time)
            
    except Exception as e:
        # Provide default graph_info even when there's an exception
        try:
            graph_info = {
                'nodes': len(graph_data['nodes']),
                'edges': len(graph_data['edges'])
            }
        except Exception:
            # If graph_data is not available, provide empty defaults
            graph_info = {
                'nodes': 0,
                'edges': 0
            }
        return render_template('custom_results.html',
                               error_message=f"An error occurred during processing: {str(e)}",
                               graph_info=graph_info)
    finally:
        # Always release both locks when done (success or failure)
        custom_search_lock.release()
        fcntl.flock(lock_fd, fcntl.LOCK_UN)
        lock_fd.close()

def create_networkx_graph(graph_data):
    """Convert the JSON graph data to NetworkX graph format"""
    G = nx.Graph()
    
    # Add nodes
    for node in graph_data['nodes']:
        G.add_node(node['id'])
    
    # Add edges with attributes
    for edge in graph_data['edges']:
        G.add_edge(edge['from'], edge['to'], attribute=edge['attribute'])
    
    return G

def run_find_hairpin():
    """Execute find_hairpin.py and return the result"""
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        find_hairpin_path = os.path.join(base_dir, 'find_hairpin.py')
        
        # Check if find_hairpin.py exists
        if not os.path.exists(find_hairpin_path):
            return {'success': False, 'error': 'find_hairpin.py not found'}
        
        # Run find_hairpin.py
        result = subprocess.run(['python', find_hairpin_path],
                                cwd=base_dir,
                                capture_output=True,
                                text=True,
                                timeout=1200)  # 20 minute timeout
        
        if result.returncode == 0:
            return {'success': True, 'output': result.stdout}
        else:
            return {'success': False, 'error': f"find_hairpin.py failed: {result.stderr}"}
            
    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'Search timed out. Please try with a simpler motif structure.'}
    except Exception as e:
        return {'success': False, 'error': f"Error running find_hairpin.py: {str(e)}"}

def read_custom_search_results():
    """Read the results from the database created by find_hairpin.py"""
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        # find_hairpin.py creates database directly in the website deploy root directory
        db_path = os.path.join(base_dir, 'hairpin.db')
        
        if not os.path.exists(db_path):
            return []
        
        conn = sqlite3.connect(db_path)
        try:
            cursor = conn.cursor()

            # Read all results from the files table
            cursor.execute("SELECT id, motif_type, pdbid, paired_nt_number, nt_number FROM files")
            results = cursor.fetchall()
        finally:
            conn.close()
        return results
        
    except Exception as e:
        print(f"Error reading custom search results: {e}")
        return []

def store_custom_results(search_results):
    """Store custom search results for download"""
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Save CSV for download
        if search_results:
            df = pd.DataFrame(search_results, columns=['ID', 'Motif Type', 'PDB ID', 'Nucleotide Range', 'Full Nucleotide List'])
            csv_path = os.path.join(base_dir, 'custom_search_results.csv')
            df.to_csv(csv_path, index=False, quoting=1)  # QUOTE_ALL
            
            # Create PDB files for download
            create_custom_pdb_files(search_results)
            
    except Exception as e:
        print(f"Error storing custom results: {e}")

def create_custom_pdb_files(search_results):
    """Create PDB files from the filecontent in the database"""
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        # find_hairpin.py creates database directly in the website deploy root directory
        db_path = os.path.join(base_dir, 'hairpin.db')
        
        # Create custom_pdb_files directory
        pdb_dir = os.path.join(base_dir, 'custom_pdb_files')
        if os.path.exists(pdb_dir):
            # Clear existing files
            for file in os.listdir(pdb_dir):
                os.remove(os.path.join(pdb_dir, file))
        else:
            os.makedirs(pdb_dir)
        
        # Read filecontent from database and create PDB files
        conn = sqlite3.connect(db_path)
        try:
            cursor = conn.cursor()
            for result in search_results:
                motif_id = result[0]
                cursor.execute("SELECT filecontent FROM files WHERE id = ?", (motif_id,))
                row = cursor.fetchone()
                if row and row[0]:
                    pdb_filename = os.path.join(pdb_dir, f'custom_match_{motif_id}.pdb')
                    with open(pdb_filename, 'w') as pdb_file:
                        pdb_file.write(row[0])
        finally:
            conn.close()
        
    except Exception as e:
        print(f"Error creating custom PDB files: {e}")

# Route to download custom search CSV results
@app.route('/download_custom_csv', methods=['GET'])
def download_custom_csv():
    """Download CSV file of custom search results"""
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        csv_path = os.path.join(base_dir, 'custom_search_results.csv')
        
        if not os.path.exists(csv_path):
            return "Custom search results not found. Please run a custom search first.", 404
        
        return send_file(csv_path, 
                         mimetype='text/csv', 
                         as_attachment=True, 
                         download_name='custom_motif_search_results.csv')
    except Exception as e:
        return f"An error occurred during CSV download: {str(e)}", 500

# Route to download custom search PDB files as ZIP
@app.route('/download_custom_pdb_zip', methods=['GET'])
def download_custom_pdb_zip():
    """Download ZIP file containing all matched PDB files from custom search"""
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        pdb_dir = os.path.join(base_dir, 'custom_pdb_files')

        if not os.path.exists(pdb_dir) or not os.listdir(pdb_dir):
            return "Custom PDB files not found. Please run a custom search first.", 404

        # Use temporary file instead of memory buffer to reduce memory usage
        temp_zip = tempfile.NamedTemporaryFile(delete=False, suffix='.zip')

        try:
            with zipfile.ZipFile(temp_zip.name, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                for root, dirs, files in os.walk(pdb_dir):
                    for file in files:
                        file_path = os.path.join(root, file)
                        arcname = os.path.relpath(file_path, start=pdb_dir)
                        zip_file.write(file_path, arcname=arcname)

            response = send_file(
                temp_zip.name,
                mimetype='application/zip',
                as_attachment=True,
                download_name='custom_motif_pdb_files.zip'
            )
            return response
        except Exception as inner_e:
            raise inner_e
        finally:
            if os.path.exists(temp_zip.name):
                os.unlink(temp_zip.name)
    except Exception as e:
        return f"An error occurred during ZIP download: {str(e)}", 500

# Route to get PDB content for 3D visualization
@app.route('/get_pdb_content/<int:motif_id>', methods=['GET'])
def get_pdb_content(motif_id):
    """Return PDB file content for a specific motif ID"""
    try:
        motif_type = request.args.get('type', 'standard')  # 'standard' or 'pseudoknot'

        conn = get_db_connection()
        try:
            cursor = conn.cursor()

            if motif_type == 'pseudoknot':
                # Query from PK table using rowid
                cursor.execute("SELECT file_content, pdbid, motif_type FROM PK WHERE rowid = ?", (motif_id,))
            else:
                # Query from data table
                cursor.execute("SELECT filecontent, pdbid, motif_type FROM data WHERE id = ?", (motif_id,))

            result = cursor.fetchone()
        finally:
            conn.close()

        if result:
            pdb_content = result[0]
            pdb_id = result[1]
            motif_type_name = result[2]

            return jsonify({
                'success': True,
                'pdb_content': pdb_content,
                'pdb_id': pdb_id,
                'motif_type': motif_type_name
            })
        else:
            return jsonify({
                'success': False,
                'error': 'Motif not found'
            }), 404

    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

# Route to download single PDB file
@app.route('/download_single_pdb/<int:motif_id>', methods=['GET'])
def download_single_pdb(motif_id):
    """Download a single PDB file for a specific motif"""
    try:
        motif_type = request.args.get('type', 'standard')  # 'standard' or 'pseudoknot'
        pdb_id = request.args.get('pdb_id', 'unknown')
        motif_name = request.args.get('motif_name', 'motif')

        conn = get_db_connection()
        try:
            cursor = conn.cursor()

            if motif_type == 'pseudoknot':
                # Query from PK table using rowid
                cursor.execute("SELECT file_content, pdbid, motif_type FROM PK WHERE rowid = ?", (motif_id,))
            else:
                # Query from data table
                cursor.execute("SELECT filecontent, pdbid, motif_type FROM data WHERE id = ?", (motif_id,))

            result = cursor.fetchone()
        finally:
            conn.close()

        if result:
            pdb_content = result[0]
            if not pdb_content:
                return "PDB content not available for this motif.", 404
            pdb_id_from_db = result[1]
            motif_type_from_db = result[2]

            # Generate filename: {PDB_ID}_{motif_type}_{ID}.pdb
            filename = f"{pdb_id_from_db}_{motif_type_from_db}_{motif_id}.pdb"

            # Return file as download
            return send_file(
                BytesIO(pdb_content.encode('utf-8')),
                mimetype='chemical/x-pdb',
                as_attachment=True,
                download_name=filename
            )
        else:
            return "Motif not found", 404

    except Exception as e:
        return f"An error occurred during download: {str(e)}", 500

if __name__ == '__main__':
    # Disable template caching in debug mode
    app.config['TEMPLATES_AUTO_RELOAD'] = True
    app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
    app.run(host='0.0.0.0', port=5000, debug=False)
    #app.run(debug=True)
