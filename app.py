# -*- coding: utf-8 -*-
"""
Created on Nov 25 11:40:07 2024
@author: JingyiLi
"""

from flask import Flask, render_template, request, send_file, jsonify
import sqlite3
import pandas as pd
import os
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
import tempfile
# import logging

# # Configure logging
# logging.basicConfig(
#     filename='app.log',        # Log file name
#     level=logging.DEBUG,       # Log level (DEBUG shows all details)
#     format='%(asctime)s [%(levelname)s] %(message)s'  # Log format
# )

# logger = logging.getLogger(__name__)

app = Flask(__name__)
app = Flask(__name__, static_folder='static', static_url_path='/static')

# Create lock for custom search to prevent concurrent executions (avoid memory issues on 4GB systems)
custom_search_lock = Lock()

# Connect to SQLite database
def get_db_connection():
    base_dir = os.path.dirname(__file__)
    db_path = os.path.join(base_dir, 'ATLAS.db')
    conn = sqlite3.connect(db_path)
    return conn

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
    lines = filecontent.strip().split('\n')

    for nt in nt_list:
        # Count atoms for this residue
        atom_count = sum(
            1 for line in lines
            if line.startswith('ATOM') and f'{line[21]}{line[22:27].strip()}' == nt
        )

        # Standard RNA nucleotides should have 15-23 atoms
        # Use conservative threshold of 10 to avoid false positives
        if atom_count < min_atoms_per_residue:
            return False

    return True

# Check C4' atom distances for coordinate quality
def check_c4_distances(filecontent, nt_list, min_dist=4.0, max_dist=8.0):
    """
    Check if C4' atom distances between adjacent nucleotides are within normal range.

    Args:
        filecontent: PDB file content string
        nt_list: List of nucleotide identifiers
        min_dist: Minimum acceptable distance in Angstroms (default 4.0)
        max_dist: Maximum acceptable distance in Angstroms (default 8.0)

    Returns:
        tuple: (is_normal: bool, reason: str)
            - (True, "Normal") if all distances are within range
            - (False, "C4' too close") if any distance < min_dist
            - (False, "Poor geometry") if any distance > max_dist
    """
    import numpy as np

    try:
        lines = filecontent.strip().split('\n')
        c4_atoms = []

        # Extract C4' atoms
        for line in lines:
            if not line.startswith('ATOM'):
                continue

            try:
                atom_name = line[12:16].strip()
                if atom_name != "C4'":
                    continue

                chain = line[21]
                res_num = line[22:27].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                full_id = f"{chain}{res_num}"
                c4_atoms.append({
                    'id': full_id,
                    'coords': np.array([x, y, z])
                })
            except (ValueError, IndexError):
                # Skip malformed lines
                continue

        # Need at least 2 C4' atoms to calculate distance
        if len(c4_atoms) < 2:
            return (True, "Normal")  # Cannot check, assume normal

        # Sort by nucleotide list order (to maintain sequence)
        id_to_atom = {atom['id']: atom for atom in c4_atoms}
        ordered_c4_atoms = []
        for nt in nt_list:
            if nt in id_to_atom:
                ordered_c4_atoms.append(id_to_atom[nt])

        # Calculate distances between adjacent C4' atoms
        for i in range(len(ordered_c4_atoms) - 1):
            dist = np.linalg.norm(ordered_c4_atoms[i]['coords'] - ordered_c4_atoms[i+1]['coords'])

            if dist < min_dist:
                return (False, "C4' too close")
            elif dist > max_dist:
                return (False, "Poor geometry")

        return (True, "Normal")

    except Exception as e:
        # On any error, assume normal to avoid false positives
        return (True, "Normal")

# Search RNA motifs based on motif_type and nt_number or subtype
def search_motif(motif_type, nt_number=None, max_results=10000):
    """
    Search for RNA motifs with memory-efficient cursor iteration.

    Args:
        motif_type: Type of motif to search for
        nt_number: Number of nucleotides (for hairpin/internal/bulge)
        max_results: Maximum number of results to return (default 5000)

    Returns:
        tuple: (filtered_rows, csv_path, zip_buffer, metadata)
        metadata includes 'truncated' flag and 'total_scanned'
    """
    try:
        conn = get_db_connection()
        cursor = conn.cursor()

        filtered_rows = []
        was_truncated = False
        total_scanned = 0

        if motif_type in ["hairpin", "internal", "bulge"]:
            # Query for hairpin, internal, and bulge from "data" table
            # Use cursor iteration instead of fetchall() to avoid loading all rows into memory
            query = "SELECT id, motif_type, pdbid, nt_number, filecontent FROM data WHERE motif_type = ?"
            cursor.execute(query, (motif_type,))

            # Iterate row by row instead of fetchall()
            for row in cursor:
                total_scanned += 1
                # Filter rows by nt_number (length of nt_number list)
                # AND exclude records with invalid filecontent
                if len(row[3].split(',')) == int(nt_number):
                    # Skip records with "No matching nodes found" error
                    if row[4] and 'No matching nodes found' not in row[4]:
                        # Check residue completeness
                        is_atom_complete = is_residue_data_complete(row[4], row[3].split(','))
                        # Check C4' distance quality
                        is_coord_normal, coord_reason = check_c4_distances(row[4], row[3].split(','))
                        # Combine quality checks
                        is_quality_ok = is_atom_complete and is_coord_normal
                        warning_reason = ""
                        if not is_atom_complete:
                            warning_reason = "Incomplete atom data"
                        if not is_coord_normal:
                            warning_reason += (" + " if warning_reason else "") + coord_reason
                        # If both checks pass, set reason to empty string
                        if is_quality_ok:
                            warning_reason = ""
                        # Add quality flags to the row tuple
                        row_with_flag = list(row) + [is_quality_ok, warning_reason]
                        filtered_rows.append(tuple(row_with_flag))
                        # Stop if we reach max_results
                        if len(filtered_rows) >= max_results:
                            was_truncated = True
                            break

        elif 'junction' in motif_type:
            # Query for multiway junctions
            motif_type_db = motif_type.replace("-way", "").replace(" ", "_")
            query = "SELECT id, motif_type, pdbid, nt_number, filecontent FROM data WHERE motif_type = ?"
            cursor.execute(query, (motif_type_db,))

            # Use cursor iteration with limit
            for row in cursor:
                total_scanned += 1
                # Skip records with "No matching nodes found" error
                if row[4] and 'No matching nodes found' not in row[4]:
                    # Check residue completeness
                    is_atom_complete = is_residue_data_complete(row[4], row[3].split(','))
                    # Check C4' distance quality
                    is_coord_normal, coord_reason = check_c4_distances(row[4], row[3].split(','))
                    # Combine quality checks
                    is_quality_ok = is_atom_complete and is_coord_normal
                    warning_reason = ""
                    if not is_atom_complete:
                        warning_reason = "Incomplete atom data"
                    if not is_coord_normal:
                        warning_reason += (" + " if warning_reason else "") + coord_reason
                    if is_quality_ok:
                        warning_reason = ""
                    # Add quality flags to the row tuple
                    row_with_flag = list(row) + [is_quality_ok, warning_reason]
                    filtered_rows.append(tuple(row_with_flag))
                    if len(filtered_rows) >= max_results:
                        was_truncated = True
                        break

        elif 'pseudoknot' in motif_type:
            # Query for pseudoknot subtypes from "PK" table
            # Note: PK table's id column is NULL, so we use rowid instead
            motif_type_db = motif_type
            query = "SELECT rowid, motif_type, pdbid, nt_number, file_content FROM PK WHERE motif_type = ?"
            cursor.execute(query, (motif_type_db.strip(),))

            # Use cursor iteration with limit
            for row in cursor:
                total_scanned += 1
                # Skip records with "No matching nodes found" error
                # Note: PK table uses file_content (with underscore)
                if row[4] and 'No matching nodes found' not in row[4]:
                    # Check residue completeness
                    is_atom_complete = is_residue_data_complete(row[4], row[3].split(','))
                    # Check C4' distance quality
                    is_coord_normal, coord_reason = check_c4_distances(row[4], row[3].split(','))
                    # Combine quality checks
                    is_quality_ok = is_atom_complete and is_coord_normal
                    warning_reason = ""
                    if not is_atom_complete:
                        warning_reason = "Incomplete atom data"
                    if not is_coord_normal:
                        warning_reason += (" + " if warning_reason else "") + coord_reason
                    if is_quality_ok:
                        warning_reason = ""
                    # Add quality flags to the row tuple
                    row_with_flag = list(row) + [is_quality_ok, warning_reason]
                    filtered_rows.append(tuple(row_with_flag))
                    if len(filtered_rows) >= max_results:
                        was_truncated = True
                        break

        else:
            # If motif_type is unknown
            filtered_rows = []

        conn.close()

        # If results are found, save them and return
        if filtered_rows:
            # Save results to CSV (now includes quality checks)
            df = pd.DataFrame(filtered_rows, columns=['ID', 'Motif Type', 'PDB ID', 'NT Number', 'Filecontent', 'Quality_OK', 'Warning_Reason'])
            base_dir = os.path.dirname(os.path.abspath(__file__))
            csv_path = os.path.join(base_dir, 'search_results.csv')
            df[['ID', 'Motif Type', 'PDB ID', 'NT Number', 'Quality_OK', 'Warning_Reason']].to_csv(csv_path, index=False)

            # Save PDB files
            pdb_dir = os.path.join(base_dir, 'pdb_files')
            if os.path.exists(pdb_dir):
                for file in os.listdir(pdb_dir):
                    os.remove(os.path.join(pdb_dir, file))
            else:
                os.makedirs(pdb_dir)

            for row in filtered_rows:
                motif_id = row[0]
                pdb_content = row[4]  # Filecontent column
                pdb_filename = os.path.join(pdb_dir, f'{motif_id}.pdb')
                with open(pdb_filename, 'w') as pdb_file:
                    pdb_file.write(pdb_content)

            # Create a ZIP archive
            zip_buffer = BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                for file in os.listdir(pdb_dir):
                    file_path = os.path.join(pdb_dir, file)
                    zip_file.write(file_path, arcname=file)
            zip_buffer.seek(0)

            # Return with metadata
            metadata = {
                'truncated': was_truncated,
                'total_scanned': total_scanned,
                'result_count': len(filtered_rows),
                'max_results': max_results
            }
            return filtered_rows, csv_path, zip_buffer, metadata

        # No results found
        return [], None, None, {'truncated': False, 'total_scanned': total_scanned, 'result_count': 0}

    except sqlite3.Error as e:
        print(f"Error during search: {e}")
        return [], None, None, {'error': str(e)}

# Route to the home page
@app.route('/')
def home():
    return render_template('index.html')

# User Guide page
@app.route('/user-guide', methods=['GET'])
def user_guide():
    return render_template('user_guide.html')

# Route to handle motif search requests
@app.route('/search', methods=['POST'])
def search():
    motif_type = request.form.get('motif_type')
    nt_number = request.form.get('nt_number')
    if motif_type == 'multiway-junction':
        # Here nt_number will be a string like '3-way junction'
        motif_type = nt_number
        nt_number_count = None
    elif motif_type == 'pseudoknot':
        motif_type = f'pseudoknot_{nt_number}'
        nt_number_count = None
    else:
        # For other types, it's a number, so convert to integer
        try:
            nt_number_count = int(nt_number)
        except ValueError:
            nt_number_count = None

    search_results, csv_path, pdb_zip_buffer, metadata = search_motif(motif_type, nt_number_count)

    if search_results:
        # Prepare warning message if results were truncated
        warning_message = None
        if metadata.get('truncated'):
            warning_message = f"Results limited to {metadata['result_count']:,} out of {metadata['total_scanned']:,} total matches (maximum {metadata['max_results']:,} results shown)."

        return render_template('results.html',
                               search_results=search_results,
                               csv_path=csv_path,
                               warning=warning_message,
                               metadata=metadata)
    else:
        return render_template('no_results.html')
    
# Route to download CSV file
@app.route('/download_csv', methods=['GET'])
def download_csv():
    csv_path = request.args.get('csv_path')

    if not csv_path or not os.path.exists(csv_path):
        return "CSV file not found.", 404

    try:
        return send_file(csv_path, mimetype='text/csv', as_attachment=True, download_name='search_results.csv')
    except Exception as e:
        return f"An error occurred during CSV download: {str(e)}", 500

# Route to download ZIP of PDB file
@app.route('/download_zip', methods=['GET'])
def download_zip():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_dir = os.path.join(base_dir, 'pdb_files')

    if not os.path.exists(pdb_dir) or not os.listdir(pdb_dir):
        return "PDB files not found.", 404

    # Use temporary file instead of memory buffer to reduce memory usage
    temp_zip = tempfile.NamedTemporaryFile(delete=False, suffix='.zip')

    try:
        # Zip the files to temporary file
        with zipfile.ZipFile(temp_zip.name, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            for root, dirs, files in os.walk(pdb_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, start=pdb_dir)
                    zip_file.write(file_path, arcname=arcname)

        # Send the temporary file
        return send_file(
            temp_zip.name,
            mimetype='application/zip',
            as_attachment=True,
            download_name='pdb_files.zip'
        )
    except Exception as e:
        # Clean up temp file on error
        if os.path.exists(temp_zip.name):
            os.unlink(temp_zip.name)
        return f"An error occurred during ZIP creation: {str(e)}", 500

# Route to download the complete database
@app.route('/download-database', methods=['GET'])
def download_database():
    """Download the entire ATLAS database file"""
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
        db_path = os.path.join(base_dir, 'ATLAS.db')
        
        if not os.path.exists(db_path):
            return "Database file not found.", 404
        
        return send_file(
            db_path,
            as_attachment=True,
            download_name='ATLAS.db',
            mimetype='application/x-sqlite3'
        )
    except Exception as e:
        return f"An error occurred during database download: {str(e)}", 500

# Route to custom motif search drawing page
@app.route('/custom-motif-search', methods=['GET'])
def custom_motif_search():
    """Render the custom motif drawing interface"""
    return render_template('custom_motif_draw.html')

# Route to process custom motif search
@app.route('/process-custom-motif', methods=['POST'])
def process_custom_motif():
    """Process the user-drawn motif and search for similar patterns"""
    # Try to acquire the lock (non-blocking)
    if not custom_search_lock.acquire(blocking=False):
        # Another custom search is already running
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
            
            return render_template('custom_results.html',
                                   search_results=search_results,
                                   graph_info=graph_info,
                                   processing_time=processing_time,
                                   unique_pdbs=unique_pdbs,
                                   total_structures=result.get('total_structures', 'Unknown'))
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
        except:
            # If graph_data is not available, provide empty defaults
            graph_info = {
                'nodes': 0,
                'edges': 0
            }
        return render_template('custom_results.html',
                               error_message=f"An error occurred during processing: {str(e)}",
                               graph_info=graph_info)
    finally:
        # Always release the lock when done (success or failure)
        custom_search_lock.release()

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
        cursor = conn.cursor()
        
        # Read all results from the files table
        cursor.execute("SELECT id, motif_type, pdbid, paired_nt_number, nt_number FROM files")
        results = cursor.fetchall()
        
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
            df = pd.DataFrame(search_results, columns=['ID', 'Motif Type', 'PDB ID', 'Nucleotide Range', 'Full NT Sequence'])
            csv_path = os.path.join(base_dir, 'custom_search_results.csv')
            df.to_csv(csv_path, index=False)
            
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
        cursor = conn.cursor()
        
        for result in search_results:
            motif_id = result[0]
            cursor.execute("SELECT filecontent FROM files WHERE id = ?", (motif_id,))
            row = cursor.fetchone()
            
            if row and row[0]:
                pdb_filename = os.path.join(pdb_dir, f'custom_match_{motif_id}.pdb')
                with open(pdb_filename, 'w') as pdb_file:
                    pdb_file.write(row[0])
        
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

            return send_file(
                temp_zip.name,
                mimetype='application/zip',
                as_attachment=True,
                download_name='custom_motif_pdb_files.zip'
            )
        except Exception as inner_e:
            # Clean up temp file on error
            if os.path.exists(temp_zip.name):
                os.unlink(temp_zip.name)
            raise inner_e
    except Exception as e:
        return f"An error occurred during ZIP download: {str(e)}", 500

# Route to get PDB content for 3D visualization
@app.route('/get_pdb_content/<int:motif_id>', methods=['GET'])
def get_pdb_content(motif_id):
    """Return PDB file content for a specific motif ID"""
    try:
        motif_type = request.args.get('type', 'standard')  # 'standard' or 'pseudoknot'

        conn = get_db_connection()
        cursor = conn.cursor()

        if motif_type == 'pseudoknot':
            # Query from PK table using rowid
            cursor.execute("SELECT file_content, pdbid, motif_type FROM PK WHERE rowid = ?", (motif_id,))
        else:
            # Query from data table
            cursor.execute("SELECT filecontent, pdbid, motif_type FROM data WHERE id = ?", (motif_id,))

        result = cursor.fetchone()
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
        cursor = conn.cursor()

        if motif_type == 'pseudoknot':
            # Query from PK table using rowid
            cursor.execute("SELECT file_content, pdbid, motif_type FROM PK WHERE rowid = ?", (motif_id,))
        else:
            # Query from data table
            cursor.execute("SELECT filecontent, pdbid, motif_type FROM data WHERE id = ?", (motif_id,))

        result = cursor.fetchone()
        conn.close()

        if result:
            pdb_content = result[0]
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
    app.run(host='0.0.0.0', port=5000, debug=False)
