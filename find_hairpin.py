# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 16:48:20 2024
I corrected the subgraph searching algorithm for hairpin. It's fast now. I stored the data to database. 
I dinnot use the compress_graph and decompress_graph to search for hairpinloops. 
Finding matching graphs call is subgraph isomorphic and save to data base, which call extract node info, which call read pdb and find nt.
I fix the noncontinues chains problem. A1 A2 A4 A5 is considered as two chains.  
@author: JingyiLi

Finding matching graphs
    /                         \
save to database     is subgraph isomorphic
    /
extract node info
    /
read pdb and find nt
"""

from networkx.algorithms import isomorphism
import time
import pickle
import os
import sqlite3
import re
import networkx as nx

def get_chain_type(node_name):
    # 1. Match uppercase letter followed by a positive integer, e.g., 'A1'
    if re.match(r'^[A-Z]\d+$', node_name):
        return 'A'  # Return 'A' as chain type

    # 2. Match uppercase letter followed by a hyphen and a positive integer, e.g., 'A-1'
    elif re.match(r'^[A-Z]-\d+$', node_name):
        return 'A-'  # Return 'A-' as chain type

    # 3. Match a positive integer enclosed in single quotes followed by a hyphen and another positive integer, e.g., "'0'-1"
    elif re.match(r"^'\d+'-\d+$", node_name):
        return "0'-"  # Return "0'-" as chain type

    # 4. Match a positive integer enclosed in single quotes followed directly by another positive integer, e.g., "'0'1"
    elif re.match(r"^'\d+'\d+$", node_name):
        return "0'"  # Return "0'" as chain type

    else:
        return None

def natural_key(node_name):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', node_name)]

def is_paired_node(G, node):
    for neighbor in G.neighbors(node):
        if G[node][neighbor]['attribute'] != [0, 0, 1]:
            return True
    return False

def is_unpaired_node(G, node):
    for neighbor in G.neighbors(node):
        if G[node][neighbor]['attribute'] != [0, 0, 1]:
            return False
    return True

def get_chains(lst):
    chains_dict = {}
    current_chain = [lst[0]] 
    def parse_item(item):
        match = re.match(r"([A-Za-z]*)'?-?(\d+)", item)
        if match:
            letter_part = match.group(1)
            number_part = int(match.group(2)) 
            return letter_part, number_part
        return None, None
    chain_index = 1 
    for i in range(1, len(lst)):
        prev_letter, prev_number = parse_item(lst[i-1])
        curr_letter, curr_number = parse_item(lst[i])

        if prev_letter == curr_letter and curr_number == prev_number + 1:
            current_chain.append(lst[i])
        else:
            chains_dict[f"Chain {chain_index}"] = current_chain  
            current_chain = [lst[i]]  
            chain_index += 1  
    chains_dict[f"Chain {chain_index}"] = current_chain

    #print(f"Generated chains: {chains_dict}")
    return chains_dict

def read_pdb_and_find_nt(pdb_filepath, node_id):
    '''input the pdb filename, and the node ids. Open the pdb file to get the lines representing 
    the atom info.
    node_id A1 A-1 '0'10 '0'-10
    pdb file name: pdb4wfn.ent.pdb (example)
    ATOM     53  H5'   C A   2       9.545  18.397   7.203  1.00  0.00           H  
    '''
    node_id = node_id.replace("'", "")
    content = []
    with open(pdb_filepath, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                current_chain = line[21:22].strip()  # Adjusted to get correct chain ID
                current_nt = line[22:27].strip()  # Adjusted to get correct nt number
                
                # Combine current chain and nt number, e.g., 'A1', 'A-1', '01'
                node_id_in_file = current_chain + current_nt
                
                # Compare with the input node_id
                if node_id_in_file == node_id:
                    content.append(line.strip())
    
    return content

def extract_node_info(node_id, pdb_filename):
    '''get the pdb filepath with pdb name
    Then call read pdb and output matching atomic info for nodes'''
    # Updated to use the correct PDB path for our project
    # Extract PDB ID from filename (e.g., "157D_annotation.txt" -> "157D")
    pdb_id = pdb_filename.replace("_annotation.txt", "")
    
    # Use the correct path to PDB files
    # Use rna_chain_corrected_2 directory within deployment folder
    current_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_filepath = os.path.join(current_dir, 'data', 'rna_chain_corrected_2', f'{pdb_id}.pdb')
    
    if not os.path.exists(pdb_filepath):
        return f"File {pdb_id}.pdb not found in rna_chain_corrected_2 folder."
    
    content = read_pdb_and_find_nt(pdb_filepath, node_id)
    return "\n".join(content) if content else "No matching nodes found"

def save_to_database(motif_type, graph_id, reduced_nt_numbers, combined_nt_numbers, combined_filecontent, conn):
    cursor = conn.cursor()
    pdbid = graph_id[:4]  # Fixed PDB ID extraction, get PDB ID from first 4 characters
    cursor.execute('''
    INSERT INTO files (motif_type, pdbid, paired_nt_number, nt_number, filecontent)
    VALUES (?, ?, ?, ?,?)
    ''', (motif_type, pdbid, reduced_nt_numbers, combined_nt_numbers, combined_filecontent))
    conn.commit()
    print(f"Saved to database: {motif_type}, {pdbid}")

def is_subgraph_isomorphic(graph, target_graph):
    """
    Check if the graph has isomorphic subgraphs. Output the mapping between 
    the node numbers of subgraphs and target graphs. Also, return the matched subgraphs from the original graph.
    The edges with different attributes are treated as different edges.
    """
    GM = isomorphism.GraphMatcher(graph, target_graph, 
                                  node_match=lambda n1,n2: True,
                                  edge_match=lambda e1, e2: e1['attribute'] == e2['attribute'])
    
    matches = []
    matched_subgraphs = []  # Use a set to avoid duplicates

    for mapping in GM.subgraph_isomorphisms_iter():
        # Standardize the mapping to avoid duplicates
        sorted_mapping = tuple(sorted(mapping.items()))  # Sort the mapping items

        if sorted_mapping not in matches:
            matches.append(sorted_mapping)
            # Extract the matched subgraph based on the current mapping
            subgraph_nodes = mapping.keys()
            matched_subgraph = graph.subgraph(subgraph_nodes).copy()
            
            # Convert edges with attributes to a hashable format
            def make_hashable(attr):
                """ Convert lists in attribute values to tuples to make them hashable """
                return {k: tuple(v) if isinstance(v, list) else v for k, v in attr.items()}
            
            sorted_edges = tuple(sorted((u, v) for u, v, attr in matched_subgraph.edges(data=True)))
            
            matched_subgraphs.append(sorted_edges)  # Add to set
    
    if matches:
        # Convert matched_subgraphs back to list of subgraphs
        unique_subgraphs = [graph.edge_subgraph(edges) for edges in matched_subgraphs]
        return matches, unique_subgraphs
    return None, None

def sort_key(node):
    # Library node ids are strings ('A1', 'A-1', "'0'1"), but a browser-drawn
    # target graph numbers its nodes with plain ints, and those ints reach here
    # through the mapping sort below. Coerce so both node id kinds are accepted.
    node = str(node).replace('-', '')
    match = re.match(r"([A-Za-z]+)(\d+)", node)
    if not match:
        match = re.match(r"'(\d+)'(\d+)", node)  
    if match:
        letter = match.group(1)
        number = int(match.group(2))  
        return (letter, number)  
    return ('', 0) 

def find_matching_subgraphs(graphs, target_graphs, conn):
    matching_counts = {}
    target_match_counts = {tg_id: 0 for tg_id in target_graphs.keys()}
    total_graphs = len(graphs)
    graph_processed = 0
    matches_info = []
    problematic_graphs = []
    for graph_id, graph in graphs.items():   
        graph_processed += 1
        count = 0
        start_time = time.time()        
        for target_graph_id, target_graph in target_graphs.items():
            mappings, matched_subgraphs = is_subgraph_isomorphic(graph, target_graph)
            # Decompress the first matched subgraph (as an example)
            if matched_subgraphs:
                for mapping, matched_graph in zip(mappings, matched_subgraphs):
                    combined_filecontent = []
                    reduced_nt_numbers = []
                    combined_nt_numbers = []
                    mapping = sorted(mapping, key=lambda x: sort_key(x[1]))
                    for a_node, t_node in mapping:
                        reduced_nt_numbers.append(a_node[:2])
                    for nodes in matched_graph:
                        nt_number = nodes
                        filecontent = extract_node_info(nt_number, graph_id)  # return atom data to file content
                        combined_nt_numbers.append(nt_number)
                        combined_filecontent.append(filecontent)
                    combined_nt_numbers = sorted(combined_nt_numbers, key=sort_key)
                    # according to the number, get a new sequence from little to big    
                    reduced_nt_numbers = [combined_nt_numbers[0], combined_nt_numbers[-1]]
                    combined_nt_numbers_str = ','.join(combined_nt_numbers)
                    reduced_nt_numbers_str = ','.join(reduced_nt_numbers)
                    combined_filecontent_str = '\n'.join(combined_filecontent)
                    
                    save_to_database('user defined', graph_id, reduced_nt_numbers_str, combined_nt_numbers_str, combined_filecontent_str, conn)
                    
        matching_counts[graph_id] = count
        elapsed_time = time.time() - start_time

        #print(f"Finished searching graph {graph_id} ({graph_processed}/{total_graphs}). Time spent: {elapsed_time:.2f} seconds.")
   
    print(f"Total graphs processed: {graph_processed}/{total_graphs}")
    print(f"Problematic graphs: {problematic_graphs}")
    
    return matching_counts, matches_info, target_match_counts, problematic_graphs
# Setup directories
current_dir = os.path.dirname(os.path.abspath(__file__))
db_path = os.path.join(current_dir, 'hairpin.db')


def main():
    """Run a custom motif search. app.py invokes this file as a subprocess.

    Kept behind a __main__ guard so the helpers above can be imported by tests
    without dropping the results table and launching a full search on import.
    """
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute('DROP TABLE IF EXISTS files')
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS files (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            motif_type TEXT,
            pdbid TEXT,
            paired_nt_number TEXT,
            nt_number TEXT,
            filecontent TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
        ''')
        # Use batch_0000_graphs.pickle as database. Resolve inputs against the
        # script directory so the search does not depend on the caller's cwd.
        file_path = os.path.join(current_dir, "batch_0000_graphs.pickle")
        with open(file_path, 'rb') as f:
            graphs = pickle.load(f)

        # Load target graphs from temporary file created by web interface
        try:
            with open(os.path.join(current_dir, "temp_target_graphs.pickle"), 'rb') as f:
                target_graphs = pickle.load(f)
            print(f"Loaded {len(target_graphs)} target graphs from temp file")
        except FileNotFoundError:
            target_graphs = {}  # fallback if no temp file exists
            print("No temp target graphs file found, using empty target graphs")
        print(f"Total number of graphs: {len(graphs)}")
        time1 = time.time()
        matching_counts, matches_info, target_match_counts, problematic_graphs = find_matching_subgraphs(graphs, target_graphs, conn)
        time2 = time.time()

        print("Matching counts per graph:", matching_counts)
        print("Matching counts per target graph:", target_match_counts)

        file_path_3 = os.path.join(current_dir, "Matched_hairpin_test.txt")
        with open(file_path_3, 'w') as file:
            for graph_id, count in matching_counts.items():
                file.write(f'Graph ID: {graph_id}, Count: {count}\n')

            for target_graph_id, match_count in target_match_counts.items():
                file.write(f"Target Graph ID: {target_graph_id}, Total Matches: {match_count}\n")

            file.write(f"Total time is: {time2 - time1} seconds\n")
            for match_info in matches_info:
                file.write(match_info)


if __name__ == '__main__':
    main()
