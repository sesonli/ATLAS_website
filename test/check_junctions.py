#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test script to check which multiway junctions exist in the ATLAS database
"""

import sqlite3
import os

def check_junctions():
    """Check all junction types in the database"""

    # Get database path
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    db_path = os.path.join(base_dir, 'ATLAS.db')

    if not os.path.exists(db_path):
        print(f"Error: Database not found at {db_path}")
        return

    print(f"Connecting to database: {db_path}\n")
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Check all junction types from 3-way to 8-way
    junction_types = [
        '3_way_junction',
        '4_way_junction',
        '5_way_junction',
        '6_way_junction',
        '7_way_junction',
        '8_way_junction'
    ]

    print("=" * 70)
    print("MULTIWAY JUNCTION ANALYSIS")
    print("=" * 70)

    total_junctions = 0

    for junction_type in junction_types:
        # Query the database
        query = "SELECT COUNT(*) FROM data WHERE motif_type = ?"
        cursor.execute(query, (junction_type,))
        count = cursor.fetchone()[0]

        total_junctions += count

        # Display results
        status = "[EXISTS]" if count > 0 else "[NOT FOUND]"
        print(f"{junction_type:20s}: {count:8,} entries  {status}")

        # If exists, show sample IDs
        if count > 0:
            cursor.execute("SELECT id, pdbid FROM data WHERE motif_type = ? LIMIT 3", (junction_type,))
            samples = cursor.fetchall()
            print(f"  Sample entries: {', '.join([f'{s[0]} (PDB: {s[1]})' for s in samples])}")

    print("=" * 70)
    print(f"TOTAL MULTIWAY JUNCTIONS: {total_junctions:,}")
    print("=" * 70)

    # Also check what's in the database with LIKE query
    print("\nSearching for any junction-related entries:")
    cursor.execute("SELECT DISTINCT motif_type FROM data WHERE motif_type LIKE '%junction%'")
    all_junction_types = cursor.fetchall()

    if all_junction_types:
        print(f"Found {len(all_junction_types)} junction type(s) in database:")
        for jtype in all_junction_types:
            cursor.execute("SELECT COUNT(*) FROM data WHERE motif_type = ?", (jtype[0],))
            count = cursor.fetchone()[0]
            print(f"  - {jtype[0]}: {count:,} entries")
    else:
        print("  No junction entries found")

    conn.close()

    print("\n" + "=" * 70)
    print("RECOMMENDATION:")
    print("=" * 70)

    # Provide recommendations
    existing_junctions = []
    for junction_type in junction_types:
        cursor = sqlite3.connect(db_path).cursor()
        cursor.execute("SELECT COUNT(*) FROM data WHERE motif_type = ?", (junction_type,))
        if cursor.fetchone()[0] > 0:
            junction_name = junction_type.replace('_', '-')
            existing_junctions.append(junction_name)
        cursor.close()

    if len(existing_junctions) > 0:
        print(f"Update search interface to include: {', '.join(existing_junctions)}")
        print(f"Update User Guide to reflect available junction types: {', '.join(existing_junctions)}")
    else:
        print("No multiway junctions found in database")

if __name__ == '__main__':
    check_junctions()
