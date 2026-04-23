import sqlite3

conn = sqlite3.connect('ATLAS.db')
cursor = conn.cursor()

# Check hairpin count
cursor.execute('SELECT COUNT(*) FROM data WHERE motif_type="hairpin"')
print('Total hairpin records:', cursor.fetchone()[0])

# Check sample records
cursor.execute('SELECT id, motif_type, nt_number FROM data WHERE motif_type="hairpin" LIMIT 5')
print('\nSample hairpin records (id, motif_type, nt_number):')
for row in cursor.fetchall():
    nt_count = len(row[2].split(','))
    print(f"  ID: {row[0]}, Type: {row[1]}, NT: {row[2][:50]}..., Count: {nt_count}")

# Test search for hairpin with 5 nucleotides (display value 3 + 2)
cursor.execute('SELECT id, motif_type, nt_number FROM data WHERE motif_type="hairpin"')
print('\nSearching for hairpin with 5 nucleotides:')
count = 0
for row in cursor:
    if len(row[2].split(',')) == 5:
        count += 1
        if count <= 3:
            print(f"  Found: ID {row[0]}, NT: {row[2][:50]}...")
    if count >= 3:
        break

print(f'\nTotal found with 5 nucleotides: checking first 100...')
cursor.execute('SELECT id, motif_type, nt_number FROM data WHERE motif_type="hairpin" LIMIT 100')
count = 0
for row in cursor.fetchall():
    if len(row[2].split(',')) == 5:
        count += 1
print(f'Found {count} hairpins with 5 nucleotides in first 100 records')

conn.close()
