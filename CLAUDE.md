# CLAUDE.md

Guidance for Claude Code when working with this repository.

## Project Overview

ATLAS is a Flask web app for searching 3D RNA structural motifs in the PDB. Two search modes: standard (type-based) and custom (user-drawn graph).

## Commands

```bash
python app.py                  # Development
pip install -r requirements.txt
```

Production: use gunicorn/uwsgi with timeouts ≥ 3700 seconds (see IT-Deployment-Notes.md).

## Architecture

**app.py** — Flask server, two routes:
- `/search` [POST]: queries ATLAS.db, returns motif matches
- `/process-custom-motif` [POST]: accepts NetworkX graph JSON, runs `find_hairpin.py` via subprocess

**find_hairpin.py** — subgraph isomorphism engine using NetworkX GraphMatcher; reads `batch_0000_graphs.pickle`, writes matches to `hairpin.db`

### Database Schema

**ATLAS.db** (8.3 GB, read-only):
- `data` table: `id, motif_type, pdbid, nt_number, filecontent`
- `PK` table: `id, motif_type, pdbid, nt_number, file_content` ← note underscore (inconsistency)

**hairpin.db** (runtime-generated):
- `files` table: `id, motif_type, pdbid, paired_nt_number, nt_number, filecontent, created_at`

### Key Routes

- `/` — index.html
- `/user-guide` — user_guide.html
- `/custom-motif-search` [GET] — custom_motif_draw.html
- `/download_csv`, `/download_zip` — standard search downloads
- `/download_custom_csv`, `/download_custom_pdb_zip` — custom search downloads
- `/download-database` — full ATLAS.db download

### Important Constraints

- **Do not rename/move files**: all paths are relative to deployment root
- **DB column inconsistency**: `data.filecontent` vs `PK.file_content`
- **Graph pickle format**: must be `{graph_id: networkx.Graph}` with edge `attribute` keys
- **PDB filename convention**: `{pdbid}.pdb` in `data/rna_chain_corrected_2/`
- **Max results**: 5000 per standard search (memory guard)
- **Custom search**: thread-locked, only one at a time; can take up to 60 min

### Node ID Formats

Handles: `A1`, `A-1`, `'0'1`, `'0'-1` (chain + residue, including negative numbering and quoted chain IDs).

### Runtime Write Permissions Required

`hairpin.db`, `custom_pdb_files/`, `pdb_files/`, `*.csv`, `temp_target_graphs.pickle`

### Memory Notes

Standard search uses cursor iteration (not `fetchall()`) to avoid loading all records into RAM. ZIP generation uses temp files, not BytesIO. Custom search is still memory-intensive (NP-complete subgraph isomorphism) — recommend 8 GB RAM.
