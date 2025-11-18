# ATLAS: 3D RNA Motif Library Website

A comprehensive web-based database for searching and visualizing RNA structural motifs from the Protein Data Bank.

## Features

- **Standard Motif Search**: Query pre-indexed RNA motifs
  - Motif types: Hairpin, Internal Loop, Bulge, Junction, Pseudoknot
  - Filter by nucleotide count
  - Quality checking based on C4' atom distances

- **Custom Motif Search**: Draw and search for custom structural patterns
  - Interactive graph-based motif drawing
  - Subgraph isomorphism matching
  - Networkx-powered pattern recognition

- **3D Structure Visualization**
  - Interactive 3Dmol.js viewer
  - PDB file downloads (single and batch)
  - Real-time structure rendering

- **Performance Optimizations**
  - Memory-efficient database queries for 4GB RAM systems
  - Result pagination (100 items per page)
  - Cursor iteration to prevent OOM crashes

## Database

**ATLAS.db** (8.30 GB, not included in repository)
- 449,195 RNA motif records with quality-checked coordinates
- Pseudoknot classification: LR, HHH, H, Hlout, LL, Hlin
- Must be downloaded separately and placed in project root

## Installation

```bash
# Clone repository
git clone https://github.com/sesonli/ATLAS_website.git
cd ATLAS_website

# Install dependencies
pip install -r requirements.txt

# Download ATLAS.db (8.3 GB) - contact maintainer for access

# Run application
python app.py
```

Visit http://localhost:5000

## Requirements

- Python 3.8+
- Flask 2.3.3
- pandas 2.3.1
- networkx 3.1
- 4GB RAM minimum (8GB recommended for custom searches)

## Project Structure

```
.
├── app.py                      # Flask web server
├── find_hairpin.py             # Subgraph isomorphism engine
├── generate_distribution_figures.py  # Data visualization
├── ATLAS.db                    # Main database (download separately)
├── data/rna_chain_corrected_2/ # Source PDB files
├── templates/                  # HTML templates
├── static/                     # CSS/JS/images
└── requirements.txt
```

## Documentation

- **CLAUDE.md**: Developer documentation and architecture guide
- **IT-Deployment-Notes.md**: Production deployment instructions

## Contact

**Nikolay V. Dokholyan, PhD**
Professor, Neurology
UVA School of Medicine
Email: dokh@virginia.edu

## License

Copyright (c) 2024. All rights reserved.

---

Generated with [Claude Code](https://claude.com/claude-code)
