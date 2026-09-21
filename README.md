# RNAdex: 3D RNA Motif Library Website

A comprehensive web-based database for searching and visualizing RNA structural motifs from the Protein Data Bank.

**Last Updated:** July 2026

**Website:** [rnadex.dokhlab.org](https://rnadex.dokhlab.org)

## Features

- **Standard Motif Search**: Query pre-indexed RNA motifs
  - Motif types: Hairpin, Internal Loop, Bulge, Junction, Pseudoknot
  - Filter by nucleotide count
  - Quality checking based on C4' atom distances

- **Custom Motif Search**: Draw and search for custom structural patterns
  - Interactive graph-based motif drawing
  - Subgraph isomorphism matching
  - Networkx-powered pattern recognition

- **Reproducible Worked Examples**
  - Kink-turn candidate retrieval and geometry confirmation
  - Sarcin-ricin and canonical GNRA workflows with pinned FR3D validation
  - Downloadable candidates, audit samples, validation JSON, and structures

- **3D Structure Visualization**
  - Interactive 3Dmol.js viewer
  - PDB file downloads (single and batch)
  - Real-time structure rendering

- **Performance Optimizations**
  - Memory-efficient database queries for 4GB RAM systems
  - Result pagination (100 items per page)
  - Cursor iteration to prevent OOM crashes

## Database

**Public download:** `RNAdex.db` (8.6 GB, not included in repository)
- 462,755 RNA motif entries from 5,846 PDB structures
- SHA-256: `3d2356b6f82016d3d5a927a0d141a5803cb1ec45266c1d26e2607d075aecc9ed`
- Pseudoknot classification: LR, HHH, H, Hlout, LL, Hlin
- Download it separately and place it in the project root as `ATLAS.db`; this
  deployment filename is retained for internal compatibility.

## Installation

```bash
# Clone repository
git clone https://github.com/sesonli/ATLAS_website.git
cd ATLAS_website

# Install dependencies
pip install -r requirements.txt

# Development and browser-level responsive tests
pip install -r requirements-dev.txt

# Download RNAdex.db from https://rnadex.dokhlab.org/download-database

# Run application
python app.py
```

Visit http://localhost:5000

The worked-example index is available at
http://localhost:5000/examples.

## Requirements

- Python 3.8+
- Flask 2.3.3
- pandas 2.3.1
- networkx 3.1
- numpy 2.3.2
- 4GB RAM minimum (8GB recommended for custom searches)

## Project Structure

```
.
├── app.py                      # Flask web server
├── find_hairpin.py             # Subgraph isomorphism engine
├── generate_distribution_figures.py  # Data visualization
├── config.py                   # Configuration management
├── create_db_indexes.sql       # Database optimization script
├── ATLAS.db                    # Internal deployment name for RNAdex.db
├── data/rna_chain_corrected_2/ # Source PDB files
├── templates/                  # HTML templates
├── static/                     # CSS/JS/images
├── requirements.txt            # Python dependencies
└── app.log                     # Application log file
```

## Recent Improvements (November 2024)

### Code Quality
- Fixed Flask application duplicate initialization issue
- Enabled logging functionality with INFO level
- Added numpy dependency to requirements.txt
- Created centralized configuration management (config.py)

### Performance Optimizations
- Database index recommendations (create_db_indexes.sql)
- Memory-efficient cursor iteration for large result sets
- Result pagination (100 items per page)
- Thread lock for custom searches to prevent memory issues

### User Experience
- Loading indicators for search operations
- Home navigation buttons on all pages
- Improved error handling and user feedback
- Contact information section

## Configuration

The application can be configured via `config.py`:
- Set `FLASK_ENV=production` for production deployment
- Configure `SECRET_KEY` environment variable for security
- Adjust `MAX_SEARCH_RESULTS` for result limits
- Modify timeout settings for custom searches
- Set `ATLAS_DB_PATH` to stage a release database without replacing `ATLAS.db`.

## Database Optimization

To improve search performance, run the provided index creation script:
```bash
sqlite3 ATLAS.db < create_db_indexes.sql
```

## Documentation

- **CLAUDE.md**: Developer documentation and architecture guide
- **IT-Deployment-Notes.md**: Production deployment instructions
- **config.py**: Configuration settings and options

## Contact

**Nikolay V. Dokholyan, PhD**
Professor
- Department of Neurology
- Department of Neuroscience
- Department of Biomedical Engineering
- Department of Microbiology, Immunology, & Cancer Biology

University of Virginia
Email: dokh@virginia.edu

## License

RNAdex source code is available under the
[MIT License](https://opensource.org/license/mit). RNAdex database contents
and downloads are available under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

---

Generated with [Claude Code](https://claude.com/claude-code)
