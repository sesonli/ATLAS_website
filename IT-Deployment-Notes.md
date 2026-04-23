# ATLAS Website — Deployment Requirements (IT)

**Last Updated:** November 2024

## Target URL and scope
- Public URL: [atlas.dokhlab.org](http://atlas.dokhlab.org)
- Scope: Deploy the ATLAS web interface as a production service. This document states requirements and must‑haves; IT may choose the exact WSGI server and reverse proxy stack.

## Environment
- Python: 3.9–3.11 (virtual environment recommended)
- Dependencies: install from `requirements.txt` (includes Flask, pandas, networkx, numpy)
- OS: Linux server recommended

## Resources
- RAM: minimum 4 GB; recommended 8 GB (or more)
- CPU: ≥ 2 vCPU
- Disk: `ATLAS.db` is ~8.3 GB. Provide ≥ 20 GB free space for DB + outputs/logs

## Deployment package
- Deploy the entire folder as‑is (do not delete, move, or rename any files or subfolders). The application uses relative paths and expects the full layout intact.
- Key components (non‑exhaustive, for reference only):
  - `app.py`
  - `templates/` (HTML templates)
  - `static/` (CSS/JS/images)
  - `requirements.txt`
  - `ATLAS.db`
  - `data/rna_chain_corrected_2/`
- Runtime outputs (created/overwritten at runtime; ensure write permission in the deployment root):
  - `hairpin.db`
  - `custom_pdb_files/`
  - `pdb_files/`
  - `search_results.csv`, `custom_search_results.csv`

## Networking and TLS
- Expose the service at `atlas.dokhlab.org` with HTTPS (TLS certificate)
- Place behind a reverse proxy (Nginx/Apache). IT may choose the app server (e.g., gunicorn/uwsgi)

## Timeouts and long jobs
- Custom motif searches may take up to 60 minutes
- Ensure upstream/proxy and app‑server timeouts ≥ 3700 seconds end‑to‑end

## Runtime and paths
- Start the service with the deployment root as the current working directory; all paths are relative
- Run in production mode (no debug/reloader)

## Backups and persistence
- Back up `ATLAS.db` before initial deployment and prior to upgrades
- Runtime outputs (`hairpin.db`, CSVs, ZIPs, `custom_pdb_files/`, `pdb_files/`) are ephemeral and may be cleared

## Security and permissions
- Run under a non‑root service account
- Restrict write access to only the runtime output paths listed above
- No external secrets or environment tokens are required

## Basic verification (post‑deploy)
- `GET /` serves the homepage (“3D RNA Motif Library”)
- `GET /custom-motif-search` loads the drawing canvas (includes “Example 1”)
- Standard search returns a results table and enables CSV/ZIP download
- Custom search (using “Example 1”) completes and exposes CSV/ZIP within the configured timeout

## Performance Optimization (Optional)

### Database Indexing
- Run `create_db_indexes.sql` on `ATLAS.db` to create performance indexes
- This will significantly improve search query speeds
- Example: `sqlite3 ATLAS.db < create_db_indexes.sql`

### Configuration Management
- A `config.py` file is now available for centralized configuration
- Set environment variable `FLASK_ENV=production` for production settings
- Set `SECRET_KEY` environment variable for production deployments

### Logging
- Application now logs to `app.log` (INFO level by default)
- Monitor this file for application events and errors
- Consider log rotation for long-running deployments

### Cleanup Script (Optional)
- Temporary files in `pdb_files/` and `custom_pdb_files/` can be safely deleted periodically
- These directories are regenerated as needed from database content

## Reference
- Website: [atlas.dokhlab.org](http://atlas.dokhlab.org)
- Contact: Nikolay V. Dokholyan, PhD (dokh@virginia.edu)
