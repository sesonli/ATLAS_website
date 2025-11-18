# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ATLAS is a 3D RNA Motif Library web application that enables researchers to search for RNA structural motifs in the Protein Data Bank. The application provides both standard motif searches and custom motif drawing/searching capabilities.

## Commands

### Running the Application
```bash
# Development mode (Flask debug server)
python app.py

# Production deployment (see IT-Deployment-Notes.md for production requirements)
# Use gunicorn or uwsgi with timeouts ≥ 3700 seconds
```

### Dependencies
```bash
pip install -r requirements.txt
```

Required packages: Flask 2.3.3, pandas 2.3.1, networkx 3.1

## Architecture

### Core Components

**app.py** - Flask web server with two primary search modes:
1. **Standard Search** (`/search`): Queries pre-indexed motifs from ATLAS.db
   - Motif types: hairpin, internal, bulge, multiway junctions, pseudoknots
   - Returns results from the `data` table (or `PK` table for pseudoknots)
   - Generates CSV and ZIP downloads of matching PDB files

2. **Custom Motif Search** (`/process-custom-motif`): User-drawn graph matching
   - Accepts NetworkX graph JSON from the web interface
   - Saves to `temp_target_graphs.pickle`
   - Invokes `find_hairpin.py` via subprocess (5-minute timeout)
   - Results stored in `hairpin.db` (separate from ATLAS.db)

**find_hairpin.py** - Subgraph isomorphism engine:
- Loads `batch_0000_graphs.pickle` (pre-computed RNA structure graphs)
- Loads `temp_target_graphs.pickle` (user-drawn target)
- Uses NetworkX GraphMatcher with edge attribute matching
- Extracts PDB atom data from `data/rna_chain_corrected_2/*.pdb`
- Saves matches to `hairpin.db` with full atom content

### Database Schema

**ATLAS.db** (8.3 GB, pre-built, read-only):
- `data` table: id, motif_type, pdbid, nt_number, filecontent
- `PK` table: id, motif_type, pdbid, nt_number, file_content (note different column name)

**hairpin.db** (runtime-generated, ephemeral):
- `files` table: id, motif_type, pdbid, paired_nt_number, nt_number, filecontent, created_at

### File Layout

The application expects the entire deployment folder structure to remain intact:
- `data/rna_chain_corrected_2/` - Source PDB files for motif extraction
- `templates/` - Jinja2 HTML templates
- `static/` - CSS/JS/images
- `batch_0000_graphs.pickle` - Pre-computed graph database (~433 KB)
- `ATLAS.db` - Primary SQLite database

Runtime outputs (require write permissions):
- `hairpin.db` - Custom search results database
- `custom_pdb_files/` - Matched PDB files from custom searches
- `pdb_files/` - Matched PDB files from standard searches
- `search_results.csv`, `custom_search_results.csv` - Downloadable result tables
- `temp_target_graphs.pickle` - Temporary user input graph

### Key Algorithms

**Node ID Formats**: The code handles multiple PDB chain/residue formats:
- `A1` - Chain A, residue 1
- `A-1` - Chain A, residue -1 (negative numbering)
- `'0'1` - Chain 0, residue 1
- `'0'-1` - Chain 0, residue -1

**Chain Detection**: `get_chains()` splits node lists into separate RNA chains when residue numbers are non-consecutive.

**Graph Matching**: Edge attributes distinguish bond types (Watson-Crick pairs, backbone connections, etc.). The matcher requires exact attribute equality for edges to match.

### Deployment Considerations

- **Timeouts**: Custom motif searches can take up to 60 minutes. Set proxy/WSGI timeouts to ≥ 3700 seconds.
- **Memory**: Recommend 8 GB RAM for graph operations on large structures.
- **Working Directory**: Must run from deployment root; all paths are relative.
- **File Permissions**: Service account needs write access to: `hairpin.db`, `*_pdb_files/`, `*.csv`, `*.pickle` (temp files).

### Templates and Routes

Key routes:
- `/` - Homepage (index.html)
- `/user-guide` - Documentation
- `/search` [POST] - Standard motif search
- `/custom-motif-search` [GET] - Drawing interface (custom_motif_draw.html)
- `/process-custom-motif` [POST] - Custom search processing
- `/download_csv`, `/download_zip` - Standard search downloads
- `/download_custom_csv`, `/download_custom_pdb_zip` - Custom search downloads
- `/download-database` - Full ATLAS.db download

### Critical Path Dependencies

When modifying custom search flow:
1. Web interface → JSON graph data
2. `create_networkx_graph()` → pickle to `temp_target_graphs.pickle`
3. `run_find_hairpin()` subprocess → `find_hairpin.py` reads pickle
4. `find_hairpin.py` → writes to `hairpin.db`
5. `read_custom_search_results()` → reads `hairpin.db`
6. `store_custom_results()` → generates CSV and `custom_pdb_files/`

### Node ID to PDB Atom Extraction

The `extract_node_info()` → `read_pdb_and_find_nt()` pipeline:
- Strips quotes from node IDs (`'0'1` → `01`)
- Matches PDB ATOM lines using chain ID (column 22) + residue number (columns 23-27)
- Returns all atoms for the matching residue

### Important Constraints

- **Do not rename/move files**: The application uses relative paths throughout.
- **Database column inconsistency**: `data.filecontent` vs `PK.file_content` (underscore).
- **Graph pickle format**: Must be dict of `{graph_id: networkx.Graph}` with edge `attribute` keys.
- **PDB filename convention**: Expects `{pdbid}.pdb` in `data/rna_chain_corrected_2/`.
- **Maximum results**: Standard searches limited to 5000 results to prevent memory issues.
- **Custom search concurrency**: Only one custom motif search can run at a time (thread lock).

## Memory Optimization (December 2024)

### Problem Background

The original implementation had a critical memory issue on 4GB systems:

**Two-Stage Filtering Bug** (app.py:search_motif):
```python
# OLD CODE (BROKEN):
cursor.execute(query, (motif_type,))
rows = cursor.fetchall()  # Loads ALL 224,796 internal records = 7.7GB!
filtered_rows = [row for row in rows if len(row[3].split(',')) == int(nt_number)]
```

This caused:
- Internal searches: 224,796 records × 17.5KB = **7.7GB memory**
- Hairpin searches: 133,555 records × 17.5KB = **4.6GB memory**
- Result: **OOM (Out of Memory) crashes** on 4GB systems

### Implemented Fixes

#### 1. Cursor Iteration (P0 - Critical)
**File**: app.py:42-112 (`search_motif` function)

**Change**: Use cursor iteration instead of `fetchall()`
```python
# NEW CODE (FIXED):
for row in cursor:  # Iterate row-by-row
    if len(row[3].split(',')) == int(nt_number):
        filtered_rows.append(row)
        if len(filtered_rows) >= max_results:
            break  # Stop early
```

**Impact**: Memory reduced from 7.7GB → **<100MB** per search

#### 2. Result Count Limiting (P1)
**Files**: app.py:42, templates/results.html:13-17

**Change**: Hard limit of 5000 results per search
- Metadata tracking: `{'truncated': bool, 'total_scanned': int, 'result_count': int}`
- Warning banner displayed when results are truncated
- Prevents excessive memory usage even after fix #1

**Impact**: Maximum memory per search capped at ~150MB

#### 3. Custom Search Concurrency Control (P1)
**Files**: app.py:20, 36, 280-288

**Change**: Thread lock to prevent simultaneous custom searches
```python
custom_search_lock = Lock()  # Global lock

@app.route('/process-custom-motif')
def process_custom_motif():
    if not custom_search_lock.acquire(blocking=False):
        return error("Another search in progress")
    try:
        # ... search logic ...
    finally:
        custom_search_lock.release()
```

**Impact**: Prevents 2× memory usage from concurrent searches (4-6GB → stays under 2GB)

#### 4. ZIP Generation Optimization (P2)
**Files**: app.py:225-256 (`download_zip`), app.py:517-550 (`download_custom_pdb_zip`)

**Change**: Use temporary files instead of in-memory BytesIO buffers
```python
# OLD: zip_buffer = BytesIO()  # Entire ZIP in RAM
# NEW: temp_zip = tempfile.NamedTemporaryFile(suffix='.zip')
```

**Impact**: Memory reduction from ~500MB → ~150MB for large result sets

### Memory Budget (4GB System)

| Component | Before Fix | After Fix | Status |
|-----------|------------|-----------|--------|
| OS + Base | 1.5GB | 1.5GB | - |
| Python Runtime | 500MB | 500MB | - |
| **Single Internal Search** | 7.7GB | <100MB | ✅ **Fixed** |
| Result Processing | 500MB | 150MB | ✅ Optimized |
| ZIP Generation | 500MB | 150MB | ✅ Optimized |
| Custom Search | 2-3GB | 2-3GB | ⚠️ Still high |
| **Total Available** | ❌ Insufficient | ✅ ~2GB headroom | **Success** |

### Performance Benchmarks

On 4GB RAM systems:

| Operation | Before | After | Improvement |
|-----------|--------|-------|-------------|
| Internal 6nt search | ❌ Crash (OOM) | ✅ <30sec | **From failure to success** |
| Hairpin 10nt search | ⚠️ Unstable | ✅ <20sec | Reliable |
| 5000-result ZIP | 500MB peak | 150MB peak | 70% reduction |
| Concurrent custom | ❌ Crash | ✅ Queue (lock) | No conflicts |

### Testing

See `TESTING_GUIDE.md` for comprehensive testing checklist.

Quick validation:
```bash
python app.py
# Test: Internal, 6 nucleotides (should complete in <30s, no crash)
# Test: Download ZIP (should succeed without memory spike)
```

### Known Limitations

1. **Custom motif search** still memory-intensive (subgraph isomorphism is NP-complete)
   - Recommendation: Upgrade to 8GB RAM if heavy custom search usage expected
   - Mitigated by: Concurrency lock (prevents multiple simultaneous searches)

2. **Result truncation** at 5000 records
   - Adjustable via `max_results` parameter in `search_motif()`
   - Increasing beyond 10,000 not recommended on 4GB systems

## UX Improvements (November 2024)

A comprehensive set of user experience enhancements were implemented to improve navigation, provide better user feedback, and enhance result readability.

### Implemented Features

#### 1. Home Button Navigation
**Files Modified**:
- `templates/results.html`
- `templates/custom_results.html`
- `templates/custom_motif_draw.html`
- `templates/user_guide.html`
- `templates/no_results.html`

**Implementation**:
```html
<a href="/" class="home-button">🏠 Home</a>

<style>
  .home-button {
    position: fixed;
    top: 20px;
    left: 20px;
    background-color: #478ac9;
    color: white;
    z-index: 1000;
    box-shadow: 0 2px 8px rgba(0,0,0,0.2);
  }
</style>
```

**Purpose**: Provides consistent navigation back to homepage from all pages, fixed at top-left corner for easy access.

#### 2. Search Processing Indicators

**Type-Based Search (index.html)**:
- Blue loading spinner with "Searching database..." message
- Displays during standard motif searches
- Form integration: `onsubmit="showLoadingSpinner(); return convertNucleotideNumber();"`

**Graph-Based Search (custom_motif_draw.html)**:
- Red loading spinner with detailed progress message
- Warns users: "This may take up to 5 minutes"
- Prevents window closure during search
- Explains subgraph isomorphism algorithm processing

**Critical Bug Fix**: Form action changed from relative `action="search"` to absolute `action="/search"` to ensure proper navigation after search completion.

#### 3. Frontend Pagination System

**Files Modified**:
- `templates/results.html`
- `templates/custom_results.html`

**Implementation**:
- 100 results per page
- JavaScript-based pagination (no backend changes)
- Controls: First | Previous | Page X of Y | Next | Last
- Dynamic page range display: "Showing X-Y of Z results"
- Button state management (disabled on first/last page)
- Smooth scroll to top on page navigation

**Code Structure**:
```javascript
const rowsPerPage = 100;
function displayPage(page) {
  const start = (page - 1) * rowsPerPage;
  const end = start + rowsPerPage;
  rows.forEach((row, index) => {
    row.style.display = (index >= start && index < end) ? '' : 'none';
  });
  // Update pagination controls and scroll to top
}
```

**Impact**: Dramatically improved readability for large result sets (previously displayed all results on single page).

#### 4. User Guide Integration

**File Modified**: `templates/index.html`

**Changes**:
- User Guide condensed and integrated into homepage as section
- Navigation changed from `/user-guide` route to `#user-guide` anchor link
- Two-part structure:
  - Quick Start guide on homepage (simplified)
  - Full detailed guide accessible via link to `/user-guide`
- Content simplified: Removed technical one-hot encoding notation
- Clarified UI flow: "Option 1/2" instead of "Step 2/3" for drawing vs loading examples

#### 5. Contact Information Section

**File Modified**: `templates/index.html`

**Implementation**:
```html
<section class="u-clearfix u-section-contact" id="contact">
  <div style="background-color: #e8f4f8; padding: 25px; text-align: center;">
    <p><strong>Nikolay V. Dokholyan, PhD</strong></p>
    <p>Professor, Neurology</p>
    <p>UVA School of Medicine</p>
    <p><strong>Email:</strong> <a href="mailto:dokh@virginia.edu">dokh@virginia.edu</a></p>
  </div>
</section>
```

**Purpose**: Provides clear contact information for user inquiries and feedback.

**Design Decision**: Removed mailto feedback form (poor UX in many email clients) in favor of direct email display.

#### 6. Typography Adjustments

**File Modified**: `templates/index.html` (line 35)

**Homepage Title Optimization**:
```html
<h1 class="u-title" style="font-size: 3.78rem !important;">
  RNA ATLAS: 3D RNA Motif Library
</h1>
```

**Changes**:
- Reduced from default 4.5rem to 3.78rem (84% of original size)
- Removed line break to ensure single-line display
- Used absolute rem units with `!important` flag to override CSS cascade
- Maintains responsive design while improving visual hierarchy

### Technical Implementation Notes

**CSS Specificity**: Required `!important` flag to override external `nicepage.css` stylesheet (h1.u-title default: 4.5rem).

**Flask Debug Mode**: Running with `debug=False` requires manual Flask restart for template changes to take effect.

**Browser Caching**: Users may need hard refresh (Ctrl+Shift+R) to see CSS/HTML changes immediately.

**Form Submission**: Absolute paths (`action="/search"`) used instead of relative paths to prevent navigation issues after loading spinner integration.

### Files Summary

**Modified Templates** (6 files):
1. `templates/index.html` - Home buttons, loading spinner, user guide integration, contact section, title font size
2. `templates/results.html` - Home button, pagination system
3. `templates/custom_results.html` - Home button, pagination system
4. `templates/custom_motif_draw.html` - Home button, graph search loading spinner
5. `templates/user_guide.html` - Home button
6. `templates/no_results.html` - Home button

**No Backend Changes**: All improvements are frontend-only (HTML/CSS/JavaScript).

### Testing Checklist

After UX improvements:
- ✅ Home buttons visible and functional on all pages
- ✅ Type-based search shows loading spinner until results load
- ✅ Graph-based search shows loading spinner with 5-minute warning
- ✅ Results pagination works (100 items/page, all controls functional)
- ✅ User Guide accessible from homepage anchor link
- ✅ Contact information displays correctly
- ✅ Title displays on single line at correct font size
- ✅ Form submissions navigate to results page correctly

### Known Issues Resolved

**Issue 1**: Search function failure after loading spinner integration
- **Root Cause**: Relative path `action="search"` instead of `action="/search"`
- **Solution**: Changed to absolute path in form action attribute

**Issue 2**: Database corruption (all searches returning no results)
- **Root Cause**: ATLAS.db file moved/corrupted (0 bytes)
- **Solution**: Restored original 8.4GB database file
- **Prevention**: Document importance of maintaining file structure

**Issue 3**: Title font size not changing
- **Root Causes**:
  - Flask not reloading templates (debug=False)
  - `<br>` tag forcing line break
  - Relative % sizing instead of absolute rem
  - CSS specificity (needed !important)
- **Solution**: Absolute rem units + !important + remove br tag + Flask restart

### Temporary Files Cleanup

**Runtime Generated Directories**:
- `pdb_files/` - Type-based search temporary PDB files (safe to delete)
- `custom_pdb_files/` - Custom search temporary PDB files (safe to delete)
- These directories are regenerated on each download request
- Original data stored in ATLAS.db (filecontent column)

**Cleanup Strategy**:
```bash
# Manual cleanup
rm -rf pdb_files/*
rm -rf custom_pdb_files/*

# Or add to app.py startup/shutdown
import atexit
def cleanup_temp_files():
    shutil.rmtree('pdb_files', ignore_errors=True)
    os.makedirs('pdb_files', exist_ok=True)
atexit.register(cleanup_temp_files)
```
