"""Gunicorn configuration for RNAdex.

Run from the deployment root (all application paths are relative to it):

    gunicorn -c gunicorn.conf.py app:app

See IT-Deployment-Notes.md for the full deployment requirements.
"""

import multiprocessing
import os

# --- socket ---------------------------------------------------------------
bind = os.environ.get("RNADEX_BIND", "127.0.0.1:8000")

# --- workers --------------------------------------------------------------
# Threaded workers so one long custom search does not block a whole process.
# Custom searches are serialised globally anyway (fcntl lock + threading lock
# in app.py), so extra workers only add capacity for ordinary requests.
worker_class = "gthread"
workers = int(os.environ.get("RNADEX_WORKERS", min(4, multiprocessing.cpu_count())))
threads = int(os.environ.get("RNADEX_THREADS", 2))

# A standard search holds up to 10,000 rows of PDB text in memory and renders
# a large page, so give each worker room rather than packing them tightly.
max_requests = 200
max_requests_jitter = 40

# --- timeouts -------------------------------------------------------------
# A custom motif search may run for the better part of an hour. Any reverse
# proxy in front of this must use a matching (or larger) timeout, otherwise the
# proxy will cut the connection while the search is still running.
timeout = 3700
graceful_timeout = 120
keepalive = 5

# --- startup --------------------------------------------------------------
# Import the app once in the master before forking. app.py runs
# cleanup_on_startup() at import, which rmtree's search_sessions/, pdb_files/
# and custom_pdb_files/; letting every worker do that concurrently makes them
# race and log "Directory not empty" failures.
preload_app = True

# --- logging --------------------------------------------------------------
accesslog = os.environ.get("RNADEX_ACCESS_LOG", "-")
errorlog = os.environ.get("RNADEX_ERROR_LOG", "-")
loglevel = os.environ.get("RNADEX_LOG_LEVEL", "info")
access_log_format = '%(h)s "%(r)s" %(s)s %(b)s %(M)sms'

# --- misc -----------------------------------------------------------------
proc_name = "rnadex"
# Large result pages and long nucleotide lists make for long request lines.
limit_request_line = 8190
forwarded_allow_ips = os.environ.get("RNADEX_FORWARDED_ALLOW_IPS", "127.0.0.1")
