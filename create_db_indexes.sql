-- Database Index Recommendations for ATLAS.db
-- These indexes will improve query performance significantly
-- Run these commands in SQLite to create indexes

-- Index for motif_type queries (most common filter)
CREATE INDEX IF NOT EXISTS idx_data_motif_type ON data(motif_type);
CREATE INDEX IF NOT EXISTS idx_pk_motif_type ON PK(motif_type);

-- Composite index for motif_type and pdbid (frequently queried together)
CREATE INDEX IF NOT EXISTS idx_data_motif_pdbid ON data(motif_type, pdbid);
CREATE INDEX IF NOT EXISTS idx_pk_motif_pdbid ON PK(motif_type, pdbid);

-- Index for pdbid (used in filtering and joins)
CREATE INDEX IF NOT EXISTS idx_data_pdbid ON data(pdbid);
CREATE INDEX IF NOT EXISTS idx_pk_pdbid ON PK(pdbid);

-- For hairpin.db (custom search results)
-- Note: Run these on hairpin.db separately
-- CREATE INDEX IF NOT EXISTS idx_files_motif_type ON files(motif_type);
-- CREATE INDEX IF NOT EXISTS idx_files_pdbid ON files(pdbid);
-- CREATE INDEX IF NOT EXISTS idx_files_created_at ON files(created_at DESC);

-- To check existing indexes:
-- .indexes

-- To analyze database for optimization:
-- ANALYZE;

-- To check index usage:
-- EXPLAIN QUERY PLAN SELECT * FROM data WHERE motif_type = 'internal' LIMIT 100;