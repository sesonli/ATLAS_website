"""
Configuration file for ATLAS RNA 3D Motif Library
Created: November 2024
"""

import os

class Config:
    """Base configuration"""

    # Flask Configuration
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'dev-secret-key-change-in-production'
    DEBUG = False

    # Database Configuration
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    ATLAS_DB_PATH = os.path.join(BASE_DIR, 'ATLAS.db')
    HAIRPIN_DB_PATH = os.path.join(BASE_DIR, 'hairpin.db')

    # File Paths
    BATCH_GRAPHS_PATH = os.path.join(BASE_DIR, 'batch_0000_graphs.pickle')
    TEMP_GRAPHS_PATH = os.path.join(BASE_DIR, 'temp_target_graphs.pickle')
    PDB_SOURCE_DIR = os.path.join(BASE_DIR, 'data', 'rna_chain_corrected_2')
    PDB_OUTPUT_DIR = os.path.join(BASE_DIR, 'pdb_files')
    CUSTOM_PDB_OUTPUT_DIR = os.path.join(BASE_DIR, 'custom_pdb_files')

    # Search Configuration
    MAX_SEARCH_RESULTS = 5000  # Maximum results per search
    CUSTOM_SEARCH_TIMEOUT = 300  # 5 minutes in seconds
    MIN_ATOMS_PER_RESIDUE = 10  # Minimum atoms for complete residue data

    # C4' Distance Thresholds (Angstroms)
    C4_MIN_DISTANCE = 4.0
    C4_MAX_DISTANCE = 8.0

    # Logging Configuration
    LOG_FILE = 'app.log'
    LOG_LEVEL = 'INFO'

    # Server Configuration
    HOST = '0.0.0.0'
    PORT = 5000
    THREADED = True

    # Memory Management
    ENABLE_CUSTOM_SEARCH_LOCK = True  # Prevent concurrent custom searches

class DevelopmentConfig(Config):
    """Development configuration"""
    DEBUG = True
    LOG_LEVEL = 'DEBUG'

class ProductionConfig(Config):
    """Production configuration"""
    DEBUG = False
    LOG_LEVEL = 'WARNING'

    # Production timeouts (for gunicorn/uwsgi)
    WORKER_TIMEOUT = 3700  # > 1 hour for long searches

    def __init__(self):
        """Initialize production config with environment variable validation"""
        super().__init__()
        # Production should use environment variables
        self.SECRET_KEY = os.environ.get('SECRET_KEY')
        if not self.SECRET_KEY:
            raise ValueError("SECRET_KEY environment variable must be set in production")

class TestingConfig(Config):
    """Testing configuration"""
    TESTING = True
    MAX_SEARCH_RESULTS = 100  # Smaller for testing
    CUSTOM_SEARCH_TIMEOUT = 60  # 1 minute for tests

# Configuration dictionary
config = {
    'development': DevelopmentConfig,
    'production': ProductionConfig,
    'testing': TestingConfig,
    'default': DevelopmentConfig
}

def get_config(config_name=None):
    """Get configuration object based on environment"""
    if config_name is None:
        config_name = os.environ.get('FLASK_ENV', 'default')
    return config[config_name]