# Umbrellas

import os
import logging
logger = logging.getLogger('umbrellas')

_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
_handler = logging.StreamHandler()
_handler.setFormatter(_formatter)
logger.addHandler(_handler)

import MDAnalysis

from config import setup_config, setup_replicadb
from reaction import Distance, Angle, Dihedral

class Ensemble:
    """ The Ensemble class contains and manages Replicas """
     
    def __init__(self, config_path='config.ini', replicadb_path='replicas.db'):
        # Load the config file
        if not os.path.exists(config_path):
            logger.warning('Config file %s not found so default will be created' % config)
        self.config = setup_config(config_path, create=True)
        
        # Load the replica DB file
        self.replicadb = setup_replicadb(replicadb_path, create=True)
        
        if self.config['debug']:
            logger.setLevel(logging.DEBUG)
            logger.debug('Debug mode is on!')
        else:
            logger.setLevel(logging.WARNING)
        
        # Instantiate the Reaction class and validate it against the config
        self.reaction = globals()[self.config['reaction']['type']](self.config)
        
    def replicas(self):
        """ Load the replica DB """
        return self.replicadb['replicas']
        
    def save(self):
        """ Save the replica DB """
        self.replicadb.write()

class Replica:
    """ The Replica class defines a single replica in a system.
    """
    
    def __init__(self, ensemble, config):
        self.ensemble = ensemble
        self.config = config
    
    def load(self):
        """ Generate a MDanalysis Universe object"""
        
        pass
        
    def save(self, filename):
        """ Save the coordinates, optionally run through VMD"""
        pass
    
    def properties(self, properties={}):
        """ Properties are a key/value set of extra properties that a replica may have.
            They are simply written to the config file on save.
        """
     
    def coordinate(self):
        """ Return the coordinate for this replica, calculated automatically from the coordinates"""
        pass
    
    def generate(self, step=1.0):
        """ Generate another replica from self with the given step size.
            The proper step size depends on the current reaction type in use.
        """
        pass

