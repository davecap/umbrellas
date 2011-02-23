# Umbrellas

import os
import logging
logger = logging.getLogger('umbrellas')
f = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
h = logging.StreamHandler()
h.setFormatter(f)
logger.addHandler(h)

import MDAnalysis

from config import setup_config, setup_replicadb
from reaction import Distance, Angle, Dihedral

class Ensemble:
    """ The Ensemble class contains and manages Replicas """
     
    def __init__(self, config='config.ini', replicadb='replicas.db'):
        """ Load the config file"""
        if not os.path.exists(config):
            logger.warning('Config file %s not found so default will be created' % config)
        self.config = setup_config(config, create=True)
        
        if self.config['debug']:
            logger.setLevel(logging.DEBUG)
            logger.debug('Debug mode is on!')
        else:
            logger.setLevel(logging.WARNING)
        
        # Instantiate the Reaction class and validate it against the config
        self.reaction = globals()[self.config['reaction']['type']]()
        
        # Load the replica DB file, assert that it exists
        self.replicadb = setup_replicadb(replicadb, create=True)
        self.load()
        
        print self.replicas
        
    def reaction(self):
        """ Returns the current reaction type """
        pass
        
    def load(self):
        """ Load the replica DB """
        self.replicas = self.replicadb['replicas']
        
    def save(self):
        """ Save the replica DB """
        self.replicadb.write()

class Replica:
    """ The Replica class defines a single replica in a system.
    """
    
    def __init__(self, config):
        """ Load the config file"""
        pass
    
    def load(self, structure, coordinates):
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

