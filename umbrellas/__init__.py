# Umbrellas

__version__ = "1.0"

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
     
    def __init__(self, config_path='config.ini'):
        # Load the config file
        if not os.path.exists(config_path):
            logger.warning('Config file %s not found so default will be created' % config)
        self.config = setup_config(config_path, create=True)
        
        # Load the replica DB file (path from config)
        self.replicadb = setup_replicadb(self._replicadb(), create=True)
        
        if self.config['debug']:
            logger.setLevel(logging.DEBUG)
            logger.debug('Debug mode is on!')
        else:
            logger.setLevel(logging.WARNING)
        
        # Instantiate the Reaction class and validate it against the config
        self.reaction = globals()[self.config['reaction']['type']](self.config)
    
    def _basedir(self):
        return os.path.dirname(self.config.filename)
    
    def _replicadb(self):
        return os.path.join(self._basedir(), self.config['replicadb'])
    
    def replicas(self):
        """ Load the replica DB """
        return self.replicadb['replicas']
        
    def add_replica(self, name, **kwargs):
        # if not name:
        #     name = self._generate_name()
        r = Replica(self, name, **kwargs)
        replicas = self.replicas()
        
        if name in replicas:
            logging.error('Cannot create replica with name %s, already taken!' % name)
            return None
        
        replicas[name] = kwargs
        self.save()
        return r
        
    def get_replica(self, name):
        replicas = self.replicas()
        if name in replicas:
            return Replica(self, name, **replicas[name])
        else:
            logging.error('Replica with name %s not found!' % name)
            return None
        
    def save(self):
        """ Save the replica DB """
        self.replicadb.write()

class Replica:
    """ The Replica class defines a single replica in a system.
    """
    
    def __init__(self, ensemble, name, **kwargs):
        self.ensemble = ensemble
        self.name = name
        self.parameters = kwargs
    
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

