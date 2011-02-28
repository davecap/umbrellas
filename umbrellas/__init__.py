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
    STRUCTURE_PATH_KEY = 'structure'
    COORDINATES_PATH_KEY = 'coordinates'
    
    def __init__(self, ensemble, name, **kwargs):
        self.ensemble = ensemble
        self.name = name
        self.universe = None
        self._parameters = kwargs
        
    def structure(self):
        """ OPTIONAL Path to the structure file (PSF). """
        return self.parameter(Replica.STRUCTURE_PATH_KEY)
        
    def coordinates(self):
        """ REQUIRED Path to the coordinates file (PDB, DCD, CRD). """
        return self.parameter(Replica.COORDINATES_PATH_KEY)
    
    def u_structure(self):
        # self.universe.filename
        # OR
        # None if structure not set
        if not self.universe or not self.structure():
            return None
        else:
            return self.universe.filename
    
    def u_coordinates(self):
        # self.universe.filename IF STRUCTURE (topology) IS NOT SET
        # OR
        # self.universe.trajectory.filename IF STRUCTURE IS SET
        if not self.universe:
            return None
        if self.structure():
            return self.universe.trajectory.filename
        else:
            return self.universe.filename
    
    def parameter(self, k):
        if k in self._parameters:
            return self._parameters[k]
        else:
            return None
    
    def load(self):
        """ Generate a MDanalysis Universe object"""
        if not os.path.exists(self.coordinates()):
            raise Exception ('Coordinates path for replica %s not found: %s' % (self.name, self.coordinates()))
        
        if self.structure() and not os.path.exists(self.structure()):
            raise Exception ('Structure path for replica %s not found: %s' % (self.name, self.structure()))
        
        if self.structure() and self.coordinates():
            # first see if we have both structure and coordinate files
            self.universe = MDAnalysis.Universe(self.structure(), self.coordinates())
        else:
            # fallback to coordinate only
            self.universe = MDAnalysis.Universe(self.coordinates())
        return self.universe
        
    def save(self, overwrite=False):
        """ Save the coordinates to a new file (and save the new file name).
            This function will overwrite the current coordinates file if overwite is True """
        self.ensemble.save()
     
    def coordinate(self):
        """ Return the coordinate for this replica, calculated automatically from the coordinates"""
        self.load()
        pass
    
    def generate(self, step=1.0):
        """ Generate another replica from self with the given step size.
            The proper step size depends on the current reaction type in use.
        """
        pass

