import os
import logging
import time
import MDAnalysis

logger = logging.getLogger('umbrellas')
_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
_handler = logging.StreamHandler()
_handler.setFormatter(_formatter)
logger.addHandler(_handler)

from config import setup_config, setup_replicadb
from reaction import Distance, Angle, Dihedral

class Ensemble:
    """ The Ensemble class contains and manages Replicas """
     
    def __init__(self, config_path='config.ini'):
        # Load the config file
        if not os.path.exists(config_path):
            logger.warning('Config file %s not found so default will be created' % config_path)
        self.config = setup_config(config_path, create=True)
        
        # Load the replica DB file (path from config)
        self.replicadb = setup_replicadb(self._replicadb(), create=True)
        
        if self.config['debug']:
            logger.setLevel(logging.DEBUG)
            logger.debug('Debug mode is on!')
        else:
            logger.setLevel(logging.WARNING)
        
        # Instantiate the Reaction class and validate it against the config
        # TODO: improve this!
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
        replicas = self.replicas()
        if name in replicas:
            logging.error('Cannot create replica with name %s, already taken!' % name)
            return None
        r = Replica(self, name, **kwargs)
        replicas[name] = r.parameters
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
    COORDINATE_KEY = 'coordinate'
    FORCE_KEY = 'force'
    
    DEFAULTS = { FORCE_KEY: 0.0, COORDINATE_KEY: '', COORDINATES_PATH_KEY: '' }
    
    def __init__(self, ensemble, name, **kwargs):
        self.ensemble = ensemble # Ensemble object reference
        self.name = name # name of this Replica
        self._universe = None # local MDAnalysis Universe object, see self.universe()
        self.parameters = Replica.DEFAULTS.copy() # copy a set of defaults first
        
        if not self.parameter(Replica.COORDINATES_PATH_KEY):
            logging.warning('No coordinate path for replica %s' % self.name)
        
        self.parameters.update(kwargs)
        
        
    def structure(self):
        """ OPTIONAL Path to the structure file (PSF). """
        return self.parameter(Replica.STRUCTURE_PATH_KEY)
        
    def coordinates(self):
        """ REQUIRED Path to the coordinates file (PDB, DCD, CRD). """
        return self.parameter(Replica.COORDINATES_PATH_KEY)
    
    def u_structure(self):
        """ Structure/Topology file from the MDAnalysis Universe perspective """
        if not self._universe or not self.structure():
            return None
        else:
            return self._universe.filename
    
    def u_coordinates(self):
        """ Coordinate file from the MDAnalysis Universe perspective """
        if not self._universe:
            return None
        if self.structure():
            return self._universe.trajectory.filename
        else:
            return self._universe.filename
    
    def parameter(self, k):
        """ Returns a replica parameter, or None """
        if k in self.parameters:
            return self.parameters[k]
        else:
            return None
        
    def universe(self):
        """ Generate/Return a MDanalysis Universe object"""
        # check to make sure we haven't already loaded the universe
        if self._universe and self.u_structure() == self.structure() and self.u_coordinates() == self.coordinates():
            return self._universe
        
        if self.coordinates() is None or not os.path.exists(self.coordinates()):
            raise Exception('Coordinates path for replica %s not found: %s' % (self.name, self.coordinates()))
        
        if self.structure() and not os.path.exists(self.structure()):
            raise Exception('Structure file for replica %s not found: %s' % (self.name, self.structure()))
        
        if self.structure() and self.coordinates():
            # first see if we have both structure and coordinate files
            self._universe = MDAnalysis.Universe(self.structure(), self.coordinates())
        else:
            # fallback to coordinate only
            self._universe = MDAnalysis.Universe(self.coordinates())
        return self._universe
        
    def save(self, path=None, overwrite=False):
        """ Save the coordinates to a new file (and save the new file name).
            This function will overwrite the current coordinates file if overwite is True """
        # if no path is specified, use the current path
        if path is None:
            # just overwrite the current path
            path = self.coordinates()
        
        # write to a new file in the same dir as path, dont overwrite path
        if os.path.exists(path) and not overwrite:
            basename, ext = os.path.splitext(path)
            # get time in seconds as a string
            t = str(time.time()).split('.')[0]
            path = os.path.join(os.path.dirname(path), basename+'_'+t+ext)
            if os.path.exists(path):
                # this shouldn't really ever happen...
                # TODO: figure out how to generate a filename that doesn't exist
                raise Exception('Generated file name already exists (%s). Will not overwrite so try again.' % path)
            
        # write the universe
        w = MDAnalysis.Writer(path)
        w.write(self.universe())
        
        # set the new file as this coordinates file
        self.parameters[Replica.COORDINATES_PATH_KEY] = path
        # reset the coordinate
        
        # save the ensemble
        self.ensemble.save()
     
    def coordinate(self):
        """ Return the coordinate for this replica, calculated automatically from the coordinates"""
        if not self.parameter(Replica.COORDINATE_KEY):
            # save this value to the replicas.db
            self.parameters[Replica.COORDINATE_KEY] = self.ensemble.reaction.coordinate(self.universe())
            self.ensemble.save()
        return self.parameter(Replica.COORDINATE_KEY)
    
    def mutate(self, step=1.0):
        """ Mutate this replica's coordinate by step units.
            The proper step size depends on the current reaction type in use.
        """
        self.ensemble.reaction.mutate(self.universe(), step)
        # reset the coordinate after a mutation
        self.parameters[Replica.COORDINATE_KEY] = None
