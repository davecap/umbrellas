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
from reaction import *

class Ensemble:
    """ The Ensemble class contains and manages Replicas """
    
    def __init__(self, config_path='config.ini'):
        # Load the config file
        if not os.path.exists(config_path):
            logger.warning('Config file %s not found so default will be created' % config_path)
        self.config = setup_config(config_path, create=True)
        
        if self.config['debug']:
            logger.setLevel(logging.DEBUG)
            logger.debug('Debug mode is on!')
        else:
            logger.setLevel(logging.WARNING)
            
        # Load the replica DB file (path from config)
        self.replicadb = setup_replicadb(self._replicadb_path(), create=True)
        self.load()
        
        # Instantiate the Reaction class and validate it against the config
        # TODO: improve this!
        self.reaction = globals()[self.config['reaction']['type']](self.config)
    
    def _replicadb_path(self):
        return os.path.join(self.basedir(), self.config['replicadb'])
    
    def basedir(self):
        """ Returns the directory containing the current config file. """
        return os.path.dirname(self.config.filename)
    
    def add_replica(self, name, **kwargs):
        """ Add a new replica (name must be unique) with the kwargs as parameters. """
        if name in self.replicas:
            raise Exception('Cannot create replica with name %s, already taken!' % name)
        new_replica = Replica(self, name, **kwargs)
        self.replicas[name] = new_replica
        # TODO: autosave?
        return new_replica
        
    def get_replica(self, name):
        """ Get a replica by name, raises exception if not found. """
        try:
            return self.replicas[name]
        except KeyError:
            # TODO: subclassed exceptions
            raise Exception('Replica with name %s not found!' % name)
    
    def load(self):
        """ Load the replica DB from the file. Overwrites the current state. """
        self.replicas = {}
        for name in self.replicadb['replicas']:
            self.replicas[name] = Replica(self, name, **self.replicadb['replicas'][name])
                
    def save(self):
        """ Save the replica DB. """
        for name,r in self.replicas.items():
            self.replicadb['replicas'][name] = r.parameters
        self.replicadb.write()
        self.config.write()

class Replica:
    """ The Replica class defines a single replica in a system. """
    
    # reserved parameters
    # TODO: should the topology file be shared from the Ensemble?
    TOPOLOGY_PATH = 'topology' # topology file path
    # TODO: are these too similar? should one be renamed?
    COORDINATES_PATH = 'coordinates' # coordinates file path
    COORDINATE = 'coordinate' # reaction coordinate
    FORCE = 'force'
    
    DEFAULTS = { FORCE: 0.0, COORDINATES_PATH: '' }
    
    def __init__(self, ensemble, name, **kwargs):
        self.ensemble = ensemble # Ensemble object reference
        self.name = name # name of this Replica
        self._universe = None # local MDAnalysis Universe object, see self.universe()
        self.parameters = Replica.DEFAULTS.copy() # copy a set of defaults first
        self.parameters.update(kwargs) # update the defaults with the passed kwargs
        
        if not self.parameter(Replica.COORDINATES_PATH):
            raise Exception('No coordinate path for replica %s' % self.name)
            
        if not self.parameter(Replica.COORDINATE):
            logging.warning('No coordinate for replica %s, calculating now...' % self.name)
            self.coordinate()
        
    def topology(self):
        """ OPTIONAL Path to the topology file (PSF). """
        return self.parameter(Replica.TOPOLOGY_PATH)
        
    def coordinates(self):
        """ REQUIRED Path to the coordinates file (PDB, DCD, CRD). """
        return self.parameter(Replica.COORDINATES_PATH)
    
    def u_topology(self):
        """ Structure/Topology file from the MDAnalysis Universe perspective """
        if not self._universe or not self.topology():
            return None
        else:
            return self._universe.filename
    
    def u_coordinates(self):
        """ Coordinate file from the MDAnalysis Universe perspective """
        if not self._universe:
            return None
        if self.topology():
            return self._universe.trajectory.filename
        else:
            return self._universe.filename
    
    def parameter(self, k):
        """ Returns a replica parameter, or None """
        if k in self.parameters:
            return self.parameters[k]
        else:
            return None
        
    def universe(self, force=False):
        """ Generate/Return a MDanalysis Universe object"""
        # check to make sure we haven't already loaded the universe
        if not force and self._universe and self.u_topology() == self.topology() and self.u_coordinates() == self.coordinates():
            return self._universe
        
        if self.coordinates() is None or not os.path.exists(self.coordinates()):
            raise Exception('Coordinates path for replica %s not found: %s' % (self.name, self.coordinates()))
        
        if self.topology() and not os.path.exists(self.topology()):
            raise Exception('Topology file for replica %s not found: %s' % (self.name, self.topology()))
        
        if self.topology() and self.coordinates():
            # first see if we have both topology and coordinate files
            self._universe = MDAnalysis.Universe(self.topology(), self.coordinates())
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
        self.parameters[Replica.COORDINATES_PATH] = path
     
    def coordinate(self, force=False):
        """ Return the coordinate for this replica, calculated automatically from the coordinates"""
        if force or not self.parameter(Replica.COORDINATE):
            # save this value to the replicas.db
            self.parameters[Replica.COORDINATE] = self.ensemble.reaction.coordinate(self.universe(force))
        return self.parameter(Replica.COORDINATE)
    
    def mutate(self, step=1.0):
        """ Mutate this replica's coordinate by step units.
            The proper step size depends on the current reaction type in use.
        """
        self.ensemble.reaction.mutate(self.universe(), step)
        # reset the coordinate after a mutation
        self.parameters[Replica.COORDINATE] = None
