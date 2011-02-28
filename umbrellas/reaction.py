class BaseReaction:
    """ The BaseReaction defines the type of reaction for the given system.
        A reaction can be a change of distance, angle, dihedral or any other property.
    """
    def __init__(self, config):
        self.config = config['reaction']
        
    def validate(self):
        raise NotImplementedError('Subclassess of BaseReaction must define a validation method.')
        
    def coordinate(self, universe):
        raise NotImplementedError('Subclassess of BaseReaction must define a coordinate calculation.')
            
class Distance(BaseReaction):
    """ Distance reactions require 2 parameters:
        reference and components, both strings
    """
    def validate(self):
        if not self.config['target']:
            raise Exception('No target specified in reaction config')
        if not self.config['reference']:
            raise Exception('No reference specified in reaction config')
        if not ('x' in self.config['components'] or 'y' in self.config['components'] or 'z' in self.config['components']):
            raise Exception('Invalid distance components specified in reaction config, must be a combination of xyz')
        
    def coordinate(self, universe):
        """ Somehow, the MDAnalysis Universe must be passed here so we can extract the coordinate """
        target = universe.selectAtoms(self.config['target']).centerOfMass()
        reference = universe.selectAtoms(self.config['reference']).centerOfMass()
        
        t = [0.0,0.0,0.0]
        r = [0.0,0.0,0.0]
        total = 0
        
        if 'x' in self.config['components']:
            t[0] = target[0]
            r[0] = reference[0]
            total += 1
        if 'y' in self.config['components']:
            t[1] = target[1]
            r[1] = reference[1]
            total += 1
        if 'z' in self.config['components']:
            t[2] = target[2]
            r[2] = reference[2]
            total += 1
            
        distance = ((t[0]-r[0])**2 + (t[1]-r[1])**2 + (t[2]-r[2])**2)**(0.5)
        
        # if we are in 1D and the target is < than the reference, make the distance negative
        # TODO: make this optional
        if total == 1 and t < r:
            distance = distance * -1
        
        return distance
        
    
class Angle(BaseReaction):
    def validate(self):
        return True
        
    def coordinate(self, universe):
        return 0.0

class Dihedral(BaseReaction):
    def validate(self):
        return True
        
    def coordinate(self, universe):
        return 0.0
