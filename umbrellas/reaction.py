class BaseReaction:
    """ The BaseReaction defines the type of reaction for the given system.
        A reaction can be a change of distance, angle, dihedral or any other property.
    """
    def __init__(self, config):
        self.config = config['reaction']
        
    def is_valid(self):
        raise NotImplementedError('Subclassess of BaseReaction must define a validation method.')
        
    def coordinate(self):
        raise NotImplementedError('Subclassess of BaseReaction must define a coordinate calculation.')
            
class Distance(BaseReaction):
    def is_valid(self):
        return True
        
    def coordinate(self):
        return 0.0
    
class Angle(BaseReaction):
    def is_valid(self):
        return True
        
    def coordinate(self):
        return 0.0

class Dihedral(BaseReaction):
    def is_valid(self):
        return True
        
    def coordinate(self):
        return 0.0
