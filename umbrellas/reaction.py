class BaseReaction:
    """ The BaseReaction defines the type of reaction for the given system.
        A reaction can be a change of distance, angle, dihedral or any other property.
    """
    def __init__(self, config):
        pass
        
    def is_valid(self):
        raise NotImplementedError('Subclassess of BaseReaction must define a validation method.')
        
    def coordinate(self):
        raise NotImplementedError('Subclassess of BaseReaction must define a coordinate calculation.')
            
class Distance(BaseReaction):
    pass
    
class Angle(BaseReaction):
    pass

class Dihedral(BaseReaction):
    pass