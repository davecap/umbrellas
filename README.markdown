Umbrellas
=========

#### What is Umbrellas?

Umbrellas is a small framework for setting up umbrella replicas in a molecular dynamics umbrella sampling simulation.
It is based on the MDAnalysis package and optionally VMD.

#### How do I use Umbrellas?

Umbrella supports simple US reaction coordinates "out of the box":

* Distance (any xyz component combination)
* Angle
* Dihedral

All that is required is a config file in the following format:
    
    title = <optional text title for this system>
    
    [reaction]
        type = [Distance, Angle, Dihedral] # Choose one. Not case sensitive
        target = <MDAnalysis selection: target residue>
        
        # For Distance
        reference = <MDAnalysis selection: reference residue> # Distance type only
        components = xyz
        
        # Angle and Dihedral
        atoms = <list of atoms> # Comma separated. Required for Angle and Dihedral types
        
#### How does it work?

The system will generate a replicas.db file, which is just a text config file containing the replica information.
You provide the structure file and Umbrella will help you do the rest.
