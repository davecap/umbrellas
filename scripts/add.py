"""
    Add one or more PDB files to the replicas DB
"""

import optparse
import sys
import os

import umbrellas

def main():    
    usage = """
        usage: %prog [options] <config.ini> <PDB>
    """
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    
    if not args or not os.path.exists(args[0]):
        parser.error('No config.ini file found!')
        
    ensemble = umbrellas.Ensemble(args[0])
    for a in args[1:]:
        if not os.path.exists(a):
            print 'PDB file not found at %s' % a
        else:
            path, ext = os.path.splitext(a)
            name = os.path.basename(path)
            try:
                ensemble.add_replica(name=name, coordinates=a, force=2.0)
            except Exception, e:
                print str(e)
            else:
                ensemble.save()

if __name__ == "__main__":
    main()
