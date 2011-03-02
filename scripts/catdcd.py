import optparse
import sys
import os

import umbrellas

def main():    
    usage = """
        usage: %prog [options] <config.ini> <spacing>
    """
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    
    if not args or not os.path.exists(args[0]):
        parser.error('No config.ini file found!')
        
    ensemble = umbrellas.Ensemble(args[0])
    
    sorted_replicas = sorted(ensemble.replicas.values(), key=lambda r: float(r.coordinate()), reverse=False)
    
    ret = 'catdcd -o poop.dcd -otype dcd '
    for r in sorted_replicas:
        ret += '-pdb %s ' % r.coordinates()
    print ret
            
if __name__ == "__main__":
    main()
