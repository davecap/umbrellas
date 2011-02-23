import optparse
from MDAnalysis.tests.datafiles import PDB_small

import umbrellas

def main():
    usage = """
        usage: %prog [options]
    """
    
    parser = optparse.OptionParser(usage)
    #parser.add_option("-c", "--config", dest="config_file", default="config.ini", help="Config file [default: %default]")
    #parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False, help="Log to the console [default: %default]")    
    (options, args) = parser.parse_args()
    
    e = umbrellas.Ensemble()
    print e.replicas()
    print e.reaction
    
    # add a few replicas with the PDB small
    
    
if __name__=='__main__':
    main()
