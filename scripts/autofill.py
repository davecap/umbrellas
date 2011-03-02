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
    if len(args) < 2:
        parser.error('Please provide spacing')
    spacing = float(args[1])
    ensemble = umbrellas.Ensemble(args[0])
    
    done = False
    while not done:
        done = True
        # sort replicas by coordinate
        sorted_replicas = sorted(ensemble.replicas.values(), key=lambda r: float(r.coordinate()), reverse=False)
        prev = sorted_replicas[0]
        for i,r in enumerate(sorted_replicas[1:]):
            diff = abs(r.coordinate()-prev.coordinate())
            # if the difference is less than spacing
            if diff > spacing:
                done = False
                # mutate prev towards r by spacing or half the distance, whatever is smaller
                step = min(diff/2.0, spacing)
                step = -step
                print " Replicas %s (%f) and %s (%f) differ by %f, will create new replica" % (r.name, r.coordinate(), prev.name, prev.coordinate(), diff)
                print " Mutating %s (%f) by %f" % (prev.name, prev.coordinate(), step)
                new_replica = ensemble.copy_replica(prev.name)
                new_replica.mutate(step) # mutate and calc new coordinate
                new_replica.save() # save mutation
                print " New replica %s at %s (%f)" % (new_replica.name, new_replica.coordinates(), new_replica.coordinate())
            prev = r
        
        print "Reverse!"
        
        # sort replicas by coordinate REVERSE
        sorted_replicas = sorted(ensemble.replicas.values(), key=lambda r: float(r.coordinate()), reverse=True)
        prev = sorted_replicas[0]
        for i,r in enumerate(sorted_replicas[1:]):
            diff = abs(r.coordinate()-prev.coordinate())
            # if the difference is less than spacing
            if diff > spacing:
                done = False
                # mutate prev towards r by spacing or half the distance, whatever is smaller
                step = min(diff/2.0, spacing)
                print " Replicas %s (%f) and %s (%f) differ by %f, will create new replica" % (r.name, r.coordinate(), prev.name, prev.coordinate(), diff)
                print " Mutating %s (%f) by %f" % (prev.name, prev.coordinate(), step)
                new_replica = ensemble.copy_replica(prev.name)
                new_replica.mutate(step) # mutate and calc new coordinate
                new_replica.save() # save mutation
                print " New replica %s at %s (%f)" % (new_replica.name, new_replica.coordinates(), new_replica.coordinate())
            prev = r
    ensemble.save()
            
if __name__ == "__main__":
    main()
