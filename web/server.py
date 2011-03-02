import optparse
import sys
import os
import json

from flask import Flask, g, request, redirect, render_template, render_template_string, url_for, send_file, abort, jsonify, json
app = Flask(__name__)

import umbrellas
_config = None

def main():
    global _config
    
    usage = """
        usage: %prog [options] <config.ini>
    """
    parser = optparse.OptionParser(usage)
    # parser.add_option("-x", dest="x_column", default=None, help="X Column REQUIRED")    
    # parser.add_option("-y", dest="y_column", default=None, help="Y Column REQUIRED")
    options, args = parser.parse_args()
    
    if not args or not os.path.exists(args[0]):
        parser.error('No config.ini file found!')
    
    _config = args[0]
    app.run(debug=True)

@app.before_request
def before_request():
    global _config
    g.ensemble = umbrellas.Ensemble(_config)

@app.route('/')
def index():
    replicas = [ r.export() for r in g.ensemble.get_replicas() ]
    return render_template('index.html', config=g.ensemble.config, replicas=replicas)

@app.route('/replicas/<name>', methods=['GET','POST'])
def replicas(name):
    if request.method == 'POST':
        print request.form
        return jsonify(replica={})
    else:
        # GET
        return jsonify(g.ensemble.get_replica(name).export())
        
@app.route('/replicas/<name>/parameters', methods=['GET','POST'])
def replica_parameters(name):
    if request.method == 'POST':
        print request.form
        return jsonify(replica={})
    else:
        # GET
        return jsonify(g.ensemble.get_replica(name).export()['parameters'])

if __name__ == "__main__":
    main()
