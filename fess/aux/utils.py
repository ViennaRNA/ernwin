#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
import subprocess
import fess
import forgi
import forgi.graph.bulge_graph as fgb
from fess import __version__
__metaclass__=type

def get_all_subclasses(cls, include_base = False):
    """
    Thanks to fletom at http://stackoverflow.com/a/17246726/5069869
    """        
    all_subclasses = []
    if include_base:
        all_subclasses.append(cls)

    for subclass in cls.__subclasses__():
        all_subclasses.append(subclass)
        all_subclasses.extend(get_all_subclasses(subclass))
    return all_subclasses
    
def get_version_string():
    """
    If installed from github, print the version obtained by `git describe`
    On my local machine, when run within a git directory, get the commit 
    hash directly using `git describe`
    """
    try:
        #Installed with setup.py from a gitrepo
        label = "ernwin {}, forgi {}".format(fess.__complete_version__, forgi.__complete_version__)
    except:
        try: 
            #On my local machine, run from git directory. This script issues `git describe` in the ernwin and forgi directory.
            label = subprocess.check_output(["get_ernwin_version"])
        except OSError: 
            #In production, use the version variable
            label = "ernwin {}, forgi {}".format(__version__, fgb.__version__)
    return label
