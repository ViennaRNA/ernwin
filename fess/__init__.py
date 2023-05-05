from __future__ import absolute_import
import os

__author__ = "Peter Kerpedjiev, Bernhard Thiel"
__copyright__ = "Copyright 2012 - 2016"
__license__ = "GNU Affero GPL v 3.0"
__version__ = "1.2.0"
__maintainer__ = "Bernhard Thiel, Peter Kerpedjiev"
__email__ = "thiel@tbi.univie.ac.at, pkerp@tbi.univie.ac.at"

def data_file(fname):
    """Return the path to a data file of ours."""
    return os.path.join(os.path.split(__file__)[0], fname)
