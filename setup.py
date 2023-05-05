from __future__ import absolute_import
from __future__ import print_function
from distutils.core import setup
import os


setup(
      name='ernwin',
      version='1.2.0',
      description='Coarse Grain 3D RNA Structure Modelling',
      author='Bernhard Thiel, Peter Kerpedjiev',
      author_email='thiel@tbi.univie.ac.at, pkerp@tbi.univie.ac.at',
      url='http://www.tbi.univie.ac.at/~thiel/ernwin/',
      packages = ['fess', 'fess.builder', 'fess.ccd', 'fess.motif'],
      install_requires = [
         "pandas>=0.19",
         "numpy",
         "scipy",
         "biopython",
         "deepdiff",
         "networkx",
         "future",
         "forgi>=2.2.1",
         "logging_exceptions",
         "commandline_parsable>=0.2",
         "scikit-learn",
         "matplotlib",
         "contextlib2",
         "bitarray"
      ],
      package_data={'fess': ['stats/all_nr3.36.stats',
                             'stats/cylinder_intersections*.csv',
                             'stats/cde_reference_dist_nr2.110.csv',
                             'stats/sld_reference_dist_nr2.110.csv',
                             'stats/rog_reference_dist_nr2.110.csv',
                             'stats/rog_target_dist_nr2.110.csv',
                             'stats/sld_target_dist_nr2.110.csv',
                             'stats/sld_reference_dist_1S72_0.csv',
                             'stats/rog_reference_dist_1S72_0.csv',
                             'stats/rog_target_dist_1S72_0.csv',
                             'stats/sld_target_dist_1S72_0.csv',
                             'stats/residue_template.pdb',
                             'stats/AME_distributions.csv']},
      scripts=['fess/scripts/ernwin.py', 'fess/scripts/reconstruct.py', 'fess/scripts/plot_ernwin_pdd.py']
     )
