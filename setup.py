from distutils.core import setup

setup(name='ernwin',
      version='0.1',
      description='Coarse Grain 3D RNA Structure Modelling',
      author='Peter Kerpedjiev',
      author_email='pkerp@tbi.univie.ac.at',
      url='http://www.tbi.univie.ac.at/~pkerp/ernwin/',
      packages = ['fess', 'fess.builder'],
      package_data={'fess': ['stats/temp.stats', 'stats/cylinder_intersections*.csv',
                            'stats/subgraph_radius_of_gyration_sampled.csv',
                            'stats/subgraph_radius_of_gyration.csv',
                            'stats/temp.longrange.stats', 'temp.longrange.random_radius_of_gyration_beta_16.stats']},
      scripts=['fess/scripts/gibbs.py']
     )
