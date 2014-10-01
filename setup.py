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
                            'stats/temp.longrange.stats', 
                             'stats/temp.longrange.random_radius_of_gyration_beta_16.stats',
                             'stats/loop_loop3_distances_sampled.csv', 
                             'stats/loop_loop3_distances_native.csv',
                             'stats/tall.csv',
                             'stats/all_elements.csv',
                             'stats/aminors_1s72.csv',
                             'stats/aminors_1jj2_sampled.csv']},
      scripts=['fess/scripts/ernwin_go.py']
     )
