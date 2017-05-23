from distutils.core import setup
from distutils.command.build_py import build_py as _build_py
import subprocess
import os
        

try: #If we are in a git-repo, get git-describe version.
    path = os.path.abspath(os.path.dirname(__file__))    
    ernwin_version = subprocess.check_output([os.path.join(path, "git"), "describe", "--always"]).strip()
    try:
        subprocess.check_call([os.path.join(path, "git"), 'diff-index', '--quiet', 'HEAD', '--'])
    except subprocess.CalledProcessError:
        ernwin_version+="+uncommited_changes"
    #Use a subclass of build_py from distutils to costumize the build.
    class build_py(_build_py):
        def run(self):
            """
            During building, adds a variable with the complete version (from git describe)
            to fess/__init__.py.
            """
            try: 
                outfile = self.get_module_outfile(self.build_lib, ["fess"], "__init__")
                os.remove(outfile) #If we have an old version, delete it, so _build_py will copy the original version into the build directory.
            except: 
                pass
            # Superclass build
            print("Running _build_py")
            _build_py.run(self)
            print("Postprocessing")
            outfile = self.get_module_outfile(self.build_lib, ["fess"], "__init__")
            # Apped the version number to init.py
            with open(outfile, "a") as of:
                of.write('\n__complete_version__ = "{}"'.format(ernwin_version))
except OSError: #Outside of a git repo, do nothing.
    build_py = _build_py




setup(cmdclass={'build_py': build_py},
      name='ernwin',
      version='0.1',
      description='Coarse Grain 3D RNA Structure Modelling',
      author='Peter Kerpedjiev, Bernhard Thiel',
      author_email='pkerp@tbi.univie.ac.at, thiel@tbi.univie.ac.at',
      url='http://www.tbi.univie.ac.at/~thiel/ernwin/',
      packages = ['fess', 'fess.builder', 'fess.aux', 'fess.motif'],
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
                             'stats/aminors_1jj2_sampled.csv',
			     'stats/combined.stats',
                             'stats/all.stats',
                             'stats/all_nr2.110.stats']},
      scripts=['fess/scripts/ernwin_new.py']
     )
