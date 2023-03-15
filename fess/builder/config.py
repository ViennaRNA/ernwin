from __future__ import absolute_import
import os

class Configuration:
    mids_method="template"
    
    sampling_output_dir = 'sampling_out'

    default_stats_file = "stats/all_nr3.36.stats" # Can be overridden by commandline arguments.

    #The file where ideal bases are stored for base-replacement during reconstruction
    template_residue_fn = 'stats/residue_template.pdb'
    template_residues = None # This is changed programatically. It should always be None in this file.

    # The JAR3D dir should contain the jar3d-related files defined in the following 3 configuration values.
    # Further more, ernwin will create the subdirectories 'cgs' and 'pdbs'
    jar3d_dir = '~/JAR3D'
    jar3d_jar = 'jar3d_2014-12-11.jar' #Download from http://rna.bgsu.edu/data/jar3d/models/ and strore in jar3d_dir
    jar3d_motif = 'Internal Loop Motif Atlas Release 1.18.json' #Download from http://rna.bgsu.edu/rna3dhub/motifs/release/il/current and strore in jar3d_dir
    jar3d_IL = 'IL/1.18/lib/all.txt' #Download from http://rna.bgsu.edu/data/jar3d/models/ and strore in jar3d_dir
