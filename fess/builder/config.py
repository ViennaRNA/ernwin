import os

class Configuration:
    mids_method="template"
    #mids_method="basenormals"
    #base_dir = os.path.expanduser('~/projects/ernwin')
    #data_base_dir = os.path.expanduser('~/doarse/processed/')
    #pdb_base_dir = os.path.expanduser('~/data/ernwin/pdb')
    #stats_file = os.path.join(base_dir, 'fess/stats/temp.stats')
    #stem_fragment_dir = os.path.join(base_dir, 'fess/stats/stems')
    #lric_stats_fn = os.path.join(base_dir, 'fess/stats/temp.energy')
    #template_residue_fn = os.path.join(base_dir, 'fess/stats/residue_template.pdb')
    #longrange_contact_stats_fn = os.path.join(base_dir, 'fess/stats/temp.longrange.contact')
    #test_input_dir = os.path.expanduser('~/data/ernwin/processed/')
    #test_output_dir = os.path.join(base_dir, "test_output")
    sampling_output_dir = 'sampling_out'
    #barnacle_dir = '/scr/plastilin/pkerp/apps/Barnacle'
    #template_residues = None
    #stem_library = dict()
    
    # The JAR3D dir should contain the jar3d-related files defined in the following 3 configuration values.
    # Further more, ernwin will create the directories 'cgs' and 'pdbs'
    jar3d_dir = '~/JAR3D'
    jar3d_jar = 'jar3d_2014-12-11.jar' #Download from http://rna.bgsu.edu/data/jar3d/models/ and strore in jar3d_dir
    jar3d_motif = 'Internal Loop Motif Atlas Release 1.18.json' #Download from http://rna.bgsu.edu/rna3dhub/motifs/release/il/current and strore in jar3d_dir
    jar3d_IL = 'IL/1.18/lib/all.txt' #Download from http://rna.bgsu.edu/data/jar3d/models/ and strore in jar3d_dir

