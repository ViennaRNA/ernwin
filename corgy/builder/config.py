import os

class Configuration:
    mids_method = "template"
    #mids_method = "fit"
    #mids_method = "superimpose"
    #mids_method = "estimate"
    base_dir = os.path.expanduser('~/projects/ernwin')
    data_base_dir = os.path.expanduser('~/data/ernwin/processed')
    pdb_base_dir = os.path.expanduser('~/data/ernwin/pdb')
    stats_file = os.path.join(base_dir, 'fess/stats/temp.stats')
    stem_fragment_dir = os.path.join(base_dir, 'fess/stats/stems')
    lric_stats_fn = os.path.join(base_dir, 'fess/stats/temp.energy')
    template_residue_fn = os.path.join(base_dir, 'fess/stats/residue_template.pdb')
    longrange_contact_stats_fn = os.path.join(base_dir, 'fess/stats/temp.longrange.contact')

    test_input_dir = os.path.expanduser('~/data/ernwin/processed/')
    test_output_dir = os.path.join(base_dir, "test_output")
    sampling_output_dir = 'best'
    barnacle_dir = '/scr/plastilin/pkerp/apps/Barnacle'
    stem_library = dict()
