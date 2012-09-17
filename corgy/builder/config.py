import os

class Configuration:
    base_dir = os.path.expanduser('~/projects/ernwin')
    data_base_dir = os.path.expanduser('~/projects/ernwin/fess/output')
    stats_file = os.path.join(base_dir, 'fess/stats/temp.stats')
    stem_fragment_dir = os.path.join(base_dir, 'fess/stats/stems')
    lric_stats_fn = os.path.join(base_dir, 'fess/stats/temp.energy')
    template_residue_fn = os.path.join(base_dir, 'fess/stats/residue_template.pdb')
    longrange_contact_stats_fn = os.path.join(base_dir, 'fess/stats/temp.longrange.contact')

    test_input_dir = os.path.join(base_dir, "test_input")
    test_output_dir = os.path.join(base_dir, "test_output")
    barnacle_dir = '/scr/plastilin/pkerp/apps/Barnacle'
