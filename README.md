ernwin 1.1
==========

Since release 1.1 , ernwin supports python3

Installation
------------

Ernwin can be installed from the python package index using pip:
    
    pip install ernwin

On most Linux distribution this will work out of the box. 
If you get an error complaining about an unknown extension "pyx", try 
installing Cython first and then re-try the pip installation of ernwin.
Cython can either be installed using `pip install Cython`, 
or via your operating system's package manager.

If you experience any problems with the installation of ernwin, 
please open an issue at github.

Optional dependency
-------------------

To visualize sampled structures in the ernwin format with the `visualize_rna.py` script, 
PyMOL <https://pymol.org/> is required (tested with PyMOL 2.5.0).

Post-installation set-up
------------------------

If you just want to sample coarse-grained structures in the ernwin format, 
ernwin can be used directly after installation.
But in order to perform all-atom reconstruction 
and generate PDB files which can be used by other tools, 
a knowledge base of existing pdb structures has to be downloaded (can take up to an hour).
This can be done using the following steps (tested using a bash terminal on Linux -
please use the equivalent commands if you are on a different operating system):

1. Create a folder where you will put your ernwin data:

    ```mkdir ernwin_data```

2. Download  and extract the coarse grained representation of the PDB structure
   knowledge base into this folder (as an alternative you could clone the repository):

    ```cd ernwin_data
    wget https://github.com/ViennaRNA/ernwin/raw/master/RESOURCES/CGS.tar.gz
    tar -xzf CGS.tar.gz```

3. To download the corresponding PDB files in MMCIF format directly from RCSB,
   you can use a bash script available in the ernwin git repository:

    ```wget https://github.com/ViennaRNA/ernwin/raw/master/RESOURCES/download_pdb_files_for_cg_files.sh
    chmod +x download_pdb_files_for_cg_files.sh
    # In the following command, "CGS" is the folder you have downloaded and extracted before
    # This assumes you are still in the ernwin_data folder.
    # This can take several minutes up to an hour, as it downloads over 2000 PDB files.
    # You will see the output of wget displayed while it is in progress.
    ./download_pdb_files_for_cg_files.sh CGS
    # Finally, unzip all downloaded files:
    gunzip *.gz```

After this set-up, you should have over 2000 PDB files in the folder `ernwin_data` and corresponding files
in the ernwin (and forgi) coarse grain format inside `ernwin_data/CGS`.

Running your first simulation
-----------------------------

Without all-atom resonstruction:
 
    ernwin.py INPUT.fa

To enable all-atom resonstruction, you need to give the paths to the ernwin_data folder 
and the ernwin_data/CGS folder (see post-installation set-up):

    ernwin.py INPUT.fa --reconstruct-every-n 10 --source-pdb-dir <PATH>/ernwin_data --source-cg-dir <PATH>/ernwin_data/CGS

To speed-up all-atom reconstruction, you can allow ernwin to store fragments extracted from the
PDB files in ernwin_data into a new folder. 
This uses disk space but provides a significant speed-up 
as it avoids re-opening huge PDB files multiple times.
Do this by adding the option:

    --reconstruction-cache-dir ~/.cache/ernwin

To specify the number of iterations, add:

    --iter 100

Running a simulation with SAXS data
-----------------------------------

To guide the simulation using SAXS data, you need to give the pair distance distribution function 
as CSV file or as GNOM output file using the following option::

    --pdd-file PDDFILE.out 

And you *must* specify an energy that uses this PDD file:

Simple PDD energy (residue resolution)::

    --energy PDD[R]

Simple PDD energy plus long range interactions::

    --energy PDD[R],LLI,AME

Ensemble based energy (publications describing the difference to the PDD energy are in preparation) plus long range interactions::

    --energy EPD[R],LLI,AME



Other documentations
--------------------

Also see the tutorial at 
http://www.tbi.univie.ac.at/~thiel/ernwin
