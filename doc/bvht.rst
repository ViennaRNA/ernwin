
Modelling a long noncoding RNA with ERNWIN using experimental SAXS data
=======================================================================

As a last, more complex example, this section will show how to sample 3D structures of the long 
non-coding RNA Braveheart (bvht) using experimental SAXS data, as we did in the paper 
"Zinc-finger protein CNBP alters the 3-D structure of lncRNA Braveheart in solution" :ref:`[1]<ref_bvht>`.

In this section you will learn how to use artificial structures as fallback,
if the knowledge base of ernwin is insufficient, 
and how to use SAXS data with ernwin.

A secondary structure of braveheart which we can use as input to our sampling is included in the
RESOURCES folder of the ernwin git repository. 


Challenges with lncRNAs and braveheart in particular
----------------------------------------------------

Running an ernwin simulation as described in the previous section with the secondary 
structure of braveheart as input will work, 
but ernwin will issue an error-level warning at the beginning for the following reason:
Braveheart contains long unpaired regions, such as an interiot loop with 14 unpaired nucleotides on one strand. 
Ernwin's fragment library does not contain any fragment for such a big interior loop, 
so ernwin falls back to using fragments corresponding to shorter interior loops.

.. code-block:: text

  WARNING:fess.builder.stat_container._possible_stats_inner[247]: Only 0 angle-stats found for interior loop with 14 and 3 unpaired nucleotides
  ERROR:fess.builder.stat_container._possible_stats_inner[282]: Trying key interior loop with 13 and 3 unpaired nucleotides instead of interior loop with 14 and 3 unpaired nucleotides for angle-stat

This fallback mechanism works for sampling, but will cause the all-atom reconstructions to have a missing residue.

If you do not want to accept missing residues in the sampled structure, you can extend ernwin's 
fragment library with fragments from artificial structures. 
For braveheart, we have already created artificial fragments and made them available for download.
Over 50000 short artificial structures were generated  :ref:`[2] <ref2>` mostly with 
Rosetta :ref:`[3] <ref3>` (with some additional structures generated with SimRNA :ref:`[4] <ref4>`  
or iFoldRNA :ref:`[5] <ref5>`). See below for the download links.

If you want to create your own artificial fallback structures 
(which goes beyond the scope of this tutorial), you can do so as follows:

#. Create PDB structures of short fragments (individual loops) with a software of your choice.
#. Create a representation in the forgi format for each of them, using forgi's `rnaConvert.py` script.
   Either let it auto-detect the secondary structure, or enforce a secondary structure using the `--pdb-secondary-structure` option.
#. Create a stats file from the forgi files, using the `create_stats.py` script from the ernwin repository. Assuming all forgi files from the previous step are in the current directory and have the `.cg` file ending, the command would look something like this::

    python3 fess/scripts/create_stats.py *.cg 2>/dev/null >../my_fallback.stats

In this tutorial, you can continue with just downloading these fallback files, 
as described in the next section.

Download the necessary files
----------------------------

The bvht secondary structure is included in the RESOURCES folder of the ernwin repository::

  wget https://raw.githubusercontent.com/ViennaRNA/ernwin/master/RESOURCES/bvht.fa
  
To use our fallback fragments to enrich the library of known structures with artificial structures, 
you have to download 3 files:
The artificial structures `in PDB format <www.tbi.univie.ac.at/~thiel/fallback_pdbs.tar.gz>`_ 
and `in cg format <www.tbi.univie.ac.at/~thiel/fallback_cgs.tar.gz>`_ 
and a `stats file <www.tbi.univie.ac.at/~thiel/fallback.stats.gz>`_ where the fragments of the 
structures are summarized for ernwin.
You can download and unzip these files either manually or via the command line::


  wget www.tbi.univie.ac.at/~thiel/fallback_cgs.tar.gz
  tar -xzf fallback_cgs.tar.gz  
  mv fallback_cgs/*.cg /path/to/ernwin_data/CGS/
  wget www.tbi.univie.ac.at/~thiel/fallback_pdbs.tar.gz
  tar -xzf fallback_pdbs.tar.gz
  mv fallback_pdbs/*.pdb /path/to/ernwin_data/PDBs/
  wget www.tbi.univie.ac.at/~thiel/fallback.stats.gz
  gunzip fallback.stats.gz

To use these fallback structures in sampling, give the path to the fallback.stats file 
via the comman line::

  --fallback-stats-files /path/to/fallback.stats
  
  
Additionally, for all atom reconstruction, you have to place the cg and pdb files in the directories
which ernwin uses for reconstruction (the `--source-cg-dir` and `--source-pdb-dir` options). 
In these directories, the artificial structures will be next to the real structures 
that you put there when following the quick start section of the tutorial.


Finally, to sample against experimental SAXS data, obtain the SAXS file from SASBDB:

wget https://www.sasbdb.org/media/p_of_R_files/SASDHE3.out --no-check-certificate

Note that the p(r) in this file has a value of 260, which is clearly in angstrom. 
As ernwin by default assumes nm for GNOM output files, we will have to tell ernwin 
the correct unit via a separate command line argument (`--gnom-unit A`).

Modelling with experimental SAXS data
-------------------------------------

Experimental SAXS data can be given either as CSV file or as GNOM output file. 
Additionally, the pair distance distribution has to be referenced in the energy option::


  --energy PDD6[R],SLD,AME --pdd-file FILE.csv

If the file is in csv format, it should look like this:

  distance,count,error,
  0.0000E+00,0.0000E+00,0.0000E+00
  0.2000E+01,0.1665E-05,0.3661E-07
  0.4000E+01,0.3586E-05,0.5092E-07
  0.6000E+01,0.5714E-05,0.4946E-07
  0.8000E+01,0.7981E-05,0.4338E-07

Slow building of the initial structure
--------------------------------------

Due to the large size of bvht, finding an initial structure as starting point for sampling 
can take a lot of time (minutes to hours). 
It is often possible to speed this up using an experimental building mechanism (`--relaxation-builder`).

Putting it all together
-----------------------

Assuming you have obtained all needed files as described above, the following command can be used::

  ernwin.py bvht.fa --iter 5000 --reconstruct-every-n 100 --source-pdb-dir ernwin_data/PDBs --source-cg-dir ernwin_data/CGS --fallback-stats-files ffallback/fallback.stats --pdd-file SASDHE3.out --gnom-unit A --energy PDD[R],LLI,AME --constraint-energy-per-ml JDIST --relaxation-builder --seed 999 --reconstruction-cache-dir ~/.cache/ernwin 

With this command and seed 999, on our Linux machine running python 3.11, the command took 7 hours for sampling 5000 steps, 
a large part of which is spent on all-atom reconstruction. 
The frequency of all-atom reconstruction can be controlled with the `--reconstruct-every-n` option. 
Additionally, construction of the initial model takes several minutes. 
Thus, the first all-atom structure after 100 sampling steps was available after roughly 15 minutes runtime.

Some structures generated with this run are available in 
the folder RESOURCES/bvht_output of the ernwin git repository.

Using Crysol from the ATSAS package, you can evaluate the Chi^2 of the sampled structures. 
In this run, the first reconstructed structure (after 100 steps) had a Chi^2 of 12.8, 
the structure after 2700 steps had the lowest Chi^2 of the trajectory with a value of 6.05, 
and the last structure (after 5000) steps had a value of 9.750.
If you reach Chi^2 values in this range for RNA molecules with several 100 nucleotides, 
then this shows that Ernwin is installed correctly, correctly uses the pair distance distribution
as energy potential (correct unit etc) 
and the secondary structure is at least roughly compatible with the SAXS profile.
By running multiple simulations, tweaking the parameters of the energy function 
and increasing the number of steps per simulation, it is often possible to find even 
better fitting structures (e.g. with Chi^2 of 1.7 to 2.6 for bvht, 
as we reported in :ref:`[1]<ref_bvht>`).



References
----------

.. _ref_bvht:

[1] *Kim, D.N., Thiel, B.C., Mrozowich, T. et al.*
**Zinc-finger protein CNBP alters the 3-D structure of lncRNA Braveheart in solution.**
Nat Commun 11, 148 (2020). https://doi.org/10.1038/s41467-019-13942-4

.. _ref2:

[2] *Peter Kerpedjiev, Christian Höner zu Siederdissen and Ivo L. Hofacker*.
**Predicting RNA 3D structure using a coarse-grain helix-centered model.**
RNA (2015) 21:1110-1121.

.. _ref3:

[3] *R. Das and D. Baker*.
**Automated de novo prediction of native-like RNA tertiary structures.**
Proc Natl Acad Sci (2007) 104:14664-14669

.. _ref4:

[4] *Boniecki MJ, Lach G, Dawson WK, Tomala K, Lukasz P, Soltysinski T, Rother KM, Bujnicki JM*
**SimRNA: a coarse-grained method for RNA folding simulations and 3D structure prediction**
Nucleic Acids Res 2015 [doi: 10.1093/nar/gkv1479]

.. _ref5:

[5] *S. Sharma, F. Ding, and N. V. Dokholyan*
**iFoldRNA:Three-dimensional RNA structure prediction and folding**
Bioinformatics 2008, 24: 1951-1952
