Setup of the system
===================

After installation of ernwin and its dependencies, you can start sampling right away with the default parameters and generate forgi-compatible coarse grained structures. However, to reconstruct all-atom structures, you need download additional files. 

Ernwin comes with the file all_nr3.36.stats which contains angles and distances between secoondary structure elements of solved RNA structures. This is enough to sample coarse grained structures, but for all-atom reconstruction, you in addition need the pdb files which were used to create this stats file.
Create a directory, lets call it PDB_DIR, and put all pdb files referenced in the representative set of RNA structures, Version 3.36, into it (in cif format)

PDB-files do not store the secondary structure of the RNA chains and different structure annotation tools have different results. To ensure consistency, we have deposited the secondary structures and corresponding residue numbers in forgi format in this repository. When you extract the file CGS.tar.gz, this will create the directory CGS.


You can now run simulations with the following command::

    ernwin.py INPUT_FILE --source-cg-dir CGS --source-pdb-dir PDB_DIR --reconstruct-every-n 1


Modelling a long noncoding RNA with ERNWIN
==========================================

To reproduce the results shown in the paper about bvht by Kim, Thiel, Mrozowich, Hennelly, Hofacker, Patel and Sanbonmatsu (submitted), additional resources are needed: Over 50000 short artificial structures were generated  :ref:`[1] <ref1>` mostly with Rosetta :ref:`[2] <ref2>` (with some additional structures generated with SimRNA :ref:`[3] <ref3>`  or iFoldRNA :ref:`[4] <ref4>`) and used as fallback fragments, when less than 100 examples of a secondary structure element could be found in the PDB. This was especially the case for larger exterior loop segments and other large unpaired loops, which are underrepresented in the PDB due to their flexibility. These artificial PDB structures can be downloaded from www.tbi.univie.ac.at/~thiel/fallback_pdbs.tar.gz (Only needed for all-atom reconstruction). The corresponding forgi files are at www.tbi.univie.ac.at/~thiel/fallback_cgs.tar.gz . For reconstruction to work, put the files into your folders PDB_DIR and CGS respectively.

To use them in sampling, the file fallback.stats from the folder RESOURCES in this repository has to be used with the following commandline argument::

  --fallback-stat /path/to/fallback.stats



Modelling with experimental SAXS data
-------------------------------------

To give an experimental SAXS data file to ernwin, use the following commandline arguments::

  --energy PDD[R],SLD,AME --pdd-file FILE.csv

Where the file FILE.csv should be of the following format::

  distance,count,error,
  0.0000E+00,0.0000E+00,0.0000E+00
  0.2000E+01,0.1665E-05,0.3661E-07
  0.4000E+01,0.3586E-05,0.5092E-07
  0.6000E+01,0.5714E-05,0.4946E-07
  0.8000E+01,0.7981E-05,0.4338E-07


References
----------

.. _ref1:

[1] *Peter Kerpedjiev, Christian Höner zu Siederdissen and Ivo L. Hofacker*.
**Predicting RNA 3D structure using a coarse-grain helix-centered model.**
RNA (2015) 21:1110-1121.

.. _ref2:

[2] *R. Das and D. Baker*.
**Automated de novo prediction of native-like RNA tertiary structures.**
Proc Natl Acad Sci (2007) 104:14664-14669

.. _ref3:

[3] *Boniecki MJ, Lach G, Dawson WK, Tomala K, Lukasz P, Soltysinski T, Rother KM, Bujnicki JM*
**SimRNA: a coarse-grained method for RNA folding simulations and 3D structure prediction**
Nucleic Acids Res 2015 [doi: 10.1093/nar/gkv1479]

.. _ref4:

[4] *S. Sharma, F. Ding, and N. V. Dokholyan*
**iFoldRNA:Three-dimensional RNA structure prediction and folding**
Bioinformatics 2008, 24: 1951-1952
