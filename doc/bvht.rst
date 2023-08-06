
Modelling a long noncoding RNA with ERNWIN
==========================================

To reproduce the results shown in the paper "Zinc-finger protein CNBP alters the 3-D structure of lncRNA Braveheart in solution" :ref:`[1]<_ref_bvht>`, additional resources are needed: 
Over 50000 short artificial structures were generated  :ref:`[2] <ref2>` mostly with Rosetta :ref:`[3] <ref3>` (with some additional structures generated with SimRNA :ref:`[4] <ref4>`  or iFoldRNA :ref:`[5] <ref5>`) and used as fallback fragments, when less than 100 examples of a secondary structure element could be found in the PDB. This was especially the case for larger exterior loop segments and other large unpaired loops, which are underrepresented in the PDB due to their flexibility. These artificial PDB structures can be downloaded from www.tbi.univie.ac.at/~thiel/fallback_pdbs.tar.gz (Only needed for all-atom reconstruction). The corresponding forgi files are at www.tbi.univie.ac.at/~thiel/fallback_cgs.tar.gz . For reconstruction to work, put the files into your folders PDB_DIR and CGS respectively.

To use them in sampling, the file fallback.stats from the folder RESOURCES in this repository has to be used with the following commandline argument::

  --fallback-stat /path/to/fallback.stats



Modelling with experimental SAXS data
-------------------------------------

To give an experimental SAXS data file to ernwin, use the following commandline arguments::

  --energy PDD6[R],SLD,AME --pdd-file FILE.csv

Where the file FILE.csv should be of the following format::

  distance,count,error,
  0.0000E+00,0.0000E+00,0.0000E+00
  0.2000E+01,0.1665E-05,0.3661E-07
  0.4000E+01,0.3586E-05,0.5092E-07
  0.6000E+01,0.5714E-05,0.4946E-07
  0.8000E+01,0.7981E-05,0.4338E-07

The number 6 after the string PDD specifies how closely the predicted ensemble should fit the SAXS data file. Values between 1 and 10 are reasonable.

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
