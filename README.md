# Graphon_estimation_reproducibility

To simulate the 100 sample paths of 7 different graphon Hawkes processes (the unperturbed and 3 perturbed processes referred to in the document), run parallel_estimation_legendre_complete.m
For each of those 7 sets of 100 realized paths, run lines 213—241 of the same code to reproduce the contents of Appendix C, as well as the numbers in Table 2 of the main text.
Table 1 can be reproduced by running  parallel_estimation_legendre.m for the different values of T. Note that the user can alter commands like
%delete(gcp('nocreate'));
%parpool(7);
To prevent an overflow of RAM. If there is no issue with RAM overflow, you could simply surpress.

The data for the session 715093703 of the Neuropixels extracellular recording collection maintained by the Allen Institute for Brain Science can be found through the website of Allen Institute, yielding ecephys_session_715093703.nwb. We used Python to extract an .m file from this. The resulting .m file, consisting of the events in the time interval [5000;5250] is found in Allen_neuropixel.zip (unzip to get Allen_neuropixel.m).

Figure 1 can be reproduced by running allen_brain_legendre.m (for Legendre basis) and allen_brain.m (for Haar basis). Lines 271—288  and 172—189, respectively, recreate the required figures. 
Finally, Figure A3 can be reproduced by running allen_brain.m with eta=0.
