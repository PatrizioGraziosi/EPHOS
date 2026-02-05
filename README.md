# EPHOS
electron-phonon coupling in organic semiconductors

Pool of scripts to compute and parametrize the e-ph coupling in OSCs, under-testing version.

Below I describe the workflow, I assume that the user can run independently DFT (here VASP is assumed, but not mandatory) and phonons calculations with phonopy. 

- After having used at least a 2x2x2 supercell, $ phonopy dos.conf computes the phonon DOS that allow you to compute the DOS and see up to which frequency you consider the phonons. Rule of thumb is to take all the modes up to the firdt drop-to-zero of the DOS. These will include also some more modes than the nominal Z*6, for this reason I used the "pseudo_inter" expression.
- After having identified the modes and q-points, modulate the structure. The file setting_mod_pseudo_inter.conf represent an exaple, while the script write_lines_mod.m allow to print in the matlab prompt all the lines necessary for the modulation. It is a good practice to study the convegence choosing for example 10 (q-point,mode) pairs and observing how the EPC varies with the displacment tag, that means to run all the sequence with only a few random pairs for several displacement, then run again for all the pairs.
- use the make_n_run_modulated_complete_pseudo_inter script to run DFT over all the modulated structure; mind the instruction on the nuber of files
- use the script extract_freq.sh get the frequencies, you MUST set correctly the number of modes, points, and so...
- Finally, compute the EPC with the launch script analyze_extract_OSC_DefPot_v5.sh, mind the information in the first lines.
- After this, the script overall_EPC_OSC_separate_branches_v1.m allows to make the scattering parameter file for ELECTRA, which will must be run in the 'multivalley' mode.



