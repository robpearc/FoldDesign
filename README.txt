################################################################################
              ______    _     _______          _             
              |  ___|  | |   | |  _  \        (_)            
              | |_ ___ | | __| | | | |___  ___ _  __ _ _ __  
              |  _/ _ \| |/ _` | | | / _ \/ __| |/ _` | '_ \ 
              | || (_) | | (_| | |/ /  __/\__ \ | (_| | | | |
              \_| \___/|_|\__,_|___/ \___||___/_|\__, |_| |_|
                                                  __/ |      
                                                 |___/    


                      (Version 1.0, 02/14/2022)

(Copyrighted by the Regents of the University of Michigan, All rights reserved)

FoldDesign is a program for de novo protein fold design

Author: Robin Pearce

To report bugs and questions email: robpearc@umich.edu

If you use this program, please cite:
Pearce R., Huang X., Omenn G.S., and Zhang Y. De Novo Protein Fold Design
by Assembling Structural Fragments Under the Guidance of an Optimized
Force Field. Submitted, 2022.

Note, this is the stand-alone program, users may also submit jobs to 
https://zhanggroup.org/FoldDesign/

The source code is freely available to academic/non-profit users under the 
PolyForm Noncommercial License

################################################################################



######################### SYSTEM REQUIREMENTS ##################################

    System requirements: x86_64 machine, Linux kernel OS with perl and python 
    interpreters installed. 

################################################################################



########################### INSTALLATION #######################################

    Before you run FoldDesign, you must download the library by running the command 
    scripts/install_library.sh. The library should be installed in the main package 
    directory.

################################################################################



####################### RUNNING INSTRUCTIONS ###################################

    To run FoldDesign use the following command:


    perl runFoldDesign.pl -datadir=<data directory> 


    The only required input is the path to the directory that contains 
    the following file:
        
        input.txt     This file contains the input sequence and secondary 
                      structure information.


    The format of the input.txt file should be as follows:

        V H
        V H
        V E
        V E
        V C

    where the first character is the one-letter amino acid code and the second
    character is the desired secondary structure (H->helix, E->beta-strand, C->coil)
    at the given position. The program accepts any of the 20 naturally occuring amino
    acids, which will primarly be used to evaluate steric clashes between the side-chain
    centers of mass. By default, Valine may be used as a placeholder residue.


    Additionally, the following files may be present in the input directory:
        
        contact_restr.txt     This file contains any inter-residue contact restraints.
                              A contact is defined as a pair of atoms from two residues
                              that should be within 8 angstroms of each other.
       
        distance_restr.txt     This file contains any inter-residue distance restraints.
                               A distance restraint specifies the distance between a pair
                               of atoms from two residues.


    The format of the contact restraints should be as follows:
        
        9 25 1.0 4 4
        10 24 2.0 0 4
        14 30 0.5 0 0

    where the first 2 characters are the indices for residue i and j and the 3rd
    character is the weight of the contact restraint. Assigning a higher weight to 
    a contact restraint will increase its importance during the design simulations. 
    Additionally, the last two characters specify the atom type the restraints will be
    enforced for, where the 4th character is the atom type for residue i and the
    5th character is the atom type for residue j. A value of 0 means the contact
    will be enforced for the CA atom, while a value of 4 means the contact will
    be enforced for the CB atom.


    The format of the distance restraints should be as follows:
        
        9 25 3.0 10.0 4 4 0
        10 24 0.75 6.9 0 4 0
        14 30 1.0 7.1 0 0 1

    where the first 2 characters are the indices for residue i and j, the 3rd
    character is the weight of the distance restraint, and the 4th character is
    the distance between the two atoms from the residue pair. Again, assigning a 
    higher weight to a distance restraint will increase its importance during the 
    design simulations. The 5th and 6th characters specify the atom type the restraint 
    should be enforced for, where the 5th character is the atom type for residue i 
    and the 6th character is the atom type for residue j. A value of 0 means the 
    distance restraint will be enforced for the CA atom, while a value of 4 means 
    the distance restraint will be enforced for the CB atom. The final character
    is the mathematical function that will be used to enforce the distance restraint.
    Two options are available, either a strong harmonic restraint (0), or a weak
    reciprocal square restraint (1). If you specify few distance restraints, it
    is probably better to use harmonic restraints, while if you specify many restraints,
    it is probably better to use reciprocal square restraints.



    Optional arguments:

    -random_num=<NUM>     Set a specific random number to generate different 
                          designed structures. Default: 102

    -num_remc_cycles=<NUM>     Number of REMC cycles to perform during the 
                               design simulations. Default: 500

    -design_all_clusters=<(True/False)>     Whether or not to perform sequence design 
                                            on all the clusters generated by the FoldDesign 
                                            simulations. By default, sequence design will 
                                            only be performed for the lowest energy design.

    -simulation_timeout=<NUM>     Terminate the FoldDesign simulations after a specified 
                                  number of hours. By default the simulations will terminate 
                                  after 72 hours or after the given number of REMC cycles
                                  have been completed, whichever occurs first.


    The final designed structure(s) will be saved in the <datadir>/final_designs folder,
    where all designed structures will be saved as all_decoys.pdb.tar.bz2 in the 
    same folder. Note, the structures in the all_decoys.pdb.tar.bz2 file will not have 
    designed sequences or refined structures. To perform sequence design and refinement, 
    EvoDesign and ModRefiner, located in the bin directory, may be used on any of the 
    selected designs.

################################################################################
