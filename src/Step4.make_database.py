from dssp_handling import read_dssp, CHECK_HBONDS, SS_SEGMENT, AA_SEGMENT
from os.path import exists
import pickle
import copy
from path_setup import DB_PATH, DSSP_PATH, DERIVED_PATH

"""

This script reads the DSSP files and makes dictionaries with this format: 

D[pdbc][segment] = [A, B]:
     pdbc    = pdb ID + chain code
     segment = unique 11mer sequence
     A       = number of hairpins identified for this sequence
     B       = number of times this sequence is observed

Two dictionaries are saved as pickles

HAIRPINMOTIFS_NONR_PDB_SequenceStats.p4
    -> All 11mers with at least one hairpin observed
   
STPXG_NONR_PDB_SequenceStats.p4
    -> All 11mers with a [TS]PXG motif starting from residue i+3
    
"""

# These files are the specific dssp files used in preparing the manuscript
dfiles = open(DB_PATH+"dssp.files").readlines()

if DSSP_PATH == "":
    print("ERROR: You need to set the DSSP_PATH variable")
    exit()

if not exists(DSSP_PATH):
    print("ERROR: You need to download the DSSP database (~20gb)")
    print("     rsync -avz rsync://rsync.cmbi.umcn.nl/dssp/ "+DSSP_PATH)
    exit()

nonredundant_pdbs = {}
pisces_cull_file = open(DB_PATH+"cullpdb_pc90_res3.0_R1.0_inclNOTXRAY_d191216_chains49938.22503").readlines()
for p in pisces_cull_file:
    pdbc = p.split()[0]
    if len(pdbc) == 5:
        pdb = pdbc[:4].lower()
        chain = pdbc[4]
        if pdb not in nonredundant_pdbs:
            nonredundant_pdbs[pdb] = {}
        nonredundant_pdbs[pdb][chain] = 0

redundancy_map = {}
pisces_cull_file = open(DB_PATH+"cullpdb_pc90_res3.0_R1.0_inclNOTXRAY_d191216_chains49938.log.22503").readlines()
for p in pisces_cull_file:
    """
    The pisces server log file is 573216 lines of:
    "reject 3NIRA 3U7TA  96"
    
    which specifically means that 3NIR chain A has NOT been rejected
    (it is the high quality representative that has been selected)
    
    But 3U7T chain A has been rejected, because it has 96% sequence identity to 3NIR chain A
    
    """
    nr_pdb = p.split()[1][:4].lower()
    if nr_pdb in nonredundant_pdbs:
        nr_chain = p.split()[1][4]
        if nr_chain in nonredundant_pdbs[nr_pdb]:

            redundant_pdbc = p.split()[2][:4].lower()+'.'+p.split()[2][4]
            redundancy_map[redundant_pdbc] = nr_pdb+"."+nr_chain
            # This is what most map lines look like, redundant pdb+chain maps to the selected representative

            if nr_pdb+'.'+nr_chain not in redundancy_map:
                redundancy_map[nr_pdb+'.'+nr_chain] = nr_pdb+"."+nr_chain
                # Might as well also just put map the selected representatives to themselves

# STEP1: Pull out hairpin data from the DSSP database

ALLSEGS = {}
NONRSEGS = {}
for di, dfile in enumerate(dfiles):

    pdb = dfile.split('.')[0].lower()

    if di % 1000 == 0:
        print("Working on %i of %i" %(di+1, len(dfiles)))#len(nonredundant_pdbs)))

    if not exists(DSSP_PATH+dfile.split()[0]):
        print("WARNING, YOU ARE MISSING A DSSP FILE", dfile.split()[0])
    else:

        dssp = read_dssp(DSSP_PATH+dfile.split()[0])

        for c in dssp:

            pdbc = dfile.split('.')[0]+"."+c
            for d in dssp[c]:

                SSSEG = SS_SEGMENT( dssp, c, d, 11 )
                AASEG = AA_SEGMENT( dssp, c, d, 11 )

                if len(AASEG) == 11:

                    HasMotif = False
                    if SSSEG[1:-1] == "EELTTLLEE":
                        if CHECK_HBONDS( dssp, c, d ):
                            HasMotif = True

                    if pdbc in redundancy_map:
                        """
                        
                        Note 1: Statistics are collected for non-redundant PDB chains...
                        (at 90% sequence identity and resolution <= 3.0 A from pisces culling server)
                        
                        but the frequency of hairpin observation is based on all PDB chains that map to the
                        non-redundant PDB chain in question
                        
                        Ie: If a given 11 residue motif is observed in 6 PDBs with the same sequence but only
                            4 of those have the right hydrogen bond patten and DSSP assignment to count as hairpin 
                            turns then for the one nonredundant PDB in the set we store the motif as [4, 6]
                        
                        and then we filter out hairpins that aren't observed in at least 50% of the structures.
                        
                        Note 2: motifs remain mapped to nonredundant pdbs in order to facilitate bootstrapping
                        
                        """

                        nonr_pdbc = redundancy_map[pdbc]

                        # Database Initialization
                        if nonr_pdbc not in NONRSEGS:
                            NONRSEGS[nonr_pdbc] = {}
                        if AASEG not in NONRSEGS[nonr_pdbc]:
                            NONRSEGS[nonr_pdbc][AASEG] = [0,0]

                        # Counting
                        NONRSEGS[nonr_pdbc][AASEG][1] += 1
                        if HasMotif:
                            NONRSEGS[nonr_pdbc][AASEG][0] += 1

###
# All 11 mers observed in the nonredundant PDB set
pickle.dump(NONRSEGS, open(DERIVED_PATH+"ALL_NONR_PDB_SequenceStats.p4", 'wb'))

"""

For doing statistics on just the observed hairpins we make a new data set by just removing
the 11-mers that did not have hairpin observations.

"""

HAIRPINS = copy.deepcopy(NONRSEGS)
for pdbc in NONRSEGS.keys():
    for motif in NONRSEGS[pdbc].keys():
        if float(NONRSEGS[pdbc][motif][0])/float(NONRSEGS[pdbc][motif][1]) < 0.5:
            del HAIRPINS[pdbc][motif]
    if len(HAIRPINS[pdbc].keys()) == 0:
        del HAIRPINS[pdbc]
pickle.dump(HAIRPINS, open(DERIVED_PATH+"HAIRPINMOTIFS_NONR_PDB_SequenceStats.p4", 'wb'))

"""

For doing statistics on just the [ST]PXG containing 11-mers we make that dataset too

"""

STPXG = copy.deepcopy(NONRSEGS)
for pdbc in NONRSEGS.keys():
    for motif in NONRSEGS[pdbc].keys():
        MSEG = motif[3] + motif[4]
        MSEG += "X"
        MSEG += motif[6]
        if MSEG not in ("SPXG", "TPXG"):
            del STPXG[pdbc][motif]
    if len(STPXG[pdbc].keys()) == 0:
        del STPXG[pdbc]
pickle.dump(STPXG, open(DERIVED_PATH+"STPXG_NONR_PDB_SequenceStats.p4", 'wb'))
