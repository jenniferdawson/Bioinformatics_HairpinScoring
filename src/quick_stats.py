import pickle
import numpy as np
import os
def allin20( seg ):
    for n in range(11):
        if seg[n] not in ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                          'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'):
            return False
    return True


nonredundant_pdbs = {}
pisces_cull_file = open("../database/cullpdb_pc90_res3.0_R1.0_inclNOTXRAY_d191216_chains49938.22503").readlines()
for p in pisces_cull_file:
    pdbc = p.split()[0]
    if len(pdbc) == 5:
        pdb = pdbc[:4].lower()
        chain = pdbc[4]
        if pdb not in nonredundant_pdbs:
            nonredundant_pdbs[pdb] = {}
        nonredundant_pdbs[pdb][chain] = 0

redundancy_map = {}
pisces_cull_file = open("../database/cullpdb_pc90_res3.0_R1.0_inclNOTXRAY_d191216_chains49938.log.22503").readlines()
for p in pisces_cull_file:
    nr_pdb = p.split()[1][:4].lower()
    if nr_pdb in nonredundant_pdbs:
        nr_chain = p.split()[1][4]
        if nr_chain in nonredundant_pdbs[nr_pdb]:

            redundant_pdbc = p.split()[2][:4].lower()+'.'+p.split()[2][4]

            if os.path.exists("/home/rvernon/DSSP_DOWNLOAD/"+p.split()[2][:4].lower()+".dssp"):
                redundancy_map[redundant_pdbc] = nr_pdb+"."+nr_chain

            if nr_pdb+'.'+nr_chain not in redundancy_map:
                redundancy_map[nr_pdb+'.'+nr_chain] = nr_pdb+"."+nr_chain

representative_pdbcs = {}
for pdbc in redundancy_map:
    nonr = redundancy_map[pdbc]
    if nonr not in representative_pdbcs:
        representative_pdbcs[nonr] = {}
    representative_pdbcs[nonr][pdbc] = 1

DERIVED_PATH = "../derived/"

NONRSEGS = pickle.load(open(DERIVED_PATH+"HAIRPINMOTIFS_NONR_PDB_SequenceStats.p4", 'rb'))
#NONRSEGS = pickle.load(open(DERIVED_PATH+"STPXG_NONR_PDB_SequenceStats.p4", 'rb'))

N_unique = 0
N_100percent = 0
N_majority = 0
N_minority = 0
N_none = 0
N_proteins = {}
for pdbc in NONRSEGS.keys():
    for motif in NONRSEGS[pdbc].keys():
        hitif = 0
        if NONRSEGS[pdbc][motif][1] == 1:
            N_unique += 1
            N_proteins[pdbc] = 1
            hitif += 1
        elif NONRSEGS[pdbc][motif][0] > 1 and NONRSEGS[pdbc][motif][0] == NONRSEGS[pdbc][motif][1]:
            N_100percent += 1
            N_proteins[pdbc] = 1
            hitif += 1
        elif NONRSEGS[pdbc][motif][0] > 0 and float(NONRSEGS[pdbc][motif][0])/float(NONRSEGS[pdbc][motif][1]) >= 0.5:
            N_majority += 1
            N_proteins[pdbc] = 1
            hitif += 1
        elif NONRSEGS[pdbc][motif][0] > 0 and float(NONRSEGS[pdbc][motif][0])/float(NONRSEGS[pdbc][motif][1]) < 0.5:
            N_minority += 1
            hitif += 1

        if hitif == 0:
            N_none += 1
            hitif += 1

        assert(hitif == 1)

print("#N unique: %5i" % (N_unique))
print("#N 100percent: %5i" % (N_100percent))
print("#N majority: %5i" % (N_majority))
print("#N minority: %5i" % (N_minority))
print("#N none: %5i" % (N_none))
print("#N motifs taken: %5i" % (N_unique+N_100percent+N_majority))
print("#N proteins: %5i" % (len(N_proteins)))
print("\n")

N_onechain = 0
N_onepdb = 0
N_twopdbs = 0
N_morepdbs = 0

totalchains = []
totalpdbs = []
for nonr in N_proteins.keys():

    pdbs = {}
    for pdbC in representative_pdbcs[nonr].keys():
        pdb = pdbC.split('.')[0]
        pdbs[pdb] = 0

    if len(representative_pdbcs[nonr]) == 1:
        N_onechain += 1
    elif len(pdbs) == 1:
        N_onepdb += 1
    elif len(pdbs) == 2:
        N_twopdbs += 1
    elif len(pdbs) >= 3:
        N_morepdbs += 1

    totalchains.append(len(representative_pdbcs[nonr]))
    totalpdbs.append(len(pdbs))
print("#N onechain: %5i" % (N_onechain))
print("#N onepdb, multiple chains: %5i" % (N_onepdb))
print("#N two pdbs: %5i" % (N_twopdbs))
print("#N more than two pdbs: %5i" % (N_morepdbs))
print("# Average number of chains: %8.2f +/- %8.2f" % (np.mean(totalchains), np.std(totalchains)))
print("# Median number of chains: %8.2f" % (np.median(totalchains)))
print("# Average number of pdbs: %8.2f +/- %8.2f" % (np.mean(totalpdbs), np.std(totalpdbs)))
print("# Median number of pdbs: %8.2f" % (np.median(totalpdbs)))

print(np.quantile(totalchains, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]))
print(np.quantile(totalpdbs, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]))
#N unique:   832
#N 100percent:  4559
#N majority:  2390
#N minority:     0
#N none:     0
#N motifs taken:  7781
#N proteins:  5611

exit()
NONRSEGS = pickle.load(open(DERIVED_PATH+"ALL_NONR_PDB_SequenceStats.p4", 'rb'))

N_unique = 0
N_100percent = 0
N_majority = 0
N_minority = 0
N_none = 0
N_proteins = {}
for pdbc in NONRSEGS.keys():
    for motif in NONRSEGS[pdbc].keys():
        hitif = 0
        if NONRSEGS[pdbc][motif][0] == 1 and NONRSEGS[pdbc][motif][1] == 1:
            N_unique += 1
            N_proteins[pdbc] = 1
            hitif += 1
        elif NONRSEGS[pdbc][motif][0] > 1 and NONRSEGS[pdbc][motif][0] == NONRSEGS[pdbc][motif][1]:
            N_100percent += 1
            N_proteins[pdbc] = 1
            hitif += 1
        elif NONRSEGS[pdbc][motif][0] > 0 and float(NONRSEGS[pdbc][motif][0])/float(NONRSEGS[pdbc][motif][1]) >= 0.5:
            N_majority += 1
            N_proteins[pdbc] = 1
            hitif += 1
        elif NONRSEGS[pdbc][motif][0] > 0 and float(NONRSEGS[pdbc][motif][0])/float(NONRSEGS[pdbc][motif][1]) < 0.5:
            N_minority += 1
            hitif += 1

        if hitif == 0:
            N_none += 1
            hitif += 1

        assert(hitif == 1)

print("#N unique: %5i" % (N_unique))
print("#N 100percent: %5i" % (N_100percent))
print("#N majority: %5i" % (N_majority))
print("#N minority: %5i" % (N_minority))
print("#N none: %5i" % (N_none))
print("#N motifs taken: %5i" % (N_unique+N_100percent+N_majority))
print("#N proteins: %5i" % (len(N_proteins)))
#N unique:   832
#N 100percent:  4559
#N majority:  2390
#N minority:  1136
#N none: 10408820
#N motifs taken:  7781
#N proteins:  5611
