from path_setup import DERIVED_PATH
import os
import pickle
from PSSM import make_PSSM, pssm_scoreseg, allin20

from sys import argv


if os.path.exists(DERIVED_PATH+"PSSM_Matrix.pickle4"):
    scoremat = pickle.load(open(DERIVED_PATH+"PSSM_Matrix.pickle4", 'rb'))
else:

    if os.path.exists(DERIVED_PATH+"HAIRPINMOTIFS_NONR_PDB_SequenceStats.p4"):
        hairpin_data = pickle.load(open(DERIVED_PATH+"HAIRPINMOTIFS_NONR_PDB_SequenceStats.p4", 'rb'))

        scoremat = make_PSSM(hairpin_data, "segment_probability", set(), set())

        pickle.dump(scoremat, open(DERIVED_PATH+"PSSM_Matrix.pickle4", 'wb'), 4)
    else:

        print("\n\nERROR: Cannot find %s\nERROR: Cannot find %s\n\n"%
              (DERIVED_PATH+"PSSM_Matrix.pickle4",
               DERIVED_PATH+"HAIRPINMOTIFS_NONR_PDB_SequenceStats.p4"))

        print("Either the derived databath path has not been set or the files have been deleted\n\n")

        print('First: Set paths in "./src/path_setup.py"\n')
        print('If that does not work:\n\nGenerate "'+DERIVED_PATH+'HAIRPINMOTIFS_NONR_PDB_SequenceStats.p4"')
        print('by running "./src/Step4.make_database.py"\n\n')
        exit()

if len(argv) != 2:
    print("\n\nUsage:\n\tpython ./src/score11mer.py PSKQSNNKYAA")
    print("\tor\n\tpython ./src/score11mer.py ./derived/allmotifs.txt\n\n")

    print("This is the PSSM that will be used: (position counts i=0 at [T/S] in XXX[T/S]PXGXXXX\n\n")

    ostr = "Position "
    for AA in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
        ostr += "%6s " % (AA)
    print(ostr)
    for pos in range(len(scoremat)):
        ostr = "%8i " % (pos-3)
        for AA in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
            ostr += "%6.2f " % (scoremat[pos][AA])
        print(ostr)
    print("\n\n")
    exit()

seqs = []
if os.path.exists(argv[1]):
    ifile = open(argv[1]).readlines()

    for i in ifile:
        seq = i.split()[0].upper()
        if len(seq) == 11 and allin20(seq):
            seqs.append(seq)
else:
    seq = argv[1].upper()
    if len(seq) == 11 and allin20(seq):
        seqs.append(seq)

if len(seqs) == 0:
    print("\n\nError: Sequence must be 11 residues long and consist of the 20 canonical amino acids")
    exit()

#
# All scoring happens in these two lines
#
for seq in seqs:
    print("%s %10.6f" % (seq, pssm_scoreseg(seq, scoremat)))