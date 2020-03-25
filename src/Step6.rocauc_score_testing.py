import pickle
import random
import sklearn
from sklearn import metrics
import numpy as np
from PSSM import allin20, make_PSSM, pssm_scoreseg
import matplotlib.pyplot as plt
from path_setup import DERIVED_PATH, PLOT_PATH

hairpin_data = pickle.load(open(DERIVED_PATH+"HAIRPINMOTIFS_NONR_PDB_SequenceStats.p4", 'rb'))
STPXG_data = pickle.load(open(DERIVED_PATH+"STPXG_NONR_PDB_SequenceStats.p4", 'rb'))

def train_and_test( tag, normtype, leave_out_pdbc_list, leave_out_elevenmer_list ):

    scoremat = make_PSSM(hairpin_data, normtype, leave_out_pdbc_list, leave_out_elevenmer_list)

    test_scores = []
    test_truth = []
    train_scores = []
    train_truth = []

    for pdbc in STPXG_data: # all XXX[TS]PXGXXXX elevenmers, with or without hairpins
        leave_out = False

        for elevenmer in STPXG_data[pdbc].keys():
            if allin20(elevenmer):
                MSEG = elevenmer[3]+elevenmer[4]
                MSEG += "X"
                MSEG += elevenmer[6]
                assert( MSEG in ("SPXG","TPXG") )

                if pdbc in leave_out_pdbc_list:
                    leave_out = True
                if elevenmer in leave_out_elevenmer_list:
                    leave_out = True

                truth_value = 0.0
                if float(STPXG_data[pdbc][elevenmer][0])/(STPXG_data[pdbc][elevenmer][1]) >= 0.5:
                    truth_value = 1.0 # at least half of the 11mers for this pdb cluster are hairpins

                score = pssm_scoreseg(elevenmer, scoremat)

                if leave_out:
                    test_scores.append(score)
                    test_truth.append(truth_value)
                else:
                    train_scores.append(score)
                    train_truth.append(truth_value)


    test_scores = np.asarray(test_scores)
    test_truth = np.asarray(test_truth)
    train_scores = np.asarray(train_scores)
    train_truth = np.asarray(train_truth)

    auc_train = sklearn.metrics.roc_auc_score(train_truth, train_scores)

    auc_test = sklearn.metrics.roc_auc_score(test_truth, test_scores)
    fpr_test, tpr_test, thresholds = sklearn.metrics.roc_curve(test_truth, test_scores)

    return fpr_test, tpr_test, auc_test, auc_train

fig, ax = plt.subplots(figsize=(3, 3))
plt.subplots_adjust(left=0.2, bottom=0.2, right=0.95, top=0.95)
plt.xlabel("False positive rate")
plt.ylabel("True positive rate")
plt.xlim([0,1])
plt.ylim([0,1])
ax.plot([0,1],[0,1], transform=ax.transAxes, color="black")

auc_tests = []
auc_trains = []
for rseed in range(100):
    print("Running 5-fold cross validation with random seed %i of 100" %(rseed+1))
    random.seed(rseed)
    pdbc_list = list(STPXG_data.keys())
    random.shuffle(pdbc_list)
    leave_out_sets = np.array_split(pdbc_list, 5)
    for i in leave_out_sets:
        fpr_test, tpr_test, auc_test, auc_train = train_and_test("PDBC"+str(rseed), "segment_probability", set(i), set())
        plt.plot(fpr_test, tpr_test, linewidth=2, alpha=0.02, color="black")

        auc_tests.append(auc_test)
        auc_trains.append(auc_train)

plt.text(0.5, 0.1, " Average AUC\n%5.2f +/- %5.2f" % (np.mean(auc_tests), np.std(auc_tests)))
plt.savefig(PLOT_PATH+"S6.ROCAUC.svg", format="svg")
#plt.show()
