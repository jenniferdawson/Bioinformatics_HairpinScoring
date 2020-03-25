import pickle
import numpy as np
from path_setup import DERIVED_PATH, PLOT_PATH

dbdata = pickle.load(open(DERIVED_PATH+"HAIRPINMOTIFS_NONR_PDB_SequenceStats.p4", 'rb'))

aa_order = ['G', 'P', 'T', 'S', 'D', 'N', 'E', 'Q', 'K', 'R', 'A', 'H', 'C', 'W', 'Y', 'F', 'I', 'L', 'V', 'M']

def allin20( seg ):
    for n in range(11):
        if seg[n] not in ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                          'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'):
            return False
    return True


#####
#
# Writes out one motif sequence for every time a hairpin is observed in a cluster of nonredundant pdbs
#
#####

allmotifs = []

ofile_all = open(DERIVED_PATH+"allmotifs.txt", 'w')
ofile_TP = open(DERIVED_PATH+"TPmotifs.txt", 'w')
ofile_SP = open(DERIVED_PATH+"SPmotifs.txt", 'w')
for pdb in dbdata:
    for motif in dbdata[pdb]:
        if allin20(motif):
            ofile_all.write(motif+"\n")
            allmotifs.append(motif)
            MSEG = motif[3]+motif[4]
            if (MSEG == "TP"):
                ofile_TP.write(motif+"\n")
            if (MSEG == "SP"):
                ofile_SP.write(motif+"\n")
ofile_all.close()
ofile_TP.close()
ofile_SP.close()

#####
#
# Writes a grid of amino acid frequencies for each residue position across all nonredundant observations
#
#####

pdata = {}
for n in range(11):
    pdata[n] = {}

for segment in allmotifs:
    for n in range(11):
        aa = segment[n]

        if aa not in pdata[n]:
            pdata[n][aa] = 0
        pdata[n][aa] += 1


minmax_frequencies = np.zeros((20,11))

for x in range(11):
    allfreqs = []
    for aa in pdata[x]:
        freq = 100.0*float(pdata[x][aa])/float(len(allmotifs))
        print(n, aa, pdata[x][aa], 100.0*float(pdata[x][aa])/float(len(allmotifs)))
        allfreqs.append(freq)
    minf = np.min(allfreqs)
    maxf = np.max(allfreqs)

    for y, aa in enumerate(aa_order):
        if aa in pdata[x]:
            freq = 100.0*float(pdata[x][aa])/float(len(allmotifs))
            minmax_frequencies[y][x] = (freq - minf) / (maxf - minf)


import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(3, 4))

ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')

dat = ax.imshow(minmax_frequencies, cmap='viridis', interpolation='nearest', aspect="equal")
cb = fig.colorbar(dat, ax=ax)
cb.set_label("(Frequency at Position - Min.) / (Max. - Min.)")

plt.yticks(range(20), aa_order, rotation='horizontal')
plt.xticks(range(11), [-3,-2,-1,0,1,2,3,4,5,6,7], rotation='horizontal')

plt.xlabel("Motif Position")
plt.ylabel("Amino Acid")

#plt.show()
#plt.savefig(PLOT_PATH+"S2.Heatmap.svg", format="svg")
plt.savefig(PLOT_PATH+"S5.Heatmap.eps", format="eps")