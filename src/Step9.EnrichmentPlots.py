import matplotlib.pyplot as plt
import pickle
from path_setup import PLOT_PATH, DERIVED_PATH

human_lys = pickle.load(open(DERIVED_PATH+"human_log2OE_lys_basedata.p4", 'rb'))
human_ubi = pickle.load(open(DERIVED_PATH+"human_log2OE_ubi_basedata.p4", 'rb'))
human_deg = pickle.load(open(DERIVED_PATH+"human_log2OE_deg_basedata.p4", 'rb'))

mouse_lys = pickle.load(open(DERIVED_PATH+"mouse_log2OE_lys_basedata.p4", 'rb'))
mouse_ubi = pickle.load(open(DERIVED_PATH+"mouse_log2OE_ubi_basedata.p4", 'rb'))
mouse_deg = pickle.load(open(DERIVED_PATH+"mouse_log2OE_deg_basedata.p4", 'rb'))

human_lys_boot = pickle.load(open(DERIVED_PATH+"human_log2OE_lys_bootdata.p4", 'rb'))
human_ubi_boot = pickle.load(open(DERIVED_PATH+"human_log2OE_ubi_bootdata.p4", 'rb'))
human_deg_boot = pickle.load(open(DERIVED_PATH+"human_log2OE_deg_bootdata.p4", 'rb'))

mouse_lys_boot = pickle.load(open(DERIVED_PATH+"mouse_log2OE_lys_bootdata.p4", 'rb'))
mouse_ubi_boot = pickle.load(open(DERIVED_PATH+"mouse_log2OE_ubi_bootdata.p4", 'rb'))
mouse_deg_boot = pickle.load(open(DERIVED_PATH+"mouse_log2OE_deg_bootdata.p4", 'rb'))

fig, ax = plt.subplots(3,2,figsize=(7, 8))
fig.tight_layout(pad=2.0, h_pad=4.0, w_pad=3.0)

for x in range(3):
    for y in range(2):
        ax[x][y].set_xlabel("Fraction by highest score (%)")
        ax[x][y].hlines(0,0,11,linestyles='dashed',color='gray')

x = [1,2,3,4,5,6,7,8,9,10]
ax[0][0].set_yticks([0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14])
ax[0][0].boxplot(human_lys_boot, showfliers=False, whis=[0.5,99.5])
ax[0][0].scatter(x, human_lys, s=8)
ax[0][0].set_ylabel("log2(O/E)")
ax[0][0].set_title("Human: Lysine")

ax[1][0].set_yticks([0.0, 0.5, 1.0, 1.5])
ax[1][0].boxplot(human_deg_boot, showfliers=False, whis=[0.5,99.5])
ax[1][0].scatter(x, human_deg, s=8)
ax[1][0].set_ylabel("log2(O/E)")
ax[1][0].set_title("Human: Phosphodegron")

ax[2][0].set_yticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax[2][0].boxplot(human_ubi_boot, showfliers=False, whis=[0.5,99.5])
ax[2][0].scatter(x, human_ubi, s=8)
ax[2][0].set_ylabel("log2(O/E)")
ax[2][0].set_title("Human: Ubiquitination")

ax[0][1].set_yticks([0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14])
ax[0][1].boxplot(mouse_lys_boot, showfliers=False, whis=[0.5,99.5])
ax[0][1].scatter(x, mouse_lys, s=8)
ax[0][1].set_ylabel("log2(O/E)")
ax[0][1].set_title("Mouse: Lysine")

ax[1][1].set_yticks([0.0, 0.5, 1.0, 1.5])
ax[1][1].boxplot(mouse_deg_boot, showfliers=False, whis=[0.5,99.5])
ax[1][1].scatter(x, mouse_deg, s=8)
ax[1][1].set_ylabel("log2(O/E)")
ax[1][1].set_title("Mouse: Phosphodegron")

ax[2][1].set_yticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax[2][1].boxplot(mouse_ubi_boot, showfliers=False, whis=[0.5,99.5])
ax[2][1].scatter(x, mouse_ubi, s=8)
ax[2][1].set_ylabel("log2(O/E)")
ax[2][1].set_title("Mouse: Ubiquitination")

for x in range(3):
    for y in range(2):
        ax[x][y].set_xticklabels([100, 90, 80, 70, 60, 50, 40, 30, 20, 10])

#plt.savefig(PLOT_PATH+"S5.Enrichment.svg", format="svg")
plt.savefig(PLOT_PATH+"S9.Enrichment.eps", format="eps")
#plt.show()
