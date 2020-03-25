
# Step 1: Installation

In base directory of repo (~/HairpinScoring/)

```
python3 -m venv env
source ./env/bin/activate
pip install -r requirements.in
```

(works with python 3.7.3 and pip version 18.0)


To run the standalone hairpin PSSM score function run:

`python ./src/score11mer.py [11mer sequence]`

Running `python ./src/score11mer.py` without arguments will
print out the PSSM being used for scoring.


# Step 2: Set filesystem paths

`./src/path_setup.py` contains path information

Set the paths to the HairpinScoring repo
ie:
```
BASE_PATH = "/home/user/HairpinScoring/"
DSSP_PATH = "/home/user/DSSP_FILES/"
```

# Step 3 (Optional): Download DSSP data

The data that gets derived from the DSSP files is included in the repo, but if you want to recalculate it
then you'll need to download the set of files we used (listed in "./database/dssp.files") from the DSSP database.

`rsync -avz rsync://rsync.cmbi.umcn.nl/dssp/ ./DSSP_DATA/`
is one way to do this, but you may need to use their FTP site or reach out directly.

To test how the code itself works a small number of DSSP files are provided in DSSP_TEST
(just set DSSP_PATH = "/home/user/HairpinScoring/database/DSSP_TEST/" in ./src/path_setup.py)

(but this will overwrite the derived database included in the repo)

# Step 4 (Optional): Generate 11mer databases from DSSP data

`python ./src/Step4.make_database.py`

These databases have been included in the repo, and are in the "derived" folder:
HAIRPINMOTIFS_NONR_PDB_SequenceStats.p4
STPXG_NONR_PDB_SequenceStats.p4


# Step 5: Hairpin frequency statistics

`python ./src/Step5.frequency_statistics.py`

This will generate 3 text files in "derived" folder for use with "https://weblogo.berkeley.edu/logo.cgi"
```
allmotifs.txt
TPmotifs.txt
SPmotifs.txt
```

And will make a frequency plot in the "plots" folder
S5.Heatmap.eps


# Step 6: ROCAUC Testing

`python ./src/Step6.rocauc_score_testing.py`

will run 5 fold cross validation 100 times and save a summary plot "S6.ROCAUC.svg"


# Step 7: Extract 11mers from PhosphositePlus

`python ./src/Step7.MakePhosphositePlus11mers.py`

will pull out 11mers from PhosphositePlus and annotate them for closest ubiquitination site, closest lysine,
and overlap with phosphodegron motifs

# Step 8: Bootstrap

`python ./src/Step8.decile_with_boot.py`

Calculated enrichment values over 10000 random samples of the pdb list used (with replacement) and saves the data
for bootstrap analysis.

Files are saved as:
```
./derived/human_log2OE_deg_basedata.p4
./derived/human_log2OE_ubi_basedata.p4
./derived/mouse_log2OE_lys_basedata.p4
./derived/human_log2OE_deg_bootdata.p4
./derived/human_log2OE_ubi_bootdata.p4
./derived/mouse_log2OE_lys_bootdata.p4
./derived/human_log2OE_lys_basedata.p4
./derived/mouse_log2OE_deg_basedata.p4
./derived/mouse_log2OE_ubi_basedata.p4
./derived/human_log2OE_lys_bootdata.p4
./derived/mouse_log2OE_deg_bootdata.p4
./derived/mouse_log2OE_ubi_bootdata.p4
```

# Step 8: Plot Enrichment

python ./src/Step9.EnrichmentPlots.py

Calculates bootstrap enrichment distributions and makes "./plots/S9.Enrichment.eps"