#!/usr/bin/python
import numpy as np
import random
from collections import namedtuple
from PSSM import allin20, make_PSSM, pssm_scoreseg
import pickle
from math import log2
import os
from path_setup import DERIVED_PATH

phosphosite_plus_elevenmers = open(DERIVED_PATH+"phosphosite.11mers.txt").readlines()
hairpin_data = pickle.load(open(DERIVED_PATH+"HAIRPINMOTIFS_NONR_PDB_SequenceStats.p4", 'rb'))

Segdata = namedtuple('segdata', 'score ubi lys degron original_line')
scoremat = make_PSSM(hairpin_data, "segment_probability", set(), set())

pickle.dump(scoremat, open(DERIVED_PATH+"PSSM_Matrix.pickle4", 'wb'), 4)


def get_protein_blocked_data( species, silent=True ):
    """
    Extract species of interest from our full list of phosphosite plus 11mers

    :param species: "human" or "mouse"
    :return: observation data dictionary (D[protein] = [Segdata])
    """

    if not silent:
        motifs_used = {}
        pmotifs_used = {}
        proteins_used = {}

    ProteinBlocks = {}
    used_elevenmer_for_protein = {}
    for i in phosphosite_plus_elevenmers:
        l = i.split()

        l_species = i.split('|')[len(i.split('|'))-2]
        
        if l_species == species:
    
            elevenmer = l[0]
            if allin20(elevenmer.upper()):
                protkey = i.split('|')[len(i.split('|'))-1].split()[0]
    
                score = pssm_scoreseg(elevenmer.upper(), scoremat)
    
                degron = True
                if l[4].startswith("NON"):
                    degron = False
    
                dat = Segdata(score, int(l[1]), int(l[2]), degron, i)
    
                if (protkey in ProteinBlocks) == False:
                    ProteinBlocks[protkey] = []
                    used_elevenmer_for_protein[protkey] = set()
    
                # We don't want to treat multiple instances of a highly repetitive sequence
                # as independent, given that they usually can't be measured independently in mass-spec
                # so we're just using the first observation of each unique elevenmer
                #
                # For human this throws out 201 out of 9965 11mers
                # (over 100 of them coming from a single protein: MUCIN)
                #
                if l[0] not in used_elevenmer_for_protein[protkey]:
                    ProteinBlocks[protkey].append(dat)
                    used_elevenmer_for_protein[protkey].add(l[0])

                    if not silent:
                        motifs_used[l[0]] = 1
                        pmotifs_used[protkey+l[0]] = 1
                        proteins_used[protkey] = 1

    if not silent:
        print("Unique %s motifs used: %10i" % (species, len(motifs_used.keys())))
        print("Unique (to a %s protein) motifs used: %10i" % (species, len(pmotifs_used.keys())))
        print("%s proteins used: %10i" % (species, len(proteins_used.keys())))

    return ProteinBlocks


def protein_list_OEstats(ProteinBlocks, distance, decile_floor, protein_list, dumpdata=""):
    """
    Calculates log2(O/E) values for a list of proteins
    
    for a) observation of 1+ lysine within [distance:int] residues of 11mer
        b) observation of 1+ ubiquitination site within [distance:int] residues of 11mer
        c) observation of a phosphodegron motif overlapping with the [T/S]PXG motif of the 11mer
        
    :param ProteinBlocks:           observation data dictionary (D[protein] = [Segdata])
    :param distance:                measure lys/ubi sites within this distance of the 11mer
    :param decile_floor:            calculate enrichment of this this score window
    :param KEYLIST:                 the list of proteins to use (facilitates bootstrapping)

    :return: log2(OE_lys), log2(OE_ubi), log2(OE_deg)
    """

    # Pull out the list of 11mers to use
    ifile = []
    for k in protein_list:
        for i in ProteinBlocks[k]:
            ifile.append(i)

    # Calculate score range of decile bin
    SL = []
    for i in ifile:
        SL.append(i.score)

    quantile_floor = np.quantile(SL, [decile_floor])[0]

    # Collect statistics
    n_all = 0.0
    n_dbin = 0.0

    lys_all = 0.0
    lys_dbin = 0.0

    ubi_all = 0.0
    ubi_dbin = 0.0

    degron_all = 0.0
    degron_dbin = 0.0

    if dumpdata != "":
        ofile = open(DERIVED_PATH+"S4.Scored."+dumpdata+".txt", 'w')

    for i in ifile:

        SCORE = i.score
        if dumpdata != "":
            ofile.write(str(SCORE)+" "+i.original_line)

        UBI = 0.0
        if i.ubi > 0 and i.ubi <= distance:
            UBI = 1.0

        LYS = 0.0
        if i.lys > 0 and i.lys <= distance:
            LYS = 1.0

        DEGRON = 0.0
        if i.degron:
            DEGRON = 1.0

        if SCORE >= quantile_floor:
            lys_dbin += LYS
            ubi_dbin += UBI
            degron_dbin += DEGRON
            n_dbin += 1.0

        lys_all += LYS
        ubi_all += UBI
        degron_all += DEGRON
        n_all += 1.0

    if dumpdata != "": ofile.close()

    OE_lys = (lys_dbin / n_dbin) / (lys_all / n_all)
    OE_ubi = (ubi_dbin / n_dbin) / (ubi_all / n_all)
    OE_deg = (degron_dbin / n_dbin) / (degron_all / n_all)

    return log2(OE_lys), log2(OE_ubi), log2(OE_deg)

def bootstrap_species(species):
    """

    Runs full bootstrap run against the specified species

    :param species: a string, "human" or "mouse"
    :return: (saves bootstrap files in DERIVED_PATH)
    """

    ProteinBlocks = get_protein_blocked_data(species)

    base_protein_list = list(ProteinBlocks.keys())

    oe_lys_bootdata = []
    oe_ubi_bootdata = []
    oe_deg_bootdata = []
    for B in range(10000):
        if B % 100 == 0:
            print('Running %s Bootstrap # %i' %(species, B+1))

        boot_protein_list = []
        for k in range(len(base_protein_list)):
            boot_protein_list.append(random.sample(base_protein_list, 1)[0])

        oe_lys_dat = []
        oe_ubi_dat = []
        oe_deg_dat = []
        for i, floor in enumerate([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]):
            boot_oe_lys, boot_oe_ubi, boot_oe_deg = protein_list_OEstats(ProteinBlocks=ProteinBlocks,
                                                                         distance=30,
                                                                         decile_floor=floor,
                                                                         protein_list=boot_protein_list)
            oe_lys_dat.append(boot_oe_lys)
            oe_ubi_dat.append(boot_oe_ubi)
            oe_deg_dat.append(boot_oe_deg)

        oe_lys_bootdata.append(oe_lys_dat)
        oe_ubi_bootdata.append(oe_ubi_dat)
        oe_deg_bootdata.append(oe_deg_dat)

        if B % 100 == 0:
            print(oe_lys_dat)
            print(oe_deg_dat)

    oe_lys_bootdata = np.asarray(oe_lys_bootdata)
    pickle.dump(oe_lys_bootdata, open(DERIVED_PATH+species+"_log2OE_lys_bootdata.p4", 'wb'), 4)
    oe_ubi_bootdata = np.asarray(oe_ubi_bootdata)
    pickle.dump(oe_ubi_bootdata, open(DERIVED_PATH+species+"_log2OE_ubi_bootdata.p4", 'wb'), 4)
    oe_deg_bootdata = np.asarray(oe_deg_bootdata)
    pickle.dump(oe_deg_bootdata, open(DERIVED_PATH+species+"_log2OE_deg_bootdata.p4", 'wb'), 4)

def basestats_species(species):
    """
    Collects enrichment details for the specified species

    :param species: a string, "human" or "mouse"
    :return: (saves files in DERIVED_PATH)
    """
    ProteinBlocks = get_protein_blocked_data(species)
    base_protein_list = list(ProteinBlocks.keys())

    dumpdata = species

    oe_lys_dat = []
    oe_ubi_dat = []
    oe_deg_dat = []
    for i, floor in enumerate([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]):
        boot_oe_lys, boot_oe_ubi, boot_oe_deg = protein_list_OEstats(ProteinBlocks=ProteinBlocks,
                                                                     distance=30,
                                                                     decile_floor=floor,
                                                                     protein_list=base_protein_list,
                                                                     dumpdata=dumpdata)
        dumpdata = "" # only dump scores the first time protein_list_OEstats is called

        oe_lys_dat.append(boot_oe_lys)
        oe_ubi_dat.append(boot_oe_ubi)
        oe_deg_dat.append(boot_oe_deg)

    oe_lys_dat = np.asarray(oe_lys_dat)
    pickle.dump(oe_lys_dat, open(DERIVED_PATH+species+"_log2OE_lys_basedata.p4", 'wb'), 4)
    oe_ubi_dat = np.asarray(oe_ubi_dat)
    pickle.dump(oe_ubi_dat, open(DERIVED_PATH+species+"_log2OE_ubi_basedata.p4", 'wb'), 4)
    oe_deg_dat = np.asarray(oe_deg_dat)
    pickle.dump(oe_deg_dat, open(DERIVED_PATH+species+"_log2OE_deg_basedata.p4", 'wb'), 4)


def main():

    if not os.path.exists(DERIVED_PATH+"human_log2OE_lys_bootdata.p4"):
        bootstrap_species("human")
    else:
        print(DERIVED_PATH+"human_log2OE_lys_bootdata.p4 already exists, delete it to rerun")

    if not os.path.exists(DERIVED_PATH+"mouse_log2OE_lys_bootdata.p4"):
        bootstrap_species("mouse")
    else:
        print(DERIVED_PATH+"mouse_log2OE_lys_bootdata.p4 already exists, delete it to rerun")

    basestats_species("human")
    basestats_species("mouse")

    get_protein_blocked_data( "human", silent=False )
    get_protein_blocked_data( "mouse", silent=False )

if __name__ == '__main__':
    main()