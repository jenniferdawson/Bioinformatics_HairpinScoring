def read_dssp( Dfile ):
    """
    Reads a DSSP file and returns a dictionary, format = ddict[chain][res] = [aa, ss, sslong, HN2O, O2HN]
    :param Dfile: full filesystem path to a DSSP format file
    :return: dictionary ddict[chain][res] = [aa, ss, sslong, HN2O, O2HN]
                                            aa = amino acid
                                            ss = 1 letter secondary structure code: H, G, E, T or L (for all others)
                                            sslong = full DSSP assignment (columns [16:25])
                                            HN2O = HN hydrogen bonds
                                            O2HN = O hydrogen bonds
    """
    dssp = open(Dfile).readlines()
    ddict = {}

    started = False
    for d in dssp:
        dl = d.split()
        if len(dl) > 5:
            if started == False:
                if dl[0] == "#":
                    started = True
            else:
                if dl[1] != "!*" and dl[1] != "!":
                    res = int(d[6:10])
                    chain = d[11:12]
                    aa = d[13:14]
                    ss = d[16:17]

                    HN2O = {}
                    O2HN = {}

                    if float(d[42:50].split(',')[1].split()[0]) <= -1.0:
                        HN2Oa = d[42:50].split(',')[0].split()[0]
                        HN2O[HN2Oa] = 0

                    if float(d[64:72].split(',')[1].split()[0]) <= -1.0:
                        HN2Ob = d[64:72].split(',')[0].split()[0]
                        HN2O[HN2Ob] = 0

                    if float(d[52:61].split(',')[1].split()[0]) <= -1.0:
                        O2HNa = d[52:61].split(',')[0].split()[0]
                        O2HN[O2HNa] = 0

                    if float(d[75:83].split(',')[1].split()[0]) <= -1.0:
                        O2HNb = d[75:83].split(',')[0].split()[0]
                        O2HN[O2HNb] = 0

                    if ss != "G" and ss != "H" and ss != "E" and ss != "T":
                        ss = "L"

                    sslong = d[16:25]

                    if chain not in ddict:
                        ddict[chain] = {}
                    ddict[chain][res] = [aa, ss, sslong, HN2O, O2HN]
    return ddict

def CHECK_HBONDS(dssp, c, d):
    """
    Checks to see if the hydrogen bonding pattern matches a TPGGT hairpin

    dssp: a dssp dictionary (ddict[chain][res] = [aa, ss, sslong, HN2O, O2HN])
    c: chain for sequence motif
    d: starting residue number for 11 residue motif being inspected
    return: True if all hydrogen bonds match, false if one or more don't
    """

    if ("8" in dssp[c][d+1][3]) == False: return False
    if ("8" in dssp[c][d+1][4]) == False: return False
    if ("-8" in dssp[c][d+9][3]) == False: return False
    if ("-8" in dssp[c][d+9][4]) == False: return False

    if ("4" in dssp[c][d+3][3]) == False: return False
    if ("3" in dssp[c][d+3][4]) == False: return False

    if ("-3" in dssp[c][d+6][3]) == False: return False
    if ("-4" in dssp[c][d+7][4]) == False: return False

    return True

def SS_SEGMENT( dssp, chain, res, tlen ):
    """
    Grabs a tlen length secondary structure sequence from the DSSP dictionary

    Sequences use 1 letter per residue, with: "G", "H", "E", & "T" kept from dssp and "L" used for all others

    dssp: a dssp dictionary (ddict[chain][res] = [aa, ss, sslong, HN2O, O2HN])
    chain: chain to use
    res: starting residue number
    tlen: length of motif
    return: a tlen length sequence if one exists, "" if any of the residues requested do not exist
    """
    hasN = True
    sseq = ""
    for iN in range(tlen):
        if (res+iN) not in dssp[chain]:
            hasN = False

    if hasN == True:
        for iN in range(tlen):
            sseq += dssp[chain][res+iN][1]
    return sseq

def AA_SEGMENT( dssp, chain, res, tlen ):
    """
    Grabs a tlen length amino acid sequence from the DSSP dictionary
    dssp: a dssp dictionary (ddict[chain][res] = [aa, ss, sslong, HN2O, O2HN])
    chain: chain to use
    res: starting residue number
    tlen: length of motif
    return: a tlen length sequence if one exists, "" if any of the residues requested do not exist
    """
    hasN = True
    sseq = ""
    for iN in range(tlen):
        if (res+iN) not in dssp[chain]:
            hasN = False

    if hasN == True:
        for iN in range(tlen):
            sseq += dssp[chain][res+iN][0]
    return sseq.upper()
