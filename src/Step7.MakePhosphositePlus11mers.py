from path_setup import DB_PATH, DERIVED_PATH

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

aas = []
for aa in longer_names:
    aas.append(longer_names[aa])


"""

Read fasta format file containing phosphosite annotation data

(lowercase K, S & T denote ubiquitination and phosphorylation sites)

"""
inputseqs = open(DB_PATH+"Phosphosite_seq.K_ST_annotated.fasta").readlines()

FULLNAME = {}
PROTSEQS = {}
onfasta = ""
seq = ""
for i in inputseqs:
    if len(i.split()) > 0:
        if i[0] == ">":
            if onfasta != "":
                PROTSEQS[onfasta] = seq
            seq = ""
            onfasta = i.split()[0]
        else:
            seq += i.split()[0]
if onfasta != "":
    PROTSEQS[onfasta] = seq

"""

Now iterate through all 11mers looking for sequences that match XXX[T/S]PXGXXXX

For each one write a line to "" containing

11mer sequence
distance to closest ubiquitination site (0 if inside 11mer)
distance to closest lysine (has to be outside of 11mer)
tag for phosphorylation status of [T/S] residue (PHOSPHO or NONPHOS)
tag for overlap with phosphodegron (NONGRON, DEGRON1, DEGRON2, and DEGRON3)
protein name from PhosphositePlus

"""

ofile = open(DERIVED_PATH+"phosphosite.11mers.txt", 'wt')

def distance( pn, n ):
    if pn < n: return n - pn
    if pn > n+11: return pn - (n+11)
    return 0
    
for P in PROTSEQS.keys():
    PSEQ = PROTSEQS[P]
    for n in range(len(PSEQ)-11):
        SEG = PSEQ[n:n+11]
        USEG = SEG.upper()

        if USEG[4] == "P" and USEG[6] == "G":
            if USEG[3] == "T" or USEG[3] == "S":
                
                TAG = "NONPHOS"
                if SEG[3] == "t" or SEG[3] == "s":
                    TAG = "PHOSPHO"

                minlys = -1
                minub = -1
                for pn in range(len(PSEQ)):

                    dist = distance(pn, n)
                    
                    if PSEQ[pn] == "k":
                        if minub == -1 or dist < minub:
                            minub = dist
                    if PSEQ[pn] == "k" or PSEQ[pn] == "K":
                        if dist > 0:
                            if minlys == -1 or dist < minlys:
                                minlys = dist

                DTAG = "NONGRON"
                if USEG[3] == "T":
                    LIVMP = {"L":0,"I":0,"V":0,"M":0,"P":0}
                    if (USEG[0] in LIVMP) or (USEG[1] in LIVMP) or (USEG[2] in LIVMP):
                        if USEG[7] == "S" or USEG[7] == "T":
                            DTAG = "DEGRON1"
                        if USEG[7] == "E":
                            DTAG = "DEGRON2"

                    if USEG[1] == "D" or USEG[1] == "E":
                        if USEG[6] == "K":
                            STAG = "DEGRON3"

                ofile.write("%s %i %i %s %s %s\n" %(SEG, minub, minlys, TAG, DTAG, P))

ofile.close()