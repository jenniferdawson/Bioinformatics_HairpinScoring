import math
import copy

def allin20( seg ):
    for n in range(11):
        if seg[n] not in ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                          'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'):
            return False
    return True

def total( dict ):
    tot = 0
    for AA in dict:
        tot += dict[AA]
    return float(tot)


def make_PSSM( hairpin_data, normtype, leave_out_pdbc_list, leave_out_elevenmer_list ):

    all_hairpin_position_counts = {}
    for AA in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
        all_hairpin_position_counts[AA] = 0.1

    hairpin_position_counts = []
    for n in range(11):
        hairpin_position_counts.append(copy.deepcopy(all_hairpin_position_counts))

    for pdbc in hairpin_data:
        if pdbc not in leave_out_pdbc_list:                     # Remove test set pdb chains
            for elevenmer in hairpin_data[pdbc]:
                if elevenmer not in leave_out_elevenmer_list:   # Remove test set 11mers
                    if allin20(elevenmer):                      # Only canonical amino acids
                        for n in range(11):
                            AA = elevenmer[n]
                            all_hairpin_position_counts[AA] += 1
                            hairpin_position_counts[n][AA] += 1

    scores = []
    for n in range(11):
        scores.append({})

        for AA in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
            position_prob = hairpin_position_counts[n][AA] / total(hairpin_position_counts[n])

            if normtype == "segment_probability":
                amino_prob = all_hairpin_position_counts[AA] / total(all_hairpin_position_counts)

            if normtype == "equal_probability":
                amino_prob = 0.05

            l2prob = math.log( position_prob / amino_prob , 2)
            scores[n][AA] = l2prob

    return scores

def pssm_scoreseg( seg, scoremat ):
    score = 0.0
    for n in range(11):
        score += scoremat[n][seg[n]]
    return score