import argparse
import numpy as np
from Bio import Entrez
import re

parser = argparse.ArgumentParser()
filePathParser = parser.add_mutually_exclusive_group(required=True)
filePathParser.add_argument("-a", "--accessionNumber",
                            help="Gene accession number in the NCBI nucleotide DB")
filePathParser.add_argument("-s", "--geneSequence",
                            help="DNA sequence of the target gene")
parser.add_argument("-e", "--entrezEmail",
                    help="Required when accessionNumber is provided, In case of\
                    excessive usage of the E-utilities, NCBI will attempt to contact\
                    a user at the email address provided before blocking access to the\
                    E-utilities")
parser.add_argument("-o", "--outputPrefix", default="",
                    help="Output file prefix path")
args = parser.parse_args()


accessionNumber = args.accessionNumber
entrezEmail = args.entrezEmail
geneSequence = args.geneSequence
outputPrefix = args.outputPrefix


class basicAligner():

    def __init__(self):
        self.enodeDict = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
        self.scoringMatrix = np.array([[2, -1, -1, -1],
                                       [-1, 2, -1, -1],
                                       [-1, -1, 2, -1],
                                       [-1, -1, -1, 2]])
        self.mismatches = None
        self.seq1 = None
        self.seq2 = None
        self.score = None
        self.pairs = None
        self.best_match = None

    def align(self, seq1, seq2):
        self.mismatches = 0
        self.seq1 = self.str2np_arr(seq1.upper())
        self.seq2 = self.str2np_arr(seq2.upper())
        self.pairs = self.find_alignment()
        self.best_match = [
            ''.join(self._pairs2sequences(self.pairs)[0]), self.mismatches]
        return self.__str__()

    def str2np_arr(self, seq):
        np_arr = []
        for nucleotide in seq:
            if nucleotide == self.enodeDict[0]:
                np_arr.append(0)
            elif nucleotide == self.enodeDict[1]:
                np_arr.append(1)
            elif nucleotide == self.enodeDict[2]:
                np_arr.append(2)
            elif nucleotide == self.enodeDict[3]:
                np_arr.append(3)
        return np.array(np_arr, dtype=np.uint8)

    def find_alignment(self):
        self.score = np.zeros(
            (self.seq1.size+1, self.seq2.size+1), dtype=np.int16)
        self._compute_alignmentMatrix()
        return self._trackback()

    def _compute_alignmentMatrix(self):
        for i in range(1, self.seq1.size+1):
            for j in range(1, self.seq2.size+1):
                self.score[i, j] = self.score[i-1, j-1] + self._get_score(i, j)

    def _get_score(self, i, j):
        return self.scoringMatrix[self.seq1[i-1], self.seq2[j-1]]

    def _get_aligned_pair(self, i, j):
        if i > 0:
            self.nucl1 = self.enodeDict[self.seq1[i-1]]
        else:
            self.nucl1 = '_'

        if j > 0:
            self.nucl2 = self.enodeDict[self.seq2[j-1]]
        else:
            self.nucl2 = '_'

        return (self.nucl1, self.nucl2)

    def _trackback(self):
        alignmentPairs = []
        self.max_ind = np.where(self.score[-1, :] == np.max(self.score[-1, :]))
        i = self.seq1.size
        j = self.max_ind[0][0]
        while self.score[i, j] > 0:
            if self._get_score(i, j) == -1:
                self.mismatches += 1
            alignmentPairs.append(self._get_aligned_pair(i, j))
            i -= 1
            j -= 1
        while i > 0:
            alignmentPairs.append(self._get_aligned_pair(i, 0))
            i -= 1
        alignmentPairs.reverse()
        return alignmentPairs

    def _pairs2sequences(self, pairs):
        top_seq = []
        bottom_seq = []
        for (b, t) in pairs:
            bottom_seq.append(b)
            top_seq.append(t)
        return [top_seq, bottom_seq]

    def __str__(self):
        if self.pairs is None:
            return ""
        return str(self._pairs2sequences(self.pairs)[0]) + '\n' + str(self._pairs2sequences(self.pairs)[1]) + '\n' + 'Number of mismatches= '+str(self.mismatches)

    def get_best(self):
        return self.best_match


def get_exons(Acc):
    handle = Entrez.efetch(db="nucleotide", id=Acc, retmode="xml")
    records = Entrez.read(handle)
    str_records = str(records[0])
    exons = []

    while str_records.find("'GBFeature_key': 'exon'") != -1:
        exons.append([])
        exon_indx = str_records.find("'GBFeature_key': 'exon'")
        narrow = ''
        for i in range(20):
            narrow += str_records[exon_indx+47+i]
        narrow = narrow.split(',')[0][1:-1]
        min1 = int(narrow[0:narrow.find('.')])
        exons[-1].append(min1)
        max2 = int(narrow[narrow.find('.')+2::])
        exons[-1].append(max2)
        str_records = str_records[0:str_records.find(
            "'GBFeature_key': 'exon'")]+'0'+str_records[str_records.find("'GBFeature_key': 'exon'")+1::]

    exons_seq = []
    for i in exons:
        exons_seq.append(records[0]['GBSeq_sequence'][i[0]-1:i[1]-1])

    exons_seq = ''.join(exons_seq)

    return exons_seq


if not (accessionNumber is None):
    if not (entrezEmail is None):
        Entrez.email = entrezEmail
    gene = get_exons(accessionNumber).upper()
    gene = re.sub(r"\s+", "", gene)
else:
    gene = geneSequence.upper()
minStart = 0
maxEnd = len(gene)-23
candidates = dict()
# Phase 1:
# Test each 23 nucleotide for some features
for i in range(minStart, maxEnd):
    candidate = gene[i:i+23]
    Gcount = candidate.count('Gcount')
    Ccount = candidate.count('Ccount')
    GCcount = Gcount + Ccount
    GCPercent = (GCcount/23)*100
    if GCPercent > 32 and GCPercent < 55:  # 1st filter (GCPercent)
        shortRepeatsFree = True
        for j in range(23):  # 2nd filter (internal repeats)
            if candidate[j:j+5] in candidate[j+5:23]:
                shortRepeatsFree = False
                break
        if shortRepeatsFree:
            # 3rd filter (GCPercent stretches)
            match = re.match(r'GC{10}', candidate)
            if match:
                continue
            else:
                if candidate[20] in 'AT':  # 4th filter 5' end of guide A/U
                    if candidate[2] in 'CG':  # 5th filter 5' end of passenger G/C
                        Acount = candidate[13:21].count('A')
                        Ucount = candidate[13:21].count('T')
                        AUcount = Acount + Ucount
                        # 6th filter at least 4 A/U residues in 5' guide 7bp
                        if AUcount >= 4:
                            # 7th filter No G at position 13 of passenger
                            if candidate[14] in 'ATC':
                                # 8th filter  A/U at position 19 in passenger
                                if candidate[20] in 'AT':
                                    # 9th filter Gcount/Ccount at position 19 in guide
                                    if candidate[2] in 'GCPercent':
                                        candidates[candidate] = [
                                            i+1, i+23, AUcount]

# phase 2:
compDict = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
for cand in candidates:
    guide = ''
    passenger = ''
    for nucleotide in range(20, -1, -1):
        guide = guide + cand[nucleotide]
    candidates[cand].append(guide)
    for nucleotide in range(2, 23):
        if cand[nucleotide] == 'T':
            passenger = passenger + 'U'
        else:
            passenger = passenger + cand[nucleotide]
    candidates[cand].append(passenger)
    if candidates[cand][3][0] == 'U':
        candidates[cand][2] += 2
    if candidates[cand][3][0] == 'A':
        candidates[cand][2] += 1
    match = re.match(r'[A,U]{1,2}', candidates[cand][4][1:5])
    if match:
        candidates[cand][2] += 1
    else:
        candidates[cand][2] += -1
