import argparse
import numpy as np
from Bio import Entrez, SeqIO
import re
import math

parser = argparse.ArgumentParser()
inputParser = parser.add_mutually_exclusive_group(required=True)
inputParser.add_argument("-a", "--accessionNumber",
                         help="Gene accession number in the NCBI nucleotide DB")
inputParser.add_argument("-s", "--geneSequence",
                         help="DNA sequence of the target gene")
parser.add_argument("-e", "--entrezEmail", default="example@gmail.com",
                    help="optional when accessionNumber is provided, In case of\
                    excessive usage of the E-utilities, NCBI will attempt to contact\
                    a user at the email address provided before blocking access to the\
                    E-utilities")
offTargerParser = parser.add_argument_group()
parser.add_argument("-t", "--reduceOffTargets", default=0,
                    help="Reduce off targets based on the melting temperatures of both\
                        the guide and the passenger RNA. In addition, it ensure that\
                        neither has a match in the provided transcriptome")
parser.add_argument("-f", "--transcriptomePath",
                    help="fasta file, optional when reduceOffTargets equals 1")
parser.add_argument("-o", "--outputPrefix", default="",
                    help="Output file prefix path")
args = parser.parse_args()


accessionNumber = args.accessionNumber
entrezEmail = args.entrezEmail
geneSequence = args.geneSequence
outputPrefix = args.outputPrefix if args.outputPrefix != "" else "siRNAcandidates"
reduceOffTargets = bool(int(args.reduceOffTargets))
transcriptomePath = args.transcriptomePath


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


def meltingTemp(seq, compSeq, Na=0.1, CT=1e-4, A=-10.8, R=1.987):
    # SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
    nearest_neighbour = {
        'init': (0.2, -5.7), 'init_A/T': (2.2, 6.9), 'init_G/C': (0, 0),
        'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
        'sym': (0, -1.4),
        'AA/TT': (-7.6, -21.3), 'AT/TA': (-7.2, -20.4), 'TA/AT': (-7.2, -20.4),
        'CA/GT': (-8.5, -22.7), 'GT/CA': (-8.4, -22.4), 'CT/GA': (-7.8, -21.0),
        'GA/CT': (-8.2, -22.2), 'CG/GC': (-10.6, -27.2), 'GC/CG': (-9.8, -24.4),
        'GG/CC': (-8.0, -19.0)}

    # Watkins & SantaLucia (2005), Nucl Acids Res 33: 6258-6267
    mis_match = {
        'AG/TT': (1.0, 0.9), 'AT/TG': (-2.5, -8.3), 'CG/GT': (-4.1, -11.7),
        'CT/GG': (-2.8, -8.0), 'GG/CT': (3.3, 10.4), 'GG/TT': (5.8, 16.3),
        'GT/CG': (-4.4, -12.3), 'GT/TG': (4.1, 9.5), 'TG/AT': (-0.1, -1.7),
        'TG/GT': (-1.4, -6.2), 'TT/AG': (-1.3, -5.3), 'AA/TG': (-0.6, -2.3),
        'AG/TA': (-0.7, -2.3), 'CA/GG': (-0.7, -2.3), 'CG/GA': (-4.0, -13.2),
        'GA/CG': (-0.6, -1.0), 'GG/CA': (0.5, 3.2), 'TA/AG': (0.7, 0.7),
        'TG/AA': (3.0, 7.4),
        'AC/TT': (0.7, 0.2), 'AT/TC': (-1.2, -6.2), 'CC/GT': (-0.8, -4.5),
        'CT/GC': (-1.5, -6.1), 'GC/CT': (2.3, 5.4), 'GT/CC': (5.2, 13.5),
        'TC/AT': (1.2, 0.7), 'TT/AC': (1.0, 0.7),
        'AA/TC': (2.3, 4.6), 'AC/TA': (5.3, 14.6), 'CA/GC': (1.9, 3.7),
        'CC/GA': (0.6, -0.6), 'GA/CC': (5.2, 14.2), 'GC/CA': (-0.7, -3.8),
        'TA/AC': (3.4, 8.0), 'TC/AA': (7.6, 20.2),
        'AA/TA': (1.2, 1.7), 'CA/GA': (-0.9, -4.2), 'GA/CA': (-2.9, -9.8),
        'TA/AA': (4.7, 12.9), 'AC/TC': (0.0, -4.4), 'CC/GC': (-1.5, -7.2),
        'GC/CC': (3.6, 8.9), 'TC/AC': (6.1, 16.4), 'AG/TG': (-3.1, -9.5),
        'CG/GG': (-4.9, -15.3), 'GG/CG': (-6.0, -15.8), 'TG/AG': (1.6, 3.6),
        'AT/TT': (-2.7, -10.8), 'CT/GT': (-5.0, -15.8), 'GT/CT': (-2.2, -8.4),
        'TT/AT': (0.2, -1.5),
        'AI/TC': (-8.9, -25.5), 'TI/AC': (-5.9, -17.4), 'AC/TI': (-8.8, -25.4),
        'TC/AI': (-4.9, -13.9), 'CI/GC': (-5.4, -13.7), 'GI/CC': (-6.8, -19.1),
        'CC/GI': (-8.3, -23.8), 'GC/CI': (-5.0, -12.6),
        'AI/TA': (-8.3, -25.0), 'TI/AA': (-3.4, -11.2), 'AA/TI': (-0.7, -2.6),
        'TA/AI': (-1.3, -4.6), 'CI/GA': (2.6, 8.9), 'GI/CA': (-7.8, -21.1),
        'CA/GI': (-7.0, -20.0), 'GA/CI': (-7.6, -20.2),
        'AI/TT': (0.49, -0.7), 'TI/AT': (-6.5, -22.0), 'AT/TI': (-5.6, -18.7),
        'TT/AI': (-0.8, -4.3), 'CI/GT': (-1.0, -2.4), 'GI/CT': (-3.5, -10.6),
        'CT/GI': (0.1, -1.0), 'GT/CI': (-4.3, -12.1),
        'AI/TG': (-4.9, -15.8), 'TI/AG': (-1.9, -8.5), 'AG/TI': (0.1, -1.8),
        'TG/AI': (1.0, 1.0), 'CI/GG': (7.1, 21.3), 'GI/CG': (-1.1, -3.2),
        'CG/GI': (5.8, 16.9), 'GG/CI': (-7.6, -22.0),
        'AI/TI': (-3.3, -11.9), 'TI/AI': (0.1, -2.3), 'CI/GI': (1.3, 3.0),
        'GI/CI': (-0.5, -1.3)}

    # SantaLucia & Peyret (2001) Patent Application WO 01/94611
    terminal_mm = {
        'AA/TA': (-3.1, -7.8), 'TA/AA': (-2.5, -6.3), 'CA/GA': (-4.3, -10.7),
        'GA/CA': (-8.0, -22.5),
        'AC/TC': (-0.1, 0.5), 'TC/AC': (-0.7, -1.3), ' CC/GC': (-2.1, -5.1),
        'GC/CC': (-3.9, -10.6),
        'AG/TG': (-1.1, -2.1), 'TG/AG': (-1.1, -2.7), 'CG/GG': (-3.8, -9.5),
        'GG/CG': (-0.7, -19.2),
        'AT/TT': (-2.4, -6.5), 'TT/AT': (-3.2, -8.9), 'CT/GT': (-6.1, -16.9),
        'GT/CT': (-7.4, -21.2),
        'AA/TC': (-1.6, -4.0), 'AC/TA': (-1.8, -3.8), 'CA/GC': (-2.6, -5.9),
        'CC/GA': (-2.7, -6.0), 'GA/CC': (-5.0, -13.8), 'GC/CA': (-3.2, -7.1),
        'TA/AC': (-2.3, -5.9), 'TC/AA': (-2.7, -7.0),
        'AC/TT': (-0.9, -1.7), 'AT/TC': (-2.3, -6.3), 'CC/GT': (-3.2, -8.0),
        'CT/GC': (-3.9, -10.6), 'GC/CT': (-4.9, -13.5), 'GT/CC': (-3.0, -7.8),
        'TC/AT': (-2.5, -6.3), 'TT/AC': (-0.7, -1.2),
        'AA/TG': (-1.9, -4.4), 'AG/TA': (-2.5, -5.9), 'CA/GG': (-3.9, -9.6),
        'CG/GA': (-6.0, -15.5), 'GA/CG': (-4.3, -11.1), ' GG/CA': (-4.6, -11.4),
        'TA/AG': (-2.0, -4.7), 'TG/AA': (-2.4, -5.8),
        'AG/TT': (-3.2, -8.7), 'AT/TG': (-3.5, -9.4), 'CG/GT': (-3.8, -9.0),
        'CT/GG': (-6.6, -18.7), 'GG/CT': (-5.7, -15.9), 'GT/CG': (-5.9, -16.1),
        'TG/AT': (-3.9, -10.5), 'TT/AG': (-3.6, -9.8)}

    tmpSeq = seq
    tmpCompSeq = compSeq
    deltaH = 0
    deltaS = 0
    d_h = 0  # Names for indexes
    d_s = 1  # 0 and 1

    left_tmm = tmpCompSeq[:2][::-1] + '/' + tmpSeq[:2][::-1]
    if left_tmm in terminal_mm:
        deltaH += terminal_mm[left_tmm][d_h]
        deltaS += terminal_mm[left_tmm][d_s]
        tmpSeq = tmpSeq[1:]
        tmpCompSeq = tmpCompSeq[1:]

    right_tmm = tmpSeq[-2:] + '/' + tmpCompSeq[-2:]
    if right_tmm in terminal_mm:
        deltaH += terminal_mm[right_tmm][d_h]
        deltaS += terminal_mm[right_tmm][d_s]
        tmpSeq = tmpSeq[:-1]
        tmpCompSeq = tmpCompSeq[:-1]

    deltaH += nearest_neighbour['init_oneG/C'][d_h]
    deltaS += nearest_neighbour['init_oneG/C'][d_s]

    # Finally, the 'zipping'
    for basenumber in range(len(tmpSeq) - 1):
        neighbors = tmpSeq[basenumber:basenumber + 2] + \
            '/' + tmpCompSeq[basenumber:basenumber + 2]
        if neighbors in mis_match:
            deltaH += mis_match[neighbors][d_h]
            deltaS += mis_match[neighbors][d_s]
        elif neighbors[::-1] in mis_match:
            deltaH += mis_match[neighbors[::-1]][d_h]
            deltaS += mis_match[neighbors[::-1]][d_s]
        elif neighbors in nearest_neighbour:
            deltaH += nearest_neighbour[neighbors][d_h]
            deltaS += nearest_neighbour[neighbors][d_s]
        elif neighbors[::-1] in nearest_neighbour:
            deltaH += nearest_neighbour[neighbors[::-1]][d_h]
            deltaS += nearest_neighbour[neighbors[::-1]][d_s]

    Tm = ((1000 * deltaH)/(A + deltaS + R * math.log(CT / 4))) - \
        273.15 + 16.6 * math.log10(Na)
    return Tm


if not (accessionNumber is None):
    if not (entrezEmail is None):
        Entrez.email = entrezEmail
    gene = get_exons(accessionNumber).upper()
    gene = re.sub(r"\s+", "", gene)
else:
    gene = geneSequence.upper()
minStart = 0
maxEnd = len(gene)-22
candidates = dict()
# Phase 1:
# Test each 23 nucleotide for some features
for i in range(minStart, maxEnd):
    candidate = gene[i:i+23]
    if candidate in candidates:
        continue
    Gcount = candidate.count('G')
    Ccount = candidate.count('C')
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

compDict = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
alignment = basicAligner()
success = dict()
for cand in candidates:
    guide = ''
    passenger = ''
    for nucleotide in range(20, -1, -1):
        guide = guide + compDict[cand[nucleotide]]
    candidates[cand].append(guide)
    for nucleotide in range(2, 23):
        if cand[nucleotide] == 'T':
            passenger = passenger + 'U'
        else:
            passenger = passenger + cand[nucleotide]
    candidates[cand].append(passenger)
    # phase 2:
    if candidates[cand][3][0] == 'U':
        candidates[cand][2] += 2
    if candidates[cand][3][0] == 'A':
        candidates[cand][2] += 1
    match = re.match(r'[A,U]{1,2}', candidates[cand][4][1:5])
    if match:
        candidates[cand][2] += 1
    else:
        candidates[cand][2] += -1
    # phase 3:
    if candidates[cand][4][0:2] == 'AA':
        candidates[cand][2] += 1
    if candidates[cand][4][2] == 'A':
        candidates[cand][2] += 1
    if candidates[cand][4][9] == 'U':
        candidates[cand][2] += 1
    # phase 4:
    seqGuide = ''
    compSeqGuide = ''
    for nucleotide in candidates[cand][3]:
        if nucleotide == 'U':
            seqGuide += 'T'
            compSeqGuide += 'A'
        else:
            seqGuide += nucleotide
            compSeqGuide += compDict[nucleotide]

    seqPassenger = ''
    compSeqPassenger = ''
    for nucleotide in candidates[cand][4]:
        if nucleotide == 'U':
            seqPassenger += 'T'
            compSeqPassenger += 'A'
        else:
            seqPassenger += nucleotide
            compSeqPassenger += compDict[nucleotide]
    if reduceOffTargets:
        TmGuide = meltingTemp(seqGuide[1:8], compSeqGuide[1:8])
        TmPass = meltingTemp(seqPassenger[1:8], compSeqPassenger[1:8])
        if TmPass >= 21.5 or TmGuide >= 21.5:
            continue
        candidates[cand].append(TmGuide)
        candidates[cand].append(TmPass)
        if not (transcriptomePath is None):
            passengerMisM = 2
            guideMisM = 2
            valid = True
            fasta_sequences = SeqIO.parse(transcriptomePath, 'fasta')
            for fasta in fasta_sequences:
                seq = str(fasta.seq)
                if "*" in seq:
                    seq = seq.replace("*", "")
                alignment.align(compSeqGuide, seq)
                guideMisM = min(alignment.mismatches, guideMisM)
                if guideMisM < 2:
                    valid = False
                    break
                alignment.align(compSeqPassenger, seq)
                passengerMisM = min(alignment.mismatches, passengerMisM)
                if passengerMisM < 2:
                    valid = False
                    break
            if valid:
                candidates[cand].append(guideMisM)
                candidates[cand].append(passengerMisM)
    success[cand] = candidates[cand]

resultsFile = open(outputPrefix+".csv", 'w')
for cand in success:
    resultsFile.write(",".join([str(success[cand][0]), str(success[cand][1]),
                                cand, *[str(x) for x in success[cand][2:]]]) + '\n')
