import argparse
import numpy as np
from Bio import Entrez, SeqIO
import re
import math

parser = argparse.ArgumentParser()
inputParser = parser.add_mutually_exclusive_group(required=True)
inputParser.add_argument("-a", "--accessionNumber",
                         help="Gene accession number in the NCBI nucleotide DB.\
                             Multiple accession numbers should be separated by ','")
inputParser.add_argument("-s", "--geneSequence",
                         help="DNA sequence of the target gene.\
                             Multiple genes sequences should be separated by ','")
inputParser.add_argument("-i", "--targetFasta",
                         help="Fasta file of one target or more")
parser.add_argument("-e", "--entrezEmail",
                    help="optional when accessionNumber is provided, In case of\
                    excessive usage of the E-utilities, NCBI will attempt to contact\
                    a user at the email address provided before blocking access to the\
                    E-utilities")
parser.add_argument("-t", "--reduceOffTargets", default=0,
                    help="Reduce off targets based on the melting temperatures of both\
                        the guide and the passenger RNA. In addition, it ensure that\
                        neither has a match in the provided transcriptome")
parser.add_argument("-f", "--transcriptomePath",
                    help="fasta file, optional when reduceOffTargets equals 1")
parser.add_argument("-o", "--outputPrefix", default="",
                    help="Output file prefix path")
args = parser.parse_args()


accessionNumber = args.accessionNumber if args.accessionNumber is None else args.accessionNumber.split(
    ",")
entrezEmail = args.entrezEmail
geneSequence = args.geneSequence if args.geneSequence is None else args.geneSequence.split(
    ",")
targetFasta = args.targetFasta
outputPrefix = [args.outputPrefix] if args.outputPrefix != "" else [
    "siRNAcandidates"]
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


def meltingTemp(seq):
    return (seq.count('A')+seq.count('T')) * 2 + (seq.count('G')+seq.count('C')) * 4


if not (accessionNumber is None):
    if not (entrezEmail is None):
        Entrez.email = entrezEmail
    gene = [get_exons(x).upper() for x in accessionNumber]
    gene = [re.sub(r"\s+", "", x) for x in gene]
    if len(gene) != 1:
        outputPrefix = [outputPrefix[0]+"_"+str(x+1) for x in range(len(gene))]
elif not (targetFasta is None):
    geneFasta = SeqIO.parse(targetFasta, 'fasta')
    geneFastaUnpacked = [(x.seq.upper(), x.id) for x in geneFasta]
    gene = [x[0] for x in geneFastaUnpacked]
    outputPrefix = [outputPrefix[0]+"_"+x[1] for x in geneFastaUnpacked]
    del geneFastaUnpacked
else:
    gene = [x.upper() for x in geneSequence]
    if len(gene) != 1:
        outputPrefix = [outputPrefix[0]+"_"+str(x+1) for x in range(len(gene))]
for geneIndx in range(len(gene)):
    minStart = 0
    maxEnd = len(gene[geneIndx])-22
    candidates = dict()
    # Phase 1:
    # Test each 23 nucleotide for some features
    for i in range(minStart, maxEnd):
        candidate = gene[geneIndx][i:i+24]
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
                # 3rd filter (GC stretches)
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
                                        # 9th filter G/C at position 19 in guide
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
            TmGuide = meltingTemp(seqGuide[1:8])
            TmPass = meltingTemp(seqPassenger[1:8])
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

    resultsFile = open(outputPrefix[geneIndx]+".csv", 'w')
    for cand in success:
        resultsFile.write(",".join([str(success[cand][0]), str(success[cand][1]),
                                    cand, *[str(x) for x in success[cand][2:]]]) + '\n')
