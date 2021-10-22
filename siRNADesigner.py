import argparse
import numpy as np

parser = argparse.ArgumentParser()
filePathParser = parser.add_mutually_exclusive_group(required=True)
filePathParser.add_argument("-a", "--accessionNumber",
                            help="Gene accession number in the NCBI nucleotide DB")
filePathParser.add_argument("-s", "--geneSequence",
                            help="DNA sequence of the target gene")
parser.add_argument("-o", "--outputPrefix", default="",
                    help="Output file prefix path")
args = parser.parse_args()


accessionNumber = args.accessionNumber
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
