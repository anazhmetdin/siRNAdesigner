import argparse

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
