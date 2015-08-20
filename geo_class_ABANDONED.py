__author__ = 'Matteo'
__doc__='''The geo_shift.py methods were opening and closing files, which was getting irksome. Here an object is passed.
It is list of seqRecords. with the exception that the iteration is across features.
'''

from Bio import SeqIO
from collections import Sequence


class Genome(Sequence):  #got sick of double loops for chr and gene
    def __init__(self, filepath):
        self.data=list(SeqIO.parse(open(filepath, "rU"), "genbank"))
    def __getitem__(self, index):
        return self.data[index]
    def __len__(self):
        return len(self.data)
    def __iter__(self):
        return iter([gene for element in self.data for gene in element.features])


if __name__ == "__main__":
    #ref = Genome("Eco NC_000913.gb")


