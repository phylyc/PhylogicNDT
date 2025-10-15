# from enum import Enum, unique #requires >2.7.6 or 3.x
import os
import pandas as pd


def Enum(**enums):
    return type('Enum', (), enums)


Cluster = Enum(**dict([('C1', 'green'), ('C2', 'blue'), ('C3', 'red')]))

cn_state_whitelist = frozenset({(1., 2.), (0., 2.), (2., 2.)})


MutStatus = Enum(OK="OK",
                 REMOVED="REMOVED",  # blacklisted
                 GRAYLIST="GRAYLIST")  # not used in clustering

MutType = Enum(INS="INS",
               DEL="DEL",
               SNV="SNV")


class Genome(object):
    def __init__(self, build="hg19", prefix=""):
        self.build = build
        prefix = "chr" if build == "hg38" else ""
        self.CHROMS = tuple(prefix + str(c) for c in list(range(1, 23)) + ['X', 'Y'])
        self.ARMS = ("p", "q")
        self.CSIZE = {
            "hg19": [
                249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
                135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                59128983, 63025520, 48129895, 51304566, 156040895, 57227415
            ],
            "hg38": [
                248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
                133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                58617616, 64444167, 46709983, 50818468, 156040895, 57227415
            ]
        }.get(build, [])
        self.CHROM_DICT = {c: s for c, s in zip(self.CHROMS, self.CSIZE)}
        self.CENTS = {
            "hg19": [
                125000000, 93300000, 91000000, 50400000, 48400000, 61000000, 59900000, 45600000,
                49000000, 40200000, 53700000, 35800000, 17900000, 17600000, 19000000, 36600000,
                24000000, 17200000, 26500000, 27500000, 13200000, 14700000, 60600000, 12500000
            ],
            "hg38": [
                123400000, 93900000, 90900000, 50000000, 48800000, 59800000, 60100000, 45200000,
                43000000, 39800000, 53400000, 35500000, 17700000, 17200000, 19000000, 36800000,
                25100000, 18500000, 26200000, 28100000, 12000000, 15000000, 61000000, 10400000
            ]
        }.get(build, [])
        self.CENT_DICT = {c: s for c, s in zip(self.CHROMS, self.CENTS)}

        file = os.path.dirname(__file__) + '/supplement_data/' + {
            "hg19": "cytoBand.hg19.txt",
            "hg38": "cytoBand.hg38.txt"
        }.get(build, "WRONG_BUILD_FILE")
        self.cytoband_table = (
            pd.read_csv(file, sep="\t", header=None, names=["chromosome", "start", "end", "band", "stain"])
            .astype({"chromosome": str, "start": int, "end": int, "band": str})
        )
