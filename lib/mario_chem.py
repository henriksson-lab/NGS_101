"""Custom oligos printed in the defining MARIO paper."""
from chemdraw import Segment, feature, oligo


RNA_LINKER_LEFT = "CUAG"
RNA_LINKER_RIGHT = "AGCCCAUGCAAUGCGAGGA"
LINKER_CDNA = "TCGCATTGCATGGGCTACTAGCAT"
RT_ADAPTER = "AGATCGGAAGAGCGGTTCAG"
RT_BODY = "AGATCGGAAGAGCGTCGTGGATCCTGAACCGCTCTTCCGATCT"
CUT_OLIGO = "GTTCAGGATCCACGACGCTCTTCAAAA"
DP5 = "CACGACGCTCTTCCGATCT"
DP3 = "CTGAACCGCTCTTCCGATCT"
PE1 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
PE2 = "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"


def _rt_primer(name, index, trailing_n):
    key = name.lower().replace("-", "_").replace(" ", "_").replace("/", "_")
    umi = feature(key + "_umi", "umi", "random")
    sample = feature(key + "_index", "sample_index",
                     "whitelist", whitelist="three MARIO RT-primer sample barcodes")
    return oligo(name, [
        Segment("random barcode part 1", "NN", "umi", feature=umi),
        Segment("sample barcode", index, "cbc", feature=sample),
        Segment("random barcode part 2", "N" * trailing_n, "umi", feature=umi),
        Segment("RT-primer body including BamHI", RT_BODY),
    ], five="5′-/5Phos/-", three="-3′")


def oligos():
    return (
        oligo("biotinylated RNA linker", [
            Segment("RNA-linker left arm", RNA_LINKER_LEFT),
            Segment("internal biotinylated dT", "/iBiodT/", placeholder=True,
                    length_bp=1),
            Segment("RNA-linker right arm", RNA_LINKER_RIGHT),
        ]),
        oligo("cDNA strand with RNA linker", [
            Segment("four phosphorothioate-protected bases", "T*C*G*C*",
                    placeholder=True, length_bp=4),
            Segment("linker-complement arm", LINKER_CDNA[4:]),
        ]),
        oligo("pre-adenylated RT adapter", [Segment("RT adapter", RT_ADAPTER)],
              five="5′-/5rApp/-", three="-/3ddC/-3′"),
        _rt_primer("RT primer ES-1", "AGGT", 3),
        _rt_primer("RT primer ES-2 / MEF", "CGCC", 4),
        _rt_primer("RT primer ES-indirect", "CATT", 4),
        oligo("Cut_oligo", [Segment("BamHI cleavage oligo", CUT_OLIGO)],
              three="-/3InvdT/-3′"),
        oligo("truncated PCR forward primer DP5", [Segment("DP5", DP5)]),
        oligo("truncated PCR reverse primer DP3", [Segment("DP3", DP3)]),
        oligo("Illumina PE PCR forward primer 1.0", [Segment("PE1", PE1)]),
        oligo("Illumina PE PCR reverse primer 2.0", [Segment("PE2", PE2)]),
    )
