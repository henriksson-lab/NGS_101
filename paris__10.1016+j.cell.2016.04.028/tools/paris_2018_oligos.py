"""Published2018 comparator inventory only; not2016 molecular graph defaults.
Source: Lu2018 MethodsMolBiol section2.2 originaltext154–190.
"""
import illumina as il
from chemdraw import Segment, feature, oligo

RNA_ADAPTER = "AGATCGGAAGAGCGGTTCAG"
RT_RANDOM_W = "W" * 3
RT_RANDOM_N3 = "N" * 3
RT_INDEX = "ATCACG"
RT_RANDOM_N5 = "N" * 5
RT_BODY = "TACCCTTCGCTTCACACACAAG"
RT_SPACER_1 = "GGATCC"
RT_SPACER_2 = "TACTGAACCGC"
P3_TALL_V4 = "GGCATTCCTGCTGAACCGCTCTTCCGATCT"
P6_TALL_V4 = "CTCTTTCCCCTTGTGTGTGAAGCGAAGGGT"
# P3-Solexa is Illumina's PE PCR Primer 2.0; P3Tall v4 is a 3' part of the same
# historical PE Read 2 primer.
P3_SOLEXA = il.PE_PCR_PRIMER_2
if not il.PE_READ2_PRIMER.endswith(P3_TALL_V4):
    raise AssertionError("P3Tall v4 is not a 3' part of the PE Read 2 primer")
P6_SOLEXA = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
P6_CUSTOM_SEQ = "CACTCTTTCCCCTTGTGTGTGAAGCGAAGGGTA"


def oligos():
    return (
        oligo("preadenylated 3′ RNA adapter",
              [Segment("3-prime RNA adapter", RNA_ADAPTER)],
              five="5′-/5rApp/-", three="-/3Biotin/-3′"),
        oligo("RT-primer-1 representative (ATCACG)", [
            Segment("A/T molecular tag", RT_RANDOM_W, "umi", placeholder=True,
                    length_bp=3, feature=feature("rt_w3", "umi", "random")),
            Segment("random molecular tag part 1", RT_RANDOM_N3, "umi",
                    placeholder=True, length_bp=3,
                    feature=feature("rt_n8", "umi", "random", group="rt_n8", part=1)),
            Segment("sample barcode", RT_INDEX, "cbc",
                    feature=feature("rt_index", "sample_index", "whitelist",
                                    whitelist="24 PARIS RT indices")),
            Segment("random molecular tag part 2", RT_RANDOM_N5, "umi",
                    placeholder=True, length_bp=5,
                    feature=feature("rt_n8", "umi", "random", group="rt_n8", part=2)),
            Segment("RT primer body", RT_BODY),
            Segment("internal spacer 1", "[iSp18]", placeholder=True),
            Segment("spacer-flanked BamHI tract", RT_SPACER_1),
            Segment("internal spacer 2", "[iSp18]", placeholder=True),
            Segment("3-prime adapter-binding arm", RT_SPACER_2),
        ], five="5′-/5Phos/-", three="-3′"),
        oligo("P3Tall v4", [Segment("short P3 PCR primer", P3_TALL_V4)]),
        oligo("P6Tall v4", [Segment("short P6 PCR primer", P6_TALL_V4)]),
        oligo("P3-Solexa", [Segment("P7-side flow-cell PCR primer", P3_SOLEXA)]),
        oligo("P6-Solexa", [Segment("P5-side flow-cell PCR primer", P6_SOLEXA)]),
        oligo("P6_custom_SeqPrimer",
              [Segment("custom sequencing primer", P6_CUSTOM_SEQ)]),
    )
