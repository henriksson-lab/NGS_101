"""Source-exact BLISS break-tag adapter chemistry.

Yan et al. 2017 Supplementary Data 1 prints eleven classic adapters.  A1 is the
representative molecule here; Supplementary Fig. 1a establishes which end ligates and
that the bottom oligo is 5-prime phosphorylated at that blunt end.
"""
from chemdraw import Scene, Segment, feature, oligo, revcomp

T7_END = "CGATTGAGGCCGGTAATACGACTCACTATAGGGG"
RA5 = "GTTCAGAGTTCTACAGTCCGACGATC"
A1_BARCODE = "CATCACGC"
A1_TOP = T7_END + RA5 + "N" * 8 + A1_BARCODE
A1_BOTTOM = revcomp(A1_BARCODE) + "N" * 8 + revcomp(RA5) + revcomp(T7_END) + "AA"


def adapter_a1() -> Scene:
    """Checked A1 duplex, including its ligating and concatemer-blocking ends."""
    barcode = feature("sample_barcode", "sample_index", "whitelist",
                      whitelist="BLISS Supplementary Data 1 adapter set A1-A11")
    umi = feature("umi", "umi", "random")
    top = [
        Segment("T7 promoter-bearing end", T7_END, "t7"),
        Segment("RA5", RA5, "r1"),
        Segment("8-nt UMI", "N" * 8, "umi", feature=umi),
        Segment("A1 sample barcode", A1_BARCODE, "cbc", feature=barcode),
    ]
    bottom = [
        Segment("A1 sample barcode complement", revcomp(A1_BARCODE), "cbc",
                feature=barcode),
        Segment("8-nt UMI complement", "N" * 8, "umi", feature=umi),
        Segment("RA5 complement", revcomp(RA5), "r1"),
        Segment("T7 end complement", revcomp(T7_END), "t7"),
        Segment("3-prime blocking overhang", "AA"),
    ]
    sc = Scene()
    sc.strand("top", top, label="BLISS A1 top")
    sc.anneal("bottom", bottom, to="top",
              pair=("A1 sample barcode complement", "A1 sample barcode"),
              label="BLISS A1 bottom", unpaired=("3-prime blocking overhang",),
              mod5="P")
    sc.ligation_end("top", "A1 sample barcode", kind="blunt",
                    owner="adapter", text="blunt break ligates here")
    sc.mark("bottom", "3-prime blocking overhang",
            "3′-AA suppresses head-to-tail adapter concatemers")
    sc.labels("top")
    return sc


def oligos():
    return (
        oligo("BLISS A1 top", [Segment("A1 top", A1_TOP)]),
        oligo("BLISS A1 bottom", [Segment("A1 bottom", A1_BOTTOM)],
              five="5′-/5Phos/-", three="-3′"),
    )
