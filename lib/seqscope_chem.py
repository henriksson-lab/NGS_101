"""Oligos printed in the defining Seq-Scope paper."""
from chemdraw import Segment, feature, oligo


P5 = "CAAGCAGAAGACGGCATACGAGAT"
TR1 = "TCTTTCCCTACACGACGCTCTTCCGATCT"
HDMI20 = "NNVNNVNNVNNVNNVNNNNN"
HDMI32 = "NNVNBVNNVNNVNNVNNVNNVNNVNNVNNNNN"
HR1 = "TCTTGTGACTACAGCACCCTCGACTCTCGC"
POLYT = "T" * 27
DRAI = "TTTAAA"
DRAI_ADAPTER = "GACTTTCACCAGTCCATGAT"
P7 = "GTGTAGATCTCGGTGGTCGCCGTATCATT"
READ1_DRAI = "ATCATGGACTGGTGAAAGTCTTTAAA" + "A" * 27 + "GCGAGAGTCGAGGGTGCTGTAGTCACAAGA"
TR2_RANDOM = "TCAGACGTGTGCTCTTCCGATCT" + "N" * 9

HDMI_ECORI = "HNNBNBNBNBNBNBNBNNNN"
HR1B = "CCCGTTCGCAACATGTCTGGCGTCATA"
ECORI = "GAATTC"
ECORI_ADAPTER = "CGCAGTCCAG"
READ1_ECORI = "CTGGACTGCGGAATTCTATGACGCCAGACATGTTGCGAACGGG"
UMI_OLIGO = "A" * 30 + "N" * 8 + "C" + "TATGACGCCAGACATGTTGCGAACGGG"


def _hdmi(identifier, bases):
    return Segment("HDMI spatial barcode", bases, "cbc", placeholder=True,
                   length_bp=len(bases),
                   feature=feature(identifier, "spatial_barcode", "random"))


def _seed(name, identifier, bases, hr):
    return oligo(name, [
        Segment("P5-side sequence", P5), Segment("TruSeq Read 1", TR1),
        _hdmi(identifier, bases), Segment("HDMI read arm", hr),
        Segment("oligo-dT capture", POLYT), Segment("DraI site", DRAI),
        Segment("DraI adapter", DRAI_ADAPTER), Segment("P7-side sequence", P7),
    ])


def oligos():
    sample_index = feature("seqscope_i7", "sample_index", "unknown")
    return (
        _seed("HDMI-DraI", "seqscope_hdmi20", HDMI20, HR1),
        _seed("HDMI32-DraI", "seqscope_hdmi32", HDMI32, HR1),
        oligo("Read1-DraI custom sequencing primer",
              [Segment("Read1-DraI", READ1_DRAI)]),
        oligo("TruSeq Read2-conjugated random primer", [
            Segment("TruSeq Read 2", TR2_RANDOM[:-9]),
            Segment("9-nt random UMI", TR2_RANDOM[-9:], "umi",
                    feature=feature("seqscope_random_primer_umi", "umi", "random")),
        ]),
        oligo("first-round PCR forward", [
            Segment("forward-primer body", "TCTTTCCCTACACGACGC"),
            Segment("phosphorothioate-protected 3-prime bases", "*T*C",
                    placeholder=True, length_bp=2),
        ]),
        oligo("first-round PCR reverse", [
            Segment("reverse-primer body", "TCAGACGTGTGCTCTTCC"),
            Segment("phosphorothioate-protected 3-prime bases", "*G*A",
                    placeholder=True, length_bp=2),
        ]),
        oligo("second-round PCR forward", [
            Segment("full forward-primer body", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCT"),
            Segment("phosphorothioate-protected 3-prime bases", "*T*C",
                    placeholder=True, length_bp=2),
        ]),
        oligo("second-round PCR reverse", [
            Segment("P7-side primer body", "CAAGCAGAAGACGGCATACGAGAT"),
            Segment("8-mer sample index", "N" * 8, "cbc", placeholder=True,
                    feature=sample_index),
            Segment("Read-2-side primer body", "GTGACTGGAGTTCAGACGTGTGCTCTTCC"),
            Segment("phosphorothioate-protected 3-prime bases", "*G*A",
                    placeholder=True, length_bp=2),
        ]),
        oligo("HDMI-TruEcoRI", [
            Segment("P5-side sequence", P5), Segment("TruSeq Read 1", TR1),
            _hdmi("seqscope_hdmi_ecori", HDMI_ECORI),
            Segment("HDMI read arm B", HR1B), Segment("EcoRI site", ECORI),
            Segment("EcoRI adapter", ECORI_ADAPTER), Segment("P7-side sequence", P7),
        ]),
        oligo("Read1-EcoRI custom sequencing primer",
              [Segment("Read1-EcoRI", READ1_ECORI)]),
        oligo("UMI-oligo for the EcoRI array", [
            Segment("oligo-dA", "A" * 30),
            Segment("8-nt UMI", "N" * 8, "umi",
                    feature=feature("seqscope_ecori_umi", "umi", "random")),
            Segment("junction C", "C"), Segment("HR1B-binding arm", UMI_OLIGO[-28:]),
        ]),
    )
