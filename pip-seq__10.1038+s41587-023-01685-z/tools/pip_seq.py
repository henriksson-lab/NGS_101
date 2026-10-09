"""Molecular construct model for the paper-described PIP-seq workflow."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature

CELL_CASSETTE = feature("cell_barcode_cassette", "cell_barcode", "combinatorial",
                        note="four barcode blocks; internal boundaries unpublished")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"ligated barcode/linker cassette": CELL_CASSETTE, "UMI": UMI_FEATURE,
            "i7 index read": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


BEAD_PREFIX = "T" * 7
BEAD_JUNCTION = "ACGACTC"
BEAD_READ1_PART = il.TRUSEQ_READ1[3:28]
P5_SPACER = "GCCTGTCCGCGG"
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
P5_LIBRARY = il.P5 + P5_SPACER + rt.SMART_HANDLE + "AC"
I7_INDEX_READ = "TAAGGCGA"


def inferred_barcode_cassette() -> Segment:
    """Unpublished split-pool block: keep it visibly unresolved for every caller."""
    return seg("ligated barcode/linker cassette", "X" * 36, "cbc",
               placeholder=True, inferred=True,
               note="The paper states four ligated barcodes but does not print their blocks or splints.")


def inferred_read1_completion() -> Segment:
    """Unprinted ligation bases needed to complete the named Read-1 primer site."""
    return seg("Read-1-site completion", il.TRUSEQ_READ1[28:], "r1", inferred=True,
               note="These five bases are absent from the printed bead primer.")


def bead_oligo() -> Construct:
    return Construct([
        seg("bead spacer", BEAD_PREFIX), seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("bead junction", BEAD_JUNCTION),
        seg("partial TruSeq Read 1", BEAD_READ1_PART, "r1"),
        inferred_read1_completion(), inferred_barcode_cassette(),
        seg("UMI", "N" * 12, "umi", placeholder=True),
        seg("poly(T)", "T" * 19), seg("V anchor", "V", placeholder=True),
    ], name="PIP-seq bead oligo")


def tso_oligo() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrGrG", "rGrGrG", placeholder=True)]


def rt_scene(insert_nt: int = 28) -> Scene:
    mrna = [seg("transcript", "X" * (insert_nt - 1), placeholder=True),
            seg("anchor partner", "V", placeholder=True), seg("poly(A)", "A" * 19)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead", list(bead_oligo()), to="mRNA", pair=("poly(T)", "poly(A)"),
              label="released bead primer")
    sc.arrow("bead", "bulk reverse transcription and template switching")
    return sc


SEQ_PRIMERS = (
    sp.custom("Read 1", "TruSeq Read 1", il.TRUSEQ_READ1,
              "Illumina; the library omits the site's first three bases",
              "Its 5' three bases form an unpaired flap; the 3' 30 nt land exactly."),
    sp.NEXTERA["I1"], sp.NEXTERA["R2"],
)
RUN_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")


def final_library(insert_nt: int = 28) -> Construct:
    bead = list(bead_oligo())
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("P5 custom spacer", P5_SPACER),
        *bead[1:], seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 index read", I7_INDEX_READ, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="paper-described PIP-seq library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid PIP-seq library: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")
