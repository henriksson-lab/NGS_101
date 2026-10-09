"""Li et al. 2017 long-read ChIA-PET molecular model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from batch_ngs import nextera_library, seg
from chemdraw import Construct, Row, Scene, Segment
from proximity_linker import ThreePrimeBridgeLinker

TITLE = "Long-read ChIA-PET — bridge-linked chromatin contacts"
NOTES = "01_long-read-chia-pet.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/nprot.2017.012">'
          'Li et al., <i>Nature Protocols</i> (2017)</a>.')
SUMMARY = ("Protein-associated chromatin fragments are repaired and dA-tailed, joined in "
           "proximity through a double-ended biotin bridge, then tagmented; streptavidin "
           "retains only fragments that still contain the contact junction.")
CAVEAT = ("The bridge-linker bases and modification are published. The final sequencing "
          "ends use the canonical Nextera components supplied by the named Illumina kit.")

BRIDGE_F = "CGCGATATCTTATCTGACT"
BRIDGE_R = "GTCAGATAAGATATCGCGT"
BRIDGE = ThreePrimeBridgeLinker("biotin bridge linker", BRIDGE_F, BRIDGE_R,
                                overhang="T", modified_top_position=9)

def contact_segments() -> list[Segment]:
    core = BRIDGE.top_core
    return [
        seg("left contact tag", "X" * 28, placeholder=True),
        seg("left A:T ligation junction", "A", bottom="T"),
        seg("bridge core left", core[:9], "cbc"),
        seg("internal biotin-dT", core[9], "umi"),
        seg("bridge core right", core[10:], "cbc"),
        seg("right T:A ligation junction", "T", bottom="A"),
        seg("right contact tag", "X" * 28, placeholder=True),
    ]

CONTACT = Construct(contact_segments(), name="DNA–bridge–DNA proximity product")
FINAL_LIBRARY, SEQ_PRIMERS = nextera_library(contact_segments(),
                                              "long-read ChIA-PET library")
FINAL_CAPTION = ("A bridge-containing proximity product after Nextera tagmentation, "
                 "streptavidin selection and indexed PCR.")
SEQUENCING_INTRO = ("Paired reads enter from the two Tn5 ends. A read reaching either "
                    "published bridge-linker orientation proves that the selected fragment "
                    "contains a proximity-ligation junction.")
TAGMENTATION = (55, 5)
MAX_PCR_CYCLES = 13


def contact_scene() -> Scene:
    sc = Scene.duplex(list(CONTACT), label="proximity product")
    sc.junction("top", "left contact tag", "left A:T ligation junction", "ligation")
    sc.junction("top", "right T:A ligation junction", "right contact tag", "ligation")
    sc.mark("top", "internal biotin-dT", "streptavidin handle")
    return sc


def sections():
    return [
        ("Enrich protein-associated chromatin", [
            Row(chunks=[("fragment A — target protein / antibody — fragment B", None, False)])
        ], "Dual-cross-linked chromatin is sonicated and immunoprecipitated with an "
           "antibody to the protein or histone mark of interest."),
        ("Anneal the double-ended bridge linker", BRIDGE.scene().rows(),
         "Both oligos are 5′ phosphorylated. Their reverse-complementary 18-nt cores leave "
         "one 3′ T overhang at each end; the highlighted internal dT carries biotin."),
        ("Join two dA-tailed chromatin fragments", contact_scene().rows(),
         "Each terminal T pairs with a genomic 3′ dA. ** marks the two covalent genomic-"
         "fragment/bridge ligations made in one proximity-ligation reaction."),
        ("Tagment, select the junction and amplify", [
            Row(chunks=[("S5—ME ** genomic tag — biotin bridge — genomic tag ** ME—S7",
                         "me", False)]),
            Row(chunks=[("                         ↓ streptavidin capture", "umi", False)]),
        ], "Nextera Tn5 cuts the ligation product and transfers adapters. Streptavidin "
           "retains only tagmented pieces containing the biotin bridge before indexed PCR."),
    ]
