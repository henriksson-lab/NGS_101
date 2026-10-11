"""Takahashi et al. cap-trapper CAGE model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from captrap import CapTrappedHybrid
from chemdraw import Construct, Row, Scene, Segment, feature, revcomp

TITLE = "CAGE — cap analysis of gene expression"
NOTES = "01_cage.html"
SOURCE = ('Defining protocol: <a href="https://doi.org/10.1038/nprot.2012.005">'
          'Takahashi et al., <i>Nature Protocols</i> (2012)</a>.')
SUMMARY = ("Select cDNAs that reach a biotinylated RNA cap, ligate a barcoded 5′ linker, "
           "make the second strand, release a 27-nt 5′ tag with EcoP15I, and ligate the "
           "opposite Illumina linker.")

RT_ECOP = "AAGGTCTATCAGCAG" + "N" * 15
LINKER_CORE = "CCACCGACAGGTTCAGAGTTCTACAG"
BARCODE = "AGA"
ECOP_SITE = "CAGCAG"
FWD = "AATGATACGGCGACCACCGACAGGTTCAGAGTTC"
REV = "CAAGCAGAAGACGGCATACGA"
SEQUENCING = "CGGCGACCACCGACAGGTTCAGAGTTCTACAG"
THREE_UPPER = "NN" + revcomp(REV)
CAPTURED = CapTrappedHybrid(cdna_reaches_cap=True, cap_biotinylated=True)


def final_library() -> tuple[Construct, tuple]:
    # The 27-nt sequencing-primer prefix overlaps the forward PCR primer; the final
    # seven bases are supplied by the ligated 5′ linker.
    lib = Construct([
        Segment("P5-side PCR prefix", FWD[:7], "p5"),
        Segment("custom CAGE sequencing-primer site", SEQUENCING, "r1"),
        Segment("3-nt sample barcode", BARCODE, "cbc",
                feature=feature("sample_inline", "sample_index", "whitelist", whitelist="published CAGE linker set")),
        Segment("EcoP15I recognition site", ECOP_SITE),
        Segment("27-nt capped 5′ tag", "X" * 27, placeholder=True),
        Segment("two-base EcoP15I ligation overhang", "N" * 2, placeholder=True),
        Segment("P7-side linker / reverse-primer complement", revcomp(REV), "p7"),
    ], name="96-bp CAGE tag library")
    primer = sp.custom("Read 1", "CAGE custom sequencing primer", SEQUENCING,
                       "Takahashi et al. 2012, Table 1")
    problems = sp.verify(lib, (primer,), required_roles=("Read 1",))
    if problems:
        raise ValueError("invalid CAGE endpoint: " + "; ".join(problems))
    return lib, (primer,)


FINAL_LIBRARY, SEQ_PRIMERS = final_library()
FINAL_CAPTION = ("The 96-bp single-read product. The custom primer reads the three-base "
                 "sample barcode and EcoP15I site before the 27-nt transcript 5′ tag.")
SEQUENCING_INTRO = ("This protocol explicitly uses its Table 1 custom CAGE sequencing "
                    "primer, not the standard Illumina Read 1 primer.")


def sections():
    hybrid = [
        Row(chunks=[("5′ cap—RNA ================================ 3′", "w1", False)]),
        Row(chunks=[("          cDNA <========================= random N15–EcoP15I", "r1", False)]),
        Row(chunks=[("oxidize cap diol → biotin hydrazide → RNase I removes unprotected RNA → streptavidin capture", None, False)]),
    ]
    lib = FINAL_LIBRARY
    sc = Scene.duplex(list(lib), label="CAGE tag")
    sc.junction("top", "EcoP15I recognition site", "27-nt capped 5′ tag", "5′ linker ligation")
    sc.junction("top", "two-base EcoP15I ligation overhang", "P7-side linker / reverse-primer complement", "3′ linker ligation")
    sc.labels("top")
    five = AdapterLigationGeometry.unavailable("The CAGE protocol prints the 5′ linker composition and ligated product but does not provide a strand-resolved cut-facing end from which an overhang can be verified.")
    three_sc=Scene();three_sc.strand("3′ linker",[Segment("linker duplex arm",revcomp(REV),"p7"),Segment("two-base protrusion","NN",placeholder=True)],label="3′ CAGE linker")
    three_sc.overhang("3′ linker","two-base protrusion",polarity="3-prime",owner="adapter",text="published two-base protrusion; bases variable")
    three_sc.strand("EcoP15I tag",[Segment("two-base EcoP15I end","NN",placeholder=True),Segment("27-nt cap tag","X"*27,placeholder=True)],label="EcoP15I-cut tag")
    three_sc.overhang("EcoP15I tag","two-base EcoP15I end",polarity="5-prime",owner="substrate",text="two-base EcoP15I end")
    three=AdapterLigationGeometry.known("cohesive",three_sc,statement="The source specifies a two-base protrusion on the 3′ linker; its variable bases pair with the two-base EcoP15I tag end.")
    return [
        ("Cap trap full-length first strands", hybrid,
         "Only RNA/cDNA hybrids protected through the oxidized, biotinylated cap survive capture."),
        ("Install the 5′ linker and make the second strand", five.rows()+[Row(chunks=[
            ("captured cDNA → barcoded partial-duplex linker ligation → second-strand synthesis", None, False)])],
         "The linker supplies the custom sequencing site, a three-base sample barcode and an EcoP15I site."),
        ("Cut a 27-nt tag and ligate the 3′ linker", three.rows()+sc.rows(),
         "EcoP15I cuts 27 nt into the cap-derived cDNA; a linker with a two-base protrusion closes the other end. ** marks ligations."),
    ]
