"""Telomere Profiling enrichment workflow (Karimian et al. 2024)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import (Construct, Row, Scene, Segment, complement_segments, feature,
                      feature_rows, oligo, revcomp, strand_row)
from adapter_ligation import AdapterLigationGeometry
from telomere import AffinityCapture, HUMAN_C_REPEAT, TelomereEnd

TELOTAG_TOPS = {
    "NB50": "CCCTCCGATAATGGACTTTGGTAACTTCCTGCGTTGCGTACAGCAATCAGGCGCGCCGAATTCCGATCGTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB65": "CCCTCCGATATTCTCAGTCTTCCTCCAGACAAGGTGCGTACAGCAATCAGGCGCGCCGAATTCCGATCGTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB68": "CCCTCCGATAGAATCTAAGCAAACACGAAGGTGGTGCGTACAGCAATCAGGCGCGCCGAATTCCGATCGTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB88": "CCCTCCGATATCTTCTACTACCGATCCGAAGCAGTGCGTACAGCAATCAGGCGCGCCGAATTCCGATCGTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB01": "CCCTCCGATACACAAAGACACCGACAACTTTCTTTGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
    "NB70": "CCCTCCGATAACCGAGATCCTACGAATGGAGTGTTGCGTACAGCAATCAGGCGCGCCCGATCGGAATTCTTAATTAATGGTGGACTCTATGACCGTTCA",
}
TELOMERE_ARMS = tuple(
    (HUMAN_C_REPEAT * 4)[offset:offset + 18]
    for offset in (5, 4, 3, 2, 1, 0)
)


def top_parts(name="NB50"):
    sequence = TELOTAG_TOPS[name]
    return (
        Segment("5-prime TeloTag ligating flank", sequence[:10], "r3"),
        Segment("24-nt TeloTag sample barcode", sequence[10:34], "cbc",
                feature=feature("telotag_sample_index", "sample_index", "whitelist",
                                whitelist="Karimian et al. Extended Data Table 2")),
        Segment("barcode flank", sequence[34:49]),
        Segment("AscI", sequence[49:57], "r3"),
        Segment("EcoRI", sequence[57:63], "r3"),
        Segment("PvuI-HF", sequence[63:69], "r3"),
        Segment("PacI", sequence[69:77], "r3"),
        Segment("paired post-PacI flank", sequence[77:87]),
        Segment("3-prime biotin tail", sequence[87:]),
    )


def splint_family(name="NB50"):
    """The six exact phase splints, derived from their paired TeloTag strand."""
    paired_top = TELOTAG_TOPS[name][:-12]
    return tuple(revcomp(paired_top) + arm for arm in TELOMERE_ARMS)


TELOTAG_CONSTRUCT=Construct(top_parts(), name="NB50 barcoded TeloTag")
CAPTURE=AffinityCapture("biotin","adapter-encoded restriction site")


def oligos():
    top_rows = [
        oligo(name + " TeloTag", top_parts(name), mods="/5Phos/",
              three="/3BioTEG/-3′")
        for name in TELOTAG_TOPS
    ]
    splint_rows = [
        oligo(f"NB50.C{i}", (Segment(f"NB50.C{i}", sequence),))
        for i, sequence in enumerate(splint_family(), 1)
    ]
    return (*top_rows, *splint_rows)


def _telotag_scene():
    paired = complement_segments(top_parts()[:-1])
    sc = Scene()
    sc.strand("splint", [*paired, Segment("telomere-complement arm", TELOMERE_ARMS[0], "r2")],
              label="NB50.C1 splint")
    sc.anneal("TeloTag", top_parts(), to="splint",
              pair=("5-prime TeloTag ligating flank", "5-prime TeloTag ligating flank'"),
              unpaired=("3-prime biotin tail",), label="NB50 TeloTag",
              above=True, mod5="p", mod3="BioTEG")
    end = TelomereEnd(overhang_copies=3)
    chromosome = [Segment("subtelomere", "X" * end.subtelomere_bases, placeholder=True),
                  Segment("duplex telomere", end.duplex, "r1"),
                  # The six splints cover the six cyclic telomere phases.  Show
                  # the chromosome phase that actually pairs to NB50.C1 rather
                  # than forcing the canonical repeat origin onto this example.
                  Segment("G-rich overhang", revcomp(TELOMERE_ARMS[0]), "r2")]
    sc.anneal("G-rich strand", chromosome, to="splint",
              pair=("G-rich overhang", "telomere-complement arm"),
              unpaired=("subtelomere", "duplex telomere"), label="chromosome")
    bottom = complement_segments(chromosome[:2])
    sc.anneal("C-rich strand", bottom, to="G-rich strand",
              pair=("duplex telomere'", "duplex telomere"), label="chromosome")
    sc.mark("TeloTag", "24-nt TeloTag sample barcode", "published 24-nt barcode")
    sc.ligation_end("TeloTag", "5-prime TeloTag ligating flank", kind="splinted",
                    polarity="5-prime", owner="adapter",
                    text="splinted terminal TeloTag ligation")
    sc.footer("** HiFi Taq DNA ligase seals the TeloTag 5′ phosphate at the chromosome end",
              "TeloTag", "5-prime TeloTag ligating flank")
    return sc

def _ont_geometry():
    return AdapterLigationGeometry.unavailable("The Telomere Profiling protocol invokes the historical ONT ligation workflow, but the cited source does not disclose its complete motor-bearing adapter assembly.")

def panels(): return [("Digest high-molecular-weight genomic DNA away from telomere ends",[Row(chunks=[("internal restriction cut | subtelomere -- telomere -- native 3′ overhang",None,False)])],"Restriction reduces non-telomeric length while leaving the native chromosome terminus available for tagging."),("Anneal six telomeric splints and ligate the biotinylated TeloTag",_telotag_scene(),"NB50 and its six phase splints are source-exact. The splint bridges the TeloTag to the native G-rich overhang; ** marks the terminal ligation."),("Select tagged ends on streptavidin and release at the adapter restriction site",[Row(chunks=[("biotin-TeloTag -> streptavidin capture -> wash -> restriction release",None,False)])],"Affinity selection removes most untagged genomic fragments; release retains the tagged telomere molecule."),("ONT adapter terminal geometry",_ont_geometry().rows(),_ont_geometry().statement),("Attach Oxford Nanopore sequencing adapters",[Row(chunks=[("[motor-loaded ONT adapter] ** [enriched TeloTag–telomere fragment]",None,False)])],"The selected material enters the historical ligation library workflow; ** marks adapter ligation.")]
