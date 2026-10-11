"""Open-chromatin-anchored chromosome-conformation library chemistry.

HiCAR deliberately joins two unlike sequencing entries: a TruSeq Read 1 side reports the
restriction-fragment partner and a Nextera Read 2 side reports the Tn5-tagged accessible
anchor.  This module owns that hybrid construct and validates all four run primers.
"""
from __future__ import annotations

from types import SimpleNamespace

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, inline_sequence, revcomp
from protocol_model import chain, closed_circle, molecule, role


HICAR7_I0_WRITTEN = (
    "/5Phos/taGGTCTGAGTTCGCACAAAGGGAAGAGCGTCGTGTGGAGATGTGTATAAGAGACAG"
)
HICAR7_I0 = "TAGGTCTGAGTTCGCACAAAGGGAAGAGCGTCGTGTGG" + nx.ME
TRUSEQ_R1_SPLINT = "ACACGACGCTCTTCCGATCT"
I5 = il.NEBNEXT_I5_SET1["i501"]
I7 = "CAGCCTCG"
I501_PCR = il.nebnext_i5_primer(I5)
NEXTERA_I7_10_L = il.P7 + I7 + nx.S7 + nx.ME

SAMPLE_I5 = feature("hicar.i5", "sample_index", "fixed",
                    note="NEB i501 index carried by the defining library primer")
SAMPLE_I7 = feature("hicar.i7", "sample_index", "fixed",
                    note="Nextera i7-10-L index carried by the defining library primer")


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, top, tag, **kw)


def hicar_library() -> tuple[Construct, tuple]:
    """Finished HiCAR amplicon with its asymmetric primer geometry enforced."""
    lib = Construct([
        _seg("P5", il.P5, "p5"),
        _seg("i5 (i501)", I5, "cbc", feature=SAMPLE_I5),
        _seg("TruSeq Read 1 arm", il.TRUSEQ_READ1, "r1"),
        role("CviQI-cut interacting genomic DNA"),
        role("HiCAR7–splint proximity-ligation junction"),
        role("Tn5-tagged accessible genomic DNA"),
        _seg("mosaic end reverse complement", nx.ME_RC, "me"),
        _seg("S7 reverse complement", nx.S7_RC, "s7"),
        _seg("i7 reverse complement", revcomp(I7), "cbc", feature=SAMPLE_I7),
        _seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="HiCAR open-to-all contact library")
    primers = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"],
               sp.TRUSEQ["I2"], sp.NEXTERA["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("HiCAR hybrid library: " + "; ".join(problems))
    return lib, primers


def _contact_scene(label: str) -> Scene:
    return molecule(label, (
        role("CviQI-cut interacting genomic DNA"),
        role("HiCAR7 adapter and splint"),
        role("Tn5-tagged accessible genomic DNA"),
    ), junctions=(("CviQI-cut interacting genomic DNA", "HiCAR7 adapter and splint",
                  "proximity ligation"),
                 ("HiCAR7 adapter and splint", "Tn5-tagged accessible genomic DNA",
                  "Tn5 strand transfer")))


def hicar_workflow():
    initial = molecule("formaldehyde-crosslinked nuclei", (
        role("spatially proximal chromatin"), role("accessible regulatory DNA")))
    circle = closed_circle("NlaIII-trimmed HiCAR circle", (
        role("CviQI-cut interacting genomic DNA"),
        role("HiCAR7–splint junction"),
        role("Tn5-tagged accessible genomic DNA"),
    ), closure_label="dilute intramolecular T4 ligation")
    return chain(initial, (
        ("Tagment accessible chromatin with Tn5 loaded with engineered HiCAR7",
         molecule("tagged accessible chromatin", (
             role("chromatin outside the accessible site"),
             role("HiCAR7 transferred strand"),
             role("accessible genomic DNA")),
             junctions=(("HiCAR7 transferred strand", "accessible genomic DNA",
                         "Tn5 strand transfer"),)),
         "HiCAR7 carries the Tn5 mosaic end and a single-stranded 5-prime TA-compatible flank; the defining protocol starts with 30,000–100,000 fixed cells."),
        ("Digest nuclei with CviQI and splint-ligate spatial partners in situ",
         _contact_scene("open-to-all proximity-ligation product"),
         "The TruseqR1 splint bridges the engineered adapter to a CviQI-cut partner. This is the contact-recording covalent junction; Read 1 and Read 2 later report opposite sides."),
        ("Reverse crosslinks and fill the staggered Tn5 gaps with T4 DNA polymerase",
         _contact_scene("repaired linear HiCAR contact"),
         "Gap repair completes the adapter-bearing duplex before the secondary restriction digest."),
        ("Digest purified DNA with NlaIII and circularize at low concentration",
         circle,
         "NlaIII shortens the contact product; dilute T4 ligation closes individual molecules in the 4C-like library-preparation step."),
        ("Inverse-PCR between the splint side and mosaic-end side for 12 cycles",
         molecule("linear dual-index HiCAR amplicon", (
             role("P5, i501 and TruSeq Read 1 arm"),
             role("CviQI-cut interacting genomic DNA"),
             role("HiCAR7–splint contact junction"),
             role("Tn5-tagged accessible genomic DNA"),
             role("mosaic end, S7, i7-10-L and P7"))),
         "Outward primers open the selected circle into the sequencing library. The TruSeq side produces R1 from the interacting partner and the Nextera side produces R2 from accessible chromatin."),
    ))


def oligos():
    return (
        f"HiCAR7-i0 (as written, including 5-prime phosphate and lower-case TA flank) — {HICAR7_I0_WRITTEN}",
        f"TruseqR1 splint — {inline_sequence(TRUSEQ_R1_SPLINT, 'TruseqR1 splint')}",
        f"NEB primer i501 — {inline_sequence(I501_PCR, 'NEB i501 PCR primer')}",
        f"Nextera-pcr-i7-10-L — {inline_sequence(NEXTERA_I7_10_L, 'Nextera i7-10-L PCR primer')}",
    )


def get():
    final, primers = hicar_library()
    return SimpleNamespace(
        TITLE="HiCAR — open-chromatin-anchored chromosome conformation",
        NOTES="01_hicar.html",
        SOURCE='<a href="https://doi.org/10.1016/j.molcel.2022.01.023">Wei et al., <i>Molecular Cell</i> (2022)</a>; defining CviQI/NlaIII implementation.',
        SUMMARY="An engineered Tn5 adapter marks accessible chromatin and is splint-ligated in fixed nuclei to a spatially proximal CviQI fragment. NlaIII digestion, dilute recircularization and inverse PCR then recover an asymmetric contact library.",
        CAVEAT="The later authors’ Current Protocols version substitutes MseI at the first digest. This page preserves the defining paper’s CviQI implementation; the public HiCAR7-i0 oligo is shown exactly as written.",
        oligos=oligos,
        workflow=hicar_workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="The defining hybrid library has a TruSeq entry on the interacting-partner side and a Nextera entry on the accessible-anchor side. The HiCAR7–splint junction remains between those genomic regions; bracketed genomic regions vary in length.",
        SEQ_PRIMERS=primers,
        READ_LENGTHS={"Read 1": 150, "Index 1 (i7)": 8,
                      "Index 2 (i5)": 8, "Read 2": 150},
        SEQUENCING_INTRO="Read 1 enters the CviQI-cut interacting genomic partner from the TruSeq side. Read 2 enters the Tn5-tagged accessible anchor from the Nextera mosaic end. The two index reads identify the library; there is no molecule-level UMI.",
        SEQUENCING_ENDING="",
        SEQUENCING_UNAVAILABLE="",
    )
