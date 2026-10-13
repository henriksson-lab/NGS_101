"""PacBio Iso-Seq v2 / SMRTbell prep kit 3.0 workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from chemdraw import Construct, Row, Segment, feature, feature_rows, strand_row
from dumbbell import (smrtbell_circuit, smrtbell_hairpin_geometry,
                      smrtbell_sequencing_plan)
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from rna_special import isoseq_smrtbell, template_switch_scene

TITLE = "PacBio Iso-Seq v2 — SMRTbell prep kit 3.0"
NOTES = "01_isoseq.html"
SOURCE = 'Commercial protocol: PacBio <a href="https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-Iso-Seq-libraries-using-SMRTbell-prep-kit-3.0.pdf">Preparing Iso-Seq v2 libraries using SMRTbell prep kit 3.0, 102-396-000 Rev07</a>.'
SUMMARY = "Template-switch reverse transcription captures full-length poly(A) transcripts, cDNA is amplified, repaired and A-tailed, then hairpin adapters close both ends into a SMRTbell for circular-consensus sequencing."
CAVEAT = "Current cDNA-synthesis oligos and SMRTbell adapter bases are proprietary. Bracketed/unknown regions preserve molecular roles and topology without inventing sequences."
INSERT = Construct([
    seg("Iso-Seq primer barcode", "[Iso-Seq barcode]", "cbc", placeholder=True,
        inferred=True,
        feature=feature("isoseq_sample_index", "sample_index", "whitelist",
                        whitelist="PacBio Iso-Seq primer barcodes bc01–12")),
    seg("5-prime cDNA handle", "X"*14, "r1", placeholder=True, inferred=True),
    seg("full-length transcript cDNA", "X"*44, placeholder=True),
    seg("poly(A/T) boundary", "X"*12, placeholder=True),
    seg("3-prime cDNA handle", "X"*14, "r2", placeholder=True, inferred=True),
], name="barcoded full-length Iso-Seq cDNA")
SMRTBELL, SMRTBELL_ROWS = isoseq_smrtbell(INSERT)
FINAL_LIBRARY = smrtbell_circuit(SMRTBELL)
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "The SMRTbell sequencing primer anneals within the proprietary hairpin adapter. A bound polymerase repeatedly traverses forward cDNA, the opposite hairpin, reverse-complement cDNA and the starting hairpin; those passes form one circular-consensus (HiFi) read."
SEQUENCING_UNAVAILABLE = ("The procedure names the PacBio sequencing primer and hairpin "
                          "adapter but does not disclose their current bases.")
FINAL_CAPTION = "The full-length cDNA is closed by two hairpins for repeated polymerase passes."
HAIRPIN_END = smrtbell_hairpin_geometry(
    SMRTBELL,
    statement="SMRTbell prep kit 3.0 ligates T-overhang hairpins to dA-tailed cDNA.")


def final_rows():
    return SMRTBELL_ROWS

def sections():
    return [
        ("Make full-length cDNA by template switching", template_switch_scene(inferred=True).rows(), "INFERRED — the vendor discloses oligo roles but not current bases. Poly(A)-primed RT and template switching place PCR handles at both cDNA ends."),
        ("Amplify, barcode and repair full-length cDNA",
         [strand_row(INSERT), *feature_rows(INSERT, prefix_width=4),
          Row(chunks=[("PCR → damage repair → end repair → dA tailing", None, False)])],
         "The required bc01–12 primer installs the sample index during cDNA PCR. Repair and A-tailing make both duplex ends adapter-compatible."),
        ("Ligate hairpins to make the SMRTbell", SMRTBELL_ROWS, "Each ** boundary is a hairpin-adapter ligation. The Dumbbell constructor rejects missing ends or incompatible A:T overhangs."),
        ("Polymerase repeatedly traverses the insert", [Row(chunks=[("primer → insert forward → hairpin → insert reverse complement → hairpin ↻", "r1", False)])], "Multiple passes over the same full-length cDNA molecule yield the circular-consensus sequence."),
    ]


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    rna = g.input("rna", "poly(A) RNA", MaterialKind.MOLECULAR, primary=True)
    rt_primer = g.input("rt_primer", "barcoded Iso-Seq RT primer", MaterialKind.OLIGO,
                        inferred=True, source=SOURCE + " (current primer bases proprietary)")
    rt = g.input("rt", "template-switching reverse transcriptase", MaterialKind.CATALYST)
    first = g.product("first", "first-strand full-length cDNA", MaterialKind.MOLECULAR)
    g.react("reverse_transcribe", "Reverse-transcribe poly(A) RNA",
            StepKind.REVERSE_TRANSCRIPTION,
            ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (first,))
    tso = g.input("tso", "Iso-Seq template-switch oligo", MaterialKind.OLIGO,
                  inferred=True, source=SOURCE + " (current oligo bases proprietary)")
    full = g.product("full", "full-length cDNA with terminal PCR handles",
                     MaterialKind.MOLECULAR)
    g.react("template_switch", "Capture the completed cDNA end by template switching",
            StepKind.TEMPLATE_SWITCHING,
            ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (full,))
    fwd = g.input("forward", "Iso-Seq PCR forward primer", MaterialKind.OLIGO,
                  inferred=True, source=SOURCE + " (primer bases proprietary)")
    rev = g.input("reverse", "barcoded Iso-Seq PCR reverse primer", MaterialKind.OLIGO,
                  inferred=True, source=SOURCE + " (primer bases proprietary)")
    pol = g.input("polymerase", "cDNA PCR polymerase", MaterialKind.CATALYST)
    amplified = g.product("amplified", "barcoded full-length duplex cDNA",
                          MaterialKind.MOLECULAR,
                          state=tuple([strand_row(INSERT), *feature_rows(INSERT)]), inferred=True,
                          source=SOURCE + " (terminal handle bases proprietary)")
    g.react("pcr", "Amplify and barcode full-length cDNA", StepKind.PCR,
            ((full, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
             (rev, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (amplified,),
            geometry=PCRGeometry.unavailable(
                INSERT, reason="The current Iso-Seq PCR-primer bases are proprietary."))
    end_mix = g.input("end_mix", "damage-repair, end-repair and dA-tailing enzymes",
                      MaterialKind.CATALYST)
    tailed = g.product("tailed", "repaired dA-tailed full-length cDNA",
                       MaterialKind.MOLECULAR)
    g.react("end_prep", "Repair and dA-tail amplified cDNA", StepKind.END_REPAIR,
            ((amplified, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,))
    hairpins = g.input("hairpins", "T-overhang SMRTbell hairpin adapters",
                       MaterialKind.ADAPTER, inferred=True,
                       source=SOURCE + " (adapter bases proprietary)")
    ligase = g.input("ligase", "SMRTbell ligase", MaterialKind.CATALYST)
    library = g.product("library", "closed full-length Iso-Seq SMRTbell library",
                        MaterialKind.LIBRARY, state=tuple(SMRTBELL_ROWS), inferred=True,
                        source=SOURCE + " (hairpin bases proprietary)")
    g.react("ligate", "Ligate hairpins to both cDNA ends", StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (hairpins, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (library,), geometry=HAIRPIN_END)
    g.sequence(library, smrtbell_sequencing_plan(
        SMRTBELL, FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("rna", "rt_primer", "rt", "reverse_transcribe", "tso",
                              "template_switch", "forward", "reverse", "polymerase",
                              "pcr", "end_mix", "end_prep", "hairpins", "ligase",
                              "ligate", "library"), SOURCE)
    return g.require_publishable()
