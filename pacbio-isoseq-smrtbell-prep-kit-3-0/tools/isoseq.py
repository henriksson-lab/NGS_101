"""PacBio Iso-Seq v2 / SMRTbell prep kit 3.0 workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from chemdraw import Construct, Row, feature, feature_rows, revcomp, strand_row
from dumbbell import (smrtbell_circuit, smrtbell_hairpin_geometry,
                      smrtbell_sequencing_plan)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from rna_special import isoseq_smrtbell, template_switch_scene

TITLE = "PacBio Iso-Seq v2 — SMRTbell prep kit 3.0"
NOTES = "01_isoseq.html"
SOURCE = 'Commercial protocol: PacBio <a href="https://www.pacb.com/wp-content/uploads/Procedure-checklist-%E2%80%93-Preparing-Kinnex-compatible-cDNA-using-Iso-Seq-express-2.0-kit-1.pdf">Preparing Kinnex-compatible cDNA using Iso-Seq express 2.0 kit, 102-396-000 Rev10</a>; <a href="https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/IsoSeq_v2_primers_12.fasta">official primer FASTA</a>.'
SUMMARY = "Template-switch reverse transcription captures full-length poly(A) transcripts, cDNA is amplified, repaired and A-tailed, then hairpin adapters close both ends into a SMRTbell for circular-consensus sequencing."
CAVEAT = "The RT primer mix, TSO 2.0, SMRTbell adapter, and sequencing-primer bases remain proprietary. PacBio does publish the complete bc01–bc12 PCR primers and opposite cDNA amplification primer; those exact sequences and binding sites are shown."
ISOSEQ_5P_PREFIX = "CTACACGACGCTCTTCCGATCT"
ISOSEQ_5P_SUFFIX = "GCAATGAAGTCGCAGGGTTGGG"
ISOSEQ_BARCODES = (
    "ACTACAC", "ACTAGTA", "AGTGTAC", "ATCACTA", "CAGCTGT", "CAGTCAC",
    "CATGTAT", "CGTATGT", "GACATGT", "GAGTCTA", "GTAGATA", "GTATGAC",
)
ISOSEQ_3P = "AAGCAGTGGTATCAACGCAGAGTAC"


def barcoded_primer(number: int) -> str:
    if not 1 <= number <= len(ISOSEQ_BARCODES):
        raise ValueError("Iso-Seq barcode number must be 1..12")
    return ISOSEQ_5P_PREFIX + ISOSEQ_BARCODES[number - 1] + ISOSEQ_5P_SUFFIX


INSERT = Construct([
    seg("bc01 5-prime common arm", ISOSEQ_5P_PREFIX, "r1"),
    seg("bc01 sample barcode", ISOSEQ_BARCODES[0], "cbc",
        feature=feature("isoseq_sample_index", "sample_index", "whitelist",
                        whitelist="PacBio Iso-Seq primer barcodes bc01–12")),
    seg("bc01 3-prime common arm", ISOSEQ_5P_SUFFIX, "r1"),
    seg("full-length transcript cDNA", "X"*44, placeholder=True),
    seg("poly(A/T) boundary", "X"*12, placeholder=True),
    seg("cDNA amplification-primer reverse complement", revcomp(ISOSEQ_3P), "r2"),
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


def oligos():
    return [
        *(f"<b>Iso-Seq primer bc{i:02d}</b> — 5′-{barcoded_primer(i)}-3′"
          for i in range(1, len(ISOSEQ_BARCODES) + 1)),
        f"<b>Iso-Seq cDNA amplification primer (IsoSeqX_3p)</b> — 5′-{ISOSEQ_3P}-3′",
    ]


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
    rt_primer = g.input("rt_primer", "Iso-Seq RT primer mix", MaterialKind.OLIGO,
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
    fwd = g.input("forward", "published Iso-Seq bc01–bc12 primer family",
                  MaterialKind.OLIGO,
                  source=SOURCE + " (official IsoSeq_v2_primers_12.fasta)")
    rev = g.input("reverse", "published Iso-Seq cDNA amplification primer",
                  MaterialKind.OLIGO,
                  source=SOURCE + " (official IsoSeq_v2_primers_12.fasta)")
    pol = g.input("polymerase", "cDNA PCR polymerase", MaterialKind.CATALYST)
    amplified = g.product("amplified", "barcoded full-length duplex cDNA",
                          MaterialKind.MOLECULAR,
                          state=tuple([strand_row(INSERT), *feature_rows(INSERT)]),
                          source=SOURCE + " (exact terminal PCR-primer bases; transcript schematic)")
    g.react("pcr", "Amplify and barcode full-length cDNA", StepKind.PCR,
            ((full, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
             (rev, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (amplified,),
            geometry=PCRGeometry.exact(
                INSERT,
                PCRPrimer("Iso-Seq bc01", barcoded_primer(1)),
                PCRPrimer("Iso-Seq cDNA amplification primer", ISOSEQ_3P),
                statement=("PacBio's official FASTA publishes both complete primer "
                           "sequences; bc01 is the representative member of the "
                           "twelve-primer barcode family.")))
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
