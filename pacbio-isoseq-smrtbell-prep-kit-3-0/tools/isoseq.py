"""PacBio Iso-Seq v2 / SMRTbell prep kit 3.0 workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from chemdraw import Construct, Row, feature, feature_rows, revcomp, strand_row
from dumbbell import (Dumbbell, dumbbell_rows, smrtbell_circuit,
                      smrtbell_sequencing_plan)
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from rna_special import template_switch_scene

TITLE = "PacBio Iso-Seq v2 — SMRTbell prep kit 3.0"
NOTES = "01_isoseq.html"
SOURCE = 'Commercial protocol: PacBio <a href="https://www.pacb.com/wp-content/uploads/Procedure-checklist-%E2%80%93-Preparing-Kinnex-compatible-cDNA-using-Iso-Seq-express-2.0-kit-1.pdf">Preparing Kinnex-compatible cDNA using Iso-Seq express 2.0 kit, 102-396-000 Rev10</a>; <a href="https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/IsoSeq_v2_primers_12.fasta">official primer FASTA</a>; the monomer SMRTbell procedure that Rev10 links, <a href="https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-multiplexed-amplicon-libraries-using-SMRTbell-prep-kit-3.0.pdf">Preparing multiplexed amplicon libraries using SMRTbell prep kit 3.0, 102-359-000 Rev08</a>.'
SPK3 = ('PacBio 102-359-000 Rev08 sections 2–4 (repair & A-tailing, adapter ligation '
        'and 1.3X cleanup, nuclease treatment and 1.3X cleanup), the monomer route '
        'named by 102-396-000 Rev10 pp. 1 and 13')
SUMMARY = "Template-switch reverse transcription captures full-length poly(A) transcripts, cDNA is amplified with a barcoded primer, repaired and A-tailed, then hairpin adapters close both ends into a SMRTbell; a nuclease removes molecules that remain open before circular-consensus sequencing."
CAVEAT = "PacBio 102-396-000 Rev10 makes the Iso-Seq express 2.0 cDNA (12 PCR cycles for the SMRTbell prep kit 3.0 monomer route) and, since Rev09, sends monomer libraries to the multiplexed amplicon procedure 102-359-000 Rev08 for repair/A-tailing, adapter ligation and nuclease treatment. Earlier 102-396-000 revisions (Rev07, June 2025) contained their own SMRTbell 3.0 section; that revision currently returns 404 and was not read. Neither inspected procedure prints the SMRTbell adapter bases or its ligating end, so the hairpin overhang and full adapter sequence remain unverified, as do the RT primer mix, TSO 2.0 and standard sequencing-primer bases. PacBio does publish the complete bc01–bc12 PCR primers and opposite cDNA amplification primer; those exact sequences and binding sites are shown."
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
SMRTBELL = Dumbbell(INSERT, "N" * 18, "N" * 18, name="Iso-Seq SMRTbell",
                    insert_overhang="", adapter_overhang="")
SMRTBELL_ROWS = dumbbell_rows(SMRTBELL)
FINAL_LIBRARY = smrtbell_circuit(SMRTBELL)
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "The Standard sequencing primer anneals within the hairpin adapter whose bases remain unverified. A bound polymerase repeatedly traverses forward cDNA, the opposite hairpin, reverse-complement cDNA and the starting hairpin; those passes form one circular-consensus (HiFi) read."
SEQUENCING_UNAVAILABLE = ("PacBio 102-359-000 Rev08 section 5 anneals the Standard sequencing "
                          "primer and names the SMRTbell adapter, but neither inspected procedure "
                          "prints their bases; the earlier Iso-Seq-specific 102-396-000 Rev07 "
                          "SMRTbell section remains unread (404).")
FINAL_CAPTION = "The full-length cDNA is closed by two hairpins for repeated polymerase passes."
HAIRPIN_END = AdapterLigationGeometry.unavailable(
    "PacBio 102-359-000 Rev08 section 3 ligates the kit's SMRTbell adapter (ligation "
    "mix + enhancer, 20 °C 30 min) after repair & A-tailing, but prints neither its "
    "bases nor its ligating end; the earlier Iso-Seq-specific revision is unavailable.")


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
         "Twelve PCR cycles with the cDNA amplification primer and one of bc01–12 install the sample index; a 1.1X (or 0.9X) bead cleanup follows. Repair & A-tailing (37 °C 30 min, 65 °C 5 min) prepares both duplex ends for adapter ligation."),
        ("Ligate hairpins to make the SMRTbell", SMRTBELL_ROWS, "INFERRED — each ** boundary is a hairpin-adapter ligation; the exact ligating overhang is unresolved. The closed topology is checked, but its junction bases are not established. A nuclease treatment (37 °C 15 min) then digests molecules that remain open, and 1.3X bead cleanups follow ligation and nuclease treatment."),
        ("Polymerase repeatedly traverses the insert", [Row(chunks=[("standard sequencing primer → insert forward → hairpin → insert reverse complement → hairpin ↻", "r1", False)])], "Multiple passes over the same full-length cDNA molecule yield the circular-consensus sequence."),
    ]


def _bead_cleanup(g, identifier, label, template, beads, kept, reason):
    """Declare one SMRTbell bead cleanup whose supernatant is discarded."""
    waste = g.product(identifier + "_supernatant",
                      f"supernatant not retained by {beads.label}", MaterialKind.WASTE)
    g.react(identifier, label, StepKind.CLEANUP,
            ((template, InputRole.TEMPLATE), (beads, InputRole.CONDITION)),
            (kept, waste), changes_molecule=False, reason=reason)
    g.discard(waste, f"supernatant discarded after {beads.label} binding")
    return kept


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    rna = g.input("rna", "300 ng total RNA (RIN ≥7) containing poly(A) RNA",
                  MaterialKind.MOLECULAR, primary=True)
    rt_primer = g.input("rt_primer", "Iso-Seq RT primer mix", MaterialKind.OLIGO,
                        inferred=True, source=SOURCE + " (current primer bases not verified in inspected sources)")
    rt = g.input("rt", "Iso-Seq RT enzyme mix (template-switching reverse transcriptase)",
                 MaterialKind.CATALYST)
    first = g.product("first", "first-strand full-length cDNA", MaterialKind.MOLECULAR)
    g.react("reverse_transcribe",
            "Anneal RT primer mix (70 °C 5 min), reverse-transcribe 42 °C 45 min",
            StepKind.REVERSE_TRANSCRIPTION,
            ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (first,),
            source=SOURCE + " (102-396-000 Rev10 sections 2.2–2.3)")
    tso = g.input("tso", "Iso-Seq Express TSO 2.0", MaterialKind.OLIGO,
                  inferred=True, source=SOURCE + " (current oligo bases not verified in inspected sources)")
    full = g.product("full", "full-length cDNA with terminal PCR handles",
                     MaterialKind.MOLECULAR)
    g.react("template_switch", "Add TSO 2.0 and template-switch 42 °C 15 min",
            StepKind.TEMPLATE_SWITCHING,
            ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (full,),
            source=SOURCE + " (102-396-000 Rev10 steps 2.3.5–2.3.8)")
    rt_beads = g.input("rt_beads", "1.3X SMRTbell cleanup beads", MaterialKind.CONDITION)
    full_clean = g.product("full_clean", "purified full-length cDNA", MaterialKind.MOLECULAR)
    _bead_cleanup(g, "rt_cleanup", "Purify RT/template-switch reaction with 1.3X beads",
                  full, rt_beads, full_clean,
                  "102-396-000 Rev10 section 2.4: 29 µL elution buffer + 65 µL beads, "
                  "two 80% ethanol washes, 21 µL elution.")
    fwd = g.input("forward", "published Iso-Seq bc01–bc12 primer family",
                  MaterialKind.OLIGO,
                  state=tuple(barcoded_primer(i) for i in range(1, len(ISOSEQ_BARCODES) + 1)),
                  source=SOURCE + " (official IsoSeq_v2_primers_12.fasta)")
    rev = g.input("reverse", "published Iso-Seq cDNA amplification primer",
                  MaterialKind.OLIGO, state=ISOSEQ_3P,
                  source=SOURCE + " (official IsoSeq_v2_primers_12.fasta)")
    pol = g.input("polymerase", "Iso-Seq cDNA PCR mix", MaterialKind.CATALYST)
    amplified = g.product("amplified", "barcoded full-length duplex cDNA",
                          MaterialKind.MOLECULAR,
                          state=tuple([strand_row(INSERT), *feature_rows(INSERT)]),
                          source=SOURCE + " (exact terminal PCR-primer bases; transcript schematic)")
    g.react("pcr", "Amplify and barcode full-length cDNA (12 cycles for SMRTbell prep kit 3.0)",
            StepKind.PCR,
            ((full_clean, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
             (rev, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (amplified,),
            source=SOURCE + " (102-396-000 Rev10 section 3.1: 98 °C 45 s; 12 × "
                            "[98 °C 10 s, 60 °C 15 s, 72 °C 3 min]; 72 °C 5 min)",
            geometry=PCRGeometry.exact(
                INSERT,
                PCRPrimer("Iso-Seq bc01", barcoded_primer(1)),
                PCRPrimer("Iso-Seq cDNA amplification primer", ISOSEQ_3P),
                statement=("PacBio's official FASTA publishes both complete primer "
                           "sequences; bc01 is the representative member of the "
                           "twelve-primer barcode family.")))
    pcr_beads = g.input("pcr_beads", "1.1X (or 0.9X) SMRTbell cleanup beads",
                        MaterialKind.CONDITION)
    amplified_clean = g.product("amplified_clean", "size-enriched barcoded cDNA (pooled when multiplexed)",
                                MaterialKind.MOLECULAR)
    _bead_cleanup(g, "pcr_cleanup", "Purify amplified cDNA with 1.1X (or 0.9X) beads",
                  amplified, pcr_beads, amplified_clean,
                  "102-396-000 Rev10 sections 3.2–3.3: 45 µL PCR + 49.5 µL (1.1X) or "
                  "40.5 µL (0.9X, longer cDNA) beads, 47 µL elution for SMRTbell prep "
                  "kit 3.0; ≥100 ng per sample; equal-mass pooling up to 500 ng is "
                  "not drawn as a separate step.")
    end_mix = g.input("end_mix", "Repair buffer, End repair mix and DNA repair mix",
                      MaterialKind.CATALYST, source=SOURCE + " (" + SPK3 + ")")
    tailed = g.product("tailed", "repaired dA-tailed full-length cDNA",
                       MaterialKind.MOLECULAR)
    g.react("end_prep", "Repair and A-tail amplified cDNA (37 °C 30 min, 65 °C 5 min)",
            StepKind.END_REPAIR,
            ((amplified_clean, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,),
            source=SOURCE + " (" + SPK3 + ")")
    hairpins = g.input("hairpins", "SMRTbell adapter; ligating end unresolved",
                       MaterialKind.ADAPTER, inferred=True,
                       source=SOURCE + " (adapter bases not printed in inspected sources)")
    ligase = g.input("ligase", "Ligation mix + Ligation enhancer", MaterialKind.CATALYST,
                     source=SOURCE + " (" + SPK3 + ")")
    closed = g.product("closed", "hairpin-ligated cDNA (closed SMRTbells and incompletely ligated molecules)",
                       MaterialKind.MOLECULAR, state=tuple(SMRTBELL_ROWS), inferred=True,
                       source=SOURCE + " (hairpin bases not verified in inspected sources)")
    g.react("ligate", "Ligate hairpins to both cDNA ends (20 °C 30 min)", StepKind.LIGATION,
            ((tailed, InputRole.TEMPLATE), (hairpins, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (closed,), geometry=HAIRPIN_END,
            source=SOURCE + " (" + SPK3 + ")")
    lig_beads = g.input("lig_beads", "1.3X SMRTbell cleanup beads (post-ligation)",
                        MaterialKind.CONDITION)
    closed_clean = g.product("closed_clean", "purified ligation products",
                             MaterialKind.MOLECULAR, state=tuple(SMRTBELL_ROWS), inferred=True,
                             source=SOURCE + " (hairpin bases not verified in inspected sources)")
    _bead_cleanup(g, "lig_cleanup", "Purify ligation with 1.3X beads", closed, lig_beads,
                  closed_clean, "102-359-000 Rev08 section 3: 124 µL beads to 95 µL, "
                  "40 µL elution.")
    nuclease = g.input("nuclease", "Nuclease buffer + Nuclease mix", MaterialKind.CATALYST,
                       source=SOURCE + " (" + SPK3 + ")")
    digested = g.product("digested", "nuclease-resistant closed SMRTbells and digested open molecules",
                         MaterialKind.MOLECULAR, state=tuple(SMRTBELL_ROWS), inferred=True,
                         source=SOURCE + " (hairpin bases not verified in inspected sources)")
    g.react("remove_open", "Nuclease treatment of molecules retaining a free end (37 °C 15 min)",
            StepKind.DIGESTION,
            ((closed_clean, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
            (digested,), changes_molecule=False,
            source=SOURCE + " (" + SPK3 + "; the substrate specificity of the nuclease "
                            "mix is not stated in the inspected procedures)")
    nuc_beads = g.input("nuc_beads", "1.3X SMRTbell cleanup beads (post-nuclease)",
                        MaterialKind.CONDITION)
    library = g.product("library", "closed full-length Iso-Seq SMRTbell library",
                        MaterialKind.LIBRARY, state=tuple(SMRTBELL_ROWS), inferred=True,
                        source=SOURCE + " (hairpin bases not verified in inspected sources)")
    _bead_cleanup(g, "nuc_cleanup", "Purify nuclease-treated library with 1.3X beads",
                  digested, nuc_beads, library,
                  "102-359-000 Rev08 section 4: 65 µL beads to 50 µL, 26 µL elution; "
                  "digestion fragments are removed with the supernatant.")
    g.sequence(library, smrtbell_sequencing_plan(
        SMRTBELL, FINAL_LIBRARY, limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("rna", "rt_primer", "rt", "reverse_transcribe", "tso",
                              "template_switch", "rt_beads", "rt_cleanup", "full_clean",
                              "forward", "reverse", "polymerase", "pcr", "pcr_beads",
                              "pcr_cleanup", "amplified_clean", "end_mix", "end_prep",
                              "hairpins", "ligase", "ligate", "lig_beads", "lig_cleanup",
                              "closed_clean", "nuclease", "remove_open", "nuc_beads",
                              "nuc_cleanup", "library"), SOURCE)
    return g.require_publishable()
