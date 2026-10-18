"""PacBio Iso-Seq v2 / SMRTbell prep kit 3.0 workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from chemdraw import Construct, Scene, Segment, complement_segments, feature, feature_rows, revcomp, strand_row
from dumbbell import Dumbbell, dumbbell_rows
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import role, molecule

TITLE = "PacBio Iso-Seq v2 — SMRTbell prep kit 3.0"
NOTES = "01_isoseq.html"
SOURCE = 'Commercial protocol: PacBio <a href="https://www.pacb.com/wp-content/uploads/Procedure-checklist-%E2%80%93-Preparing-Kinnex-compatible-cDNA-using-Iso-Seq-express-2.0-kit-1.pdf">Preparing Kinnex-compatible cDNA using Iso-Seq express 2.0 kit, 102-396-000 Rev10</a>; <a href="https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/IsoSeq_v2_primers_12.fasta">official primer FASTA</a>; the monomer SMRTbell procedure that Rev10 links, <a href="https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-multiplexed-amplicon-libraries-using-SMRTbell-prep-kit-3.0.pdf">Preparing multiplexed amplicon libraries using SMRTbell prep kit 3.0, 102-359-000 Rev08</a>.'
SPK3 = ('PacBio 102-359-000 Rev08 sections 2–4 (repair & A-tailing, adapter ligation '
        'and 1.3X cleanup, nuclease treatment and 1.3X cleanup), the monomer route '
        'named by 102-396-000 Rev10 pp. 1 and 13')
SUMMARY = "Template-switch reverse transcription captures full-length poly(A) transcripts, cDNA is amplified with a barcoded primer, repaired and A-tailed, then hairpin adapters close both ends into a SMRTbell; a nuclease removes molecules that remain open before circular-consensus sequencing."
CAVEAT = "PacBio 102-396-000 Rev10 makes the Iso-Seq express 2.0 cDNA (12 PCR cycles for the SMRTbell prep kit 3.0 monomer route) and, since Rev09, sends monomer libraries to the multiplexed amplicon procedure 102-359-000 Rev08 for repair/A-tailing, adapter ligation and nuclease treatment. Earlier 102-396-000 revisions (Rev07, June 2025) contained their own SMRTbell 3.0 section; that revision currently returns 404 and was not read. Neither inspected procedure prints the SMRTbell adapter bases or its ligating end, so the hairpin overhang and full adapter sequence remain unverified, as do the RT primer mix, TSO 2.0 and standard sequencing-primer bases. PacBio does publish the complete bc01–bc12 PCR primers and opposite cDNA amplification primer; those exact sequences are shown. bc01 is the explicit representative used here, while the other eleven remain alternatives. Their initial sites on pre-PCR cDNA have not yet been reconstructed; PCR geometry uses the actual unbarcoded precursor without inventing a previously installed barcode."
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
                    insert_overhang="", adapter_overhang="", junction_unavailable=True)
SMRTBELL_ROWS = dumbbell_rows(SMRTBELL)
# Keep source-owned barcode features and unknown molecular lengths through the
# traversal; the 18-character hairpin bodies are drawing widths, never 18 nt claims.
FINAL_LIBRARY = Construct([
    *list(INSERT),
    Segment("right hairpin; bases and length unresolved", SMRTBELL.right_adapter, "r1", placeholder=True, inferred=True),
    *complement_segments(list(INSERT)),
    Segment("left hairpin; bases and length unresolved", SMRTBELL.left_adapter, "r1", placeholder=True, inferred=True),
], name="Iso-Seq monomer polymerase circuit; schematic unknown regions")
if FINAL_LIBRARY.top().upper() != SMRTBELL.pass_sequence.upper():
    raise ValueError("local Iso-Seq traversal differs from its declared dumbbell")
PRE_PCR = molecule("unbarcoded post-template-switch cDNA", (
    role("TSO-2.0-derived cDNA end; not yet reconstructed"),
    role("full-length transcript cDNA"),
    role("RT-primer-derived cDNA end; not yet reconstructed"),
))
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
                  MaterialKind.MOLECULAR, state=molecule("poly(A) RNA substrate", (role("RNA body"), role("native poly(A) tail; variable length"))), primary=True)
    rt_primer = g.input("rt_primer", "Iso-Seq RT primer mix", MaterialKind.OLIGO, state=molecule("kit RT primer mix", (role("source-named RT primer species; sequence/length not yet assigned"),)),
                        inferred=True, source=SOURCE + " (current primer bases not verified in inspected sources)")
    rt = g.input("rt", "Iso-Seq RT enzyme mix (template-switching reverse transcriptase)",
                 MaterialKind.CATALYST)
    first = g.product("first", "first-strand full-length cDNA", MaterialKind.MOLECULAR, state=molecule("first-strand cDNA", (role("RNA-derived cDNA"), role("RT-primer-derived end; not yet reconstructed"))))
    g.react("reverse_transcribe",
            "Anneal RT primer mix (70 °C 5 min), reverse-transcribe 42 °C 45 min",
            StepKind.REVERSE_TRANSCRIPTION,
            ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (first,),
            source=SOURCE + " (102-396-000 Rev10 sections 2.2–2.3)")
    tso = g.input("tso", "Iso-Seq Express TSO 2.0", MaterialKind.OLIGO,
                  inferred=True, source=SOURCE + " (current oligo bases not verified in inspected sources)")
    full = g.product("full", "full-length cDNA with terminal PCR handles",
                     MaterialKind.MOLECULAR, state=PRE_PCR)
    g.react("template_switch", "Add TSO 2.0 and template-switch 42 °C 15 min",
            StepKind.TEMPLATE_SWITCHING,
            ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (full,),
            source=SOURCE + " (102-396-000 Rev10 steps 2.3.5–2.3.8)")
    rt_beads = g.input("rt_beads", "1.3X SMRTbell cleanup beads", MaterialKind.CONDITION)
    full_clean = g.product("full_clean", "purified full-length cDNA", MaterialKind.MOLECULAR, state=PRE_PCR)
    _bead_cleanup(g, "rt_cleanup", "Purify RT/template-switch reaction with 1.3X beads",
                  full, rt_beads, full_clean,
                  "102-396-000 Rev10 section 2.4: 29 µL elution buffer + 65 µL beads, "
                  "two 80% ethanol washes, 21 µL elution.")
    fwd = g.input("forward", "published Iso-Seq bc01 primer (representative choice, one barcode primer per sample)",
                  MaterialKind.OLIGO,
                  state=Construct(list(INSERT)[:3], name="source-selected representative bc01 primer"),
                  source=SOURCE + " (official IsoSeq_v2_primers_12.fasta)")
    rev = g.input("reverse", "published Iso-Seq cDNA amplification primer",
                  MaterialKind.OLIGO, state=ISOSEQ_3P,
                  source=SOURCE + " (official IsoSeq_v2_primers_12.fasta)")
    pol = g.input("polymerase", "Iso-Seq cDNA PCR mix", MaterialKind.CATALYST)
    amplified = g.product("amplified", "barcoded full-length duplex cDNA",
                          MaterialKind.MOLECULAR,
                          state=INSERT,
                          source=SOURCE + " (exact terminal PCR-primer bases; transcript schematic)")
    g.react("pcr", "Amplify and barcode full-length cDNA (12 cycles for SMRTbell prep kit 3.0)",
            StepKind.PCR,
            ((full_clean, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
             (rev, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (amplified,),
            source=SOURCE + " (102-396-000 Rev10 section 3.1: 98 °C 45 s; 12 × "
                            "[98 °C 10 s, 60 °C 15 s, 72 °C 3 min]; 72 °C 5 min)",
            geometry=PCRGeometry.unavailable(
                full_clean.state,
                reason=("The official FASTA supplies complete PCR primers, while the actual starting cDNA precedes introduction of the barcode primer's 5′ portion. "
                        "The current model has not reconstructed the RT-primer/TSO-derived precursor ends or the initial strand-copying phase. "
                        "An already barcoded duplex product cannot establish those starting binding sites.")))
    pcr_beads = g.input("pcr_beads", "1.1X (or 0.9X) SMRTbell cleanup beads",
                        MaterialKind.CONDITION)
    amplified_clean = g.product("amplified_clean", "size-enriched barcoded cDNA (pooled when multiplexed)",
                                MaterialKind.MOLECULAR, state=INSERT)
    _bead_cleanup(g, "pcr_cleanup", "Purify amplified cDNA with 1.1X (or 0.9X) beads",
                  amplified, pcr_beads, amplified_clean,
                  "102-396-000 Rev10 sections 3.2–3.3: 45 µL PCR + 49.5 µL (1.1X) or "
                  "40.5 µL (0.9X, longer cDNA) beads, 47 µL elution for SMRTbell prep "
                  "kit 3.0; ≥100 ng per sample; equal-mass pooling up to 500 ng is "
                  "not drawn as a separate step.")
    input_beads = g.input("spk3_input_beads", "1.3X SMRTbell beads for Rev08 amplicon input cleanup", MaterialKind.CONDITION)
    spk3_input = g.product("spk3_input", "barcoded cDNA after Rev08 input cleanup", MaterialKind.MOLECULAR, state=INSERT)
    _bead_cleanup(g, "spk3_input_cleanup", "Clean up amplicon input before Rev08 repair", amplified_clean, input_beads, spk3_input,
                  "102-359-000 Rev08 section1.5: 1.3X bead cleanup before section2 repair/A-tailing.")
    end_mix = g.input("end_mix", "Repair buffer, End repair mix and DNA repair mix",
                      MaterialKind.CATALYST, source=SOURCE + " (" + SPK3 + ")")
    tailed = g.product("tailed", "repaired dA-tailed full-length cDNA",
                       MaterialKind.MOLECULAR)
    g.react("end_prep", "Repair and A-tail amplified cDNA (37 °C 30 min, 65 °C 5 min)",
            StepKind.END_REPAIR,
            ((spk3_input, InputRole.TEMPLATE), (end_mix, InputRole.CATALYST)), (tailed,),
            source=SOURCE + " (" + SPK3 + ")")
    hairpins = g.input("hairpins", "SMRTbell adapter; ligating end unresolved",
                       MaterialKind.ADAPTER, inferred=True,
                       source=SOURCE + " (adapter bases not printed in inspected sources)")
    ligase = g.input("ligase", "Ligation mix + Ligation enhancer", MaterialKind.CATALYST,
                     source=SOURCE + " (" + SPK3 + ")")
    closed = g.product("closed", "ligation mixture; a closed SMRTbell candidate is shown",
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
    digested = g.product("digested", "nuclease-treated ligation mixture; a closed candidate is shown",
                         MaterialKind.MOLECULAR, state=tuple(SMRTBELL_ROWS), inferred=True,
                         source=SOURCE + " (hairpin bases not verified in inspected sources)")
    other_nuclease_products = g.product("nuclease_other_products", "other nuclease-treatment products; identities and yields unresolved", MaterialKind.MOLECULAR)
    g.react("remove_open", "Treat ligation products with the source nuclease mix (37 °C 15 min)",
            StepKind.DIGESTION,
            ((closed_clean, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
            (digested, other_nuclease_products),
            source=SOURCE + " (" + SPK3 + "; the substrate specificity of the nuclease "
                            "mix is not stated in the inspected procedures)")
    g.mark_terminal(Terminal(other_nuclease_products.id, Disposition.UNRESOLVED,
        "The inspected guide names Nuclease treatment without specifying substrate selectivity or individual products. "
        "Their later 1.3X bead partition is not separately reconstructed; hydrolysis is not a physical discard."))
    nuc_beads = g.input("nuc_beads", "1.3X SMRTbell cleanup beads (post-nuclease)",
                        MaterialKind.CONDITION)
    library = g.product("library", "closed full-length Iso-Seq SMRTbell library",
                        MaterialKind.LIBRARY, state=tuple(SMRTBELL_ROWS), inferred=True,
                        source=SOURCE + " (hairpin bases not verified in inspected sources)")
    _bead_cleanup(g, "nuc_cleanup", "Purify nuclease-treated library with 1.3X beads",
                  digested, nuc_beads, library,
                  "102-359-000 Rev08 section 4: 65 µL beads to 50 µL, 26 µL elution; "
                  "bulk supernatant is discarded; no species-specific recovery is assigned to individual nuclease products.")
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, acquired_roles=("Polymerase read",),
        read_uses={"Polymerase read": "include HiFi consensus from repeated circuits"},
        limitation=SEQUENCING_UNAVAILABLE))
    g.declare_reviewed_scope(("rna", "rt_primer", "rt", "reverse_transcribe", "tso",
                              "template_switch", "rt_beads", "rt_cleanup", "full_clean",
                              "forward", "reverse", "polymerase", "pcr", "pcr_beads",
                              "pcr_cleanup", "amplified_clean", "spk3_input_beads", "spk3_input_cleanup", "end_mix", "end_prep",
                              "hairpins", "ligase", "ligate", "lig_beads", "lig_cleanup",
                              "closed_clean", "nuclease", "remove_open", "nuc_beads",
                              "nuc_cleanup", "library"), SOURCE)
    return g.require_publishable()
