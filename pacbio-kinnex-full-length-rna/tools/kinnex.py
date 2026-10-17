"""PacBio Kinnex full-length RNA ordered-array workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature
from concatemer import SegmentedArray
from adapter_ligation import AdapterLigationGeometry
from dumbbell import (Dumbbell, dumbbell_rows, smrtbell_circuit,
                      smrtbell_sequencing_plan)
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, StepKind)


def role(name: str, width: int, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("inferred", True)
    kw.setdefault("note", "Schematic placeholder width; physical sequence length is not assigned.")
    return Segment(name, "X" * width, tag, placeholder=True, **kw)


TITLE = "PacBio Kinnex full-length RNA — eight-insert arrays"
NOTES = "01_kinnex-full-length-rna.html"
SOURCE = ('Commercial defining source: PacBio <a href="https://www.pacb.com/wp-content/'
          'uploads/Procedure-checklist-Preparing-Kinnex-libraries-using-the-Kinnex-full-'
          'length-RNA-kit.pdf">Preparing Kinnex libraries using the Kinnex full-length '
          'RNA kit, 103-238-700 Rev09</a>.')
SUMMARY = ("Full-length cDNA is copied in eight parallel orientation-specific PCRs; the "
           "products assemble into ordered eight-insert arrays whose terminal adapters "
           "form nuclease-resistant SMRTbell templates for HiFi sequencing.")
CAVEAT = ("PacBio's Kinnex index FASTA confirms the four six-base terminal-adapter barcodes. "
          "The official MAS8 reference FASTA gives nine short segmentation motifs (A–H and Q), "
          "and the SMRT Link guide prints Iso-Seq v2 cDNA primer-trimming references. "
          "These analysis references do not establish complete Kinnex PCR primer constructs, "
          "ligation-junction geometry or terminal-adapter arms. Their order and topology remain "
          "shown as named unknown DNA where the physical construct is unresolved. "
          "The 3-prime trimming reference excludes poly(T) and is reverse-complemented; "
          "it is not an oligo ordering sequence. Placeholder widths are schematic, not physical lengths.")

KINNEX_TERMINAL_BARCODES = {"bcM0001": "ACAGTC", "bcM0002": "ATGACG",
                            "bcM0003": "CACGTG", "bcM0004": "CATCGC"}
# Analysis-reference sequences, not a physical PCR-primer ordering declaration.
# SMRT Link v13.1 p100 explicitly excludes Ts and reverse-complements the 3p entry.
ISOSEQ_V2_BC01_5P_REFERENCE = "CTACACGACGCTCTTCCGATCTACTACACGCAATGAAGTCGCAGGGTTGGG"
ISOSEQ_V2_3P_REFERENCE = "AAGCAGTGGTATCAACGCAGAGTAC"

# PacBio MAS-Seq Adapter v3 (MAS8) analysis reference, mas8_primers.fasta.
# These are segmentation recognition motifs, not complete physical PCR oligos.
MAS8_REFERENCE_MOTIFS = {
    "A": "AGCTTACTTGTGAAGAT", "B": "ACTTGTAAGCTGTCTAT",
    "C": "ACTCTGTCAGGTCCGAT", "D": "ACCTCCTCCTCCAGAAT",
    "E": "AACCGGACACACTTAGT", "F": "AGAGTCCAATTCGCAGT",
    "G": "AATCAAGGCTTAACGGT", "H": "ATGTTGAATCCTAGCGT",
    "Q": "AGTAGCTGTGTGCA",
}

def oligos():
    return [f"<b>Kinnex terminal-adapter barcode {name}</b> — 5′-{bases}-3′ "
            "(barcode only; adapter body unresolved)"
            for name, bases in KINNEX_TERMINAL_BARCODES.items()]

ISOSEQ_INDEX = feature("isoseq_sample_index", "sample_index", "whitelist",
                       whitelist="PacBio Iso-Seq primer barcodes bc01–12")
INSERTS = tuple(Construct([
                    Segment("Iso-Seq primer barcode", "[Iso-Seq barcode]", "cbc",
                            placeholder=True, inferred=True, feature=ISOSEQ_INDEX),
                    role("full-length cDNA", 18, "r1")],
                          name=f"orientation-specific Kinnex product {i}")
                for i in range(1, 9))
JUNCTIONS = tuple(role(f"segmentation junction {i}|{i + 1}", 8, "r3")
                  for i in range(1, 8))
ARRAY = SegmentedArray(INSERTS, JUNCTIONS,
                       role("left terminal adapter", 12, "r2"),
                       role("right terminal adapter", 12, "r2"),
                       name="eight-insert Kinnex array")
SMRTBELL = Dumbbell(ARRAY.core(), "N" * 18, "N" * 18,
                    name="Kinnex array SMRTbell", insert_overhang="",
                    adapter_overhang="", junction_unavailable=True)
KINNEX_PCR_REACTIONS = 8
KINNEX_PRIMER_MIXES = ("A", "B", "C", "D", "E", "F", "G", "HQ")  # 103-238-700 Rev09 p6, step 4.1
ISOSEQ_PCR_CYCLES = 10   # 103-238-700 Rev09 step 3.1.5
KINNEX_PCR_CYCLES = 9    # 103-238-700 Rev09 step 4.1 PCR program
ARRAY_FORMATION = (45, 60)
DNA_REPAIR = (45, 30)
NUCLEASE = (37, 15)


def polymerase_cycle() -> Construct:
    core = ARRAY.core()
    return Construct([
        role("eight cDNA inserts and seven segmentation junctions, forward", len(core), "r1"),
        role("right terminal adapter / primer site", 18, "r2"),
        role("eight cDNA inserts and seven segmentation junctions, reverse", len(core), "r1"),
        role("left terminal adapter / primer site", 18, "r2"),
    ], name="one Kinnex SMRTbell polymerase circuit")


def _bead_cleanup(graph, identifier, label, template, beads, kept, reason):
    """One SMRTbell bead cleanup: the bound fraction is kept, the supernatant discarded."""
    waste = graph.product(identifier + "_supernatant",
                          f"supernatant not retained by {beads.label}", MaterialKind.WASTE)
    graph.react(identifier, label, StepKind.CLEANUP,
                ((template, InputRole.TEMPLATE), (beads, InputRole.CONDITION)),
                (kept, waste), changes_molecule=False, reason=reason)
    graph.discard(waste, f"supernatant discarded after {beads.label} binding")
    return kept


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("rna", "300 ng total RNA (RIN ≥7) containing poly(A) RNA",
                      MaterialKind.BIOLOGICAL, primary=True)
    rt_primer = graph.input(
        "rt_primer", "poly(A)-primed Iso-Seq RT primer mix", MaterialKind.OLIGO,
        inferred=True, source=SOURCE + " (primer bases unverified)")
    reverse_transcriptase = graph.input(
        "reverse_transcriptase", "Iso-Seq RT enzyme mix (template-switching reverse transcriptase)",
        MaterialKind.CATALYST)
    first = graph.product("first_strand", "full-length first-strand cDNA",
                          MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe",
                "Anneal RT primer mix (70 °C 5 min), reverse-transcribe 42 °C 45 min",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first,))
    tso = graph.input(
        "tso", "Iso-Seq Express TSO 2.0", MaterialKind.OLIGO,
        inferred=True, source=SOURCE + " (oligo bases unverified)")
    full_length = graph.product(
        "full_length_cdna", "full-length cDNA with terminal PCR handles",
        MaterialKind.MOLECULAR)
    graph.react("template_switch", "Add TSO 2.0 and template-switch 42 °C 15 min",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (full_length,))
    rt_beads = graph.input("rt_beads", "1.3X SMRTbell cleanup beads", MaterialKind.CONDITION)
    full_clean = graph.product("full_length_clean", "purified full-length cDNA",
                               MaterialKind.MOLECULAR)
    _bead_cleanup(graph, "rt_cleanup", "Purify RT/template-switch reaction with 1.3X beads",
                  full_length, rt_beads, full_clean,
                  "103-238-700 Rev09 section 2.4: 29 µL elution buffer + 65 µL beads, 21 µL eluate.")
    isoseq_forward = graph.input(
        "isoseq_forward", "Iso-Seq barcoded cDNA amplification primer (selected member unassigned)",
        MaterialKind.OLIGO, source=SOURCE + " (kit reagent role; analysis FASTA is not a full physical oligo specification)")
    isoseq_reverse = graph.input(
        "isoseq_reverse", "Iso-Seq cDNA amplification primer (physical oligo sequence unassigned)",
        MaterialKind.OLIGO, source=SOURCE + " (kit reagent role; analysis FASTA is not a full physical oligo specification)")
    cdna_pcr_mix = graph.input("cdna_pcr_mix", "Iso-Seq cDNA PCR mix", MaterialKind.CATALYST)
    polymerase = graph.input("pcr_polymerase", "Kinnex PCR mix", MaterialKind.CATALYST)
    amplified = graph.product(
        "amplified_cdna", "barcoded full-length cDNA pool", MaterialKind.MOLECULAR,
        state=INSERTS[0], inferred=True, source=SOURCE + " (selected barcode and completed cDNA ends not specified)")
    graph.react("isoseq_pcr", f"Amplify and barcode full-length cDNA ({ISOSEQ_PCR_CYCLES} cycles)",
                StepKind.PCR,
                ((full_clean, InputRole.TEMPLATE),
                 (isoseq_forward, InputRole.STRUCTURAL),
                 (isoseq_reverse, InputRole.STRUCTURAL),
                 (cdna_pcr_mix, InputRole.CATALYST)), (amplified,),
                source=SOURCE + " (step 3.1: 98 °C 45 s; 10 × [98 °C 10 s, 60 °C 15 s, 72 °C 3 min]; 72 °C 5 min)",
                geometry=PCRGeometry.unavailable(
                    INSERTS[0], reason="The source names the cDNA amplification reagents. Primer-trimming references do not establish complete physical PCR oligos or their binding geometry in this declaration."))

    cdna_beads = graph.input("cdna_beads", "0.9X (or 0.86X) SMRTbell cleanup beads",
                             MaterialKind.CONDITION)
    cdna_clean = graph.product("amplified_cdna_clean",
                               "size-enriched barcoded cDNA (55 ng, equal-mass pool when multiplexed)",
                               MaterialKind.MOLECULAR, state=INSERTS[0], inferred=True,
                               source=amplified.source)
    _bead_cleanup(graph, "cdna_cleanup", "Purify amplified cDNA with 0.9X beads",
                  amplified, cdna_beads, cdna_clean,
                  "103-238-700 Rev09 sections 3.2–3.3: 45 µL PCR + 55 µL elution buffer + "
                  "90 µL beads (86 µL for longer cDNA), 24 µL eluate; ≥25 ng (55 ng "
                  "recommended) enters Kinnex PCR; equal-mass pooling of barcoded cDNA "
                  "to 55 ng is not drawn as a separate step.")
    aliquots = tuple(graph.product(
        f"aliquot_{i}", f"Kinnex orientation-{i} PCR aliquot",
        MaterialKind.MOLECULAR, state=INSERTS[i - 1]) for i in range(1, 9))
    graph.react("split_eight",
                "Make a Kinnex PCR master mix with 55 ng cDNA and aliquot 22.5 µL into eight tubes",
                StepKind.SPLIT, ((cdna_clean, InputRole.TEMPLATE),), aliquots,
                changes_molecule=False, source=SOURCE + " (steps 4.1.4–4.1.5)")
    # Each PCR receives one named primer mix; the mixes' composition is not printed.
    # Splitting a mix into upstream/downstream primers is an inferred PCR requirement.
    segmentation_forward = graph.input(
        "segmentation_forward",
        "Kinnex primer mixes A–G, HQ: upstream-priming members (inferred; composition not printed)",
        MaterialKind.OLIGO, inferred=True,
        source=SOURCE + " (p6 reagent list and step 4.1 name primer mixes A, B, C, D, E, F, G and HQ; primer bases unverified)")
    segmentation_reverse = graph.input(
        "segmentation_reverse",
        "Kinnex primer mixes A–G, HQ: downstream-priming members (inferred; composition not printed)",
        MaterialKind.OLIGO, inferred=True,
        source=SOURCE + " (p6 reagent list and step 4.1 name primer mixes A, B, C, D, E, F, G and HQ; primer bases unverified)")
    oriented = []
    reviewed = [
        "rna", "rt_primer", "reverse_transcriptase", "reverse_transcribe", "tso",
        "template_switch", "isoseq_forward", "isoseq_reverse", "pcr_polymerase",
        "isoseq_pcr", "split_eight", "rt_beads", "rt_cleanup", "full_length_clean",
        "cdna_pcr_mix", "cdna_beads", "cdna_cleanup", "amplified_cdna_clean",
        "segmentation_forward", "segmentation_reverse",
    ]
    for i, aliquot in enumerate(aliquots, 1):
        product = graph.product(
            f"oriented_{i}", f"orientation-specific Kinnex PCR product {i}",
            MaterialKind.MOLECULAR, state=INSERTS[i - 1], inferred=True,
            source=SOURCE + " (segmentation-junction bases unverified)")
        step_id = f"orientation_pcr_{i}"
        mix = KINNEX_PRIMER_MIXES[i - 1]
        graph.react(
            step_id, f"INFERRED — install orientation-{i} segmentation ends by PCR with primer mix {mix} "
                     f"({KINNEX_PCR_CYCLES} cycles)", StepKind.PCR,
            ((aliquot, InputRole.TEMPLATE),
             (segmentation_forward, InputRole.STRUCTURAL),
             (segmentation_reverse, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (product,),
            source=SOURCE + " (step 4.1: 98 °C 3 min; 9 × [98 °C 20 s, 68 °C 30 s, 72 °C 4 min]; 72 °C 5 min)",
            geometry=PCRGeometry.unavailable(
                INSERTS[i - 1],
                reason="The kit fixes the eight primer-mix identities (A–G, HQ); their primer composition and sequences are not printed."))
        oriented.append(product)
        reviewed.append(step_id)
    pooled = graph.product("pooled_products", "equal-volume pool of eight oriented products",
                           MaterialKind.MOLECULAR, state=tuple(ARRAY.rows()))
    graph.react("pool_eight", "Pool exactly 23 µL of each of the eight oriented products",
                StepKind.POOL, tuple((x, InputRole.TEMPLATE) for x in oriented), (pooled,),
                changes_molecule=False, source=SOURCE + " (step 4.2)")
    pool_beads = graph.input("pool_beads", "1.05X SMRTbell cleanup beads", MaterialKind.CONDITION)
    pooled_clean = graph.product("pooled_clean", "purified pooled Kinnex PCR products",
                                 MaterialKind.MOLECULAR, state=tuple(ARRAY.rows()))
    _bead_cleanup(graph, "pool_cleanup", "Purify the pool with 1.05X beads", pooled,
                  pool_beads, pooled_clean,
                  "103-238-700 Rev09 step 4.2: 193 µL beads to 184 µL, 40 µL elution; "
                  "typical yield 6–18 µg, 5 µg enters array formation.")
    reviewed.extend(("pool_eight", "pool_beads", "pool_cleanup", "pooled_clean"))

    assembly_components = graph.input(
        "assembly_components", "Kinnex adapter bc01–04 mix (one per library); junction ends from the primer mixes",
        MaterialKind.ADAPTER, inferred=True,
        source=SOURCE + " (step 5.1; junction and terminal-adapter bases unverified)")
    assembly_enzyme = graph.input("assembly_enzyme", "Kinnex enzyme in Kinnex array and repair buffer",
                                  MaterialKind.CATALYST)
    ligase = graph.input("array_ligase", "Kinnex ligase", MaterialKind.CATALYST)
    closed = graph.product(
        "closed_array", "Kinnex arrays with terminal adapters (full and partial arrays)",
        MaterialKind.MOLECULAR, state=tuple(dumbbell_rows(SMRTBELL)), inferred=True,
        source=SOURCE + " (junction and terminal-adapter bases unverified)")
    graph.react(
        "assemble_array", "INFERRED — 'Kinnex primer digestion/ligation' (45 °C 60 min) assembles the ordered array and ligates terminal adapters",
        StepKind.LIGATION,
        ((pooled_clean, InputRole.TEMPLATE), (assembly_components, InputRole.STRUCTURAL),
         (assembly_enzyme, InputRole.CATALYST), (ligase, InputRole.CATALYST)), (closed,),
        geometry=AdapterLigationGeometry.unavailable(
            "Step 5.1 names a 'Kinnex primer digestion/ligation' program that assembles cDNA segments into a linear array with a barcoded terminal adapter; the enzyme's activity, junction and terminal-adapter geometry are not printed. No blunt or cohesive end is assigned and no research-paper USER mechanism is imported."))
    repair = graph.input("repair", "DNA repair mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "repaired Kinnex arrays",
                             MaterialKind.MOLECULAR, state=tuple(dumbbell_rows(SMRTBELL)),
                             inferred=True, source=closed.source)
    graph.react("repair_array", "Add DNA repair mix directly and repair 45 °C 30 min",
                StepKind.END_REPAIR,
                ((closed, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,),
                source=SOURCE + " (step 5.1 DNA damage repair program)")
    array_beads = graph.input("array_beads", "1X SMRTbell cleanup beads (post-array)",
                              MaterialKind.CONDITION)
    repaired_clean = graph.product("repaired_clean", "purified Kinnex arrays",
                                   MaterialKind.MOLECULAR, state=tuple(dumbbell_rows(SMRTBELL)),
                                   inferred=True, source=closed.source)
    _bead_cleanup(graph, "array_cleanup", "Purify arrays with 1X beads", repaired,
                  array_beads, repaired_clean, "103-238-700 Rev09 section 5.2: 60 µL beads.")
    nuclease = graph.input("nuclease", "Nuclease buffer + Nuclease mix", MaterialKind.CATALYST)
    digested = graph.product("digested", "nuclease-resistant full arrays and digested partial arrays",
                             MaterialKind.MOLECULAR, state=tuple(dumbbell_rows(SMRTBELL)),
                             inferred=True, source=closed.source)
    graph.react("remove_incomplete", "Nuclease treatment 37 °C 15 min digests partial arrays",
                StepKind.DIGESTION,
                ((repaired_clean, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
                (digested,), changes_molecule=False,
                source=SOURCE + " (section 6.1); application note 102-326-591 Rev03 p3: "
                                "partial arrays are removed during nuclease digestion")
    final_beads = graph.input("final_beads", "1X SMRTbell cleanup beads (final)",
                              MaterialKind.CONDITION)
    library = graph.product("library", SMRTBELL.name, MaterialKind.LIBRARY,
                            state=tuple(dumbbell_rows(SMRTBELL)), inferred=True,
                            source=closed.source)
    _bead_cleanup(graph, "final_cleanup", "Final 1X bead cleanup of the nuclease-treated library",
                  digested, final_beads, library,
                  "103-238-700 Rev09 section 6.2: 50 µL beads to 50 µL; expect 10–25% "
                  "recovery of the Kinnex PCR product; digestion fragments leave with the supernatant.")
    graph.sequence(library, smrtbell_sequencing_plan(
        SMRTBELL, smrtbell_circuit(SMRTBELL),
        limitation="Section 7 anneals the Kinnex sequencing primer, but no inspected source prints its bases or its terminal-adapter binding site."))
    reviewed.extend((
        "assembly_components", "assembly_enzyme", "array_ligase", "assemble_array",
        "repair", "repair_array", "array_beads", "array_cleanup", "repaired_clean",
        "nuclease", "remove_incomplete", "digested", "final_beads", "final_cleanup",
        "library"))
    graph.declare_reviewed_scope(reviewed, SOURCE)
    return graph.require_publishable()
