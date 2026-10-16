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


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("rna", "poly(A) RNA", MaterialKind.BIOLOGICAL,
                      primary=True)
    rt_primer = graph.input(
        "rt_primer", "poly(A)-primed Iso-Seq RT primer mix", MaterialKind.OLIGO,
        inferred=True, source=SOURCE + " (primer bases unverified)")
    reverse_transcriptase = graph.input(
        "reverse_transcriptase", "template-switching reverse transcriptase",
        MaterialKind.CATALYST)
    first = graph.product("first_strand", "full-length first-strand cDNA",
                          MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe", "Reverse-transcribe poly(A) RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first,))
    tso = graph.input(
        "tso", "Kinnex full-length-RNA template-switch oligo", MaterialKind.OLIGO,
        inferred=True, source=SOURCE + " (oligo bases unverified)")
    full_length = graph.product(
        "full_length_cdna", "full-length cDNA with terminal PCR handles",
        MaterialKind.MOLECULAR)
    graph.react("template_switch", "Capture the completed cDNA 5-prime end",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (full_length,))
    isoseq_forward = graph.input(
        "isoseq_forward", "Iso-Seq barcoded cDNA amplification primer (selected member unassigned)",
        MaterialKind.OLIGO, source=SOURCE + " (kit reagent role; analysis FASTA is not a full physical oligo specification)")
    isoseq_reverse = graph.input(
        "isoseq_reverse", "Iso-Seq cDNA amplification primer (physical oligo sequence unassigned)",
        MaterialKind.OLIGO, source=SOURCE + " (kit reagent role; analysis FASTA is not a full physical oligo specification)")
    polymerase = graph.input("pcr_polymerase", "Kinnex PCR polymerase",
                             MaterialKind.CATALYST)
    amplified = graph.product(
        "amplified_cdna", "barcoded full-length cDNA pool", MaterialKind.MOLECULAR,
        state=INSERTS[0], inferred=True, source=SOURCE + " (selected barcode and completed cDNA ends not specified)")
    graph.react("isoseq_pcr", "Amplify and barcode full-length cDNA", StepKind.PCR,
                ((full_length, InputRole.TEMPLATE),
                 (isoseq_forward, InputRole.STRUCTURAL),
                 (isoseq_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    INSERTS[0], reason="The source names the cDNA amplification reagents. Primer-trimming references do not establish complete physical PCR oligos or their binding geometry in this declaration."))

    aliquots = tuple(graph.product(
        f"aliquot_{i}", f"Kinnex orientation-{i} PCR aliquot",
        MaterialKind.MOLECULAR, state=INSERTS[i - 1]) for i in range(1, 9))
    graph.react("split_eight", "Divide cDNA into eight orientation-specific PCR reactions",
                StepKind.SPLIT, ((amplified, InputRole.TEMPLATE),), aliquots,
                changes_molecule=False)
    segmentation_forward = graph.input(
        "segmentation_forward", "eight orientation-specific Kinnex forward primers",
        MaterialKind.OLIGO, inferred=True,
        source=SOURCE + " (segmentation-primer bases unverified)")
    segmentation_reverse = graph.input(
        "segmentation_reverse", "eight orientation-specific Kinnex reverse primers",
        MaterialKind.OLIGO, inferred=True,
        source=SOURCE + " (segmentation-primer bases unverified)")
    oriented = []
    reviewed = [
        "rna", "rt_primer", "reverse_transcriptase", "reverse_transcribe", "tso",
        "template_switch", "isoseq_forward", "isoseq_reverse", "pcr_polymerase",
        "isoseq_pcr", "split_eight", "segmentation_forward", "segmentation_reverse",
    ]
    for i, aliquot in enumerate(aliquots, 1):
        product = graph.product(
            f"oriented_{i}", f"orientation-specific Kinnex PCR product {i}",
            MaterialKind.MOLECULAR, state=INSERTS[i - 1], inferred=True,
            source=SOURCE + " (segmentation-junction bases unverified)")
        step_id = f"orientation_pcr_{i}"
        graph.react(
            step_id, f"INFERRED — install orientation-{i} segmentation ends by PCR", StepKind.PCR,
            ((aliquot, InputRole.TEMPLATE),
             (segmentation_forward, InputRole.STRUCTURAL),
             (segmentation_reverse, InputRole.STRUCTURAL),
             (polymerase, InputRole.CATALYST)), (product,),
            geometry=PCRGeometry.unavailable(
                INSERTS[i - 1],
                reason="The kit fixes the eight reaction identities; their full physical primer sequences are not established in this declaration."))
        oriented.append(product)
        reviewed.append(step_id)
    pooled = graph.product("pooled_products", "equal-volume pool of eight oriented products",
                           MaterialKind.MOLECULAR, state=tuple(ARRAY.rows()))
    graph.react("pool_eight", "Pool equal volumes of all eight oriented products",
                StepKind.POOL, tuple((x, InputRole.TEMPLATE) for x in oriented), (pooled,),
                changes_molecule=False)
    reviewed.append("pool_eight")

    assembly_components = graph.input(
        "assembly_components", "segmentation junctions and barcoded terminal adapters",
        MaterialKind.ADAPTER, inferred=True,
        source=SOURCE + " (junction and terminal-adapter bases unverified)")
    assembly_enzyme = graph.input("assembly_enzyme", "Kinnex array-formation enzyme",
                                  MaterialKind.CATALYST)
    ligase = graph.input("array_ligase", "Kinnex array ligase", MaterialKind.CATALYST)
    closed = graph.product(
        "closed_array", "closed ordered eight-insert Kinnex SMRTbell",
        MaterialKind.MOLECULAR, state=tuple(dumbbell_rows(SMRTBELL)), inferred=True,
        source=SOURCE + " (junction and terminal-adapter bases unverified)")
    graph.react(
        "assemble_array", "INFERRED — assemble the ordered array and close its termini",
        StepKind.LIGATION,
        ((pooled, InputRole.TEMPLATE), (assembly_components, InputRole.STRUCTURAL),
         (assembly_enzyme, InputRole.CATALYST), (ligase, InputRole.CATALYST)), (closed,),
        geometry=AdapterLigationGeometry.unavailable(
            "The kit specifies array order and covalent closure. Full physical junction and terminal-adapter geometry is not established in this declaration; no blunt or cohesive end is assigned."))
    repair = graph.input("repair", "Kinnex DNA-repair mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "repaired closed Kinnex arrays",
                             MaterialKind.MOLECULAR, state=tuple(dumbbell_rows(SMRTBELL)),
                             inferred=True, source=closed.source)
    graph.react("repair_array", "Repair assembled Kinnex arrays", StepKind.END_REPAIR,
                ((closed, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    nuclease = graph.input("nuclease", "free-end nuclease", MaterialKind.CATALYST)
    library = graph.product("library", SMRTBELL.name, MaterialKind.LIBRARY,
                            state=tuple(dumbbell_rows(SMRTBELL)), inferred=True,
                            source=closed.source)
    incomplete = graph.product("incomplete_arrays", "arrays retaining exposed DNA ends",
                               MaterialKind.WASTE)
    graph.react("remove_incomplete", "Digest arrays retaining exposed ends",
                StepKind.DIGESTION,
                ((repaired, InputRole.TEMPLATE), (nuclease, InputRole.CATALYST)),
                (library, incomplete), changes_molecule=False)
    graph.discard(incomplete, "does not form a completely closed Kinnex SMRTbell")
    graph.sequence(library, smrtbell_sequencing_plan(
        SMRTBELL, smrtbell_circuit(SMRTBELL),
        limitation="The Kinnex kit names the PacBio sequencing primer but does not print its terminal-adapter binding sequence in the checked guide."))
    reviewed.extend((
        "assembly_components", "assembly_enzyme", "array_ligase", "assemble_array",
        "repair", "repair_array", "nuclease", "remove_incomplete", "library"))
    graph.declare_reviewed_scope(reviewed, SOURCE)
    return graph.require_publishable()
