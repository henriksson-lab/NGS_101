from chemdraw import Workflow
from types import SimpleNamespace
from rna_cell_assays import Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, single, sp, token
from rna_cell_assays import nextera_library

from rna_cell_assays import chain
from protocol_graph import AnnotationSlot, SourceClaim, SourceDiscrepancy

def find_seq_workflow() -> Workflow:
    return chain("heterogeneous cell suspension", single("one cell", (token("genome"),
                                                                        token("poly(A) RNA"))), (
        ("Co-encapsulate one cell, molten allyl agarose, and acrydite oligo(dT)",
         single("hydrogel bead precursor", (token("cell genome"), token("cell RNA"),
                                             token("covalent oligo(dT) matrix"))),
         "APS/TEMED links acrydite oligo(dT) to allyl agarose."),
        ("Lyse strongly, wash cold, then reverse-transcribe and template-switch on the hydrogel",
         single("genome-capture bead", (token("retained genome"), token("matrix-bound cDNA"))),
         "LiDS/Proteinase K lysis at 55 °C; cold high-salt washes; Maxima H Minus with the TSO "
         "(Smart-seq3 TSO in the demonstrated run), 30 min RT then 90 min 42 °C, EDTA stop."),
        ("Re-inject beads with multiplex TaqMan PCR, thermocycle, and sort double-positive droplets",
         single("fluorescent selected bead", (token("retained genome"), token("matrix-bound cDNA"),
                                               token("TaqMan marker signal"))),
         "Demonstrated assay: HIV IPDA Ψ (FAM) and env (VIC) probes, 55 cycles with 88 °C denaturation. "
         "Marker PCR is a selection measurement; it is not the transcriptome library."),
        ("Melt sorted droplets at 60 °C, amplify cDNA with Smart-seq3 WTA primers, then Nextera XT",
         single("FIND-seq library", (token("Nextera left arm"), token("selected-cell cDNA"),
                                      token("Nextera right arm"))),
         "Bulk bins or single sorted droplets; KAPA HiFi WTA, 0.8× AMPure XP, Nextera XT per manufacturer."),
    ))

DEFINITION = Definition("FIND-seq — nucleic-acid cytometry for rare cells", "10.1038/s41596-024-01021-y",
        "Shin et al., <i>Nature Protocols</i> (2024)", "01_find-seq.html",
        "A hydrogel bead retains one cell's genome and oligo(dT)-captured transcriptome while marker-specific TaqMan PCR drives fluorescent droplet sorting before whole-transcriptome library preparation.",
        "Sorted-droplet cDNA is amplified with the published Smart-seq3 WTA primers (Smart-seq2 oligos are also listed) and tagmented with the Nextera XT kit; the protocol does not state a sequencing run or read layout, so the standard Nextera ends are shown without acquired cycles.")

d = DEFINITION
final, primers = nextera_library([token('selected-cell transcript cDNA', length=44)], 'FIND-seq transcriptome library')
reads = None
OLIGOS = (
            "<b>Smart-seq2 bead capture</b> — 5′-acrydite-TTTTTTAAGCAGTGGTATCAACGCAGAGTAC(T)₃₀VN-3′",
            "<b>Smart-seq2 template switch</b> — 5′-AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG-3′",
            "<b>Smart-seq3 bead capture</b> — 5′-acrydite-TTTTTTACGAGCATCAGCAGCATACGA(T)₃₀V-3′",
            "<b>Smart-seq2 Smart PCR primer</b> — 5′-AAGCAGTGGTATCAACGCAGAGT-3′",
            "<b>Smart-seq3 template switch</b> — 5′-biotin-AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG-3′; N₈ is the UMI",
            "<b>Smart-seq3 forward (WTA) primer</b> — 5′-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G-3′; * phosphorothioate",
            "<b>Smart-seq3 reverse (WTA) primer</b> — 5′-ACGAGCATCAGCAGCATAC*G*-3′; * phosphorothioate",
            "<b>env TaqMan (selection only)</b> — F 5′-AGTGGTGCAGAGAGAAAAAAGAGC-3′; R 5′-GTCTGGCCTGTACCGTCAGC-3′; probe 5′-VIC-CCTTGGGTTCTTGGGA-MGBNHQ-3′ (quencher label as printed)",
            "<b>Ψ TaqMan (selection only)</b> — F 5′-CAGGACTCGGCTTGCTGAAG-3′; R 5′-GCACCCATCTCTCTCCTTCTAGC-3′; probe 5′-FAM-TTTTGGCGTACTCACCAGT-MGBNHQ-3′ (quencher label as printed)",
        )
workflow_fn = find_seq_workflow
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
    SEQ_PRIMERS=primers, FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=reads, SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("" if primers else "Exact sequencing-primer sequences are not fixed by the defining protocol or were not exposed in its accessible primary text; read entry and cycle roles are stated explicitly."),
    workflow=workflow_fn, oligos=lambda: OLIGOS,
)


def find_seq_graph(definition, final: Construct,
                   primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """FIND-seq graph from hydrogel partition through selected-cell library."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    cells = graph.input("cell_suspension", "heterogeneous cell suspension",
                        MaterialKind.BIOLOGICAL, state=find_seq_workflow().initial, primary=True)
    matrix = graph.input("allyl_agarose", "molten allyl agarose matrix",
                         MaterialKind.CONDITION)
    capture_oligo = graph.input("acrydite_oligodt",
                                "acrydited oligo(dT) primer (Smart-seq3 or Smart-seq2 form)",
                                MaterialKind.OLIGO)
    bead = graph.product("cell_hydrogel", "one-cell hydrogel bead precursor",
                         MaterialKind.MOLECULAR,
                         state=single("hydrogel bead precursor", (
                             token("cell genome"), token("cell RNA"),
                             token("covalent oligo(dT) matrix"))))
    empty = graph.product("empty_partition_material", "empty or multiplet partitions",
                          MaterialKind.WASTE)
    graph.react("encapsulate", "Encapsulate cells in oligo(dT) hydrogels",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (matrix, InputRole.CONDITION),
                 (capture_oligo, InputRole.STRUCTURAL)), (bead, empty))
    graph.discard(empty, "partitions without one usable cell are not selected for sequencing")
    proteinase = graph.input("proteinase_k", "proteinase K", MaterialKind.CATALYST)
    lids = graph.input("lids", "LiDS lysis buffer, 55 °C", MaterialKind.CONDITION)
    lysed = graph.product("lysed_bead", "lysed genome-retaining hydrogel bead",
                          MaterialKind.MOLECULAR, state=bead.state)
    graph.react("lyse", "Lyse cells strongly inside the hydrogel",
                StepKind.DEPROTEINIZATION,
                ((bead, InputRole.TEMPLATE), (proteinase, InputRole.CATALYST),
                 (lids, InputRole.CONDITION)), (lysed,))
    washes = graph.input("cold_high_salt_washes", "4 °C high-salt and 5× RT buffer washes",
                         MaterialKind.CONDITION)
    lysis_waste = graph.product("lysis_reagents", "removed oil, PFO, LiDS and Proteinase K",
                                MaterialKind.WASTE)
    washed = graph.product("washed_bead", "washed genome-capture bead with hybridized mRNA",
                           MaterialKind.MOLECULAR, state=bead.state)
    graph.react("break_and_wash", "Gel at 4 °C, break the emulsion and wash the beads cold",
                StepKind.CLEANUP,
                ((lysed, InputRole.TEMPLATE), (washes, InputRole.CONDITION)),
                (washed, lysis_waste))
    graph.discard(lysis_waste, "aspirated during the PFO drop-breaking and buffer washes")
    rtase = graph.input("reverse_transcriptase", "Maxima H Minus reverse transcriptase",
                        MaterialKind.CATALYST)
    tso = graph.input("tso", "template-switching oligo (Smart-seq3 biotin-TSO with N8 UMI "
                      "in the demonstrated run)", MaterialKind.OLIGO)
    copied = graph.product(
        "genome_capture_bead", "genome-retaining bead with matrix-bound cDNA",
        MaterialKind.MOLECULAR,
        state=single("genome-capture bead", (
            token("retained genome"), token("matrix-bound cDNA"))))
    graph.react("reverse_transcription",
                "Reverse-transcribe and template-switch on the hydrogel (30 min RT, 90 min 42 °C)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((washed, InputRole.TEMPLATE),
                 (capture_oligo, InputRole.STRUCTURAL),
                 (tso, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (copied,))
    taqman = graph.input("taqman_assay",
                         "multiplex TaqMan assays (demonstrated: HIV IPDA Ψ-FAM and env-VIC) "
                         "with TaqPath ProAmp master mix", MaterialKind.SELECTION)
    selected = graph.product(
        "selected_bead", "fluorescence-selected genome-capture bead",
        MaterialKind.MOLECULAR,
        state=single("fluorescent selected bead", (
            token("retained genome"), token("matrix-bound cDNA"),
            token("TaqMan marker signal"))))
    rejected = graph.product("rejected_beads", "droplets outside the double-positive gate",
                             MaterialKind.WASTE)
    graph.react("marker_sort",
                "Re-inject beads with TaqMan PCR, thermocycle (55 cycles, 88 °C denaturation) "
                "and sort double-positive droplets",
                StepKind.AFFINITY_CAPTURE,
                ((copied, InputRole.TEMPLATE), (taqman, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "droplets outside the gate are not deflected by the droplet sorter")
    recovery = graph.input("cdna_recovery_condition",
                           "freeze at −80 °C, then 10 min at 60 °C to melt agarose",
                           MaterialKind.CONDITION)
    recovered = graph.product("recovered_cdna", "aqueous phase with agarose-attached selected-cell cDNA",
                              MaterialKind.MOLECULAR,
                              state=single("selected-cell cDNA", (token("matrix-attached cDNA"),)))
    graph.react("recover_cdna", "Melt sorted droplets and transfer the aqueous phase",
                StepKind.EXTRACTION,
                ((selected, InputRole.TEMPLATE),
                 (recovery, InputRole.CONDITION)), (recovered,))
    preamp_forward = graph.input("smart_forward_primer", "Smart-seq3 forward (WTA) primer",
                                 MaterialKind.OLIGO)
    preamp_reverse = graph.input("smart_reverse_primer", "Smart-seq3 reverse (WTA) primer",
                                 MaterialKind.OLIGO)
    polymerase = graph.input("library_polymerase", "KAPA HiFi HotStart polymerase",
                             MaterialKind.CATALYST)
    amplified = graph.product("amplified_cdna", "free double-stranded WTA product",
                              MaterialKind.MOLECULAR, state=recovered.state)
    graph.react("smartseq_pcr", "Whole-transcriptome amplification with Smart-seq3 primers",
                StepKind.PCR,
                ((recovered, InputRole.TEMPLATE),
                 (preamp_forward, InputRole.STRUCTURAL),
                 (preamp_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Primer sequences are published, but the cDNA-end states are "
                    "not drawn; Smart-seq2 oligos are listed as an alternative."))
    ampure = graph.input("ampure_08x", "0.8× AMPure XP", MaterialKind.CONDITION)
    wta_clean = graph.product("purified_wta", "purified WTA product", MaterialKind.MOLECULAR,
                              state=amplified.state)
    wta_waste = graph.product("wta_cleanup_waste",
                              "unbound primers, agarose and agarose-attached template",
                              MaterialKind.WASTE)
    graph.react("wta_cleanup", "Purify WTA product with 0.8× AMPure XP", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (ampure, InputRole.CONDITION)),
                (wta_clean, wta_waste))
    graph.discard(wta_waste, "not retained by the 0.8× bead cleanup")
    loaded = graph.input("loaded_nextera_adapters", "Nextera XT kit tagmentation reagent",
                         MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cdna", "tagmented selected-cell cDNA",
                           MaterialKind.MOLECULAR, state=wta_clean.state)
    tn5_waste = graph.product("non_library_tn5_products", "non-library Tn5 products",
                              MaterialKind.WASTE)
    graph.react("tagment", "Tagment selected-cell cDNA", StepKind.TAGMENTATION,
                ((wta_clean, InputRole.TEMPLATE), (loaded, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "not the productive asymmetric Nextera product")
    gap_fill = graph.input("gap_fill_polymerase", "Nextera gap-fill polymerase",
                           MaterialKind.CATALYST)
    filled = graph.product("gap_filled_cdna", "gap-filled tagmented selected-cell cDNA",
                           MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("fill_tagmentation_gaps", "Fill the staggered tagmentation gaps",
                StepKind.EXTENSION,
                ((tagged, InputRole.TEMPLATE),
                 (gap_fill, InputRole.CATALYST)), (filled,))
    index_forward = graph.input("index_forward_primer", "Nextera indexed forward primer",
                                MaterialKind.OLIGO)
    index_reverse = graph.input("index_reverse_primer", "Nextera indexed reverse primer",
                                MaterialKind.OLIGO)
    library = graph.product("sequencing_library", "FIND-seq transcriptome library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "Index and amplify the FIND-seq library", StepKind.PCR,
                ((filled, InputRole.TEMPLATE),
                 (index_forward, InputRole.STRUCTURAL),
                 (index_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Index choice varies; the final Nextera endpoint is validated."))
    graph.sequence(library, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers)))
    graph.annotate(SourceDiscrepancy(
        "acrydite_oligodt", "Name of the agarose-conjugated oligo",
        (SourceClaim("Materials, Oligonucleotides",
                     "Acrydited oligo dT primer (Smart-seq3) is listed separately from the "
                     "Smart-seq3 forward primer (Nextera-handle WTA primer without acrydite or dT)."),
         SourceClaim("Anticipated results, first demonstration paragraph",
                     "agarose conjugated with the oligo dT primer (Smart-seq3 forward primer)")),
        "The bead capture oligo is the acrydited Smart-seq3 oligo(dT) primer.",
        "Only an acrydite-modified oligo(dT) can be linked to allyl agarose by APS/TEMED and "
        "capture poly(A) RNA; the listed forward primer has neither feature.",
        "The parenthetical label is treated as a naming slip; the protocol does not state "
        "otherwise.", AnnotationSlot.BEFORE_WORKFLOW))
    graph.declare_reviewed_scope((
        "cell_suspension", "allyl_agarose", "acrydite_oligodt", "encapsulate",
        "proteinase_k", "lids", "lyse", "cold_high_salt_washes", "break_and_wash",
        "reverse_transcriptase", "tso",
        "reverse_transcription", "taqman_assay", "marker_sort",
        "cdna_recovery_condition", "recover_cdna", "smart_forward_primer",
        "smart_reverse_primer", "library_polymerase", "smartseq_pcr",
        "ampure_08x", "wta_cleanup", "purified_wta",
        "loaded_nextera_adapters", "tn5", "tagment", "gap_fill_polymerase",
        "fill_tagmentation_gaps", "index_forward_primer",
        "index_reverse_primer", "index_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = 'find-seq'
    return find_seq_graph(d, final, primers)
SPEC.protocol_graph = protocol_graph
