"""Local source and reaction model for truseq-dna-pcr-free."""
from chemdraw import Construct, Segment, Workflow
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from endprep import dA_tailed_scene, repair_and_dA_tail
from batch_ngs import truseq_full_length_adapter, truseq_pcr_free_rows
from illumina_vendor_workflows import Definition, _duplex, _sequence, _truseq_pcr_free_final, _workflow
from types import SimpleNamespace

TRUSEQ_FULL_LENGTH = truseq_full_length_adapter()
PCR_FREE_LIGATION_STATEMENT = (
    "The reference guide's ligation control joins “3′-T overhang adapters to 3′-A "
    "overhang inserts”, and Illumina Adapter Sequences prints both full-length arms, so "
    "the ligated adaptor already carries P5, i5, i7 and P7.")


def _pcr_free_workflow(final: Construct) -> Workflow:
    sheared = _duplex("Covaris-sheared genomic insert")
    repaired = _duplex("blunt size-selected genomic insert")
    end_prep = repair_and_dA_tail(Construct([
        Segment("size-selected genomic insert", "N" * 36, placeholder=True)]))
    tailed = dA_tailed_scene(end_prep, "dA-tailed genomic insert")
    w = _workflow("high-quality genomic DNA", _duplex("high-molecular-weight genomic DNA"), (
        ("Covaris shear", "fragmented DNA", sheared,
         "The guide specifies 350-bp or 550-bp insert workflows."),
        ("End-repair and two-sided SPB size-selection", "blunt size-selected DNA", repaired,
         "End Repair Mix 2 fills 5-prime overhangs and chews back 3-prime overhangs; diluted "
         "SPB then removes oversized fragments and undiluted SPB removes small ones."),
        ("Add one 3-prime dA to each end", "dA-tailed DNA", tailed,
         "The complementary adapter overhang is a single 3-prime T."),
    ))
    w.ligate_adapter(
        "Ligate full-length indexed TruSeq adapters",
        TRUSEQ_FULL_LENGTH.geometry(end_prep, statement=PCR_FREE_LIGATION_STATEMENT),
        truseq_pcr_free_rows(final, "genomic DNA insert",
                             label="PCR-free sequencing library"),
        name="PCR-free sequencing library",
        note="Sealing the adaptor finishes the library: Stop Ligation Buffer ends the 30 \u00b0C "
             "10-minute reaction, two SPB rounds clean the product, no enrichment PCR follows, "
             "and qPCR measures the cluster-forming molecules.")
    return w

def _pcr_free_graph(d, library, primers):
    source = d.source
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "high-quality high-molecular-weight genomic DNA",
                      MaterialKind.BIOLOGICAL,
                      state=_duplex("high-molecular-weight genomic DNA"), primary=True)
    shear = graph.input("shear_condition", "Covaris shearing condition",
                        MaterialKind.CONDITION)
    fragments = graph.product("fragments", "Covaris-sheared genomic DNA",
                              MaterialKind.MOLECULAR,
                              state=_duplex("Covaris-sheared genomic insert"))
    graph.react("shear", "Covaris-shear genomic DNA", StepKind.FRAGMENTATION,
                ((dna, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    repair = graph.input("repair_enzyme", "DNA end-repair enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "blunt repaired genomic fragments",
                             MaterialKind.MOLECULAR,
                             state=_duplex("blunt size-selected genomic insert"))
    graph.react("end_repair", "End-repair sheared DNA", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    # The guide size-selects in two sided rounds (Repair Ends and Select Library Size):
    # diluted SPB first binds the large fragments and the supernatant is kept, then
    # undiluted SPB binds the target inserts and that supernatant is discarded.
    diluted_beads = graph.input("diluted_spb", "diluted SPB (Sample Purification Beads)",
                                MaterialKind.SELECTION)
    retained_small = graph.product("below_upper_cutoff",
                                   "supernatant fragments below the upper size cutoff",
                                   MaterialKind.MOLECULAR, state=repaired.state)
    too_large = graph.product("large_fragments", "bead-bound oversized DNA fragments",
                              MaterialKind.WASTE)
    graph.react("remove_large", "Bind oversized fragments to diluted SPB and keep the supernatant",
                StepKind.SELECTION,
                ((repaired, InputRole.TEMPLATE), (diluted_beads, InputRole.SELECTION)),
                (retained_small, too_large))
    graph.discard(too_large, "stays on the diluted SPB beads when the supernatant is transferred")
    beads = graph.input("size_selection_beads", "undiluted SPB (Sample Purification Beads)",
                        MaterialKind.SELECTION)
    selected = graph.product("selected", "size-selected blunt genomic inserts",
                             MaterialKind.MOLECULAR, state=repaired.state)
    rejected = graph.product("rejected_sizes", "undersized DNA fragments", MaterialKind.WASTE)
    graph.react("size_select", "Bind the insert-size range to undiluted SPB and discard the supernatant",
                StepKind.SELECTION,
                ((retained_small, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "remains in the discarded supernatant of the undiluted-SPB round")
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme",
                          MaterialKind.CATALYST)
    end_prep = repair_and_dA_tail(Construct([
        Segment("size-selected genomic insert", "N" * 36, placeholder=True)]))
    tailed = graph.product("tailed", "single-dA-tailed genomic inserts",
                           MaterialKind.MOLECULAR,
                           state=dA_tailed_scene(end_prep, "dA-tailed genomic insert"))
    graph.react("tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((selected, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("adapter", "full-length indexed TruSeq 3-prime-T adapters",
                          MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "TruSeq adapter DNA ligase", MaterialKind.CATALYST)
    stop = graph.input("stop_ligation", "STL (Stop Ligation Buffer)", MaterialKind.CONDITION)
    ligated = graph.product("ligated", "adapter-ligated genomic fragments",
                            MaterialKind.MOLECULAR, state=library)
    geometry = TRUSEQ_FULL_LENGTH.geometry(
        end_prep, statement=PCR_FREE_LIGATION_STATEMENT)
    graph.react("ligate", "Ligate full-length indexed TruSeq adapters without enrichment PCR",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST), (stop, InputRole.CONDITION)),
                (ligated,), geometry=geometry)
    # "Clean Up Ligated Fragments": two SPB rounds before quantification; no PCR follows.
    cleanup_beads = graph.input("cleanup_spb", "SPB for the two ligation-cleanup rounds",
                                MaterialKind.SELECTION)
    final_product = graph.product("library", "PCR-free indexed TruSeq library",
                                  MaterialKind.LIBRARY, state=library)
    ligation_waste = graph.product("ligation_waste",
                                   "unbound adapters and ligation reagents",
                                   MaterialKind.WASTE)
    graph.react("clean_up_ligation", "Clean up the ligated fragments in two SPB rounds",
                StepKind.SELECTION,
                ((ligated, InputRole.TEMPLATE), (cleanup_beads, InputRole.SELECTION)),
                (final_product, ligation_waste))
    graph.discard(ligation_waste, "removed in the SPB supernatant and ethanol washes")
    _sequence(graph, final_product, library, primers, d.sequencing_intro)
    reviewed = ["dna", "shear_condition", "shear", "repair_enzyme", "end_repair",
                "diluted_spb", "remove_large", "size_selection_beads", "size_select",
                "tailing_enzyme", "tail", "adapter", "ligase", "stop_ligation", "ligate",
                "cleanup_spb", "clean_up_ligation", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

DEFINITION = Definition(
        "TruSeq DNA PCR-Free", "01_truseq_dna_pcr_free.html",
        'Illumina <a href="https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseq-dna-pcr-free-workflow/truseq-dna-pcr-free-workflow-reference-1000000039279-00.pdf">TruSeq DNA PCR-Free Reference Guide 1000000039279 v00</a>.',
        "Covaris-sheared genomic DNA is end-repaired, bead-size-selected, dA-tailed and ligated directly to indexed TruSeq adapters without library enrichment PCR.",
        "The guide lists separate single, combinatorial-dual and unique-dual index components. The schematic uses the combinatorial-dual adapter family; i5/i7 identities are placeholders for a chosen kit well.",
        "A cluster-capable dual-index TruSeq molecule made entirely by adapter ligation; there is no library PCR.",
        "Standard TruSeq primers read both genomic-insert ends and the two sample indexes.",
        "truseq_pcr_free", "genomic DNA insert", _pcr_free_workflow)


def protocol_graph():
    return _pcr_free_graph(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS)

FINAL_LIBRARY, SEQ_PRIMERS = _truseq_pcr_free_final(DEFINITION.insert_name, DEFINITION.title + " library", adapter=TRUSEQ_FULL_LENGTH)
d = DEFINITION
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes, SOURCE=d.source, SUMMARY=d.summary,
    CAVEAT=d.caveat, FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=SEQ_PRIMERS,
    FINAL_CAPTION=d.final_caption, SEQUENCING_INTRO=d.sequencing_intro,
    workflow=lambda: _pcr_free_workflow(FINAL_LIBRARY),
    final_rows=lambda: truseq_pcr_free_rows(FINAL_LIBRARY, d.insert_name, highlight_components=False),
    protocol_graph=protocol_graph,
)
globals().update(vars(SPEC))
