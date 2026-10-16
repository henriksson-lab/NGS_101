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
        ("End-repair and bead size-selection", "blunt size-selected DNA", repaired,
         "End repair removes the 3-prime and 5-prime overhangs made by shearing."),
        ("Add one 3-prime dA to each end", "dA-tailed DNA", tailed,
         "The complementary adapter overhang is a single 3-prime T."),
    ))
    w.ligate_adapter(
        "Ligate full-length indexed TruSeq adapters",
        TRUSEQ_FULL_LENGTH.geometry(end_prep, statement=PCR_FREE_LIGATION_STATEMENT),
        truseq_pcr_free_rows(final, "genomic DNA insert",
                             label="PCR-free sequencing library"),
        name="PCR-free sequencing library",
        note="Sealing the adaptor finishes the library: no enrichment PCR follows, and qPCR "
             "measures the cluster-forming molecules.")
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
    beads = graph.input("size_selection_beads", "size-selection magnetic beads",
                        MaterialKind.SELECTION)
    selected = graph.product("selected", "size-selected blunt genomic inserts",
                             MaterialKind.MOLECULAR, state=repaired.state)
    rejected = graph.product("rejected_sizes", "off-size DNA fragments", MaterialKind.WASTE)
    graph.react("size_select", "Select the requested insert-size range",
                StepKind.SELECTION,
                ((repaired, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "falls outside the selected library insert range")
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
    final_product = graph.product("library", "PCR-free indexed TruSeq library",
                                  MaterialKind.LIBRARY, state=library)
    geometry = TRUSEQ_FULL_LENGTH.geometry(
        end_prep, statement=PCR_FREE_LIGATION_STATEMENT)
    graph.react("ligate", "Ligate full-length indexed TruSeq adapters without enrichment PCR",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (final_product,), geometry=geometry)
    _sequence(graph, final_product, library, primers, d.sequencing_intro)
    reviewed = ["dna", "shear_condition", "shear", "repair_enzyme", "end_repair",
                "size_selection_beads", "size_select", "tailing_enzyme", "tail",
                "adapter", "ligase", "ligate", "library"]
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
