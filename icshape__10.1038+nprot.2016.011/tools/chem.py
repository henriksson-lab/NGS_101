from chemdraw import Workflow, Segment
from protocol_model import role
import illumina as il
from types import SimpleNamespace
from rna_cell_assays import CircularizationGeometry, Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, reported_adapter_end, single, token
from rna_cell_assays import custom_library, feature

from rna_cell_assays import chain

def icshape_workflow() -> Workflow:
    modified = token("NAI-N3 2-prime-O adduct", tag="umi", length=6)
    return chain("RNA in living cells", single("folded cellular RNA", (token("RNA"),)), (
        ("Acylate flexible RNA 2′-hydroxyls with NAI-N3",
         single("chemically probed RNA", (token("RNA upstream"), modified,
                                           token("RNA downstream"))),
         "NAI-N3 supplies both the structure-dependent adduct and an azide enrichment handle."),
        ("Click DIBO-biotin to the azide, fragment, and repair RNA 3′ ends",
         single("biotinylated RNA fragment", (token("RNA fragment"), modified,
                                               token("DIBO-biotin", tag="cbc", length=6))),
         "Copper-free click chemistry makes only NAI-N3-modified fragments selectable."),
        ("Ligate the preadenylated 3′ RNA adapter and enrich on streptavidin",
         single("adapter-ligated selected RNA", (token("RNA fragment"),
                                                  token("3-prime RNA adapter")),
                junctions=(("RNA fragment", "3-prime RNA adapter", "RNA ligation"),)),
         "T4 RNA ligase 1 joins the repaired RNA end to the blocked adapter."),
        ("Reverse-transcribe until the adduct stops synthesis",
         single("barcoded truncated cDNA", (token("13-nt experimental barcode + UMI", tag="cbc", length=13),
                                             token("cDNA ending one base before the adduct"))),
         "The RT primer carries 9 random bases and a 4-base experimental barcode."),
        ("Circularize cDNA with CircLigase II and PCR",
         single("icSHAPE library", (token("left Solexa arm"),
                                     token("13-nt experimental barcode + UMI", tag="cbc", length=13),
                                     token("RT-stop cDNA"), token("right Solexa arm"))),
         "A read's first transcript base marks the RT stop; analysis shifts one nucleotide to the adduct."),
    ))

DEFINITION = Definition("icSHAPE — in vivo click SHAPE", "10.1038/nprot.2016.011",
        "Flynn et al., <i>Nature Protocols</i> (2016)", "01_icshape.html",
        "NAI-N3 marks flexible RNA 2′-hydroxyls in cells, DIBO-biotin enables enrichment, and reverse-transcriptase stops encode nucleotide-resolution RNA structure.",
        "A 75-cycle single-end Illumina read begins with the 13-nt inline experimental barcode/UMI and then reports cDNA whose first mapped base locates the RT stop.")

d = DEFINITION
final = custom_library('icSHAPE library', (token('experimental barcode', tag='cbc', length=4, molecular_feature=feature('experimental_barcode', 'sample_index', 'fixed')), token('UMI', tag='umi', length=9, molecular_feature=feature('umi', 'umi', 'random')), token('RT-stop cDNA', length=40)))
primers = ()
reads = None
OLIGOS = (
            "<b>3′ adapter</b> — 5′-/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/-3′",
            "<b>RT-primer-1 architecture</b> — 5′-/5Phos/DDDNN AACC NNNN AGATCGGAAGAGCGTCGTGGA /iSp18/ GGATCC /iSp18/ TACTGAACCGC-3′",
            "<b>P5-Solexa PCR primer</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′",
            "<b>P3-Solexa PCR primer</b> — 5′-CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT-3′",
        )
workflow_fn = icshape_workflow
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
    SEQ_PRIMERS=primers, FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=reads, SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("" if primers else "The defining protocol prints both Solexa library-PCR primers, but the physical run-primer binding geometry has not been reconstructed from the final circularized and amplified library; read entry and cycle roles are stated explicitly."),
    workflow=workflow_fn, oligos=lambda: OLIGOS,
)


def icshape_graph(definition, final: Construct) -> ProtocolGraph:
    """Typed icSHAPE graph through click enrichment, RT stop and circularization."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    rna = graph.input("cellular_rna", "folded RNA in living cells",
                      MaterialKind.BIOLOGICAL,
                      state=single("folded cellular RNA", (token("RNA"),)),
                      primary=True)
    nai = graph.input("nai_n3", "NAI-N3 acylation reagent", MaterialKind.CONDITION)
    modified_state = single("chemically probed RNA", (
        role("RNA bearing 2-prime-O NAI-N3 substituents; not nucleotide identifiers"),))
    modified = graph.product("modified_rna", "NAI-N3-modified RNA",
                             MaterialKind.MOLECULAR, state=modified_state)
    graph.react("acylate", "Acylate flexible RNA 2-prime-hydroxyls",
                StepKind.CHEMICAL_CONVERSION,
                ((rna, InputRole.TEMPLATE), (nai, InputRole.CONDITION)), (modified,))

    dibo = graph.input("dibo_biotin", "DIBO-biotin click reagent",
                       MaterialKind.MOLECULAR)
    click_condition = graph.input("click_condition", "copper-free click condition",
                                  MaterialKind.CONDITION)
    biotinylated = graph.product(
        "biotinylated_rna", "DIBO-biotin-tagged modified RNA",
        MaterialKind.MOLECULAR,
        state=single("biotinylated RNA", (
            role("RNA bearing biotin-clicked NAI-N3 substituents; not a barcode"),)))
    graph.react("click_biotin", "Click DIBO-biotin to the NAI-N3 azide",
                StepKind.CHEMICAL_CONVERSION,
                ((modified, InputRole.TEMPLATE), (dibo, InputRole.STRUCTURAL),
                 (click_condition, InputRole.CONDITION)), (biotinylated,))
    fragmentation = graph.input("rna_fragmentation_reagent", "RNA fragmentation reagent",
                                MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_rna", "fragmented biotinylated RNA",
                               MaterialKind.MOLECULAR, state=biotinylated.state)
    graph.react("fragment_rna", "Fragment the biotinylated RNA",
                StepKind.FRAGMENTATION,
                ((biotinylated, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragmented,))
    repair = graph.input("rna_end_repair_enzyme", "RNA 3-prime-end repair enzyme",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired_rna", "3-prime-end-repaired RNA fragments",
                             MaterialKind.MOLECULAR, state=fragmented.state)
    graph.react("repair_rna_end", "Repair RNA 3-prime ends", StepKind.END_REPAIR,
                ((fragmented, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))

    adapter = graph.input("preadenylated_adapter",
                          "preadenylated blocked 3-prime RNA adapter",
                          MaterialKind.ADAPTER, state=Construct([Segment("published RNA adapter", "AGATCGGAAGAGCGGTTCAG")]))
    ligase = graph.input("t4_rnl1", "T4 RNA ligase 1",
                         MaterialKind.CATALYST)
    ligated = graph.product(
        "adapter_ligated_rna", "3-prime-adapter-ligated RNA",
        MaterialKind.MOLECULAR,
        state=single("adapter-ligated RNA", (
            token("RNA fragment"), token("3-prime RNA adapter")),
            junctions=(("RNA fragment", "3-prime RNA adapter", "RNA ligation"),)))
    adapter_geometry = reported_adapter_end(
        "icSHAPE 3-prime RNA adapter", kind="single-strand", polarity="5-prime",
        terminal_sequence="AGATCGGAAGAGCGGTTCAG",
        statement="The protocol prints the 5′ preadenylated, 3′-blocked RNA adapter sequence used by truncated T4 RNL2.")
    graph.react("adapter_ligation", "Ligate the preadenylated 3-prime RNA adapter",
                StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=adapter_geometry)
    rt_primer = graph.input("rt_primer", "barcoded RT-primer-1",
                            MaterialKind.OLIGO)
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    truncated = graph.product(
        "truncated_cdna", "RT-stop cDNA with its biotinylated RNA template",
        MaterialKind.MOLECULAR,
        state=single("barcoded truncated cDNA", (
            token("13-nt experimental barcode + UMI", tag="cbc", length=13),
            token("cDNA ending one base before the adduct"))))
    graph.react("reverse_transcription", "Reverse-transcribe until the adduct stops synthesis",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (truncated,))
    streptavidin = graph.input("streptavidin", "streptavidin capture reagent",
                               MaterialKind.SELECTION)
    selected = graph.product("selected_rna", "streptavidin-selected RNA/cDNA hybrids",
                             MaterialKind.MOLECULAR, state=truncated.state)
    unselected = graph.product("unselected_rna", "unmodified or uncaptured RNA",
                               MaterialKind.WASTE)
    graph.react("streptavidin_capture", "Enrich biotinylated RNA",
                StepKind.AFFINITY_CAPTURE,
                ((truncated, InputRole.TEMPLATE),
                 (streptavidin, InputRole.SELECTION)), (selected, unselected))
    graph.discard(unselected, "not retained by streptavidin enrichment")

    rnases = graph.input("rna_release", "RNase A/T1/H followed by DMSO/heat elution", MaterialKind.CATALYST)
    released = graph.product("released_cdna", "cDNA released from selected hybrids", MaterialKind.MOLECULAR, state=truncated.state)
    graph.react("release_cdna", "Digest RNA and elute cDNA from the beads", StepKind.DIGESTION,
                ((selected, InputRole.TEMPLATE), (rnases, InputRole.CATALYST)), (released,))
    circligase = graph.input("circligase_ii", "CircLigase II", MaterialKind.CATALYST)
    circular = graph.product("circular_cdna", "circularized RT-stop cDNA",
                             MaterialKind.MOLECULAR, state=truncated.state)
    graph.react("circularize", "Circularize cDNA", StepKind.CIRCULARIZATION,
                ((released, InputRole.TEMPLATE),
                 (circligase, InputRole.CATALYST)), (circular,),
                geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="The icSHAPE workflow circularizes RT-stop cDNA with CircLigase II."))
    forward = graph.input("pcr_forward_primer", "Solexa PCR forward primer",
                          MaterialKind.OLIGO, state=Construct([Segment("published P5-Solexa", il.P5 + il.TRUSEQ_READ1[4:])]))
    reverse = graph.input("pcr_reverse_primer", "Solexa PCR reverse primer",
                          MaterialKind.OLIGO, state=Construct([Segment("published P3-Solexa", il.P7 + "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT")]))
    polymerase = graph.input("pcr_polymerase", "PCR polymerase", MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "icSHAPE sequencing library",
                            MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("library_pcr", "PCR-amplify the circular cDNA", StepKind.PCR,
                ((circular, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    Construct([token("unassembled circular cDNA template")]), reason="Both ordered Solexa PCR primers are printed in the original reagent list. Exact binding geometry on the circular cDNA remains unassembled in this model; this is not a source-absence claim."))
    graph.sequence(library, SequencingPlan(limitation=definition.endpoint))
    graph.declare_reviewed_scope((
        "cellular_rna", "nai_n3", "acylate", "dibo_biotin", "click_condition",
        "click_biotin", "rna_fragmentation_reagent", "fragment_rna",
        "rna_end_repair_enzyme", "repair_rna_end", "preadenylated_adapter",
        "t4_rnl1", "adapter_ligation", "streptavidin",
        "streptavidin_capture", "rt_primer", "reverse_transcriptase",
        "reverse_transcription", "rna_release", "release_cdna", "circligase_ii", "circularize",
        "pcr_forward_primer", "pcr_reverse_primer", "pcr_polymerase",
        "library_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = 'icshape'
    return icshape_graph(d, final)
SPEC.protocol_graph = protocol_graph
SPEC.CAVEAT += (" This bounded update follows the NAI-N3-treated branch: RT precedes streptavidin capture, followed by RNase/DMSO-heat cDNA release. "
    "The graph still summarizes gel/column recovery stages and leaves full circular-template/final-library assembly unresolved. "
    "Published unit discrepancy: Step49 lists 1 µl of 1 µM P3 primer in 50 µl total, but its final-concentration column says 20 mM. "
    "Those listed volumes imply 20 nM by dilution; the intended value is unresolved and the model assigns no corrected concentration.")
