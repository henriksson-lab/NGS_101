"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import seg, truseq_library
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, da_tailed, molecule, role
from dsb_capture import BreakCaptureDesign, _illumina_break_library, _finish_pcr


def workflow():
    return chain(molecule("purified genomic DNA", (role("intact genomic DNA"),)), (
        ("Cleave naked genomic DNA with Cas9–sgRNA in vitro",
         molecule("Cas9-digested genome", (role("uncut genomic DNA"), role("programmable Cas9 cleavage ends"))),
         "Cleavage is performed on cell-free DNA. All resulting DNA—not an enriched subset—continues into whole-genome library preparation."),
        ("Remove sgRNA and protein; purify the digested DNA",
         molecule("clean digenome", (role("genomic DNA retaining Cas9-created ends"),)),
         "RNase A and protein removal stop the reaction without adding a molecular tag to the cleavage ends."),
        ("Covaris-shear the whole digenome",
         molecule("mixed-origin WGS fragments", (role("fragments carrying one Cas9 end and one shear end"), role("fragments carrying two random shear ends"))),
         "Cas9 ends survive among a much larger background of ordinary random WGS fragment ends."),
        ("End-repair, size-select and add one 3-prime dA to each end",
         da_tailed("dA-tailed Digenome-seq WGS fragment", "unenriched genomic fragment"),
         "Both the Cas9-derived and randomly sheared endpoint classes pass through this end-prepared state."),
        ("INFERRED — ligate indexed TruSeq adapters to the dA-tailed fragments",
         molecule("Digenome-seq library", (role("TruSeq P5 arm"), role("unenriched genomic fragment"), role("TruSeq indexed P7 arm"))),
         "PCR-free preparation preserves read-start counts; there is no cut-end affinity capture or locus-specific amplification."),
        ("Paired-end whole-genome sequencing and cleavage-score analysis",
         molecule("aligned digenome reads", (role("forward and reverse 5-prime read-start pileups"), role("genome-wide depth background"))),
         "A true nuclease cut produces coincident read starts on both strands; ordinary shear endpoints remain dispersed."),
    ))


def _digenome_library() -> tuple[Construct, tuple]:
    return truseq_library(
        [seg("unenriched genomic DNA; Cas9 or random-shear boundary at either end",
             "X" * 42, placeholder=True)],
        "Digenome-seq TruSeq WGS library", dual_index=False, i7_nt=6,
        inferred_adapters=True)


TITLE = 'Digenome-seq — digested-genome whole-genome sequencing'
NOTES = '01_digenome_seq.html'
SOURCE = '<a href="https://doi.org/10.1038/nmeth.3284">Kim et al., Nature Methods (2015)</a>; library architecture follows the specified Illumina TruSeq DNA PCR-free WGS chemistry.'
SUMMARY = 'Cas9–sgRNA cleaves purified genomic DNA in vitro, after which the entire unenriched genome receives an ordinary PCR-free whole-genome library preparation. Coincident forward- and reverse-strand read starts reveal cleavage sites.'
CAVEAT = "Digenome-seq has no cut-site adapter, UMI or enrichment barcode; the assay-specific information is the genomic fragment boundary itself. INFERRED — no source available here gives one adapter base for this library: the defining paper's Methods is paywalled and absent, and its supplement contains no library chemistry. The TruSeq shell, the six-base single index and the index read are inherited from the authors' later Genome Research report and Illumina documentation, so the whole adapter shell is drawn as inferred."
FINAL_CAPTION = 'One ordinary TruSeq PCR-free WGS molecule. Its genomic insert may carry a Cas9-created boundary or only random shearing boundaries; the molecule itself has no label distinguishing them.'
READ_LENGTHS = {'Read 1': 150, 'Index 1 (i7)': 6, 'Read 2': 150}
SEQUENCING_INTRO = 'The defining study used Illumina HiSeq X Ten whole-genome sequencing. Digenome analysis depends on paired genomic alignments and the concentration of identical 5′ endpoints, not an inline molecular identifier.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
FINAL_LIBRARY, SEQ_PRIMERS = _digenome_library()
def oligos():
    return ()

def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    dna = graph.input("dna", "purified genomic DNA", MaterialKind.BIOLOGICAL,
                      state=molecule("purified genomic DNA", (role("intact genomic DNA"),)),
                      primary=True)
    sgrna = graph.input("sgrna", "programmable sgRNA", MaterialKind.OLIGO)
    cas9 = graph.input("cas9", "Cas9 nuclease", MaterialKind.CATALYST)
    cut = graph.product("cut", "Cas9-digested whole genome", MaterialKind.MOLECULAR,
                        state=dna.state)
    graph.react("cleave", "Cleave naked genomic DNA with Cas9–sgRNA",
                StepKind.GENOME_EDITING,
                ((dna, InputRole.TEMPLATE), (sgrna, InputRole.STRUCTURAL),
                 (cas9, InputRole.CATALYST)), (cut,))
    cleanup = graph.input("cleanup_condition", "RNase and protein-removal condition",
                          MaterialKind.CONDITION)
    clean = graph.product("clean", "purified digenome retaining Cas9-created ends",
                          MaterialKind.MOLECULAR, state=cut.state)
    graph.react("cleanup", "Remove sgRNA and protein", StepKind.DEPROTEINIZATION,
                ((cut, InputRole.TEMPLATE), (cleanup, InputRole.CONDITION)), (clean,))
    shear = graph.input("shear_condition", "Covaris shearing condition",
                        MaterialKind.CONDITION)
    fragments = graph.product("fragments", "mixed Cas9/shear and shear/shear WGS fragments",
                              MaterialKind.MOLECULAR, state=clean.state)
    graph.react("shear", "Covaris-shear the whole digenome", StepKind.FRAGMENTATION,
                ((clean, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    repair = graph.input("repair_enzyme", "WGS end-repair enzyme mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "blunt digenome WGS fragments", MaterialKind.MOLECULAR,
                             state=fragments.state)
    graph.react("repair", "End-repair all genomic fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    beads = graph.input("size_selection_beads", "WGS size-selection beads",
                        MaterialKind.SELECTION)
    selected = graph.product("selected", "size-selected WGS fragments",
                             MaterialKind.MOLECULAR, state=repaired.state)
    rejected = graph.product("rejected", "off-size genomic fragments", MaterialKind.WASTE)
    graph.react("size_select", "Size-select repaired fragments", StepKind.SELECTION,
                ((repaired, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "falls outside the whole-genome library insert range")
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "dA-tailed digenome WGS fragments",
                           MaterialKind.MOLECULAR, state=selected.state)
    graph.react("tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((selected, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("truseq_adapters", "TruSeq 3-prime-T adapters",
                           MaterialKind.ADAPTER)
    ligase = graph.input("ligase", "TruSeq adapter DNA ligase", MaterialKind.CATALYST)
    library = graph.product("library", "PCR-free Digenome-seq WGS library",
                            MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("ligate", "Ligate TruSeq adapters without PCR", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (library,),
                geometry=reported_adapter_end(
                    "TruSeq WGS adapter", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The PCR-free TruSeq adapter T overhang ligates to the dA-tailed WGS fragment."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=tuple(SEQ_PRIMERS)))
    reviewed = ["dna", "sgrna", "cas9", "cleave", "cleanup_condition", "cleanup",
                "shear_condition", "shear", "repair_enzyme", "repair",
                "size_selection_beads", "size_select", "tailing_enzyme", "tail",
                "truseq_adapters", "ligase", "ligate", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
