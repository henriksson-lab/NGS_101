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

SITE_A1_F = "GTTGACATGCTGGATTGAGACTTCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"

SITE_A1_R = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGAAGTCTCAATCCAGCATGTCAAC"

SITE_A2_R = "ACTACTAATACGACT"

SITE_A2_F7 = "GTCGTATTAGTAGTANNNNNNNAGATCGGAAGAGCACACGTCTGAACTCC"

SITE_A2_F6 = "GTCGTATTAGTAGTANNNNNNAGATCGGAAGAGCACACGTCTGAACTCC"

SITE_A2_F5 = "GTCGTATTAGTAGTANNNNNAGATCGGAAGAGCACACGTCTGAACTCC"

SITE_RECOVERY_F = "GGAGTTCAGACGTGTGCTC"

SITE_RECOVERY_R = "GTTGACATGCTGGATTGAGACTTC"

SITE_INDEX_F = "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACG"

SITE_INDEX_R = "CAAGCAGAAGACGGCATACGAGATTANNNNNNNNGACTGGAGTTCAGACGTGTGCTC"


def workflow():
    return chain(molecule("purified HMW genomic DNA", (role("intact genomic DNA"),)), (
        ("Cleave naked genomic DNA with Cas9–sgRNA in vitro",
         molecule("Cas9-digested genome", (role("genomic DNA"), role("programmable Cas9 cleavage ends"))),
         "A concentration series of each sgRNP exposes biochemically accessible on- and off-target sites; this is not a measurement of cleavage inside cells."),
        ("Add one 3-prime dA to each Cas9 cleavage end",
         da_tailed("dA-tailed Cas9-cleaved genome", "Cas9-cleavage-adjacent genomic DNA"),
         "The Cas9-created ends are shown in their cohesive ligation-ready state."),
        ("Ligate biotinylated Adapter 1 to the dA-tailed cleavage ends",
         molecule("cut-end-labelled genome", (role("5-prime-biotin Adapter 1"), role("Cas9-cleavage-adjacent genomic DNA")),
                  junctions=(("5-prime-biotin Adapter 1", "Cas9-cleavage-adjacent genomic DNA", "T/A ligation"),)),
         "Adapter 1 is installed before fragmentation, so its biotin identifies the Cas9-created end."),
        ("Fragment with dsDNA Fragmentase to about 200–800 bp",
         molecule("labelled genomic fragments", (role("biotinylated Cas9 end"), role("cut-to-random genomic fragment"), role("Fragmentase-created distal end"))),
         "The new distal end is physically distinct from the already labelled cleavage end."),
        ("End-repair and ligate the pooled N5/N6/N7 Adapter 2",
         molecule("dual-adapter SITE-Seq fragment", (role("biotinylated Adapter 1"), role("cut-to-random genomic DNA"), role("5/6/7-nt diversity spacer + Adapter 2")),
                  junctions=(("biotinylated Adapter 1", "cut-to-random genomic DNA", "cleavage-end ligation"),
                             ("cut-to-random genomic DNA", "5/6/7-nt diversity spacer + Adapter 2", "distal ligation"))),
         "Pooling three adapter lengths creates base diversity; the paper does not define this tract as a molecule-counting UMI."),
        ("Capture biotinylated fragments on streptavidin beads",
         molecule("enriched SITE-Seq fragments", (role("Adapter 1-labelled Cas9 end"), role("adjacent genomic DNA"), role("Adapter 2 distal end"))),
         "Only molecules retaining the Adapter 1-labelled cleavage end survive the affinity selection."),
        ("Recovery PCR, then dual-index PCR",
         molecule("SITE-Seq sequencing library", (role("P5"), role("i5 sample index"), role("Read 1 site"), role("DNA beginning at a Cas9 cut"), role("5/6/7-nt diversity spacer"), role("Index 1 site"), role("i7 sample index"), role("P7"))),
         "The two published PCR rounds first recover captured molecules and then install both flow-cell ends and eight-base sample indexes."),
    ))


def _site_seq_library() -> tuple[Construct, tuple]:
    """Indexed product, with the N7 Adapter 2 branch shown explicitly."""
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("i5 sample index", "I" * 8, "cbc", placeholder=True,
                feature=feature("sample_index_i5", "sample_index", "unknown")),
        Segment("Read 1 arm", il.TRUSEQ_READ1, "r1"),
        role("genomic DNA beginning at the in-vitro Cas9 cleavage end"),
        Segment("Adapter 2 length-diversity spacer (N7 branch)", "N" * 7, "cbc", placeholder=True),
        Segment("Index 1 + Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("i7 sample index reverse complement", "I" * 8, "cbc", placeholder=True,
                feature=feature("sample_index_i7", "sample_index", "unknown")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="SITE-Seq library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("SITE-Seq final library: " + "; ".join(problems))
    return lib, primers


TITLE = 'SITE-Seq — in-vitro Cas9 cleavage-site sequencing'
NOTES = '01_site_seq.html'
SOURCE = '<a href="https://doi.org/10.1038/nmeth.4284">Cameron et al., Nature Methods (2017)</a>; exact chemistry from the published supplementary protocol and Supplementary Table 26.'
SUMMARY = 'Cas9–sgRNA cleaves purified genomic DNA in vitro; a biotinylated adapter labels those ends before enzymatic fragmentation, allowing streptavidin enrichment and Read 1 entry directly beside each candidate cut site.'
CAVEAT = ''
FINAL_CAPTION = 'N7 Adapter 2 branch; N5 and N6 differ only in diversity-spacer length. Adapter 1 defines the Cas9-cleavage end, while Adapter 2 occupies the later Fragmentase end.'
READ_LENGTHS = {'Read 1': 150, 'Index 1 (i7)': 8, 'Index 2 (i5)': 8}
SEQUENCING_INTRO = 'The supplementary protocol recommends 2–3 million 150-base single-end reads per sample. Both published index primers carry eight variable bases; the standard index primers land on the completed dual-index library.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
FINAL_LIBRARY, SEQ_PRIMERS = _site_seq_library()
def oligos():
    return (f'Adapter 1 Forward — 5′-/5Biosg/{SITE_A1_F[:-1]}*{SITE_A1_F[-1]}-3′', f'Adapter 1 Reverse — 5′-{SITE_A1_R}-3′', f'Adapter 2 N7 Forward — 5′-/5Phos/{SITE_A2_F7}-3′', f'Adapter 2 N6 Forward — 5′-/5Phos/{SITE_A2_F6}-3′', f'Adapter 2 N5 Forward — 5′-/5Phos/{SITE_A2_F5}-3′', f'Adapter 2 Reverse — 5′-{SITE_A2_R[:-1]}*{SITE_A2_R[-1]}-3′', f'Recovery PCR Forward — 5′-{SITE_RECOVERY_F}-3′', f'Recovery PCR Reverse — 5′-{SITE_RECOVERY_R}-3′', f'Index Forward — 5′-{SITE_INDEX_F}-3′', f'Index Reverse — 5′-{SITE_INDEX_R}-3′')

def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    dna = graph.input("dna", "purified high-molecular-weight genomic DNA",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("purified HMW genomic DNA", (role("intact genomic DNA"),)),
                      primary=True)
    sgrna = graph.input("sgrna", "programmable sgRNA", MaterialKind.OLIGO)
    cas9 = graph.input("cas9", "Cas9 nuclease", MaterialKind.CATALYST)
    cut = graph.product("cut", "Cas9-cleaved genomic DNA", MaterialKind.MOLECULAR,
                        state=dna.state)
    graph.react("cleave", "Cleave naked genomic DNA with Cas9–sgRNA",
                StepKind.GENOME_EDITING,
                ((dna, InputRole.TEMPLATE), (sgrna, InputRole.STRUCTURAL),
                 (cas9, InputRole.CATALYST)), (cut,))
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "dA-tailed Cas9 cleavage ends",
                           MaterialKind.MOLECULAR, state=cut.state)
    graph.react("tail", "Add one 3-prime dA at Cas9 cleavage ends", StepKind.TAILING,
                ((cut, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter1 = graph.input("adapter1", "biotinylated T-overhang Adapter 1",
                           MaterialKind.ADAPTER)
    ligase1 = graph.input("ligase1", "Adapter-1 DNA ligase", MaterialKind.CATALYST)
    labelled = graph.product("labelled", "Adapter-1-labelled Cas9 ends",
                             MaterialKind.MOLECULAR, state=tailed.state)
    graph.react("ligate_adapter1", "Ligate biotinylated Adapter 1 at cleavage ends",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter1, InputRole.STRUCTURAL),
                 (ligase1, InputRole.CATALYST)), (labelled,),
                geometry=reported_adapter_end(
                    "SITE-Seq Adapter 1", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="Adapter 1 labels the dA-tailed Cas9 end before fragmentation."))
    fragmentase = graph.input("fragmentase", "dsDNA Fragmentase", MaterialKind.CATALYST)
    fragments = graph.product("fragments", "Adapter-1-labelled 200–800-bp fragments",
                              MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("fragment", "Fragment labelled genomic DNA", StepKind.FRAGMENTATION,
                ((labelled, InputRole.TEMPLATE), (fragmentase, InputRole.CATALYST)),
                (fragments,))
    repair = graph.input("distal_repair", "distal-end repair enzyme mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired distal fragments", MaterialKind.MOLECULAR,
                             state=fragments.state)
    graph.react("repair_distal", "End-repair Fragmentase-created distal ends",
                StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    adapter2 = graph.input("adapter2", "pooled N5/N6/N7 Adapter 2", MaterialKind.ADAPTER)
    ligase2 = graph.input("ligase2", "Adapter-2 DNA ligase", MaterialKind.CATALYST)
    dual = graph.product("dual", "dual-adapter SITE-Seq fragments", MaterialKind.MOLECULAR,
                         state=repaired.state)
    graph.react("ligate_adapter2", "Ligate the length-diverse distal Adapter 2",
                StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (adapter2, InputRole.STRUCTURAL),
                 (ligase2, InputRole.CATALYST)), (dual,),
                geometry=reported_adapter_end(
                    "SITE-Seq Adapter 2", kind="blunt",
                    statement="The Fragmentase-created distal end is repaired before the pooled N5/N6/N7 Adapter 2 is ligated."))
    streptavidin = graph.input("streptavidin", "streptavidin capture beads",
                               MaterialKind.SELECTION)
    captured = graph.product("captured", "biotin-enriched SITE-Seq fragments",
                             MaterialKind.MOLECULAR, state=dual.state)
    rejected = graph.product("rejected", "Adapter-1-negative fragments", MaterialKind.WASTE)
    graph.react("capture", "Capture Adapter-1-labelled fragments",
                StepKind.AFFINITY_CAPTURE,
                ((dual, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not retain the cleavage-end Adapter 1")
    rec_f = graph.input("recovery_forward", "recovery PCR forward primer", MaterialKind.OLIGO)
    rec_r = graph.input("recovery_reverse", "recovery PCR reverse primer", MaterialKind.OLIGO)
    polymerase = graph.input("polymerase", "SITE-Seq PCR polymerase", MaterialKind.CATALYST)
    recovered = graph.product("recovered", "recovery-amplified SITE-Seq fragments",
                              MaterialKind.MOLECULAR, state=captured.state)
    graph.react("recovery_pcr", "Recovery-PCR captured fragments", StepKind.PCR,
                ((captured, InputRole.TEMPLATE), (rec_f, InputRole.STRUCTURAL),
                 (rec_r, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (recovered,), geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="Recovery PCR precedes the separately indexed final construct."))
    reviewed = ["dna", "sgrna", "cas9", "cleave", "tailing_enzyme", "tail",
                "adapter1", "ligase1", "ligate_adapter1", "fragmentase", "fragment",
                "distal_repair", "repair_distal", "adapter2", "ligase2", "ligate_adapter2",
                "streptavidin", "capture", "recovery_forward", "recovery_reverse",
                "polymerase", "recovery_pcr"]
    _finish_pcr(graph, "index", recovered, FINAL_LIBRARY, SEQ_PRIMERS,
                SEQUENCING_INTRO, reviewed, label="SITE-Seq indexed")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
