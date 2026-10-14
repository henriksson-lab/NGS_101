"""Direct double-strand-break capture chemistries.

The shared model keeps the two molecular ends distinct: the proximal end is the
pre-existing cellular break labelled in situ, whereas the distal end is created later
by fragmentation.  A constructor cannot silently exchange those origins.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import seg, truseq_library
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, da_tailed, molecule, role


@dataclass(frozen=True)
class BreakCaptureDesign:
    name: str
    proximal_adapter: str
    distal_adapter: str
    immobilization: str
    direct_break_read: bool

    def __post_init__(self):
        if not all((self.name, self.proximal_adapter, self.distal_adapter,
                    self.immobilization)):
            raise ValueError("break capture requires named proximal/distal ends and immobilization")
        if self.proximal_adapter == self.distal_adapter:
            raise ValueError("the native-break and fragmentation ends must remain distinguishable")


# Crosetto et al. Supplementary Table 5; the hairpins are written continuously 5'->3'.
BLESS_PROXIMAL = "TACTACCTCGAGAGTTACGCTAGGGATAACAGGGTAATATAGTTTTTTTCTATATTACCCTGTTATCCCTAGCGTAACTCTCGAGGTAGTA"
BLESS_DISTAL = "CGTCGTCTCGAGAGTTACGCTAGGGATAACAGGGTAATATAGTTTTTTTCTATATTACCCTGTTATCCCTAGCGTAACTCTCGAGACGACG"
BLESS_PCR_PROXIMAL = "CCCTAGCGTAACTCTCGAGGTAGTA"
BLESS_PCR_DISTAL = "CTAGCGTAACTCTCGAGACGACG"
BLESS_PROXIMAL_READ_TAG = "TCGAGGTAGTA"
BLESS_DISTAL_READ_TAG = "TCGAGACGACG"

# Lensing et al. publish the DSBCapture oligo set in Fig. 1a and Supplementary Fig. 1c
# *as raster figure artwork only*: no base of it occurs as text anywhere in the author
# manuscript, Supplementary Figures or Supplementary Tables, so the transcription cannot
# be re-checked against extracted source text (``tools/selftest_source_oligos.py`` records
# that as a declared exception below).  Every drawn base is a canonical Illumina TruSeq
# part, so the set is assembled here from ``illumina`` and from its own strands rather than
# retyped — which also means a figure transcription and a reconstruction out of canonical
# constants are indistinguishable.  See 01_dsbcapture.md "Oligo provenance".
DSB_P5 = il.TRUSEQ_P5_FULL                 # modified P5 top strand as drawn, 58 nt
DSB_P7 = il.INDEX1_PRIMER                  # modified P7 top strand as drawn, 33 nt
# Both duplexes are drawn with one unpaired 3′ T on the strand that T/A-ligates to the
# dA-tailed genomic end, so each bottom strand follows from its own top strand.
DSB_P5_BOTTOM = revcomp(DSB_P5)[1:]        # the 3′ T of the P5 top strand stays unpaired
DSB_P7_BOTTOM = revcomp(DSB_P7) + "T"      # the P7 duplex carries its T on the bottom strand
DSB_PCR_F = DSB_P5[:44]                    # drawn as the first 44 bases of the P5 adapter
DSB_AD006 = "GCCAAT"                       # TruSeq sample index AD006, read as the i7 index
# PCR R is drawn as P7, the index reverse complement, then 21 bases complementary to the
# modified P7 top strand; the paper prints three such example primers (AD006 is shown).
DSB_PCR_R_AD06 = il.P7 + revcomp(DSB_AD006) + revcomp(DSB_P7)[:21]

# Cameron et al. SITE-Seq supplementary protocol / Supplementary Table 26.
# Modifications are kept in the ordering strings returned by ``get``; these constants
# are the molecular bases used to derive the final library.
SITE_A1_F = "GTTGACATGCTGGATTGAGACTTCCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
SITE_A1_R = "GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGAAGTCTCAATCCAGCATGTCAAC"
# The 15th base is the 3'-phosphorothioate T overhang that T/A-ligates to the
# end-repaired, dA-tailed Fragmentase end; the first 14 bases pair SITE_A2_F*[:14].
SITE_A2_R = "ACTACTAATACGACT"
SITE_A2_F7 = "GTCGTATTAGTAGTANNNNNNNAGATCGGAAGAGCACACGTCTGAACTCC"
SITE_A2_F6 = "GTCGTATTAGTAGTANNNNNNAGATCGGAAGAGCACACGTCTGAACTCC"
SITE_A2_F5 = "GTCGTATTAGTAGTANNNNNAGATCGGAAGAGCACACGTCTGAACTCC"
SITE_RECOVERY_F = "GGAGTTCAGACGTGTGCTC"
SITE_RECOVERY_R = "GTTGACATGCTGGATTGAGACTTC"
SITE_INDEX_F = "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACG"
SITE_INDEX_R = "CAAGCAGAAGACGGCATACGAGATTANNNNNNNNGACTGGAGTTCAGACGTGTGCTC"


DESIGNS = {
    "bless": BreakCaptureDesign("BLESS", "biotinylated proximal hairpin",
                                 "distal hairpin", "intact fixed nuclei", False),
    "i-bless": BreakCaptureDesign("i-BLESS", "biotinylated proximal hairpin",
                                   "distal hairpin", "cells embedded in agarose beads", False),
    "dsbcapture": BreakCaptureDesign("DSBCapture", "biotinylated modified P5 adapter",
                                      "modified P7 adapter", "intact fixed nuclei", True),
    "site-seq": BreakCaptureDesign("SITE-Seq", "biotinylated Adapter 1",
                                    "degenerate Adapter 2", "purified HMW genomic DNA", True),
}


def _bless_workflow(key: str):
    d = DESIGNS[key]
    initial = molecule("yeast cells" if key == "i-bless" else d.immobilization,
                       (role("genomic DNA with a native DSB"),))
    prefix = ((
        "Embed, spheroplast and deproteinate cells in agarose beads",
        molecule("agarose-protected genomic DNA", (role("native DSB in deproteinated DNA"),)),
        "Agarose-bead immobilization is the defining i-BLESS change and protects yeast DNA from handling-induced shearing.",
    ),) if key == "i-bless" else ()
    return chain(initial, prefix + (
        ("Blunt and 5′-phosphorylate native break ends",
         molecule("ligation-ready native breaks", (role("blunt 5-prime-phosphorylated DSB"),)),
         "End repair makes 5′ overhangs, 3′ overhangs and blunt breaks enter the same ligation state."),
        ("Ligate the biotinylated proximal hairpin in situ",
         molecule("proximal-labelled break", (role("biotinylated proximal hairpin"), role("break-adjacent genomic DNA")),
                  junctions=(("biotinylated proximal hairpin", "break-adjacent genomic DNA", "T4 DNA ligase"),)),
         "The proximal barcode marks the end that existed in the cell before extraction."),
        ("Fragment DNA and capture proximal-labelled molecules on streptavidin",
         molecule("captured break fragment", (role("biotinylated proximal hairpin"), role("break-to-random genomic fragment"), role("fragmentation-created free end"))),
         "Fragmentation creates the distal end only after the native break has already been labelled."),
        ("Blunt the free end and ligate the distal hairpin",
         molecule("dual-hairpin break fragment", (role("proximal hairpin"), role("break-to-random genomic fragment"), role("distal hairpin")),
                  junctions=(("proximal hairpin", "break-to-random genomic fragment", "native-break ligation"),
                             ("break-to-random genomic fragment", "distal hairpin", "distal ligation"))),
         "Different proximal and distal barcodes preserve which end is the native DSB."),
        ("Linearize with I-SceI, PCR-amplify, then cut with XhoI",
         molecule("barcoded BLESS insert", (role("proximal 11-nt tag"), role("break-to-random genomic DNA"), role("distal 11-nt tag"))),
         "I-SceI releases the captured molecule from its hairpins; XhoI leaves the short end-identifying tags."),
        ("INFERRED — install the sequencing-kit ends",
         molecule("Illumina-compatible break library", (role("Illumina left end", inferred=True), role("proximal or distal tag"), role("genomic insert"), role("Illumina right end", inferred=True))),
         "The paper names the final library kit but does not print its complete adapter construct."),
    ))


def _dsbcapture_workflow():
    return chain(molecule("fixed nuclei", (role("genomic DNA with a native DSB"),)), (
        ("Blunt native break ends and add one 3′ dA",
         molecule("A-tailed native break", (role("repaired DSB"), role("3-prime dA"))),
         "End repair normalizes the break; Klenow exo− then creates a cohesive ligation end."),
        ("Ligate the biotinylated, T-overhang modified P5 adapter in situ",
         molecule("P5-labelled break", (role("biotinylated modified P5 adapter"), role("break-adjacent genomic DNA")),
                  junctions=(("biotinylated modified P5 adapter", "break-adjacent genomic DNA", "T/A ligation"),)),
         "Because the Read 1 site is already at the native break, R1 starts with the first genomic base."),
        ("Digest free adapter with lambda exonuclease, extract and sonicate",
         molecule("labelled sonication fragment", (role("modified P5 at native break"), role("break-to-random genomic fragment"), role("sonication-created free end"))),
         "Covalently break-labelled fragments retain biotin; free adapters are removed before extraction."),
        ("Streptavidin-capture, blunt and add one 3-prime dA to the distal end",
         da_tailed("distal-end dA-tailed DSBCapture fragment", "modified-P5-labelled break-to-random genomic DNA"),
         "The native-break P5 end is already labelled; only the sonication-created distal end is prepared."),
        ("Ligate modified P7 to the dA-tailed distal end",
         molecule("dual-adapter DSBCapture fragment", (role("modified P5"), role("break-to-random genomic DNA"), role("modified P7")),
                  junctions=(("modified P5", "break-to-random genomic DNA", "native-break ligation"),
                             ("break-to-random genomic DNA", "modified P7", "distal T/A ligation"))),
         "The second adapter is added only to the random sonication end."),
        ("INFERRED — PCR-add P5, i7 and P7; size-select the library",
         molecule("DSBCapture sequencing library", (role("P5 and Read 1 site"), role("DNA beginning at the DSB"), role("Read 2 and Index 1 site"), role("six-base sample index", inferred=True), role("P7"))),
         "The drawn PCR primers complete a single-index Illumina library; the index window's six-base width is measured from Supplementary Fig. 1c artwork and never stated in the text."),
    ))


def _site_seq_workflow():
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


def _digenome_workflow():
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


def _illumina_break_library(name: str, *, inferred: bool, tagged: bool,
                            index_length: int = 6,
                            index_inferred: bool = False) -> tuple[Construct, tuple]:
    """`index_inferred` marks only the sample-index window as an inference.

    DSBCapture publishes its adapters and primers but never states an index length, an
    index read or demultiplexing; the index window is read off figure artwork while the
    adapter arms around it are drawn bases.  Keeping that as a constructor argument stops
    a caller from presenting a figure-measured window as a source-stated one.
    """
    index_inferred = inferred or index_inferred
    left = Segment("P5 + Read 1 arm", il.TRUSEQ_P5_FULL, "p5", inferred=inferred)
    payload = ([Segment("proximal read tag", BLESS_PROXIMAL_READ_TAG, "cbc",
                        feature=feature("proximal_break_end", "feature_barcode", "fixed")),
                role("break-to-random genomic DNA"),
                Segment("distal read tag", BLESS_DISTAL_READ_TAG, "cbc",
                        feature=feature("distal_fragment_end", "feature_barcode", "fixed"))]
               if tagged else [role("DNA beginning at repaired DSB")])
    lib = Construct([left, *payload, Segment("dA junction", "A", inferred=inferred),
                     Segment("Index 1 + Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=inferred),
                     Segment(f"{index_length}-base i7 reverse complement", "I" * index_length, "cbc", placeholder=True,
                             feature=feature("sample_index_i7", "sample_index", "unknown"),
                             inferred=index_inferred),
                     Segment("P7 reverse complement", il.P7_RC, "p7", inferred=inferred)], name=name)
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib, primers


def _finish_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label,
                forward_state=None, reverse_state=None, step_label=None):
    forward = graph.input(prefix + "_forward", label + " PCR forward primer",
                          MaterialKind.OLIGO, state=forward_state)
    reverse = graph.input(prefix + "_reverse", label + " PCR reverse primer",
                          MaterialKind.OLIGO, state=reverse_state)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr",
                step_label or ("Amplify and index the " + label + " library"),
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The final sequencing construct is source-backed, while a single universal PCR-binding geometry is not available for every kit revision."))
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers)))
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library


def _bless_graph(key, spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    sample = graph.input("sample", "yeast cells" if key == "i-bless" else "intact fixed nuclei",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("DSB-bearing sample", (role("genomic DNA with native DSBs"),)),
                         primary=True)
    current = sample
    reviewed = ["sample"]
    if key == "i-bless":
        deproteinization = graph.input("deproteinization_condition",
                                       "agarose immobilization, spheroplasting and deproteinization condition",
                                       MaterialKind.CONDITION)
        protected = graph.product("protected_dna", "agarose-protected deproteinated yeast DNA",
                                  MaterialKind.MOLECULAR, state=sample.state)
        graph.react("immobilize", "Immobilize and deproteinate cells in agarose beads",
                    StepKind.DEPROTEINIZATION,
                    ((sample, InputRole.TEMPLATE),
                     (deproteinization, InputRole.CONDITION)), (protected,))
        current = protected
        reviewed.extend(("deproteinization_condition", "immobilize"))
    repair = graph.input("proximal_repair", "native-break blunting and phosphorylation enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired_breaks", "blunt 5-prime-phosphorylated native breaks",
                             MaterialKind.MOLECULAR, state=current.state)
    graph.react("repair_breaks", "Blunt and 5-prime-phosphorylate native breaks",
                StepKind.END_REPAIR,
                ((current, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    proximal = graph.input("proximal_adapter", "biotinylated proximal hairpin",
                           MaterialKind.ADAPTER)
    proximal_ligase = graph.input("proximal_ligase", "proximal-hairpin DNA ligase",
                                  MaterialKind.CATALYST)
    labelled = graph.product("proximal_labelled", "proximal-hairpin-labelled native breaks",
                             MaterialKind.MOLECULAR,
                             state=molecule("proximal-labelled break", (
                                 role("biotinylated proximal hairpin"),
                                 role("break-adjacent genomic DNA")),
                                 junctions=(("biotinylated proximal hairpin",
                                             "break-adjacent genomic DNA",
                                             "proximal ligation"),)))
    graph.react("ligate_proximal", "Ligate the proximal hairpin in situ",
                StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (proximal, InputRole.STRUCTURAL),
                 (proximal_ligase, InputRole.CATALYST)), (labelled,),
                geometry=reported_adapter_end(
                    "BLESS proximal hairpin", kind="blunt",
                    statement="The repaired native DSB is blunt and 5′ phosphorylated before proximal-hairpin ligation."))
    fragmentation = graph.input("fragmentation_condition", "genomic DNA fragmentation condition",
                                MaterialKind.CONDITION)
    fragments = graph.product("fragments", "proximal-labelled break-to-random fragments",
                              MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("fragment", "Fragment DNA after proximal break labelling",
                StepKind.FRAGMENTATION,
                ((labelled, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)),
                (fragments,))
    streptavidin = graph.input("streptavidin", "streptavidin capture beads",
                               MaterialKind.SELECTION)
    captured = graph.product("captured", "captured proximal-labelled break fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("uncaptured", "unlabelled genomic fragments", MaterialKind.WASTE)
    graph.react("capture", "Capture proximal-labelled molecules on streptavidin",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not retain the biotinylated native-break adapter")
    distal_repair = graph.input("distal_repair", "distal-end blunting enzyme mix",
                                MaterialKind.CATALYST)
    distal_ready = graph.product("distal_ready", "blunt captured distal fragment ends",
                                 MaterialKind.MOLECULAR, state=captured.state)
    graph.react("repair_distal", "Blunt the fragmentation-created distal end",
                StepKind.END_REPAIR,
                ((captured, InputRole.TEMPLATE), (distal_repair, InputRole.CATALYST)),
                (distal_ready,))
    distal = graph.input("distal_adapter", "distinct distal hairpin",
                         MaterialKind.ADAPTER)
    distal_ligase = graph.input("distal_ligase", "distal-hairpin DNA ligase",
                                MaterialKind.CATALYST)
    dual = graph.product("dual", "dual-hairpin break fragments", MaterialKind.MOLECULAR,
                         state=molecule("dual-hairpin break fragment", (
                             role("proximal hairpin"), role("break-to-random genomic DNA"),
                             role("distal hairpin")),
                             junctions=(("proximal hairpin", "break-to-random genomic DNA",
                                         "native-break ligation"),
                                        ("break-to-random genomic DNA", "distal hairpin",
                                         "distal ligation"))))
    graph.react("ligate_distal", "Ligate the distinct distal hairpin",
                StepKind.LIGATION,
                ((distal_ready, InputRole.TEMPLATE), (distal, InputRole.STRUCTURAL),
                 (distal_ligase, InputRole.CATALYST)), (dual,),
                geometry=reported_adapter_end(
                    "BLESS distal hairpin", kind="blunt",
                    statement="The fragmentation-created distal end is blunted before the distinct distal hairpin is ligated."))
    iscei = graph.input("iscei", "I-SceI hairpin-release nuclease", MaterialKind.CATALYST)
    released = graph.product("released", "linearized dual-tag break inserts",
                             MaterialKind.MOLECULAR, state=dual.state)
    graph.react("linearize", "Linearize the captured hairpin product with I-SceI",
                StepKind.DIGESTION,
                ((dual, InputRole.TEMPLATE), (iscei, InputRole.CATALYST)), (released,))
    xhoi = graph.input("xhoi", "XhoI tag-trimming nuclease", MaterialKind.CATALYST)
    trimmed = graph.product("trimmed", "break inserts retaining proximal/distal tags",
                            MaterialKind.MOLECULAR, state=released.state)
    graph.react("trim_tags", "Cut with XhoI while retaining end-identifying tags",
                StepKind.DIGESTION,
                ((released, InputRole.TEMPLATE), (xhoi, InputRole.CATALYST)), (trimmed,))
    reviewed.extend(("proximal_repair", "repair_breaks", "proximal_adapter",
                     "proximal_ligase", "ligate_proximal", "fragmentation_condition",
                     "fragment", "streptavidin", "capture", "distal_repair",
                     "repair_distal", "distal_adapter", "distal_ligase",
                     "ligate_distal", "iscei", "linearize", "xhoi", "trim_tags"))
    _finish_pcr(graph, "library", trimmed, spec.FINAL_LIBRARY, spec.SEQ_PRIMERS,
                spec.SEQUENCING_INTRO, reviewed, label=spec.TITLE)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _dsbcapture_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
    nuclei = graph.input("nuclei", "fixed nuclei containing native DNA breaks",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("fixed nuclei", (role("genomic DNA with native DSBs"),)),
                         primary=True)
    repair = graph.input("repair_enzyme", "native-break blunting enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired", "blunt repaired native breaks",
                             MaterialKind.MOLECULAR, state=nuclei.state)
    graph.react("repair", "Blunt native break ends", StepKind.END_REPAIR,
                ((nuclei, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_enzyme", "single-dA tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "dA-tailed repaired native breaks",
                           MaterialKind.MOLECULAR, state=repaired.state)
    graph.react("tail", "Add one 3-prime dA to native break ends", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    p5 = graph.input("modified_p5", "biotinylated T-overhang modified P5 adapter",
                     MaterialKind.ADAPTER, state=(DSB_P5, DSB_P5_BOTTOM))
    p5_ligase = graph.input("p5_ligase", "modified-P5 adapter ligase", MaterialKind.CATALYST)
    labelled = graph.product("labelled", "modified-P5-labelled native breaks",
                             MaterialKind.MOLECULAR, state=tailed.state)
    graph.react("ligate_p5", "Ligate modified P5 at the native break",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p5_ligase, InputRole.CATALYST)), (labelled,),
                geometry=reported_adapter_end(
                    "DSBCapture modified P5", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The modified P5 adapter carries a single 3′ T complementary to the dA-tailed native break."))
    exonuclease = graph.input("lambda_exonuclease", "lambda exonuclease",
                              MaterialKind.CATALYST)
    cleaned = graph.product("cleaned", "covalently P5-labelled break molecules",
                            MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("remove_free_adapter", "Digest free modified-P5 adapter",
                StepKind.DIGESTION,
                ((labelled, InputRole.TEMPLATE), (exonuclease, InputRole.CATALYST)),
                (cleaned,))
    shear = graph.input("sonication_condition", "DNA sonication condition",
                        MaterialKind.CONDITION)
    fragments = graph.product("fragments", "P5-labelled break-to-sonication fragments",
                              MaterialKind.MOLECULAR, state=cleaned.state)
    graph.react("sonicate", "Sonicate extracted break-labelled DNA", StepKind.FRAGMENTATION,
                ((cleaned, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    streptavidin = graph.input("streptavidin", "streptavidin capture beads",
                               MaterialKind.SELECTION)
    captured = graph.product("captured", "captured modified-P5 break fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("uncaptured", "unlabelled sonication fragments", MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated break fragments",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not retain the biotinylated native-break end")
    distal_repair = graph.input("distal_repair", "distal-end repair enzyme mix",
                                MaterialKind.CATALYST)
    distal_blunt = graph.product("distal_blunt", "blunt sonication-created distal ends",
                                 MaterialKind.MOLECULAR, state=captured.state)
    graph.react("repair_distal", "Blunt the sonication-created distal end",
                StepKind.END_REPAIR,
                ((captured, InputRole.TEMPLATE), (distal_repair, InputRole.CATALYST)),
                (distal_blunt,))
    distal_tailing = graph.input("distal_tailing", "distal single-dA tailing enzyme",
                                 MaterialKind.CATALYST)
    distal_tailed = graph.product("distal_tailed", "dA-tailed distal break fragments",
                                  MaterialKind.MOLECULAR, state=distal_blunt.state)
    graph.react("tail_distal", "Add one 3-prime dA to the distal end", StepKind.TAILING,
                ((distal_blunt, InputRole.TEMPLATE),
                 (distal_tailing, InputRole.CATALYST)), (distal_tailed,))
    p7 = graph.input("modified_p7", "T-overhang modified P7 adapter", MaterialKind.ADAPTER,
                     state=(DSB_P7, DSB_P7_BOTTOM))
    p7_ligase = graph.input("p7_ligase", "modified-P7 adapter ligase", MaterialKind.CATALYST)
    dual = graph.product("dual", "modified-P5/P7 break fragments", MaterialKind.MOLECULAR,
                         state=distal_tailed.state)
    graph.react("ligate_p7", "Ligate modified P7 only at the distal end",
                StepKind.LIGATION,
                ((distal_tailed, InputRole.TEMPLATE), (p7, InputRole.STRUCTURAL),
                 (p7_ligase, InputRole.CATALYST)), (dual,),
                geometry=reported_adapter_end(
                    "DSBCapture modified P7", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The modified P7 adapter carries the T overhang complementary to the dA-tailed sonication end."))
    reviewed = ["nuclei", "repair_enzyme", "repair", "tailing_enzyme", "tail",
                "modified_p5", "p5_ligase", "ligate_p5", "lambda_exonuclease",
                "remove_free_adapter", "sonication_condition", "sonicate", "streptavidin",
                "capture", "distal_repair", "repair_distal", "distal_tailing", "tail_distal",
                "modified_p7", "p7_ligase", "ligate_p7"]
    _finish_pcr(graph, "library", dual, spec.FINAL_LIBRARY, spec.SEQ_PRIMERS,
                spec.SEQUENCING_INTRO, reviewed, label="DSBCapture",
                forward_state=DSB_PCR_F, reverse_state=DSB_PCR_R_AD06,
                step_label="INFERRED — amplify with PCR F and an AD-indexed PCR R; the"
                           " sample-index window is measured from Supplementary Fig. 1c"
                           " artwork, not stated in the text")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


# Every DSBCapture oligo is published only as raster figure artwork (Fig. 1a and
# Supplementary Fig. 1c); the article XML, Supplementary Figures and Supplementary Tables
# contain no base of it, so a verbatim text comparison is impossible rather than failing.
_DSB_FIGURE_ONLY = ("Lensing et al. draw this oligo in Fig. 1a / Supplementary Fig. 1c as"
                    " raster artwork; no base of the DSBCapture oligo set occurs as text in"
                    " the author manuscript, Supplementary Figures or Supplementary Tables,"
                    " so there is no extracted source text to compare against.")
SOURCE_OLIGO_EXCEPTIONS = {
    "modified_p5": _DSB_FIGURE_ONLY,
    "modified_p7": _DSB_FIGURE_ONLY,
    "library_forward": _DSB_FIGURE_ONLY,
    "library_reverse": _DSB_FIGURE_ONLY + " AD006 is one of three drawn example primers.",
}


def _site_seq_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
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
                    spec.FINAL_LIBRARY, reason="Recovery PCR precedes the separately indexed final construct."))
    reviewed = ["dna", "sgrna", "cas9", "cleave", "tailing_enzyme", "tail",
                "adapter1", "ligase1", "ligate_adapter1", "fragmentase", "fragment",
                "distal_repair", "repair_distal", "adapter2", "ligase2", "ligate_adapter2",
                "streptavidin", "capture", "recovery_forward", "recovery_reverse",
                "polymerase", "recovery_pcr"]
    _finish_pcr(graph, "index", recovered, spec.FINAL_LIBRARY, spec.SEQ_PRIMERS,
                spec.SEQUENCING_INTRO, reviewed, label="SITE-Seq indexed")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _digenome_graph(spec):
    source = spec.SOURCE
    graph = ProtocolGraph(spec.TITLE, authority=source)
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
                            MaterialKind.LIBRARY, state=spec.FINAL_LIBRARY)
    graph.react("ligate", "Ligate TruSeq adapters without PCR", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (library,),
                geometry=reported_adapter_end(
                    "TruSeq WGS adapter", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The PCR-free TruSeq adapter T overhang ligates to the dA-tailed WGS fragment."))
    graph.sequence(library, SequencingPlan(
        library=spec.FINAL_LIBRARY, primers=tuple(spec.SEQ_PRIMERS)))
    reviewed = ["dna", "sgrna", "cas9", "cleave", "cleanup_condition", "cleanup",
                "shear_condition", "shear", "repair_enzyme", "repair",
                "size_selection_beads", "size_select", "tailing_enzyme", "tail",
                "truseq_adapters", "ligase", "ligate", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _attach_graph(key, spec):
    if key in {"bless", "i-bless"}:
        spec.protocol_graph = lambda: _bless_graph(key, spec)
    elif key == "dsbcapture":
        spec.protocol_graph = lambda: _dsbcapture_graph(spec)
    elif key == "site-seq":
        spec.protocol_graph = lambda: _site_seq_graph(spec)
    elif key == "digenome-seq":
        spec.protocol_graph = lambda: _digenome_graph(spec)
    return spec


def _spec(key, **fields):
    return _attach_graph(key, SimpleNamespace(**fields))


def get(key: str):
    if key not in DESIGNS and key != "digenome-seq":
        raise KeyError(key)
    if key == "digenome-seq":
        lib, primers = _digenome_library()
        return _spec(key,
            TITLE="Digenome-seq — digested-genome whole-genome sequencing",
            NOTES="01_digenome_seq.html",
            SOURCE='<a href="https://doi.org/10.1038/nmeth.3284">Kim et al., Nature Methods (2015)</a>; library architecture follows the specified Illumina TruSeq DNA PCR-free WGS chemistry.',
            SUMMARY="Cas9–sgRNA cleaves purified genomic DNA in vitro, after which the entire unenriched genome receives an ordinary PCR-free whole-genome library preparation. Coincident forward- and reverse-strand read starts reveal cleavage sites.",
            CAVEAT="Digenome-seq has no cut-site adapter, UMI or enrichment barcode; the assay-specific information is the genomic fragment boundary itself. INFERRED — no source available here gives one adapter base for this library: the defining paper's Methods is paywalled and absent, and its supplement contains no library chemistry. The TruSeq shell, the six-base single index and the index read are inherited from the authors' later Genome Research report and Illumina documentation, so the whole adapter shell is drawn as inferred.",
            oligos=lambda: (), workflow=_digenome_workflow,
            FINAL_LIBRARY=lib,
            FINAL_CAPTION="One ordinary TruSeq PCR-free WGS molecule. Its genomic insert may carry a Cas9-created boundary or only random shearing boundaries; the molecule itself has no label distinguishing them.",
            SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": 150, "Index 1 (i7)": 6, "Read 2": 150},
            SEQUENCING_INTRO="The defining study used Illumina HiSeq X Ten whole-genome sequencing. Digenome analysis depends on paired genomic alignments and the concentration of identical 5′ endpoints, not an inline molecular identifier.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")
    if key == "site-seq":
        lib, primers = _site_seq_library()
        return _spec(key,
            TITLE="SITE-Seq — in-vitro Cas9 cleavage-site sequencing", NOTES="01_site_seq.html",
            SOURCE='<a href="https://doi.org/10.1038/nmeth.4284">Cameron et al., Nature Methods (2017)</a>; exact chemistry from the published supplementary protocol and Supplementary Table 26.',
            SUMMARY="Cas9–sgRNA cleaves purified genomic DNA in vitro; a biotinylated adapter labels those ends before enzymatic fragmentation, allowing streptavidin enrichment and Read 1 entry directly beside each candidate cut site.", CAVEAT="",
            oligos=lambda: (
                f"Adapter 1 Forward — 5′-/5Biosg/{SITE_A1_F[:-1]}*{SITE_A1_F[-1]}-3′",
                f"Adapter 1 Reverse — 5′-{SITE_A1_R}-3′",
                f"Adapter 2 N7 Forward — 5′-/5Phos/{SITE_A2_F7}-3′",
                f"Adapter 2 N6 Forward — 5′-/5Phos/{SITE_A2_F6}-3′",
                f"Adapter 2 N5 Forward — 5′-/5Phos/{SITE_A2_F5}-3′",
                f"Adapter 2 Reverse — 5′-{SITE_A2_R[:-1]}*{SITE_A2_R[-1]}-3′",
                f"Recovery PCR Forward — 5′-{SITE_RECOVERY_F}-3′",
                f"Recovery PCR Reverse — 5′-{SITE_RECOVERY_R}-3′",
                f"Index Forward — 5′-{SITE_INDEX_F}-3′",
                f"Index Reverse — 5′-{SITE_INDEX_R}-3′"),
            workflow=_site_seq_workflow, FINAL_LIBRARY=lib,
            FINAL_CAPTION="N7 Adapter 2 branch; N5 and N6 differ only in diversity-spacer length. Adapter 1 defines the Cas9-cleavage end, while Adapter 2 occupies the later Fragmentase end.",
            SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": 150, "Index 1 (i7)": 8, "Index 2 (i5)": 8},
            SEQUENCING_INTRO="The supplementary protocol recommends 2–3 million 150-base single-end reads per sample. Both published index primers carry eight variable bases; the standard index primers land on the completed dual-index library.",
            SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")

    if key == "dsbcapture":
        lib, primers = _illumina_break_library("DSBCapture library", inferred=False,
                                               tagged=False, index_inferred=True)
        return _spec(key,
            TITLE="DSBCapture — direct in-situ sequencing of DNA breaks", NOTES="01_dsbcapture.html",
            SOURCE='Lensing et al., Nature Methods (2016); <a href="https://doi.org/10.1038/nmeth.3960">doi:10.1038/nmeth.3960</a>. The oligo set below is published only as figure artwork (Fig. 1a, Supplementary Fig. 1c).',
            SUMMARY="A biotinylated modified P5 adapter is ligated to repaired breaks inside fixed nuclei; after sonication and capture, modified P7 is installed only at the random distal end, so Read 1 begins at the DSB.",
            CAVEAT="<b>Oligo provenance.</b> Lensing et al. print the modified P5, modified P7 and control modified P5 duplexes and the PCR primers only as raster artwork in Fig. 1a and Supplementary Fig. 1c: no base of them occurs as text in the author manuscript, Supplementary Figures or Supplementary Tables, so the bases shown here cannot be re-checked against extracted source text. Every drawn base is a canonical Illumina TruSeq part, so this page derives the set from the shared <code>illumina</code> constants and from each duplex's own top strand; a figure transcription and a reconstruction out of those constants are therefore indistinguishable. The sample-index window is measured from the same artwork — the paper never states an index length, an index read or demultiplexing — and is shaded as inferred.",
            oligos=lambda: (
                f"modified P5 top — 5′-/5Biosg/{DSB_P5[0]}*{DSB_P5[1:-1]}*{DSB_P5[-1]}-3′",
                f"modified P5 bottom — 5′-/5Phos/{DSB_P5_BOTTOM}-3′",
                f"modified P7 top — 5′-/5Phos/{DSB_P7}-3′",
                f"modified P7 bottom — 5′-{DSB_P7_BOTTOM[:-1]}*{DSB_P7_BOTTOM[-1]}-3′",
                f"PCR F — 5′-{DSB_PCR_F}-3′",
                f"PCR R (AD006 example, one of three drawn) — 5′-{DSB_PCR_R_AD06}-3′",
                "control modified P5 — Supplementary Fig. 1c also publishes a non-biotinylated"
                " control modified P5 adapter, the negative control of Supplementary Fig. 1e"
                " (lane C). Its bases are not transcribed here: the artwork is the only source,"
                " and whether it differs from the modified P5 above by nothing but the omitted"
                " 5′ biotin is not stated anywhere in the text."),
            workflow=_dsbcapture_workflow, FINAL_LIBRARY=lib,
            FINAL_CAPTION="Published modified-P5/modified-P7 library; the proximal orientation is shown. Read 1 starts at the repaired native break. The shaded sample-index window is the figure-measured inference described above.",
            SEQ_PRIMERS=primers, READ_LENGTHS={"Read 1": 75, "Index 1 (i7)": 6, "Read 2": 75},
            READ_USES={"Index 1 (i7)": "undetermined"},
            SEQUENCING_INTRO="The paper used 2×75-bp NextSeq runs (and 2×50-bp MiSeq for the EcoRV pilot); downstream break calling retains the proximal Read 1. The indexed PCR R implies an index read, but the paper never describes one, so the Index 1 length is taken from the figure's index window and its handling is left undetermined.", SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")

    is_i = key == "i-bless"
    lib, primers = _illumina_break_library(f"{DESIGNS[key].name} library", inferred=True,
                                           tagged=True, index_length=8 if is_i else 6)
    return _spec(key,
        TITLE=("i-BLESS — immobilized Breaks Labeling, Enrichment on Streptavidin and Sequencing" if is_i else
               "BLESS — Breaks Labeling, Enrichment on Streptavidin and Sequencing"),
        NOTES=("01_i_bless.html" if is_i else "01_bless.html"),
        SOURCE=(
            'Biernacka et al., Communications Biology (2018); <a href="https://doi.org/10.1038/s42003-018-0165-9">doi:10.1038/s42003-018-0165-9</a>.' if is_i else
            'Crosetto et al., Nature Methods (2013); <a href="https://doi.org/10.1038/nmeth.2408">doi:10.1038/nmeth.2408</a>.'),
        SUMMARY=(("Agarose beads immobilize yeast cells through lysis and direct break labelling; the otherwise BLESS-like proximal/distal hairpins retain which sequenced end is the native DSB.") if is_i else
                 "A barcoded biotin hairpin labels repaired DSBs in fixed nuclei before fragmentation; streptavidin capture and a differently barcoded distal hairpin preserve the native-break end."),
        CAVEAT=("The proximal/distal tags are exact. A 2016 ThruPLEX guide establishes an eight-base index, but the study does not identify its index tube and the kit's stem-loop adapter oligos are undisclosed. The final adapter geometry is inferred." if is_i else
                "The proximal/distal tags are exact. The paper names the final Illumina library-preparation reagents but does not print that kit's complete adapter geometry; those outer regions are marked inferred."),
        oligos=lambda: (
            f"proximal hairpin — 5′-/5Phos/{BLESS_PROXIMAL[:45]}/iBiodT/{BLESS_PROXIMAL[46:]}-3′",
            f"distal hairpin — 5′-/5Phos/{BLESS_DISTAL}-3′",
            f"proximal PCR primer — 5′-{BLESS_PCR_PROXIMAL}-3′",
            f"distal PCR primer — 5′-{BLESS_PCR_DISTAL}-3′"),
        workflow=lambda: _bless_workflow(key), FINAL_LIBRARY=lib,
        FINAL_CAPTION="One insert orientation is shown. Proximal and distal 11-nt tags identify which genomic end was the native break; inferred outer sequencing-kit regions are shaded.",
        SEQ_PRIMERS=primers, READ_LENGTHS={"Read 1": 75, "Index 1 (i7)": 8 if is_i else 6, "Read 2": 75},
        SEQUENCING_INTRO=("The defining i-BLESS study used 2×75-bp HiSeq 2500/4000 sequencing." if is_i else
                          "BLESS used paired-end Illumina sequencing; either end can carry the proximal tag, which is selected computationally."),
        SEQUENCING_ENDING="", SEQUENCING_UNAVAILABLE="")
