from chemdraw import Construct, Scene, Segment, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from protocol_model import role
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind,Highlight
import seqprimers as sp
from pcr import PCRGeometry, PCRPrimer
from rage_seq_chem import R1_PRIMER, TSO_PRIMER, oligos as rage_seq_oligos

# Identifier lengths from Methods "Demultiplexing nanopore sequencing data":
# 16-nt cell barcode; the 13 nt downstream are removed to drop the 10-nt UMI.
CELL_BARCODE = feature("chromium_3p_v2_cell_barcode", "cell_barcode", "whitelist",
                       whitelist="matched CellRanger 2.0 cell-barcode list")
UMI = feature("chromium_3p_v2_umi", "umi", "random")


def _amplicon_segments() -> list[Segment]:
    """Full-length 10x 3′ v2 cDNA as amplified by the published TSO and R1 primers.

    The R1 primer is the 5′ end of the strand carrying the barcode; the TSO primer
    anneals to the opposite end. Barcode then UMI order and lengths come from the
    demultiplexing Methods; oligo(dT) length and the remaining TSO-derived bases are
    not printed in the paper and stay role-level.
    """
    return [Segment("R1 primer site", R1_PRIMER, "r1"),
            Segment("16-nt cell barcode", "B" * 16, "cbc", placeholder=True, feature=CELL_BARCODE),
            Segment("10-nt UMI", "U" * 10, "umi", placeholder=True, feature=UMI),
            role("oligo(dT) from the 10x v2 bead primer"),
            Segment("cDNA interval (receptor or non-receptor; antisense to mRNA)", "X" * 24, "insert", placeholder=True),
            role("remaining template-switch-derived bases unassigned"),
            Segment("TSO primer site", revcomp(TSO_PRIMER), "tso")]


def _amplicon(name: str) -> Construct:
    return Construct(_amplicon_segments(), name=name)


def _long_library() -> Construct:
    return Construct([role("baselineSQK-LSK108 adapter; native assembly unassigned",inferred=True), *_amplicon_segments(),
                      role("baselineSQK-LSK108 opposite end; native assembly unassigned",inferred=True)],
                     name="RAGE-Seq long-read receptor library")


def _duplex(c: Construct) -> Scene:
    sc = Scene.duplex(list(c), label=c.name,
                      unpaired=tuple(seg.name + "'" for seg in c if seg.is_role_token()))
    sc.strand("carried background",[role("co-retained off-target/unconsumed cDNA or probe species; quantitative recovery unassigned")],label="independent carried molecular background")
    return sc


TITLE="RAGE-Seq — repertoire and gene expression by sequencing"
NOTES="01_rage_seq.html"
SOURCE='Singh et al., NatureCommunications(2019), <a href="https://doi.org/10.1038/s41467-019-11049-4">doi:10.1038/s41467-019-11049-4</a>; Methods droplet RNA, capture design, targeted capture, Nanopore sequencing and demultiplexing.'
SUMMARY="Source-amplified Chromium3′v2 cDNA is aliquoted to matched Illumina expression and double-hybrid-capture long-read receptor libraries."
CAVEAT="This graph begins with the source-amplified3′v2 cDNA pool, whose initial commercial RT/amplification is a declared upstream boundary. Exact R1/TSO reamplification primers and16ntbarcode/10ntUMI come from the defining paper; oligo(dT), remaining native TSO, full capture-bound first-cycle geometry and historicalONT end preparation/adapter assembly remain unassigned. Physical hybrid capture enrichment does not establish pure receptor recovery. LSK1081D is baseline; the one LSK3081D2 sample was analysed as1D, not merged into one adapter chemistry."
LONG_LIBRARY=_long_library()
SHORT_LIBRARY=Construct([role("source Chromium3′v2 Illumina adapter/index assembly unassigned"),Segment("source16ntcellbarcode","N"*16,placeholder=True,feature=CELL_BARCODE),Segment("source10ntUMI","N"*10,placeholder=True,feature=UMI),role("matched transcriptome fragment; full end geometry unassigned")],name="RAGE matched expression library")
FINAL_LIBRARIES=(("Long-read receptor library",LONG_LIBRARY,(),"Native adapter assembly bounded; source R1/TSO and identifier lineage retained.","BaselineLSK108 onR9.4.1; statedLSK308/R9.5.1 sample alternatives remain separate."),("Matched expression library",SHORT_LIBRARY,(),"Matched3′v2 transcriptome aliquot; historical library finishing remains source-named boundary.","SourceNextSeq500 “150bps,paired-end”; per-read/index cycle allocation unassigned."))
SEQUENCING_UNAVAILABLE="Historical kit/end geometry and unassigned primer sites are not supplied by generic endpoints."
READ_SPECS={"Long-read receptor library":(sp.read_spec("Sequencing output",None,"MinION1D receptor cDNA with16ntcellbarcode/10ntUMI; identifiers extracted computationally","extract_then_include"),),"Matched expression library":(sp.read_spec("Read 1","unassigned","source NextSeq500150bp paired-end statement; exact per-read allocation/native geometry unassigned","undetermined"),sp.read_spec("Read 2","unassigned","source NextSeq500150bp paired-end statement; exact per-read allocation/native geometry unassigned","undetermined"),sp.read_spec("Index 1 (i7)","unassigned","physical sample index; acquisition/cycles unassigned","undetermined"),sp.read_spec("Index 2 (i5)","unassigned","physical index architecture and acquisition unassigned","undetermined"))}

def oligos():return rage_seq_oligos()

def _pcr(graph, ident, label, template, cycles, conc, reviewed):
    fwd = graph.input(ident + "_r1", f"R1 primer ({conc})", MaterialKind.OLIGO,
                      state=Construct([Segment("R1 primer", R1_PRIMER)], name="R1 primer"))
    rev = graph.input(ident + "_tso", f"TSO primer ({conc})", MaterialKind.OLIGO,
                      state=Construct([Segment("TSO primer", TSO_PRIMER)], name="TSO primer"))
    pol = graph.input(ident + "_kapa", "KAPA HotStart HiFi ReadyMix", MaterialKind.CATALYST)
    amp = _amplicon(label)
    out = graph.product(ident + "_product", label, MaterialKind.MOLECULAR, state=_duplex(amp))
    graph.react(ident, f"PCR with TSO and R1 primers ({cycles})", StepKind.PCR,
                ((template, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
                 (rev, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (out,),
                geometry=PCRGeometry.unavailable(template.state,reason="Actual incoming source-amplified or bead-captured cDNA is retained. Printed reamplification primers define product-end references, not complete native input/first-cycle geometry."))
    return out


def _capture(graph, ident, label, template, reviewed):
    probes = graph.input(ident + "_probes", "SeqCap RNA Choice probes (644 exons) + SeqCap EZ hybridization/wash reagents, no HE oligos",
                         MaterialKind.SELECTION)
    captured = graph.product(ident + "_captured", label, MaterialKind.MOLECULAR,
                             state=captured_state(label))
    rejected = graph.product(ident + "_washed", "unbound capture fractions and washes; exact specificity unassigned", MaterialKind.WASTE)
    graph.react(ident, "Hybridize overnight at 47 °C and capture on streptavidin beads",
                StepKind.AFFINITY_CAPTURE,
                ((template, InputRole.TEMPLATE), (probes, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "not retained by the biotinylated probes; removed in the SeqCap EZ washes")
    return captured


def protocol_graph():
    source=SOURCE
    final=LONG_LIBRARY
    graph = ProtocolGraph(TITLE, authority=source)
    cdna_state = _duplex(_amplicon("amplified 10x 3′ v2 full-length cDNA"))
    cdna = graph.input("cdna", "Chromium 3′ v2 full-length cDNA (2 extra amplification cycles)",
                       MaterialKind.BIOLOGICAL, state=cdna_state, primary=True)
    expression = graph.product("expression_branch", "30 µL (tumor 35 µL) cDNA for the Illumina library",
                               MaterialKind.MOLECULAR, state=cdna_state)
    receptor = graph.product("receptor_branch", "remaining cDNA for targeted capture",
                             MaterialKind.MOLECULAR, state=cdna_state)
    graph.react("split", "Split expression and immune-receptor branches",
                StepKind.SPLIT, ((cdna, InputRole.TEMPLATE),), (expression, receptor),
                changes_molecule=False)
    reviewed=()
    pre = _pcr(graph, "pre_capture_pcr", "re-amplified cDNA", receptor,
               "5 cycles cell lines / 20 primary cells", "3 µM", reviewed)
    beads = graph.input("ampure", "AMPure XP beads", MaterialKind.SELECTION)
    cleaned = graph.product("pre_capture_clean", "purified cDNA (500 ng–1 µg into capture)",
                            MaterialKind.MOLECULAR, state=pre.state)
    spent = graph.product("pre_capture_spent", "primers, nucleotides and reagents not bound by beads",
                          MaterialKind.WASTE)
    graph.react("ampure_cleanup", "AMPure XP purification", StepKind.CLEANUP,
                ((pre, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (cleaned, spent),
                changes_molecule=False)
    graph.discard(spent, "left in the AMPure XP supernatant and washes")
    c1 = _capture(graph, "capture1", "first-round captured receptor cDNA", cleaned, reviewed)
    p1 = _pcr(graph, "post_capture1_pcr", "first-round enriched cDNA", c1, "5 cycles", "1 µM", reviewed)
    c2 = _capture(graph, "capture2", "second-round captured receptor cDNA", p1, reviewed)
    p2 = _pcr(graph, "post_capture2_pcr", "double-captured receptor cDNA (0.6–2 kb)", c2, "20 cycles", "1 µM", reviewed)
    motor = graph.input("motor_adapter", "baselineONT SQK-LSK1081D adapter; sourceLSK308samplealternative not substituted",
                        MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "kit adapter ligation reagents", MaterialKind.CATALYST)
    long_library = graph.product("long_library", "RAGE-Seq long-read receptor library",
                                 MaterialKind.LIBRARY, state=final)
    graph.react("ligate_motor", "Prepare the ONT 1D ligation library (kit protocol)",
                StepKind.LIGATION,
                ((p2, InputRole.TEMPLATE), (motor, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (long_library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "BaselineLSK108 kit end-preparation/motor joining remain unassigned; original kit protocols were not followed within this capped source scope. The separateLSK3081D2 sample is not treated as equivalent chemistry."))
    graph.sequence(long_library, SequencingPlan(library=LONG_LIBRARY,acquisition="unspecified",limitation="SourceMinION1D:R9.4.1FLO-MIN106,threecell-lineR9.5.1FLO-MIN107; Albacore2.2.7. HistoricalONTadapter/motor/endassembly bounded."))
    short=graph.product("expression_library",SHORT_LIBRARY.name,MaterialKind.LIBRARY,state=SHORT_LIBRARY)
    graph.react("expression_library_prep", "Standard Chromium 3′ v2 shearing and Illumina library construction",
                StepKind.MAGIC, ((expression, InputRole.TEMPLATE),), (short,),
                reason="Methods states the remainder of the Chromium 3′ v2 protocol was followed (shearing, Illumina library construction); this page does not redraw that vendor chemistry.")
    graph.sequence(short, SequencingPlan(library=SHORT_LIBRARY,acquired_roles=("Read 1","Read 2"),limitation="NextSeq500 source150bps paired-end, >50000 raw reads/cell; exact per-read/index cycle allocation and native read geometry unassigned."))
    graph.declare_reviewed_scope(REVIEWED_TARGETS, source)
    return graph.require_publishable()


def captured_state(label):
    sc=Scene();sc.strand("captured cDNA",_amplicon_segments(),label=label+" bead-selected cDNA strand class; native pairing unassigned")
    sc.strand("capture probes",[role("biotinylated SeqCap probes; exact bound-strand geometry unassigned")],label="independent capture oligo population")
    sc.strand("background",[role("carried off-target cDNA; quantitative capture specificity unassigned")],label="independent retained background")
    return sc

REVIEWED_TARGETS=("cdna","split","pre_capture_pcr_r1","pre_capture_pcr_tso","pre_capture_pcr_kapa","pre_capture_pcr","ampure","ampure_cleanup","capture1_probes","capture1","post_capture1_pcr_r1","post_capture1_pcr_tso","post_capture1_pcr_kapa","post_capture1_pcr","capture2_probes","capture2","post_capture2_pcr_r1","post_capture2_pcr_tso","post_capture2_pcr_kapa","post_capture2_pcr","motor_adapter","adapter_ligase","ligate_motor","long_library","expression_library_prep","expression_library")
