"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import SequencingPlan
from protocol_graph import StepKind
from protocol_specs import role

DEFINITION = Definition("RAGE-Seq — repertoire and gene expression by sequencing", "10.1038/s41467-019-11049-4", "Singh et al., Nature Communications (2019)", "Hybrid capture enriches full-length immune-receptor cDNAs from barcoded droplet libraries for long-read sequencing alongside short-read expression.", "nanopore", ("motor adapter", "cell barcode", "full-length receptor cDNA"), "MinION 1D reads (R9.4.1 FLO-MIN106; three cell-line samples R9.5.1 FLO-MIN107) retain the 16-nt cell barcode, 10-nt UMI and full receptor isoform; barcodes are matched in the first or last 200 nt of reads >250 nt against the short-read barcode list. The Oxford Nanopore adapter sequences are not printed in the paper.", "RNA", "Exact TSO and R1 amplification primers are printed in Methods. The 644-exon Roche–NimbleGen capture panel is specified by target genes rather than probe sequences, and the SQK-LSK108/SQK-LSK308 adapter sequences are not printed in the paper.")


from chemdraw import Construct, Scene, Segment, feature, revcomp
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
            Segment("receptor cDNA (antisense to mRNA)", "X" * 24, "insert", placeholder=True),
            role("template-switch-derived bases not printed"),
            Segment("TSO primer site", revcomp(TSO_PRIMER), "tso")]


def _amplicon(name: str) -> Construct:
    return Construct(_amplicon_segments(), name=name)


def _long_library() -> Construct:
    return Construct([role("SQK-LSK108 / SQK-LSK308 adapter"), *_amplicon_segments(),
                      role("SQK-LSK108 / SQK-LSK308 adapter (other end)")],
                     name="RAGE-Seq long-read receptor library")


def _duplex(c: Construct) -> Scene:
    sc = Scene.duplex(list(c), label=c.name,
                      unpaired=tuple(seg.name + "'" for seg in c if seg.is_role_token()))
    return sc


SPEC = get('rage-seq', definition=DEFINITION,
           final_override=_long_library(),
           oligos_fn=rage_seq_oligos,
           sequencing_mechanism="nanopore_motor",
           motor_adapter_name="SQK-LSK108 / SQK-LSK308 motor adapter")
SPEC.CAVEAT = ("Exact bases: the TSO and R1 amplification primers printed in Methods. "
               "Cell-barcode/UMI lengths and order come from the demultiplexing Methods; "
               "the oligo(dT), remaining template-switch bases, the Roche–NimbleGen probe "
               "pool and the Oxford Nanopore adapters are role-level.")


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
                geometry=PCRGeometry.exact(amp, PCRPrimer("R1 primer", R1_PRIMER),
                                           PCRPrimer("TSO primer", TSO_PRIMER),
                                           statement="Methods, Targeted capture: TSO and R1 primers printed; they flank the full-length 10x 3′ v2 cDNA."))
    reviewed.extend((fwd.id, rev.id, pol.id, ident, out.id))
    return out


def _capture(graph, ident, label, template, reviewed):
    probes = graph.input(ident + "_probes", "SeqCap RNA Choice probes (644 exons) + SeqCap EZ hybridization/wash reagents, no HE oligos",
                         MaterialKind.SELECTION)
    captured = graph.product(ident + "_captured", label, MaterialKind.MOLECULAR,
                             state=_duplex(_amplicon(label)))
    rejected = graph.product(ident + "_washed", "non-receptor cDNA washed off", MaterialKind.WASTE)
    graph.react(ident, "Hybridize overnight at 47 °C and capture on streptavidin beads",
                StepKind.AFFINITY_CAPTURE,
                ((template, InputRole.TEMPLATE), (probes, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "not retained by the biotinylated probes; removed in the SeqCap EZ washes")
    reviewed.extend((probes.id, ident, captured.id, rejected.id))
    return captured


def _local_graph(d, final):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
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
    reviewed = ["cdna", "split"]
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
    reviewed.extend(("ampure", "ampure_cleanup", "pre_capture_clean", "pre_capture_spent"))
    c1 = _capture(graph, "capture1", "first-round captured receptor cDNA", cleaned, reviewed)
    p1 = _pcr(graph, "post_capture1_pcr", "first-round enriched cDNA", c1, "5 cycles", "1 µM", reviewed)
    c2 = _capture(graph, "capture2", "second-round captured receptor cDNA", p1, reviewed)
    p2 = _pcr(graph, "post_capture2_pcr", "double-captured receptor cDNA (0.6–2 kb)", c2, "20 cycles", "1 µM", reviewed)
    motor = graph.input("motor_adapter", "ONT SQK-LSK108 1D (one sample SQK-LSK308) adapter",
                        MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "kit adapter ligation reagents", MaterialKind.CATALYST)
    long_library = graph.product("long_library", "RAGE-Seq long-read receptor library",
                                 MaterialKind.LIBRARY, state=final)
    graph.react("ligate_motor", "Prepare the ONT 1D ligation library (kit protocol)",
                StepKind.LIGATION,
                ((p2, InputRole.TEMPLATE), (motor, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (long_library,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The paper names SQK-LSK108/SQK-LSK308 but not their end-preparation or adapter-end chemistry; the kit protocols were not read in this pass."))
    graph.sequence(long_library, SequencingPlan(limitation=d.sequencing))
    reviewed.extend(("motor_adapter", "adapter_ligase", "ligate_motor", "long_library"))
    short = graph.product("expression_library", "10x 3′ v2 Illumina expression library (not drawn)",
                          MaterialKind.LIBRARY,
                          state=Construct([role("Chromium 3′ v2 Illumina library")],
                                          name="RAGE-Seq matched expression library"))
    graph.react("expression_library_prep", "Standard Chromium 3′ v2 shearing and Illumina library construction",
                StepKind.MAGIC, ((expression, InputRole.TEMPLATE),), (short,),
                reason="Methods states the remainder of the Chromium 3′ v2 protocol was followed (shearing, Illumina library construction); this page does not redraw that vendor chemistry.")
    graph.sequence(short, SequencingPlan(limitation="NextSeq 500, 150-bp paired-end, >50,000 raw reads per cell (Methods); read structure follows the 10x 3′ v2 guide, not redrawn here."))
    reviewed.extend(("expression_library_prep", "expression_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY)

SPEC.protocol_graph = protocol_graph
