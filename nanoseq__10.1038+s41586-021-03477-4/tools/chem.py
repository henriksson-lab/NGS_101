"""Source-owned protocol definition and reaction graph.

Authority: Abascal et al., Nature 593, 405-410 (2021), doi:10.1038/s41586-021-03477-4.
The published article text is paywalled here; the accepted author manuscript
(Cambridge repository 1810/319175) Methods and the published Supplementary Information
(MOESM1, Notes 2 and 4) were read. Exact passages are cited beside each choice.
"""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, feature, oligo, revcomp
from nanoseq_chem import (oligo_scene, a_tailed_fragment, ddbtp_blocked_nick, duplex_adapter_end,
                          hpych4v_fragment, umi, INSERT_N,
                          NANOQPCR1, NANOQPCR2)
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from pcr import PCRGeometry
from protocol_specs import Definition, get
from rna_library_operations import molecule

AM = "Abascal et al. 2021 author manuscript"
SI = "Abascal et al. 2021 published Supplementary Information"

DEFINITION = Definition(
    "NanoSeq — nanorate duplex sequencing", "10.1038/s41586-021-03477-4",
    "Abascal et al., Nature (2021)",
    "HpyCH4V restriction fragmentation avoids end repair, and A-tailing with dATP plus "
    "dideoxy non-A nucleotides (ddBTPs) makes nick-extended strands unamplifiable before "
    "duplex-barcoded adapter ligation and duplex consensus sequencing.",
    "truseq", ("3-nt random barcode", "genomic insert", "3-nt random barcode"),
    "NanoSeq and BotSeqS libraries were sequenced as 150 bp paired-end reads on HiSeq2500, "
    "HiSeqX and NovaSeq. Each read begins NNNTCA (or NNNXTCA): a 3-nt random barcode, an "
    "optional spacer base, the adapter T and the HpyCH4V CA end; barcodes and genomic "
    "coordinates group Watson and Crick read families into duplex consensus.",
    "DNA",
    "INFERRED adapter and PCR shell: the paper names xGen Duplex Seq Adapters (IDT 1080799) "
    "and dried UDI-containing PCR primers but prints neither. The TruSeq-type Read 1/Read 2 "
    "arms, P5/P7 and 8-nt index windows drawn in the final library are inferred from the "
    "printed qPCR primers NanoqPCR1/2 (which match TruSeq Read 1 and Read 2 arm bases) and "
    "from the read layout; the optional spacer base is not drawn. The published article "
    "text itself was not available here; the accepted author manuscript was read.")


def nanoseq_oligos():
    """qPCR quantification primers printed in the author-manuscript Methods."""
    return (
        oligo("NanoqPCR1 (HPLC)", list(NANOQPCR1)),
        oligo("NanoqPCR2 (HPLC)", list(NANOQPCR2)),
    )


def final_library():
    """HpyCH4V NanoSeq library; adapter/PCR shell inferred (see DEFINITION caveat)."""
    lib = Construct([
        Segment("P5", il.P5, "p5", inferred=True),
        Segment("i5", "J" * 8, "cbc", placeholder=True, inferred=True,
                feature=feature("sample_index_i5", "sample_index", "unknown")),
        Segment("Read 1 adapter arm", il.TRUSEQ_READ1, "r1", inferred=True),
        umi("random barcode A", "nanoseq_umi_a"),
        Segment("adapter T", "T"),
        Segment("HpyCH4V 5' end", "CA"),
        Segment("genomic insert", "X" * INSERT_N, placeholder=True),
        Segment("HpyCH4V 3' end", "TG"),
        Segment("dA tail", "A"),
        Segment("random barcode B (complement)", "nnn", "umi", placeholder=True, length_bp=3,
                feature=feature("nanoseq_umi_b", "umi", "random")),
        Segment("Read 2 adapter arm (rc)", revcomp(il.TRUSEQ_READ2), "r2", inferred=True),
        Segment("i7 reverse complement", "I" * 8, "cbc", placeholder=True, inferred=True,
                feature=feature("sample_index_i7", "sample_index", "unknown")),
        Segment("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name="NanoSeq final library (HpyCH4V)")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("; ".join(problems))
    return lib, primers


FINAL, PRIMERS = final_library()


SPEC = get('nanoseq', definition=DEFINITION,
           final_override=(FINAL, PRIMERS), oligos_fn=nanoseq_oligos)


def _cleanup(graph, ident, label, template, beads, state, kind=MaterialKind.MOLECULAR):
    kept = graph.product(ident, label, kind, state=state)
    waste = graph.product(ident + "_supernatant", label + ": bead supernatant/washes",
                          MaterialKind.WASTE)
    graph.react(ident + "_cleanup", "AMPure XP recovery: " + label, StepKind.CLEANUP,
                ((template, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (kept, waste))
    graph.discard(waste, "unbound reaction components removed with the bead supernatant")
    return kept


def _local_graph(d):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("dna", "10 ng genomic DNA or LCM lysate in 20 µl", MaterialKind.BIOLOGICAL,
                      state=molecule("genomic DNA", ("genomic DNA",)), primary=True,
                      source=AM + " Methods, NanoSeq library preparation")
    beads = graph.input("ampure", "AMPure XP beads (50:50 water mix; later 2.5x, 1x, 2x 0.7x)",
                        MaterialKind.SELECTION)
    bound = _cleanup(graph, "bead_bound", "bead-bound input DNA (eluted, beads retained)",
                     dna, beads, None)
    hpy = graph.input("hpych4v", "HpyCH4V (5 U/µl, 0.5 µl) in CutSmart, 37 °C 15 min",
                      MaterialKind.CATALYST, source=AM + " Methods; " + SI + " Note 2")
    fragments = graph.product("fragments", "blunt HpyCH4V fragments (5'-CA…TG-3')",
                              MaterialKind.MOLECULAR, state=hpych4v_fragment())
    graph.react("digest", "On-bead HpyCH4V fragmentation", StepKind.DIGESTION,
                ((bound, InputRole.TEMPLATE), (hpy, InputRole.CATALYST)), (fragments,))
    frag_clean = _cleanup(graph, "fragments_clean", "2.5x AMPure-purified fragments",
                          fragments, beads, fragments.state)
    klenow = graph.input("klenow", "Klenow fragment (3'→5' exo-) in NEBuffer 4",
                         MaterialKind.CATALYST)
    nts = graph.input("datp_ddbtp", "1 mM equimolar dATP/ddBTPs (ddCTP, ddGTP, ddTTP)",
                      MaterialKind.CONDITION, source=AM + " Methods; " + SI + " Note 4")
    tailed = graph.product("a_tailed", "3'-dA-tailed HpyCH4V fragments", MaterialKind.MOLECULAR,
                           state=a_tailed_fragment())
    blocked = graph.product("ddbtp_blocked", "candidate: nick-extended strand ending in a ddBTP",
                            MaterialKind.MOLECULAR, state=ddbtp_blocked_nick(),
                            source=SI + " Note 4 (proposed mechanism)")
    graph.react("a_tail", "A-tail with dATP plus ddBTPs (37 °C 30 min)", StepKind.TAILING,
                ((frag_clean, InputRole.TEMPLATE), (klenow, InputRole.CATALYST),
                 (nts, InputRole.CONDITION)), (tailed, blocked))
    graph.mark_terminal(Terminal(blocked.id, Disposition.UNRESOLVED,
                                 "Not physically removed: the ddBTP-terminated strand is carried "
                                 "into ligation but is unamplifiable (SI Note 4); its other "
                                 "products are not modelled."))
    adapters = graph.input("duplex_adapters", "xGen Duplex Seq Adapters (IDT 1080799), 15 µM",
                           MaterialKind.ADAPTER, state=duplex_adapter_end(),
                           source=AM + " Methods (adapter strands not printed)")
    ligase = graph.input("ligase", "T4 DNA ligase (400 U/µl) + ATP, NEBuffer 4",
                         MaterialKind.CATALYST)
    adapted = graph.product("adapted", "adapter-ligated duplex-barcoded fragments",
                            MaterialKind.MOLECULAR,
                            state=molecule("adapter-ligated fragment", (
                                "adapter + NNN + T", "CA…genomic insert…TG", "A + nnn + adapter"),
                                joined=True), inferred=True)
    graph.react("ligate", "Ligate duplex adapters (20 °C 20 min)", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.known(
                    "cohesive", a_tailed_fragment(), duplex_adapter_end(),
                    statement="3' dA-tailed fragment ends meet the adapter 3' T overhang that "
                              "follows the 3-nt random barcode (AM read-processing Methods: "
                              "'T is the adapter overhang'); remaining adapter strands unprinted."))
    ligated = _cleanup(graph, "ligated_clean", "1x AMPure-purified ligation product",
                       adapted, beads, adapted.state)
    qpcr_mix = graph.input("kapa_qpcr", "KAPA library quantification kit (KK4835) SYBR FAST mix",
                           MaterialKind.CATALYST)
    qpcr1 = graph.input("nanoqpcr1", "NanoqPCR1 (HPLC), added to the KAPA mix", MaterialKind.OLIGO,
                        state=oligo_scene(NANOQPCR1),
                        source=AM + " Methods, DNA quantification")
    qpcr2 = graph.input("nanoqpcr2", "NanoqPCR2 (HPLC), added to the KAPA mix", MaterialKind.OLIGO,
                        state=oligo_scene(NANOQPCR2),
                        source=AM + " Methods, DNA quantification")
    readout = graph.product("qpcr_readout", "qPCR concentration (fmol/µl, 1:500 dilution)",
                            MaterialKind.QC)
    quantified = graph.product("quantified", "quantified adapter-ligated library",
                               MaterialKind.MOLECULAR, state=ligated.state)
    graph.react("qpcr", "qPCR-quantify a 1:500 aliquot", StepKind.PCR,
                ((ligated, InputRole.TEMPLATE), (qpcr1, InputRole.STRUCTURAL),
                 (qpcr2, InputRole.STRUCTURAL), (qpcr_mix, InputRole.CATALYST)),
                (quantified, readout), changes_molecule=False,
                geometry=PCRGeometry.unavailable(
                    ligated.state, reason="Adapter strands are not printed, so qPCR primer "
                                          "landing on the ligated molecule cannot be drawn."))
    graph.mark_terminal(Terminal(readout.id, Disposition.QC))
    diluted = graph.product("diluted", "library diluted to ~0.3 fmol in 25 µl (bottleneck)",
                            MaterialKind.MOLECULAR, state=ligated.state)
    target = graph.input("dilution_target", "volume from qPCR fmol/µl and r_opt ≈ 5.1 read pairs "
                         "per fragment (typically 0.3 fmol for a 15x run)", MaterialKind.CONDITION,
                         source=AM + " Methods, library dilution and sequencing efficiency")
    remainder = graph.product("undiluted_remainder", "ligated library not taken into the dilution",
                              MaterialKind.MOLECULAR, state=ligated.state)
    graph.react("bottleneck", "Dilute to the target fmol (typically 0.3 fmol)",
                StepKind.SAMPLING, ((quantified, InputRole.TEMPLATE), (target, InputRole.CONDITION)),
                (diluted, remainder), changes_molecule=False)
    graph.mark_terminal(Terminal(remainder.id, Disposition.UNRESOLVED,
                                 "Library not taken into the diluted PCR input; its storage or "
                                 "reuse (e.g. undiluted matched-normal libraries) is not assigned "
                                 "per sample."))
    fwd = graph.input("udi_forward", "UDI-containing PCR primer (i5 side; not printed)",
                      MaterialKind.OLIGO, inferred=True)
    rev = graph.input("udi_reverse", "UDI-containing PCR primer (i7 side; not printed)",
                      MaterialKind.OLIGO, inferred=True)
    q5 = graph.input("q5", "NEBNext Ultra II Q5 Master Mix; 98/65 °C cycling, 10–16 cycles by input",
                     MaterialKind.CATALYST)
    amplified = graph.product("amplified", "UDI-indexed NanoSeq library", MaterialKind.LIBRARY,
                              state=FINAL, inferred=True)
    graph.react("library_pcr", "Amplify with UDI-containing PCR primers", StepKind.PCR,
                ((diluted, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
                 (rev, InputRole.STRUCTURAL), (q5, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    FINAL, reason="The UDI-containing PCR primers are named but not printed."))
    final = _cleanup(graph, "library", "two consecutive 0.7x AMPure clean-ups",
                     amplified, beads, FINAL, kind=MaterialKind.LIBRARY)
    graph.sequence(final, SequencingPlan(library=FINAL, primers=PRIMERS,
                                         read_lengths={"Read 1": 150, "Read 2": 150}))
    graph.declare_reviewed_scope(tuple(graph.materials) + tuple(graph.steps), source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION)


SPEC.protocol_graph = protocol_graph
