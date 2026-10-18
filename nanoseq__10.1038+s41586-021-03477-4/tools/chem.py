"""Source-owned protocol definition and reaction graph.

Authority: Abascal et al., Nature 593, 405-410 (2021), doi:10.1038/s41586-021-03477-4.
The published article text is paywalled here; the accepted author manuscript
(Cambridge repository 1810/319175) Methods and the published Supplementary Information
(MOESM1, Notes 2 and 4) were read. Exact passages are cited beside each choice.
"""
from __future__ import annotations

from types import SimpleNamespace
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, feature, oligo
from nanoseq_chem import (oligo_scene, a_tailed_fragment, ddbtp_blocked_nick, duplex_adapter_end,
                          hpych4v_fragment, umi,
                          NANOQPCR1, NANOQPCR2)
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from pcr import PCRGeometry
from protocol_model import molecule, role

AM = "Abascal et al. 2021 author manuscript"
SI = "Abascal et al. 2021 published Supplementary Information"

TITLE="NanoSeq — nanorate duplex sequencing"
SOURCE='Abascal et al., <a href="https://doi.org/10.1038/s41586-021-03477-4">Nature (2021)</a>; accepted author manuscript Methods and published SI Notes2/4.'
SUMMARY="HpyCH4V fragmentation, modified dATP/ddBTP A-tailing, duplex-barcoded adapter ligation and a molecular bottleneck precede UDI PCR and duplex consensus analysis."
CAVEAT="Selected HpyCH4V plus ddBTP route, without the optional diversity spacer. Complete historical xGen1080799 adapter strands, UDI primer ends and index widths remain unassigned. The exact qPCR oligos do not establish a full TruSeq/P5/P7 shell. Proposed nick-termination chemistry remains a candidate mechanism (SI Note4)."
SEQUENCING="150bp paired-end HiSeq2500/HiSeqX/NovaSeq. Source read starts NNNTCA or NNNXTCA: 3 random barcode bases, optional spacer, adapter T and HpyCH4V CA. The no-spacer representative is drawn; index acquisition and cycles remain unassigned."
READ_SPECS=(sp.read_spec("Read 1",150,"NNN barcode, adapter T, HpyCH4V CA then genomic sequence; optional spacer source alternative","extract_then_include"),sp.read_spec("Read 2",150,"Opposite duplex barcode and genomic end; optional spacer source alternative","extract_then_include"),sp.read_spec("Index 1 (i7)","unassigned","UDI index; acquired cycles unassigned","undetermined"),sp.read_spec("Index 2 (i5)","unassigned","UDI index; acquired cycles unassigned","undetermined"))


def core_parts():
 return (umi("random barcode A","nanoseq_umi_a"),Segment("adapter T","T"),Segment("HpyCH4V 5′ end","CA"),role("genomic insert; variable length"),Segment("HpyCH4V 3′ end","TG"),Segment("dA tail","A"),Segment("random barcode B complement","nnn","umi",placeholder=True,length_bp=3,feature=feature("nanoseq_umi_b","umi","random")))


def adapted_precursor():
 return molecule("actual ligated native xGen fragment; strand/end occupancy unassigned",(role("native xGen left arm"),*core_parts(),role("native xGen right arm")))


def nanoseq_oligos():
    """qPCR quantification primers printed in the author-manuscript Methods."""
    return (
        oligo("NanoqPCR1 (HPLC)", list(NANOQPCR1)),
        oligo("NanoqPCR2 (HPLC)", list(NANOQPCR2)),
    )


def final_library():
 return Construct([role("UDI PCR left end; native arrangement and index width unassigned",molecular_feature=feature("nanoseq_udi_left","sample_index","unknown",note="Index carried within this unresolved primer-derived end; position, width and bases unassigned")),*core_parts(),role("UDI PCR right end; native arrangement and index width unassigned",molecular_feature=feature("nanoseq_udi_right","sample_index","unknown",note="Index carried within this unresolved primer-derived end; position, width and bases unassigned"))],name="NanoSeq HpyCH4V final target library"),()

FINAL,PRIMERS=final_library()


def _cleanup(graph, ident, label, template, beads, state, kind=MaterialKind.MOLECULAR):
    kept = graph.product(ident, label, kind, state=state)
    waste = graph.product(ident + "_supernatant", label + ": bead supernatant/washes",
                          MaterialKind.WASTE)
    graph.react(ident + "_cleanup", "AMPure XP recovery: " + label, StepKind.CLEANUP,
                ((template, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (kept, waste))
    graph.discard(waste, "unbound reaction components removed with the bead supernatant")
    return kept


def protocol_graph():
    source=SOURCE
    graph=ProtocolGraph(TITLE,authority=source)
    dna = graph.input("dna", "10 ng genomic DNA or LCM lysate in 20 µl", MaterialKind.BIOLOGICAL,
                      state=molecule("genomic DNA", (role("genomic DNA"),)), primary=True,
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
                            source=SI + " Note 4 (proposed mechanism)",inferred=True)
    graph.react("a_tail", "A-tail with dATP/ddBTPs; INFERRED nick-termination candidate (37 °C 30 min)", StepKind.TAILING,
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
                            state=adapted_precursor())
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
    qpcr_mix = graph.input("kapa_qpcr", "KAPA library quantification KK4835 SYBR FAST mix plus supplied primer premix",
                           MaterialKind.CATALYST)
    qpcr1 = graph.input("nanoqpcr1", "NanoqPCR1 (HPLC), added to the KAPA mix", MaterialKind.OLIGO,
                        state=oligo_scene(NANOQPCR1),
                        source=AM + " Methods, DNA quantification")
    qpcr2 = graph.input("nanoqpcr2", "NanoqPCR2 (HPLC), added to the KAPA mix", MaterialKind.OLIGO,
                        state=oligo_scene(NANOQPCR2),
                        source=AM + " Methods, DNA quantification")
    dilution_500=graph.input("qpcr_dilution","1:500 aliquot dilution; 2µl sample in10µl qPCR",MaterialKind.CONDITION)
    aliquot=graph.product("qpcr_aliquot","diluted qPCR aliquot",MaterialKind.MOLECULAR,state=ligated.state)
    quantified=graph.product("quantified","retained main library stock; concentration inferred from aliquot",MaterialKind.MOLECULAR,state=ligated.state)
    graph.react("sample_qpcr","Sample and dilute aliquot; retain main stock",StepKind.SAMPLING,((ligated,InputRole.TEMPLATE),(dilution_500,InputRole.CONDITION)),(aliquot,quantified),changes_molecule=False)
    readout = graph.product("qpcr_readout", "qPCR concentration (fmol/µl, 1:500 dilution)",
                            MaterialKind.QC)
    qpcr_products=graph.product("qpcr_products","quantification PCR amplicons in assay aliquot",MaterialKind.MOLECULAR,state=molecule("qPCR amplicons; native complete ends unassigned",(role("NanoqPCR/KAPA assay amplicon"),)))
    graph.react("qpcr", "qPCR-quantify a 1:500 aliquot", StepKind.PCR,
                ((aliquot, InputRole.TEMPLATE), (qpcr1, InputRole.STRUCTURAL),
                 (qpcr2, InputRole.STRUCTURAL), (qpcr_mix, InputRole.CATALYST)),
                (qpcr_products, readout),
                geometry=PCRGeometry.unavailable(
                    aliquot.state, reason="Adapter strands are not printed, so qPCR primer "
                                          "landing on the ligated molecule cannot be drawn."))
    graph.mark_terminal(Terminal(readout.id,Disposition.QC))
    graph.mark_terminal(Terminal(qpcr_products.id,Disposition.UNRESOLVED,"Analytical aliquot; not the main library stock used for production PCR."))
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
                              state=FINAL)
    graph.react("library_pcr", "Amplify with UDI-containing PCR primers", StepKind.PCR,
                ((diluted, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
                 (rev, InputRole.STRUCTURAL), (q5, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    diluted.state, reason="Actual preceding diluted native adapter-ligated fragment; UDI-containing PCR primers and exact strand ends remain unassigned."))
    final = _cleanup(graph, "library", "two consecutive 0.7x AMPure clean-ups",
                     amplified, beads, FINAL, kind=MaterialKind.LIBRARY)
    graph.sequence(final, SequencingPlan(library=FINAL, acquired_roles=("Read 1","Read 2"),limitation="Source-positive paired150 reads; exact native sequencing-primer sites remain unassigned."))
    graph.declare_reviewed_scope(("dna","ampure","bead_bound_cleanup","hpych4v","digest","fragments_clean_cleanup","klenow","datp_ddbtp","a_tail","duplex_adapters","ligase","ligate","ligated_clean_cleanup","kapa_qpcr","qpcr_dilution","sample_qpcr","nanoqpcr1","nanoqpcr2","qpcr","dilution_target","bottleneck","udi_forward","udi_reverse","q5","library_pcr","library_cleanup","library"),source)
    return graph.require_publishable()


SPEC=SimpleNamespace(TITLE=TITLE,NOTES="01_nanoseq.html",SOURCE=SOURCE,SUMMARY=SUMMARY,CAVEAT=CAVEAT,FINAL_LIBRARY=FINAL,SEQ_PRIMERS=PRIMERS,READ_LENGTHS=None,ACQUIRED_ROLES=("Read 1","Read 2"),SHOW_READ_CONTENT=False,READ_SPECS=READ_SPECS,SEQUENCING_INTRO=SEQUENCING,SEQUENCING_ENDING=SEQUENCING,SEQUENCING_UNAVAILABLE="Native sequencing-primer geometry and UDI index acquisition remain unassigned.",FINAL_CAPTION="Source-reported HpyCH4V/barcode core; native UDI ends are role-level.",oligos=nanoseq_oligos,protocol_graph=protocol_graph)
