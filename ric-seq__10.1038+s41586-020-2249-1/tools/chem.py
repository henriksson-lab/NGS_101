"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from interaction_capture import FixedRNPContact
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from pcr import PCRGeometry
from protocol_graph import SequencingPlan
from rna_library_operations import molecule
from adapter_ligation import AdapterLigationGeometry

DEFINITION = Definition("RIC-seq — RNA in situ conformation sequencing", "10.1038/s41586-020-2249-1", "Cai et al., Nature (2020)", "Crosslinked RNA–protein assemblies are fragmented and proximity-ligated in situ to capture global spatial RNA–RNA contacts.", "custom", ("RNA A cDNA", "RNA B cDNA"), "Paired-end reads identify the two RNA fragments joined at each in-situ ligation junction; the 2020 Methods state HiSeq X Ten paired-end sequencing at Novogene; read lengths are not in the inspected texts.", "RNA", caveat=("Reaction steps follow the Cai 2020 Nature Methods ('Construction of RIC-seq libraries'): 1% formaldehyde fixation with 0.125 M glycine quench, permeabilization, in-situ micrococcal nuclease fragmentation, FastAP removal of 3′ phosphates, T4 RNA ligase addition of pCp–biotin (overnight, 16 °C, PEG 20000), a second FastAP treatment to remove the pCp 3′ phosphate, T4 PNK 5′ phosphorylation, in-situ proximity ligation with T4 RNA ligase (overnight, 16 °C), proteinase K/TRIzol LS recovery, DNase I, optional Ribo-off rRNA depletion, 94 °C Mg2+ fragmentation, MyOne C1 enrichment, N6/SuperScript II RT, dUTP second strand, end repair, A-tailing, ligation of Illumina Y-shaped adaptors, USER digestion with Pfx PCR (11–13 cycles) and MinElute gel purification of 200–450 bp products. The Methods print no bases for the Y-shaped adaptors, and no read lengths. The PEI adapter A/B printed in the authors' later 2021 Nature Protocols Table 1 are not assumed for the 2020 libraries, so the final library stays role-level. Supplementary Table 9 prints PE1.0 and twelve P7 index primers whose index windows are 8 nt (primers 1–7) and 6 nt (primers 8–12) as printed; the 2021 protocol's index primer instead prints a 6-N window."))


from protocol_specs import workflow
from interaction_capture import FixedRNPContact


def _source_workflow():
    contact = FixedRNPContact("formaldehyde crosslinking")
    return workflow(molecule("fixed cells", ("spatially proximal RNAs",)), (
        ("Cross-link RNA–protein assemblies in situ and fragment RNA",
         contact.complex(),
         "Formaldehyde preserves spatial contacts, including protein-mediated ones."),
        ("Label RNA ends and proximity-ligate them in situ",
         contact.chimera(),
         "The ligation junction records a spatial RNA–RNA contact."),
        ("Reverse crosslinks, enrich chimeras and make a paired-end library",
         molecule("RIC-seq library", ("left adapter", "RNA A cDNA", "RNA B cDNA", "right adapter"), joined=True,
                  junctions=(("left adapter", "RNA A cDNA", "adapter ligation"),
                             ("RNA A cDNA", "RNA B cDNA", "copied proximity-ligation junction"),
                             ("RNA B cDNA", "right adapter", "adapter ligation"))),
         "Paired sequence identifies both sides of the RNA contact."),
    ))

SPEC = get('ric-seq', definition=DEFINITION,
           workflow_fn=_source_workflow)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    topology = FixedRNPContact("formaldehyde crosslinking")
    cells = graph.input("cells", "cells containing spatially proximal RNA–protein assemblies",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("fixed-cell input", ("spatially proximal RNAs",)), primary=True)
    formaldehyde = graph.input("crosslinker", "1% formaldehyde in PBS, 10 min room temperature; 0.125 M glycine quench "
                               "(Cai 2020 Methods; GEO GSE127188 GSM3629915 treatment protocol)",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed", "crosslinked RNA–protein neighbourhoods",
                          MaterialKind.MOLECULAR, state=topology.complex())
    graph.react("crosslink", "Crosslink RNA–protein assemblies in situ",
                StepKind.CROSS_LINKING,
                ((cells, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (fixed,))
    permeabilization = graph.input("permeabilization", "permeabilization buffer (NP-40, Triton X-100, "
                                   "Tween 20, SUPERase In), 15 min on ice; PNK-buffer washes "
                                   "(GEO GSE127188 GSM3629915 treatment protocol)", MaterialKind.CONDITION)
    permeable = graph.product("permeabilized", "permeabilized crosslinked cells",
                              MaterialKind.MOLECULAR, state=topology.complex())
    graph.react("permeabilize", "Permeabilize the crosslinked cells", StepKind.PERMEABILIZATION,
                ((fixed, InputRole.TEMPLATE), (permeabilization, InputRole.CONDITION)), (permeable,))
    fragmentation = graph.input("mnase", "micrococcal nuclease 0.03 U/µl, 10 min 37 °C; stopped by "
                                "EGTA washes (GEO GSE127188 GSM3629915 treatment protocol)",
                                MaterialKind.CATALYST)
    fragments = graph.product("fragments", "fragmented fixed RNA neighbourhoods",
                              MaterialKind.MOLECULAR, state=topology.complex())
    graph.react("fragment", "Fragment RNA in situ with micrococcal nuclease", StepKind.FRAGMENTATION,
                ((permeable, InputRole.TEMPLATE), (fragmentation, InputRole.CATALYST)), (fragments,))
    meth = "Cai 2020 Methods, 'pCp–biotin labelling and proximity ligation'"
    fastap1 = graph.input("fastap_1", f"FastAP alkaline phosphatase 10 U, 37 °C 15 min; PNK+EGTA, high-salt "
                          f"and PNK washes ({meth})", MaterialKind.CATALYST)
    dephos = graph.product("dephosphorylated", "fixed RNA fragments with 3′-OH ends (3′ phosphates removed)",
                           MaterialKind.MOLECULAR, state=topology.complex())
    graph.react("dephosphorylate_3p", "Remove MNase-generated 3′ phosphates", StepKind.DEPHOSPHORYLATION,
                ((fragments, InputRole.TEMPLATE), (fastap1, InputRole.CATALYST)), (dephos,))
    pcp = graph.input("pcp_biotin", "pCp–biotin 4 µl of 1 mM (Thermo Fisher 20160)", MaterialKind.OLIGO)
    pcp_ligase = graph.input("pcp_ligase", f"T4 RNA ligase 100 U, RNA ligase buffer, RNase inhibitor, "
                             f"15% PEG 20000, overnight 16 °C; PNK washes ({meth})", MaterialKind.CATALYST)
    labelled = graph.product("labelled", "fixed RNA fragments with 3′-terminal pCp–biotin (3′ phosphate)",
                             MaterialKind.MOLECULAR, state=topology.complex())
    graph.react("label_ends", "Ligate pCp–biotin to RNA 3′ ends in situ (omitted in −pCp controls, ED Fig. 1b,c)",
                StepKind.LIGATION,
                ((dephos, InputRole.TEMPLATE), (pcp, InputRole.STRUCTURAL), (pcp_ligase, InputRole.CATALYST)),
                (labelled,), geometry=AdapterLigationGeometry.unavailable(
                    "Single-nucleotide 3′-end label: the Methods give T4 RNA ligase joining pCp–biotin "
                    "to dephosphorylated RNA 3′-OH ends, but the RNA ends are random MNase fragments "
                    "with no defined bases, so no base-level junction is drawn."))
    fastap2 = graph.input("fastap_2", f"FastAP alkaline phosphatase 10 U, 15 min ({meth})", MaterialKind.CATALYST)
    cp_oh = graph.product("cp_biotin_oh", "fixed RNA fragments ending in Cp–biotin with a 3′-OH",
                          MaterialKind.MOLECULAR, state=topology.complex())
    graph.react("dephosphorylate_pcp", "Remove the 3′ phosphate of the added pCp–biotin",
                StepKind.DEPHOSPHORYLATION,
                ((labelled, InputRole.TEMPLATE), (fastap2, InputRole.CATALYST)), (cp_oh,))
    pnk = graph.input("pnk", f"T4 polynucleotide kinase (EK0032) with 1.5 mM ATP in PNK buffer A, 37 °C 45 min; "
                      f"PNK+EGTA and PNK washes ({meth})", MaterialKind.CATALYST)
    kinased = graph.product("kinased", "fixed RNA fragments with 5′-phosphate and 3′ Cp–biotin-OH ends",
                            MaterialKind.MOLECULAR, state=topology.complex())
    graph.react("phosphorylate_5p", "Phosphorylate RNA 5′ ends", StepKind.PHOSPHORYLATION,
                ((cp_oh, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (kinased,))
    ligase = graph.input("proximity_ligase", f"T4 RNA ligase 100 U (Thermo Fisher EL0021), RNA ligase buffer, "
                         f"RNase inhibitor, 0.1 mg/ml BSA in 200 µl, overnight 16 °C ({meth})", MaterialKind.CATALYST)
    chimera = graph.product("chimera", "in-situ proximity-ligated RNA chimera with internal C–biotin at the junction",
                            MaterialKind.MOLECULAR, state=topology.chimera())
    graph.react("ligate", "Proximity-ligate labelled RNA ends in situ",
                StepKind.PROXIMITY_LIGATION,
                ((kinased, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (chimera,))
    reversal = graph.input("reversal", "proteinase K in 0.5% SDS buffer, 37 °C 60 min then 56 °C 15 min; "
                           "TRIzol LS/chloroform extraction (GEO GSE127188 GSM3629915 extract protocol)",
                           MaterialKind.CONDITION)
    purified = graph.product("purified", "de-crosslinked RNA chimeras",
                             MaterialKind.MOLECULAR, state=chimera.state)
    graph.react("reverse_crosslinks", "Reverse crosslinks and recover RNA chimeras",
                StepKind.DEPROTEINIZATION,
                ((chimera, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (purified,))
    dnase = graph.input("dnase", "RQ1 DNase I, 37 °C 20 min; acid phenol:chloroform, ethanol "
                        "precipitation (GEO GSE127188 GSM3629915 extract protocol)", MaterialKind.CATALYST)
    dnased = graph.product("dnase_treated", "DNase-treated total RNA including chimeras",
                           MaterialKind.MOLECULAR, state=chimera.state)
    dna_waste = graph.product("dna_digested", "digested DNA, including RNA–DNA ligation by-products",
                              MaterialKind.WASTE)
    graph.react("dnase_digest", "Digest DNA and RNA–DNA ligation by-products", StepKind.DIGESTION,
                ((purified, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)), (dnased, dna_waste))
    graph.discard(dna_waste, "digested by DNase I")
    depletion = graph.input("rrna_depletion", "Ribo-off rRNA Depletion Kit (Vazyme N406-01); optional "
                            "per GEO — GSE127188 holds rRNA_depleted and Total libraries",
                            MaterialKind.SELECTION)
    depleted = graph.product("depleted", "rRNA-depleted total RNA including chimeras",
                             MaterialKind.MOLECULAR, state=chimera.state)
    rrna = graph.product("rrna_removed", "depleted rRNA", MaterialKind.WASTE)
    graph.react("deplete_rrna", "Optionally deplete rRNA (ED Fig. 1a, top; GEO extract protocol)",
                StepKind.SELECTION,
                ((dnased, InputRole.TEMPLATE), (depletion, InputRole.SELECTION)), (depleted, rrna))
    graph.discard(rrna, "removed by rRNA depletion")
    refrag = graph.input("rna_fragmentation", "94 °C for 5 min in 1× first-strand buffer (Mg2+), then ice "
                         "(GEO GSE127188 GSM3629915 extract protocol)", MaterialKind.CONDITION)
    sized = graph.product("sized", "fragmented rRNA-depleted RNA", MaterialKind.MOLECULAR,
                          state=chimera.state)
    graph.react("fragment_rna", "Fragment the recovered RNA (ED Fig. 1a, middle)",
                StepKind.FRAGMENTATION,
                ((depleted, InputRole.TEMPLATE), (refrag, InputRole.CONDITION)), (sized,))
    beads = graph.input("c1_beads", "MyOne Streptavidin C1 beads in TWB, 30 min room temperature; "
                        "eluted 3× at 95 °C in PK buffer (GEO GSE127188 GSM3629915 extract protocol)",
                        MaterialKind.SELECTION)
    enriched = graph.product("enriched", "pCp–biotin-labelled RNA fragments after C1 pull-down",
                             MaterialKind.MOLECULAR, state=chimera.state)
    unbound = graph.product("unbound", "unlabelled RNA fragments in C1 flow-through",
                            MaterialKind.WASTE)
    graph.react("c1_pulldown", "Enrich biotin-labelled RNA on C1 beads (ED Fig. 1a, bottom)",
                StepKind.AFFINITY_CAPTURE,
                ((sized, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (enriched, unbound))
    graph.discard(unbound, "not retained by streptavidin C1 beads")
    reviewed = ["cells", "crosslinker", "crosslink", "permeabilization", "permeabilize",
                "mnase", "fragment", "fastap_1", "dephosphorylated", "dephosphorylate_3p",
                "pcp_ligase", "fastap_2", "cp_biotin_oh", "dephosphorylate_pcp", "pnk", "kinased",
                "phosphorylate_5p", "dnase", "dnase_treated", "dnase_digest",
                "pcp_biotin", "label_ends", "proximity_ligase", "ligate", "reversal",
                "reverse_crosslinks", "rrna_depletion", "deplete_rrna", "rna_fragmentation",
                "fragment_rna", "c1_beads", "c1_pulldown"]
    rt_primer = graph.input("library_rt_primer", "N6 random primer (Cai 2020 Methods; GEO GSE127188 GSM3629915)",
                            MaterialKind.OLIGO)
    rt = graph.input("library_rt_enzyme", "SuperScript II reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "RIC-seq chimera cDNA", MaterialKind.MOLECULAR,
                         state=molecule("RIC-seq chimera cDNA", ("cDNA insert",)))
    graph.react("library_rt", "Reverse-transcribe the enriched RNA", StepKind.REVERSE_TRANSCRIPTION,
                ((enriched, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    second = graph.input("library_second_strand", "RNase H, E. coli DNA polymerase I and dNTPs with "
                         "80% of dTTP replaced by dUTP, 16 °C 2 h (GEO extract protocol)",
                         MaterialKind.CATALYST)
    ds = graph.product("library_ds_cdna", "double-stranded cDNA with dUTP-marked second strand",
                       MaterialKind.MOLECULAR,
                       state=molecule("ds cDNA (dUTP second strand)", ("cDNA insert",)))
    graph.react("library_second_strand_synthesis", "Synthesize the dUTP-marked second strand",
                StepKind.EXTENSION, ((cdna, InputRole.TEMPLATE), (second, InputRole.CATALYST)), (ds,))
    repair = graph.input("library_end_repair", "T4 DNA polymerase, Klenow DNA polymerase and T4 PNK "
                         "after 1.8× AMPure (GEO extract protocol)", MaterialKind.CATALYST)
    repaired = graph.product("library_repaired", "end-repaired ds cDNA", MaterialKind.MOLECULAR,
                             state=ds.state)
    graph.react("library_end_repair_step", "End-repair the ds cDNA", StepKind.END_REPAIR,
                ((ds, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    klenow = graph.input("library_a_tailing", "Klenow exo− (GEO extract protocol)", MaterialKind.CATALYST)
    tailed = graph.product("library_a_tailed", "3′-A-tailed ds cDNA", MaterialKind.MOLECULAR,
                           state=ds.state)
    graph.react("library_a_tail", "Add 3′ adenines", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (klenow, InputRole.CATALYST)), (tailed,))
    adaptors = graph.input("library_adapters", "Illumina Y-shaped adaptors, 2 pmol (Cai 2020 Methods; "
                           "bases and ends not printed there)", MaterialKind.ADAPTER)
    ligase_dna = graph.input("library_adapter_ligase", "Quick T4 DNA ligase, 20 °C 15 min; AMPure ×2 "
                             "(GEO extract protocol)", MaterialKind.CATALYST)
    adapted = graph.product("library_adapted", "adapter-flanked chimera cDNA", MaterialKind.MOLECULAR,
                            state=molecule("adapter-flanked cDNA", ("left adapter", "cDNA insert", "right adapter")))
    graph.react("library_adapter_install", "Ligate paired-end adapters to the A-tailed cDNA",
                StepKind.LIGATION, ((tailed, InputRole.TEMPLATE), (adaptors, InputRole.STRUCTURAL),
                                    (ligase_dna, InputRole.CATALYST)),
                (adapted,), geometry=AdapterLigationGeometry.unavailable(
                    "The Cai 2020 Methods name 'Illumina Y-shaped adaptors (2 pmol)' after "
                    "Klenow exo- A-tailing but print no adaptor bases or ends. The PEI adapter A/B "
                    "printed in Cao 2021 Nat Protoc Table 1 belong to that later protocol and are "
                    "not assumed for the 2020 libraries."))
    user = graph.input("library_user", "USER enzyme, 37 °C 15 min in the PCR mix to digest the dUTP "
                       "strand (GEO extract protocol)", MaterialKind.CATALYST)
    forward = graph.input("library_forward", "Illumina PE1.0 primer (Supplementary Table 9)",
                          MaterialKind.OLIGO)
    reverse = graph.input("library_reverse", "Index primer 1–12 (Supplementary Table 9)",
                          MaterialKind.OLIGO)
    polymerase = graph.input("library_polymerase", "Platinum Pfx DNA polymerase, 11–13 cycles (GEO extract protocol)", MaterialKind.CATALYST)
    amplified = graph.product("library_amplified", "PCR-amplified RIC-seq library",
                              MaterialKind.MOLECULAR, state=final)
    graph.react("library_pcr", "Amplify and index the RIC-seq library", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST),
                 (user, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    final, reason="Supplementary Table 9 prints PE1.0 and the index primers, but "
                                  "the Cai 2020 Methods print no bases for the Y-shaped adaptors "
                                  "they amplify."))
    gel = graph.input("gel", "agarose gel, 200–450 bp products, MinElute Gel Extraction (Qiagen 28606) "
                      "(Cai 2020 Methods; ED Fig. 1b,c)",
                      MaterialKind.SELECTION)
    library = graph.product("library_library", "RIC-seq sequencing library", MaterialKind.LIBRARY,
                            state=final)
    off_size = graph.product("gel_excluded", "library fragments outside the excised gel region",
                             MaterialKind.WASTE)
    graph.react("gel_purify", "Gel-purify the library (ED Fig. 1b,c)", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (library, off_size))
    graph.discard(off_size, "outside excised gel region")
    graph.sequence(library, SequencingPlan(library=final, primers=tuple(primers)) if primers
                   else SequencingPlan(limitation=d.sequencing))
    reviewed.extend(("library_rt_primer", "library_rt_enzyme", "library_rt",
                     "library_second_strand", "library_ds_cdna", "library_second_strand_synthesis",
                     "library_end_repair", "library_repaired", "library_end_repair_step",
                     "library_a_tailing", "library_a_tailed", "library_a_tail",
                     "library_adapter_ligase", "library_user", "library_adapters",
                     "library_adapter_install", "library_forward", "library_reverse",
                     "library_polymerase", "library_pcr", "library_amplified", "gel",
                     "gel_purify", "library_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
