"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from chil_chem import AD1_NOMX as CHIL_AD1_NOMX
from chil_chem import AD2_INDEX as CHIL_AD2_INDEX
from chil_chem import CHIL_V1_FW, CHIL_V1_RV, MEDSB_FW, MEDSB_RV
from pcr import PCRGeometry
from pcr import PCRPrimer
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import SequencingPlan
from protocol_graph import StepKind
from rna_library_operations import molecule

DEFINITION = Definition("ChIL-seq — chromatin integration labeling", "10.1038/s41596-020-0375-8", "Handa et al., Nature Protocols (2020)", "An oligo-conjugated antibody targets chromatin; nearby Tn5 integration and T7 transcription linearly amplify target-proximal DNA before sequencing.", "custom", ("ChIL-probe end", "target-proximal DNA"), "Illumina sequencing starts from the ChIL-probe-derived end and reports nearby genomic sequence.", "DNA",
    caveat=("Source scope: the defining Nature Protocols article (Handa et al. 2020, Procedure Steps 11\u201384, "
            "Reagent setup, Figs. 1\u20132) and its Supplementary Information (Table 4 oligos) were read. Drawn: "
            "sequential in situ Tn5 (90 ng/well) then 20 nM MEDS-B loading, 37\u00b0C transposition in TAPS-DMF, "
            "SDS stop, in situ fill-in with T4 DNA polymerase + T4 DNA ligase + dNTPs (RT, 30 min), overnight "
            "in situ T7 transcription, DNase I and RNeasy MinElute RNA recovery, Read2+N11-primed SMARTScribe "
            "first strand (42\u00b0C 90 min, SMART-Seq v4 template-switching oligonucleotide present), SeqAmp PCR "
            "with Ad1_noMX + Ad2_index (12\u201317 cycles; 17 for \u22641k cells), two 1:1 AMPure XP cleanups, "
            "E-Gel 250\u2013350-bp size selection and MinElute, then 51-cycle single-end sequencing. The "
            "SMART-Seq v4 oligonucleotide is a proprietary kit component (bases not printed); its template-switched "
            "cDNA end lies outside the Ad1_noMX priming site and is not drawn. The 8-nt Ad2 index is read for "
            "demultiplexing, but the index-read cycle count is not stated. The optional MinElute pooling of "
            "libraries and the multitarget (mtChIL-seq, v3 probes, Ad1_mt, dual-index) branch are not drawn. "
            "Supplement caption slips that do not change the drawn library: Supplementary Fig. 4e cites "
            "\u2018Supplementary Table 4\u2019 for the Gapdh/Neurod6 promoter primers, which are printed in "
            "Supplementary Table 3; Supplementary Fig. 2 says \u2018in situ transcription\u2019 for a genomic-PCR "
            "transposition assay."))


from protocol_specs import workflow, _payload
from chil_chem import oligos as chil_oligos
from batch_ngs import nextera_library
from chemdraw import Segment


def _source_workflow():
    return workflow(molecule("fixed chromatin", ("antibody target", "nearby DNA")), (
        ("Bind an antibody carrying a ChIL oligonucleotide",
         molecule("targeted chromatin", ("antibody–ChIL probe", "nearby genomic DNA")),
         "The probe contains a T7 promoter and sequence used to recover nearby DNA."),
        ("Tagment nearby genomic DNA in situ",
         molecule("tagged target DNA", ("ChIL probe", "target-proximal genomic DNA", "Tn5 end")),
         "Tn5 strand transfer joins the ChIL-DNA mosaic end to nearby genomic DNA."),
        ("Fill in, transcribe in situ with T7, prime cDNA with Read2+N11, then PCR",
         molecule("ChIL-seq cDNA", ("ChIL-probe end", "target-proximal DNA copy"), joined=True),
         "T4 DNA polymerase/ligase fill-in completes the integrated probe; T7 transcription "
         "from its promoter amplifies the neighbourhood linearly (Procedure Steps 18\u201348)."),
    ))

SPEC = get('chil-seq', definition=DEFINITION,
           workflow_fn=_source_workflow,
           final_override=nextera_library(
               [Segment("target-proximal genomic DNA", "X" * 28, placeholder=True),
                Segment("N11 random-priming site (Read2+N11)", "N" * 11,
                        placeholder=True, length_bp=11)],
               DEFINITION.title + " final library", dual_index=False),
           oligos_fn=chil_oligos)


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    proc = source + "; Procedure Step "
    graph = ProtocolGraph(d.title, authority=source)
    chromatin = graph.input("chromatin", "fixed, permeabilized cells in a 96-well plate",
                            MaterialKind.BIOLOGICAL,
                            state=molecule("fixed chromatin", ("antibody target", "nearby DNA")),
                            primary=True)
    antibody = graph.input("antibody",
                           "primary antibody, then secondary ChIL probe (azido ChIL-DNA on DBCO-IgG) "
                           "at 1 \u00b5g/ml in 0.5 M NaCl, 4\u00b0C overnight",
                           MaterialKind.SELECTION,
                           state=(CHIL_V1_FW, CHIL_V1_RV),
                           source=proc + "7A; Box 3; Supplementary Table 4")
    bound = graph.product("bound", "ChIL-probe-targeted chromatin", MaterialKind.MOLECULAR,
                          state=molecule("targeted chromatin", (
                              "antibody\u2013ChIL probe", "nearby genomic DNA")))
    unbound = graph.product("unbound", "unbound antibody and ChIL probe", MaterialKind.WASTE)
    graph.react("bind", "Bind the ChIL-oligo-conjugated antibody, then wash",
                StepKind.BINDING,
                ((chromatin, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (bound, unbound), source=proc + "7\u20138")
    graph.discard(unbound, "removed with the antibody solution and three PBS washes (Step 8)")
    tn5 = graph.input("tn5", "Tn5 transposase, 90 ng/well, added first (10 min RT)",
                      MaterialKind.CATALYST, source=proc + "11")
    transposome = graph.input("transposome",
                              "20 nM annealed MEDS-B (MEDSB-Fw/Rv), added after Tn5 (1 h RT); "
                              "the bound ChIL-DNA supplies the other mosaic end",
                              MaterialKind.ADAPTER,
                              state=(MEDSB_FW, MEDSB_RV),
                              source=proc + "12; Reagent setup; Supplementary Table 4")
    taps = graph.input("taps_dmf", "1\u00d7 TAPS-DMF buffer (Mg2+), 37\u00b0C 1 h; then 0.2% SDS 10 min",
                       MaterialKind.CONDITION, source=proc + "13\u201317")
    tagged = graph.product("tagged", "ChIL-DNA integrated into target-proximal genomic DNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("tagged target DNA", (
                               "ChIL probe", "target-proximal genomic DNA", "Tn5 end")))
    waste = graph.product("tn5_waste", "unbound Tn5 and MEDS-B, removed by washes",
                          MaterialKind.WASTE)
    graph.react("tagment", "Load Tn5 and MEDS-B in situ, then transpose at 37\u00b0C",
                StepKind.TAGMENTATION,
                ((bound, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST), (taps, InputRole.CONDITION)), (tagged, waste),
                source=proc + "11\u201317")
    graph.discard(waste, "supernatant and PBS/dialysis-buffer washes discarded (Steps 13\u201314, 17)")
    # deferred-followup-03: the defining Procedure (Steps 18-22) replaces the earlier
    # MAGIC boundary with an in situ fill-in by T4 DNA polymerase + T4 DNA ligase.
    # The molecular states stay role-level: the gap/overhang lengths at the
    # integration junction are not printed, so no base geometry is asserted.
    fill_mix = graph.input("fill_in_enzymes",
                           "T4 DNA polymerase (0.015 U/\u00b5l) + T4 DNA ligase (2 U/\u00b5l), "
                           "0.1 mM dNTPs in 1\u00d7 T4 DNA ligase buffer, RT 30 min; then 0.2% SDS",
                           MaterialKind.CATALYST, source=proc + "18\u201322")
    extended = graph.product("transcription_ready",
                             "filled-in, sealed integrated ChIL probe ready for in situ T7 transcription",
                             MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("post_transposition", "Fill in and seal the integration junctions in situ",
                StepKind.EXTENSION, ((tagged, InputRole.TEMPLATE),
                                     (fill_mix, InputRole.CATALYST)), (extended,),
                source=proc + "18\u201322")
    t7 = graph.input("t7_polymerase",
                     "T7 RNA polymerase (2 U/\u00b5l), 2 mM each NTP, RNase inhibitor, 37\u00b0C overnight",
                     MaterialKind.CATALYST, source=proc + "23\u201325")
    rna = graph.product("rna", "linearly amplified target-proximal RNA copy",
                        MaterialKind.MOLECULAR, state=extended.state)
    graph.react("transcribe", "Linearly amplify with in situ T7 transcription",
                StepKind.TRANSCRIPTION,
                ((extended, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (rna,),
                source=proc + "23\u201325")
    dnase = graph.input("dnase", "recombinant DNase I, 37\u00b0C 30 min in the transcription well",
                        MaterialKind.CATALYST, source=proc + "26")
    rna_treated = graph.product("rna_dnase", "DNase-treated ChIL RNA in well supernatant",
                                MaterialKind.MOLECULAR, state=rna.state)
    dna_digest = graph.product("dna_digest", "DNase I digestion products of template DNA",
                               MaterialKind.MOLECULAR)
    graph.react("dnase_digest", "Digest DNA with DNase I and recover the supernatant",
                StepKind.DIGESTION, ((rna, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)),
                (rna_treated, dna_digest), source=proc + "26\u201327")
    rneasy = graph.input("rneasy", "RNeasy MinElute column (RLT + ethanol, RPE, 80% ethanol; "
                         "10 \u00b5l water, re-eluted)", MaterialKind.SELECTION,
                         source=proc + "28\u201341")
    rna_pure = graph.product("rna_purified", "purified ChIL RNA", MaterialKind.MOLECULAR,
                             state=rna.state)
    rna_flow = graph.product("rna_column_flow", "RNeasy flow-through and washes, including DNase I digestion products",
                             MaterialKind.WASTE)
    graph.react("rna_cleanup", "Purify ChIL RNA on an RNeasy MinElute column",
                StepKind.CLEANUP, ((rna_treated, InputRole.TEMPLATE),
                                   (dna_digest, InputRole.TEMPLATE),
                                   (rneasy, InputRole.SELECTION)),
                (rna_pure, rna_flow), source=proc + "28\u201341")
    graph.discard(rna_flow, "column flow-through and washes discarded (Steps 31\u201336)")
    rt_primer = graph.input("library_rt_primer",
                            "Read2+N11 first-strand primer (0.96 \u00b5M in annealing mix; 72\u00b0C 3 min)",
                            MaterialKind.OLIGO,
                            source=proc + "43\u201345; Supplementary Table 4")
    tso = graph.input("library_tso",
                      "SMART-Seq v4 oligonucleotide (template-switching oligo; proprietary, bases not printed)",
                      MaterialKind.OLIGO, source=proc + "46")
    rt = graph.input("library_rt_enzyme",
                     "SMARTScribe reverse transcriptase, 42\u00b0C 90 min then 70\u00b0C 10 min",
                     MaterialKind.CATALYST, source=proc + "46\u201348")
    cdna = graph.product("library_cdna", "ChIL target-proximal cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("ChIL cDNA", ("target-proximal cDNA",)))
    graph.react("library_rt", "Reverse-transcribe ChIL RNA with Read2+N11",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna_pure, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (tso, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (cdna,),
                source=proc + "43\u201348")
    forward = graph.input("library_forward", "Ad1_noMX PCR primer (0.24 \u00b5M)", MaterialKind.OLIGO,
                          state=CHIL_AD1_NOMX,
                          source=proc + "49A; Supplementary Table 4")
    reverse = graph.input("library_reverse", "Ad2_index PCR primer (8-nt index; 0.24 \u00b5M)",
                          MaterialKind.OLIGO,
                          source=proc + "49A; Supplementary Table 4")
    pcr_polymerase = graph.input("library_polymerase",
                                 "SeqAmp DNA polymerase; 95\u00b0C 1 min, 12\u201317 \u00d7 (98\u00b0C 10 s, "
                                 "65\u00b0C 30 s, 68\u00b0C 3 min), 72\u00b0C 10 min",
                                 MaterialKind.CATALYST, source=proc + "49A")
    library = graph.product("library_library", "ChIL target-proximal sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Amplify and index the ChIL-seq library", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (pcr_polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    final, PCRPrimer("Ad1_noMX", CHIL_AD1_NOMX),
                    PCRPrimer("Ad2_index", CHIL_AD2_INDEX),
                    statement="Supplementary Table 4 prints both complete PCR primers."),
                source=proc + "49A")
    ampure = graph.input("ampure", "AMPure XP beads 1:1 (with 10\u00d7 lysis buffer), twice; "
                         "27 \u00b5l EB", MaterialKind.SELECTION, source=proc + "50\u201361")
    cleaned = graph.product("library_cleaned", "AMPure-purified ChIL library",
                            MaterialKind.LIBRARY, state=final)
    bead_waste = graph.product("ampure_supernatant", "AMPure supernatants and ethanol washes",
                               MaterialKind.WASTE)
    graph.react("ampure_cleanup", "Purify the library twice with 1:1 AMPure XP",
                StepKind.CLEANUP, ((library, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
                (cleaned, bead_waste), source=proc + "50\u201361")
    graph.discard(bead_waste, "bead supernatants and washes discarded (Steps 52\u201355)")
    egel = graph.input("egel", "2% E-Gel Size Select II, recover 250\u2013350 bp; MinElute, 12 \u00b5l EB",
                       MaterialKind.SELECTION, source=proc + "72\u201381")
    selected = graph.product("library_selected", "size-selected ChIL library",
                             MaterialKind.LIBRARY, state=final)
    off_size = graph.product("off_size", "fragments outside 250\u2013350 bp and column washes",
                             MaterialKind.WASTE)
    graph.react("size_selection", "Size-select 250\u2013350 bp on E-Gel and purify on MinElute",
                StepKind.SELECTION, ((cleaned, InputRole.TEMPLATE), (egel, InputRole.SELECTION)),
                (selected, off_size), source=proc + "72\u201381")
    graph.discard(off_size, "left in the gel or column flow-through (Steps 72\u201378)")
    graph.sequence(selected, SequencingPlan(
        library=final, primers=tuple(primers),
        acquired_roles=("Read 1", "Index 1 (i7)"),
        read_lengths={"Read 1": 51}))
    reviewed = ["chromatin", "antibody", "bind", "tn5", "transposome", "taps_dmf", "tagment",
                "fill_in_enzymes", "post_transposition", "t7_polymerase", "transcribe",
                "dnase", "dnase_digest", "rna_dnase", "rneasy", "rna_cleanup", "library_rt_primer",
                "library_tso", "library_rt_enzyme", "library_rt",
                "library_forward", "library_reverse", "library_polymerase",
                "library_pcr", "library_library", "ampure", "ampure_cleanup",
                "egel", "size_selection", "library_selected"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
