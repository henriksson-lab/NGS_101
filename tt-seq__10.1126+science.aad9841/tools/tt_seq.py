"""Schwalb et al. 2016 TT-seq molecular workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "TT-seq — transient transcriptome sequencing"
NOTES = "01_tt-seq.html"
SOURCE = ('Defining paper: <a href="https://doi.org/10.1126/science.aad9841">'
          'Schwalb et al., <i>Science</i> (2016)</a>.')
SUMMARY = ("A five-minute 4-thiouridine pulse marks newly made RNA. Fragmenting RNA before "
           "thiol capture makes only the newly synthesized pieces recoverable for sequencing.")
CAVEAT = ("The paper names the NuGEN Ovation Human Blood RNA-seq kit but does not print its "
          "adapter oligos. The Illumina-compatible library below is therefore an inferred "
          "run geometry; the TT-seq label, fragmentation and capture steps are source-defined.")

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5", inferred=True),
    seg("Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=True),
    seg("strand-specific cDNA from captured RNA", "X" * 38, placeholder=True),
    seg("dA junction", "A", inferred=True),
    seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=True),
    seg("six-base i7 index reverse complement", "I" * 6, "cbc",
        placeholder=True, inferred=True,
        feature=feature("sample_i7", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
], name="inferred single-index TT-seq library")
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                         required_roles=tuple(p.role for p in SEQ_PRIMERS)):
    raise ValueError("TT-seq final library: " + "; ".join(problems))

FINAL_CAPTION = ("INFERRED — Illumina-compatible paired-end, six-base single-index geometry; "
                 "the defining supplement reports paired 50-base reads plus a six-base barcode.")
SEQUENCING_INTRO = ("The inferred standard primer sites show how the named HiSeq run enters the "
                    "captured cDNA. Adapter bases are not asserted as NuGEN kit disclosures.")


def sections():
    return [
        ("Pulse-label nascent RNA", [
            Row(chunks=[("pre-existing RNA     --------------------------", None, False)]),
            Row(chunks=[("new RNA (5 min)      ----4sU---------4sU-------", "umi", False)]),
        ], "Cells receive 500 µM 4sU for five minutes; labelled spike-ins are added during extraction."),
        ("Fragment before selection", [
            Row(chunks=[("long labelled RNA -- BioRuptor 30 s ON / 30 s OFF --> shorter labelled RNA pieces", None, False)]),
        ], "Fragmentation precedes labelled-RNA purification; that ordering is the distinguishing TT-seq step."),
        ("Biotinylate and capture labelled pieces", [
            Row(chunks=[("RNA—4sU  + thiol-reactive biotin  ->  RNA—4sU—biotin", "umi", False)]),
            Row(chunks=[("RNA—4sU—biotin  -> streptavidin capture -> DTT elution", None, False)]),
        ], "Only fragments carrying a pulse-labelled thiol are enriched; unlabeled RNA is washed away."),
        ("INFERRED — convert captured RNA to the sequencing library", [
            Row(chunks=[("DNase treatment -> strand-specific NuGEN Ovation library preparation -> paired-end library", None, True)]),
        ], "The defining supplement delegates this construct-changing step to the commercial kit."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    cells=graph.input("cells","cells during the five-minute pulse",MaterialKind.BIOLOGICAL,
                      primary=True)
    s4u=graph.input("s4u","4-thiouridine pulse label",MaterialKind.MOLECULAR)
    rnap=graph.input("rna_polymerase","cellular RNA polymerases",MaterialKind.CATALYST)
    labelled=graph.product("labelled_rna","4sU-labelled newly synthesized RNA",
                           MaterialKind.MOLECULAR)
    graph.react("pulse_label","Transcribe RNA during the five-minute 4sU pulse",
                StepKind.TRANSCRIPTION,
                ((cells,InputRole.TEMPLATE),(s4u,InputRole.STRUCTURAL),
                 (rnap,InputRole.CATALYST)),(labelled,))
    sonication=graph.input("fragmentation_condition","BioRuptor fragmentation condition",
                           MaterialKind.CONDITION)
    fragmented=graph.product("fragmented","fragmented 4sU-labelled RNA",
                             MaterialKind.MOLECULAR)
    graph.react("fragment","Fragment RNA before selection",StepKind.FRAGMENTATION,
                ((labelled,InputRole.TEMPLATE),(sonication,InputRole.CONDITION)),
                (fragmented,))
    biotin=graph.input("biotin_reagent","thiol-reactive biotinylation condition",
                       MaterialKind.CONDITION)
    biotinylated=graph.product("biotinylated","biotinylated 4sU-containing RNA pieces",
                               MaterialKind.MOLECULAR)
    graph.react("biotinylate","Covalently biotinylate 4sU thiols",
                StepKind.CHEMICAL_CONVERSION,
                ((fragmented,InputRole.TEMPLATE),(biotin,InputRole.CONDITION)),
                (biotinylated,))
    streptavidin=graph.input("streptavidin","streptavidin selection reagent",
                            MaterialKind.SELECTION)
    captured=graph.product("captured","captured newly synthesized RNA fragments",
                           MaterialKind.MOLECULAR)
    unlabelled=graph.product("unlabelled","unlabelled RNA fragments",MaterialKind.WASTE)
    graph.react("capture","Capture biotinylated RNA fragments",StepKind.AFFINITY_CAPTURE,
                ((biotinylated,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,unlabelled))
    graph.discard(unlabelled,"unlabelled RNA is removed by streptavidin washing")
    rt_primer=graph.input("rt_primer","NuGEN strand-specific RT primer",MaterialKind.OLIGO,
                          inferred=True,source=SOURCE+" (commercial kit oligo undisclosed)")
    rt=graph.input("reverse_transcriptase","NuGEN reverse transcriptase",
                   MaterialKind.CATALYST)
    cdna=graph.product("cdna","strand-specific cDNA from captured RNA",
                       MaterialKind.MOLECULAR,inferred=True,
                       source=SOURCE+" (commercial kit construct)")
    graph.react("reverse_transcribe","Convert captured RNA to strand-specific cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    adapter=graph.input("adapter","NuGEN Illumina-compatible library adapter",
                        MaterialKind.ADAPTER,inferred=True,
                        source=SOURCE+" (commercial kit adapter undisclosed)")
    adapted=graph.product("adapted","adapter-completed TT-seq cDNA",
                          MaterialKind.MOLECULAR,inferred=True,
                          source=SOURCE+" (commercial kit construct)")
    graph.react("attach_adapter","Install the commercial sequencing adapters",
                StepKind.ADAPTER_ATTACHMENT,
                ((cdna,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL)),(adapted,))
    fwd=graph.input("pcr_forward","NuGEN PCR forward primer",MaterialKind.OLIGO,
                    inferred=True,source=SOURCE+" (commercial kit oligo undisclosed)")
    rev=graph.input("pcr_reverse","NuGEN indexed PCR reverse primer",MaterialKind.OLIGO,
                    inferred=True,source=SOURCE+" (commercial kit oligo undisclosed)")
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY,inferred=True,
                          source=SOURCE+" (inferred Illumina-compatible endpoint)")
    graph.react("library_pcr","Amplify and index the captured cDNA library",StepKind.PCR,
                ((adapted,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
                 (rev,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The defining supplement names the NuGEN kit without printing its PCR oligos."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "cells","s4u","rna_polymerase","pulse_label","fragmentation_condition",
        "fragment","biotin_reagent","biotinylate","streptavidin","capture",
        "rt_primer","reverse_transcriptase","reverse_transcribe","adapter",
        "attach_adapter","pcr_forward","pcr_reverse","pcr_polymerase",
        "library_pcr","library"),SOURCE)
    return graph.require_publishable()
