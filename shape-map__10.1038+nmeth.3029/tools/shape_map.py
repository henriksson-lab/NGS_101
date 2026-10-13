"""SHAPE-MaP chemical probing and mutational profiling workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, truseq_library
from chemdraw import Row
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_special import mutational_read_rows

TITLE = "SHAPE-MaP — RNA structure by mutational profiling"
NOTES = "01_shape-map.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1038/nmeth.3029">Siegfried et al. (2014)</a>; implementation guidance: <a href="https://doi.org/10.1038/nprot.2018.044">Busan &amp; Weeks (2018)</a>.'
SUMMARY = "A SHAPE reagent acylates conformationally flexible RNA 2′-hydroxyls. MaP reverse transcription records those adducts as mutations in cDNA, which is converted to an indexed sequencing library."
CAVEAT = "SHAPE-MaP was selected over DMS-MaPseq because it is the defining general-purpose MaP workflow and has an authoritative full protocol. The chemical adduct is a mark, not an added sequencing adapter."
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([seg("MaP cDNA amplicon", "X"*44, placeholder=True)], "SHAPE-MaP library", inferred_adapters=True)
FINAL_CAPTION = "Indexed amplicon containing MaP-encoded substitutions at SHAPE-reactive positions."
SEQUENCING_INTRO = "Standard sequencing primers read the cDNA; structure is inferred from mutation frequencies rather than termination positions."

def sections():
    return [
        ("Modify flexible RNA nucleotides", [Row(chunks=[("folded RNA     paired stem ║ loop / flexible region ~~~~~", None, False)]), Row(chunks=[("SHAPE reagent                         ●  ●●   ●", "umi", False)])], "Electrophilic SHAPE reagent forms covalent 2′-O-adducts preferentially at flexible nucleotides; untreated and denatured controls accompany the modified sample."),
        ("Encode adducts as cDNA mutations", mutational_read_rows(), "MaP conditions encourage reverse transcriptase to traverse adducts while misincorporating, rather than stopping."),
        ("Amplify target and install sequencing arms", [Row(chunks=[("gene-specific PCR / library PCR → indexed double-stranded cDNA", None, False)]), Row(chunks=[("INFERRED — canonical flow-cell arms flank the assayed RNA-derived insert", None, True)])], "INFERRED — library-arm bases depend on the chosen amplicon implementation; the panel below states the modeled standard configuration."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("rna","folded RNA",MaterialKind.MOLECULAR,primary=True)
    shape=graph.input("shape_reagent","electrophilic SHAPE modification condition",
                      MaterialKind.CONDITION)
    modified=graph.product("modified","2-prime-O-adduct-modified RNA",
                           MaterialKind.MOLECULAR)
    graph.react("modify","Acylate flexible RNA 2-prime hydroxyls",
                StepKind.CHEMICAL_CONVERSION,
                ((rna,InputRole.TEMPLATE),(shape,InputRole.CONDITION)),(modified,))
    rt_primer=graph.input("rt_primer","target-specific MaP RT primer",MaterialKind.OLIGO,
                          inferred=True,source=SOURCE+" (target-specific identity varies)")
    rt=graph.input("map_rt","mutational-profiling reverse transcriptase",
                   MaterialKind.CATALYST)
    cdna=graph.product("map_cdna","mutation-encoded MaP cDNA",MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe","Traverse SHAPE adducts under MaP conditions",
                StepKind.REVERSE_TRANSCRIPTION,
                ((modified,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    forward=graph.input("pcr_forward","target/library PCR forward primer",
                        MaterialKind.OLIGO,inferred=True,
                        source=SOURCE+" (implementation-specific primer)")
    reverse=graph.input("pcr_reverse","indexed target/library PCR reverse primer",
                        MaterialKind.OLIGO,inferred=True,
                        source=SOURCE+" (implementation-specific primer)")
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY,inferred=True,
                          source=SOURCE+" (modelled standard amplicon endpoint)")
    graph.react("library_pcr","Amplify the target and install sequencing arms",
                StepKind.PCR,
                ((cdna,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="SHAPE-MaP target primers and library-arm implementation vary by assayed RNA."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "rna","shape_reagent","modify","rt_primer","map_rt",
        "reverse_transcribe","pcr_forward","pcr_reverse","pcr_polymerase",
        "library_pcr","library"),SOURCE)
    return graph.require_publishable()
