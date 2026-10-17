"""SHAPE-MaP chemical probing and mutational profiling workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg, truseq_library
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Row
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_special import mutational_read_rows

TITLE = "SHAPE-MaP — RNA structure by mutational profiling"
NOTES = "01_shape-map.html"
SOURCE = 'Defining paper: <a href="https://doi.org/10.1038/nmeth.3029">Siegfried et al. (2014)</a>; implementation guidance: <a href="https://doi.org/10.1038/nprot.2018.044">Busan &amp; Weeks (2018)</a>.'
SUMMARY = "A SHAPE reagent acylates conformationally flexible RNA 2′-hydroxyls. MaP reverse transcription records those adducts as mutations in cDNA, which is converted to an indexed sequencing library."
CAVEAT = "The 2014 paper gives exact fixed primer arms for its targeted two-round PCR route, plus a separate fragmented-cDNA TruSeq ligation route. This schematic represents the latter endpoint; the source does not print its ordered adapter pair or the sample barcode identities."
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library([seg("fragment-derived MaP cDNA", "X"*44, placeholder=True)], "SHAPE-MaP library", inferred_adapters=True)
FINAL_CAPTION = "Indexed fragmented-cDNA library containing MaP-encoded substitutions at SHAPE-reactive positions."
SEQUENCING_INTRO = "Standard sequencing primers read the cDNA; structure is inferred from mutation frequencies rather than termination positions."

def sections():
    return [
        ("Modify flexible RNA nucleotides", [Row(chunks=[("folded RNA     paired stem ║ loop / flexible region ~~~~~", None, False)]), Row(chunks=[("SHAPE reagent                         ●  ●●   ●", "umi", False)])], "Electrophilic SHAPE reagent forms covalent 2′-O-adducts preferentially at flexible nucleotides; untreated and denatured controls accompany the modified sample."),
        ("Encode adducts as cDNA mutations", mutational_read_rows(), "MaP conditions encourage reverse transcriptase to traverse adducts while misincorporating, rather than stopping."),
        ("Build the fragmented-cDNA sequencing library", [Row(chunks=[("fragment → RT → second strand → end repair → dA tail → adapter ligation → PCR", None, False)]), Row(chunks=[("INFERRED — canonical flow-cell arms flank the assayed RNA-derived insert", None, True)])], "INFERRED — the selected forked-adapter bases remain unverified; this view follows the fragmented-cDNA library route."),
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
    fragment_condition=graph.input("fragment_condition","RNA fragmentation condition",MaterialKind.CONDITION)
    fragmented=graph.product("fragmented_rna","fragmented SHAPE-modified RNA",MaterialKind.MOLECULAR)
    graph.react("fragment","Fragment modified RNA for the library-ligation branch",StepKind.FRAGMENTATION,
                ((modified,InputRole.TEMPLATE),(fragment_condition,InputRole.CONDITION)),(fragmented,))
    rt_primer=graph.input("rt_primer","random nonamer RT primers (fragmented ribosome/intron branch)",MaterialKind.OLIGO,
                          inferred=True,source=SOURCE+" (target-specific identity varies)")
    rt=graph.input("map_rt","mutational-profiling reverse transcriptase",
                   MaterialKind.CATALYST)
    cdna=graph.product("map_cdna","mutation-encoded MaP cDNA",MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe","Traverse SHAPE adducts under MaP conditions",
                StepKind.REVERSE_TRANSCRIPTION,
                ((fragmented,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    second_mix=graph.input("second_mix","NEBNext second-strand synthesis module",MaterialKind.CATALYST)
    duplex=graph.product("duplex","double-stranded MaP cDNA",MaterialKind.MOLECULAR)
    graph.react("second_strand","Synthesize the second cDNA strand",StepKind.EXTENSION,((cdna,InputRole.TEMPLATE),(second_mix,InputRole.CATALYST)),(duplex,))
    repair_mix=graph.input("repair_mix","NEBNext end-repair module",MaterialKind.CATALYST)
    repaired=graph.product("repaired","end-repaired MaP cDNA",MaterialKind.MOLECULAR)
    graph.react("end_repair","Repair double-stranded cDNA ends",StepKind.END_REPAIR,((duplex,InputRole.TEMPLATE),(repair_mix,InputRole.CATALYST)),(repaired,))
    tail_mix=graph.input("tail_mix","NEBNext dA-tailing module",MaterialKind.CATALYST)
    tailed=graph.product("tailed","dA-tailed MaP cDNA",MaterialKind.MOLECULAR)
    graph.react("da_tail","Add ligation-ready dA tails",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(tail_mix,InputRole.CATALYST)),(tailed,))
    adapters=graph.input("truseq_adapters","TruSeq-compatible forked adapters",MaterialKind.ADAPTER,inferred=True)
    ligase=graph.input("ligation_mix","quick DNA ligation module",MaterialKind.CATALYST)
    adapted=graph.product("adapted","adapter-ligated MaP cDNA",MaterialKind.MOLECULAR,inferred=True)
    graph.react("ligate_adapters","Ligate forked adapters to fragmented-cDNA library",StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapters,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),
                geometry=AdapterLigationGeometry.unavailable("Selected ordered forked-adapter end geometry has not been verified in this model."))
    forward=graph.input("pcr_forward","fragmented-library PCR forward primer",
                        MaterialKind.OLIGO,inferred=True,
                        source=SOURCE+" (implementation-specific primer)")
    reverse=graph.input("pcr_reverse","indexed fragmented-library PCR reverse primer",
                        MaterialKind.OLIGO,inferred=True,
                        source=SOURCE+" (implementation-specific primer)")
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY,inferred=True,
                          source=SOURCE+" (modelled fragmented-cDNA endpoint)")
    graph.react("library_pcr","Amplify the fragmented-cDNA library (30-cycle Q5 emulsion PCR)",
                StepKind.PCR,
                ((adapted,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="Selected library PCR pair and binding geometry remain unassigned for this fragmented-cDNA branch."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "fragment_condition","fragment","second_mix","second_strand","repair_mix","end_repair","tail_mix","da_tail","truseq_adapters","ligation_mix","ligate_adapters",
        "rna","shape_reagent","modify","rt_primer","map_rt",
        "reverse_transcribe","pcr_forward","pcr_reverse","pcr_polymerase",
        "library_pcr","library"),SOURCE)
    return graph.require_publishable()
