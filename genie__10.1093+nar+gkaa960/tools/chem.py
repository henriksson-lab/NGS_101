"""GenIE targeted cDNA readout."""
from batch_ngs import seg, truseq_library
from chemdraw import Row
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "GenIE — endogenous variant editing with targeted cDNA sequencing"
NOTES = "01_genie.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1093/nar/gkaa960">Cooper et al., <i>Nucleic Acids Research</i> (2020)</a>.'
SUMMARY = "Arrayed CRISPR editing creates endogenous alleles; locus-specific RT-PCR measures each allele's effect on transcript abundance or splicing."
FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(
    [seg("locus-specific cDNA amplicon", "X" * 56, placeholder=True)],
    "GenIE targeted cDNA library")
FINAL_CAPTION = "The insert is the locus-specific cDNA amplicon spanning the edited variant and the transcript feature being measured."
SEQUENCING_INTRO = "Standard paired-end Illumina primers read the indexed targeted amplicon."

def sections():
    return [
        ("Create an allelic series at the endogenous locus", [Row(chunks=[("arrayed guide + HDR donor → edited cell population", "cbc", False)])], "Each well contains one designed edit or control."),
        ("Reverse-transcribe cellular RNA", [Row(chunks=[("RNA → locus-containing cDNA", "r1", False)])], "The molecular readout is RNA, rather than genomic editing frequency alone."),
        ("Amplify the locus from cDNA", [Row(chunks=[("outer locus PCR → nested allele-spanning amplicon", "r2", False)])], "Primer pairs are locus-specific and are tabulated in the paper's supplement."),
        ("Install indexes and sequence", [Row(chunks=[("P5—index—read arm—target cDNA—read arm—index—P7", "cbc", False)])], "Indexed amplicons from the arrayed experiment are pooled for Illumina sequencing."),
    ]

def protocol_graph():
    source="Defining source: https://doi.org/10.1093/nar/gkaa960"
    graph=ProtocolGraph("GenIE",authority=source)
    cells=graph.input("cells","arrayed cell population",MaterialKind.BIOLOGICAL,primary=True)
    guide=graph.input("guide","arrayed CRISPR guide RNA",MaterialKind.OLIGO)
    donor=graph.input("hdr_donor","designed HDR donor",MaterialKind.OLIGO)
    editor=graph.input("editor","CRISPR genome-editing nuclease",MaterialKind.CATALYST)
    edited=graph.product("edited_cells","endogenous allelic-series cell population",
                         MaterialKind.BIOLOGICAL)
    graph.react("edit","Create the endogenous allelic series",StepKind.GENOME_EDITING,
                ((cells,InputRole.TEMPLATE),(guide,InputRole.STRUCTURAL),
                 (donor,InputRole.STRUCTURAL),(editor,InputRole.CATALYST)),(edited,))
    extraction=graph.input("rna_extraction","RNA extraction conditions",MaterialKind.CONDITION)
    rna=graph.product("rna","RNA from the edited cell population",MaterialKind.MOLECULAR)
    graph.react("extract_rna","Extract cellular RNA",StepKind.EXTRACTION,
                ((edited,InputRole.TEMPLATE),(extraction,InputRole.CONDITION)),(rna,))
    rt_primer=graph.input("rt_primer","locus-compatible reverse-transcription primer",
                          MaterialKind.OLIGO,inferred=True,
                          source=source+" (experiment-specific primer choice)")
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    cdna=graph.product("cdna","locus-containing cDNA",MaterialKind.MOLECULAR)
    graph.react("reverse_transcription","Reverse-transcribe cellular RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    outer_f=graph.input("outer_forward","outer locus PCR forward primer",MaterialKind.OLIGO)
    outer_r=graph.input("outer_reverse","outer locus PCR reverse primer",MaterialKind.OLIGO)
    polymerase=graph.input("polymerase","amplicon PCR polymerase",MaterialKind.CATALYST)
    outer=graph.product("outer_amplicon","outer locus cDNA amplicon",MaterialKind.MOLECULAR)
    graph.react("outer_pcr","Amplify the edited locus from cDNA",StepKind.PCR,
                ((cdna,InputRole.TEMPLATE),(outer_f,InputRole.STRUCTURAL),
                 (outer_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(outer,),
                geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                    reason="Primer pairs are locus-specific; the final page represents their amplicon generically."))
    nested_f=graph.input("nested_forward","nested allele-spanning forward primer",MaterialKind.OLIGO)
    nested_r=graph.input("nested_reverse","nested allele-spanning reverse primer",MaterialKind.OLIGO)
    nested=graph.product("nested_amplicon","nested allele-spanning cDNA amplicon",MaterialKind.MOLECULAR)
    graph.react("nested_pcr","Generate the nested allele-spanning amplicon",StepKind.PCR,
                ((outer,InputRole.TEMPLATE),(nested_f,InputRole.STRUCTURAL),
                 (nested_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(nested,),
                geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                    reason="The represented insert is a generic locus-specific nested amplicon."))
    i5=graph.input("index_forward","indexed Illumina forward primer",MaterialKind.OLIGO)
    i7=graph.input("index_reverse","indexed Illumina reverse primer",MaterialKind.OLIGO)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("index_pcr","Install indexes and flow-cell arms",StepKind.PCR,
                ((nested,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                    reason="The final standard Illumina shell is exact; the locus-specific pre-index amplicon is generic."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
    graph.declare_reviewed_scope((
        "cells","guide","hdr_donor","editor","edit","rna_extraction","extract_rna",
        "rt_primer","reverse_transcriptase","reverse_transcription","outer_forward",
        "outer_reverse","polymerase","outer_pcr","nested_forward","nested_reverse",
        "nested_pcr","index_forward","index_reverse","index_pcr","library"),source)
    return graph.require_publishable()
