"""scGRO-seq clickable nascent-RNA chemistry."""
from chemdraw import Row
from multibranch import illumina_branch
from protocol_extensions import clicked_nascent_rna
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
TITLE="scGRO-seq — single-cell nascent RNA by click chemistry"
NOTES="01_scgro-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41586-024-07517-7">Mahat et al., <i>Nature</i> (2024)</a>.'
SUMMARY="Nuclear run-on incorporates chain-terminating O-propargyl nucleotides; click chemistry attaches a barcode/UMI oligo. Fragmentation, reverse transcription, pre-amplification, DASH depletion and final indexed PCR recover the library."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("scGRO-seq","nascent-RNA cDNA",cell=10,umi=5,
                                          inferred_adapters=True)
FINAL_CAPTION="INFERRED — the 10-nt cell barcode and 5-nt UMI represent the third 43-nt azide adaptor in Supplementary Table 7; the complete final adapter junction is not strand-resolved here."
SEQUENCING_INTRO="The paper reports a 64/43/8/8 NovaSeq run and trims Read 1 motif CCCCTGTCTCTTATACACAT and Read 2 motif AGATCGGAAGAGCGTCGTGT before extracting cell barcodes and UMIs from Read 1. The displayed primer binding remains an inferred final-library model."
CAVEAT=("Supplementary Table 7 gives three alternative 5′-azide adaptor designs with "
        "different barcode/UMI dimensions (10/6, 12/4 and 10/5). The final-library "
        "drawing uses the third design as a representative; its Illumina shell remains inferred.")
READ_LENGTHS={"Read 1":64,"Read 2":43,"Index 1 (i7)":8,"Index 2 (i5)":8}

def oligos():
    return [
        "<b>5′-AzScBc design 1, 43 nt</b> — 5′-[azideT]JJJJJJJJJJNNNNNNGATCGTCGGACTGTAGAACTCTGAAC-3′ (J = cell barcode; N = UMI)",
        "<b>5′-AzScBc design 2, 86 nt hairpin</b> — 5′-[azideT]NNNNJJJJJJJJJJJJGATCGTCGGACTGTAGAACTCTGCGGCCGTGCTCGTTTTCGAGCACGGCCGCAGAGTTCTACAGTCCGA-3′",
        "<b>5′-AzScBc design 3, 43 nt</b> — 5′-[azideT]ACAGGNNNNNJJJJJJJJJJAGATCGGAAGAGCGTCGTGTAG[SpC3]-3′",
    ]
def sections():
    return [("Run on engaged polymerases with 3-prime O-propargyl NTPs",[Row(chunks=[("nascent RNA → RNA—O-propargyl-NMP (chain terminated)", "w1", False)])],"The terminal alkyne marks the active RNA 3-prime end."),("Click the azide barcode oligo to the alkyne",clicked_nascent_rna().rows(),"Copper-catalysed azide–alkyne cycloaddition creates a triazole linkage that RT can traverse."),("Reverse-transcribe and add sequencing adapters",[Row(chunks=[("clicked RNA → barcoded cDNA → PCR/index adapters", "r1", False)])],"The clicked oligo supplies molecule and cell identity before amplification."),]

def protocol_graph():
    source="Defining source: https://doi.org/10.1038/s41586-024-07517-7"
    graph=ProtocolGraph("scGRO-seq",authority=source)
    complexes=graph.input("engaged_complexes","engaged nuclear transcription complexes",
                          MaterialKind.BIOLOGICAL,primary=True)
    nucleotide=graph.input("propargyl_nucleotide","3′ O-propargyl run-on nucleotide",
                           MaterialKind.MOLECULAR)
    polymerase=graph.input("engaged_polymerase","engaged RNA polymerase",MaterialKind.CATALYST)
    labelled=graph.product("labelled_rna","alkyne-terminated nascent RNA",MaterialKind.MOLECULAR)
    graph.react("run_on","Run on engaged polymerases with 3′ O-propargyl NTPs",StepKind.RUN_ON,
                ((complexes,InputRole.TEMPLATE),(nucleotide,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),(labelled,))
    click_oligo=graph.input("click_oligo","5′-azide barcode/UMI adaptor (design 3)",
                            MaterialKind.ADAPTER,
                            source=source+" (Supplementary Table 7, design 3)")
    clicked=graph.product("clicked_rna","barcoded clicked nascent RNA",MaterialKind.MOLECULAR,
                          state=clicked_nascent_rna(barcode_length=10,umi_length=5))
    unclicked=graph.product("unclicked_azide","unreacted azide-barcode oligo",MaterialKind.MOLECULAR)
    graph.react("click_attachment","Click the azide barcode oligo to the alkyne",
                StepKind.ADAPTER_ATTACHMENT,
                ((labelled,InputRole.TEMPLATE),(click_oligo,InputRole.STRUCTURAL)),(clicked,unclicked))
    zinc=graph.input("zinc_chloride","ZnCl2 fragmentation condition",MaterialKind.CONDITION)
    fragments=graph.product("fragmented_clicked_rna","fragmented clicked RNA–DNA conjugates",MaterialKind.MOLECULAR)
    graph.react("fragment_clicked_rna","Fragment clicked RNA with ZnCl2",StepKind.FRAGMENTATION,
                ((clicked,InputRole.TEMPLATE),(zinc,InputRole.CONDITION)),(fragments,))
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    tso=graph.input("template_switch_oligo","template-switching oligo",MaterialKind.OLIGO,
                    inferred=True,source=source+" (ordering sequence unresolved)")
    cdna=graph.product("cdna","barcoded nascent-RNA cDNA",MaterialKind.MOLECULAR)
    empty_cdna=graph.product("empty_cdna","cDNA from unclicked azide adaptor and TSO",MaterialKind.MOLECULAR)
    graph.react("reverse_transcription","Reverse-transcribe through the clicked junction",
                StepKind.REVERSE_TRANSCRIPTION,
                ((fragments,InputRole.TEMPLATE),(unclicked,InputRole.TEMPLATE),
                 (tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,empty_cdna))
    p5=graph.input("p5_primer","P5 library PCR primer",MaterialKind.OLIGO,
                   inferred=True,source=source+" (specific ordering primer unresolved)")
    p7=graph.input("p7_primer","indexed P7 library PCR primer",MaterialKind.OLIGO,
                   inferred=True,source=source+" (specific index primer unresolved)")
    pcr=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    pre_library=graph.product("pre_library","pre-amplified barcoded cDNA",MaterialKind.LIBRARY)
    empty_amplicon=graph.product("empty_amplicon","pre-amplified empty-adaptor artifact",MaterialKind.LIBRARY)
    graph.react("pre_amplify","Pre-amplify desired cDNA and empty-adaptor artifacts (6 cycles)",StepKind.PCR,
                ((cdna,InputRole.TEMPLATE),(empty_cdna,InputRole.TEMPLATE),
                 (p5,InputRole.STRUCTURAL),(p7,InputRole.STRUCTURAL),
                 (pcr,InputRole.CATALYST)),(pre_library,empty_amplicon),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The pre-amplified scGRO-seq amplicon and empty-adaptor junctions are not strand-resolved."))
    cas9=graph.input("cas9","SpCas9",MaterialKind.CATALYST)
    guide=graph.input("dash_guide","guide RNA targeting the ACAGG-dependent empty-adaptor amplicon",
                      MaterialKind.OLIGO,inferred=True,source=source+" (guide ordering sequence unresolved)")
    dash_library=graph.product("dash_library","barcoded cDNA surviving DASH",MaterialKind.LIBRARY)
    cut_empty=graph.product("cut_empty","Cas9-cleaved empty-adaptor amplicons",MaterialKind.WASTE)
    graph.react("dash","DASH: cut empty-adaptor amplicons with Cas9–guide complex",StepKind.DIGESTION,
                ((pre_library,InputRole.TEMPLATE),(empty_amplicon,InputRole.TEMPLATE),
                 (guide,InputRole.STRUCTURAL),(cas9,InputRole.CATALYST)),(dash_library,cut_empty))
    graph.discard(cut_empty,"Cas9 cleavage prevents amplification of the empty-adaptor artifact")
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr","Final indexed PCR of DASH-surviving molecules",StepKind.PCR,
                ((dash_library,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(pcr,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The compact branch model records the final construct but not separate PCR-primer landing objects."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
    graph.declare_reviewed_scope((
        "engaged_complexes","propargyl_nucleotide","engaged_polymerase","run_on",
        "click_oligo","click_attachment","zinc_chloride","fragment_clicked_rna",
        "reverse_transcriptase","template_switch_oligo","reverse_transcription",
        "p5_primer","p7_primer","pcr_polymerase","pre_amplify","cas9","dash_guide",
        "dash","library_pcr","library"),source)
    return graph.require_publishable()
