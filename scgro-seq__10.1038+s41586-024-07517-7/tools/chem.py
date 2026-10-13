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
SUMMARY="Nuclear run-on incorporates a chain-terminating O-propargyl nucleotide; click ligation attaches the barcode/UMI/RT handle used to recover nascent RNA from single cells."
FINAL_LIBRARY,SEQ_PRIMERS=illumina_branch("scGRO-seq","nascent-RNA cDNA",cell=8,umi=8)
FINAL_CAPTION="Cell barcode and UMI precede cDNA copied from the run-on-labelled RNA 3-prime end."
SEQUENCING_INTRO="Declared Illumina primers bind the completed scGRO-seq library."
def sections():
    return [("Run on engaged polymerases with 3-prime O-propargyl NTPs",[Row(chunks=[("nascent RNA → RNA—O-propargyl-NMP (chain terminated)", "w1", False)])],"The terminal alkyne marks the active RNA 3-prime end."),("Click the azide barcode oligo to the alkyne",clicked_nascent_rna().rows(),"Copper-catalysed azide–alkyne cycloaddition creates a triazole linkage that RT can traverse."),("Reverse-transcribe and add sequencing adapters",[Row(chunks=[("clicked RNA → barcoded cDNA → PCR/index adapters", "r1", False)])],"The clicked oligo supplies molecule and cell identity before amplification."),]

def protocol_graph():
    source="Defining source: https://doi.org/10.1038/s41586-024-07517-7"
    graph=ProtocolGraph("scGRO-seq",authority=source)
    complexes=graph.input("engaged_complexes","engaged nuclear transcription complexes",
                          MaterialKind.BIOLOGICAL,primary=True)
    nucleotide=graph.input("propargyl_nucleotide","3′ O-propargyl run-on nucleotide",
                           MaterialKind.OLIGO)
    polymerase=graph.input("engaged_polymerase","engaged RNA polymerase",MaterialKind.CATALYST)
    labelled=graph.product("labelled_rna","alkyne-terminated nascent RNA",MaterialKind.MOLECULAR)
    graph.react("run_on","Run on engaged polymerases with 3′ O-propargyl NTPs",StepKind.RUN_ON,
                ((complexes,InputRole.TEMPLATE),(nucleotide,InputRole.STRUCTURAL),
                 (polymerase,InputRole.CATALYST)),(labelled,))
    click_oligo=graph.input("click_oligo","azide barcode/UMI/RT oligo",MaterialKind.ADAPTER,
                            state=clicked_nascent_rna())
    clicked=graph.product("clicked_rna","barcoded clicked nascent RNA",MaterialKind.MOLECULAR,
                          state=clicked_nascent_rna())
    graph.react("click_attachment","Click the azide barcode oligo to the alkyne",
                StepKind.ADAPTER_ATTACHMENT,
                ((labelled,InputRole.TEMPLATE),(click_oligo,InputRole.STRUCTURAL)),(clicked,))
    rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    cdna=graph.product("cdna","barcoded nascent-RNA cDNA",MaterialKind.MOLECULAR)
    graph.react("reverse_transcription","Reverse-transcribe through the clicked junction",
                StepKind.REVERSE_TRANSCRIPTION,
                ((clicked,InputRole.TEMPLATE),(click_oligo,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(cdna,))
    p5=graph.input("p5_primer","P5 library PCR primer",MaterialKind.OLIGO)
    p7=graph.input("p7_primer","indexed P7 library PCR primer",MaterialKind.OLIGO)
    pcr=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("library_pcr","Amplify and add the Illumina sequencing shell",StepKind.PCR,
                ((cdna,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(pcr,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The compact branch model records the final construct but not separate PCR-primer landing objects."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
        read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
    graph.declare_reviewed_scope((
        "engaged_complexes","propargyl_nucleotide","engaged_polymerase","run_on",
        "click_oligo","click_attachment","reverse_transcriptase","reverse_transcription",
        "p5_primer","p7_primer","pcr_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
