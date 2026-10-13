"""Zubradt et al. genome-wide DMS-MaPseq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Scene, Segment, circle_rows, revcomp
from circular import circularize_ssdna
from circular import CircularizationGeometry
from rna_special import adapter_ligation_scene, mutational_read_rows
from adapter_ligation import reported_adapter_end
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "DMS-MaPseq — mutations report RNA structure"
NOTES = "01_dms-mapseq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.4057">Zubradt et al., <i>Nature Methods</i> (2017)</a>.'
SUMMARY = "DMS modifies exposed A and C bases in living cells; TGIRT reads through those lesions and records them as substitutions before cDNA circularization and custom-primer Illumina sequencing."
CAVEAT = "This page follows the paper's genome-wide branch. Its targeted RT-PCR and UMI-targeted branches are related protocols with different library conversion."
LINKER = "CACTCGGGCACCAAGGA"
RT_HANDLE = "GATCGTCGGACTGTAGAACTCTGAACCTGTCG"
READ1 = "GCAGCGACAGGTTCAGAGTTCTACAGTCCGACGATC"
CDNA = Construct([Segment("RT-primer handle", RT_HANDLE, "r1"),
                  Segment("mutation-bearing cDNA", "X"*42, placeholder=True)],
                 name="DMS-MaPseq cDNA")
CIRCLE = circularize_ssdna(CDNA, five_prime_phosphate=True)
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = "Custom Read 1 primer oNTI202 anneals to the PCR-added handle and reads into the mutation-bearing cDNA; the paper used 50-nt single-end HiSeq 4000 reads."

def final_primer_scene():
    template = [Segment("oNTI202 binding site", revcomp(READ1), "r1"),
                Segment("mutation-bearing cDNA", "X"*42, placeholder=True)]
    primer = [Segment("oNTI202", READ1, "r1")]
    sc=Scene(); sc.strand("library",template,label="PCR-completed library")
    sc.anneal("Read 1 primer",primer,to="library",pair=("oNTI202","oNTI202 binding site"),label="oNTI202")
    sc.arrow("Read 1 primer","Read 1 enters cDNA")
    sc.labels("library")
    return sc

def sections():
    return [
        ("Modify exposed A and C bases with DMS", mutational_read_rows(),
         "DMS reacts preferentially at unpaired Watson-Crick faces; the genome-wide treatment averages about one modification per 50 nt."),
        ("Fragment RNA, dephosphorylate and ligate the 3-prime linker", adapter_ligation_scene(adapter_name="5-prime-rApp linker-2").rows(),
         "Zn2+ fragmentation is followed by size selection and rSAP treatment; truncated T4 RNA ligase 2 joins the preadenylated, 3-prime-blocked linker."),
        ("Reverse-transcribe through lesions with TGIRT", [*mutational_read_rows(), Row(chunks=[("TGIRT cDNA carries substitutions at DMS adducts", "r1", False)])],
         "The long 57-degree reaction favors processive read-through so multiple modifications can be encoded on one cDNA."),
        ("Remove RNA and circularize the phosphorylated cDNA", circle_rows(CDNA,"CircLigase closure"),
         "The RT primer supplies the 5-prime phosphate required for CircLigase; the cDNA becomes a single-stranded circle."),
        ("PCR-add Illumina arms and bind the custom sequencing primer", final_primer_scene().rows(),
         "Nine to thirteen PCR cycles add the indexed flow-cell arms. oNTI202—not a stock TruSeq Read 1 primer—starts the single-end read."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("rna","cellular RNA in vivo",MaterialKind.MOLECULAR,primary=True)
    dms=graph.input("dms","dimethyl sulfate modification condition",
                    MaterialKind.CONDITION)
    modified=graph.product("modified","DMS-modified cellular RNA",MaterialKind.MOLECULAR)
    graph.react("modify","Modify exposed RNA A and C bases with DMS",
                StepKind.CHEMICAL_CONVERSION,
                ((rna,InputRole.TEMPLATE),(dms,InputRole.CONDITION)),(modified,))
    zinc=graph.input("zinc","zinc fragmentation condition",MaterialKind.CONDITION)
    fragmented=graph.product("fragmented","zinc-fragmented DMS-modified RNA",
                             MaterialKind.MOLECULAR)
    graph.react("fragment","Fragment modified RNA",StepKind.FRAGMENTATION,
                ((modified,InputRole.TEMPLATE),(zinc,InputRole.CONDITION)),(fragmented,))
    size_rule=graph.input("size_rule","RNA size-selection rule",MaterialKind.SELECTION)
    selected=graph.product("selected","size-selected RNA fragments",MaterialKind.MOLECULAR)
    offsize=graph.product("offsize","off-size RNA fragments",MaterialKind.WASTE)
    graph.react("size_select","Select the intended RNA fragment sizes",StepKind.SELECTION,
                ((fragmented,InputRole.TEMPLATE),(size_rule,InputRole.SELECTION)),
                (selected,offsize))
    graph.discard(offsize,"off-size fragments are excluded")
    phosphatase=graph.input("phosphatase","rSAP",MaterialKind.CATALYST)
    dephosphorylated=graph.product("dephosphorylated","3-prime-dephosphorylated RNA fragments",
                                  MaterialKind.MOLECULAR)
    graph.react("dephosphorylate","Prepare RNA 3-prime ends for linker ligation",
                StepKind.DEPHOSPHORYLATION,
                ((selected,InputRole.TEMPLATE),(phosphatase,InputRole.CATALYST)),
                (dephosphorylated,))
    linker=graph.input("linker","5-prime-rApp linker-2 with blocked 3-prime end",
                       MaterialKind.ADAPTER)
    ligase=graph.input("ligase","truncated T4 RNA ligase 2",MaterialKind.CATALYST)
    linked=graph.product("linked","3-prime-linker-ligated modified RNA",
                         MaterialKind.MOLECULAR)
    geometry=reported_adapter_end(
        "linker-2",kind="single-strand",polarity="5-prime",terminal_sequence=LINKER,
        statement="The defining protocol prints preadenylated linker-2 and ligates it to dephosphorylated RNA 3-prime ends.")
    graph.react("ligate_linker","Ligate preadenylated linker-2",StepKind.LIGATION,
                ((dephosphorylated,InputRole.TEMPLATE),(linker,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(linked,),geometry=geometry)
    rt_primer=graph.input("rt_primer","phosphorylated TGIRT reverse-transcription primer",
                          MaterialKind.OLIGO)
    tgirt=graph.input("tgirt","TGIRT reverse transcriptase",MaterialKind.CATALYST)
    cdna=graph.product("cdna","mutation-bearing phosphorylated cDNA",
                       MaterialKind.MOLECULAR,state=CDNA)
    graph.react("reverse_transcribe","Read through DMS lesions with TGIRT",
                StepKind.REVERSE_TRANSCRIPTION,
                ((linked,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (tgirt,InputRole.CATALYST)),(cdna,))
    rnase=graph.input("rnase","RNA-removal enzyme",MaterialKind.CATALYST)
    clean_cdna=graph.product("clean_cdna","RNA-free phosphorylated cDNA",
                             MaterialKind.MOLECULAR,state=CDNA)
    rna_waste=graph.product("rna_waste","digested RNA template",MaterialKind.WASTE)
    graph.react("remove_rna","Remove the RNA template",StepKind.DIGESTION,
                ((cdna,InputRole.TEMPLATE),(rnase,InputRole.CATALYST)),
                (clean_cdna,rna_waste))
    graph.discard(rna_waste,"the copied RNA template is enzymatically removed")
    circligase=graph.input("circligase","CircLigase",MaterialKind.CATALYST)
    circle=graph.product("circle","circular mutation-bearing cDNA",MaterialKind.MOLECULAR,
                         state=CIRCLE)
    graph.react("circularize","Circularize phosphorylated cDNA",
                StepKind.CIRCULARIZATION,
                ((clean_cdna,InputRole.TEMPLATE),(circligase,InputRole.CATALYST)),
                (circle,),geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="CircLigase closes the RT-primer-supplied 5-prime phosphate to the cDNA 3-prime hydroxyl."))
    forward=graph.input("pcr_forward","circle PCR forward primer",MaterialKind.OLIGO)
    reverse=graph.input("pcr_reverse","indexed circle PCR reverse primer",MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library","PCR-completed DMS-MaPseq library",
                          MaterialKind.LIBRARY,state=final_primer_scene())
    graph.react("library_pcr","Add indexed Illumina arms",StepKind.PCR,
                ((circle,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    CDNA,reason="The complete indexed outer PCR oligos are not reconstructed in this compact module."))
    graph.sequence(library,SequencingPlan(
        limitation="Custom primer oNTI202 binds the modelled PCR-added site and reads mutation-bearing cDNA; the complete indexed outer-arm construct is not available."))
    graph.declare_reviewed_scope((
        "rna","dms","modify","zinc","fragment","size_rule","size_select",
        "phosphatase","dephosphorylate","linker","ligase","ligate_linker",
        "rt_primer","tgirt","reverse_transcribe","rnase","remove_rna",
        "circligase","circularize","pcr_forward","pcr_reverse","pcr_polymerase",
        "library_pcr","library"),SOURCE)
    return graph.require_publishable()
