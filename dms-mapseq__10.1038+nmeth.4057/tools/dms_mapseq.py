"""Zubradt et al. genome-wide DMS-MaPseq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import (Construct, Row, Scene, Segment, circle_rows, feature,
                      oligo, revcomp)
from circular import circularize_ssdna
from circular import CircularizationGeometry
import illumina as il
from rna_special import adapter_ligation_scene, mutational_read_rows
from adapter_ligation import reported_adapter_end
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
import seqprimers as sp

TITLE = "DMS-MaPseq — mutations report RNA structure"
NOTES = "01_dms-mapseq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.4057">Zubradt et al., <i>Nature Methods</i> (2017)</a>.'
SUMMARY = "DMS modifies exposed A and C bases in living cells; TGIRT reads through those lesions and records them as substitutions before cDNA circularization and custom-primer Illumina sequencing."
CAVEAT = "This page follows the paper's genome-wide branch. Its targeted RT-PCR and UMI-targeted branches are related protocols with different library conversion."
LINKER = "CACTCGGGCACCAAGGA"
RT_HANDLE = "GATCGTCGGACTGTAGAACTCTGAACCTGTCG"
RT_PRIMER = (RT_HANDLE + il.P7 + revcomp(LINKER))
PCR_FORWARD = "CAAGCAGAAGACGGCATACGA"
INDEX_PCR_PREFIX = il.P5 + il.INDEX1_PRIMER
INDEX_PCR_ANNEAL = "CGACAGGTTCAGAGTTC"
INDEX_PCR = INDEX_PCR_PREFIX + "N" * 6 + INDEX_PCR_ANNEAL
READ1 = revcomp(RT_HANDLE)
CDNA = Construct([Segment("RT-primer handle", RT_HANDLE, "r1"),
                  Segment("P7 in RT primer", il.P7, "p7"),
                  Segment("linker-2 complement", revcomp(LINKER), "r2"),
                  Segment("mutation-bearing cDNA", "X"*42, placeholder=True)],
                 name="DMS-MaPseq cDNA")
CIRCLE = circularize_ssdna(CDNA, five_prime_phosphate=True)
PCR_TEMPLATE = Construct([
    Segment("P7 in RT primer", il.P7, "p7"),
    Segment("linker-2 complement", revcomp(LINKER), "r2"),
    Segment("mutation-bearing cDNA", "X"*42, placeholder=True),
    Segment("RT-primer handle", RT_HANDLE, "r1"),
], name="cDNA circle linearized between RT handle and P7")
FINAL_LIBRARY = Construct([
    Segment("P7", il.P7, "p7"),
    Segment("linker-2 complement", revcomp(LINKER), "r2"),
    Segment("mutation-bearing cDNA", "X"*42, placeholder=True),
    Segment("oNTI202 binding site", RT_HANDLE, "r1"),
    Segment("sample index", "N"*6, "i7", placeholder=True,
            feature=feature("dms_mapseq_sample_index", "sample_index", "unknown",
                            note="six-base TruSeq index selected in the PCR primer")),
    Segment("Index 1 primer site", revcomp(il.INDEX1_PRIMER), "r2"),
    Segment("P5 reverse complement", il.P5_RC, "p5"),
], name="PCR-completed DMS-MaPseq library")
SEQ_PRIMERS = (
    sp.custom("Read 1", "oNTI202", READ1, "Zubradt et al. Supplementary Table 1"),
    sp.TRUSEQ["I1"],
)
ACQUIRED_ROLES = ("Read 1", "Index 1 (i7)")
READ_LENGTHS = {"Read 1": 50, "Index 1 (i7)": 6}
READ_USES = {"Read 1": "include", "Index 1 (i7)": "demultiplex"}
SEQUENCING_INTRO = "The source-exact custom oNTI202 primer reads into mutation-bearing cDNA; the standard Index 1 primer reads the six-base sample index introduced by the indexing PCR primer."
FINAL_CAPTION = "Source-exact PCR endpoint, linearized from P7 to P5. The insert is schematic; every oligo-derived region is sequence exact."


def oligos():
    return (
        oligo("Linker 2", (Segment("linker-2", LINKER),),
              mods="/5rApp/", three="/3ddC/-3′"),
        oligo("oCJ200-link2 RT primer",
              (Segment("custom handle", RT_HANDLE),
               Segment("internal spacer", "/iSp18/", placeholder=True),
               Segment("P7", il.P7),
               Segment("linker-2 annealing region", revcomp(LINKER))),
              mods="/5Phos/", three="-3′"),
        oligo("oNTI231 library PCR primer", (Segment("oNTI231", PCR_FORWARD),)),
        oligo("indexing PCR primer",
              (Segment("P5", il.P5), Segment("Index 1 primer site", il.INDEX1_PRIMER),
               Segment("six-base TruSeq index", "N"*6, placeholder=True,
                       feature=feature("dms_mapseq_sample_index", "sample_index", "unknown")),
               Segment("circle annealing region", INDEX_PCR_ANNEAL))),
        oligo("oNTI202 custom Read 1 primer", (Segment("oNTI202", READ1),)),
    )

def final_primer_scene():
    primer = [Segment("oNTI202", READ1, "r1")]
    sc=Scene(); sc.strand("library",list(FINAL_LIBRARY),label="PCR-completed library")
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
    deadenylase=graph.input("deadenylase","5-prime deadenylase",MaterialKind.CATALYST)
    recj=graph.input("recj","RecJ exonuclease",MaterialKind.CATALYST)
    cleaned=graph.product("linker_depleted","RNA after excess-linker digestion",MaterialKind.MOLECULAR)
    linker_fragments=graph.product("linker_fragments","digested excess linker",MaterialKind.MOLECULAR)
    graph.react("digest_linker","Digest excess linker with deadenylase/RecJ (30°C, 1 h)",StepKind.DIGESTION,
                ((linked,InputRole.TEMPLATE),(deadenylase,InputRole.CATALYST),(recj,InputRole.CATALYST)),(cleaned,linker_fragments))
    graph.mark_terminal(Terminal(linker_fragments.id,Disposition.UNRESOLVED,reason="Recovery of each digested linker fragment is not represented."))
    ribozero=graph.input("ribozero","RiboZero rRNA-depletion reagent",MaterialKind.SELECTION)
    depleted=graph.product("rrna_depleted","rRNA-depleted linker-ligated RNA",MaterialKind.MOLECULAR)
    rrna=graph.product("rrna_fraction","rRNA-depletion fraction",MaterialKind.MOLECULAR)
    graph.react("deplete_rrna","Deplete rRNA with RiboZero before RT",StepKind.SELECTION,
                ((cleaned,InputRole.TEMPLATE),(ribozero,InputRole.SELECTION)),(depleted,rrna))
    graph.mark_terminal(Terminal(rrna.id,Disposition.UNRESOLVED,reason="Exact capture and recovery fractions are not modeled."))
    rt_primer=graph.input("rt_primer","phosphorylated TGIRT reverse-transcription primer",
                          MaterialKind.OLIGO)
    tgirt=graph.input("tgirt","TGIRT reverse transcriptase",MaterialKind.CATALYST)
    cdna=graph.product("cdna","mutation-bearing phosphorylated cDNA",
                       MaterialKind.MOLECULAR,state=CDNA)
    graph.react("reverse_transcribe","Read through DMS lesions with TGIRT (57°C, 1.5 h)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((depleted,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
                 (tgirt,InputRole.CATALYST)),(cdna,))
    rnase=graph.input("rnase","NaOH RNA hydrolysis condition",MaterialKind.CONDITION)
    clean_cdna=graph.product("clean_cdna","RNA-free phosphorylated cDNA",
                             MaterialKind.MOLECULAR,state=CDNA)
    rna_waste=graph.product("rna_waste","hydrolyzed RNA template",MaterialKind.MOLECULAR)
    graph.react("remove_rna","Hydrolyze RNA with NaOH (95°C, 3 min)",StepKind.CHEMICAL_CONVERSION,
                ((cdna,InputRole.TEMPLATE),(rnase,InputRole.CONDITION)),
                (clean_cdna,rna_waste))
    graph.mark_terminal(Terminal(rna_waste.id,Disposition.UNRESOLVED,reason="RNA is chemically hydrolyzed; subsequent precipitation/gel recovery is not represented."))
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
    graph.react("library_pcr","Add indexed Illumina arms (9–13 Phusion HF cycles)",StepKind.PCR,
                ((circle,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.exact(
                    PCR_TEMPLATE,
                    PCRPrimer("oNTI231", PCR_FORWARD),
                    PCRPrimer("six-base TruSeq indexing primer", INDEX_PCR,
                              anneal_nt=len(INDEX_PCR_ANNEAL)),
                    statement=("Supplementary Table 1 prints both PCR primers; the "
                               "circle is linearized between the RT handle and P7 so "
                               "their inward-facing sites are checked directly.")))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS,
        acquired_roles=ACQUIRED_ROLES, read_lengths=READ_LENGTHS,
        read_uses=READ_USES))
    graph.declare_reviewed_scope((
        "deadenylase","recj","digest_linker","ribozero","deplete_rrna",
        "rna","dms","modify","zinc","fragment","size_rule","size_select",
        "phosphatase","dephosphorylate","linker","ligase","ligate_linker",
        "rt_primer","tgirt","reverse_transcribe","rnase","remove_rna",
        "circligase","circularize","pcr_forward","pcr_reverse","pcr_polymerase",
        "library_pcr","library"),SOURCE)
    return graph.require_publishable()
