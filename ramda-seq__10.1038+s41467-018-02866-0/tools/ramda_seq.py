"""Hayashi et al. 2018 RamDA-seq chemistry."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]; sys.path[:0]=[str(ROOT/"lib")]
from batch_ngs import nextera_library, seg
from chemdraw import MolecularState, Scene, Segment, Workflow, oligo, revcomp
from rna_prep import rna_template, rna_dna_hybrid, single_strand
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE="RamDA-seq — single-cell full-length total RNA"
NOTES="01_ramda-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-018-02866-0">Hayashi et al., <i>Nature Communications</i> (2018)</a>.'
SUMMARY="Not-so-random and oligo-dT primers initiate RT; DNase I nicks cDNA while RNase-H-minus RT displaces and amplifies it, then second-strand synthesis and Nextera XT make the sequencing library."
INITIAL_NAME="single-cell total RNA"
INITIAL_ROWS=rna_template("total RNA",38)
FINAL_LIBRARY, SEQ_PRIMERS=nextera_library([seg("RamDA cDNA insert","X"*38,placeholder=True)],"RamDA-seq library")
FINAL_CAPTION="Nextera XT library made from double-stranded RT-RamDA product."
SEQUENCING_INTRO="The defining study used single-read 76-cycle NextSeq sequencing; the page also locates the complete standard Nextera primer set present on the finished library."
READ_LENGTHS={"Read 1":76}
ACQUIRED_ROLES=("Read 1",)
SHOW_READ_CONTENT=True
FIRST_NSR_1 = "AAAGCA"  # Supplementary Data 3, 1st_MmNSR_1
SECOND_NSR_1 = revcomp(FIRST_NSR_1)  # 2nd_MmNSR_1 in the same table

def oligos():
    yield oligo("1st_MmNSR_1", [Segment("first-strand NSR", FIRST_NSR_1)])
    yield oligo("2nd_MmNSR_1", [Segment("second-strand NSR", SECOND_NSR_1)])
    yield oligo("oligo(dT)18", [Segment("poly-T capture primer", "T" * 18)])

def workflow():
    amplified=single_strand("displaced cDNA",[seg("gp32-protected displaced cDNA","X"*46,placeholder=True)],"displaced cDNA products")
    ds=Scene.duplex([seg("amplified cDNA","X"*42,placeholder=True)],label="double-stranded cDNA").rows()
    w=Workflow(MolecularState(INITIAL_NAME,tuple(INITIAL_ROWS)))
    for action,rows,note in [
      ("NSR + oligo-dT reverse transcription",rna_dna_hybrid(length=38),"A pool of 408 first-strand NSR hexamers avoids exact rRNA matches while oligo(dT) captures poly(A) RNA."),
      ("DNase-I nicking and RT displacement",amplified,"DNase I nicks cDNA; RNase-H-minus PrimeScript extends from nicked 3' ends, while T4 gp32 promotes displacement and protects products."),
      ("Second-strand NSR synthesis",ds,"Complementary second-strand NSRs prime Klenow Fragment (3'→5' exo−)."),
      ("Nextera XT tagmentation and PCR",Scene.duplex(list(FINAL_LIBRARY),label="sequencing library").rows(),"The study used scaled Nextera XT and 13 PCR cycles (14 in its advanced method)."),
    ]: w.react(action,rows,note=note)
    return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("rna","single-cell total RNA",MaterialKind.MOLECULAR,primary=True)
    nsr=graph.input("nsr_primers","408-member not-so-random first-strand primer pool",
                    MaterialKind.OLIGO)
    oligo_dt=graph.input("oligo_dt","oligo-dT primer",MaterialKind.OLIGO)
    rt=graph.input("reverse_transcriptase","RNase-H-minus PrimeScript reverse transcriptase",
                   MaterialKind.CATALYST)
    hybrid=graph.product("hybrid","random- and oligo-dT-primed RNA/cDNA hybrids",
                         MaterialKind.MOLECULAR,state=tuple(rna_dna_hybrid(length=38)))
    graph.react("reverse_transcribe","Prime total RNA with NSR and oligo-dT pools",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna,InputRole.TEMPLATE),(nsr,InputRole.STRUCTURAL),
                 (oligo_dt,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(hybrid,))
    dnase=graph.input("dnase","DNase I",MaterialKind.CATALYST)
    nicked=graph.product("nicked","DNase-I-nicked cDNA hybrids",MaterialKind.MOLECULAR)
    graph.react("nick","Nick first-strand cDNA",StepKind.DIGESTION,
                ((hybrid,InputRole.TEMPLATE),(dnase,InputRole.CATALYST)),(nicked,))
    gp32=graph.input("gp32","T4 gene-32 single-strand-binding protein",
                     MaterialKind.CONDITION)
    displaced=graph.product("displaced","displacement-amplified cDNA products",
                            MaterialKind.MOLECULAR)
    graph.react("displace","Extend from cDNA nicks with strand displacement",
                StepKind.EXTENSION,
                ((nicked,InputRole.TEMPLATE),(rt,InputRole.CATALYST),
                 (gp32,InputRole.CONDITION)),(displaced,))
    second_nsr=graph.input("second_nsr","second-strand NSR primer pool",
                           MaterialKind.OLIGO)
    klenow=graph.input("klenow","Klenow Fragment exo-minus",MaterialKind.CATALYST)
    duplex=graph.product("duplex","double-stranded RT-RamDA product",
                         MaterialKind.MOLECULAR)
    graph.react("second_strand","Synthesize the second RamDA strand",StepKind.EXTENSION,
                ((displaced,InputRole.TEMPLATE),(second_nsr,InputRole.STRUCTURAL),
                 (klenow,InputRole.CATALYST)),(duplex,))
    transposome=graph.input("transposome","Nextera XT loaded transposomes",
                            MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagged","tagmented RamDA cDNA",MaterialKind.MOLECULAR)
    waste=graph.product("tagment_waste","non-amplifiable tagmentation products",
                        MaterialKind.WASTE)
    graph.react("tagment","Tagment double-stranded RamDA product",StepKind.TAGMENTATION,
                ((duplex,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.discard(waste,"non-amplifiable Tn5 outcomes are not sequenced")
    i5=graph.input("i5_primer","Nextera i5 indexing primer",MaterialKind.OLIGO)
    i7=graph.input("i7_primer","Nextera i7 indexing primer",MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("index_pcr","Complete and index the Nextera library",StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="Index identities vary while the final Nextera sites are exact."))
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),read_lengths=READ_LENGTHS,
        acquired_roles=ACQUIRED_ROLES))
    graph.declare_reviewed_scope((
        "rna","nsr_primers","oligo_dt","reverse_transcriptase","reverse_transcribe",
        "dnase","nick","gp32","displace","second_nsr","klenow",
        "second_strand","transposome","tn5","tagment","i5_primer","i7_primer",
        "pcr_polymerase","index_pcr","library"),SOURCE)
    return graph.require_publishable()
