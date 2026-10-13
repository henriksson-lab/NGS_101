"""Di et al. 2020 SHERRY chemistry."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]; sys.path[:0]=[str(ROOT/"lib")]

from batch_ngs import nextera_library, seg
from chemdraw import MolecularState, Scene, Workflow
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_prep import rna_template, rna_dna_hybrid
import nextera as nx

TITLE="SHERRY — direct tagmentation of RNA/DNA hybrids"
NOTES="01_sherry.html"
SOURCE='Defining source: <a href="https://doi.org/10.1073/pnas.1919800117">Di et al., <i>PNAS</i> (2020)</a>.'
SUMMARY="Oligo-dT reverse transcription makes an RNA/cDNA hybrid; Tn5 directly fragments and tags that hybrid before gap extension and indexed PCR."

INITIAL_NAME="poly(A) RNA"
INITIAL_ROWS=rna_template("poly(A) RNA", 36)
FINAL_LIBRARY, SEQ_PRIMERS=nextera_library([seg("transcript-derived insert","X"*38,placeholder=True)],"SHERRY library")
FINAL_CAPTION="The amplifiable S5–insert–S7 product after gap extension and indexed PCR."
SEQUENCING_INTRO="The paper used Illumina NextSeq 500 or HiSeq 4000; standard Nextera primer sites are installed by the published transposomes and indexed PCR primers."
def _tagged():
    sc=Scene()
    sc.strand("tagged RNA",[seg("transferred S5–ME",nx.ADAPTOR_S5,"s5+me"),seg("RNA fragment","X"*30,placeholder=True)],label="tagged RNA product")
    sc.strand("tagged cDNA",[seg("transferred S7–ME",nx.ADAPTOR_S7,"s7+me"),seg("cDNA fragment","X"*30,placeholder=True)],label="tagged cDNA product")
    sc.labels("tagged RNA"); sc.labels("tagged cDNA")
    return sc.rows()

def workflow():
    w=Workflow(MolecularState(INITIAL_NAME,tuple(INITIAL_ROWS)))
    for action,rows,note in [
      ("Oligo-dT reverse transcription",rna_dna_hybrid(length=36),"SuperScript II copies poly(A) RNA; no second-strand synthesis is performed."),
      ("Tn5 tagment RNA/DNA hybrid",_tagged(),"Strand tests demonstrated adapter transfer to both fragmented RNA and cDNA; these are transferred-strand products, not a gap-filled duplex."),
      ("Gap extension and indexed PCR",Scene.duplex(list(FINAL_LIBRARY),label="sequencing library").rows(),"Polymerase extension makes the tagmented cDNA amplifiable; indexed common primers complete P5 and P7."),
    ]: w.react(action,rows,note=note)
    return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("rna","poly(A) RNA",MaterialKind.MOLECULAR,primary=True)
    oligo_dt=graph.input("oligo_dt","poly(T)30VN reverse-transcription primer",
                         MaterialKind.OLIGO)
    rt=graph.input("reverse_transcriptase","SuperScript II",MaterialKind.CATALYST)
    hybrid=graph.product("rna_cdna_hybrid","RNA/cDNA hybrid",MaterialKind.MOLECULAR,
                         state=tuple(rna_dna_hybrid(length=36)))
    graph.react("reverse_transcribe","Make an RNA/cDNA hybrid without second-strand synthesis",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna,InputRole.TEMPLATE),(oligo_dt,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(hybrid,))
    transposome=graph.input("loaded_transposome","S5/S7-loaded Tn5 transposomes",
                            MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagged_hybrid","adapter-transferred RNA/cDNA fragments",
                         MaterialKind.MOLECULAR,state=tuple(_tagged()))
    nonlibrary=graph.product("nonlibrary_tagmentation","non-amplifiable tagmentation products",
                             MaterialKind.WASTE)
    graph.react("tagment","Fragment and tag the RNA/cDNA hybrid",StepKind.TAGMENTATION,
                ((hybrid,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,nonlibrary))
    graph.discard(nonlibrary,"homotypic or otherwise non-amplifiable tagmentation products")
    extension_polymerase=graph.input("extension_polymerase",
                                     "SuperScript II or Bst 2.0 gap-extension enzyme",
                                     MaterialKind.CATALYST)
    extended=graph.product("extended","gap-extended tagmented cDNA",
                           MaterialKind.MOLECULAR)
    graph.react("gap_extend","Extend the transferred adapter junctions",StepKind.EXTENSION,
                ((tagged,InputRole.TEMPLATE),(extension_polymerase,InputRole.CATALYST)),
                (extended,))
    i5=graph.input("i5_primer","indexed Nextera i5 PCR primer",MaterialKind.OLIGO)
    i7=graph.input("i7_primer","indexed Nextera i7 PCR primer",MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    graph.react("index_pcr","Select S5/S7 products and complete P5/P7 arms",
                StepKind.PCR,
                ((extended,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The final Nextera sites are exact, but one experiment-specific i5/i7 primer pair is not singled out."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "rna","oligo_dt","reverse_transcriptase","reverse_transcribe",
        "loaded_transposome","tn5","tagment","extension_polymerase","gap_extend",
        "i5_primer","i7_primer","pcr_polymerase","index_pcr","library"),SOURCE)
    return graph.require_publishable()
