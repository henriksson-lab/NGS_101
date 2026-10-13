"""Original m6A-seq / MeRIP-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Row
from terminal_rna import immunoprecipitation_rows
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule,role

TITLE="m6A-seq / MeRIP-seq — antibody enrichment of methylated RNA"
NOTES="01_m6a-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nature11112">Dominissini et al., <i>Nature</i> (2012)</a>.'
SUMMARY="Poly(A)-selected RNA is fragmented, enriched with anti-m6A antibody, competitively eluted with free m6A, and converted into an Illumina cDNA library alongside an input control."
CAVEAT="The paper defines the enrichment assay but does not establish a unique adapter sequence architecture; the final conventional Illumina arms are therefore marked inferred."
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("m6A-enriched RNA-derived cDNA","X"*56,placeholder=True)],"m6A-seq IP library",inferred_adapters=True)
FINAL_CAPTION="The IP and matched input use the same cDNA-library chemistry; enrichment is biological selection, not a retained molecular barcode."
SEQUENCING_INTRO="Standard Illumina primers read the enriched cDNA insert and sample indexes."
READ_LENGTHS={"Read 1":50,"Read 2":50}
SHOW_READ_CONTENT=True

def sections(): return [
 ("Select polyadenylated RNA and fragment it",[Row(chunks=[("poly(A)+ RNA → approximately 100-nt RNA fragments",None,False)])],"Fragmentation limits immunoprecipitation enrichment to a local interval around a modified residue."),
 ("Immunoprecipitate m6A-containing fragments",immunoprecipitation_rows("m6A"),"Anti-m6A antibody retains methylated fragments while an aliquot is kept as input."),
 ("Competitively elute with free m6A",[Row(chunks=[("antibody-bound RNA + excess free m6A → released enriched RNA", "w1",False)])],"The supplementary validation identifies free-m6A competition as the specific elution."),
 ("Convert IP and input RNA into cDNA libraries",[Row(chunks=[("RNA fragment → reverse transcription → adapter-complete cDNA",None,False)])],"IP enrichment is quantified relative to the matched input library."),
]


def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 rna=graph.input("rna","total cellular RNA",MaterialKind.BIOLOGICAL,state=molecule("cellular RNA",(role("polyadenylated transcripts"),role("other RNA"))),primary=True)
 poly_dt=graph.input("poly_dt","poly(dT) selection matrix",MaterialKind.SELECTION)
 selected=graph.product("poly_a_rna","poly(A)-selected RNA",MaterialKind.MOLECULAR,state=rna.state)
 depleted=graph.product("non_poly_a","non-polyadenylated RNA",MaterialKind.WASTE)
 graph.react("select_polya","Select polyadenylated RNA",StepKind.SELECTION,((rna,InputRole.TEMPLATE),(poly_dt,InputRole.SELECTION)),(selected,depleted))
 graph.discard(depleted,"is removed by poly(A) selection")
 fragmentation=graph.input("fragmentation_condition","RNA fragmentation condition",MaterialKind.CONDITION)
 fragments=graph.product("fragments","approximately 100-nt poly(A)-RNA fragments",MaterialKind.MOLECULAR,state=selected.state)
 graph.react("fragment","Fragment poly(A)-selected RNA",StepKind.FRAGMENTATION,((selected,InputRole.TEMPLATE),(fragmentation,InputRole.CONDITION)),(fragments,))
 ip_input=graph.product("ip_input","RNA aliquot for immunoprecipitation",MaterialKind.MOLECULAR,state=fragments.state)
 control=graph.product("input_control","matched input-control RNA aliquot",MaterialKind.MOLECULAR,state=fragments.state)
 graph.react("split","Reserve matched input and IP aliquots",StepKind.SPLIT,((fragments,InputRole.TEMPLATE),),(ip_input,control),changes_molecule=False)
 antibody=graph.input("antibody","anti-m6A antibody",MaterialKind.SELECTION)
 bound=graph.product("bound","antibody-bound m6A-containing RNA fragments",MaterialKind.MOLECULAR,state=ip_input.state)
 unbound=graph.product("unbound","m6A-depleted IP flow-through",MaterialKind.WASTE)
 graph.react("immunoprecipitate","Immunoprecipitate m6A-containing fragments",StepKind.AFFINITY_CAPTURE,((ip_input,InputRole.TEMPLATE),(antibody,InputRole.SELECTION)),(bound,unbound))
 graph.discard(unbound,"is not retained by the anti-m6A antibody")
 free_m6a=graph.input("free_m6a","free m6A competitive-elution condition",MaterialKind.CONDITION)
 enriched=graph.product("enriched","competitively eluted m6A-enriched RNA",MaterialKind.MOLECULAR,state=bound.state)
 graph.react("elute","Competitively elute enriched RNA with free m6A",StepKind.EXTRACTION,((bound,InputRole.TEMPLATE),(free_m6a,InputRole.CONDITION)),(enriched,))
 reviewed=["rna","poly_dt","select_polya","fragmentation_condition","fragment","split","antibody","immunoprecipitate","free_m6a","elute"]
 for prefix,template in (("ip",enriched),("input",control)):
  rt_primer=graph.input(prefix+"_rt_primer",prefix+" library RT primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (library kit not specified)")
  rt=graph.input(prefix+"_rt",prefix+" reverse transcriptase",MaterialKind.CATALYST)
  cdna=graph.product(prefix+"_cdna",prefix+" RNA-derived cDNA",MaterialKind.MOLECULAR,state=template.state)
  graph.react(prefix+"_reverse_transcribe","Reverse-transcribe the "+prefix+" RNA",StepKind.REVERSE_TRANSCRIPTION,((template,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
  fwd=graph.input(prefix+"_forward",prefix+" library PCR forward primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (library kit not specified)")
  rev=graph.input(prefix+"_reverse",prefix+" library PCR reverse primer",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (library kit not specified)")
  pol=graph.input(prefix+"_polymerase",prefix+" library PCR polymerase",MaterialKind.CATALYST)
  library=graph.product(prefix+"_library",prefix+" m6A-seq library",MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (canonical compatible shell inferred)")
  graph.react(prefix+"_library_pcr","Build the "+prefix+" cDNA library",StepKind.PCR,((cdna,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason="The defining paper does not establish a unique adapter/primer architecture."))
  graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),read_lengths=READ_LENGTHS))
  reviewed.extend((rt_primer.id,rt.id,prefix+"_reverse_transcribe",fwd.id,rev.id,pol.id,prefix+"_library_pcr",library.id))
 graph.declare_reviewed_scope(reviewed,SOURCE)
 return graph.require_publishable()
