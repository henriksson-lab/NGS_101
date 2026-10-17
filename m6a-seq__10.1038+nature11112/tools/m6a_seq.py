"""Original m6A-seq / MeRIP-seq molecular model (Dominissini et al. 2012).

Source passages (author-hosted PDF of the published Nature article, text twin
Dominissini_Nature_2012.pdf.txt): Methods "RNA preparation" and "m6A-seq"
(lines 957-980) and "Calculation of coverage" (lines 1031-1040).
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, Row
import seqprimers as sp
from terminal_rna import immunoprecipitation_rows
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule,role

TITLE="m6A-seq / MeRIP-seq — antibody enrichment of methylated RNA"
NOTES="01_m6a-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/nature11112">Dominissini et al., <i>Nature</i> (2012)</a>.'
SUMMARY="Poly(A)-selected RNA is chemically fragmented to about 100 nt, immunoprecipitated with an anti-m6A antibody, competitively eluted with free N6-methyladenosine, and converted with an Illumina mRNA sequencing kit into a library alongside the matched input fragments."
CAVEAT=("The Methods name only an “mRNA sequencing kit (Illumina)”, with no catalogue number, adapter, "
        "reverse-transcription or PCR oligos. Library construction is therefore a warned kit boundary and the "
        "final library ends are role tokens, not a reconstructed adapter shell or sample-index design.")
# Kit geometry is not reconstructed from the historical Illumina kit original;
# do not assign canonical TruSeq arms or sample-index reads to this library.
FINAL_LIBRARY=Construct([role("Illumina mRNA-seq kit adapter end"),
                         role("IP or input RNA-fragment-derived cDNA"),
                         role("opposite kit adapter end")],name="source-limited m6A-seq library")
SEQ_PRIMERS=()
FINAL_CAPTION="The IP and matched input use the same unreconstructed kit chemistry; enrichment is biological selection, not a retained molecular barcode."
SEQUENCING_ENDING=("Illumina Genome Analyzer GAIIx, 10 pM template per sample for cluster generation, genomic sequencing kit V2 "
                   "(Methods “m6A-seq”). The coverage analysis treats each read as a sequenced 36-nucleotide stretch "
                   "extended in silico to 100 nt in the 5′→3′ direction (Methods “Calculation of coverage”).")
SEQUENCING_UNAVAILABLE="The historical Illumina mRNA sequencing kit adapters and sequencing-primer sites have not been reconstructed."
READ_SPECS=(sp.read_spec("Read 1",36,
    "Methods: 36-nucleotide sequenced stretch of the IP or input RNA-fragment cDNA; a second read is not described "
    "and strand/orientation relative to the RNA is not stated.","include"),)

def sections(): return [
 ("Select polyadenylated RNA and fragment it",[Row(chunks=[("poly(A)+ RNA → ZnCl2, 94 °C 5 min → ~100-nt RNA fragments",None,False)])],"Fragmentation limits immunoprecipitation enrichment to a local interval around a modified residue."),
 ("Immunoprecipitate m6A-containing fragments",immunoprecipitation_rows("m6A"),"Anti-m6A antibody retains methylated fragments while an aliquot is kept as input."),
 ("Competitively elute with free m6A",[Row(chunks=[("antibody-bound RNA + free N6-methyladenosine → released enriched RNA", "w1",False)])],"Methods: bound RNA is eluted with free N6-methyladenosine in IPP buffer."),
 ("Convert IP and input RNA into libraries",[Row(chunks=[("RNA fragment → Illumina mRNA sequencing kit → library",None,False)])],"IP enrichment is quantified relative to the matched input library."),
]


def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 rna=graph.input("rna","DNase-treated total RNA (PerfectPure kits, 5 Prime)",MaterialKind.BIOLOGICAL,state=molecule("cellular RNA",(role("polyadenylated transcripts"),role("other RNA"))),primary=True)
 poly_dt=graph.input("poly_dt","GenElute mRNA miniprep kit, one round (Sigma-Aldrich)",MaterialKind.SELECTION)
 selected=graph.product("poly_a_rna","poly(A)-selected RNA",MaterialKind.MOLECULAR,state=molecule("poly(A)+ RNA",(role("polyadenylated transcripts"),)))
 depleted=graph.product("non_poly_a","non-polyadenylated RNA",MaterialKind.WASTE)
 graph.react("select_polya","Select polyadenylated RNA",StepKind.SELECTION,((rna,InputRole.TEMPLATE),(poly_dt,InputRole.SELECTION)),(selected,depleted))
 graph.discard(depleted,"is removed by poly(A) selection")
 fragmentation=graph.input("fragmentation_condition","10 mM ZnCl2, 10 mM Tris-HCl pH 7, 94 °C 5 min; EDTA stop and ethanol precipitation",MaterialKind.CONDITION)
 fragments=graph.product("fragments","approximately 100-nt poly(A)-RNA fragments",MaterialKind.MOLECULAR,state=molecule("RNA fragments",(role("~100-nt RNA fragment, with or without m6A"),)))
 graph.react("fragment","Chemically fragment poly(A)-selected RNA",StepKind.FRAGMENTATION,((selected,InputRole.TEMPLATE),(fragmentation,InputRole.CONDITION)),(fragments,))
 ip_input=graph.product("ip_input","RNA aliquot for immunoprecipitation",MaterialKind.MOLECULAR,state=fragments.state)
 control=graph.product("input_control","matched input-control RNA fragments",MaterialKind.MOLECULAR,state=fragments.state)
 graph.react("split","Reserve matched input and IP aliquots",StepKind.SPLIT,((fragments,InputRole.TEMPLATE),),(ip_input,control),changes_molecule=False)
 antibody=graph.input("antibody","affinity-purified anti-m6A polyclonal antibody (Synaptic Systems) in IPP buffer, 2 h 4 °C; protein-A beads (Repligen) 2 h 4 °C",MaterialKind.SELECTION)
 bound=graph.product("bound","bead-bound m6A-containing RNA fragments",MaterialKind.MOLECULAR,state=molecule("antibody-bound RNA",(role("m6A-containing RNA fragment"),)))
 unbound=graph.product("unbound","unbound RNA and washes",MaterialKind.WASTE)
 graph.react("immunoprecipitate","Immunoprecipitate m6A-containing fragments and wash extensively",StepKind.AFFINITY_CAPTURE,((ip_input,InputRole.TEMPLATE),(antibody,InputRole.SELECTION)),(bound,unbound))
 graph.discard(unbound,"is not retained by the antibody–protein-A beads")
 free_m6a=graph.input("free_m6a","free N6-methyladenosine (Sigma-Aldrich) in IPP buffer; ethanol precipitation",MaterialKind.CONDITION)
 enriched=graph.product("enriched","competitively eluted m6A-enriched RNA",MaterialKind.MOLECULAR,state=bound.state)
 graph.react("elute","Competitively elute enriched RNA with free m6A",StepKind.EXTRACTION,((bound,InputRole.TEMPLATE),(free_m6a,InputRole.CONDITION)),(enriched,))
 reviewed=["rna","poly_dt","select_polya","fragmentation_condition","fragment","split","antibody","immunoprecipitate","free_m6a","elute"]
 for prefix,template in (("ip",enriched),("input",control)):
  kit=graph.input(prefix+"_library_kit",prefix+" mRNA sequencing kit reagents (Illumina; catalogue number not given)",MaterialKind.CONDITION)
  library=graph.product(prefix+"_library",prefix+" m6A-seq library",MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (role-only endpoint; historical kit geometry unresolved)")
  graph.react(prefix+"_library_kit_processing","Generate the "+prefix+" library with the Illumina mRNA sequencing kit",StepKind.MAGIC,
              ((template,InputRole.TEMPLATE),(kit,InputRole.CONDITION)),(library,),
              reason="The Methods name only an Illumina mRNA sequencing kit; its reverse transcription, adapter and PCR chemistry is not reconstructed from the kit original.")
  graph.sequence(library,SequencingPlan(limitation=SEQUENCING_ENDING))
  reviewed.extend((kit.id,prefix+"_library_kit_processing",library.id))
 graph.declare_reviewed_scope(reviewed,SOURCE)
 return graph.require_publishable()
