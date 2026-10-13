"""Spear-ATAC accessibility and guide-capture branches."""
from chemdraw import Row
from multibranch import illumina_branch
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
TITLE="Spear-ATAC — single-cell ATAC CRISPR screening"
NOTES="01_spear-atac.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-021-23213-w">Pierce et al., <i>Nature Communications</i> (2021)</a>.'
SUMMARY="A 10x scATAC-derived workflow captures accessible chromatin and enriches guide identity from the same barcoded nuclei."
CAVEAT="The paper defines the shared cell barcode and the two enrichment branches, but does not print each completed library as one contiguous molecule. Outer sequencing shells are marked inferred."
ATAC,AP=illumina_branch("Spear-ATAC chromatin library","accessible chromatin fragment",architecture="nextera",cell=16,inferred_adapters=True)
GUIDE,GP=illumina_branch("Spear-ATAC guide library","guide-capture amplicon",cell=16,guide=20,inferred_adapters=True)
FINAL_LIBRARIES=(("Chromatin-accessibility library",ATAC,AP,"Tagmented chromatin fragments carry the droplet cell barcode.","Nextera-compatible primers read the ATAC branch."),("Guide-identity library",GUIDE,GP,"Guide enrichment preserves the same droplet cell barcode.","Illumina primers read the guide branch."))
def sections(): return [
 ("Tagment nuclei before partitioning",[Row(chunks=[("accessible chromatin ** Tn5 mosaic ends","me",False)])],"Bulk tagmentation marks open chromatin."),
 ("Barcode nuclei in droplets",[Row(chunks=[("ATAC fragment + gel-bead cell barcode","cbc",False)])],"Gap fill transfers one droplet identity to accessible fragments."),
 ("Amplify the shared barcoded material",[Row(chunks=[("barcoded ATAC products → PCR pool","r1",False)])],"The pool contains genomic ATAC molecules and guide-associated material."),
 ("INFERRED — complete the two sequencing shells",[Row(chunks=[("ATAC library | guide-specific linear PCR → guide library","r2",False)])],"Guide enrichment is source-established; the fully drawn outer adapter completion is inferred."),]

def protocol_graph():
 source="Defining source: https://doi.org/10.1038/s41467-021-23213-w"
 inferred_source=source+" (completed outer library shells are reconstructed)"
 graph=ProtocolGraph("Spear-ATAC",authority=source)
 nuclei=graph.input("nuclei","CRISPR-screen nuclei",MaterialKind.BIOLOGICAL,primary=True)
 transposome=graph.input("transposome","10x scATAC-compatible transposome",MaterialKind.ADAPTER)
 tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=graph.product("tagged_nuclei","nuclei containing tagmented accessible chromatin",
                      MaterialKind.MOLECULAR)
 twaste=graph.product("tagmentation_waste","non-library tagmentation products",MaterialKind.WASTE)
 graph.react("tagment","Tagment nuclei before partitioning",StepKind.TAGMENTATION,
             ((nuclei,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
              (tn5,InputRole.CATALYST)),(tagged,twaste))
 graph.discard(twaste,"fragments lacking usable library ends do not amplify")
 loading=graph.input("gem_loading","10x droplet loading condition",MaterialKind.CONDITION)
 gems=graph.product("productive_gems","GEMs containing nuclei and barcode beads",MaterialKind.MOLECULAR)
 excluded=graph.product("excluded_gems","empty or multiplet GEMs",MaterialKind.WASTE)
 graph.react("partition","Partition tagmented nuclei with gel-bead barcodes",StepKind.PARTITION,
             ((tagged,InputRole.TEMPLATE),(loading,InputRole.CONDITION)),(gems,excluded))
 graph.discard(excluded,"empty and multiplet GEMs are excluded from accepted cells")
 bead=graph.input("bead_primer","gel-bead cell-barcode primer",MaterialKind.OLIGO,
                  inferred=True,source=inferred_source)
 opposite=graph.input("opposite_primer","opposite-end ATAC amplification primer",MaterialKind.OLIGO,
                      inferred=True,source=inferred_source)
 polymerase=graph.input("polymerase","gap-fill and PCR polymerase",MaterialKind.CATALYST)
 shared=graph.product("shared_pool","shared cell-barcoded ATAC/guide-associated product pool",
                      MaterialKind.MOLECULAR,inferred=True,source=inferred_source)
 graph.react("gem_pcr","Gap-fill and amplify the shared barcoded material",StepKind.PCR,
             ((gems,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
              (opposite,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(shared,),
             geometry=PCRGeometry.unavailable(ATAC,
                 reason="The paper defines the shared barcode and branches but not one contiguous pre-branch molecule."))
 atac_in=graph.product("atac_branch","chromatin-accessibility enrichment branch",MaterialKind.MOLECULAR)
 guide_in=graph.product("guide_branch","guide-identity enrichment branch",MaterialKind.MOLECULAR)
 graph.react("branch","Divide the shared pool into ATAC and guide enrichments",StepKind.SPLIT,
             ((shared,InputRole.TEMPLATE),),(atac_in,guide_in))
 atac_f=graph.input("atac_forward","ATAC branch forward primer",MaterialKind.OLIGO,
                    inferred=True,source=inferred_source)
 atac_r=graph.input("atac_reverse","ATAC branch reverse primer",MaterialKind.OLIGO,
                    inferred=True,source=inferred_source)
 atac_lib=graph.product("atac_library",ATAC.name,MaterialKind.LIBRARY,state=ATAC,
                        inferred=True,source=inferred_source)
 graph.react("atac_pcr","Complete the chromatin-accessibility sequencing library",StepKind.PCR,
             ((atac_in,InputRole.TEMPLATE),(atac_f,InputRole.STRUCTURAL),
              (atac_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(atac_lib,),
             geometry=PCRGeometry.unavailable(ATAC,
                 reason="The completed ATAC shell is reconstructed from the 10x-compatible layout."))
 guide_linear=graph.input("guide_linear_primer","guide-specific linear-extension primer",MaterialKind.OLIGO)
 guide_product=graph.product("guide_enriched","guide-specific enriched products",MaterialKind.MOLECULAR)
 graph.react("guide_extension","Linearly enrich guide identity from the shared pool",StepKind.EXTENSION,
             ((guide_in,InputRole.TEMPLATE),(guide_linear,InputRole.STRUCTURAL),
              (polymerase,InputRole.CATALYST)),(guide_product,))
 guide_f=graph.input("guide_forward","guide-library forward PCR primer",MaterialKind.OLIGO,
                     inferred=True,source=inferred_source)
 guide_r=graph.input("guide_reverse","guide-library reverse PCR primer",MaterialKind.OLIGO,
                     inferred=True,source=inferred_source)
 guide_lib=graph.product("guide_library",GUIDE.name,MaterialKind.LIBRARY,state=GUIDE,
                         inferred=True,source=inferred_source)
 graph.react("guide_pcr","Complete the guide-identity sequencing library",StepKind.PCR,
             ((guide_product,InputRole.TEMPLATE),(guide_f,InputRole.STRUCTURAL),
              (guide_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(guide_lib,),
             geometry=PCRGeometry.unavailable(GUIDE,
                 reason="Guide enrichment is source-established; the completed outer shell is reconstructed."))
 for material,library,primers in ((atac_lib,ATAC,AP),(guide_lib,GUIDE,GP)):
  graph.sequence(material,SequencingPlan(
      library=library,primers=primers,acquired_roles=tuple(p.role for p in primers),
      read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in primers},
      limitation="Completed outer sequencing shells are marked inferred."))
 graph.declare_reviewed_scope((
     "nuclei","transposome","tn5","tagment","gem_loading","partition","bead_primer",
     "opposite_primer","polymerase","gem_pcr","branch","atac_forward","atac_reverse",
     "atac_pcr","guide_linear_primer","guide_extension","guide_forward","guide_reverse",
     "guide_pcr","atac_library","guide_library"),source)
 return graph.require_publishable()
