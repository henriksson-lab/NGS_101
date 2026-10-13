"""snapTotal-seq MALBAC/template-switch total-RNA workflow."""
from batch_ngs import seg, truseq_library
from chemdraw import Construct, MolecularState, Row, Workflow, feature
from endprep import (ThreePrimeOverhangAdapter, cohesive_ligation_geometry,
                     dA_tailed_scene, ligation_product_scene, repair_and_dA_tail)
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
TITLE="snapTotal-seq — single-cell total-RNA sequencing"
NOTES="01_snaptotal-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-024-50291-3">Niu et al., <i>Nature Communications</i> (2024)</a>.'
SUMMARY="Low-temperature MALBAC priming and template switching amplify total RNA; a later double-strand conversion adds an 8-nt cell barcode to only one cDNA end."
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("captured RNA sequence","X"*42,placeholder=True),seg("cell barcode","B"*8,"cbc",placeholder=True,feature=feature("cell_barcode","cell_barcode","whitelist",whitelist="Niu et al. Supplementary Data 12")),seg("barcode-to-UMI spacer","L"*35,placeholder=True),seg("UMI","U"*5,"umi",placeholder=True,feature=feature("umi","umi","random"))],"snapTotal-seq standard-primer library")
FINAL_CAPTION="The paper's standard-Illumina option reads RNA sequence in Read 1; Read 2 reaches the 8-nt cell barcode and 5-nt UMI."
SEQUENCING_INTRO="This page shows the paper's standard-primer alternative; the original run configuration also provides custom Read 2 and index primers."
Y_ADAPTER=ThreePrimeOverhangAdapter("snapTotal-seq Y-adapter ligating stem","GCTCTTCCGATCT","GATCGGAAGAGC")
def oligos(): return ["<b>Cell-barcode DSC primer</b> — 5′-GTGTGCTCTTCCGATCTNNNNNNNNAGGAGAGTGTGAGTGATGGTTGAGGATGTGTGGAG-3′"]
def workflow():
 w=Workflow(MolecularState("total RNA",(Row(chunks=[("cellular total RNA",None,False)]),)))
 w.react("Prime total RNA at low temperature",[Row(chunks=[("MALBAC primer → dispersed RNA sites","r1",False)])],name="multiply primed RNA",note="Random multiple annealing avoids dependence on poly(A).")
 w.react("Template-switch and amplify cDNA",[Row(chunks=[("RNA → cDNA with common MALBAC-derived ends","tso",False)])],name="amplified cDNA",note="Ten annealing/extension cycles precede cDNA PCR.")
 w.react("Add one cell barcode by double-strand conversion",[Row(chunks=[("one cDNA end ** 8-nt cell barcode ** UMI-bearing handle","cbc",False)])],name="barcoded cDNA",note="Adding the barcode after PCR avoids complementary barcode arms at both ends.")
 tailed=repair_and_dA_tail(Construct([seg("tagged cDNA fragment","X"*36,placeholder=True)]))
 w.react("Pool, shear, enrich tagged fragments and add 3′ dA",dA_tailed_scene(tailed).rows(),name="dA-tailed tagged fragments",note="A second conversion enriches fragments carrying cell barcode and UMI.")
 geom=cohesive_ligation_geometry(tailed,Y_ADAPTER,statement="Supplementary Data 12 prints both Y-adapter oligos; their ligating stem leaves a terminal 3′-T complementary to each insert 3′-dA.")
 product=ligation_product_scene(FINAL_LIBRARY,"captured RNA sequence","UMI",label="snapTotal-seq library")
 w.ligate_adapter("Ligate the Y adapter",geom,product.rows(),name="indexed Illumina library",note="PCR produces either the custom-primer or standard-primer variant drawn below.")
 return w

def protocol_graph():
 source="Defining source: https://doi.org/10.1038/s41467-024-50291-3"
 graph=ProtocolGraph("snapTotal-seq",authority=source)
 rna=graph.input("total_rna","cellular total RNA",MaterialKind.BIOLOGICAL,primary=True)
 malbac=graph.input("malbac_primer","low-temperature MALBAC primer pool",MaterialKind.OLIGO)
 rt=graph.input("reverse_transcriptase","template-switching reverse transcriptase",MaterialKind.CATALYST)
 cdna=graph.product("first_strand_cdna","multiply primed total-RNA cDNA",MaterialKind.MOLECULAR)
 graph.react("reverse_transcription","Prime dispersed RNA sites and synthesize cDNA",
             StepKind.REVERSE_TRANSCRIPTION,
             ((rna,InputRole.TEMPLATE),(malbac,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
 tso=graph.input("template_switch_oligo","MALBAC template-switch oligo",MaterialKind.OLIGO)
 switched=graph.product("switched_cdna","template-switched cDNA",MaterialKind.MOLECULAR)
 graph.react("template_switch","Template-switch the nascent cDNA",StepKind.TEMPLATE_SWITCHING,
             ((cdna,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(switched,))
 polymerase=graph.input("polymerase","MALBAC, conversion and library PCR polymerase",
                        MaterialKind.CATALYST)
 amplified=graph.product("amplified_cdna","MALBAC-amplified cDNA",MaterialKind.MOLECULAR)
 graph.react("malbac_pcr","Amplify cDNA through MALBAC and PCR cycles",StepKind.PCR,
             ((switched,InputRole.TEMPLATE),(malbac,InputRole.STRUCTURAL),
              (polymerase,InputRole.CATALYST)),(amplified,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The MALBAC primer family acts through repeated common ends rather than the final library shell.",symmetric_species=True))
 barcode_primer=graph.input("barcode_dsc_primer","cell-barcode double-strand-conversion primer",
                            MaterialKind.OLIGO)
 barcoded=graph.product("barcoded_cdna","one-ended cell-barcoded cDNA",MaterialKind.MOLECULAR)
 graph.react("barcode_conversion","Add one cell barcode by double-strand conversion",
             StepKind.EXTENSION,
             ((amplified,InputRole.TEMPLATE),(barcode_primer,InputRole.STRUCTURAL),
              (polymerase,InputRole.CATALYST)),(barcoded,))
 shear=graph.input("shearing_condition","physical shearing condition",MaterialKind.CONDITION)
 sheared=graph.product("sheared_fragments","sheared barcoded cDNA fragments",MaterialKind.MOLECULAR)
 graph.react("shear","Pool and shear the barcoded cDNA",StepKind.FRAGMENTATION,
             ((barcoded,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(sheared,))
 enrichment=graph.input("tag_selection","tagged-fragment enrichment reagent",MaterialKind.SELECTION)
 enriched=graph.product("enriched_fragments","tag-bearing cDNA fragments",MaterialKind.MOLECULAR)
 rejected=graph.product("untagged_fragments","untagged cDNA fragments",MaterialKind.WASTE)
 graph.react("enrich","Enrich fragments carrying the cell-barcode/UMI tag",
             StepKind.SELECTION,
             ((sheared,InputRole.TEMPLATE),(enrichment,InputRole.SELECTION)),(enriched,rejected))
 graph.discard(rejected,"untagged sheared fragments are removed before library construction")
 end_repair=graph.input("end_repair_enzymes","end-repair enzymes",MaterialKind.CATALYST)
 repaired=graph.product("repaired_fragments","end-repaired tagged fragments",MaterialKind.MOLECULAR)
 graph.react("end_repair","End-repair enriched tagged fragments",StepKind.END_REPAIR,
             ((enriched,InputRole.TEMPLATE),(end_repair,InputRole.CATALYST)),(repaired,))
 tailing=graph.input("tailing_polymerase","dA-tailing polymerase",MaterialKind.CATALYST)
 tailed_construct=repair_and_dA_tail(Construct([seg("tagged cDNA fragment","X"*36,placeholder=True)]))
 tailed=graph.product("tailed_fragments","dA-tailed tagged fragments",MaterialKind.MOLECULAR,
                      state=dA_tailed_scene(tailed_construct))
 graph.react("da_tail","Add 3′ dA overhangs",StepKind.TAILING,
             ((repaired,InputRole.TEMPLATE),(tailing,InputRole.CATALYST)),(tailed,))
 adapter=graph.input("y_adapter","published Y adapter",MaterialKind.ADAPTER)
 ligase=graph.input("ligase","DNA ligase",MaterialKind.CATALYST)
 precursor=graph.product("adapter_ligated","Y-adapter-ligated library precursor",
                         MaterialKind.MOLECULAR)
 graph.react("adapter_ligation","Ligate the Y adapter",StepKind.LIGATION,
             ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
              (ligase,InputRole.CATALYST)),(precursor,),
             geometry=cohesive_ligation_geometry(
                 tailed_construct,Y_ADAPTER,
                 statement="Supplementary Data 12 prints both Y-adapter strands and their terminal 3′-T."))
 p5=graph.input("p5_primer","standard-variant P5 PCR primer",MaterialKind.OLIGO)
 p7=graph.input("p7_primer","standard-variant indexed P7 PCR primer",MaterialKind.OLIGO)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","PCR-complete the standard-primer library variant",StepKind.PCR,
             ((precursor,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
              (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The final standard-primer alternative is represented as a completed construct."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),
     read_uses={p.role:("demultiplex" if "Index" in p.role else "include") for p in SEQ_PRIMERS}))
 graph.declare_reviewed_scope((
     "total_rna","malbac_primer","reverse_transcriptase","reverse_transcription",
     "template_switch_oligo","template_switch","polymerase","malbac_pcr",
     "barcode_dsc_primer","barcode_conversion","shearing_condition","shear",
     "tag_selection","enrich","end_repair_enzymes","end_repair","tailing_polymerase",
     "da_tail","y_adapter","ligase","adapter_ligation","p5_primer","p7_primer",
     "library_pcr","library"),source)
 return graph.require_publishable()
