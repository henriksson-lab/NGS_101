"""M3-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from adapter_ligation import AdapterLigationGeometry
from chemdraw import MolecularState,Row,Workflow
from microbial_singlecell import bacterial_rna_library,droplet_barcode_rows,indexed_random_rt_scene
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE="M3-seq — multiplexed microbial single-cell RNA-seq"
NOTES="01_m3-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41564-023-01462-3">McNulty et al., <i>Nature Microbiology</i> (2023)</a>.'
SUMMARY="Random in-situ RT supplies a plate barcode and UMI; 10x ATAC droplets ligate a second barcode, followed by Tn5/T7 amplification and post hoc RNase-H rRNA depletion."
CAVEAT="The page preserves the important RNA→cDNA→T7 RNA→cDNA transitions. Proprietary droplet-barcode bases remain role tokens."
FINAL_LIBRARY,SEQ_PRIMERS=bacterial_rna_library("m3seq","M3-seq RNA library")
FINAL_CAPTION="The final cDNA library retains the round-one barcode, round-two barcode and UMI through T7 amplification and depletion."
SEQUENCING_INTRO="Read cycles recover the combinatorial cell identity, UMI and bacterial transcript insert."
READ_LENGTHS={"Read 1":28,"Read 2":90}
SHOW_READ_CONTENT=True
INITIAL_ROWS=(Row(chunks=[("fixed bacterial cells containing total RNA",None,False)]),)
INITIAL_NAME="Fixed bacterial cells"

def sections():
 wf=Workflow(MolecularState("Fixed bacteria",(Row(chunks=[("fixed bacterial RNA",None,False)]),)))
 wf.react("Index fixed bacteria by random-primed in-situ RT",indexed_random_rt_scene(method="m3seq").rows(),note="The RT primer carries BC1 and a UMI.")
 wf.ligate_adapter("Ligate the second index in 10x ATAC droplets",AdapterLigationGeometry.unavailable("M3-seq uses the commercial 10x ATAC droplet barcode reagent for BC2 ligation; the paper does not publish its oligo or ligation-ready overhang geometry."),droplet_barcode_rows("m3seq"),note="BC2 is ligated to cell-associated cDNA; BC1+BC2 defines a cell.")
 wf.react("Lyse cells and synthesize the second cDNA strand",[Row(chunks=[("single-stranded indexed cDNA → Klenow/random-primer dsDNA",None,False)])],note="Bulk second-strand synthesis makes the library amplifiable.")
 wf.react("Tagment, PCR-add a T7 promoter and transcribe",[Row(chunks=[("dsDNA → Tn5 fragment → T7-promoter PCR → amplified RNA","me",False)])],note="T7 transcription linearly amplifies indexed library molecules.")
 wf.react("Deplete rRNA library molecules",[Row(chunks=[("RNA library + rRNA DNA probes → RNase H cleavage",None,False)])],note="Post hoc depletion acts on amplified RNA rather than intact cells.")
 wf.react("Reverse-transcribe depleted RNA and add sequencing adapters",[Row(chunks=[("depleted RNA → cDNA → final adapter PCR",None,False)])],note="A final RT and PCR return the selected molecules to an Illumina-compatible cDNA library.")
 return [(r.action,r.after.rows,r.note) for r in wf.reactions]


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("fixed_cells","fixed permeabilized bacterial cells",MaterialKind.BIOLOGICAL,
               state=INITIAL_ROWS,primary=True)
 rt_primer=g.input("indexed_rt_primer","BC1/UMI random RT primer",MaterialKind.OLIGO)
 rt=g.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
 cdna=g.product("bc1_cdna","BC1- and UMI-tagged bacterial cDNA",MaterialKind.MOLECULAR,
                state=tuple(indexed_random_rt_scene(method="m3seq").rows()))
 g.react("in_situ_rt","Index fixed bacteria by random-primed in-situ RT",
         StepKind.REVERSE_TRANSCRIPTION,
         ((cells,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),
         (cdna,))
 droplet=g.input("droplet_condition","10x ATAC droplet partitioning condition",
                 MaterialKind.CONDITION)
 occupied=g.product("occupied_droplets","droplets containing BC1-tagged cells and beads",
                    MaterialKind.MOLECULAR)
 empty=g.product("excluded_droplets","empty or unusable droplets",MaterialKind.WASTE)
 g.react("partition","Partition indexed cells into 10x ATAC droplets",StepKind.PARTITION,
         ((cdna,InputRole.TEMPLATE),(droplet,InputRole.CONDITION)),(occupied,empty),
         changes_molecule=False)
 g.discard(empty,"empty and unusable droplets do not yield retained molecules")
 bc2=g.input("bc2_adapter","commercial 10x ATAC BC2 adapter",MaterialKind.ADAPTER,
             inferred=True,source=SOURCE+" (commercial oligo sequence proprietary)")
 ligase=g.input("dna_ligase","DNA ligase",MaterialKind.CATALYST)
 dual=g.product("dual_barcoded_cdna","BC1/BC2-tagged bacterial cDNA",MaterialKind.MOLECULAR,
                state=tuple(droplet_barcode_rows("m3seq")))
 g.react("ligate_bc2","Ligate the second cell index in droplets",StepKind.LIGATION,
         ((occupied,InputRole.TEMPLATE),(bc2,InputRole.STRUCTURAL),
          (ligase,InputRole.CATALYST)),(dual,),
         geometry=AdapterLigationGeometry.unavailable(
             "The paper does not publish the commercial barcode oligo or overhang geometry."))
 random_primer=g.input("random_second_strand_primer","random second-strand primer",
                       MaterialKind.OLIGO)
 klenow=g.input("klenow","Klenow DNA polymerase",MaterialKind.CATALYST)
 ds=g.product("double_stranded_cdna","double-stranded dual-indexed cDNA",MaterialKind.MOLECULAR)
 g.react("second_strand","Synthesize the second cDNA strand",StepKind.EXTENSION,
         ((dual,InputRole.TEMPLATE),(random_primer,InputRole.STRUCTURAL),
          (klenow,InputRole.CATALYST)),(ds,))
 tn5_adapter=g.input("tn5_adapter","Tn5 mosaic-end adapter",MaterialKind.ADAPTER)
 tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=g.product("tagmented_cdna","Tn5-tagged double-stranded cDNA",MaterialKind.MOLECULAR)
 tag_waste=g.product("tagment_waste","non-amplifiable tagmentation products",MaterialKind.WASTE)
 g.react("tagment","Tagment double-stranded indexed cDNA",StepKind.TAGMENTATION,
         ((ds,InputRole.TEMPLATE),(tn5_adapter,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
         (tagged,tag_waste))
 g.discard(tag_waste,"products lacking the required amplifiable transfer-end combination")
 t7_fwd=g.input("t7_forward","T7-promoter PCR primer",MaterialKind.OLIGO)
 t7_rev=g.input("t7_reverse","T7-amplification reverse primer",MaterialKind.OLIGO)
 pol=g.input("t7_pcr_polymerase","T7-promoter PCR polymerase",MaterialKind.CATALYST)
 t7_template=g.product("t7_template","T7-promoter-bearing indexed DNA",MaterialKind.MOLECULAR)
 g.react("t7_pcr","PCR-install the T7 promoter",StepKind.PCR,
         ((tagged,InputRole.TEMPLATE),(t7_fwd,InputRole.STRUCTURAL),
          (t7_rev,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(t7_template,),
         geometry=PCRGeometry.unavailable(
             FINAL_LIBRARY,reason="The intermediate T7-promoter template is not reconstructed base-for-base."))
 t7=g.input("t7_polymerase","T7 RNA polymerase",MaterialKind.CATALYST)
 amplified_rna=g.product("amplified_rna","linearly amplified indexed RNA",MaterialKind.MOLECULAR)
 g.react("transcribe","Linearly amplify the indexed library by T7 transcription",
         StepKind.TRANSCRIPTION,
         ((t7_template,InputRole.TEMPLATE),(t7,InputRole.CATALYST)),(amplified_rna,))
 rrna_probe=g.input("rrna_probe","rRNA depletion DNA probes",MaterialKind.OLIGO)
 hybrid=g.product("rrna_hybrid","probe-hybridized amplified RNA",MaterialKind.MOLECULAR)
 g.react("hybridize_rrna","Hybridize depletion probes to rRNA-derived molecules",
         StepKind.HYBRIDIZATION,
         ((amplified_rna,InputRole.TEMPLATE),(rrna_probe,InputRole.STRUCTURAL)),(hybrid,))
 rnaseh=g.input("rnase_h","RNase H",MaterialKind.CATALYST)
 depleted=g.product("depleted_rna","rRNA-depleted amplified RNA",MaterialKind.MOLECULAR)
 rrna_waste=g.product("rrna_waste","RNase-H-cleaved rRNA library molecules",MaterialKind.WASTE)
 g.react("deplete_rrna","Cleave rRNA-derived amplified molecules",StepKind.DIGESTION,
         ((hybrid,InputRole.TEMPLATE),(rnaseh,InputRole.CATALYST)),(depleted,rrna_waste))
 g.discard(rrna_waste,"cleaved rRNA-derived library molecules are not reverse-transcribed")
 final_rt_primer=g.input("final_rt_primer","final reverse-transcription primer",MaterialKind.OLIGO)
 final_cdna=g.product("final_cdna","depleted second-generation cDNA",MaterialKind.MOLECULAR)
 g.react("final_rt","Reverse-transcribe the depleted amplified RNA",
         StepKind.REVERSE_TRANSCRIPTION,
         ((depleted,InputRole.TEMPLATE),(final_rt_primer,InputRole.STRUCTURAL),
          (rt,InputRole.CATALYST)),(final_cdna,))
 fwd=g.input("library_forward","final library PCR forward primer",MaterialKind.OLIGO)
 rev=g.input("library_reverse","final library PCR reverse primer",MaterialKind.OLIGO)
 final_pol=g.input("library_polymerase","final library PCR polymerase",MaterialKind.CATALYST)
 library=g.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 g.react("library_pcr","Add sequencing adapters to the depleted cDNA",StepKind.PCR,
         ((final_cdna,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
          (rev,InputRole.STRUCTURAL),(final_pol,InputRole.CATALYST)),(library,),
         geometry=PCRGeometry.unavailable(
             FINAL_LIBRARY,reason="The commercial barcode arm prevents a fully exact compact PCR geometry."))
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),
                                   read_lengths=READ_LENGTHS))
 g.declare_reviewed_scope(("fixed_cells","indexed_rt_primer","reverse_transcriptase",
   "in_situ_rt","droplet_condition","partition","bc2_adapter","dna_ligase","ligate_bc2",
   "random_second_strand_primer","klenow","second_strand","tn5_adapter","tn5","tagment",
   "t7_forward","t7_reverse","t7_pcr_polymerase","t7_pcr","t7_polymerase","transcribe",
   "rrna_probe","hybridize_rrna","rnase_h","deplete_rrna","final_rt_primer","final_rt",
   "library_forward","library_reverse","library_polymerase","library_pcr","library"),SOURCE)
 return g.require_publishable()
