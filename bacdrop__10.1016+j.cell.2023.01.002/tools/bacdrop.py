"""BacDrop molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from microbial_singlecell import bacterial_rna_library, droplet_barcode_rows, indexed_random_rt_scene
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE="BacDrop — droplet bacterial single-cell RNA-seq"
NOTES="01_bacdrop.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2023.01.002">Ma et al., <i>Cell</i> (2023)</a>.'
SUMMARY="Fixed bacteria receive a plate-specific RT index, are pooled at several cells per 10x droplet, and acquire a second bead barcode after TdT poly(A)-tailing of cDNA."
CAVEAT="The two-index cell identity is published, while commercial 10x bead sequence is represented by its molecular role."
FINAL_LIBRARY,SEQ_PRIMERS=bacterial_rna_library("bacdrop","BacDrop RNA library")
FINAL_CAPTION="The plate barcode, droplet barcode and UMI remain physically linked to each bacterial cDNA insert."
SEQUENCING_INTRO="Read cycles recover the two-part cell identity and UMI before transcript sequence."
READ_LENGTHS={"Read 1":60,"Read 2":39,"Index 1 (i7)":8,"Index 2 (i5)":16}
SHOW_READ_CONTENT=True

def oligos():
 return [
  "<b>SMRT_dT second-strand primer</b> — 5′-AAGCAGTGGTATCAACGCAGAGTTTTTTTTTTTTTTTVN-3′",
  "<b>SMRT_PCR cDNA enrichment primer</b> — 5′-AAGCAGTGGTATCAACGCAGAGT-3′",
  "<b>Partial P5 cDNA/library primer</b> — 5′-AATGATACGGCGACCACCGAGA-3′",
 ]
INITIAL_ROWS=(Row(chunks=[("fixed bacterial cells containing total RNA",None,False)]),)
INITIAL_NAME="Fixed bacterial cells"

def sections(): return [
 ("Fix, permeabilize and deplete abundant rRNA in cells",[Row(chunks=[("bacterium: RNA — rRNA probes/RNase H — DNase I cleanup",None,False)])],"Universal rRNA probes and RNase H reduce ribosomal molecules before indexing."),
 ("Install the first barcode by random-primed in-situ RT",indexed_random_rt_scene(method="bacdrop").rows(),"Each plate well contributes one RT barcode and each copied molecule receives a UMI."),
 ("Pool cells and partition several bacteria per 10x droplet",[Row(chunks=[("many plate indexes + one droplet bead barcode", "cbc",False)])],"Combinatorial plate and droplet indexes permit super-loading."),
 ("TdT-tail cDNA and transfer the droplet barcode",droplet_barcode_rows("bacdrop"),"A poly(A) tail substitutes for inefficient bacterial template switching; bead-primed second-strand synthesis transfers the second index."),
 ("Enrich cDNA and complete the Illumina library",[Row(chunks=[("two-index cDNA → enrichment PCR → sequencing library",None,False)])],"The cell-associated cDNA is released and amplified after emulsion break."),
]


def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 cells=graph.input("fixed_cells","fixed bacterial cells",MaterialKind.BIOLOGICAL,
                   state=INITIAL_ROWS,primary=True)
 permeabilizer=graph.input("permeabilizer","cell-wall permeabilization condition",
                           MaterialKind.CONDITION)
 permeable=graph.product("permeable_cells","permeabilized fixed bacteria",MaterialKind.MOLECULAR)
 graph.react("permeabilize","Permeabilize fixed bacterial cells",StepKind.PERMEABILIZATION,
             ((cells,InputRole.TEMPLATE),(permeabilizer,InputRole.CONDITION)),(permeable,))
 probes=graph.input("rrna_probes","universal rRNA depletion probes",MaterialKind.OLIGO)
 hybrid=graph.product("rrna_hybrid","rRNA–depletion-probe hybrids in cells",MaterialKind.MOLECULAR)
 graph.react("hybridize_rrna","Hybridize universal probes to abundant rRNA",
             StepKind.HYBRIDIZATION,
             ((permeable,InputRole.TEMPLATE),(probes,InputRole.STRUCTURAL)),(hybrid,))
 rnaseh=graph.input("rnase_h","RNase H",MaterialKind.CATALYST)
 depleted=graph.product("depleted_cells","rRNA-depleted permeabilized bacteria",
                        MaterialKind.MOLECULAR)
 rrna_waste=graph.product("rrna_waste","cleaved rRNA and depletion probes",MaterialKind.WASTE)
 graph.react("deplete_rrna","Cleave probe-bound rRNA",StepKind.DIGESTION,
             ((hybrid,InputRole.TEMPLATE),(rnaseh,InputRole.CATALYST)),
             (depleted,rrna_waste))
 graph.discard(rrna_waste,"cleaved abundant rRNA is excluded from downstream copying")
 rt_primer=graph.input("indexed_rt_primer","well-indexed random RT primer",MaterialKind.OLIGO)
 rt=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
 indexed_a=graph.product("indexed_cells_a","one group of plate-indexed cells",MaterialKind.MOLECULAR,
                         state=tuple(indexed_random_rt_scene(method="bacdrop").rows()))
 indexed_b=graph.product("indexed_cells_b","other plate-indexed cell groups",MaterialKind.MOLECULAR)
 graph.react("in_situ_rt","Install plate barcode and UMI by in-situ RT",
             StepKind.REVERSE_TRANSCRIPTION,
             ((depleted,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
              (rt,InputRole.CATALYST)),(indexed_a,indexed_b))
 pooled=graph.product("pooled_cells","pooled plate-indexed bacteria",MaterialKind.MOLECULAR)
 graph.react("pool_wells","Pool all plate-indexed cell groups",StepKind.POOL,
             ((indexed_a,InputRole.TEMPLATE),(indexed_b,InputRole.TEMPLATE)),(pooled,),
             changes_molecule=False)
 droplets=graph.input("droplet_condition","10x droplet partitioning condition",MaterialKind.CONDITION)
 occupied=graph.product("occupied_droplets","droplets containing indexed bacteria and a bead",
                        MaterialKind.MOLECULAR)
 empty=graph.product("excluded_droplets","empty or unusable droplets",MaterialKind.WASTE)
 graph.react("partition","Partition pooled bacteria with barcode beads",StepKind.PARTITION,
             ((pooled,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,empty),
             changes_molecule=False)
 graph.discard(empty,"empty and otherwise unusable droplets yield no retained cell library")
 tdt=graph.input("tdt","terminal deoxynucleotidyl transferase",MaterialKind.CATALYST)
 tailed=graph.product("tailed_cdna","poly(A)-tailed indexed cDNA",MaterialKind.MOLECULAR)
 graph.react("tail_cdna","Add poly(A) tails to first-strand cDNA",StepKind.TAILING,
             ((occupied,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(tailed,))
 bead_primer=graph.input("bead_primer","commercial droplet-barcode bead primer",
                         MaterialKind.OLIGO,inferred=True,
                         source=SOURCE+" (commercial barcode bases proprietary)")
 extension_pol=graph.input("extension_polymerase","second-strand DNA polymerase",
                           MaterialKind.CATALYST)
 barcoded=graph.product("dual_barcoded_cdna","plate- and droplet-barcoded cDNA",
                        MaterialKind.MOLECULAR,state=tuple(droplet_barcode_rows("bacdrop")))
 graph.react("transfer_droplet_barcode","Prime the poly(A) tail and transfer the bead barcode",
             StepKind.EXTENSION,
             ((tailed,InputRole.TEMPLATE),(bead_primer,InputRole.STRUCTURAL),
              (extension_pol,InputRole.CATALYST)),(barcoded,))
 fwd=graph.input("pcr_forward","BacDrop library PCR forward primer",MaterialKind.OLIGO)
 rev=graph.input("pcr_reverse","BacDrop library PCR reverse primer",MaterialKind.OLIGO)
 pcr_pol=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","Enrich cDNA and complete the Illumina library",StepKind.PCR,
             ((barcoded,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
              (rev,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The commercial bead sequence is proprietary."))
 graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),
                                       read_lengths=READ_LENGTHS))
 graph.declare_reviewed_scope(("fixed_cells","permeabilizer","permeabilize","rrna_probes",
   "hybridize_rrna","rnase_h","deplete_rrna","indexed_rt_primer","reverse_transcriptase",
   "in_situ_rt","pool_wells","droplet_condition","partition","tdt","tail_cdna",
   "bead_primer","extension_polymerase","transfer_droplet_barcode","pcr_forward",
   "pcr_reverse","pcr_polymerase","library_pcr","library"),SOURCE)
 return graph.require_publishable()
