"""BacDrop molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Row
from batch_ngs import seg, truseq_library
from chemdraw import Scene, feature
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

def indexed_random_rt_scene(*, method: str, barcode_len: int = 13, umi_len: int = 8) -> Scene:
    """In-cell random-primed RT carrying a plate barcode and UMI."""
    rna = [seg("bacterial RNA", "X" * 34, placeholder=True)]
    primer = [
        seg("PCR handle", "H" * 16, "r1", placeholder=True),
        seg("plate barcode", "B" * barcode_len, "cbc", placeholder=True,
            feature=feature(method + "_plate", "cell_barcode", "combinatorial",
                            group=method + "_cell", part="plate well")),
        seg("UMI", "U" * umi_len, "umi", placeholder=True,
            feature=feature(method + "_umi", "umi", "random")),
        seg("random primer", "N" * 6, placeholder=True),
    ]
    sc = Scene(); sc.strand("RNA", rna, label="fixed, permeabilized bacterium")
    sc.anneal("RT primer", primer, to="RNA", pair=("random primer", "bacterial RNA"),
              label="well-indexed random RT primer",
              unpaired=("PCR handle", "plate barcode", "UMI", "random primer"))
    sc.arrow("RT primer", "in situ reverse transcription")
    sc.labels("RT primer")
    return sc


def droplet_barcode_rows(method: str) -> list[Row]:
    return [Row(chunks=[
        ("cDNA — [plate barcode] — [UMI] — [droplet barcode]", "cbc", False)]),
        Row(chunks=[
            ("plate + droplet indexes jointly identify the bacterium", None, False)])]


def bacterial_rna_library(method: str, name: str):
    plate = feature(method + "_plate", "cell_barcode", "combinatorial",
                    group=method + "_cell", part="plate well")
    droplet = feature(method + "_droplet", "cell_barcode", "combinatorial",
                      group=method + "_cell", part="droplet")
    umi = feature(method + "_umi", "umi", "random")
    return truseq_library([
        seg("plate barcode", "B" * 13, "cbc", placeholder=True, feature=plate),
        seg("UMI", "U" * 8, "umi", placeholder=True, feature=umi),
        seg("bacterial cDNA", "X" * 40, placeholder=True),
    ], name, i5_nt=16, i5_feature=droplet, i5_name="droplet barcode (i5)")



TITLE="BacDrop — droplet bacterial single-cell RNA-seq"
NOTES="01_bacdrop.html"
SOURCE='Defining source: <a href="https://doi.org/10.1016/j.cell.2023.01.002">Ma et al., <i>Cell</i> (2023)</a>.'
SUMMARY="Fixed bacteria receive a plate-specific RT index, are pooled at several cells per 10x droplet, and acquire a second bead barcode after TdT poly(A)-tailing of cDNA."
CAVEAT="The two-index cell identity and 10x bead primer's constant arms are published; individual bead barcode members and the full plate-indexed RT primer table remain unverified because the cited supplementary table is unavailable."
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
  '<b>10x PN-2000210 gel-bead oligo</b> — 5′-AATGATACGGCGACCACCGAGATCTACAC-[16-nt barcode]-TCGTCGGCAGCGTC-3′ (<a href="https://cdn.10xgenomics.com/image/upload/v1666737555/support-documents/CG000209_Chromium_NextGEM_SingleCell_ATAC_ReagentKits_v1.1_UserGuide_RevG.pdf">CG000209 Rev G, Appendix p. 51</a>; the barcode member is a placeholder).',
 ]
INITIAL_ROWS=(Row(chunks=[("fixed bacterial cells containing total RNA",None,False)]),)
INITIAL_NAME="Fixed bacterial cells"

def sections(): return [
 ("Fix, permeabilize and deplete abundant rRNA in cells",[Row(chunks=[("bacterium: RNA — rRNA probes/RNase H — DNase I cleanup",None,False)])],"Universal rRNA probes and RNase H reduce ribosomal molecules before indexing."),
 ("Install the first barcode by random-primed in-situ RT",indexed_random_rt_scene(method="bacdrop").rows(),"Each plate well contributes one RT barcode and each copied molecule receives a UMI."),
 ("Pool cells and TdT-tail their cDNA",[Row(chunks=[("pooled plate-indexed cells → TdT poly(A) tails on cDNA",None,False)])],"The poly(A) tail is added in cells before droplet generation."),
 ("Partition several bacteria per 10x droplet",[Row(chunks=[("many plate indexes + one droplet bead barcode", "cbc",False)])],"Combinatorial plate and droplet indexes permit super-loading."),
 ("Transfer the droplet barcode",droplet_barcode_rows("bacdrop"),"The bead primer copies poly(A)-tailed cDNA during second-strand synthesis in droplets."),
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
 tdt=graph.input("tdt","terminal deoxynucleotidyl transferase",MaterialKind.CATALYST)
 tailed=graph.product("tailed_cdna","poly(A)-tailed indexed cDNA in pooled cells",MaterialKind.MOLECULAR)
 graph.react("tail_cdna","Add poly(A) tails to first-strand cDNA before droplet loading",StepKind.TAILING,
             ((pooled,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(tailed,))
 droplets=graph.input("droplet_condition","10x droplet partitioning condition",MaterialKind.CONDITION)
 occupied=graph.product("occupied_droplets","droplets containing indexed bacteria and a bead",
                        MaterialKind.MOLECULAR)
 empty=graph.product("excluded_droplets","empty or unusable droplets",MaterialKind.WASTE)
 graph.react("partition","Partition pooled bacteria with barcode beads",StepKind.PARTITION,
             ((tailed,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,empty),
             changes_molecule=False)
 graph.discard(empty,"empty and otherwise unusable droplets yield no retained cell library")
 bead_primer=graph.input("bead_primer","commercial droplet-barcode bead primer",
                         MaterialKind.OLIGO,inferred=True,
                         source=SOURCE+"; 10x CG000209 Rev G, Appendix p. 51 (constant arms printed; barcode member unresolved)")
 extension_pol=graph.input("extension_polymerase","second-strand DNA polymerase",
                           MaterialKind.CATALYST)
 barcoded=graph.product("dual_barcoded_cdna","plate- and droplet-barcoded cDNA",
                        MaterialKind.MOLECULAR,state=tuple(droplet_barcode_rows("bacdrop")))
 graph.react("transfer_droplet_barcode","Prime the poly(A) tail and transfer the bead barcode",
             StepKind.EXTENSION,
             ((occupied,InputRole.TEMPLATE),(bead_primer,InputRole.STRUCTURAL),
              (extension_pol,InputRole.CATALYST)),(barcoded,))
 fwd=graph.input("pcr_forward","BacDrop library PCR forward primer",MaterialKind.OLIGO)
 rev=graph.input("pcr_reverse","BacDrop library PCR reverse primer",MaterialKind.OLIGO)
 pcr_pol=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("library_pcr","Enrich cDNA and complete the Illumina library",StepKind.PCR,
             ((barcoded,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),
              (rev,InputRole.STRUCTURAL),(pcr_pol,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The source prints bead-primer constant arms, but the plate RT primer family and exact bead barcode member remain unresolved."))
 graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),
                                       read_lengths=READ_LENGTHS))
 graph.declare_reviewed_scope(("fixed_cells","permeabilizer","permeabilize","rrna_probes",
   "hybridize_rrna","rnase_h","deplete_rrna","indexed_rt_primer","reverse_transcriptase",
   "in_situ_rt","pool_wells","droplet_condition","partition","tdt","tail_cdna",
   "bead_primer","extension_polymerase","transfer_droplet_barcode","pcr_forward",
   "pcr_reverse","pcr_polymerase","library_pcr","library"),SOURCE)
 return graph.require_publishable()
