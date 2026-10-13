"""Microbe-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import nextera_library,seg
from chemdraw import Row,feature
from microbial_singlecell import indexed_mda_rows
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE="Microbe-seq — droplet single-microbe genome sequencing"
NOTES="01_microbe-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1126/science.abm1483">Zheng et al., <i>Science</i> (2022)</a>.'
SUMMARY="Individual microbes are lysed and MDA-amplified in droplets, tagmented by reagent fusion, then fused with barcode-bead droplets for indexed PCR."
CAVEAT="The paper establishes droplet operations and library order but not every commercial bead base; those regions are structural placeholders."
CELL=feature("microbeseq_droplet","cell_barcode","whitelist",whitelist="Microbe-seq barcode bead set")
FINAL_LIBRARY,SEQ_PRIMERS=nextera_library([seg("microbe barcode","B"*16,"cbc",placeholder=True,feature=CELL),seg("single-microbe genomic DNA","X"*46,placeholder=True)],"Microbe-seq WGS library")
FINAL_CAPTION="A bead barcode identifies all tagmented MDA fragments derived from one microbial droplet."
SEQUENCING_INTRO="Paired genomic reads are associated with the bead-derived microbe barcode and sample indexes."
READ_LENGTHS={"Read 1":150,"Read 2":150}
SHOW_READ_CONTENT=True

def sections(): return [
 ("Encapsulate and lyse individual microbes",[Row(chunks=[("one microbe + lysis droplet",None,False)])],"Lysis remains physically partitioned by cell."),
 ("Fuse with amplification droplets for MDA",indexed_mda_rows()[:1],"Phi29 multiple-displacement amplification creates enough genomic DNA while retaining droplet identity."),
 ("Fuse with Tn5 droplets for tagmentation",indexed_mda_rows()[1:2],"Nextera mosaic ends are transferred to amplified genomic fragments inside droplets."),
 ("Fuse with barcode-bead/PCR droplets",indexed_mda_rows()[2:],"Bead primers install one droplet barcode on the tagmented fragments."),
 ("Break emulsion and complete sequencing adapters",[Row(chunks=[("barcoded genomic amplicons → sample-index PCR → Illumina library",None,False)])],"Bulk PCR completes the standard sequencing arms."),
]

def protocol_graph():
 source="Defining source: https://doi.org/10.1126/science.abm1483"
 graph=ProtocolGraph("Microbe-seq",authority=source)
 microbes=graph.input("microbes","microbial suspension",MaterialKind.BIOLOGICAL,primary=True)
 loading=graph.input("droplet_loading","single-microbe droplet loading condition",MaterialKind.CONDITION)
 droplets=graph.product("single_microbe_droplets","droplets containing individual microbes",MaterialKind.BIOLOGICAL)
 excluded=graph.product("excluded_droplets","empty or multiplet droplets",MaterialKind.WASTE)
 graph.react("partition","Partition microbes into lysis droplets",StepKind.PARTITION,
             ((microbes,InputRole.TEMPLATE),(loading,InputRole.CONDITION)),(droplets,excluded))
 graph.discard(excluded,"empty and multiplet droplets do not yield accepted single-microbe libraries")
 lysis=graph.input("lysis_mix","droplet lysis reagents",MaterialKind.CONDITION)
 dna=graph.product("released_genomes","lysed single-microbe genomic DNA",MaterialKind.MOLECULAR)
 debris=graph.product("cell_debris","insoluble microbial debris",MaterialKind.WASTE)
 graph.react("lysis","Lyse microbes inside droplets",StepKind.SOLUBILIZATION,
             ((droplets,InputRole.TEMPLATE),(lysis,InputRole.CONDITION)),(dna,debris))
 graph.discard(debris,"insoluble debris does not contribute to the sequencing library")
 random=graph.input("mda_primers","MDA random primers",MaterialKind.OLIGO,
                    inferred=True,source=source+" (commercial MDA primer bases unavailable)")
 phi29=graph.input("phi29","phi29 DNA polymerase",MaterialKind.CATALYST)
 mda=graph.product("mda_products","single-microbe MDA products",MaterialKind.MOLECULAR)
 graph.react("mda","Amplify each released genome by multiple displacement",StepKind.EXTENSION,
             ((dna,InputRole.TEMPLATE),(random,InputRole.STRUCTURAL),(phi29,InputRole.CATALYST)),(mda,))
 transposome=graph.input("transposome","Nextera-loaded Tn5 transposome",MaterialKind.ADAPTER)
 tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=graph.product("tagmented_dna","tagmented MDA fragments",MaterialKind.MOLECULAR)
 twaste=graph.product("tagmentation_waste","non-library tagmentation products",MaterialKind.WASTE)
 graph.react("tagment","Fuse Tn5 droplets and tagment MDA products",StepKind.TAGMENTATION,
             ((mda,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
             (tagged,twaste))
 graph.discard(twaste,"fragments lacking two usable library ends are not amplified")
 bead=graph.input("barcode_bead_primer","barcode-bead PCR primer",MaterialKind.OLIGO,
                  inferred=True,source=source+" (commercial bead sequence is only partly disclosed)")
 opposite=graph.input("opposite_primer","opposite Nextera PCR primer",MaterialKind.OLIGO)
 polymerase=graph.input("droplet_polymerase","droplet PCR polymerase",MaterialKind.CATALYST)
 barcoded=graph.product("barcoded_fragments","bead-barcoded genomic amplicons",MaterialKind.MOLECULAR)
 graph.react("barcode_pcr","Fuse barcode-bead droplets and copy the microbe barcode",StepKind.PCR,
             ((tagged,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
              (opposite,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(barcoded,),
             geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                 reason="The paper does not disclose every commercial bead-primer base."))
 i5=graph.input("index_primer_i5","sample-index PCR i5 primer",MaterialKind.OLIGO)
 i7=graph.input("index_primer_i7","sample-index PCR i7 primer",MaterialKind.OLIGO)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 graph.react("index_pcr","Break the emulsion and complete sequencing adapters",StepKind.PCR,
             ((barcoded,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
              (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                 reason="The pre-index-PCR intermediate is not reported base-for-base; final Nextera arms are canonical."))
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),
     acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS,
     read_uses={"Read 1":"include","Index 1 (i7)":"demultiplex",
                "Index 2 (i5)":"demultiplex","Read 2":"include"},
     limitation="Commercial bead-linker bases are not fully disclosed."))
 graph.declare_reviewed_scope((
     "microbes","droplet_loading","partition","lysis_mix","lysis","mda_primers",
     "phi29","mda","transposome","tn5","tagment","barcode_bead_primer",
     "opposite_primer","droplet_polymerase","barcode_pcr","index_primer_i5",
     "index_primer_i7","index_pcr","library"),source)
 return graph.require_publishable()
