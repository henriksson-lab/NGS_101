"""Illumina two-PCR V3-V4 16S library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg
from chemdraw import Construct, Scene, duplex_rows, revcomp
from targeted_ngs import two_stage_amplicon
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule, role

TITLE="Illumina 16S metagenomic sequencing"
NOTES="01_illumina-16s.html"
SOURCE='Illumina <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf">16S Metagenomic Sequencing Library Preparation, 15044223 B</a>.'
SUMMARY="A V3–V4 locus PCR installs Nextera overhangs; a limited second PCR adds dual sample indexes and complete Illumina flow-cell arms."
CAVEAT="The locus contains the guide's degenerate IUPAC positions, so the schematic uses a target placeholder. The guide prints the Nextera XT N701–N712 and S501–S508 index sequences, but does not assign one pair to a particular sample; the index segments therefore remain variable."
FWD="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG"
REV="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC"
FWD_LOCUS="CCTACGGGNGGCWGCAG"
REV_LOCUS="GACTACHVGGGTATCTAATCC"
FWD_OVERHANG=FWD[:-len(FWD_LOCUS)]
REV_OVERHANG=REV[:-len(REV_LOCUS)]

PCR1_PRODUCT=Construct([
 seg("forward Nextera overhang",FWD_OVERHANG,"s5"),
 seg("forward V3 locus primer",FWD_LOCUS,"r1",placeholder=True),
 seg("16S V3–V4 target","N"*40,placeholder=True),
 seg("reverse V4 locus primer complement",revcomp(REV_LOCUS),"r2",placeholder=True),
 seg("reverse Nextera overhang complement",revcomp(REV_OVERHANG),"s7"),
],name="tailed V3–V4 PCR1 product")

FINAL_LIBRARY,SEQ_PRIMERS=two_stage_amplicon("16S V3–V4 amplicon",40)
FINAL_CAPTION="Dual-indexed V3–V4 amplicon after the second PCR."
SEQUENCING_INTRO="All four landing sites come from the Nextera overhang/index-PCR architecture."


def pcr1_scene():
    scene=Scene.duplex(list(PCR1_PRODUCT),label="PCR1 product")
    scene.labels("top")
    return scene

def sections():
    return [
      ("PCR1 — select V3–V4 and install universal tails",pcr1_scene().rows(),
       "The locus-specific 3′ ends determine the amplicon; their 5′ overhangs become PCR2 priming sites."),
      ("PCR2 — add dual indexes and flow-cell arms",
       duplex_rows(FINAL_LIBRARY,label="PCR2 product"),
       "Limited-cycle index PCR turns the cleaned PCR1 product into a sequencing-ready library."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("community_dna","community genomic DNA",MaterialKind.BIOLOGICAL,
                    state=molecule("community DNA",(role("bacterial 16S loci"),)),primary=True)
    fwd=graph.input("v3_forward","V3 locus primer with Nextera overhang",MaterialKind.OLIGO)
    rev=graph.input("v4_reverse","V4 locus primer with Nextera overhang",MaterialKind.OLIGO)
    pol1=graph.input("pcr1_polymerase","locus-PCR polymerase",MaterialKind.CATALYST)
    product=graph.product("pcr1_product",PCR1_PRODUCT.name,MaterialKind.MOLECULAR,state=PCR1_PRODUCT)
    geometry=PCRGeometry.exact(PCR1_PRODUCT,PCRPrimer("V3 forward",FWD),PCRPrimer("V4 reverse",REV),statement="The Illumina guide prints both complete V3–V4 primers including their Nextera overhangs.")
    graph.react("pcr1","Select V3–V4 and install universal tails",StepKind.PCR,
                ((dna,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(pol1,InputRole.CATALYST)),(product,),geometry=geometry)
    i5=graph.input("i5_primer","Nextera i5 index-PCR primer",MaterialKind.OLIGO)
    i7=graph.input("i7_primer","Nextera i7 index-PCR primer",MaterialKind.OLIGO)
    pol2=graph.input("pcr2_polymerase","index-PCR polymerase",MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    graph.react("pcr2","Add dual indexes and flow-cell arms",StepKind.PCR,
                ((product,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),(i7,InputRole.STRUCTURAL),(pol2,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason="The guide lists every N701–N712 and S501–S508 index sequence but does not select a sample pair; final sequencing-primer sites are checked on the completed library."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope(["community_dna","v3_forward","v4_reverse","pcr1_polymerase","pcr1","i5_primer","i7_primer","pcr2_polymerase","pcr2","library"],SOURCE)
    return graph.require_publishable()
