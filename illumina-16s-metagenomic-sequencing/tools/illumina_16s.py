"""Protocol-owned Illumina 16S V3–V4 two-PCR declaration (15044223 Rev. B)."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, feature, iupac_complement, oligo, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
from protocol_model import role

TITLE="Illumina 16S metagenomic sequencing"
NOTES="01_illumina-16s.html"
SOURCE='Illumina <a href="https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf">16S Metagenomic Sequencing Library Preparation, 15044223 B</a>.'
SUMMARY="A V3–V4 locus PCR installs universal overhangs. After cleanup and aliquoting, index PCR extends that actual product, preserving both gene-specific primer regions. A second cleanup precedes MiSeq v3 paired-end sequencing; there is no Tn5 reaction."
CAVEAT=("The native community DNA and biological interval have no selected sequence or fixed length. "
        "The defining guide prints both full locus primers and the eight-base N701–N712/S501–S508 index families. "
        "No experimental sample or index pair is selected here. INFERRED — complete index-PCR primers and flow-cell "
        "caps use a locally declared canonical Nextera-compatible extension; a matched historical full-primer "
        "original remains unassigned. The source-known PCR1 product and its two locus-primer regions are retained. "
        "MiSeq Index 2 uses grafted P5 and seven non-imaged chemistry cycles; the displayed P5 reference is an "
        "effective post-extension read-start anchor, not an exact physical grafted oligo. The manufacturer "
        "mechanism is checked in Indexed Sequencing Overview15057455 v08 (November 2020), not assigned as the "
        "authors’ historical instrument-recipe revision.")
FWD="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG"
REV="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC"
FWD_LOCUS="CCTACGGGNGGCWGCAG"
REV_LOCUS="GACTACHVGGGTATCTAATCC"
FWD_OVERHANG=FWD[:-len(FWD_LOCUS)]
REV_OVERHANG=REV[:-len(REV_LOCUS)]
I5=feature("sixteen_s_i5","sample_index","whitelist",
           whitelist="Illumina15044223B S501–S508; member unselected")
I7=feature("sixteen_s_i7","sample_index","whitelist",
           whitelist="Illumina15044223B N701–N712; member unselected")


def locus(name: str, sequence: str, tag: str) -> Segment:
    return Segment(name, sequence, tag, placeholder=True, length_bp=len(sequence),
                   bottom=iupac_complement(sequence))


def forward_primer() -> list[Segment]:
    return [Segment("forward universal overhang",FWD_OVERHANG,"s5"),
            locus("forward V3 locus primer",FWD_LOCUS,"r1")]


def reverse_primer() -> list[Segment]:
    return [Segment("reverse universal overhang",REV_OVERHANG,"s7"),
            locus("reverse V4 locus primer",REV_LOCUS,"r2")]


PCR1_PRODUCT=Construct([
 *forward_primer(),
 role("biological 16S interval between the two primer regions"),
 locus("reverse V4 locus primer complement",iupac_complement(REV_LOCUS)[::-1],"r2"),
 Segment("reverse universal overhang complement",revcomp(REV_OVERHANG),"s7"),
],name="tailed V3–V4 PCR1 product")
GENOMIC_TEMPLATE=Construct([role("unassigned community genomic 16S template")],
                           name="native community DNA")


def index_primer(side: str, *, inferred: bool = True) -> list[Segment]:
    """Canonical whole-primer reconstruction; source-selected kit generation unassigned."""
    if side == "i5":
        return [Segment("P5",il.P5,"p5",inferred=inferred),
                Segment("i5", "J"*8,"cbc",placeholder=True,feature=I5,inferred=inferred),
                Segment("S5 annealing arm",nx.S5,"s5",inferred=inferred)]
    if side == "i7":
        return [Segment("P7",il.P7,"p7",inferred=inferred),
                Segment("i7", "I"*8,"cbc",placeholder=True,feature=I7,inferred=inferred),
                Segment("S7 annealing arm",nx.S7,"s7",inferred=inferred)]
    raise ValueError("index primer must be i5 or i7")


def index_template(side: str) -> str:
    return "".join(s.top for s in index_primer(side)).replace("J","N").replace("I","N")


def final_library(*, inferred: bool = True) -> Construct:
    """Extend the declared PCR1 product; retain its source-known molecular regions."""
    return Construct([
        Segment("P5",il.P5,"p5",inferred=inferred),
        Segment("i5","J"*8,"cbc",placeholder=True,feature=I5,inferred=inferred),
        *list(PCR1_PRODUCT),
        Segment("i7 reverse complement","I"*8,"cbc",placeholder=True,feature=I7,inferred=inferred),
        Segment("P7 reverse complement",il.P7_RC,"p7",inferred=inferred),
    ],name="dual-indexed V3–V4 amplicon")


FINAL_LIBRARY=final_library()
SEQ_PRIMERS=(sp.NEXTERA["R1"],sp.NEXTERA["I1"],
             sp.custom("Index 2 (i5)","effective P5 anchor after seven dark cycles",il.P5,
                       "Illumina Indexed Sequencing Overview15057455v08 p4",
                       "Computational post-extension read-start anchor; not the literal physical grafted primer."),
             sp.NEXTERA["R2"])
READ_LENGTHS={"Read 1":300,"Index 1 (i7)":8,"Index 2 (i5)":8,"Read 2":300}
ACQUIRED_ROLES=tuple(READ_LENGTHS)
FINAL_CAPTION="INFERRED — canonical index-PCR caps extend the actual PCR1 product. Both source locus-primer regions remain; biological interval length is unassigned."
SEQUENCING_INTRO=("The MiSeq v3 300/8/8/300 view starts Read 1 with the 17-base forward locus-primer region "
                  "and Read 2 with the 21-base reverse locus-primer region. The internal biological interval "
                  "has variable length; no fixed readthrough into the opposite end is asserted. Index 2 "
                  "has seven additional non-imaged chemistry cycles before its eight reported bases.")


def oligos():
    yield oligo("Full source forward amplicon primer",forward_primer())
    yield oligo("Full source reverse amplicon primer",reverse_primer())
    yield oligo("INFERRED — canonical S5XX index primer; member unselected",index_primer("i5"))
    yield oligo("INFERRED — canonical N7XX index primer; member unselected",index_primer("i7"))


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("community_dna","community genomic DNA",MaterialKind.BIOLOGICAL,
                    state=GENOMIC_TEMPLATE,primary=True)
    fwd=graph.input("v3_forward","full source V3 forward primer",MaterialKind.OLIGO,
                    state=FWD)
    rev=graph.input("v4_reverse","full source V4 reverse primer",MaterialKind.OLIGO,
                    state=REV)
    pol1=graph.input("pcr1_polymerase","KAPA HiFi HotStart ReadyMix",MaterialKind.CATALYST)
    product=graph.product("pcr1_product",PCR1_PRODUCT.name,MaterialKind.MOLECULAR,state=PCR1_PRODUCT)
    graph.react("pcr1","Amplify V3–V4 and install overhangs (25 cycles)",StepKind.PCR,
                ((dna,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),
                 (pol1,InputRole.CATALYST)),(product,),
                geometry=PCRGeometry.unavailable(GENOMIC_TEMPLATE,
                  reason="Both full primers are printed, but this declaration does not select a native community template or its binding coordinates. The PCR1 product is not substituted as the genomic input."))
    beads1=graph.input("cleanup1_beads","AMPure XP beads for PCR1 cleanup",MaterialKind.SELECTION)
    clean1=graph.product("clean_pcr1","50 µl purified PCR1 product",MaterialKind.MOLECULAR,state=PCR1_PRODUCT)
    waste1=graph.product("cleanup1_waste","material not recovered by cleanup1",MaterialKind.WASTE)
    graph.react("cleanup1","Purify the amplicon PCR product",StepKind.CLEANUP,
                ((product,InputRole.TEMPLATE),(beads1,InputRole.SELECTION)),(clean1,waste1),changes_molecule=False)
    graph.discard(waste1,"Cleanup wash/excluded material; exact molecular composition unassigned.")
    aliquot=graph.product("pcr2_aliquot","5 µl PCR1 aliquot for index PCR",MaterialKind.MOLECULAR,state=PCR1_PRODUCT)
    remainder=graph.product("pcr1_remainder","45 µl PCR1 remainder not used in this workflow",MaterialKind.MOLECULAR,state=PCR1_PRODUCT)
    transfer=graph.input("transfer_volume","source-specified 5 µl of 50 µl transfer",MaterialKind.CONDITION)
    graph.react("aliquot","Transfer 5 µl; retain the 45 µl remainder for possible other use",StepKind.SAMPLING,
                ((clean1,InputRole.TEMPLATE),(transfer,InputRole.CONDITION)),(aliquot,remainder),changes_molecule=False)
    graph.mark_terminal(Terminal(remainder.id,Disposition.STORED,
        "The guide states the remaining 45 µl is not used in this protocol and can be stored for other uses; storage here is optional, not a molecular rejection."))
    i5=graph.input("i5_primer","INFERRED canonical S5XX index primer",MaterialKind.OLIGO,
                   state=tuple(index_primer("i5")),inferred=True)
    i7=graph.input("i7_primer","INFERRED canonical N7XX index primer",MaterialKind.OLIGO,
                   state=tuple(index_primer("i7")),inferred=True)
    pol2=graph.input("pcr2_polymerase","KAPA HiFi HotStart ReadyMix",MaterialKind.CATALYST)
    amplified=graph.product("indexed_amplicon",FINAL_LIBRARY.name,MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    geometry=PCRGeometry.exact(PCR1_PRODUCT,
        PCRPrimer("inferred S5XX full primer",index_template("i5"),anneal_nt=len(nx.S5)),
        PCRPrimer("inferred N7XX full primer",index_template("i7"),anneal_nt=len(nx.S7)),
        statement="Conditional canonical whole-primer model: S5/S7 terminal 14/15-base arms bind the actual PCR1 overhangs. Full historical primer/cap identity and a sample pair are unassigned; this is not their source verification.")
    graph.react("pcr2","INFERRED — extend the actual PCR1 product with dual indexes and caps (8 cycles)",StepKind.PCR,
                ((aliquot,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),(i7,InputRole.STRUCTURAL),
                 (pol2,InputRole.CATALYST)),(amplified,),geometry=geometry)
    beads2=graph.input("cleanup2_beads","AMPure XP beads for index-PCR cleanup",MaterialKind.SELECTION)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    waste2=graph.product("cleanup2_waste","material not recovered by cleanup2",MaterialKind.WASTE)
    graph.react("cleanup2","Purify the indexed library before quantification/pooling",StepKind.CLEANUP,
                ((amplified,InputRole.TEMPLATE),(beads2,InputRole.SELECTION)),(library,waste2),changes_molecule=False)
    graph.discard(waste2,"Cleanup wash/excluded material; exact molecular composition unassigned.")
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
                   acquired_roles=ACQUIRED_ROLES,read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope(("community_dna","v3_forward","v4_reverse","pcr1_polymerase","pcr1",
        "cleanup1_beads","cleanup1","transfer_volume","aliquot","i5_primer","i7_primer","pcr2_polymerase","pcr2",
        "cleanup2_beads","cleanup2","library"),SOURCE)
    return graph.require_publishable()
