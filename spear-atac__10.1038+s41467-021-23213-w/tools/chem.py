"""Spear-ATAC accessibility and guide-capture branches."""
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, Segment, feature, oligo, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Highlight, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
TITLE="Spear-ATAC — single-cell ATAC CRISPR screening"
NOTES="01_spear-atac.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41467-021-23213-w">Pierce et al., <i>Nature Communications</i> (2021)</a>.'
SUMMARY="A 10x scATAC-derived workflow captures accessible chromatin and enriches guide identity from the same barcoded nuclei."
CAVEAT=("Supplementary Data 8 prints the five Spear-ATAC primers (oMCB1672, oSP1735, oSP2053, "
        "oSP1594 and the indexed P7 primer) and they are drawn verbatim. The 10x scATAC v1 "
        "bead oligo and ATAC index primers are vendor material named but not printed here, and "
        "the pSP618 vector sequence between the guide barcode and the spacer is described only in "
        "words, so that stretch is a placeholder. Variable cell, guide and sample indexes remain "
        "typed fields, not missing bases.")

CELL=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x scATAC 16-nt whitelist")
GUIDE_BC=feature("guide_barcode","guide_barcode","whitelist",whitelist="Spear-ATAC 10-nt guide-barcode set")
I7=feature("sample_index_i7","sample_index","whitelist",whitelist="selected 8-nt i7 set")

# Supplementary Data 8 (MOESM11) lists the complete primer set; the article Methods
# and Supplementary Note 3 print the same bases.
OSP1594="AATGATACGGCGACCACCGAGA"
OSP2053_ADAPTER=il.TRUSEQ_READ2
OSP2053_SPECIFIC="AAGTATCCCTTGGAGAACCACCTTG"
OSP2053=OSP2053_ADAPTER+OSP2053_SPECIFIC
OSP1735="GCTACATTTTACATGATAGGCTTGG"
OMCB1672="GCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCTTAGCTCTTAAAC"
GUIDE_PRIMER_SITE=OMCB1672
# The published indexed primer stops after the first 22 bases of TruSeq Read 2.
P7_INDEX_TAIL=il.TRUSEQ_READ2[:22]
P7_TRUSEQ_R2=il.P7+"N"*8+P7_INDEX_TAIL
# Downstream of the spacer the library must contain the complement of the
# sgRNA-specific 3' half of oSP2053; those bases are derived from that primer.
GUIDE_DOWNSTREAM=revcomp(OSP2053_SPECIFIC)

ATAC=Construct([
 seg("P5",il.P5,"p5"),seg("cell barcode","N"*16,"cbc",placeholder=True,feature=CELL),
 seg("S5",nx.S5,"s5"),seg("left mosaic end",nx.ME,"me"),
 seg("accessible chromatin fragment","X"*60,placeholder=True),
 seg("right mosaic end reverse complement",nx.ME_RC,"me"),seg("S7 reverse complement",nx.S7_RC,"s7"),
 seg("i7 reverse complement","N"*8,"cbc",placeholder=True,feature=I7),seg("P7 reverse complement",il.P7_RC,"p7")],
 name="Spear-ATAC chromatin library")
AP=(sp.NEXTERA["R1"],sp.NEXTERA["I1"],sp.NEXTERA["I2"],sp.NEXTERA["R2"])

GUIDE=Construct([
 seg("P5",il.P5,"p5"),seg("cell barcode","N"*16,"cbc",placeholder=True,feature=CELL),
 seg("S5",nx.S5,"s5"),seg("left mosaic end",nx.ME,"me"),
 seg("guide barcode","N"*10,"cbc",placeholder=True,feature=GUIDE_BC),
 seg("pSP618 constant region between the guide barcode and the oMCB1672 site","X"*20,placeholder=True),
 seg("custom Read 1 primer site",GUIDE_PRIMER_SITE,"r1"),
 seg("guide spacer","X"*20,"guide",placeholder=True),
 seg("oSP2053 annealing site (complement of its sgRNA-specific 3' half)",GUIDE_DOWNSTREAM),
 seg("TruSeq Read 2 reverse complement",revcomp(il.TRUSEQ_READ2),"r2"),
 seg("i7 reverse complement","N"*8,"cbc",placeholder=True,feature=I7),seg("P7 reverse complement",il.P7_RC,"p7")],
 name="Spear-ATAC guide library")
GUIDE_R1=sp.custom("Read 1","oMCB1672 custom guide Read 1",OMCB1672,"Pierce et al. 2021 Supplementary Note 3")
GUIDE_I1=sp.custom("Index 1 (i7)","guide sample-index primer",revcomp(il.TRUSEQ_READ2),"Pierce et al. 2021 Supplementary Note 3")
GUIDE_I2=sp.NEXTERA["I2"]
GP=(GUIDE_R1,GUIDE_I1,GUIDE_I2)
for library,primers in ((ATAC,AP),(GUIDE,GP)):
 problems=sp.verify(library,primers,required_roles=tuple(p.role for p in primers))
 if problems: raise ValueError(library.name+": "+"; ".join(problems))

FINAL_LIBRARIES=(("Chromatin-accessibility library",ATAC,AP,"Tagmented chromatin fragments carry the droplet cell barcode.","Published Nextera primers read the ATAC branch."),("Guide-identity library",GUIDE,GP,"Guide enrichment preserves the same droplet cell barcode.","Published custom primers read guide spacer, sample index and cell barcode."))

def oligos():
 yield oligo("oSP1594 P5 PCR primer",[seg("P5 PCR",OSP1594,"p5")])
 yield oligo("oSP2053 biotinylated guide linear-PCR primer",[seg("guide enrichment",OSP2053,"r2")],five="5′-/Bio/-",three="-3′")
 yield oligo("oMCB1672 custom guide Read 1 primer",[seg("guide Read 1",OMCB1672,"r1")])
 yield oligo("oSP1735 sgRNA-backbone primer spiked into the 10x GEM master mix",[seg("sgRNA backbone",OSP1735)])
 yield oligo("indexed P7 PCR primer (Supplementary Data 8)",[seg("P7",il.P7,"p7"),seg("i7","N"*8,"cbc",placeholder=True,feature=I7),seg("TruSeq Read 2 first 22 bases",P7_INDEX_TAIL,"r2")])
def sections(): return [
 ("Tagment nuclei before partitioning",[Row(chunks=[("accessible chromatin ** Tn5 mosaic ends","me",False)])],"Bulk tagmentation marks open chromatin."),
 ("Barcode nuclei in droplets",[Row(chunks=[("ATAC fragment + gel-bead cell barcode","cbc",False)])],"Gap fill transfers one droplet identity to accessible fragments."),
 ("Amplify the shared barcoded material",[Row(chunks=[("barcoded ATAC products → PCR pool","r1",False)])],"The pool contains genomic ATAC molecules and guide-associated material."),
 ("Complete the two sequencing libraries",[Row(chunks=[("ATAC library | guide-specific linear PCR → guide library","r2",False)])],"The supplementary protocol prints the terminal primers and completed library structures."),]

def protocol_graph():
 source="Defining source: https://doi.org/10.1038/s41467-021-23213-w"
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
 bead=graph.input("bead_primer","published P5–16-nt-cell-barcode–S5 gel-bead oligo",MaterialKind.OLIGO,
                  source=source)
 opposite=graph.input("opposite_primer","oSP1735 sgRNA-backbone primer spiked into the GEM master mix (1.2 µL of 50 µM)",
                      MaterialKind.OLIGO,source=source,
                      state=Construct((Segment("oSP1735",OSP1735),),name="oSP1735"))
 polymerase=graph.input("polymerase","gap-fill and PCR polymerase",MaterialKind.CATALYST)
 shared=graph.product("shared_pool","shared cell-barcoded ATAC/guide-associated product pool",
                      MaterialKind.MOLECULAR,source=source)
 graph.react("gem_pcr","Gap-fill and amplify the shared barcoded material",StepKind.PCR,
             ((gems,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
              (opposite,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(shared,),
             geometry=PCRGeometry.unavailable(ATAC,
                 reason="The paper defines the shared barcode and branches but not one contiguous pre-branch molecule."))
 atac_f=graph.input("atac_forward","SI-PCR Primer B / P5-side primer",MaterialKind.OLIGO,
                    source=source)
 atac_r=graph.input("atac_reverse","indexed N7xx ATAC primer",MaterialKind.OLIGO,
                    source=source)
 atac_full=graph.product("atac_full_library","finished scATAC-seq library (~150 nM in 20 µL)",
                         MaterialKind.LIBRARY,state=ATAC,source=source)
 graph.react("atac_pcr","Complete the chromatin-accessibility sequencing library",StepKind.PCR,
             ((shared,InputRole.TEMPLATE),(atac_f,InputRole.STRUCTURAL),
              (atac_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(atac_full,),
             geometry=PCRGeometry.exact(
                 ATAC,PCRPrimer("oSP1594 / SI-PCR Primer B",OSP1594),
                 PCRPrimer("indexed N7xx",nx.n7xx_primer("N"*8)),
                 statement="The supplementary protocol prints the P5-side primer and indexed N7xx architecture."))
 atac_lib=graph.product("atac_library","scATAC library aliquot taken to sequencing",
                        MaterialKind.LIBRARY,state=ATAC,source=source)
 guide_in=graph.product("guide_input","2.5 µL scATAC-library aliquot used for sgRNA enrichment",
                        MaterialKind.MOLECULAR,state=ATAC,source=source)
 graph.react("split_library","Take 2.5 µL of the finished scATAC library for sgRNA enrichment",
             StepKind.SPLIT,((atac_full,InputRole.TEMPLATE),),(atac_lib,guide_in))
 guide_linear=graph.input("guide_linear_primer","5′-biotinylated oSP2053 sgRNA-specific primer (0.5 µL of 100 µM per 100 µL)",MaterialKind.OLIGO,source=source,
                          state=Construct((Segment("TruSeq Read 2 adapter",OSP2053_ADAPTER,"r2"),
                                           Segment("sgRNA-specific 3′ half",OSP2053_SPECIFIC)),
                                          name="5′-biotin-oSP2053"))
 linear_enzyme=graph.input("linear_polymerase","2X NEBNext master mix, 25 linear cycles (98/63/72 °C), eight 100 µL reactions",MaterialKind.CATALYST,source=source)
 biotin_strands=graph.product("biotinylated_strands","5′-biotinylated single strands copied from sgRNA-containing library molecules",MaterialKind.MOLECULAR)
 linear_carryover=graph.product("linear_carryover","unextended scATAC library molecules carried through PCR1",MaterialKind.MOLECULAR)
 graph.react("guide_extension","Linearly amplify sgRNA fragments with biotinylated oSP2053 (25 cycles)",StepKind.EXTENSION,
             ((guide_in,InputRole.TEMPLATE),(guide_linear,InputRole.STRUCTURAL),
              (linear_enzyme,InputRole.CATALYST)),(biotin_strands,linear_carryover))
 minelute1=graph.input("minelute_pcr1","pooled 8×100 µL reactions, Qiagen MinElute purification into 40.5 µL",MaterialKind.CONDITION)
 cleaned1=graph.product("pcr1_cleaned","purified PCR1 material",MaterialKind.MOLECULAR)
 cleaned1_waste=graph.product("pcr1_waste","primers and buffer removed by MinElute",MaterialKind.WASTE)
 graph.react("pcr1_cleanup","Pool and MinElute-purify the linear-amplification product",StepKind.CLEANUP,
             ((biotin_strands,InputRole.TEMPLATE),(linear_carryover,InputRole.TEMPLATE),
              (minelute1,InputRole.CONDITION)),(cleaned1,cleaned1_waste))
 graph.discard(cleaned1_waste,"removed with the MinElute flow-through")
 beads=graph.input("streptavidin_beads","Streptavidin MyOne C1 beads in B&W buffer (2 M NaCl), washed twice with B&W and twice with water",MaterialKind.SELECTION,source=source)
 bound=graph.product("bead_bound","bead-bound biotinylated sgRNA fragments",MaterialKind.MOLECULAR)
 unbound=graph.product("unbound_scatac","unbound scATAC library fragments",MaterialKind.WASTE)
 graph.react("streptavidin_capture","Capture the biotinylated sgRNA fragments and wash away scATAC fragments",
             StepKind.AFFINITY_CAPTURE,
             ((cleaned1,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(bound,unbound))
 graph.discard(unbound,"washed off the beads before PCR2; the note warns that skipping this capture contaminates the sgRNA library with the scATAC library")
 guide_f=graph.input("guide_forward","oSP1594 P5 PCR primer",MaterialKind.OLIGO,source=source,
                     state=Construct((Segment("oSP1594",OSP1594,"p5"),),name="oSP1594"))
 guide_r=graph.input("guide_reverse","indexed P7 primer (Supplementary Data 8), unique per sample",MaterialKind.OLIGO,source=source,
                     state=Construct((Segment("P7",il.P7,"p7"),
                                      Segment("eight-base index of choice","N"*8,"cbc",length_bp=8,feature=I7),
                                      Segment("TruSeq Read 2 first 22 bases",P7_INDEX_TAIL,"r2")),
                                     name="indexed P7 primer"))
 guide_lib=graph.product("guide_amplified",GUIDE.name,MaterialKind.LIBRARY,state=GUIDE,
                         source=source)
 graph.react("guide_pcr","Exponentially amplify off the beads with oSP1594 and the indexed P7 primer (15 cycles)",StepKind.PCR,
             ((bound,InputRole.TEMPLATE),(guide_f,InputRole.STRUCTURAL),
              (guide_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(guide_lib,),
             geometry=PCRGeometry.exact(
                 GUIDE,PCRPrimer("oSP1594",OSP1594),
                 PCRPrimer("indexed P7 primer",P7_TRUSEQ_R2),
                 statement="The article Methods and Supplementary Data 8 print oSP1594 and the indexed P7 primer, whose 3′ end is the first 22 bases of TruSeq Read 2."))
 magnet=graph.input("bead_disenrichment","magnet 3 min, keep the supernatant, MinElute into 20.5 µL EB",MaterialKind.CONDITION)
 guide_final=graph.product("guide_library","purified sgRNA sequencing library (supernatant)",
                           MaterialKind.LIBRARY,state=GUIDE,source=source)
 retained_beads=graph.product("retained_beads","beads retaining the biotinylated PCR1 strands",MaterialKind.WASTE)
 graph.react("bead_release","Magnet the beads away and MinElute the amplified supernatant",StepKind.CLEANUP,
             ((guide_lib,InputRole.TEMPLATE),(magnet,InputRole.CONDITION)),(guide_final,retained_beads))
 graph.discard(retained_beads,"the biotinylated template strands stay on the beads and are discarded with them")
 uses={"Read 1":"include","Read 2":"include","Index 1 (i7)":"demultiplex",
       "Index 2 (i5)":"extract"}
 for material,library,primers in ((atac_lib,ATAC,AP),(guide_final,GUIDE,GP)):
  graph.sequence(material,SequencingPlan(
      library=library,primers=primers,acquired_roles=tuple(p.role for p in primers),
      read_uses={p.role:uses[p.role] for p in primers},
      limitation=("Variable barcode and index members are represented by their published lengths and roles. "
                  "The 16-nt droplet cell barcode is acquired as an index read on both libraries and is the "
                  "identifier that links a guide to its accessibility profile; the inspected sources give no "
                  "cycle counts, so no read lengths are asserted.")))
 graph.annotate(Highlight("guide_extension","The sgRNA library is made from the finished scATAC library",
     "Supplementary Note 3 and the Methods start the enrichment from 2.5 µL of the completed ~150 nM scATAC "
     "library, not from a parallel aliquot of the GEM pool. Linear amplification with biotinylated oSP2053 "
     "precedes streptavidin capture; the note states that exponential amplification or omitted capture "
     "contaminates the sgRNA library with the scATAC library."))
 graph.annotate(Highlight("guide_library","The sequenced sgRNA library is the PCR2 supernatant",
     "After PCR2 the beads are magnetised again and the supernatant is kept: the biotinylated PCR1 strands "
     "stay on the beads, so the sequenced molecules are the newly synthesised, indexed copies."))
 graph.annotate(Highlight("guide_amplified","Vector interior between the barcode and the spacer is not transcribed",
     "The Methods describe the pSP618 cassette in words — a mouse U6 promoter carrying a 34 bp Nextera Read 2 "
     "adapter, the Perturb-seq cr1 constant region, a constant region and a 34 bp Nextera Read 1 adapter — but "
     "no inspected original prints those bases. They are drawn as a placeholder; the Addgene pSP618 plasmid "
     "(GenBank MW852482) was not read in this pass. Only the oMCB1672 site, the oSP2053 annealing site and "
     "the two Illumina ends are source-exact."))
 graph.annotate(Highlight("bead_primer","The 10-bp sgRNA barcode is the alternative read-out",
     "Barcodes adjacent to the Nextera Read 1 adapter allow sgRNA identity to be read directly from the scATAC "
     "library, but the authors recommend the targeted amplification shown here and the oMCB1672 custom Read 1 "
     "primer; the barcode window is kept as a typed, unassigned field."))
 graph.declare_reviewed_scope((
     "nuclei","transposome","tn5","tagment","gem_loading","partition","bead_primer",
     "opposite_primer","polymerase","gem_pcr","atac_forward","atac_reverse",
     "atac_pcr","atac_full_library","split_library","atac_library","guide_input",
     "guide_linear_primer","linear_polymerase","guide_extension","minelute_pcr1",
     "pcr1_cleanup","streptavidin_beads","streptavidin_capture","guide_forward",
     "guide_reverse","guide_pcr","bead_disenrichment","bead_release","guide_library"),source)
 return graph.require_publishable()
