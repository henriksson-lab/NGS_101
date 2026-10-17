from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Construct,Segment
from circular import ClosedCircle,rolling_circle
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,SourceClaim,SourceDiscrepancy,StepKind
from protocol_model import molecule,role

TITLE="eccDNA Circle-Seq — enrichment of native DNA circles"
NOTES="01_eccdna-circle-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1073/pnas.1508825112">Møller et al., <i>PNAS</i> (2015)</a>.'
SUMMARY="Native extrachromosomal circles are column-enriched after alkaline lysis, survive NotI cutting and 130 h of Plasmid-Safe exonuclease depletion of linear DNA, are amplified by phi29 (REPLI-g Mini) rolling-circle amplification, then sonicated and indexed on an Apollo 324 robot for single-end Illumina sequencing."
CAVEAT="The protocol enriches pre-existing eccDNA; it does not circularize genomic DNA during library preparation, and the drawn circle is one strand of a representative native double-stranded circle. Spiked pBR322, pUC19 and pUG72 circles travel with the eccDNA as internal detection standards. The phi29 priming site is representative: the paper cites the REPLI-g Mini protocol without describing its primers. Møller et al. state only that the Apollo 324 system added adapters and 6-base barcode index labels; the adapter bases, the joining chemistry and any library PCR are not described, so that step is a warned boundary and the standard single-index Illumina arms (6-nt index window) are shown as inferred."
UNIT=Construct([Segment("native eccDNA unit","ACGTTGCAACGATCGT",placeholder=False)],name="representative eccDNA (one strand of a native double-stranded circle)")
# Native circles pre-exist in the cells (Results, Fig. 1B): no ligation step closes them, so the
# circle is declared closed directly rather than built by ssDNA circularisation.
CIRCLE=ClosedCircle(UNIT,"native eccDNA unit","native eccDNA unit")
# One representative priming site; the paper names only "phi29 polymerase ... REPLI-g Mini Kit
# protocol" and does not describe the kit's primers.
RCA=rolling_circle(CIRCLE,"TGCAACGT",copies=3)
FINAL_LIBRARY,SEQ_PRIMERS=truseq_library([seg("phi29-amplified eccDNA fragment","X"*48,placeholder=True)],"eccDNA Circle-Seq library",dual_index=False,inferred_adapters=True,i7_nt=6)  # Methods "Sequencing": "barcode index labels (6-base oligos)"; index position inferred
FINAL_CAPTION="INFERRED — sonicated phi29 product after Apollo 324 adapter and 6-base index addition; adapter bases, junction and index position are not printed and are drawn as a standard single-index Illumina library."
SEQUENCING_INTRO="The study used single-end HiSeq 2000 sequencing; Results report 141-nt reads, and Methods (Mapping of Reads) add that the longest reads (215 nt, samples S2, Z4 and R1) were adapter-trimmed. Reads were demultiplexed by index barcodes, but how the 6-base index was read is not stated. The inferred library also shows the opposite standard primer site without implying that Read 2 was collected."

def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 dna=graph.input("cellular_dna","cellular DNA containing native eccDNA circles and linear chromosomes",MaterialKind.BIOLOGICAL,state=molecule("cellular DNA",(role("native eccDNA circles"),role("linear chromosomes"))),primary=True)
 purification=graph.input("circle_purification","Glass-bead disruption (0.5 mm, 10 min vortex) in L1 solution, then Plasmid Mini AX alkaline treatment, protein precipitation and ion-exchange column (A&A Biotechnology); eluted DNA precipitated into TE",MaterialKind.CONDITION)
 spikes=graph.input("plasmid_spike_ins","pBR322, pUC19 and pUG72 plasmid spike-ins at 1:1, 1:50 and 1:2,500 per cell (assuming 8 × 10^9 genomes), added to the lysate before the column",MaterialKind.MOLECULAR,source=SOURCE+" (Methods, Purification and enrichment of eccDNA; Fig. S5)")
 fraction=graph.product("circle_enriched_fraction","circle-enriched cellular DNA fraction with spiked plasmid circles (6–16 µg)",MaterialKind.MOLECULAR,state=molecule("circle-enriched DNA",(role("native eccDNA circles"),role("spiked plasmid circles"),role("residual linear chromosomal DNA"))))
 graph.react("purify_circles","Lyse 10^10 cells, spike plasmids and purify the circular-DNA-enriched fraction",StepKind.EXTRACTION,((dna,InputRole.TEMPLATE),(spikes,InputRole.SPIKE_IN),(purification,InputRole.CONDITION)),(fraction,))
 noti=graph.input("noti","NotI restriction endonuclease",MaterialKind.CATALYST)
 cut=graph.product("noti_treated","NotI-treated circle-enriched DNA",MaterialKind.MOLECULAR,state=fraction.state)
 graph.react("cut_linear","Cut linear chromosomes with NotI (37 °C, 17 h; 80 °C 5 min)",StepKind.DIGESTION,((fraction,InputRole.TEMPLATE),(noti,InputRole.CATALYST)),(cut,))
 exonuclease=graph.input("plasmid_safe_dnase","Plasmid-Safe ATP-dependent DNase",MaterialKind.CATALYST)
 circles=graph.product("native_circles","exonuclease-resistant native eccDNA circles",MaterialKind.MOLECULAR,state=CIRCLE)
 linear_waste=graph.product("linear_dna_waste","degraded free-ended chromosomal DNA",MaterialKind.WASTE)
 graph.react("digest_linear","Digest free-ended DNA with Plasmid-Safe DNase, ATP and enzyme replenished every 24 h for 130 h; 70 °C 30 min",StepKind.DIGESTION,((cut,InputRole.TEMPLATE),(exonuclease,InputRole.CATALYST)),(circles,linear_waste))
 graph.discard(linear_waste,"has free ends and is depleted during native-circle enrichment")
 phi29=graph.input("phi29","phi29 strand-displacing polymerase",MaterialKind.CATALYST)
 amplified=graph.product("rca_product","rolling-circle-amplified eccDNA concatemers",MaterialKind.MOLECULAR,state=RCA)
 graph.react("rolling_circle","Amplify 1% of surviving circles with phi29 (REPLI-g Mini, 30 °C, 16 h)",StepKind.EXTENSION,((circles,InputRole.TEMPLATE),(phi29,InputRole.CATALYST)),(amplified,))
 sonication=graph.input("sonication","Covaris LE220 sonication to ~300 nt",MaterialKind.CONDITION)
 fragments=graph.product("fragments","sonicated phi29 concatemer fragments",MaterialKind.MOLECULAR,state=amplified.state)
 graph.react("fragment","Sonicate the rolling-circle products",StepKind.FRAGMENTATION,((amplified,InputRole.TEMPLATE),(sonication,InputRole.CONDITION)),(fragments,))
 apollo=graph.input("apollo_reagents","Apollo 324 (IntegenX) adapters and 6-base barcode index oligos",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (Methods, Sequencing: reagents named, bases not printed)")
 library=graph.product("sequencing_library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (compatible single-index Illumina shell inferred)")
 graph.react("apollo_library","Apollo 324 adds adapters and barcode index labels to 300 ng fragmented DNA",StepKind.MAGIC,((fragments,InputRole.TEMPLATE),(apollo,InputRole.STRUCTURAL)),(library,),reason="Methods (Sequencing) state only that 300 ng of sonicated DNA was loaded onto a robotic Apollo 324 system adding adapters and barcode index labels (6-base oligos); adapter sequences, end repair/ligation geometry, enzymes and any library PCR are not described, and the Apollo 324 kit documentation was not inspected.")
 graph.annotate(SourceDiscrepancy("sequencing_library","Sequencing depth: 100 Mb versus 40–100 million reads",(SourceClaim("Methods, Sequencing","The 12 samples were sequenced up to 100 Mb as single ends (HiSeq 2000)"),SourceClaim("Results, Circle-Seq method","40–100 million 141-nt single-end reads per population; saturation near 30–45 million reads")),"The model records single-end Read 1 acquisition only and no depth value.","100 Mb of 141-nt reads would be about 0.7 million reads, two orders of magnitude below the read counts in Results and Fig. 2; the Methods unit may mean 100 million reads.","Which figure is a typographical slip is not stated by the authors; depth does not change the molecular construct."))
 graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS),acquired_roles=("Read 1",)))
 graph.declare_reviewed_scope(("cellular_dna","circle_purification","plasmid_spike_ins","purify_circles","noti","cut_linear","plasmid_safe_dnase","digest_linear","phi29","rolling_circle","sonication","fragment","apollo_reagents","apollo_library","sequencing_library"),SOURCE)
 return graph.require_publishable()
