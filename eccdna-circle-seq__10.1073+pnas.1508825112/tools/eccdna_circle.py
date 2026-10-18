from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct,Scene,Segment,circle_rows,feature
from protocol_graph import Disposition,Terminal,InputRole,MaterialKind,ProtocolGraph,SequencingPlan,SourceClaim,SourceDiscrepancy,StepKind
from protocol_model import Drawing,molecule,role
import seqprimers as sp

TITLE="eccDNA Circle-Seq — enrichment of native DNA circles"
NOTES="01_eccdna-circle-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1073/pnas.1508825112">Møller et al., <i>PNAS</i> (2015)</a>.'
SUMMARY="Native extrachromosomal circles are column-enriched after alkaline lysis, survive NotI/exonuclease treatment as uncut-circle subsets through130h of Plasmid-Safe exonuclease depletion of linear DNA, are amplified by phi29 (REPLI-g Mini) rolling-circle amplification, then sonicated and indexed on an Apollo 324 robot for single-end Illumina sequencing."
CAVEAT="The protocol enriches pre-existing eccDNA; it does not circularize genomic DNA during library preparation, and the drawn circle is one strand of a representative native double-stranded circle. Spiked pBR322, pUC19 and pUG72 circles travel with the eccDNA as internal detection standards. Phi29 priming sites/sequences remain unassigned in this model; the source names REPLI-g Mini. Møller et al. state only that the Apollo 324 system added adapters and 6-base barcode index labels; adapter bases, joining chemistry, index placement and library PCR remain unassigned in this bounded source scope, so Apollo324 finishing remains an explicit boundary. No canonical TruSeq shell or fabricated RCA template/primer sequence is supplied."
# Native circles pre-exist: no arbitrary DNA sequence or ssDNA closure is supplied.
CIRCLE = Drawing(tuple(circle_rows(Construct([role("native eccDNA unit; sequence unassigned")]), "native circle strand (representative)")))
RCA = molecule("phi29-amplified circular templates", (role("circle-derived repeated DNA; sequence, priming sites and branch/repeat geometry unassigned"),))
FINAL_LIBRARY = Construct([role("Apollo324 library molecule; insert/adapter/index order and end geometry unassigned")], name="eccDNA Circle-Seq Apollo324 endpoint; molecular index position unassigned")
SEQ_PRIMERS = ()
FINAL_CAPTION = "Apollo324 adds adapters and six-base barcode oligos; their molecular placement, joining geometry and library PCR status are unassigned. No generic TruSeq or dA junction is inserted."
SEQUENCING_INTRO = "Source141-nt single-end HiSeq2000 reads are reported; Mapping of Reads also identifies215-nt S2/Z4/R1 reads. Barcode demultiplexing does not establish an acquired index-cycle recipe."
SEQUENCING_ENDING = SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE = "Physical sequencing-primer sequences/sites, index placement and index acquisition cycles remain unassigned."
READ_SPECS = (sp.read_spec("Read 1", "141;215 in S2/Z4/R1", "phi29-amplified circle-derived DNA; source alternative samples", "include"),)


def populations(label, *names):
 scene=Scene()
 for i,name in enumerate(names):
  key="population_"+str(i)
  scene.strand(key,[role(name)],label=label+": "+name)
 return scene


def surviving_circle_populations():
 return Drawing(tuple(CIRCLE.rows() + circle_rows(Construct([role("independent spike-in plasmid circle; sequence unassigned")]), "native spike circle strand (representative)")))


def apollo_reagent_inventory():
 scene=Scene()
 scene.strand("adapters",[role("Apollo324 adapter reagents; molecular sequences unassigned")],label="adapter inventory; no end geometry assigned")
 scene.strand("barcodes",[Segment("source six-base barcode oligos", "N"*6, "cbc", placeholder=True,
              feature=feature("sample_index", "sample_index", "unknown"))],label="source-named six-base barcode inventory; placement unassigned")
 return scene


def protocol_graph():
 graph=ProtocolGraph(TITLE,authority=SOURCE)
 dna=graph.input("cellular_dna","cellular DNA containing native eccDNA circles and linear chromosomes",MaterialKind.BIOLOGICAL,state=populations("cellular DNA","native eccDNA circles","linear chromosomes"),primary=True)
 purification=graph.input("circle_purification","Glass-bead disruption (0.5 mm, 10 min vortex) in L1 solution, then Plasmid Mini AX alkaline treatment, protein precipitation and ion-exchange column (A&A Biotechnology); eluted DNA precipitated into TE",MaterialKind.CONDITION)
 spikes=graph.input("plasmid_spike_ins","pBR322, pUC19 and pUG72 plasmid spike-ins at 1:1, 1:50 and 1:2,500 per cell (assuming 8 × 10^9 genomes), added to the lysate before the column",MaterialKind.MOLECULAR,source=SOURCE+" (Methods, Purification and enrichment of eccDNA; Fig. S5)")
 fraction=graph.product("circle_enriched_fraction","circle-enriched cellular DNA fraction with spiked plasmid circles (6–16 µg)",MaterialKind.MOLECULAR,state=populations("circle-enriched DNA","native eccDNA circles","independent spiked plasmid circles","residual linear chromosomal DNA"))
 graph.react("purify_circles","Lyse 10^10 cells, spike plasmids and purify the circular-DNA-enriched fraction",StepKind.EXTRACTION,((dna,InputRole.TEMPLATE),(spikes,InputRole.SPIKE_IN),(purification,InputRole.CONDITION)),(fraction,))
 noti=graph.input("noti","NotI restriction endonuclease",MaterialKind.CATALYST)
 cut=graph.product("noti_treated","NotI-treated circle-enriched DNA",MaterialKind.MOLECULAR,state=populations("NotI-treated mixture","surviving native eccDNA circles","surviving spiked plasmid circles","free-ended NotI/linear DNA substrates"))
 graph.react("cut_linear","Cut linear chromosomes with NotI (37 °C, 17 h; 80 °C 5 min)",StepKind.DIGESTION,((fraction,InputRole.TEMPLATE),(noti,InputRole.CATALYST)),(cut,))
 exonuclease=graph.input("plasmid_safe_dnase","Plasmid-Safe ATP-dependent DNase",MaterialKind.CATALYST)
 circles=graph.product("native_circles","exonuclease-resistant native eccDNA circles",MaterialKind.MOLECULAR,state=surviving_circle_populations())
 linear_waste=graph.product("linear_dna_waste","hydrolysis products of free-ended DNA; physical removal unassigned",MaterialKind.MOLECULAR)
 graph.react("digest_linear","Digest free-ended DNA with Plasmid-Safe DNase, ATP and enzyme replenished every 24 h for 130 h; 70 °C 30 min",StepKind.DIGESTION,((cut,InputRole.TEMPLATE),(exonuclease,InputRole.CATALYST)),(circles,linear_waste))
 graph.mark_terminal(Terminal(linear_waste.id,Disposition.UNRESOLVED,"Hydrolysis reduces intact free-ended DNA but does not establish physical removal of digestion products from the mixture."))
 used=graph.product("rca_aliquot","1% (5µl) exonuclease-treated sample for REPLI-g",MaterialKind.MOLECULAR,state=circles.state)
 remainder=graph.product("unused_sample","remaining99% exonuclease-treated sample",MaterialKind.MOLECULAR,state=circles.state)
 graph.react("take_rca_aliquot","Take1% (5µl) for phi29 amplification",StepKind.SPLIT,((circles,InputRole.TEMPLATE),),(used,remainder),changes_molecule=False)
 graph.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,"The source uses1% for amplification; subsequent disposition of the remaining sample is unassigned."))
 primers=graph.input("repli_primers","REPLI-g Mini priming reagents; physical oligos/geometry unassigned",MaterialKind.OLIGO)
 phi29=graph.input("phi29","phi29 strand-displacing polymerase",MaterialKind.CATALYST)
 amplified=graph.product("rca_product","rolling-circle-amplified eccDNA concatemers",MaterialKind.MOLECULAR,state=RCA)
 graph.react("rolling_circle","Amplify 1% of surviving circles with phi29 (REPLI-g Mini, 30 °C, 16 h)",StepKind.EXTENSION,((used,InputRole.TEMPLATE),(primers,InputRole.STRUCTURAL),(phi29,InputRole.CATALYST)),(amplified,))
 sonication=graph.input("sonication","Covaris LE220 sonication to ~300 nt",MaterialKind.CONDITION)
 fragments=graph.product("fragments","sonicated phi29 concatemer fragments",MaterialKind.MOLECULAR,state=molecule("sonicated phi29 DNA",(role("circle-derived sonication fragment; exact boundaries unassigned"),)))
 graph.react("fragment","Sonicate the rolling-circle products",StepKind.FRAGMENTATION,((amplified,InputRole.TEMPLATE),(sonication,InputRole.CONDITION)),(fragments,))
 apollo=graph.input("apollo_reagents","Apollo 324 (IntegenX) adapters and 6-base barcode index oligos",MaterialKind.ADAPTER,state=apollo_reagent_inventory(),source=SOURCE+" (Methods, Sequencing: named adapter and six-base barcode inventory; placement unassigned)")
 library=graph.product("sequencing_library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True,source=SOURCE+" (Apollo324 endpoint geometry unassigned)")
 graph.react("apollo_library","Apollo 324 adds adapters and barcode index labels to 300 ng fragmented DNA",StepKind.MAGIC,((fragments,InputRole.TEMPLATE),(apollo,InputRole.STRUCTURAL)),(library,),reason="Methods (Sequencing) state only that 300 ng of sonicated DNA was loaded onto a robotic Apollo 324 system adding adapters and barcode index labels (6-base oligos); adapter sequences, end repair/ligation geometry, enzymes and library PCR are unassigned within this inspected scope; Apollo324 kit documentation was not inspected.")
 graph.annotate(SourceDiscrepancy("sequencing_library","Sequencing depth: 100 Mb versus 40–100 million reads",(SourceClaim("Methods, Sequencing","The 12 samples were sequenced up to 100 Mb as single ends (HiSeq 2000)"),SourceClaim("Results, Circle-Seq method","40–100 million 141-nt single-end reads per population; saturation near 30–45 million reads")),"The model records single-end Read 1 acquisition only and no depth value.","100 Mb of 141-nt reads would be about 0.7 million reads, two orders of magnitude below the read counts in Results and Fig. 2; the Methods unit may mean 100 million reads.","Which figure is a typographical slip is not stated by the authors; depth does not change the molecular construct."))
 graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,acquired_roles=("Read 1",),limitation=SEQUENCING_INTRO+" "+SEQUENCING_UNAVAILABLE))
 graph.declare_reviewed_scope(("cellular_dna","circle_purification","plasmid_spike_ins","purify_circles","noti","cut_linear","plasmid_safe_dnase","digest_linear","take_rca_aliquot","rca_aliquot","unused_sample","repli_primers","phi29","rolling_circle","sonication","fragment","apollo_reagents","apollo_library","sequencing_library"),SOURCE)
 return graph.require_publishable()
