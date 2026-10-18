"""Protocol-owned native declaration; source boundaries stay local."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Scene,Segment,feature,revcomp
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as M,InputRole as I,StepKind as K,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry
import crispr as cr
TITLE='lentiCRISPR v1 screening'
NOTES='01_lenticrispr-v1.html'
SOURCE='Shalem et al., <a href="https://doi.org/10.1126/science.1247005">Science2014, published manuscriptPMC4089965 Results53–69</a>; depositedAddgene49535 map for invariantcassette components; July2014GeCKO plasmid-primer workbook as a separately scoped comparator.'
SUMMARY='The graph follows the original v1 pooled-screen lineage. Native vector DNA, particle-associated RNA and integrated genomic DNA are distinct molecular states; the original screen readout remains bounded by its unread SupplementaryMethods.'
CAVEAT='The gated original SupplementaryInformation remains an unresolved source obligation. The later Joung one-PCR library is not substituted for this experiment. The July2014 workbook explicitly describes plasmid readout and refers genomicDNA preparation back to the SupplementaryMethods; its inlinebarcodes do not establish the original screen run. Representative inserted20ntguide/extraG geometry is INFERRED from the depositedempty-vector components; guide identity, exactclonedjunctions, integration sites and native readout ends remain unassigned.'
FINAL_LIBRARY=None
SEQ_PRIMERS=()
SEQUENCING_ENDING='The original screen was deep-sequenced, but native library primer strands, barcode scheme, loaded sequencing primers and acquired run lengths remain unverified here. A separate later compatible readout has its own <a href="../joung-one-pcr-crispr-screen-readout__10.1038+nprot.2017.016/joung-one-pcr-crispr-screen-readout.html">Joung one-PCR page</a>.'
GUIDE=feature('guide','guide_barcode','unknown')
# Only these invariant hU6/scaffold components are reused; deposited original
# confirms both. The map's filler is not an experimentally cloned guide member.
U6=cr.U6_UPSTREAM+cr.U6_3PRIME
SCAFFOLD=cr.SCAFFOLD_V1

def guide_inventory(sc):
 sc.strand('guide_inventory',[role('source guide-identity inventory; native width/strand placement unassigned',molecular_feature=GUIDE)],label='unplaced identifier inventory; not another physical polymer')

def dna_state(stage):
 sc=Scene()
 if stage=='empty':
  sc.strand('vector',[role('native circular empty vector context'),Segment('hU6 component',U6),role('deposited cloning filler; native cut/end chemistry unassigned'),Segment('scaffold component',SCAFFOLD),role('cPPT/CTS downstream of sourcecassette; remainingvectorcontext')],label='partial empty-vector molecular view; no inserted guide yet')
 elif stage in ('plasmid','genome','screen'):
  prefix=[role('native genomic leftflank; integration sites/copy number unassigned')] if stage!='plasmid' else [role('native circular vector leftcontext')]
  suffix=[role('native genomic rightflank')] if stage!='plasmid' else [role('cPPT/CTS downstreamofcassette; native circular rightcontext')]
  parts=prefix+[Segment('hU6 component',U6),Segment('representative extraG','G',inferred=True),Segment('representative20ntguide','N'*20,placeholder=True,inferred=True,feature=GUIDE),Segment('scaffold component',SCAFFOLD)]+suffix
  sc.strand('DNA',parts,label='INFERRED representative inserted-guide locus; integration/circular context not completed')
  sc.strand('other',[role('possible independent unreacted/original DNA, assemblyinsert/filler, oligos, RNA andcellular polymers; class recovery unassigned')])
 else:
  sc.strand('DNA',[role(stage+' native DNA population; full sequence and molecular ends unassigned')]);sc.strand('other',[role('possible independent co-carried DNA/RNA/oligos; class recovery unassigned')])
 if stage not in ('empty','plasmid','genome','screen'):guide_inventory(sc)
 return sc

def bounded(g,x,why):g.mark_terminal(Terminal(x.id,Disposition.UNRESOLVED,why))

def oligos():return ()

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 vector=g.input('vector','empty lentiCRISPRv1 vector; deposited49535 molecular components',M.MOLECULAR,state=dna_state('empty'),primary=True)
 poolstate=Scene();poolstate.strand('oligos',[role('independent synthesized guide-pool members; full source strands and pooled assembly chemistry unreviewed')])
 guide_inventory(poolstate)
 pool=g.input('guide_pool','source GeCKO synthesized guide library',M.OLIGO,state=poolstate,primary=True)
 clone_condition=g.input('pool_construction','source pooled-guide construction; original SupplementaryMethods unread',M.CONDITION)
 cloned=g.product('cloned','INFERRED representative guide-bearing native plasmid population',M.MOLECULAR,state=dna_state('plasmid'),inferred=True)
 g.react('construct_pool','INFERRED — original pooled-vector construction boundary',K.MAGIC,((vector,I.TEMPLATE),(pool,I.STRUCTURAL),(clone_condition,I.CONDITION)),(cloned,),reason='Published main identifies the synthesized-guide pooled library. Original detailed assembly and clone junctions are in unread SupplementaryMethods; an individual-guide duplex-ligation recipe is not substituted.')
 packaging=g.input('packaging','source lentiviral delivery preparation; detailed chemistry unreviewed',M.CONDITION)
 particle_state=Scene();particle_state.strand('RNA',[role('particle-associated vectorRNA; exact genome ends and copy/assembly chemistry unassigned')]);particle_state.strand('other',[role('independent particle proteins/lipids and possible carried polymers; not a DNA/RNA concatemer')])
 guide_inventory(particle_state)
 particles=g.product('particles','native source delivery particles withRNA genome',M.MOLECULAR,state=particle_state)
 plasmid_rest=g.product('plasmid_remainder','original plasmidDNA/unincorporated classes',M.MOLECULAR,state=cloned.state)
 g.react('package','Source delivery-particle preparation boundary',K.MAGIC,((cloned,I.TEMPLATE),(packaging,I.CONDITION)),(particles,plasmid_rest),reason='Original intracellular transcription/packaging chemistry is unreviewed. ParticleRNA is not relabelled unchanged plasmidDNA; detailed strand transformations remain unresolved.');bounded(g,plasmid_rest,'Delivery-particle production doesnot chemically erase original plasmidDNA; subsequent physical recovery unassigned.')
 cellstate=Scene();cellstate.strand('genomes',[role('native A375 cellular genomes before guide delivery')]);cellstate.strand('other',[role('independent cellular RNA/proteins/other polymers')])
 cells=g.input('cells','source A375 screen-cell population',M.BIOLOGICAL,state=cellstate,primary=True)
 delivery=g.input('delivery','source low-copy delivery; exact intracellular RT/integration chemistry unreviewed',M.CONDITION)
 integrated=g.product('integrated','INFERRED representative guide-bearing cellular genomes',M.MOLECULAR,state=dna_state('genome'),inferred=True)
 unintegrated=g.product('unintegrated','nonintegrated particle/cellular molecular classes',M.MOLECULAR,state=particle_state)
 g.react('deliver','INFERRED — source transduction/integrated-guide lineage',K.MAGIC,((cells,I.TEMPLATE),(particles,I.STRUCTURAL),(delivery,I.CONDITION)),(integrated,unintegrated),reason='Main article establishes delivery and genome-integrated guide lineage; molecular RT/integration mechanism, sites and multiplicity distribution remain outside inspected Methods. No guarantee of one inserted cassette per cell.');bounded(g,unintegrated,'Unintegrated particleRNA/cDNA and untransduced cells are not equated with physical disposal; subsequent recovery unknown.')
 screen=g.input('screen_condition','published A375 screen selection/timepoint; biological cohort andguideabundance change',M.SELECTION)
 selected=g.product('selected','INFERRED guide-bearing screen cohort; abundance distribution unassigned',M.MOLECULAR,state=dna_state('screen'),inferred=True)
 nonharvested=g.product('nonharvested','other sourcecohort molecular classes',M.MOLECULAR,state=integrated.state)
 g.react('screen','INFERRED — propagate source screen cohort to harvest',K.SELECTION,((integrated,I.TEMPLATE),(screen,I.SELECTION)),(selected,nonharvested),changes_molecule=False);bounded(g,nonharvested,'Nonharvested biologicalcohort/classes; physicaldisposalnotassigned.')
 extract=g.input('extraction','source screen genomicDNA extraction; detailed method unread',M.CONDITION)
 gdna=g.product('genomic_dna','native screen genomicDNA andpossible co-extractedclasses',M.MOLECULAR,state=dna_state('genome'),inferred=True)
 g.react('extract','INFERRED — extract native screen genomicDNA',K.EXTRACTION,((selected,I.TEMPLATE),(extract,I.CONDITION)),(gdna,))
 readout=g.input('readout_reagents','original screen library preparation; source primer/reagent/run declarations unverified',M.CONDITION)
 lib=g.product('library','original screen sequencing library; native endsunassigned',M.LIBRARY,state=dna_state('source library'))
 g.react('original_readout','Bounded original screen readout',K.MAGIC,((gdna,I.TEMPLATE),(readout,I.CONDITION)),(lib,),reason='Original screen SupplementaryMethods remain unread; neither a later genomicDNA one-PCR recipe nor a plasmid-only workbook is its verified native preparation.')
 g.sequence(lib,SequencingPlan(limitation='Original screen primers, barcode scheme and acquired runlengths unverified in this cappedpass. July2014 workbook is plasmid-specific; compatibleJoung2017 is a separately cited protocol, not this terminal.'))
 g.declare_reviewed_scope(('vector','guide_pool','pool_construction','construct_pool','packaging','package','cells','delivery','deliver','screen_condition','screen','extraction','extract','readout_reagents','original_readout','library'),SOURCE)
 return g.require_publishable()

EXTRA_MOLECULAR_STATES=(('Deposited empty-vector molecular components',Construct(list(dna_state('empty').strands['vector'].segs),name='deposited emptyvector context'),'Source-deposited empty-vector view; hU6/scaffold components and cPPT ordering are source-qualified. Unknown context/filler spans have no assigned length; this is not a cloned-guide member.'),)
