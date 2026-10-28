"""Study RRredBS and matchedBS arms; source5fC→5hmC→CMS remains physical."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Construct,Scene,Segment
from protocol_graph import ProtocolGraph,InputRole as R,MaterialKind as K,StepKind as S,SequencingPlan
from pcr import PCRGeometry
from protocol_model import role,molecule
import seqprimers as sp
PAGE_FILE='redbs-seq.html';TITLE='redBS-seq — genomic5fC reduction with matchedBS comparator';NOTES='01_redbs-seq.html'
SOURCE='<a href="https://doi.org/10.1038/nchem.1893">Booth et al.(2014)</a>,associated manuscriptPMC4188980;publisherSupplementaryMethodspp2–5.'
SUMMARY='MspI/TruSeq libraries undergo NaBH4 reduction or remain untreated, then bisulfite and KAPA Uracil+ PCR. 5fC becomes 5hmC then CMS in redBS, and U in ordinary BS; the difference in C readout estimates 5fC.'
CAVEAT='The represented study branch is MspI genomic RRredBS plus matched untreatedBS;the study oxBS arm and synthetic-control library recipes are not merged. ChemicalCMS is not unchangedC. Base-class examples do not assign genomic coordinates or perfect conversion. Exact originalTruSeq adapters/primer members and VORMethods remain unreviewed. Table3 synthetic-control primers are not substituted for genomicTruSeq primers.'
INITIAL=('C','5mC','5hmC','5fC')
@dataclass(frozen=True)
class RedState:
 stage:str
 carryover:str=''
 def __post_init__(self):
  if self.stage not in {'input','reduced','red_bisulfite','bs_bisulfite'}:raise ValueError('unknownredBSchemicalstage')
 def rows(self):
  values={'input':INITIAL,'reduced':('C','5mC','5hmC','5hmC'),'red_bisulfite':('U','5mC','CMS','CMS'),'bs_bisulfite':('U','5mC','CMS','U')}[self.stage]
  sc=Scene()
  for n,(old,new) in enumerate(zip(INITIAL,values)):
   sc.strand(str(n),[role('genomicflank'),Segment('physicalbase '+new,'U' if new=='U' else 'C'),role('genomicflank')],label=old+' startingbaseclass')
   sc.note(str(n),'Physicalproduct:'+new+(';reported'+('T' if new=='U' else 'C') if 'bisulfite' in self.stage else '')+';source idealizedpath,notquantitativeyield')
  if self.carryover:sc.note('0','Separatepreparationcomponents:'+self.carryover+';notcovalentlypartofDNA.')
  return sc.rows()
FINAL_LIBRARY=Construct((role('unassembledgenomicTruSeqend'),role('genome-derivedPCRsequence;redBS/BSreadoutmappingabove'),role('unassembledoppositeTruSeqend')),name='redBS/BSgenomiclibraryroles')
FINAL_CAPTION='Source-positive genomiclibraryroles;historicaladapter/indexmembersandfullstrandassemblyunassigned.'
SEQ_PRIMERS=();SEQUENCING_UNAVAILABLE='Source: HiSeq 2000, single-end 100 bases. Exact run primers, index mask and study adapter members remain unassigned.'
SEQUENCING_ENDING='Compare C fractions between redBS and BS at matched genomic positions. Incomplete conversion and background require source controls; a C call does not uniquely identify 5fC.'
READ_SPECS=(sp.read_spec('Read 1',100,'bisulfite-derivedgenomicsequence;CMS/5mCreadC,UreadT','include'),)
def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Booth2014publisherSupplementaryMethods,genomicRRredBS/BSarms;associatedarticlechemicalrationale')
 dna=g.input('dna','purifiedmouseESgenomicDNA;representative1µgpreparations',K.MOLECULAR,state=RedState('input'),primary=True)
 mspi=g.input('mspi','MspIovernight37C',K.CATALYST)
 cut=g.product('cut','MspIdigestedDNA',K.MOLECULAR,state=RedState('input'))
 g.react('digest','Digest genomic DNA with MspI',S.DIGESTION,((dna,R.TEMPLATE),(mspi,R.CATALYST)),(cut,))
 edta=g.input('edta','EDTAquench',K.CONDITION)
 stopped=g.product('stopped','quench-reactionDNApreparation',K.MOLECULAR,state=RedState('input','EDTA/digestioncomponents'))
 g.react('stop_digest','Quench restriction digest with EDTA',S.QUENCHING,((cut,R.TEMPLATE),(edta,R.CONDITION)),(stopped,))
 def clean(i,x,label,*,outstate=None,kind=K.MOLECULAR):
  r=g.input(i+'_condition',label,K.CONDITION);y=g.product(i+'_dna','recoveredDNA',kind,state=x.state if outstate is None else outstate);w=g.product(i+'_waste','physicallyremovedfractions',K.WASTE)
  g.react(i,label,S.CLEANUP,((x,R.TEMPLATE),(r,R.CONDITION)),(y,w));g.discard(w,'Actualsourceworkupseparatesremovedfractions;nocompletepurityclaim.');return y
 current=clean('digest_cleanup',stopped,'GeneJETPCRpurification;EBelution',outstate=RedState('input'))
 controls=g.input('controls','source5fC/5hmCsyntheticDNAcontrolmixture',K.MOLECULAR,state=molecule('spike-inDNA',(role('syntheticcontrolDNAclasses;notgenomicprimers'),)))
 tru=g.input('truseq','sourceTruSeqkit;1:10adapters;historicalidentityunassembled',K.CONDITION)
 prepared=g.product('prepared','independentlypreparedmatchedTruSeqlibraries',K.MOLECULAR,state=RedState('input'))
 g.react('library_prep','Source TruSeq library preparation after control spike-in',S.MAGIC,((current,R.TEMPLATE),(controls,R.SPIKE_IN),(tru,R.CONDITION)),(prepared,),reason='SourcespecifiesgenomicTruSeqpreparation,adapterdilution/finalworkup;exacthistoricalkitendprocessingandadapterpolymersnotreconstructed. Nosynthetic-controlprimerrecipeimported.')
 current=clean('truseq_cleanup',prepared,'finalAMPureworkup;dry30min;waterelution')
 allocation=g.input('allocation','allocateindependentlypreparedlibraries:study4redBS,4BS;oxBSarmoutsidegraph',K.CONDITION)
 red=g.product('red_input','preparedlibrariesallocatedredBS',K.MOLECULAR,state=RedState('input'))
 bs=g.product('bs_input','matchedpreparedlibrariesallocateduntreatedBS',K.MOLECULAR,state=RedState('input'))
 g.react('allocate','Allocate matched prepared libraries to analytical arms',S.SPLIT,((current,R.TEMPLATE),(allocation,R.CONDITION)),(red,bs))
 bh4=g.input('borohydride','5µLfresh1MNaBH4+15µLDNA;250mMfinal,1hroomtempdark',K.CONDITION)
 reduced=g.product('reduced','5fCreducedto5hmC',K.MOLECULAR,state=RedState('reduced','reductionreagentsandproducts'))
 g.react('reduce','Reduce 5fC to 5hmC with NaBH4',S.CHEMICAL_CONVERSION,((red,R.TEMPLATE),(bh4,R.CONDITION)),(reduced,))
 acetate=g.input('acetate','10µL750mMsodiumacetate,pH5;wait≥10minoruntilgasreleaseends',K.CONDITION)
 quenched=g.product('quenched','reducedDNAafterquench',K.MOLECULAR,state=RedState('reduced','aqueousquenchcomponents;sourcegasventing,separationnotyetperformed'))
 g.react('quench_reduction','Quench reduction; vent gas without implying purification',S.QUENCHING,((reduced,R.TEMPLATE),(acetate,R.CONDITION)),(quenched,))
 forward=g.input('forward_primer','genomicTruSeqprimermixforwardpolymer;memberunassigned',K.OLIGO,state=molecule('forwardprimer',(role('exactgenomicTruSeqforwardprimerunassigned'),)))
 reverse=g.input('reverse_primer','genomicTruSeqprimermixreversepolymer;memberunassigned',K.OLIGO,state=molecule('reverseprimer',(role('exactgenomicTruSeqreverseprimerunassigned'),)))
 kapa=g.input('kapa','KAPAHiFiUracil+mastermix',K.CATALYST)
 for prefix,x,chemical in [('red',quenched,'red_bisulfite'),('bs',bs,'bs_bisulfite')]:
  reagent=g.input(prefix+'_bisulfite_condition','EpiTect:30µLDNA,80µLmix,30µLprotect;twoFFPEcycles',K.CONDITION)
  converted=g.product(prefix+'_converted','sourcebisulfiteproducts',K.MOLECULAR,state=RedState(chemical,'bisulfitereactioncomponentsuntilkitworkup'))
  g.react(prefix+'_bisulfite','Bisulfite conversion; physical products differ from readout',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(reagent,R.CONDITION)),(converted,))
  purified=clean(prefix+'_workup',converted,'EpiTectFFPEkitworkup',outstate=RedState(chemical))
  amp=g.product(prefix+'_amplified',prefix+'PCRlibrary',K.MOLECULAR,state=FINAL_LIBRARY)
  g.react(prefix+'_pcr','KAPA Uracil+ PCR, 12 cycles at source 98/60/72°C program',S.PCR,((purified,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(kapa,R.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(molecule('bisulfiteDNA',(role('sourcekitadapter-bearingconvertedDNA'),)),reason='ExacthistoricalgenomicTruSeqprimermembersandstrandgeometryunassigned;controlTable3oligosnotborrowed.'))
  library=clean(prefix+'_final',amp,'AMPure50µL;Qubit/TapeStationQC;pool',kind=K.LIBRARY)
  g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE,read_lengths={'Read 1':100},acquired_roles=('Read 1',)))
 g.declare_reviewed_scope(('dna','mspi','cut','digest','edta','stopped','stop_digest','digest_cleanup_condition','digest_cleanup_dna','digest_cleanup_waste','digest_cleanup','controls','truseq','prepared','library_prep','truseq_cleanup_condition','truseq_cleanup_dna','truseq_cleanup_waste','truseq_cleanup','allocation','red_input','bs_input','allocate','borohydride','reduced','reduce','acetate','quenched','quench_reduction','forward_primer','reverse_primer','kapa','red_bisulfite_condition','red_converted','red_bisulfite','red_workup_condition','red_workup_dna','red_workup_waste','red_workup','red_amplified','red_pcr','red_final_condition','red_final_dna','red_final_waste','red_final','bs_bisulfite_condition','bs_converted','bs_bisulfite','bs_workup_condition','bs_workup_dna','bs_workup_waste','bs_workup','bs_amplified','bs_pcr','bs_final_condition','bs_final_dna','bs_final_waste','bs_final'),g.authority)
 return g.require_publishable()
def render_page():
 from batch_page import render
 return render(sys.modules[__name__])
def selftest():
 protocol_graph();render_page();print('PASS source-qualifiedredBS/BSgraphandrendering;noexactgenomicoligoclaims')
