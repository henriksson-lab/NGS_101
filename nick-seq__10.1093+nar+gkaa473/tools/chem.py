"""Published Nick-seq positive NT/TdT branches; source-local end/backbone states."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Construct,Scene,Segment,feature
from protocol_graph import InputRole as R,MaterialKind as K,StepKind as S,ProtocolGraph,SequencingPlan
from protocol_model import molecule,role
import seqprimers as sp
PAGE_FILE='nick-seq.html'
TITLE='Nick-seq — complementary nick-translation and polyT libraries'
NOTES='01_nick-seq.html'
SOURCE='<a href="https://doi.org/10.1093/nar/gkaa473">Cao et al. (2020)</a>, published Materials and Methods, mechanism and supplementary source inventory.'
SUMMARY='Block old extendable ends, generate target nicks, then split DNA into phosphorothioate-protected nick-translation and TdT polyT-tailing branches. Separate libraries provide complementary mapping evidence.'
CAVEAT='This graph represents E.coli nickase-validation material and the two positive libraries. Source negative-control library processing and alternative iodine/EndoIV front ends remain in the notes, not silently reconstructed. Known end/backbone classes are distinguished from unassigned genomic coordinates and reaction yields. SMART ChIP-seq kit internals/oligos remain explicit unvalidated boundaries. The published AP-site application differs in H2O2 concentration: Methods lists0.1/0.5/1/2mM; Results describes0.2mM. This nickase branch depends on neither value. Supplemental DOCX graphics remain unread after a failed local render; no sequence absence is inferred.'
@dataclass(frozen=True)
class NickState:
 end_class:str
 backbone:str='native DNA; positions unassigned'
 carryover:str=''
 def rows(self):
  sc=Scene();sc.strand('dna',[role('genomic sequence; coordinate unassigned')],label=self.end_class)
  sc.note('dna','Backbone: '+self.backbone)
  if self.end_class=='TdT polyT-tailed DNA':sc.strand('tail',[role('polyT attached to eligible DNA3-prime end; length unassigned')],label='source terminal extension; no fixed tail length')
  if self.carryover:sc.strand('carryover',[role(self.carryover)],label='material still in reaction until physical cleanup')
  return sc.rows()

def library_roles(prefix):return Construct((role('unassembled SMART ChIP library end'),role(prefix+' genome-derived sequence'),Segment('source sample barcode; position/member/length unassigned','[sample index]',placeholder=True,feature=feature('nick_'+prefix.lower()+'_sample_index','sample_index','unknown')),role('unassembled opposite SMART ChIP end')),name=prefix+' source library roles')
NT_LIBRARY=library_roles('NT')
TDT_LIBRARY=library_roles('TdT')
FINAL_LIBRARY=Construct((role('branch-specific unassembled kit end'),role('NT or TdT genome-derived library sequence'),role('branch-specific unassembled opposite end')),name='Nick-seq library-role overview')
FINAL_CAPTION='Role-only overview; distinct NT and TdT terminal libraries do not imply a common chemically assembled product or exact adapter architecture.'
SEQ_PRIMERS=()
SEQUENCING_UNAVAILABLE='NextSeq500 paired-end75bp; exact historical kit adapter/primer polymers, acquired barcode mask and selected barcode members remain unassigned.'
SEQUENCING_ENDING='NT and TdT readouts are combined computationally. First3bp of Read1 are trimmed in source processing; they are not declared a UMI. GenomicT adjacency can confound tail-based calling, motivating complementary NT information.'
READ_SPECS=(sp.read_spec('Read 1',75,'NT- or TdT-derived genomic sequence; source trims first3bp','include'),sp.read_spec('Read 2',75,'opposite library end; no exact technical-cycle assignment supplied','include'))

def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Cao2020 published Materials and Methods; positive nickase-validation NT/TdT branches')
 dna=g.input('dna','E.coliDH10B genomic DNA,1µg',K.MOLECULAR,state=NickState('starting DNA end mixture',carryover='possible RNA carryover'),primary=True)
 enzymes=g.input('fragmentation_enzymes','three SEPARATE digests:NciI;HindIII+XhoI;SalI+XbaI+NdeI;RNaseA',K.CATALYST)
 fragments=g.product('fragments','restriction-fragmented DNA preparations',K.MOLECULAR,state=NickState('pre-existing3-primeOH/3-primeP end mixture',carryover='RNA digestion products'))
 g.react('fragment','Digest in three separate preparations with RNaseA',S.DIGESTION,((dna,R.TEMPLATE),(enzymes,R.CATALYST)),(fragments,))
 def recover(prefix,template,label,*,state=None,kind=K.MOLECULAR):
  cond=g.input(prefix+'_condition',label,K.CONDITION)
  retained=g.product(prefix+'_dna','retained DNA',kind,state=template.state if state is None else state)
  waste=g.product(prefix+'_waste','removed recovery fractions',K.WASTE)
  g.react(prefix,label,S.CLEANUP,((template,R.TEMPLATE),(cond,R.CONDITION)),(retained,waste));g.discard(waste,'Physical source cleanup removes separated fractions; exact composition and recovery efficiency unassigned.')
  return retained
 current=recover('fragment_cleanup',fragments,'Qiagen PCR purify each digest;pool',state=NickState('pre-existing DNA end mixture'))
 sap=g.input('sap','shrimp alkaline phosphatase',K.CATALYST)
 oldoh=g.product('old_oh','DNA preparation with3-primeP class dephosphorylated',K.MOLECULAR,state=NickState('pre-existing3-primeOH substrate class; no complete conversion guarantee',carryover='phosphatase reaction products'))
 g.react('dephosphorylate_old','SAP37C30min;heat70C10min',S.DEPHOSPHORYLATION,((current,R.TEMPLATE),(sap,R.CATALYST)),(oldoh,))
 pol=g.input('pol_i','DNA polymeraseI',K.CATALYST)
 dd=g.input('ddntp','ddNTP blocking mixture',K.MOLECULAR)
 blocked=g.product('blocked','dd-terminated pre-existing extendable end class',K.MOLECULAR,state=NickState('old extendable ends dd-blocked; no free3-primeOH at blocked class',carryover='excess ddNTP and phosphatase products'))
 g.react('block_old','PolI/ddNTP block37C40min',S.EXTENSION,((oldoh,R.TEMPLATE),(dd,R.STRUCTURAL),(pol,R.CATALYST)),(blocked,))
 processed=g.product('processed_block_mix','blocked DNA plus processed excess nucleotides',K.MOLECULAR,state=NickState('old extendable ends dd-blocked',carryover='processed excess ddNTP/small-molecule products,still in tube'))
 g.react('process_excess_ddntp','SAP37C30min;heat75C10min; excess-nucleotide processing',S.CHEMICAL_CONVERSION,((blocked,R.TEMPLATE),(sap,R.CATALYST)),(processed,))
 current=recover('blocking_cleanup',processed,'DyeEx desalting',state=NickState('old extendable ends dd-blocked'))
 nickases=g.input('nickases','Nb.BsmI andNb.BsrDI;65C1h,80C20min,slow cooling',K.CATALYST)
 nicked=g.product('nicked','new target nick3-primeOH sites with old ends blocked',K.MOLECULAR,state=NickState('new nick3-primeOH versus old dd-blocked ends'))
 g.react('create_nicks','Create target nicks after old-end blocking; no additional cleanup',S.DIGESTION,((current,R.TEMPLATE),(nickases,R.CATALYST)),(nicked,))
 split=g.input('split_condition','divide sample between NT andTdT reactions',K.CONDITION)
 nt=g.product('nt_input','NT aliquot',K.MOLECULAR,state=nicked.state)
 tdt=g.product('tdt_input','TdT aliquot',K.MOLECULAR,state=nicked.state)
 g.react('split','Split between complementary library chemistries',S.SPLIT,((nicked,R.TEMPLATE),(split,R.CONDITION)),(nt,tdt))
 alpha=g.input('alpha_thio','alpha-thio-dNTP mixture;15C90min',K.MOLECULAR)
 ntproduct=g.product('nt_product','PT-containing copied DNA within original template preparation',K.MOLECULAR,state=NickState('nick-translated DNA preparation','new copied phosphorothioate backbone plus original DNA; copied span unassigned',carryover='unincorporated nucleotide/reaction products'))
 g.react('nick_translate','PolI alpha-thio-dNTP nick translation;heat75C20min',S.EXTENSION,((nt,R.TEMPLATE),(alpha,R.STRUCTURAL),(pol,R.CATALYST)),(ntproduct,))
 ntclean=recover('nt_cleanup',ntproduct,'DyeEx cleanup before nuclease treatment',state=NickState('nick-translated DNA','PT-containing copied tract plus original DNA'))
 exo=g.input('exoiii','ExonucleaseIII200U;37C60min',K.CATALYST)
 exoproduct=g.product('exo_product','ExoIII-treated DNA preparation',K.MOLECULAR,state=NickState('ExoIII-treated preparation','PT-containing copied tract retained as source target class',carryover='DNA digestion fragments and remaining background DNA'))
 g.react('digest_exoiii','ExoIII digest; debris remains in tube',S.DIGESTION,((ntclean,R.TEMPLATE),(exo,R.CATALYST)),(exoproduct,))
 heat=g.input('nt_heat','95C3min;ice quench',K.CONDITION)
 denatured=g.product('nt_denatured','heat-denatured ExoIII preparation',K.MOLECULAR,state=exoproduct.state)
 g.react('denature_nt','Denature before RecJf digestion',S.THERMAL_TREATMENT,((exoproduct,R.TEMPLATE),(heat,R.CONDITION)),(denatured,))
 recj=g.input('recjf','RecJf60U;37C60min; optional additional60U for complex samples;80C10min stop',K.CATALYST)
 digested=g.product('nt_digested','PT-resistant oligonucleotide class plus digestion products',K.MOLECULAR,state=NickState('nuclease-treated NT preparation','PT-containing copied oligonucleotide target class',carryover='ExoIII/RecJf digestion products and possible remaining background'))
 g.react('digest_recjf','RecJf digest; retain products until physical purification',S.DIGESTION,((denatured,R.TEMPLATE),(recj,R.CATALYST)),(digested,))
 purified=recover('nt_oligo_cleanup',digested,'Zymo Oligo Clean & Concentrator',state=NickState('recovered NT target oligonucleotide class','source PT-containing copied DNA; purity/yield unassigned'))
 tdt_heat=g.input('tdt_heat','95C3min in20µL water',K.CONDITION)
 tdt_denatured=g.product('tdt_denatured','denatured nick-converted DNA',K.MOLECULAR,state=tdt.state)
 g.react('denature_tdt','Denature nick-converted DNA before tailing',S.THERMAL_TREATMENT,((tdt,R.TEMPLATE),(tdt_heat,R.CONDITION)),(tdt_denatured,))
 enzyme=g.input('tdt','Clontech TdT andDNA SMART T-Tailing Mix;37C20min,70C10min stop',K.CATALYST)
 tailed=g.product('tdt_tailed','TdT polyT-tailed DNA',K.MOLECULAR,state=NickState('TdT polyT-tailed DNA'))
 g.react('tail_t','Add source polyT tail to eligible3-prime ends',S.TAILING,((tdt_denatured,R.TEMPLATE),(enzyme,R.CATALYST)),(tailed,))
 for prefix,substrate,library_state in [('nt',purified,NT_LIBRARY),('tdt',tailed,TDT_LIBRARY)]:
  kit=g.input(prefix+'_kit','Clontech SMART ChIP-seq kit; historical revision unassigned',K.CONDITION)
  kit_oligos=g.input(prefix+'_kit_oligos','source kit adapter/primer polymers; exact members unassembled',K.OLIGO,state=molecule('kit oligos',(role('historical SMART ChIP-seq oligo polymers unassigned'),)))
  lib=g.product(prefix+'_kit_product',prefix+' indexed kit PCR library',K.MOLECULAR,state=library_state)
  g.react(prefix+'_kit_library',('NT recovered-oligo library preparation' if prefix=='nt' else 'TdT primer annealing/template switching')+';source kit12-cycle finalPCR',S.MAGIC,((substrate,R.TEMPLATE),(kit,R.CONDITION),(kit_oligos,R.STRUCTURAL)),(lib,),reason='Source positively names SMART ChIP-seq kit workflow and12PCR cycles; historical kit chemistry/oligo geometry not reconstructed here. No generic RNA/Nextera recipe imported.')
  selected=recover(prefix+'_final_cleanup',lib,'AMPure selection of represented positive class; source matched-control pooling branch deferred',kind=K.LIBRARY)
  g.sequence(selected,SequencingPlan(library=library_state,limitation=SEQUENCING_UNAVAILABLE,read_lengths={'Read 1':75,'Read 2':75},acquired_roles=('Read 1','Read 2')))
 g.declare_reviewed_scope(('dna','fragmentation_enzymes','fragments','fragment','fragment_cleanup_condition','fragment_cleanup_dna','fragment_cleanup_waste','fragment_cleanup','sap','old_oh','dephosphorylate_old','pol_i','ddntp','blocked','block_old','processed_block_mix','process_excess_ddntp','blocking_cleanup_condition','blocking_cleanup_dna','blocking_cleanup_waste','blocking_cleanup','nickases','nicked','create_nicks','split_condition','nt_input','tdt_input','split','alpha_thio','nt_product','nick_translate','nt_cleanup_condition','nt_cleanup_dna','nt_cleanup_waste','nt_cleanup','exoiii','exo_product','digest_exoiii','nt_heat','nt_denatured','denature_nt','recjf','nt_digested','digest_recjf','nt_oligo_cleanup_condition','nt_oligo_cleanup_dna','nt_oligo_cleanup_waste','nt_oligo_cleanup','tdt_heat','tdt_denatured','denature_tdt','tdt','tdt_tailed','tail_t','nt_kit','nt_kit_oligos','nt_kit_product','nt_kit_library','nt_final_cleanup_condition','nt_final_cleanup_dna','nt_final_cleanup_waste','nt_final_cleanup','tdt_kit','tdt_kit_oligos','tdt_kit_product','tdt_kit_library','tdt_final_cleanup_condition','tdt_final_cleanup_dna','tdt_final_cleanup_waste','tdt_final_cleanup'),g.authority)
 return g.require_publishable()

def render_page():
 from batch_page import render
 return render(sys.modules[__name__])

def selftest():
 # No source-exact ordering oligos asserted in this partial model. Structural
 # graph/source coverage and renderer construction are the available local checks.
 protocol_graph();render_page()
 print('PASS Nick-seq partial graph and rendering; no source-exact oligo checks asserted')
