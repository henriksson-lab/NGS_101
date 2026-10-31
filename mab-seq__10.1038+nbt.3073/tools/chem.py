"""Genome-scale MAB-seq: CpG-qualified physical products and reported bases."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, oligo, feature
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole as R, MaterialKind as K, ProtocolGraph, SequencingPlan, StepKind as S, Disposition, Terminal
from protocol_model import molecule, role
PAGE_FILE='mab-seq.html'
TITLE='MAB-seq — CpG protection before bisulfite conversion'
NOTES='01_mab-seq.html'
SOURCE='<a href="https://doi.org/10.1038/nbt.3073">Wu et al. (2014)</a>; associated accepted manuscript PMC4269366, genome-scale Methods and publisher Supplementary Information.'
SUMMARY='Methylated adapters are ligated before two-round M.SssI/SAM protection. Bisulfite converts 5fC and 5caC into products read as T; KAPA Uracil+ PCR generates the sequencing library.'
CAVEAT='CpG-qualified alternative molecular classes are shown, not assigned genomic positions or perfect conversion. Chemical5caU is kept distinct from reportedT; the5hmC bisulfite derivative remains unassigned. This graph represents the genome-scale route, whose adapter ligation precedes methylase treatment. Locus-specific four-round treatment, caMAB and redBS branches are not substituted. The complete published VOR Methods and exact historical NEB index-primer members remain unreviewed; source-associated accepted-manuscript identity is explicit.'
ADAPTER_F=il.TRUSEQ_READ1
ADAPTER_R=il.NEBNEXT_ARM_READ2
INITIAL=('C','5mC','5hmC','5fC','5caC')
PROTECTED=('5mC','5mC','5hmC','5fC','5caC')
CONVERTED=('5mC','5mC','unassigned bisulfite derivative of5hmC','U','5caU')
REPORTED=('C','C','C','T','T')
@dataclass(frozen=True)
class CpGState:
 stage:str
 adapters:bool=False
 def __post_init__(self):
  if self.stage not in {'input','protected','bisulfite'}:raise ValueError('unrecognized MAB chemical stage')
 def rows(self):
  sc=Scene();states={'input':INITIAL,'protected':PROTECTED,'bisulfite':CONVERTED}[self.stage]
  for i,(original,current) in enumerate(zip(INITIAL,states)):
   base=role(current) if current.startswith('unassigned') else Segment('physical base: '+current,'U' if current in {'U','5caU'} else 'C')
   sc.strand(str(i),[role('genomic flank'),base,Segment('CpG neighboring G','G'),role('genomic flank')],label=original+' starting CpG case')
   sc.note(str(i),'Physical state: '+current+('; sequencing reports '+REPORTED[i] if self.stage=='bisulfite' else '')+('; source methylated adapters retained; exact strand assembly unassigned' if self.adapters else ''))
  return sc.rows()
FINAL_LIBRARY=Construct((role('unassembled NEBNext primer-derived end'),role('genome-derived PCR product; CpG physical-to-readout mapping shown above'),Segment('sample index; member and length unassigned','[i7 index]',placeholder=True,feature=feature('mab_i7','sample_index','unknown')),role('opposite unassembled primer-derived end')),name='MAB library roles; exact indexed PCR polymers unresolved')
FINAL_CAPTION='Source-established library roles; exact historical indexed/universal primer members and full endpoint assembly remain unassigned.'
SEQ_PRIMERS=()
SEQUENCING_ENDING='Genome-scale libraries: HiSeq2500 single-end60or76bp. Final PCR cycles are qPCR-determined, not the five-cycle pre-amplification count.'
SEQUENCING_UNAVAILABLE='HiSeq2500 single-end60or76bp is source reported; exact run-primer molecules, acquired index mask and selected index member remain unassigned.'
READ_SPECS=(sp.read_spec('Read 1','60or76','bisulfite-derived genomic sequence; CpG T signal combines5fC/5caC and incomplete protection background','include'),sp.read_spec('Index 1 (i7)','unassigned','NEBNext Index Primers Set1 sample index; member unassigned','undetermined'))

def adapter_scene():
 sc=Scene();core=ADAPTER_F[-13:-1]
 sc.strand('forward',[Segment('forward Y arm',ADAPTER_F[:-13]),Segment('paired core',core),Segment('terminal T','T')],label='methylated custom forward',mod3='final bond phosphorothioate')
 sc.anneal('reverse',[Segment('reverse paired core',ADAPTER_R[:12]),Segment('reverse Y arm',ADAPTER_R[12:])],to='forward',pair=('reverse paired core','paired core'),label='methylated custom reverse',mod5='5Phos',unpaired=('reverse Y arm',))
 sc.note('forward','Source says5mC-methylated adapters; per-position modification map is not independently assigned. Letters preserve source base identities.')
 sc.ligation_end('forward','terminal T',kind='cohesive',polarity='3-prime',owner='adapter')
 return sc

def oligos():return (oligo('methylated adapter forward',[Segment('forward',ADAPTER_F)],mods='source5mC-methylated; per-position map unassigned; final phosphorothioate bond'),oligo('methylated adapter reverse',[Segment('reverse',ADAPTER_R)],mods='5Phos; source5mC-methylated; per-position map unassigned'))

def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Wu2014 associated manuscript, genome-scale Library generation Methods; publisher SI controls',representation='published_instance')
 dna=g.input('dna','unenriched genomic DNA,1µg',K.MOLECULAR,state=CpGState('input'),primary=True)
 spike=g.input('lambda','unmethylated lambda control DNA,1:400',K.MOLECULAR,state=molecule('control DNA',(role('unmethylated lambda DNA'),)))
 repaired=g.product('repaired','spiked-in, end-repaired DNA',K.MOLECULAR,state=CpGState('input'))
 kit=g.input('endprep','NEBNext Ultra DNA LibraryPrep end repair; exact genomic end geometry unassembled',K.CATALYST)
 g.react('repair','Spike unmethylated lambda and repair DNA ends',S.END_REPAIR,((dna,R.TEMPLATE),(spike,R.SPIKE_IN),(kit,R.CATALYST)),(repaired,))
 adapters=g.input('adapters','5mC custom Y-adapter pair',K.ADAPTER,state=adapter_scene())
 ligase=g.input('ligation_mix','source NEBNext Ultra adapter-ligation mix',K.CATALYST)
 ligated=g.product('ligated','methylated-adapter-ligated DNA',K.MOLECULAR,state=CpGState('input',True))
 g.react('ligate','Ligate methylated custom adapters before M.SssI',S.LIGATION,((repaired,R.TEMPLATE),(adapters,R.STRUCTURAL),(ligase,R.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.known('cohesive',adapter_scene(),statement='Printed forward/reverse adapters define12-bp duplex core and3-primeT end. Exact genomic repaired-end assembly not declared.'))
 def cleanup(prefix,substrate,label,*,state=None,kind=K.MOLECULAR):
  cond=g.input(prefix+'_condition',label,K.CONDITION)
  retained=g.product(prefix+'_dna','retained DNA',kind,state=state if state is not None else substrate.state)
  waste=g.product(prefix+'_waste','discarded purification fractions',K.WASTE)
  g.react(prefix,label,S.CLEANUP,((substrate,R.TEMPLATE),(cond,R.CONDITION)),(retained,waste));g.discard(waste,'Source removes purification supernatants/washes; recovery efficiency unassigned.')
  return retained
 current=cleanup('adapter_cleanup',ligated,'1.2xAMPure cleanup')
 enzyme=g.input('msssi','M.SssI CpG methyltransferase M0226M',K.CATALYST)
 for number,buffer in [(1,'Mg-free Tris/NaCl/EDTA'),(2,'NEBbuffer2 with10mMMgCl2')]:
  cond=g.input('round'+str(number)+'_conditions',buffer+';1U/µL M.SssI,640µMSAM;4h in25µL then add25µL and8h',K.CONDITION)
  sam=g.input('round'+str(number)+'_sam','Actual SAM supply for this round: initial and supplemental reaction inventories,640µM source concentration; not a buffer-only condition',K.MOLECULAR)
  product=g.product('protected'+str(number),'CpG-protected DNA,round'+str(number),K.MOLECULAR,state=CpGState('protected',True))
  g.react('methylate'+str(number),'Round'+str(number)+' compound M.SssI/SAM protection',S.CHEMICAL_CONVERSION,((current,R.TEMPLATE),(enzyme,R.CATALYST),(sam,R.STRUCTURAL),(cond,R.CONDITION)),(product,),reason='One compound source round: initial25µL for4h then fresh25µL same-concentration M.SssI/SAM for8h. Aggregate reagent inventory, not a single simultaneous50µL incubation. Source concentrations do not establish complete protection or yield.')
  current=cleanup('round'+str(number)+'_purify',product,'PCI extraction and ethanol precipitation after this round')
 bs=g.input('bisulfite_condition','EpiTect59104; kit thermal cycle repeated twice',K.CONDITION)
 converted=g.product('converted','bisulfite-treated DNA; physical products retained',K.MOLECULAR,state=CpGState('bisulfite',True))
 g.react('bisulfite','Convert source-qualified CpG substrate classes',S.CHEMICAL_CONVERSION,((current,R.TEMPLATE),(bs,R.CONDITION)),(converted,))
 forward=g.input('universal_primer','NEBNext Multiplex Oligos Set1 universal primer; exact historical polymer unassigned',K.OLIGO,state=molecule('universal PCR oligo',(role('exact universal primer unassigned'),)))
 reverse=g.input('index_primer','NEBNext Index Primers Set1 indexed primer; member unassigned',K.OLIGO,state=molecule('indexed PCR oligo',(Segment('index-bearing primer, exact polymer unassigned','[indexed primer]',placeholder=True,feature=feature('mab_i7','sample_index','unknown')),)))
 kapa_stock=g.input('kapa','Actual KAPA HiFi Uracil+ HotStart ReadyMix; aggregate preamp/final-PCR allocations unassigned',K.MOLECULAR)
 kapa=g.product('kapa_enzyme','Polymerase component of the same actual ReadyMix',K.CATALYST)
 kapa_nt=g.product('kapa_nt','Nucleotide component of the same actual ReadyMix',K.MOLECULAR)
 kapa_buffer=g.product('kapa_buffer','Buffer component of the same actual ReadyMix',K.CONDITION)
 g.react('kapa_components','Expose inventory of the source ReadyMix',S.MAGIC,((kapa_stock,R.TEMPLATE),),(kapa,kapa_nt,kapa_buffer),changes_molecule=False,reason='Same actual named ReadyMix supply, not three fresh stocks. Exact formulation and allocations among source preamplification and final-PCR reactions remain unassigned.')
 pre=g.product('preamplified','five-cycle preamplified DNA',K.MOLECULAR,state=FINAL_LIBRARY)
 geom=PCRGeometry.unavailable(molecule('bisulfite-treated template',(role('source adapter-bearing bisulfite DNA; exact binding polymers unresolved'),)),reason='Exact historical universal/index primer members and complete strand assembly unassigned; no binding geometry invented.')
 g.react('preamp','Pre-amplify five cycles with KAPA Uracil+',S.PCR,((converted,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(kapa,R.CATALYST),(kapa_nt,R.STRUCTURAL),(kapa_buffer,R.CONDITION)),(pre,),geometry=geom)
 split=g.input('split_condition','reserve analytical aliquot for qPCR cycle determination',K.CONDITION)
 qc=g.product('qpcr_aliquot','aliquot for qPCR',K.QC,state=pre.state)
 retained=g.product('retained_preamp','retained material for scaled-up final PCR',K.MOLECULAR,state=pre.state)
 g.react('split_qpcr','Reserve qPCR aliquot; retain library material',S.SPLIT,((pre,R.TEMPLATE),(split,R.CONDITION)),(qc,retained))
 g.mark_terminal(Terminal(qc.id,Disposition.QC,'Source determines optimal cycles by qPCR; exact analytical reaction setup not reconstructed here. This aliquot is not pooled into the final library.'))
 cycles=g.input('final_cycle_choice','final cycles determined by qPCR; numeric assignment unassigned',K.CONDITION)
 amp=g.product('amplified','scaled-up PCR products',K.MOLECULAR,state=FINAL_LIBRARY)
 g.react('final_pcr','Scaled-up PCR at qPCR-determined cycle count',S.PCR,((retained,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(kapa,R.CATALYST),(kapa_nt,R.STRUCTURAL),(kapa_buffer,R.CONDITION),(cycles,R.CONDITION)),(amp,),geometry=PCRGeometry.unavailable(FINAL_LIBRARY,reason='Selected historical primer polymers remain unassigned.'))
 library=cleanup('final_cleanup',amp,'1.2xAMPure final cleanup',kind=K.LIBRARY)
 loading=g.product('loading_aliquot','Selected sequencing loading aliquot; volume unassigned',K.LIBRARY,state=library.state)
 remaining=g.product('library_remainder','Untransferred purified library stock; amount unassigned',K.LIBRARY,state=library.state)
 g.react('select_loading','Select sequencing loading aliquot',S.SPLIT,((library,R.TEMPLATE),),(loading,remaining),changes_molecule=False)
 g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Loading does not imply consumption of the whole physical stock. Amounts remain unassigned; no multi-sample pooling arrangement is inferred.'))
 g.sequence(loading,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE,acquired_roles=('Read 1',)))
 g.declare_reviewed_scope(('round1_sam','round2_sam','kapa_components','select_loading','loading_aliquot','library_remainder','dna','lambda','repaired','endprep','repair','adapters','ligation_mix','ligated','ligate','adapter_cleanup_condition','adapter_cleanup_dna','adapter_cleanup_waste','adapter_cleanup','msssi','round1_conditions','protected1','methylate1','round1_purify_condition','round1_purify_dna','round1_purify_waste','round1_purify','round2_conditions','protected2','methylate2','round2_purify_condition','round2_purify_dna','round2_purify_waste','round2_purify','bisulfite_condition','converted','bisulfite','universal_primer','index_primer','kapa','preamplified','preamp','split_condition','qpcr_aliquot','retained_preamp','split_qpcr','final_cycle_choice','amplified','final_pcr','final_cleanup_condition','final_cleanup_dna','final_cleanup_waste','final_cleanup'),g.authority)
 return g.require_publishable()

def render_page():
 from page import head,info
 from chemdraw import panel
 from protocol_views import assemble_page,WorkflowSection
 g=protocol_graph()
 return assemble_page(g,(head(TITLE),'<div class="wrap" data-protocol-complete="true">',f'<h1>{TITLE}</h1>',f'<p class="research-notes"><a href="{NOTES}">Research notes</a></p>',info(SOURCE),
  '<h2>Protect ordinary CpG cytosines before conversion</h2>',
  '<p>Genome-scale MAB-seq ligates methylated adapters, then uses two rounds of M.SssI with SAM to protect unmodified CpG cytosines. Bisulfite treatment leaves a combined5fC/5caC signal reported as T. Five-cycle preamplification precedes qPCR-guided final amplification.</p>',
  info(CAVEAT),
  '<h2>Physical products and sequencing calls</h2>',
  '<p>The panels show alternative source-defined CpG classes, not particular genomic sequences.5fC gives U and5caC gives5caU; both are reported as T after amplification. A T call alone cannot distinguish those starting modifications. Incomplete protection or conversion contributes background; the diagrams establish neither perfect conversion nor a quantitative assay result. Non-CpG cytosines are outside this protection table.</p>',
  '<h2>Known methylated adapter orders</h2>','<seq>'+''.join(oligos())+'</seq>',
  panel(adapter_scene().rows(),caption='Printed methylated custom Y-adapter pair: twelve-base common core,3-primeT end, forward terminal phosphorothioate and reverse5-prime phosphate. A per-position5mC map and occupied genomic junction are not invented.'),
  '<h2>Source and reaction boundaries</h2>',
  '<p>This view uses the explicitly identified associated2014 author manuscript and published supplementary controls, without asserting full final-version equivalence. It starts with unenriched genomic DNA; optional ChIP enrichment, locus-specific four-round MAB and caMAB are separate routes. Each protection round is a compound schedule:25µL for4h, then25µL fresh same-concentration enzyme/SAM and8h; purification follows each round. ReadyMix catalyst, nucleotide and buffer classes come from the same named physical supply.</p>',
  '<h2>Read setup and aliquots</h2>',
  '<p>The source reports HiSeq2500 single-end60or76bp. Historical indexed/universal primer members, complete amplification products, exact sequencing-primer molecules and acquired index mask remain unassembled. The qPCR aliquot stays separate from retained DNA; final library loading has an explicit untransferred stock remainder. No pool membership, read landing or index sequence is inferred from the kit name.</p>',
  WorkflowSection('Full preparation and molecular validation'),'</div>'))

def selftest():
 import re
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();src=Source(sources_root()/HERE.parent.name/'nbt.3073_PMC4269366.xml.txt','associated manuscript genome-scale Methods')
 if have(c,src,label='MAB custom adapter original'):
  t=re.sub(r'\s+','',src.path.read_text())
  c('forward source polymer with final phosphorothioate',ADAPTER_F[:-1]+'*'+ADAPTER_F[-1] in t,True)
  c('reverse source polymer with5Phos','/5Phos/'+ADAPTER_R in t,True)
 g=protocol_graph();render_page()
 c('SAM remains a molecular contribution in each source round',[any(p.material=='round'+str(n)+'_sam' and p.role is R.STRUCTURAL for p in g.steps['methylate'+str(n)].inputs) for n in (1,2)],[True,True])
 c('one actual ReadyMix supply creates its component inventory',[p.material for p in g.steps['kapa_components'].inputs],['kapa'])
 c('final amplification uses produced nucleotide inventory',any(p.material=='kapa_nt' for p in g.steps['final_pcr'].inputs),True)
 c('qPCR aliquot terminates separately',g.terminals['qpcr_aliquot'].disposition,Disposition.QC)
 c('sequencing uses selected aliquot',tuple(g.sequencing),('loading_aliquot',))
 c('untransferred library survives loading',g.terminals['library_remainder'].disposition,Disposition.UNRESOLVED)
 c.report()
