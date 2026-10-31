"""Source-local AtRNL terminal chemistry; never pretends the substrate has3-OH."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import illumina as il
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,feature,oligo,circle_rows
from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import Disposition,Terminal,InputRole as R,MaterialKind as K,ProtocolGraph,SequencingPlan,StepKind as S
from protocol_model import molecule,role
PAGE_FILE='ribose-seq.html'
TITLE='Ribose-seq — cyclic-phosphate capture and ribonucleotide-bearing circles'
NOTES='01_ribose-seq.html'
SOURCE='<a href="https://doi.org/10.1038/nmeth.3259">Koh et al. (2015)</a>, associated manuscript PMC4686381; <a href="https://doi.org/10.1038/s41592-019-0505-9">2019 Addendum and corrected Supplementary Table8</a>.'
SUMMARY='Restriction fragments receive a UMI-bearing adapter before alkaline cleavage at embedded ribonucleotides. AtRNL closes compatible fragments; T5 removes unclosed DNA and Tpt1 removes the residual junction2-prime phosphate before PCR.'
CAVEAT='INFERRED — circle segment order is derived from the source adapters and reaction scheme; genomic intervals and ribonucleotide identity are unassigned. AtRNL substrate2-prime/3-prime chemistry is declared explicitly; no free3-primeOH circle helper is used. The2019 Addendum corrects upstream/downstream context plots, Table4 and reverse-strand offsets; those analysis claims are not taken from the uncorrected2015 version. Exact circular-PCR strand assembly and loaded run primers remain unassembled. Adaptor.S has20bases plus terminal amine; Am is not a21st nucleotide.'
L_FIXED='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGGGAGTTCAGACGTGTGCTCTTCCGATCTAGCCAGCGCAGACCGTGAGGT'
SHORT='CCTCACGGTCTGCGCTGGCT'
UMI=feature('ribose_umi','umi','random')
INDEX=feature('ribose_i7','sample_index','fixed')
INDEXES=('CGTGAT','ACATCG','GCCTAA','TGGTCA')
PCR2=il.NEBNEXT_UNIVERSAL_PRIMER
PCR1_ARM=il.TRUSEQ_READ2[:-1]
PCR1=il.P7+INDEXES[0]+PCR1_ARM
CIRCLE_REFERENCE=Construct((Segment('UMI','N'*8,'umi',placeholder=True,feature=UMI,inferred=True),Segment('long adapter',L_FIXED,'r1',inferred=True),role('genomic DNA upstream of captured rNMP'),role('retained terminal rNMP; base identity unassigned')),name='INFERRED — captured circle reference order')
@dataclass(frozen=True)
class RiboseState:
 topology:str
 terminal_chemistry:str
 junction_2prime_phosphate:bool=False
 linear_debris:bool=False
 def __post_init__(self):
  if self.topology not in {'linear','circle'}:raise ValueError('unknown Ribose topology')
  if self.topology=='linear' and self.junction_2prime_phosphate:raise ValueError('unclosed fragment cannot have a ligation-junction phosphate')
 def rows(self):
  if self.topology=='circle':
   rows=circle_rows(CIRCLE_REFERENCE,'AtRNL junction; '+('residual2-primeP' if self.junction_2prime_phosphate else '2-primeP removed by Tpt1'))
   if self.linear_debris:rows.extend(molecule('digestion products still in tube',(role('T5-generated linear DNA digestion fragments; not covalently incorporated in circle'),)).rows())
   return rows
  sc=Scene();sc.strand('linear',list(CIRCLE_REFERENCE),label='INFERRED — rNMP-bearing linear reference',mod5='5Phos',mod3=self.terminal_chemistry)
  sc.note('linear','Compatible rNMP end classes are distinct from nonproductive3-primeP and other linear backgrounds; exact fractions unassigned.')
  return sc.rows()
FINAL_LIBRARY=Construct((role('PCR2-derived end; circular amplification assembly unassigned'),Segment('eight UMI cycles','N'*8,'umi',placeholder=True,feature=UMI),role('captured rNMP and upstream DNA, reported as reverse-complement sequence'),Segment('selected sample index',INDEXES[0],'cbc',feature=INDEX),role('PCR1-derived opposite end; exact strand assembly unassigned')),name='Ribose library roles; printed index1 example')
FINAL_CAPTION='Role-level library with source-established eight-cycle UMI and printed Index1 example; full circular-PCR product is not declared sequence-exact.'
SEQ_PRIMERS=()
SEQUENCING_UNAVAILABLE='MiSeq50-cycle single-end run; exact acquired index mask and loaded run-primer oligos are not assigned. PCR primers are source-known, but full circular-template PCR geometry remains unassembled.'
SEQUENCING_ENDING='First8raw read cycles are UMI; remove them before alignment. The captured rNMP identity is the reverse complement of the first aligned base. The2019 Addendum changes flanking-context directionality, not this physical capture order.'
READ_SPECS=(sp.read_spec('Read 1',50,'eight UMI cycles followed by reverse-complement genomic capture sequence','extract_then_include'),sp.read_spec('Index 1 (i7)','unassigned','printed six-base PCR1 index family; acquired mask unassigned','undetermined'))
SOURCE_OLIGO_EXCEPTIONS={'short_adapter':'Rendered corrected Table8 has20-base Adaptor.S followed by-Am(amine). Shared scraper incorrectly consumes the A in Am as base21; local source test checks exact20-base-Am notation and long/PCR polymers.'}

def adapter_scene():
 sc=Scene()
 sc.strand('long',[Segment('UMI','N'*8,'umi',placeholder=True,feature=UMI),Segment('unpaired long adapter',L_FIXED[:-21]),Segment('paired adapter core',L_FIXED[-21:-1]),Segment('terminal T','T')],label='Adaptor.L',mod5='5Phos')
 sc.anneal('short',[Segment('short adapter',SHORT)],to='long',pair=('short adapter','paired adapter core'),label='Adaptor.S',mod5='5Phos',mod3='amine')
 sc.ligation_end('long','terminal T',kind='cohesive',polarity='3-prime',owner='adapter')
 return sc

def oligos():
 rows=[oligo('Adaptor.L',[Segment('UMI','N'*8,'umi',placeholder=True,feature=UMI),Segment('long fixed sequence',L_FIXED)],mods='5Phos'),oligo('Adaptor.S',[Segment('short adapter',SHORT)],mods='5Phos;3-prime amine')]
 rows.extend(oligo('PCR.1.Index'+str(i+1),[Segment('P7',il.P7),Segment('printed index',x,'cbc',feature=feature('ribose_i7','sample_index','fixed')),Segment('PCR1 arm',PCR1_ARM)]) for i,x in enumerate(INDEXES))
 return (*rows,oligo('PCR.2',[Segment('PCR2',PCR2)]))

def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Koh2015 associated manuscript library Methods and sequencing; corrected2019 Table8; Addendum',representation='published_instance')
 dna=g.input('dna','yeast genomic DNA with embedded ribonucleotides',K.MOLECULAR,state=molecule('genomic substrate',(role('DNA with embedded rNMP sites; positions unassigned'),)),primary=True)
 restriction=g.input('restriction_enzymes','SspI,DraI,EcoRV; overnight37C',K.CATALYST)
 fragments=g.product('fragments','restriction fragments500–3000bp',K.MOLECULAR,state=dna.state)
 g.react('restrict','Restriction-fragment genomic DNA before end labeling',S.DIGESTION,((dna,R.TEMPLATE),(restriction,R.CATALYST)),(fragments,))
 klenow=g.input('klenow','Klenow exo-minus;37C30min',K.CATALYST)
 datp=g.input('datp','Source dATP for A-tailing; concentration unassigned',K.MOLECULAR)
 tailed=g.product('a_tailed','A-tailed genomic fragments',K.MOLECULAR,state=molecule('A-tailed DNA',(role('genomic fragment with source A-tailed ends'),)))
 g.react('a_tail','A-tail restriction fragments',S.TAILING,((fragments,R.TEMPLATE),(klenow,R.CATALYST),(datp,R.STRUCTURAL)),(tailed,))
 def clean(prefix,template,label,*,kind=K.MOLECULAR,retained_state=None):
  reagent=g.input(prefix+'_condition',label,K.CONDITION)
  retained=g.product(prefix+'_retained','purified DNA',kind,state=template.state if retained_state is None else retained_state)
  waste=g.product(prefix+'_waste','purification removed fractions',K.WASTE)
  g.react(prefix,label,S.CLEANUP,((template,R.TEMPLATE),(reagent,R.CONDITION)),(retained,waste));g.discard(waste,'Source purification removes supernatants/washes; exact recovery unassigned.')
  return retained
 current=clean('a_tail_cleanup',tailed,'QIAGEN spin-column cleanup')
 long_input=g.input('long_adapter','5-phosphorylated Adaptor.L',K.OLIGO,state=Construct((Segment('UMI','N'*8,'umi',placeholder=True,feature=UMI),Segment('long fixed sequence',L_FIXED))))
 short_input=g.input('short_adapter','Adaptor.S;5-phosphate,3-amine',K.OLIGO,state=Construct((Segment('short20nt',SHORT),)))
 adapter=g.product('adapters','preannealed Adaptor.L/Adaptor.S',K.ADAPTER,state=adapter_scene())
 g.react('anneal_adapters','Preanneal source Adaptor.L/Adaptor.S',S.HYBRIDIZATION,((long_input,R.TEMPLATE),(short_input,R.STRUCTURAL)),(adapter,))
 ligase=g.input('t4_ligase','T4DNA ligase;overnight15C',K.CATALYST)
 ligated=g.product('adapter_ligated','UMI-adapter-ligated genomic fragments',K.MOLECULAR,state=molecule('adapter-bearing substrate',(role('source adapter-bearing genomic duplex with embedded rNMP'),)))
 g.react('adapter_ligation','Attach UMI-bearing T-overhang adapter',S.LIGATION,((current,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.known('cohesive',adapter_scene(),statement='CorrectedTable8 and library Methods define20-bp core and3-primeT opposite source A-tail.'))
 current=clean('adapter_cleanup',ligated,'AMPure purification')
 alkali=g.input('alkali','0.3MNaOH,55C2h;neutralize and purify',K.CONDITION)
 cleaved=g.product('cleaved','rNMP-ended single-stranded fragments',K.MOLECULAR,state=RiboseState('linear','mixture:2-prime,3-primecyclicP /2-primeP /nonproductive3-primeP'))
 g.react('alkaline_cleave','Denature and cleave at embedded ribonucleotides',S.CHEMICAL_CONVERSION,((current,R.TEMPLATE),(alkali,R.CONDITION)),(cleaved,))
 current=clean('alkali_cleanup',cleaved,'Neutralization and AMPure purification')
 atrnl=g.input('atrnl','AtRNL1µM;30C1h',K.CATALYST)
 atp=g.input('atp','Source ATP30µM in AtRNL reaction',K.MOLECULAR)
 circles=g.product('circles','AtRNL-closed rNMP-bearing circles with residual2-primeP',K.MOLECULAR,state=RiboseState('circle','closed',True))
 linear=g.product('unclosed','unclosed/nonproductive linear DNA',K.MOLECULAR,state=cleaved.state)
 g.react('close','AtRNL intramolecular closure of compatible end classes',S.CIRCULARIZATION,((current,R.TEMPLATE),(atrnl,R.CATALYST),(atp,R.STRUCTURAL)),(circles,linear),geometry=CircularizationGeometry.known('single-strand intramolecular ligation',statement='Source AtRNL uses ribonucleotide2-prime,3-primecyclicP/2-primeP substrates and5-primeP; no claim of a pre-existing3-primeOH terminus.'))
 cleanup_condition=g.input('post_atrnl_condition','AMPure purification following AtRNL',K.CONDITION)
 clean_circle=g.product('post_atrnl_circles','recovered circles after AtRNL cleanup',K.MOLECULAR,state=circles.state)
 clean_linear=g.product('post_atrnl_linear','remaining unclosed linear DNA after AtRNL cleanup',K.MOLECULAR,state=linear.state)
 clean_waste=g.product('post_atrnl_waste','AtRNL cleanup removed fractions',K.WASTE)
 g.react('post_atrnl_cleanup','AMPure purification after AtRNL; remaining linears retained for T5 treatment',S.CLEANUP,((circles,R.TEMPLATE),(linear,R.TEMPLATE),(cleanup_condition,R.CONDITION)),(clean_circle,clean_linear,clean_waste));g.discard(clean_waste,'Source AtRNL purification removes wash/supernatant fractions; does not establish complete removal of linear DNA.')
 t5=g.input('t5','T5exonuclease,37C2h',K.CATALYST)
 retained=g.product('exonuclease_resistant','retained closed circles',K.MOLECULAR,state=circles.state)
 degraded=g.product('degraded','T5-digested linear products remaining in reaction tube',K.MOLECULAR,state=molecule('digestion products',(role('linear DNA digestion fragments; not physically removed by catalysis'),)))
 g.react('digest_linear','Degrade unclosed linear DNA while retaining circles',S.DIGESTION,((clean_circle,R.TEMPLATE),(clean_linear,R.TEMPLATE),(t5,R.CATALYST)),(retained,degraded))
 tpt=g.input('tpt1','Tpt1 1µM;30C1h',K.CATALYST)
 nad=g.input('nad','Source NAD+10mM in Tpt1 reaction',K.MOLECULAR)
 dephos=g.product('dephosphorylated','closed circle with junction2-primeP removed',K.MOLECULAR,state=RiboseState('circle','closed',False,True))
 g.react('remove_2prime_phosphate','Tpt1 removes residual junction2-primephosphate',S.CHEMICAL_CONVERSION,((retained,R.TEMPLATE),(degraded,R.TEMPLATE),(tpt,R.CATALYST),(nad,R.STRUCTURAL)),(dephos,))
 current=clean('circle_cleanup',dephos,'AMPure purification after Tpt1; remove digestion fragments in discarded fractions',retained_state=RiboseState('circle','closed',False))
 p1=g.input('pcr1','PCR.1.Index1 example; four source index choices',K.OLIGO,state=Construct([Segment('P7',il.P7),Segment('Index1',INDEXES[0],'cbc',feature=INDEX),Segment('PCR1 arm',PCR1_ARM)]))
 p2=g.input('pcr2','PCR.2',K.OLIGO,state=Construct([Segment('PCR2',PCR2)]))
 polymerase=g.input('polymerase','source alternative:Phusion or EconoTaq',K.CATALYST)
 amplified=g.product('amplified','circle-derived PCR library',K.MOLECULAR,state=FINAL_LIBRARY)
 g.react('pcr','Amplify with PCR1/PCR2; exact circular product geometry unassembled',S.PCR,((current,R.TEMPLATE),(p1,R.STRUCTURAL),(p2,R.STRUCTURAL),(polymerase,R.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(molecule('circular template',(role('rNMP-bearing source circle with both printed primer sites; exact strand assembly unassigned'),)),reason='Printed primers known; full circular-template linear amplicon derivation not completed in capped model. No free3-primeOH constructor is substituted.'))
 page_aliquot=g.product('page_aliquot','Analytical aliquot for6%PAGE confirmation; amount unassigned',K.QC,state=amplified.state)
 retained_pcr=g.product('retained_pcr','PCR material retained for library purification',K.MOLECULAR,state=amplified.state)
 g.react('split_page','Reserve analytical PAGE aliquot and retain library material',S.SPLIT,((amplified,R.TEMPLATE),),(page_aliquot,retained_pcr),changes_molecule=False)
 g.mark_terminal(Terminal(page_aliquot.id,Disposition.QC,'Source PAGE confirms amplified DNA; analytical amount and gel recovery fraction unassigned, not asserted to be returned to the library.'))
 library=clean('final_cleanup',retained_pcr,'Post-PCR purification; source recovery unassigned',kind=K.LIBRARY)
 contribution=g.product('pool_contribution','Selected equimolar contribution of this library',K.LIBRARY,state=library.state)
 library_remainder=g.product('library_remainder','Untransferred purified library stock; amount unassigned',K.LIBRARY,state=library.state)
 g.react('select_pool_contribution','Select contribution for equimolar multi-library pool',S.SPLIT,((library,R.TEMPLATE),),(contribution,library_remainder),changes_molecule=False)
 g.mark_terminal(Terminal(library_remainder.id,Disposition.UNRESOLVED,'Untransferred library stock; selected contribution and residual amounts unassigned.'))
 others=g.input('other_libraries','Other indexed libraries for the source equimolar pool; members/counts unassigned',K.LIBRARY)
 pooled=g.product('pooled','Equimolar sequencing library pool normalized10nM',K.LIBRARY,state=FINAL_LIBRARY)
 g.react('pool','Mix distinct indexed libraries; normalize10nM',S.POOL,((contribution,R.TEMPLATE),(others,R.TEMPLATE)),(pooled,),changes_molecule=False)
 loading=g.product('loading_aliquot','Selected sequencing loading aliquot; volume unassigned',K.LIBRARY,state=pooled.state)
 pool_remainder=g.product('pool_remainder','Untransferred sequencing-pool stock; amount unassigned',K.LIBRARY,state=pooled.state)
 g.react('select_loading','Select aliquot from sequencing pool',S.SPLIT,((pooled,R.TEMPLATE),),(loading,pool_remainder),changes_molecule=False)
 g.mark_terminal(Terminal(pool_remainder.id,Disposition.UNRESOLVED,'The source reports pooled sequencing, not consumption of the entire physical stock. Loading and residual volumes unassigned.'))
 g.sequence(loading,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE,read_lengths={'Read 1':50},acquired_roles=('Read 1',)))
 g.declare_reviewed_scope(('datp','atp','nad','split_page','page_aliquot','select_pool_contribution','library_remainder','other_libraries','pool','select_loading','pool_remainder','loading_aliquot','dna','restriction_enzymes','fragments','restrict','klenow','a_tailed','a_tail','a_tail_cleanup_condition','a_tail_cleanup_retained','a_tail_cleanup_waste','a_tail_cleanup','long_adapter','short_adapter','adapters','anneal_adapters','t4_ligase','adapter_ligated','adapter_ligation','adapter_cleanup_condition','adapter_cleanup_retained','adapter_cleanup_waste','adapter_cleanup','alkali','cleaved','alkaline_cleave','alkali_cleanup_condition','alkali_cleanup_retained','alkali_cleanup_waste','alkali_cleanup','atrnl','circles','unclosed','close','post_atrnl_condition','post_atrnl_circles','post_atrnl_linear','post_atrnl_waste','post_atrnl_cleanup','t5','exonuclease_resistant','degraded','digest_linear','tpt1','dephosphorylated','remove_2prime_phosphate','circle_cleanup_condition','circle_cleanup_retained','circle_cleanup_waste','circle_cleanup','pcr1','pcr2','polymerase','amplified','pcr','final_cleanup_condition','final_cleanup_retained','final_cleanup_waste','final_cleanup'),g.authority)
 return g.require_publishable()

def render_page():
 from page import head,info
 from chemdraw import panel
 from protocol_views import assemble_page,WorkflowSection
 g=protocol_graph()
 return assemble_page(g,(head(TITLE),'<div class="wrap" data-protocol-complete="true">',f'<h1>{TITLE}</h1>',f'<p class="research-notes"><a href="{NOTES}">Research notes</a></p>',info(SOURCE),
  '<h2>Capture the DNA next to an embedded ribonucleotide</h2>',
  '<p>Ribose-seq attaches UMI-bearing adapters to restriction fragments, then cleaves DNA at embedded ribonucleotides with alkali. AtRNL closes compatible single-stranded fragments. T5 digests unclosed DNA; Tpt1 removes the residual junction phosphate before PCR.</p>',
  info(CAVEAT),
  '<h2>Source revision and molecular scope</h2>',
  '<p>The preparation is declared from the associated 2015 author manuscript and the corrected publisher supplement. The 2019 Addendum corrects upstream/downstream sequence-context plots, Table4 and reverse-strand analysis offsets. The uncorrected directional claims are excluded; complete accepted/final-version equivalence is not asserted. The addendum text rechecked in this pass is an associated manuscript.</p>',
  '<p>AtRNL-compatible cyclic-phosphate/2-prime-phosphate end classes and the residual junction2-prime phosphate remain distinct from a free3-primeOH substrate. The exact circular-PCR product and genomic coordinates are unassembled. Index1 is a printed primer reference choice, not a claim that every experimental sample used it.</p>',
  '<h2>Known adapter and primer orders</h2>','<seq>'+''.join(oligos())+'</seq>',
  panel(adapter_scene().rows(),caption='Source-exact adapter common-arm pairing and3-primeT overhang. Adaptor.S has twenty nucleotides and a non-nucleotide3-prime amine; the source5-prime phosphates are retained.'),
  '<h2>Sequencing and stock boundaries</h2>',
  '<p>MiSeq single-end50 cycles: the first eight raw cycles are the UMI. The six-base printed index family is known; the acquired index mask and exact loaded sequencing-primer molecules remain unassigned. PAGE analytical material is separated from retained library material. Equimolar pool contributions and final loading aliquots have explicit untransferred remainders; amounts and other pool members are unassigned.</p>',
  WorkflowSection('Full preparation and molecular validation'),'</div>'))

def selftest():
 import re
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();src=Source(sources_root()/'doi_10.1038+s41592-019-0505-9'/'s41592-019-0505-9_41592_2019_505_MOESM1_ESM.pdf.txt','corrected2019Table8, rendered page20 verified')
 if have(c,src,label='correctedTable8'):
  t=re.sub(r'\s+','',src.path.read_text())
  for label,seq in [('Adaptor.L','P-'+'N'*8+L_FIXED),('Adaptor.S20nt amine','P-'+SHORT+'-Am'),('PCR2',PCR2)]+[('PCR1.Index'+str(i+1),il.P7+x+PCR1_ARM) for i,x in enumerate(INDEXES)]:c('source '+label,seq in t,True)
 g=protocol_graph();render_page()
 c('source molecular cofactors remain explicit inputs',[(k,g.materials[k].kind) for k in ('datp','atp','nad')],[(k,K.MOLECULAR) for k in ('datp','atp','nad')])
 c('Tpt1 retains molecular NAD contribution',any(p.material=='nad' and p.role is R.STRUCTURAL for p in g.steps['remove_2prime_phosphate'].inputs),True)
 c('PAGE material terminates analytically',g.terminals['page_aliquot'].disposition,Disposition.QC)
 c('sequencing consumes a selected loading aliquot',tuple(g.sequencing),('loading_aliquot',))
 c('untransferred pool survives as an explicit terminal',g.terminals['pool_remainder'].disposition,Disposition.UNRESOLVED)
 c('selected pool is a physical combination, not a concatenated molecule',g.steps['pool'].kind,S.POOL)
 c.report()
