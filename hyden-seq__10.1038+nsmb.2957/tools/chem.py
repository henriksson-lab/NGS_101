"""Original alkaline HydEn-seq; 3-prime-phosphatase-minus PNK is explicit."""
from dataclasses import replace
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import illumina as il
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,feature,oligo,revcomp,complement_segments
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry,PCRPrimer
from protocol_graph import ProtocolGraph,InputRole as R,MaterialKind as K,StepKind as S,SequencingPlan,MolecularRelation
from protocol_model import molecule,role
PAGE_FILE='hyden-seq.html';TITLE='HydEn-seq — hydrolytic5-prime end capture';NOTES='01_hyden-seq.html'
SOURCE='<a href="https://doi.org/10.1038/nsmb.2957">Clausen et al.(2015)</a>, associated manuscript PMC4351163 and publisher SupplementaryTable1,renderedpage4.'
SUMMARY='Alkaline hydrolysis creates5-primeOH DNA ends.3-prime-phosphatase-minusPNK phosphorylates them; T4RNA ligase attaches a5-amine-blocked adapter. Random-hexamer-adapter/T7 synthesis and indexedPCR produce the library.'
CAVEAT='INFERRED — fragment assemblies are derived from source oligos; genomic sequence and random priming positions remain unassigned. Source Table1 labels first indexed row ARC78, whereas associated-manuscript Methods names ARC79; this graph uses unambiguous ARC84 example ACATCG. Index lengths vary across source members. ARC76random6 is priming degeneracy,notUMI. Native ends and background can contribute; an aligned end is not automatically a newly incorporated rNMP. Complete publishedVOR and source-linked original details remain unreviewed.'
A=il.TRUSEQ_READ1;B=il.TRUSEQ_READ2
ARC76=B+'N'*6;ARC77=revcomp(B);P5=il.NEBNEXT_UNIVERSAL_PRIMER;I7='ACATCG';P7=il.P7+I7+B
F=Construct([Segment('ARC49',P5)])
RPR=Construct([Segment('P7',il.P7),Segment('ARC84index',I7,'cbc',feature=feature('hyden_i7','sample_index','fixed')),Segment('ARC84arm',B)])
TEMPLATE=Construct([Segment('ARC140',A),role('DNA downstream of hydrolyzed rNMP; random-primer extent unassigned'),Segment('opposite adapter arm',revcomp(B))],name='INFERRED — second-strand template')
GEO=PCRGeometry.exact(TEMPLATE,PCRPrimer('ARC49',P5,len(A)),PCRPrimer('ARC84',P7,len(B)),statement='Source primer polymers and locally derived entry arms; entire genomic role preserved.')
assert GEO.forward_hit.end==len(A) and GEO.reverse_hit.start==len(A)+len(TEMPLATE.segments[1])
FINAL_LIBRARY=Construct(tuple(replace(s,inferred=True,top=s.top.upper() if s.feature else s.top) for s in (*F,TEMPLATE.segments[1],*complement_segments(tuple(RPR)))),name='INFERRED — ARC84 HydEn reference')
FINAL_CAPTION='INFERRED — source ARC84 example, whole genomic interval retained between exact primer arms. No source study-index assignment or universal fragment length implied.'
SEQ_PRIMERS=();SEQUENCING_UNAVAILABLE='Source HiSeq2500 paired-end50; exact loaded run primers/index mask and selected experiment assignment remain unassigned.'
SEQUENCING_ENDING='Read5-prime ends report captured DNA downstream of hydrolysis. Source shifts counts one base upstream to infer rNMP coordinates; computational shift does not change molecular bases.'
READ_SPECS=(sp.read_spec('Read 1',50,'captured downstream DNA5-prime boundary','include'),sp.read_spec('Read 2',50,'opposite fragment/priming end','include'),sp.read_spec('Index 1 (i7)','unassigned','ARC84molecular six-base example ACATCG; acquiredmaskunassigned','undetermined'))
def end_state(label,*,rna=False):
 sc=Scene();sc.strand('dna',[role('genomic DNA fragment; end position unassigned')],label=label)
 if rna:sc.strand('rna',[role('possible RNA carryover/hydrolysis fragments; fate quantified neither here nor by DNA diagram')],label='separate RNA material')
 return sc

def adapter_pool():
 sc=Scene();sc.strand('arc76',[Segment('ARC76arm',B),Segment('random6','N'*6)],label='ARC76;last2bondsphosphorothioate')
 sc.anneal('arc77',[Segment('ARC77',ARC77)],to='arc76',pair=('ARC77','ARC76arm'),label='ARC77;last2bondsphosphorothioate')
 return sc

def geometry():
 sc=Scene();sc.strand('arc140',[Segment('ARC140',A)],label='ARC140',mod5='AmMC6')
 sc.strand('dna',[role('DNA downstream of hydrolysis')],label='phosphorylated DNA',mod5='5Phos',mod3='phosphate-ended class retained')
 sc.ligation_end('arc140','ARC140',kind='single-strand',polarity='3-prime',owner='adapter')
 sc.ligation_end('dna','DNA downstream of hydrolysis',kind='single-strand',polarity='5-prime',owner='substrate')
 return AdapterLigationGeometry.known('single-strand',sc,statement='Source ARC1403-OH accepts DNA5-P after3-phosphatase-minusPNK; ARC1405-amine prevents a free5-OH reagent end.')

def oligos():return (oligo('ARC140',[Segment('ARC140',A)],mods='5AmMC6'),oligo('ARC76',[Segment('ARC76',ARC76)],mods='last2phosphorothioatebonds'),oligo('ARC77',[Segment('ARC77',ARC77)],mods='last2phosphorothioatebonds'),oligo('ARC49',[Segment('ARC49',P5)]),oligo('ARC84',[Segment('ARC84',P7)]))
def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Original HydEn-seq Methods; publisherTable1,selectedARC84')
 dna=g.input('dna','MasterPure-extracted yeast DNA;noRNaseA treatment',K.MOLECULAR,state=end_state('starting preparation',rna=True),primary=True)
 def op(i,label,kind,substrate,reagent,state,*,rk=K.CONDITION,rr=R.CONDITION):
  reagent=g.input(i+'_reagent',reagent,rk)
  out=g.product(i+'_product',label,K.MOLECULAR,state=state)
  g.react(i,label,kind,((substrate,R.TEMPLATE),(reagent,rr)),(out,));return out
 def clean(i,substrate,label,*,state=None,kind=K.MOLECULAR):
  reagent=g.input(i+'_condition',label,K.CONDITION);out=g.product(i+'_dna','recovered DNA',kind,state=substrate.state if state is None else state);w=g.product(i+'_waste','removed purification fractions',K.WASTE)
  g.react(i,label,S.CLEANUP,((substrate,R.TEMPLATE),(reagent,R.CONDITION)),(out,w));g.discard(w,'Physical source purification; no quantitative purity/recovery guarantee.');return out
 hydro=op('hydrolyze','KOHhydrolysis0.3M55C2h',S.CHEMICAL_CONVERSION,dna,'0.3MKOH,55C2h',end_state('DNA5-OH and blocked3-phosphate classes after hydrolysis',rna=True))
 current=clean('ethanol',hydro,'ethanol precipitation',state=end_state('hydrolyzed DNA target class; RNA-removal completeness unassigned'))
 current=op('heat1','Heat85C3min',S.THERMAL_TREATMENT,current,'85C3min',current.state)
 current=op('phosphorylate','Phosphorylate5-OH; preserve3-phosphate class',S.PHOSPHORYLATION,current,'10U3-prime-phosphatase-minusT4PNK,37C30min;65C20min',end_state('DNA5-P; hydrolysis3-phosphate class preserved'),rk=K.CATALYST,rr=R.CATALYST)
 current=clean('pnk_cleanup',current,'HighPrepPCRbeads')
 current=op('heat2','Heat85C3min beforeRNA-ligase capture',S.THERMAL_TREATMENT,current,'85C3min',current.state)
 a=g.input('arc140','ARC1405-amine/C6 adapter',K.ADAPTER,state=Construct([Segment('ARC140',A)]))
 ligase=g.input('rna_ligase','T4RNAligase10U;25%PEG8000,1mMCoCl3(NH3)6;overnightroomtemperature',K.CATALYST)
 joined_state=Scene();joined_state.strand('joined',[Segment('ARC140',A),role('captured DNA downstream ofrNMP')],label='source sealed adapter-DNA reference',mod5='AmMC6');joined_state.note('joined','Junction:ARC1403-OH toDNA5-P; native DNA5end captured after PNK.')
 joined=g.product('joined','ARC140-ligated target DNA',K.MOLECULAR,state=joined_state)
 g.react('capture','LigateARC140toDNA5-primeP',S.LIGATION,((current,R.TEMPLATE),(a,R.STRUCTURAL),(ligase,R.CATALYST)),(joined,),geometry=geometry())
 current=clean('ligation_cleanup',joined,'HighPrepPCRbeads')
 current=op('heat3','Heat85C3min before second-strand adapter',S.THERMAL_TREATMENT,current,'85C3min',current.state)
 pool=g.input('second_strand_adapter','preannealedARC76/ARC77',K.ADAPTER,state=adapter_pool())
 primed=g.product('primed','random-primed genomic template class;coordinatesunassigned',K.MOLECULAR,state=end_state('ARC76random6priming class;exactsiteunassigned'))
 g.react('anneal','AnnealARC76/ARC77for5minroomtemperature',S.HYBRIDIZATION,((current,R.TEMPLATE),(pool,R.STRUCTURAL)),(primed,))
 t7=g.input('t7','T7DNApolymerase4U',K.CATALYST)
 filled=g.product('filled','INFERRED — second-strand end architecture',K.MOLECULAR,state=TEMPLATE,inferred=True)
 g.react('second_strand','INFERRED — T7second-strand synthesis fromsource adapter',S.EXTENSION,((primed,R.TEMPLATE),(pool,R.STRUCTURAL),(t7,R.CATALYST)),(filled,))
 current=clean('second_cleanup',filled,'HighPrepPCRbeads')
 f=g.input('arc49','ARC49PCRprimer',K.OLIGO,state=F);r=g.input('arc84','ARC84PCRprimer example',K.OLIGO,state=RPR);kapa=g.input('kapa','KAPAHiFiHotstartReadyMix;cyclesunassigned',K.CATALYST)
 amp=g.product('amp','INFERRED — ARC84PCRreference',K.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
 g.react('pcr','INFERRED — indexedPCRwithARC49/ARC84example',S.PCR,((current,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(kapa,R.CATALYST)),(amp,),geometry=GEO,relations=(MolecularRelation.PCR_PRODUCT,))
 library=clean('final_cleanup',amp,'HighPrepPCRbeads;poollibraries',kind=K.LIBRARY)
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE,read_lengths={'Read 1':50,'Read 2':50},acquired_roles=('Read 1','Read 2')))
 g.declare_reviewed_scope(('dna','hydrolyze_reagent','hydrolyze_product','hydrolyze','ethanol_condition','ethanol_dna','ethanol_waste','ethanol','heat1_reagent','heat1_product','heat1','phosphorylate_reagent','phosphorylate_product','phosphorylate','pnk_cleanup_condition','pnk_cleanup_dna','pnk_cleanup_waste','pnk_cleanup','heat2_reagent','heat2_product','heat2','arc140','rna_ligase','joined','capture','ligation_cleanup_condition','ligation_cleanup_dna','ligation_cleanup_waste','ligation_cleanup','heat3_reagent','heat3_product','heat3','second_strand_adapter','primed','anneal','t7','filled','second_strand','second_cleanup_condition','second_cleanup_dna','second_cleanup_waste','second_cleanup','arc49','arc84','kapa','amp','pcr','final_cleanup_condition','final_cleanup_dna','final_cleanup_waste','final_cleanup'),g.authority)
 return g.require_publishable()
def render_page():
 from batch_page import render
 return render(sys.modules[__name__])
def selftest():
 import re
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();source=Source(sources_root()/HERE.parent.name/'nsmb.2957_41594_2015_BFnsmb2957_MOESM6_ESM.pdf.txt','publisherTable1renderedpage4')
 if have(c,source,label='HydEnTable1'):
  t=re.sub(r'\s+','',source.path.read_text())
  for n,x in [('ARC140','/5AmMC6/'+A),('ARC76',B+'NNNN*N*N'),('ARC77',ARC77[:-2]+'*'+ARC77[-2]+'*'+ARC77[-1]),('ARC49',P5),('ARC84',P7)]:c('source '+n,x in t,True)
 protocol_graph();render_page();c.report()
