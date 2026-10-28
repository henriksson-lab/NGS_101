"""Original XR-seq, NHF1/CS-B lesion-IP order; source-qualified split-end ligation."""
from dataclasses import replace
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import illumina as il
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,feature,oligo,revcomp,complement_segments
from pcr import PCRGeometry,PCRPrimer
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole as R,MaterialKind as K,StepKind as S,ProtocolGraph,SequencingPlan,Terminal,Disposition,MolecularRelation
from protocol_model import molecule,role
PAGE_FILE='xr-seq.html'
TITLE='XR-seq — excised UV-repair oligonucleotides'
NOTES='01_xr-seq.html'
SOURCE='<a href="https://doi.org/10.1101/gad.261271.115">Hu et al., Genes & Development (2015)</a>; published Methods, Supplemental Material and Table5.'
SUMMARY='TFIIH-associated excision products receive two splinted DNA adapters. Sequential lesion-specific immunoprecipitation separates CPD and6-4 products; appropriate photolyase repairs the selected oligo before indexed KAPA PCR.'
CAVEAT='INFERRED — representative fragment-end and PCR-reference assemblies follow the printed oligos, with genomic sequence unassigned. Five random bases in each blocked splint are hybridization degeneracy, not UMIs. This graph follows NHF1/CS-B anti-CPD first, then anti-6-4 on flowthrough; XP-C reverses that order. Exact source-linked lesion-IP buffer recipe, sample index assignments and loaded run primers remain unreviewed. A nominal30mer or mean26nt is not a fixed insert-length constraint.'
A5F=il.TRUSEQ_SMALL_RNA_RA5
A5R='N'*5+revcomp(A5F)
A3F=il.TRUSEQ_SMALL_RNA_RA3
A3R=revcomp(A3F)+'N'*5
RP1=il.P5+A5F[:21]
RPIN_ARM='GTGACTGGAGTTC'+revcomp(A3F)[1:]
INDEX=feature('xr_i7','sample_index','unknown')
FORWARD=Construct([Segment('RP1',RP1)])
REVERSE=Construct([Segment('P7',il.P7),Segment('six-base sample index','N'*6,'cbc',placeholder=True,feature=INDEX),Segment('RPIn annealing arm',RPIN_ARM)])
TEMPLATE=Construct([Segment('A5F',A5F,'r1'),role('photorepaired excised DNA; genomic sequence unassigned'),Segment('A3F',A3F,'r2')],name='INFERRED — repaired adapter-ligated template')
GEOMETRY=PCRGeometry.exact(TEMPLATE,PCRPrimer('RP1',RP1,21),PCRPrimer('RPIn',REVERSE.top(),len(A3F)),statement='Published exact PCR polymers and A5F/A3F assembly; index member unassigned.')
# Preserve the entire unbounded genomic role semantically. Shared dna_window
# intentionally refuses coordinate slicing of prose roles; neither primer cuts it.
assert GEOMETRY.forward_hit.end<=len(A5F)
assert GEOMETRY.reverse_hit.start>=len(A5F)+len(TEMPLATE.segments[1])
_AMP=Construct((*FORWARD,Segment('retained A5F suffix',A5F[GEOMETRY.forward_hit.end:]),TEMPLATE.segments[1],Segment('retained A3F prefix',A3F[:GEOMETRY.reverse_hit.start-len(A5F)-len(TEMPLATE.segments[1])]),*complement_segments(tuple(REVERSE))),name='INFERRED — XR PCR reference')
FINAL_LIBRARY=Construct(tuple(replace(s,inferred=True,top=s.top.upper() if s.feature is not None else s.top) for s in _AMP),name='INFERRED — indexed XR reference library')
FINAL_CAPTION='INFERRED — PCR reference derived from full source primers and retained template interval; source sample index member and genomic sequence remain unassigned.'
SEQ_PRIMERS=()
SEQUENCING_UNAVAILABLE='Source HiSeq2000 single-end50; acquired index mask and exact loaded sequencing primers remain unassigned. Source-known PCR primers are not automatically treated as run primers.'
SEQUENCING_ENDING='Source collects1x50 reads; adapter-trimmed inserts report the repair-excised strand. Sequencing measures repair products rather than all genomic photolesions.'
READ_SPECS=(sp.read_spec('Read 1',50,'excised repair oligonucleotide and possible flanking adapter read-through','include'),sp.read_spec('Index 1 (i7)','unassigned','six-base molecular index; selected member and acquired cycles unassigned','undetermined'))

def adapter_pool():
 sc=Scene()
 sc.strand('a5f',[Segment('A5F',A5F)],label='A5F')
 sc.anneal('a5r',[Segment('5-end random splint','N'*5),Segment('A5R core',revcomp(A5F))],to='a5f',pair=('A5R core','A5F'),label='A5R',mod3='SpacerC3')
 sc.strand('a3f',[Segment('A3F',A3F)],label='A3F',mod5='5Phos',mod3='SpacerC3')
 sc.anneal('a3r',[Segment('A3R core',revcomp(A3F)),Segment('3-end random splint','N'*5)],to='a3f',pair=('A3R core','A3F'),label='A3R',mod3='SpacerC3')
 return sc

def ligation_geometry():
 five=Scene();five.strand('a5f',[Segment('A5F',A5F)],label='5-end adapter')
 five.anneal('a5r',[Segment('random5','N'*5),Segment('A5R core',revcomp(A5F))],to='a5f',pair=('A5R core','A5F'),mod3='SpacerC3',label='blocked splint')
 five.anneal('dna',[Segment('5-prime genomic end','N'*5),role('remaining excised DNA')],to='a5r',pair=('5-prime genomic end','random5'),label='representative excised5-end')
 five.ligation_end('a5f','A5F',kind='splinted',polarity='3-prime',owner='adapter')
 five.ligation_end('dna','5-prime genomic end',kind='splinted',polarity='5-prime',owner='substrate')
 three=Scene();three.strand('a3f',[Segment('A3F',A3F)],label='3-end adapter',mod5='5Phos',mod3='SpacerC3')
 three.anneal('a3r',[Segment('A3R core',revcomp(A3F)),Segment('random3','N'*5)],to='a3f',pair=('A3R core','A3F'),label='blocked splint',mod3='SpacerC3')
 three.anneal('dna',[role('remaining excised DNA'),Segment('3-prime genomic end','N'*5)],to='a3r',pair=('3-prime genomic end','random3'),label='representative excised3-end')
 three.ligation_end('a3f','A3F',kind='splinted',polarity='5-prime',owner='adapter')
 three.ligation_end('dna','3-prime genomic end',kind='splinted',polarity='3-prime',owner='substrate')
 return AdapterLigationGeometry.known('splinted',five,three,statement='Published Methods: A5F/A5R and A3F/A3R random5nt splinted end ligations, SpacerC3 blocks; two end-detail views of one captured oligomer.')

def oligos():
 return tuple(oligo(n,[Segment(n,s)],mods=m) for n,s,m in [('A5F',A5F,''),('A5R',A5R,'3-primeSpacerC3'),('A3F',A3F,'5Phos;3-primeSpacerC3'),('A3R',A3R,'3-primeSpacerC3'),('RP1',RP1,''),('RPIn family',REVERSE.top(),'sixN=selected sample index')])

def preparation_scene(dna_label,rna_label):
 sc=Scene()
 sc.strand('dna',[role(dna_label)],label='DNA species')
 sc.strand('rna',[role(rna_label)],label='separate RNA material')
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Hu2015 original XR-seq Methods; SupplementalMaterial/Table5; NHF1/CS-B route')
 cells=g.input('cells','NHF1 or CS-B cells after254nm20J/m2 exposure and source-specific recovery',K.BIOLOGICAL,state=molecule('repair-active cells',(role('TFIIH-associated excised UV-damaged DNA oligomers and other cellular material'),)),primary=True)
 lysis=g.input('lysis_condition','bufferA,Dounce lysis,chromatin centrifugation',K.CONDITION)
 soluble=g.product('soluble','soluble fraction containing excision complexes',K.MOLECULAR,state=cells.state)
 pellet=g.product('pellet','separated chromatin pellet',K.WASTE)
 g.react('lysis','Lyse and separate chromatin pellet from soluble material',S.EXTRACTION,((cells,R.TEMPLATE),(lysis,R.CONDITION)),(soluble,pellet));g.discard(pellet,'Only soluble fraction is carried into represented TFIIH capture; no chemical destruction asserted.')
 tf=g.input('tfii h_capture'.replace(' ',''),'anti-p89/XPB,anti-p62,A/Gagarose; RNaseA; source washes',K.SELECTION)
 bound=g.product('tf_bound','TFIIH immunoprecipitated material',K.MOLECULAR,state=molecule('TFIIH-selected pool',(role('protein-associated excision oligonucleotides; purity and conversion unassigned'),)))
 unbound=g.product('tf_unbound','TFIIH capture supernatants/washes',K.WASTE)
 g.react('tf_capture','Capture TFIIH-associated excision oligomers',S.AFFINITY_CAPTURE,((soluble,R.TEMPLATE),(tf,R.SELECTION)),(bound,unbound));g.discard(unbound,'Source capture/washing physically separates unretained material.')
 def recover(prefix,template,label,*,state=None,kind=K.MOLECULAR):
  cond=g.input(prefix+'_condition',label,K.CONDITION)
  retained=g.product(prefix+'_dna','recovered DNA',kind,state=template.state if state is None else state)
  waste=g.product(prefix+'_waste','removed purification fractions',K.WASTE)
  g.react(prefix,label,S.EXTRACTION,((template,R.TEMPLATE),(cond,R.CONDITION)),(retained,waste));g.discard(waste,'Source recovery removes separated fractions; no perfect recovery asserted.')
  return retained
 current=recover('tf_elution',bound,'65C SDS elution;PCI extraction/ethanol precipitation',state=preparation_scene('excised UV-damaged DNA oligomers','potential RNA carryover'))
 rnase=g.input('rnase','RNaseA/T1,37C1h',K.CATALYST)
 treated=g.product('rnase_treated','DNA plus RNA digestion products',K.MOLECULAR,state=preparation_scene('UV-damaged excised DNA','RNA digestion fragments remaining in tube'))
 g.react('digest_rna','Digest contaminating RNA before physical recovery',S.DIGESTION,((current,R.TEMPLATE),(rnase,R.CATALYST)),(treated,))
 current=recover('rna_cleanup',treated,'PCI/G50filtration/ethanol recovery',state=molecule('purified excision products',(role('UV-damaged excised DNA,unassigned sequence'),)))
 adapters=g.input('adapters','preannealed5-end and3-end splinted DNA adapters',K.ADAPTER,state=adapter_pool())
 ligase=g.input('ligase','T4DNA ligaseHC,PEG8000,overnight16C',K.CATALYST)
 ligated=g.product('ligated','INFERRED — adapter-ligated damaged excision oligomers',K.MOLECULAR,state=Construct((Segment('A5F',A5F),role('excised DNA containing CPD or6-4PP'),Segment('A3F',A3F)),name='INFERRED — ligated damaged strand'),inferred=True)
 g.react('ligate','INFERRED — ligate two end-specific splinted adapters',S.LIGATION,((current,R.TEMPLATE),(adapters,R.STRUCTURAL),(ligase,R.CATALYST)),(ligated,),geometry=ligation_geometry())
 current=recover('ligation_cleanup',ligated,'PCI extraction and ethanol precipitation')
 cpd_ab=g.input('cpd_antibody','anti-CPD antibody; source-linked IP procedure',K.SELECTION)
 cpd=g.product('cpd_selected','purified CPD-selected adapter-ligated DNA',K.MOLECULAR,state=current.state)
 flow=g.product('cpd_flowthrough','CPD-IP flowthrough for subsequent6-4PP capture',K.MOLECULAR,state=current.state)
 g.react('cpd_capture','Capture CPD; retain flowthrough for6-4PP IP',S.AFFINITY_CAPTURE,((current,R.TEMPLATE),(cpd_ab,R.SELECTION)),(cpd,flow))
 pp_ab=g.input('pp_antibody','anti-6-4PP antibody; source-linked IP procedure',K.SELECTION)
 pp=g.product('pp_selected','purified6-4PP-selected adapter-ligated DNA',K.MOLECULAR,state=current.state)
 remaining=g.product('pp_flowthrough','material not selected by subsequent6-4PP capture',K.MOLECULAR,state=current.state)
 g.react('pp_capture','Capture6-4PP from CPD-IP flowthrough',S.AFFINITY_CAPTURE,((flow,R.TEMPLATE),(pp_ab,R.SELECTION)),(pp,remaining))
 g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'Remaining unselected flowthrough is not a sequencing branch here; no physical disposal or complete depletion is asserted.'))
 forward=g.input('rp1','RP1 PCR primer',K.OLIGO,state=FORWARD)
 reverse=g.input('rpin','RPIn six-base index family',K.OLIGO,state=REVERSE)
 kapa=g.input('kapa','KAPA HiFi HotStart ReadyMix',K.CATALYST)
 for prefix,dna,enzyme in [('cpd',cpd,'20pmol Drosophila CPD photolyase'),('pp',pp,'1.25pmol E.coli6-4PP photolyase')]:
  photo=g.input(prefix+'_photolyase',enzyme,K.CATALYST)
  light=g.input(prefix+'_light','366nm18mW/cm2,2h on ice,sourcePLR buffer',K.CONDITION)
  repaired=g.product(prefix+'_repaired','INFERRED — photorepaired adapter-ligated oligomers',K.MOLECULAR,state=TEMPLATE,inferred=True)
  g.react(prefix+'_photorepair','INFERRED — lesion-specific photorepair before PCR',S.CHEMICAL_CONVERSION,((dna,R.TEMPLATE),(photo,R.CATALYST),(light,R.CONDITION)),(repaired,))
  qc_cond=g.input(prefix+'_qc_condition','1% repaired sample for quality check',K.CONDITION)
  qc=g.product(prefix+'_qc','repaired sample QC aliquot',K.QC,state=TEMPLATE)
  retained=g.product(prefix+'_for_pcr','remaining repaired material',K.MOLECULAR,state=TEMPLATE)
  g.react(prefix+'_qc_split','Reserve1% repaired material for QC',S.SPLIT,((repaired,R.TEMPLATE),(qc_cond,R.CONDITION)),(qc,retained));g.mark_terminal(Terminal(qc.id,Disposition.QC,'Source allows unrepaired or repaired aliquotQC; graph represents repaired-aliquot choice.'))
  purified=recover(prefix+'_repair_cleanup',retained,'PCI extraction/ethanol precipitation after photorepair')
  amplified=g.product(prefix+'_amplified','INFERRED — indexed XR PCR reference',K.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
  cycles=g.input(prefix+'_pcr_cycles','Table5 sample-specific6–15cycles; exact sample unassigned',K.CONDITION)
  g.react(prefix+'_pcr','INFERRED — source RP1/RPIn KAPA amplification',S.PCR,((purified,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(kapa,R.CATALYST),(cycles,R.CONDITION)),(amplified,),geometry=GEOMETRY,relations=(MolecularRelation.PCR_PRODUCT,))
  library=recover(prefix+'_gel',amplified,'PCI/ethanol;10%nativePAGE130–155nt band excision;NaCl elution/ethanol;PicoGreen',kind=K.LIBRARY)
  g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE,read_lengths={'Read 1':50},acquired_roles=('Read 1',)))
 g.declare_reviewed_scope(('cells','lysis_condition','soluble','pellet','lysis','tfiih_capture','tf_bound','tf_unbound','tf_capture','tf_elution_condition','tf_elution_dna','tf_elution_waste','tf_elution','rnase','rnase_treated','digest_rna','rna_cleanup_condition','rna_cleanup_dna','rna_cleanup_waste','rna_cleanup','adapters','ligase','ligated','ligate','ligation_cleanup_condition','ligation_cleanup_dna','ligation_cleanup_waste','ligation_cleanup','cpd_antibody','cpd_selected','cpd_flowthrough','cpd_capture','pp_antibody','pp_selected','pp_flowthrough','pp_capture','rp1','rpin','kapa','cpd_photolyase','cpd_light','cpd_repaired','cpd_photorepair','cpd_qc_condition','cpd_qc','cpd_for_pcr','cpd_qc_split','cpd_repair_cleanup_condition','cpd_repair_cleanup_dna','cpd_repair_cleanup_waste','cpd_repair_cleanup','cpd_amplified','cpd_pcr_cycles','cpd_pcr','cpd_gel_condition','cpd_gel_dna','cpd_gel_waste','cpd_gel','pp_photolyase','pp_light','pp_repaired','pp_photorepair','pp_qc_condition','pp_qc','pp_for_pcr','pp_qc_split','pp_repair_cleanup_condition','pp_repair_cleanup_dna','pp_repair_cleanup_waste','pp_repair_cleanup','pp_amplified','pp_pcr_cycles','pp_pcr','pp_gel_condition','pp_gel_dna','pp_gel_waste','pp_gel'),g.authority)
 return g.require_publishable()

def render_page():
 from batch_page import render
 return render(sys.modules[__name__])

def selftest():
 import re
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();source=Source(sources_root()/HERE.parent.name/'gad.261271.115_PMC4421983.xml.txt','published Oligonucleotides and adaptors Methods')
 if have(c,source,label='XR source oligos'):
  t=re.sub(r'\s+','',source.path.read_text())
  for label,sequence in [('A5F',A5F),('A5R',A5R+'-SpC3'),('A3F',A3F+'-SpC3'),('A3R',A3R+'-SpC3'),('RP1',RP1),('RPIn',il.P7+'X'*6+RPIN_ARM)]:c('source '+label,sequence in t,True)
 protocol_graph();render_page();c.report()
