"""GLOE core route: native 3-prime end ligation precedes fragmentation."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import MolecularRelation, InputRole as R, MaterialKind as K, ProtocolGraph, SequencingPlan, StepKind as S
from protocol_model import molecule, role
PAGE_FILE='gloe-seq.html'
TITLE='GLOE-Seq — capture native 3-prime hydroxyl ends before fragmentation'
NOTES='01_gloe-seq.html'
SOURCE='<a href="https://doi.org/10.1016/j.molcel.2020.03.027">Sriramachandran et al. (2020)</a>, published STAR Methods, MethodsS1 and rendered FigureS1.'
SUMMARY='A random splinter bridges a native DNA 3-prime hydroxyl and a phosphorylated biotinylated adapter. Ligation precedes sonication; streptavidin enrichment, complementary-strand synthesis, distal adapter ligation and PCR form a strand-specific library.'
CAVEAT='INFERRED — strand assemblies are derived from the printed oligos and FigureS1. Genomic sequence and the selected six-base P7 index are unassigned. This page follows the native-end core route, not a lesion-specific enzyme treatment or an assigned dual-index run. Heat release from streptavidin is represented operationally; the specific covalent biotin-release chemistry is unassigned.'
SPLINT_CORE=il.TRUSEQ_READ1[-22:]
SPLINT=SPLINT_CORE+'N'*6
SECOND=il.NEBNEXT_UNIVERSAL_PRIMER[17:]
PROXIMAL=revcomp(SECOND)+'TTTT'
DISTAL1=il.TRUSEQ_READ2[2:]
DISTAL2=il.NEBNEXT_ARM_READ2
P5=il.NEBNEXT_UNIVERSAL_PRIMER
P7_ARM=il.TRUSEQ_READ2
P7=il.P7+'N'*6+P7_ARM
SAMPLE_INDEX=feature('gloe_i7','sample_index','unknown')
CAPTURED=Construct((role('original captured DNA strand; native 3-prime end at adapter junction'),Segment('proximal adapter',PROXIMAL,'r1')),name='source captured single strand')
TEMPLATE=Construct((Segment('second-strand primer',SECOND,'r1'),role('reverse complement of captured genomic fragment'),Segment('distal arm reverse complement',revcomp(DISTAL1),'r2')),name='INFERRED — A/T-ligated pre-PCR template')
FINAL_LIBRARY=Construct((Segment('P5 primer',P5,'p5',inferred=True),role('reverse complement of captured genomic fragment'),Segment('P7 binding arm reverse complement',revcomp(P7_ARM),'r2',inferred=True),Segment('six-base P7 index reverse complement','N'*6,'cbc',placeholder=True,feature=SAMPLE_INDEX,inferred=True),Segment('P7 reverse complement',il.P7_RC,'p7',inferred=True)),name='INFERRED — GLOE single-index endpoint')
FINAL_CAPTION='INFERRED — single-index family printed in MethodsS1; six molecular index bases do not imply a six-cycle instrument mask. Read sequence is the reverse complement of the original captured DNA.'
SOURCE_OLIGO_EXCEPTIONS={'p7':'MethodsS1 p4 prints (X)6 between P7 and the complete annealing arm. The shared scraper does not expand this notation; local selftest verifies the complete printed polymer with (X)6 verbatim.'}
SEQ_PRIMERS=()
SEQUENCING_UNAVAILABLE='STAR Methods specifies single-end NextSeq500 runs of 75/84/92/150/160 bases, seven index cycles for single-index runs, or 8+8 for dual-index runs. Exact acquired run-primer oligos, run assignment and the dual-index molecular architecture are unassigned here; MethodsS1 prints a six-base P7 index.'
SEQUENCING_ENDING='GLOE-Pipe assigns each mapped read to the opposite strand: direct mode reports the original 3-prime terminal nucleotide, indirect mode shifts the inferred signal to the base immediately 3-prime of the break for lesion-conversion enzymes cutting 5-prime of the lesion.'
READ_SPECS=(sp.read_spec('Read 1','75/84/92/150/160; run dependent','reverse complement of original captured fragment, starting at its 3-prime boundary','include'),sp.read_spec('Index 1 (i7)','7 for source single-index runs','six molecular bases in represented P7 family; extra acquired cycle is not another declared index base','demultiplex'),sp.read_spec('Index 2 (i5)','not assigned to this molecular branch','source also reports dual-index runs; architecture deferred','undetermined'))

def oligos():
 rows=[]
 for n,s,m in [('3898 splinter',SPLINT,'3-prime NH2; phosphorothioate between final two random bases'),('3899 proximal',PROXIMAL,'5Phos; 3-prime biotin-dT'),('3790 second-strand primer',SECOND,''),('3791 distal',DISTAL1,''),('3792 distal',DISTAL2,''),('P5',P5,''),('P7 six-base family',P7,'six N positions denote selected sample index')]:rows.append(oligo(n,[Segment(n,s,placeholder=False)],mods=m))
 return rows

def proximal_geometry():
 # Six arbitrary terminal genomic bases are complementary to six random splint
 # bases. N is a degenerate substrate class, never a UMI annotation.
 sc=Scene()
 sc.strand('proximal',[Segment('adapter annealing core',PROXIMAL[:len(SPLINT_CORE)]),Segment('adapter tail',PROXIMAL[len(SPLINT_CORE):])],label='3899 adapter',mod5='5Phos',mod3='biotin-dT')
 sc.anneal('splinter',[Segment('splinter core',SPLINT_CORE),Segment('random splinter bases','N'*6,placeholder=True)],to='proximal',pair=('splinter core','adapter annealing core'),label='3898 splinter; final random bond phosphorothioate',mod3='NH2')
 sc.anneal('substrate',[role('genomic flank'),Segment('terminal genomic bases','N'*6,placeholder=True)],to='splinter',pair=('terminal genomic bases','random splinter bases'),label='native genomic 3-OH end')
 sc.ligation_end('substrate','terminal genomic bases',kind='splinted',polarity='3-prime',owner='substrate')
 sc.ligation_end('proximal','adapter annealing core',kind='splinted',polarity='5-prime',owner='adapter')
 return AdapterLigationGeometry.known('splinted',sc,statement='MethodsS1 p4, y23-28; FigureS1A: six-base random splinter, 3-NH2 block and 5-phosphorylated proximal adapter.')

def distal_geometry():
 sc=Scene()
 sc.strand('distal1',[Segment('distal core',DISTAL1[:-1]),Segment('distal T overhang','T')],label='3791 distal adapter')
 sc.anneal('distal2',[Segment('distal complementary arm',DISTAL2)],to='distal1',pair=('distal complementary arm','distal core'),label='3792 distal adapter')
 sc.ligation_end('distal1','distal T overhang',kind='cohesive',polarity='3-prime',owner='adapter')
 return AdapterLigationGeometry.known('cohesive',sc,statement='FigureS1A shows T-overhang adapter against A-tailed genomic duplex; MethodsS1 steps41-42. Source oligos3791/3792 define the duplex.')

def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Published GLOE2020 STAR Methods; MethodsS1 pp4,6-10; FigureS1A')
 dna=g.input('dna','gently extracted DNA carrying native free 3-prime hydroxyl ends',K.MOLECULAR,state=molecule('native substrate',(role('DNA containing native nicks, gaps or DSB ends'),)),primary=True)
 heat=g.input('heat','95 C denaturation, ice quench; yeast y23 or mammalian m22',K.CONDITION)
 denatured=g.product('denatured','denatured DNA exposing native 3-prime ends',K.MOLECULAR,state=molecule('single-stranded substrate',(role('captured genomic strand ending in native 3-OH'),)))
 g.react('denature','Denature intact DNA before fragmentation',S.THERMAL_TREATMENT,((dna,R.TEMPLATE),(heat,R.CONDITION)),(denatured,))
 adapter=g.input('proximal_adapter','3899 biotinylated/5-phosphorylated adapter annealed with blocked 3898 splinter',K.ADAPTER,state=proximal_geometry().scenes[0])
 ligase=g.input('t4_ligase','T4 DNA ligase; source ligation program',K.CATALYST)
 ligated=g.product('ligated','native end covalently sealed to proximal adapter',K.MOLECULAR,state=CAPTURED)
 g.react('capture_ligation','Seal native 3-OH to proximal adapter 5-phosphate before sonication',S.LIGATION,((denatured,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)),(ligated,),geometry=proximal_geometry())
 sonic=g.input('sonication','Bioruptor Pico sonication to about 200 nt',K.CONDITION)
 fragmented=g.product('fragmented','adapter-bearing sonicated fragments',K.MOLECULAR,state=CAPTURED)
 untagged=g.product('untagged','sonicated fragments lacking captured proximal end',K.MOLECULAR,state=molecule('other sonication fragments',(role('DNA without proximal biotin tag'),)))
 g.react('fragment','Sonicate only after original-end ligation',S.FRAGMENTATION,((ligated,R.TEMPLATE),(sonic,R.CONDITION)),(fragmented,untagged))
 beads=g.input('streptavidin','MyOne C1 streptavidin beads; SSC/NaOH washes, 95 C heat elution',K.SELECTION)
 retained=g.product('retained','eluted adapter-linked DNA after splinter removal',K.MOLECULAR,state=CAPTURED)
 wash_state=Scene()
 wash_state.strand('splinter',[Segment('3898 removed splinter',SPLINT)],label='removed splinter')
 wash_state.strand('other',[role('untagged sonication fragments and other unretained material')],label='other unretained molecules')
 wash=g.product('capture_wash','capture wash material including untagged fragments and removed splinter',K.WASTE,state=wash_state)
 g.react('enrich','AMPure cleanup; streptavidin capture, NaOH splinter removal and heat release',S.AFFINITY_CAPTURE,((fragmented,R.TEMPLATE),(untagged,R.TEMPLATE),(beads,R.SELECTION)),(retained,wash));g.discard(wash,'MethodsS1 steps30-37 separate supernatants and washes; released biotin chemical product not assigned.')
 primer=g.input('second_primer','3790 complementary-strand primer',K.OLIGO,state=Construct([Segment('3790',SECOND)]))
 phusion=g.input('phusion','Phusion FLASH polymerase mix',K.CATALYST)
 duplex=g.product('duplex','INFERRED — complementary strand synthesized from 3790',K.MOLECULAR,state=Construct((Segment('3790',SECOND),role('reverse complement of genomic fragment')),name='INFERRED — complementary strand'),inferred=True)
 g.react('second_strand','INFERRED — complementary strand synthesis',S.EXTENSION,((retained,R.TEMPLATE),(primer,R.STRUCTURAL),(phusion,R.CATALYST)),(duplex,))
 second_cleanup=g.input('second_cleanup_condition','AMPure purification and water elution; MethodsS1 steps39-40',K.CONDITION)
 clean_duplex=g.product('clean_duplex','purified second-strand duplex',K.MOLECULAR,state=duplex.state,inferred=True)
 second_waste=g.product('second_waste','second-strand cleanup supernatants and washes',K.WASTE)
 g.react('second_cleanup','Purify second-strand duplex with AMPure',S.CLEANUP,((duplex,R.TEMPLATE),(second_cleanup,R.CONDITION)),(clean_duplex,second_waste));g.discard(second_waste,'MethodsS1 steps39-40 remove bead supernatants and washes before water elution.')
 endprep=g.input('endprep','NEBNext UltraII EndPrep; 20 C 30 min / 65 C 30 min',K.CATALYST)
 polished=g.product('polished','end-polished and A-tailed duplex',K.MOLECULAR,state=molecule('end-prepped DNA',(role('3790-derived genomic duplex with distal 3-prime A overhang'),)))
 g.react('end_repair','Polish and A-tail with source EndPrep',S.END_REPAIR,((clean_duplex,R.TEMPLATE),(endprep,R.CATALYST)),(polished,))
 distal=g.input('distal_adapter','annealed distal adapter3791/3792; T-overhang',K.ADAPTER,state=distal_geometry().scenes[0])
 mix=g.input('ligation_mix','NEBNext UltraII Ligation Master Mix and enhancer; 20 C 20 min',K.CATALYST)
 preamp=g.product('preamp','INFERRED — distal-adapter-ligated template',K.MOLECULAR,state=TEMPLATE,inferred=True)
 g.react('distal_ligation','INFERRED — ligate distal adapter at A/T cohesive junction',S.LIGATION,((polished,R.TEMPLATE),(distal,R.STRUCTURAL),(mix,R.CATALYST)),(preamp,),geometry=distal_geometry())
 clean_condition=g.input('pre_pcr_cleanup_condition','salt/PEG cleanup and two AMPure purification rounds (steps43-46)',K.CONDITION)
 clean=g.product('clean_preamp','purified distal-ligated template',K.MOLECULAR,state=TEMPLATE,inferred=True)
 pre_waste=g.product('pre_pcr_waste','pre-PCR bead supernatants and washes',K.WASTE)
 g.react('pre_pcr_cleanup','Purify distal-ligated DNA before amplification',S.CLEANUP,((preamp,R.TEMPLATE),(clean_condition,R.CONDITION)),(clean,pre_waste));g.discard(pre_waste,'MethodsS1 steps43-46 bead cleanup removes supernatants and washes.')
 p5=g.input('p5','P5 primer',K.OLIGO,state=Construct([Segment('P5',P5)]))
 p7=g.input('p7','P7 six-base indexed primer family',K.OLIGO,state=Construct([Segment('P7 prefix',il.P7),Segment('six-base index','N'*6,placeholder=True,feature=SAMPLE_INDEX),Segment('P7 annealing arm',P7_ARM)]))
 q5=g.input('q5','Q5 DNA polymerase; eight cycles',K.CATALYST)
 amplified=g.product('amplified','INFERRED — amplified single-index GLOE library',K.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
 g.react('pcr','INFERRED — P5/P7 PCR; source eight cycles',S.PCR,((clean,R.TEMPLATE),(p5,R.STRUCTURAL),(p7,R.STRUCTURAL),(q5,R.CATALYST)),(amplified,),geometry=PCRGeometry.exact(TEMPLATE,PCRPrimer('P5',P5,anneal_nt=len(SECOND)),PCRPrimer('P7 family',P7,anneal_nt=len(DISTAL1)),statement='MethodsS1 p4, steps47-49; endpoint derived from printed primer extensions.'),relations=(MolecularRelation.PCR_PRODUCT,))
 ampure=g.input('ampure','two AMPure cleanups; Qubit QC; 4 nM sequencing pool',K.CONDITION)
 library=g.product('library','INFERRED — purified GLOE library',K.LIBRARY,state=FINAL_LIBRARY,inferred=True)
 waste=g.product('pcr_waste','post-PCR bead supernatants and washes',K.WASTE)
 g.react('cleanup','Purify library twice with AMPure; quantify and pool',S.CLEANUP,((amplified,R.TEMPLATE),(ampure,R.CONDITION)),(library,waste));g.discard(waste,'MethodsS1 steps48-50 remove supernatants and washes; no exact composition assigned.')
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE))
 g.declare_reviewed_scope(('dna','heat','denatured','denature','proximal_adapter','t4_ligase','ligated','capture_ligation','sonication','fragmented','untagged','fragment','streptavidin','retained','capture_wash','enrich','second_primer','phusion','duplex','second_strand','second_cleanup_condition','clean_duplex','second_waste','second_cleanup','endprep','polished','end_repair','distal_adapter','ligation_mix','preamp','distal_ligation','pre_pcr_cleanup_condition','clean_preamp','pre_pcr_waste','pre_pcr_cleanup','p5','p7','q5','amplified','pcr','ampure','library','pcr_waste','cleanup'),g.authority)
 return g.require_publishable()

def render_page():
 from batch_page import render
 return render(sys.modules[__name__])

def selftest():
 import re
 from checks import Check, Source, have
 from storage import sources_root
 c=Check();source=Source(sources_root()/HERE.parent.name/'j.molcel.2020.03.027_supp_mmc2_Methods_S1.pdf.txt','published MethodsS1 p4')
 if have(c,source,label='MethodsS1 oligos'):
  text=re.sub(r'\s+','',source.path.read_text())
  for name,seq in [('3898',SPLINT_CORE+'NNNNN*N-NH2'),('3899','PO4-'+PROXIMAL+'-Bio'),('3790',SECOND),('3791',DISTAL1),('3792',DISTAL2),('P5',P5),('P7',il.P7+'(X)6'+P7_ARM)]:c('source '+name,seq in text,True)
 protocol_graph();render_page();c.report()
