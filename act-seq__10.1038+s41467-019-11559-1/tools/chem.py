"""Published bulk ACT-seq; corrected Supplementary Data2 owns all oligos."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import MolecularRelation, Disposition, InputRole as R, MaterialKind as K, ProtocolGraph, SequencingPlan, StepKind as S, Terminal
from protocol_model import molecule, role

PAGE_FILE='act-seq.html'
TITLE='ACT-seq — published bulk antibody-guided chromatin tagmentation'
NOTES='01_act-seq.html'
SOURCE='<a href="https://doi.org/10.1038/s41467-019-11559-1">Carter et al. (2019)</a>, published bulk Methods and corrected Supplementary Data2; <a href="https://doi.org/10.1038/s41467-020-18309-8">2020 Author Correction</a>.'
SUMMARY='Adapter-loaded Protein A–Tn5 binds the antibody before chromatin targeting. Magnesium activates tagmentation; DNA extraction, initial extension, indexed PCR and gel selection produce the bulk library.'
CAVEAT=('INFERRED — the A/B fragment-end assembly is derived from the printed adapters and PCR primers. Tn5 dimer composition, A/A and B/B frequencies and physical fates are unassigned. '
'This page represents bulk ACT-seq and the index pair printed in corrected Data2, not the separate iACT single-cell endpoint. '
'The published correction replaced duplicated Supplementary Data2 with the oligo workbook. The publisher workbook contains oligos; the local PMC mirror remains a stale mapping-statistics workbook. '
'The source’s single-cell Results describe 5000 cells per first well whereas Methods describe 10000; neither value is used to replace the bulk method.')
A=nx.S5+nx.ME
B=nx.S7+nx.ME
I5='TAGATCGC'
I7='CATGCCTA'
# Corrected workbook bulk library primers include the entire complex oligo.
P5=il.P5+I5+A
P7=il.P7+I7+B

def arm(name,seq,tag): return Segment(name,seq,tag,inferred=True)
def fragment(left=A,right=B):
 return Construct((arm('left adapter',left,'s5'),role('target-proximal genomic DNA'),arm('right adapter reverse complement',revcomp(right),'s7')),name='INFERRED — gap-filled fragment')
TEMPLATE=fragment()
FINAL_LIBRARY=Construct((arm('P5',il.P5,'p5'),Segment('printed i5',I5,'cbc',feature=feature('act_bulk_i5','sample_index','fixed'),inferred=True),arm('A',A,'s5'),role('target-proximal genomic DNA'),arm('B reverse complement',revcomp(B),'s7'),Segment('printed i7 reverse complement',revcomp(I7),'cbc',feature=feature('act_bulk_i7','sample_index','fixed'),inferred=True),arm('P7 reverse complement',il.P7_RC,'p7')),name='INFERRED — bulk A/B indexed library')
FINAL_CAPTION='INFERRED — productive A/B end assembly with the exact index pair printed in corrected Supplementary Data2; other index assignments and same-end product abundances are not implied.'
SEQ_PRIMERS=()
SEQUENCING_UNAVAILABLE='The paper specifies HiSeq2500 paired-end sequencing; loaded run-primer oligos and cycle lengths are unassigned in this capped model. Printed PCR primers establish molecular index sites but do not establish an acquired run mask.'
SEQUENCING_ENDING='HiSeq2500 paired-end sequencing (published Sequencing and data analysis).'
READ_SPECS=(sp.read_spec('Read 1','unassigned','one end of target-proximal genomic DNA','include'),sp.read_spec('Index 1 (i7)','unassigned','printed molecular i7 CATGCCTA; acquired cycles unassigned','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','printed molecular i5 TAGATCGC; acquired cycles unassigned','undetermined'),sp.read_spec('Read 2','unassigned','opposite end of target-proximal genomic DNA','include'))

def adapters():
 sc=Scene()
 for name,seq in [('A',A),('B',B)]:
  sc.strand(name,[Segment(name+' entry',seq[:-len(nx.ME)],'s5'),Segment(name+' ME',nx.ME,'me')],label='source '+name)
  sc.anneal(name+' complement',[Segment(name+' MErev',nx.ME_RC,'me')],to=name,pair=(name+' MErev',name+' ME'),label='pMENTS',mod5='5Phos')
 return sc

def oligos():
 return tuple(oligo(n,[Segment(n,s)],mods=m) for n,s,m in [('pMENTS',nx.ME_RC,'5Phos'),('Bulk 5-prime complex',A,''),('Bulk 3-prime complex',B,''),('Bulk 5-prime library',P5,''),('Bulk 3-prime library',P7,'')])

def protocol_graph():
 g=ProtocolGraph(TITLE,authority='Carter2019 published bulk Methods; corrected Data2; Author Correction2020')
 cells=g.input('cells','HEK293T cells; published bulk entry',K.BIOLOGICAL,state=molecule('cells',(role('chromatin DNA'),)),primary=True)
 perm=g.input('permeabilization','complex buffer; 10 min on ice',K.CONDITION)
 prepared=g.product('permeabilized','permeabilized cells',K.MOLECULAR,state=cells.state)
 g.react('permeabilize','Permeabilize and prepare 1000-cell aliquots',S.PERMEABILIZATION,((cells,R.TEMPLATE),(perm,R.CONDITION)),(prepared,))
 pool=g.input('adapters','corrected source adapter pool A/B, annealed with phosphorylated pMENTS',K.ADAPTER,state=adapters())
 enzyme=g.input('pa_tnp','Protein A–Tn5 fusion PA-Tnp; binding protein',K.MOLECULAR)
 catalytic=g.input('tn5_activity','Tn5 catalytic activity in targeted PA-Tnp',K.CATALYST)
 antibody=g.input('antibody','target-specific antibody',K.SELECTION)
 loaded=g.product('loaded','adapter-loaded PA-Tnp',K.MOLECULAR,state=adapters())
 g.react('load','Load complex adapters: 10 min at 25 C',S.BINDING,((pool,R.TEMPLATE),(enzyme,R.STRUCTURAL)),(loaded,))
 complexed=g.product('antibody_complex','preassembled antibody–PA-Tnp complex',K.MOLECULAR,state=loaded.state)
 g.react('complex_antibody','Bind antibody to loaded PA-Tnp: 60 min at 25 C',S.BINDING,((loaded,R.TEMPLATE),(antibody,R.SELECTION)),(complexed,))
 selection=g.input('target_binding','antibody-specific chromatin association and removal of unbound complexes',K.SELECTION)
 targeted_state=adapters()
 targeted_state.strand('chromatin',[role('chromatin DNA')],label='antibody-targeted chromatin; loaded adapters remain noncovalently associated')
 target=g.product('targeted','washed antibody-guided adapter-loaded PA-Tnp on chromatin',K.MOLECULAR,state=targeted_state)
 unbound=g.product('unbound','removed unbound complex',K.WASTE,state=molecule('removed complex',(role('unbound antibody and PA-Tnp adapter complexes'),)))
 g.react('target','Bind cells for 60 min at room temperature; three washes and rinse',S.AFFINITY_CAPTURE,((prepared,R.TEMPLATE),(complexed,R.STRUCTURAL),(selection,R.SELECTION)),(target,unbound));g.discard(unbound,'Published bulk washing removes unbound complex before Mg activation.')
 mg=g.input('magnesium','10 mM MgCl2; 30 min at 37 C',K.CONDITION)
 tagged=g.product('tagged','INFERRED — A/B-bearing gapped fragment class',K.MOLECULAR,state=molecule('gapped fragment',(role('A/B-bearing DNA with transposition gaps; strand geometry not assembled'),)),inferred=True)
 same_state=Scene()
 same_state.strand('aa',[role('A/A-bearing fragment')],label='possible A/A class')
 same_state.strand('bb',[role('B/B-bearing fragment')],label='possible B/B class')
 same=g.product('same_ends','INFERRED — same-end classes; fates unassigned',K.MOLECULAR,state=same_state,inferred=True)
 # The washed composite is both the genomic substrate and the carrier of loaded adapters.
 g.react('tagment','INFERRED — end classes after source-reported Mg-activated tagmentation',S.TAGMENTATION,((target,R.TEMPLATE),(target,R.STRUCTURAL),(catalytic,R.CATALYST),(mg,R.CONDITION)),(tagged,same))
 g.mark_terminal(Terminal(same.id,Disposition.UNRESOLVED,'End classes are logical possibilities of the declared pool. Their abundance and physical fate are unassigned; no absolute same-primer amplification prohibition is asserted.'))
 extract=g.input('extraction','EDTA/SDS/Proteinase K; phenol/chloroform and ethanol precipitation',K.CONDITION)
 recovered=g.product('recovered','purified tagmented DNA',K.MOLECULAR,state=tagged.state)
 waste=g.product('extraction_waste','extraction removed fractions',K.WASTE)
 g.react('extract','Stop tagmentation; digest protein and recover DNA',S.EXTRACTION,((tagged,R.TEMPLATE),(extract,R.CONDITION)),(recovered,waste));g.discard(waste,'Removed organic/aqueous fractions and precipitation supernatants.')
 polymerase=g.input('phusion','Phusion High-Fidelity PCR Master Mix',K.CATALYST)
 filled=g.product('filled','INFERRED — A/B fragment after initial extension',K.MOLECULAR,state=TEMPLATE,inferred=True)
 g.react('initial_extension','INFERRED — initial 72 C 5 min extension before PCR cycling',S.EXTENSION,((recovered,R.TEMPLATE),(polymerase,R.CATALYST)),(filled,))
 p5=g.input('p5_primer','printed bulk 5-prime library primer',K.OLIGO,state=Construct([Segment('P5 prefix',il.P5),Segment('printed i5',I5,'cbc',feature=feature('act_bulk_i5','sample_index','fixed')),Segment('A annealing arm',A)]))
 p7=g.input('p7_primer','printed bulk 3-prime library primer',K.OLIGO,state=Construct([Segment('P7 prefix',il.P7),Segment('printed i7',I7,'cbc',feature=feature('act_bulk_i7','sample_index','fixed')),Segment('B annealing arm',B)]))
 amplified=g.product('amplified','INFERRED — amplified bulk A/B library',K.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
 geometry=PCRGeometry.exact(TEMPLATE,PCRPrimer('bulk P5',P5,anneal_nt=len(A)),PCRPrimer('bulk P7',P7,anneal_nt=len(B)),statement='Corrected Data2 full primers; derived A/B template; bulk Methods PCR.')
 g.react('pcr','INFERRED — indexed A/B PCR; source 15 cycles 98 C/65 C/72 C',S.PCR,((filled,R.TEMPLATE),(p5,R.STRUCTURAL),(p7,R.STRUCTURAL),(polymerase,R.CATALYST)),(amplified,),geometry=geometry,relations=(MolecularRelation.PCR_PRODUCT,))
 gel=g.input('gel','agarose electrophoresis, band excision, QIAquick purification',K.CONDITION)
 library=g.product('library','INFERRED — gel-purified bulk A/B library',K.LIBRARY,state=FINAL_LIBRARY,inferred=True)
 outside=g.product('outside_band','material outside excised gel band and purification washes',K.WASTE)
 g.react('gel_select','Gel-excise desired fragment range; QIAquick purification',S.CLEANUP,((amplified,R.TEMPLATE),(gel,R.CONDITION)),(library,outside));g.discard(outside,'Source excises desired gel fragments; universal size range not assigned.')
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE))
 g.declare_reviewed_scope(('cells','permeabilization','permeabilized','permeabilize','adapters','pa_tnp','tn5_activity','antibody','loaded','load','antibody_complex','complex_antibody','target_binding','targeted','unbound','target','magnesium','tagged','same_ends','tagment','extraction','recovered','extraction_waste','extract','phusion','filled','initial_extension','p5_primer','p7_primer','amplified','pcr','gel','library','outside_band','gel_select'),g.authority)
 return g.require_publishable()

def render_page():
 from batch_page import render
 return render(sys.modules[__name__])

def selftest():
 import re
 from checks import Check, Source, have
 from storage import sources_root
 c=Check();source=Source(sources_root()/HERE.parent.name/'s41467-019-11559-1_41467_2019_11559_MOESM3_ESM.xlsx.txt','corrected publisher Data2; Author Correction2020')
 if have(c,source,label='corrected Data2 original'):
  text=re.sub(r'\s+','',source.path.read_text())
  for name,seq in [('pMENTS','5\'Phos-'+nx.ME_RC),('bulk A',A),('bulk B',B),('bulk P5',P5),('bulk P7',P7)]:c('source '+name,seq in text,True)
 protocol_graph();render_page();c.report()
