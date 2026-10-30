"""Source-local MAESTRO capture model, deliberately incomplete and private."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from illumina import TRUSEQ_READ1,TRUSEQ_READ2
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from pcr import PCRGeometry
SOURCE='PMC9089460 author manuscript Methods163–205; final published article not read. Probe Sequences File first438SNV row. No whole duplex-adapter recipe inherited.'
PROBE='CACATCCGGCGTTCGGCCGA'
CONVERSION_PRIMERS=(TRUSEQ_READ1[1:21],TRUSEQ_READ2[8:28])
def one(label):
 s=Scene();s.strand(label,[Segment(label,'['+label+']',placeholder=True)]);return s
def probe_scene():
 s=Scene();s.strand('first438SNVpanelprobe',[Segment('source probe body',PROBE)]);return s
def selected_probe_binding(probe,target,inferred=True):
 """Exact supplied-strand pairing only; never predicts capture or mismatch selectivity."""
 if not probe or set(probe+target)-set('ACGT') or len(probe)!=len(target):raise ValueError('Explicit equal-length DNA strands required')
 s=Scene();a=Segment('supplied probe',probe,inferred=inferred);b=Segment('supplied target',target,inferred=inferred);s.strand('probe',[a]);s.anneal('target',[b],to='probe',pair=('supplied target','supplied probe'));return s

def operation(g,key,x,label,kind=S.MAGIC,reagent='Source reagents; complete molecular geometry unassembled'):
 stock=g.input(key+'_stock',reagent,K.CONDITION);y=g.product(key+'_product',label+'; full mixture, offproducts and reagent remnants unassembled',K.MOLECULAR)
 kw={'reason':'Actual operation source-supported; full molecular inventory and geometry not assembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
 g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CONDITION)),(y,),**kw);return y

def pcr(g,key,x,cycles):
 p=g.input(key+'_primer1','Source postcapture forward amplification primer; exactidentity/sequence/landing unassigned',K.OLIGO);p2=g.input(key+'_primer2','Source postcapture reverse amplification primer; exactidentity/sequence/landing unassigned',K.OLIGO);e=g.input(key+'_enzyme','Source PCR catalyst and conditions; exact reagent formulation unreviewed',K.CATALYST);y=g.product(key+'_product','New amplified DNA plus original template, probes, primers and reaction remnants',K.MOLECULAR)
 g.react(key,f'Source postcapture PCR {cycles} cycles',S.PCR,((x,R.TEMPLATE),(p,R.STRUCTURAL),(p2,R.STRUCTURAL),(e,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(one('Full custom duplexUMI-adapter and PCR geometry unassembled'),reason='No imported duplex-library constructor or default historical adapter.'));return y

def capture(g,key,x,round_number):
 probe=g.input(key+'_probes','Source biotinylated short allele-selective probe panel; complete panel/end chemistry unassembled',K.OLIGO)
 r=g.input(key+'_reagents','Source xGen Hybridization and Wash kit/UniversalBlockers/Cot1; '+('1ug library,0.75pmol/ul probes' if round_number==1 else 'second capture: half volumes Cot1/blockers/probes'),K.SELECTION)
 y=g.product(key+'_hybridized','Probe-library hybrids plus mismatched/unbound library,free probes/blockers and all mixture',K.MOLECULAR)
 g.react(key+'_hybridize','95C30s;65→50C decrease1C/48min;50C≥4h,total16h',S.HYBRIDIZATION,((x,R.TEMPLATE),(probe,R.STRUCTURAL),(r,R.SELECTION)),(y,))
 kept=g.product(key+'_retained','Source capture-retained preparation; capture chemistry/geometry and recovery unassembled',K.MOLECULAR);waste=g.product(key+'_wash','Physical unbound/wash fractions; not computationally filtered reads',K.WASTE)
 g.react(key+'_wash_step','Source hybrid capture and heated washes50C,lid55C',S.AFFINITY_CAPTURE,((y,R.TEMPLATE),(r,R.SELECTION)),(kept,waste));g.discard(waste,'Physical source capture wash fractions; no perfect mutant/wildtype partition.');return kept

def protocol_graph():
 g=ProtocolGraph('MAESTRO — version-qualified sheared gDNA capture',authority=SOURCE,representation='published_instance')
 x=g.input('genomic','Source HapMap genomic DNA mixture; selected dilution value unassigned',K.MOLECULAR,primary=True)
 x=operation(g,'shear',x,'CovarisE220/LE220 shearing to source150bp fragments',S.FRAGMENTATION)
 retained=g.product('quantified','Retained shearedDNA after Picogreen assay allocation; actual assay fraction unassigned',K.MOLECULAR);qc=g.product('picogreen','Independent Picogreen dsDNA assay aliquot',K.QC);g.react('quantify','Source shearedDNA Quant-iTPicogreen assay; exact fraction unassigned',S.SPLIT,((x,R.TEMPLATE),),(retained,qc));g.mark_terminal(Terminal(qc.id,Disposition.QC,'Analytical DNA does not rejoin library preparation.'))
 used=g.product('library_input','Source20ng DNA library-construction allocation',K.MOLECULAR);rest=g.product('dna_remainder','Unused genomic DNA; mass and fate unassigned',K.MOLECULAR);g.react('allocate','Actual20ng library input allocation',S.SPLIT,((retained,R.TEMPLATE),),(used,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No silently consumed whole input.'))
 x=operation(g,'library',used,'KAPA HyperPrep with custom IDT dual-index duplexUMI adapters; full kit/adapter geometry unassembled')
 x=operation(g,'mass_boundary',x,'Source capture input1ug library; extra PCR or reindexing only when needed, selected choice unassigned')
 x=capture(g,'capture1',x,1);x=pcr(g,'pcr1',x,16);x=capture(g,'capture2',x,2);x=pcr(g,'pcr2',x,16)
 lib=g.product('final','Source final captured library; complete molecular geometry unassembled',K.LIBRARY);qc=g.product('final_qc','Source final quantification assay aliquot; fraction unassigned',K.QC);g.react('final_quantification','Final captured product quantified',S.SPLIT,((x,R.TEMPLATE),),(lib,qc));g.mark_terminal(Terminal(qc.id,Disposition.QC,'Quantification measurement separate from retained library.'))
 other=g.input('other_libraries','Independent source sequencing-pool libraries; identities/amounts unassigned',K.LIBRARY,primary=True);pool=g.product('pool','Source pooled final libraries',K.LIBRARY);g.react('pooling','Source final libraries pooled for sequencing',S.POOL,((lib,R.TEMPLATE),(other,R.TEMPLATE)),(pool,))
 g.sequence(pool,SequencingPlan(limitation='Source HiSeq2500 PE101 orHiSeqX PE151; selected run/library geometry unassigned. Duplex consensus/noise filtering computational, not molecular disappearance. Target depth10000x per site.'))
 g.add_source_item(SourceItem('pending','Final published main, historical KAPA/IDT customadapter and directly cited capture originals','Author manuscript source-qualified; no source absence or complete published-instance molecular validation.'))
 return g

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();s=selected_probe_binding('ACGT','ACGT');c('Exact supplied probe pairing builds',bool(s),True)
 try:selected_probe_binding('AAAA','AAAA');bad=False
 except ValueError:bad=True
 c('Mismatched supplied pairing rejected',bad,True)
 g=protocol_graph();c('Two captures preserve physical wash fractions',all(k in g.terminals for k in ('capture1_wash','capture2_wash')),True);c('PostcapturePCR consumes recovered capture only',g.steps['pcr1'].inputs[0].material,'capture1_retained');c('SecondPCR cycles remain MAESTRO16','16 cycles' in g.steps['pcr2'].label,True);c('No implicit source support for complete output',bool(g.publication_problems()),True)
 src=Source(sources_root()/'maestro__10.1038+s41551-022-00855-9'/'s41551-022-00855-9_PMC9089460.xml.txt','get_sources --slugs maestro__10.1038+s41551-022-00855-9; author manuscript Methods')
 if have(c,src,label='SourceAM'):
  t=src.path.read_text();c('Both whole conversion assay primers verified',all(p in t for p in CONVERSION_PRIMERS),True);c('Whole fluorophore/quencher probe retained in source','/5HEX/-CAACAAACATGCCATCTCCTTCTCCTGA-/ZEN//3IaBkFQ/' in t,True)
 src=Source(sources_root()/'maestro__10.1038+s41551-022-00855-9'/'s41551-022-00855-9_supp_NIHMS1775951-supplement-Probe_Sequences_File.xlsx.txt','Probe Sequences File first438SNV row')
 if have(c,src,label='SourceProbeWorkbook'):c('Whole source probe/coordinate row verified','1\t2329646\t2329666\t2329660\tA\tG\t'+PROBE+'\t+' in src.path.read_text(),True)
 c.report()
if __name__=='__main__':selftest()
