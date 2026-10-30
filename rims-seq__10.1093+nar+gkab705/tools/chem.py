"""Protocol-owned, source-qualified incomplete preparation graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from pcr import PCRGeometry

def unknown(label):
 s=Scene();s.strand('preparation',[Segment('unassembled','['+label+']',placeholder=True)]);return s

def op(g,key,x,label,kind=S.MAGIC,reagent='Source reagents; complete formulation/geometry not assembled'):
 role=R.CATALYST if kind in {S.DIGESTION,S.DEPROTEINIZATION,S.EXTENSION,S.END_REPAIR,S.TAILING,S.NICK_SEALING} else R.STRUCTURAL if kind==S.HYBRIDIZATION else R.CONDITION
 st=g.input(key+'_stock',reagent,K.CATALYST if role==R.CATALYST else K.OLIGO if role==R.STRUCTURAL else K.CONDITION);y=g.product(key+'_product',label+'; all substrates/products/companion species and reaction remnants retained until stated workup',K.MOLECULAR)
 kw={'reason':'Known source operation; exact molecular states, ends and geometry not yet assembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
 g.react(key,label,kind,((x,R.TEMPLATE),(st,role)),(y,),**kw);return y

def cleanup(g,key,x,label):
 st=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Source recovered preparation; recovery/carryover not presumed perfect',K.MOLECULAR);w=g.product(key+'_waste','Actual separated physical workup fractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(st,R.CONDITION)),(y,w));g.discard(w,'Source physical separation, not computational filtering.');return y

def split(g,key,x,labels):
 ys=[g.product(key+'_'+str(i),lab,K.MOLECULAR) for i,lab in enumerate(labels)];g.react(key,key,S.SPLIT,((x,R.TEMPLATE),),ys);return ys

def assay(g,key,x,label):
 y,a=split(g,key,x,('Retained preparation after source assay allocation; assay amount unassigned',label));g.mark_terminal(Terminal(a.id,Disposition.QC,'Separate analytical material; no implicit return to preparation.'));return y

def pcr(g,key,x,label,primers):
 ps=[g.input(key+'_p'+str(i),lab,K.OLIGO,state=state) for i,(lab,state) in enumerate(primers)];e=g.input(key+'_enzyme',label+' source catalyst',K.CATALYST);y=g.product(key+'_product','New DNA copies plus original DNA/free primers/reagents',K.MOLECULAR);geom=PCRGeometry.unavailable(unknown('Whole operational PCR geometry not assembled'),reason='Source primer identities do not by themselves establish all inward sites/products.',symmetric_species=len(ps)==1);g.react(key,label,S.PCR,((x,R.TEMPLATE),*((p,R.STRUCTURAL) for p in ps),(e,R.CATALYST)),(y,),geometry=geom);return y

def finish(g,x,limitation):
 lib=g.product('library','Source recovered sequencing preparation; complete molecular state unassembled',K.LIBRARY);g.react('library_identity','Source recovered sequencing library',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation=limitation));g.add_source_item(SourceItem('pending','Full source-specific substrate/adapter/end/read geometry and historical kit originals','Selected original scope only; ordinary unresolved work is not confirmed source absence.'));return g

SOURCE='Published RIMS-seq main PMC8565308 Methods114–120/Results188–212 and publishedSI Supplementarytext1/FigS1caption. Selected historical unmethylated-adapter RIMS route; bisulfite/methylated-adapter/Q5U comparison is separate. Exact kit revision and molecular/read geometry unassembled.'
def protocol_graph():
 g=ProtocolGraph('RIMS-seq — source limited alkaline treatment',authority=SOURCE,representation='published_instance')
 x=g.input('genomic_dna','Source100ng isolated genomic DNA; actual sample/coordinates/methylated sites unassigned',K.MOLECULAR,primary=True)
 x=op(g,'sonication',x,'Source CovarisS2 in50ul1xTE,nominal200bp insert setting; actual distribution not predicted',S.FRAGMENTATION,'Source mechanical fragmentation conditions')
 x=op(g,'historical_preligation',x,'Source NEBNextUltraII E7645 historical manufacturer preparative steps; exact state/formulation/version before ligation unassembled')
 adapter=g.input('loop_adapter','Source regular UNMETHYLATED loop-shaped adapter; exact sequence/end geometry and USER substrate positions unassembled',K.ADAPTER)
 enzyme=g.input('ligation_reagents','Historical kit ligation enzyme/reagents, exact identities and conditions not assembled from selected main',K.CATALYST)
 y=g.product('adapter_ligation_product','Source ligated/unligated DNA, free unmethylated-loop adapter and all reaction companions',K.MOLECULAR)
 from adapter_ligation import AdapterLigationGeometry
 g.react('adapter_ligation','Source standard unmethylated-loop adapter ligation',S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(enzyme,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown('Historical source geometry unassembled; source disclosure not audited'),),'Historical kit version, literal oligos/ends and exact sealed junctions not assembled; no source absence assertion or current kit substitution.'))
 x=op(g,'user',y,'Source included USER treatment BEFORE alkaline exposure; exact adapter cleavage/end inventory unassembled',S.DIGESTION,'Source kit USER reagent; not a second post-deamination removal of newly formed genomic U')
 x=op(g,'limited_deamination',x,'Selected source add1MNaOHpH13stock to0.1Mfinal,60C3h; limitedm5C→T/C→U and other original/damaged classes retained, no perfect conversion prediction',S.CHEMICAL_CONVERSION,'Source NaOH and heat; no invented enzyme-mediated conversion')
 x=op(g,'neutralization',x,'Source immediately cool onice and add1Maceticacidstock to0.1Mfinal; other0.5/1MNaOH conditions are separate source tests',S.QUENCHING,'Source acetic acid and cooling; actual final pH not assumed from nominal equal molarity')
 x=cleanup(g,'zymo',x,'Source ZymoOligoCleanandConcentratorD4060;20ul0.1xTE eluate. No skip-cleanup shortcut imported from an untested discussion suggestion')
 x=pcr(g,'library_pcr',x,'Source historicalUltraII ER7645 protocol with E7337A multiplexoligos; sample-specificPCRcycles optimized, actual enzyme/primer identities/cycles unassembled here',(('Source E7337A PCR primer1; full strand/index unassembled',None),('Source E7337A PCR primer2; full strand/index unassembled',None)))
 x=cleanup(g,'post_pcr',x,'Source0.9xNEBNextSamplePurificationBeadsE7137AA,elute25ul0.1xTE; no perfect removal of every original uracil-bearing molecule')
 x=assay(g,'tapestation',x,'Separate source TapeStationHighSensitivityDNA1000 analytical allocation; amount unassigned')
 return finish(g,x,'Source Illumina paired-end sequencing, actual instrument/runlength/indexcycles and complete molecular primer geometry unassembled. SourceSI attributes R1C→T versus R2G→A imbalance to limitedm5C deamination and U-intolerant amplification; this readout explanation does not prove quantitative perfect U-gating, conversion, strand origin or physical destruction of all U-containing originals. Bisulfite control uses methylatedadapter/Q5U and separate lambdaDNAspike; not inherited here. OptionalXP12validation control and8oxoG-LCMS oligo test are separate experiments, not automatic sample components. Variant/motif filtering and assembly are computational.')

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph();order=list(g.steps)
 c('USER precedes alkaline genomic deamination',order.index('user')<order.index('limited_deamination'),True)
 c('Acid neutralization precedes physical column',order.index('neutralization')<order.index('zymo'),True)
 c('Genomic deamination is chemical conversion',g.steps['limited_deamination'].kind,S.CHEMICAL_CONVERSION)
 c('QC is a separate analytical fraction','tapestation_1'in g.terminals,True)
 c('Historical exact geometry remains incomplete',bool(g.publication_problems()),True)
 src=Source(sources_root()/'rims-seq__10.1093+nar+gkab705'/'nar_gkab705_PMC8565308.xml.txt','PublishedmainMethods114–120 andseparatebisulfite122–130')
 if have(c,src,label='Publishedmain'):
  t=''.join(src.path.read_text().split());rims=t.split('RIMS-seqlibrarypreparation',1)[1].split('Bisulfite-seqlibrarypreparation',1)[0]
  c('Selected source adapter is unmethylated, not methylated comparator','regularunmethylatedloop-shapedadapter'in rims,True)
  c('SourceRIMS PCR not silently borrowed Q5U from comparator','Q5U'not in rims and 'E7337A'in rims,True)
  c('Source includes USER then added alkaline treatment',rims.index('AftertheUSERtreatment')<rims.index('heatalkalinedeamination'),True)
 c.report()
if __name__=='__main__':selftest()
