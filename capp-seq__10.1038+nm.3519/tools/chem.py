"""Source-local incomplete chemistry graph; no public finished marker."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(label):
 s=Scene();s.strand('material',[Segment('unassembled','['+label+']',placeholder=True)]);return s

def op(g,key,x,label,kind=S.MAGIC,reagent='Source reagents; exact chemical inventory/geometry unassembled'):
 role=R.CATALYST if kind==S.DIGESTION else R.STRUCTURAL if kind==S.HYBRIDIZATION else R.CONDITION
 st=g.input(key+'_stock',reagent,K.CATALYST if role==R.CATALYST else K.OLIGO if role==R.STRUCTURAL else K.CONDITION);y=g.product(key+'_product',label+'; retained native/copied DNA, companion species and reaction remnants unassembled',K.MOLECULAR)
 kw={'reason':'Source operation established, complete molecular inventory/geometry not yet assembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
 g.react(key,label,kind,((x,R.TEMPLATE),(st,role)),(y,),**kw);return y

def cleanup(g,key,x,label):
 st=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Source recovered preparation; carryover/recovery unquantified',K.MOLECULAR);w=g.product(key+'_waste','Actual separated cleanup fractions',K.WASTE);g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(st,R.CONDITION)),(y,w));g.discard(w,'Actual source physical cleanup; no flawless partition implied.');return y

def split(g,key,x,labels):
 ys=[g.product(key+'_'+str(i),lab,K.MOLECULAR) for i,lab in enumerate(labels)];g.react(key,key,S.SPLIT,((x,R.TEMPLATE),),ys);return ys

def qc(g,key,x,label):
 y,a=split(g,key,x,('Retained preparation after assay allocation; exact fraction unassigned',label));g.mark_terminal(Terminal(a.id,Disposition.QC,'Separate source diagnostic material; not returned to preparation.'));return y

def pcr(g,key,x,label,pair_labels):
 ps=[g.input(key+'_p'+str(i),p,K.OLIGO) for i,p in enumerate(pair_labels)];e=g.input(key+'_catalyst','Source PCR enzyme/conditions; full molecular geometry unassembled',K.CATALYST);y=g.product(key+'_product','New PCR DNA plus original DNA/free primers/reagents; no silent purification',K.MOLECULAR);g.react(key,label,S.PCR,((x,R.TEMPLATE),*((p,R.STRUCTURAL) for p in ps),(e,R.CATALYST)),(y,),geometry=PCRGeometry.unavailable(scene('Whole source-specific template and inward primer sites unassembled'),reason='Exact individual oligo inventories do not establish complete operational PCR geometry.'));return y

def library(g,x,label,limitation):
 y=g.product('library',label,K.LIBRARY);g.react('library_identity','Source recovered sequencing preparation',S.GENERIC,((x,R.TEMPLATE),),(y,),changes_molecule=False);g.sequence(y,SequencingPlan(limitation=limitation));g.add_source_item(SourceItem('pending','Historical kit/vendor/predecessor originals and full molecular/read geometry','Ordinary bounded deferred scope; no accepted absence or whole source-validation claim.'));return g

SOURCE='Version-qualified PMC4016134 author manuscript Methods128–146 and associated SI A3/A4 pp12–13; final published main/SI not yet read.'
def graph_for_replicates(replicates):
 if type(replicates)is not int or not 4<=replicates<=6:raise ValueError('Source4–6parallelpostcapturePCRreactions only')
 g=ProtocolGraph('CAPP-Seq — original source-qualified low-input plasmaDNA preparation',authority=SOURCE,representation='published_instance')
 x=g.input('plasma_dna','Source7–32ng extracted plasmaDNA; no further fragmentation; sample assignment uncurated',K.MOLECULAR,primary=True)
 x=op(g,'kapa_with_bead',x,'SourcepreligationKAPAlibrary enzymatic/with-bead boundary;full intermediate order/chemistryunassembled')
 ad=g.input('truseq_adapters','Source indexedIlluminaTruSeqadapters,100foldmolarexcess; exact historical ends/index unassigned',K.ADAPTER);e=g.input('ligase','Source ligase reaction formulation unreviewed',K.CATALYST);y=g.product('ligated','Ligated/unligatedDNA,adapters and reaction remnants; exact junction unassembled',K.MOLECULAR);g.react('ligation','Source16h16Cligation',S.LIGATION,((x,R.TEMPLATE),(ad,R.STRUCTURAL),(e,R.CATALYST)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(scene('Exact source DNA/adaptor ends unassembled'),),'Known source ligation; actual molecular geometry not yet assembled.'))
 x=cleanup(g,'peg_selection',y,'Source single-step40ulPEGbuffer0.8Xsize selection')
 x=pcr(g,'pre_capture_pcr',x,'Source4–9cycles dependent on inputDNA,500nMbackboneprimers; exact choice unassigned',('SourceIllumina forward backbone primer','SourceIllumina reverse backbone primer'))
 x=qc(g,'library_qc',x,'Separate sourceNanoDrop purity,KAPAqPCR concentration andBioanalyzerDNA1000 length assays; fractions unassigned')
 others=g.input('other_indexed','Independent source indexedlibraries; selectedpooltotal9–12unassigned',K.LIBRARY,primary=True);pool=g.product('capture_pool','Source9–12indexedlibrarypool; exactmembers/amountsunassigned',K.MOLECULAR);g.react('pool_capture','Source9–12indexedlibraries co-captured',S.POOL,((x,R.TEMPLATE),(others,R.TEMPLATE)),(pool,))
 x=op(g,'hybrid_selection',pool,'Source custom biotinylatedDNASeqCapEZChoice selector hybridselection; historical kit protocol unreviewed',S.HYBRIDIZATION,'SourceRocheNimbleGenSeqCapEZChoice,hg19panel; exact probe/end/capturegeometryunassembled')
 x=cleanup(g,'capture_recovery',x,'Source hybrid-selection physical retained/unbound separation; exactworkupunassembled')
 branches=split(g,'post_capture_split',x,tuple('Source parallel50ulPCR'+str(i+1)+' input; exactfractionunassigned' for i in range(replicates)))
 ys=[pcr(g,'post_capture_'+str(i),t,'Source12–14cycles,1xKAPAHiFiHotStartReadyMix,2uMbackboneprimers,50ul;exactcyclechoiceunassigned',('SourceIllumina forward backbone primer','SourceIllumina reverse backbone primer')) for i,t in enumerate(branches)]
 pooled=g.product('post_capture_pool','Pooled source parallelPCRproducts',K.MOLECULAR);g.react('pool_post_capture','Source4–6PCRreactions physically pooled',S.POOL,tuple((y,R.TEMPLATE) for y in ys),(pooled,));x=cleanup(g,'qiaquick',pooled,'SourceQIAquickPCRpurification')
 return library(g,x,'OriginalCAPPSeq sequencing library; fullkit/adapter/readgeometryunassembled','IlluminaHiSeq2000 PE100. OriginalindexedTruSeqsamples; no assignedmoleculeUMI/lateriDESrecipe. Barcode/variantfiltering computational,not chemical removal. PlasmaDNAunfragmented;Covaris200bp appliesothergenomicDNA route only.')

def protocol_graph():return graph_for_replicates(4)

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=graph_for_replicates(4);h=graph_for_replicates(6);c('Source4and6parallelPCRsretained',(len(g.steps['pool_post_capture'].inputs),len(h.steps['pool_post_capture'].inputs)),(4,6));c('PostPCRQIAquickconsumesphysicalpool',g.steps['qiaquick'].inputs[0].material,'post_capture_pool');c('Noadditionalplasmafragmentation','shear' not in g.steps,True);c('LigationselectionbeforePCR',g.steps['pre_capture_pcr'].inputs[0].material,'peg_selection_kept');c('Noacceptedfullgeometry',bool(g.publication_problems()),True)
 try:graph_for_replicates(3);bad=False
 except ValueError:bad=True
 c('Unpublishedparallelcountrejected',bad,True)
 src=Source(sources_root()/'capp-seq__10.1038+nm.3519'/'nm.3519_PMC4016134.xml.txt','SourceAMMethodsNGSlibraryconstruction')
 if have(c,src,label='SourceAM'):
  t=src.path.read_text();c('Source16h16Cplus100foldadaptersverified','16 h at 16 °C using 100-fold molar excess' in t,True)
 c.report()
if __name__=='__main__':selftest()
