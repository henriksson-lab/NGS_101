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

from illumina import P5,P7
CS1='ACACTGACGACATGGTTCTACA'
CS2='TACGGTAGCAGAGACTTGGTCT'
PE1=P5[:25]
PE2=P7[:24]
INDEX='GTATCGTCGT'
SOURCE='Final2012paper+SI in original PTABexhibit;SIpp2–3/TableS2/FigS1 controls selected case, published discrepancies retained in note.'
def barcode_primer(orientation='forward',indexed=False):
 if orientation not in {'forward','reverse'}:raise ValueError('Two source read orientations only')
 cs=(CS2 if orientation=='forward' else CS1) if indexed else (CS1 if orientation=='forward' else CS2)
 s=Scene();segs=[Segment('PE2' if indexed else 'PE1',PE2 if indexed else PE1)]
 if indexed:segs.append(Segment('FLD0001 sample index',INDEX,feature=feature('tam_fld0001','sample_index','fixed')))
 segs.append(Segment('source CS handle',cs));s.strand('source TableS2 primer',[*segs]);return s

def protocol_graph():
 g=ProtocolGraph('TAm-Seq — source SI bidirectional sample-indexed libraries',authority=SOURCE,representation='published_instance')
 x=g.input('input_dna','Source extracted fragmented plasmaDNA; amount/sample identity unassigned',K.MOLECULAR,primary=True)
 x=pcr(g,'preamplification',x,'Source10ul,50nM each targetprimer,FastStartHiFi;95C10min then15cycles95C15s/60C4min',('Source CS1-tagged forward target primer pool','Source CS2-tagged reverse target primer pool'))
 x=op(g,'exo_sap',x,'4ulExoSAPit,37C15min then80C15min',S.DIGESTION,'Source ExoSAPit; primer/dNTP hydrolysis products retained')
 x=op(g,'dilute5',x,'Source fivefold dilution in PCR-grade water',S.GENERIC,'PCR-gradewater')
 used,rest=split(g,'access_input',x,('Source1ul diluted preparation plus4ul pre-sample mastermix','Unused preamplified preparation; exact amount/fate unassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Unused source stock not silently consumed.'))
 st=g.input('ifc_conditions','Source IFCAXcontroller/individual targetprimerpairs; exact chamberoccupancy/unusedvolume unassigned',K.CONDITION)
 x=g.product('access_array_product','Source target-specific AccessArray33nlchambers;source0.7percentinputvolume each;full48chambergeometryunassembled',K.MOLECULAR);rem=g.product('access_array_remainder','Unallocated IFCDNA/reagents; exact amount/fate unassigned',K.MOLECULAR);g.react('access_array','SourceAccessArrayIFC spatialpartitioning',S.PARTITION,((used,R.TEMPLATE),(st,R.CONDITION)),(x,rem));g.mark_terminal(Terminal(rem.id,Disposition.UNRESOLVED,'WholeloadedIFCvolume notsilentlyassignedtoonechamber.'))
 x=pcr(g,'singleplex',x,'Source35cycleBiomark programme;FastStartHiFi/CS-taggedtargetprimers',('Source single target CS1-forward primer per chamber','Source single target CS2-reverse primer per chamber'))
 x=op(g,'harvest',x,'Source0.05percentTween20 harvest withIFCAXcontroller; recovered source sample product',S.EXTRACTION,'Source harvesting solution')
 used,rest=split(g,'harvest_input',x,('Source1ul harvested product only; water added in following dilution','Remaining harvested stock; amount/fate unassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No invented full-stock consumption.'))
 x=op(g,'dilute100',used,'Source100fold dilution',S.GENERIC,'99ulPCR-freewater')
 forward,reverse,rest=split(g,'two_plate_inputs',x,('Source1ul to first orientation20ulPCR','Source1ul to second orientation20ulPCR','Remaining100fold diluted stock; amount/fate unassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Two actual1ul transfers do not consume entire dilution.'))
 ys=[]
 for name,t in [('forward',forward),('reverse',reverse)]:
  ys.append(pcr(g,'barcode_'+name,t,'SourceSI15barcodePCRcycles95C15s/60C30s/72C4min;main10cycles remains disputed',('SourceTableS2 PE1-'+('CS1' if name=='forward' else 'CS2'),'SourceTableS2 PE2-FLD0001-'+('CS2' if name=='forward' else 'CS1'))))
  ys[-1]=qc(g,'bioanalyzer_'+name,ys[-1],'Source prepool Agilent2100 Bioanalyzer assay aliquot; expectedinsert~180bp')
 pooled=g.product('orientation_pool','Both source read orientations physically pooled',K.MOLECULAR);g.react('pool_orientations','Source all barcoded products pooled',S.POOL,tuple((y,R.TEMPLATE) for y in ys),(pooled,))
 x=cleanup(g,'ampure',pooled,'SourceAMPureXP bead:amplicon1.8:1cleanup');x=qc(g,'pool_quantification',x,'Source pooled library Bioanalyzer quantification assay')
 return library(g,x,'Source two-orientation sample-indexed library pool; whole target-specific geometry unassembled','GAIIx SE100plus10base index read,customCS1/CS2 readprimers and reversecomplements forindex. Fixed sampleindex introduced aftertargetPCR,not moleculeUMI. Main10vsSI15barcodecycles andSIp3PE1vsTableS2/FigS1PE2placement unresolved;selectedSI15/TableS2 topology.')

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph();c('Actualtwo1ultransfersretainstock','two_plate_inputs_2' in g.terminals,True);c('ExoSAPbeforefivefolddilution',g.steps['dilute5'].inputs[0].material,'exo_sap_product');c('OrientationPCRsremainindependent',g.steps['barcode_forward'].inputs[0].material!=g.steps['barcode_reverse'].inputs[0].material,True);c('Analyticalprepoolassaysseparate',all('bioanalyzer_'+n+'_1' in g.terminals for n in ('forward','reverse')),True);c('Noacceptedfullgeometry',bool(g.publication_problems()),True)
 for n in ('forward','reverse'):
  c('Actualsourceprimerbuilds '+n,bool(barcode_primer(n,True)),True)
 src=Source(sources_root()/'tam-seq__10.1126+scitranslmed.3003726'/'source_pdf_1ef17b49d5c4fba54e01.pdf.txt','ExactPTABoriginal2012paper+SI;TableS1/S2originalp26rendered')
 if have(c,src,label='PublishedSI'):
  t=src.path.read_text();c('WholeliteralCSandPEbodiesverified',all(v in t for v in (CS1,CS2,PE1,PE2,INDEX)),True);c('PublishedSIidentityverified','10.1126/scitranslmed.3003726' in t and 'Published 30 May 2012' in t,True)
 c.report()
if __name__=='__main__':selftest()
