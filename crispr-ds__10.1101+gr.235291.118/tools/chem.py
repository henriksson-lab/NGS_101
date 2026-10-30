"""Protocol-owned, source-qualified incomplete preparation graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry,reported_adapter_end

def unknown(label):
 s=Scene();s.strand('preparation',[Segment('unassembled','['+label+']',placeholder=True)]);return s

def op(g,key,x,label,kind=S.MAGIC,reagent='Source reagents; complete formulation/geometry not assembled'):
 role=R.CATALYST if kind in {S.DIGESTION,S.EXTENSION,S.END_REPAIR,S.TAILING,S.NICK_SEALING} else R.STRUCTURAL if kind==S.HYBRIDIZATION else R.CONDITION
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

def oligo_scene(name,seq,mod5='',mod3=''):
 s=Scene();s.strand(name,[Segment(name,seq)],mod5=mod5,mod3=mod3);return s

def finish(g,x,limitation):
 lib=g.product('library','Source recovered sequencing preparation; complete molecular state unassembled',K.LIBRARY);g.react('library_identity','Source recovered sequencing library',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation=limitation));g.add_source_item(SourceItem('pending','Full source-specific substrate/adapter/end/read geometry and historical kit originals','Selected original scope only; ordinary unresolved work is not confirmed source absence.'));return g

from illumina import P5,TRUSEQ_READ1,TRUSEQ_READ2
MWS13=P5[:22]
MWS20='GTGACTGGAGTTCAGACGTGTGC'
MWS60=P5+TRUSEQ_READ1[4:]+'I'*12+'TGACT'
MSW61='GTCA'+'I'*12+revcomp(TRUSEQ_READ2)
SOURCE='Published GenomeResearch2018 PMC6169890 Methods140–172 andFig1C;SI TableS1/2 andMaterialscontents inspected. Commercial10bpUMI prototype differs old12bp adapter; whole geometry unassembled.'
def protocol_graph():
 g=ProtocolGraph('CRISPR-DS — source one-capture preparation',authority=SOURCE,representation='published_instance')
 x=g.input('dna','Source high-integrity extracted genomicDNA10–250ng selected amount unassigned; no whole-genome reference sequence invented',K.MOLECULAR,primary=True)
 cr=g.input('crrna','Source12targetedcrRNAs; wholeRNA scaffold/end chemistry and genomic cut map unassembled',K.OLIGO);tr=g.input('tracrrna','SourceIDTtracrRNA; exactsequence unreviewed',K.OLIGO);heat=g.input('guide_heat','Source95C5mincomplexing; individualorpooledcrRNA choice unassigned',K.CONDITION);guide=g.product('guides','SourcecrRNA:tracrRNA pool; unpaired RNA companions unassembled',K.MOLECULAR);g.react('guide_anneal','Source crRNA/tracrRNA assembly',S.HYBRIDIZATION,((cr,R.TEMPLATE),(tr,R.STRUCTURAL),(heat,R.CONDITION)),(guide,))
 cas=g.input('cas9','SourceNEBS.pyogenesCas9 about30nM;1xbuffer25C10min then DNA37Covernight',K.CATALYST);cut=g.product('cut','Cut/uncut DNA,guideRNA,Cas9 and reaction remnants; source selectedTP53fragments plus other genome not erased',K.MOLECULAR);g.react('cas9_digest','SourcemultiplexCas9 digestion30ul',S.DIGESTION,((x,R.TEMPLATE),(guide,R.STRUCTURAL),(cas,R.CATALYST)),(cut,));x=op(g,'cas9_heat',cut,'70C10minenzymeinactivation',S.THERMAL_TREATMENT)
 x=cleanup(g,'negative_selection',x,'Source0.5xAMPure3min;KEEP SUPERNATANT,bead-boundhighMWDNA separated')
 x=cleanup(g,'buffer_exchange',x,'Source1.8xAMPurepurification;KEEP recoveredDNA,50ulTELowelution')
 x=op(g,'erat',x,'NEBNextUltraIIERAT20C30min then65C30min;end-repair reagent retained althoughsource saysunnecessaryforCas9blunts',S.TAILING,'SourceNEBNextERAT reagents; completepolarity/junction geometry unassembled')
 ad=g.input('ds_adapter','CommercialTwinStrandprototype adapter;10bpcomplementaryrandomUMI andsource3dToverhang;exactwholeoligo/endinventoryunreviewed',K.ADAPTER);e=g.input('ligase','NEBNextligationmastermix,2.5ul15uMadapter20C15min',K.CATALYST);y=g.product('ligated','Ligated/unligatedDNA,prototypeadapters andreactionremnants',K.MOLECULAR);g.react('ligation','SourcecommercialDSadapterligation',S.LIGATION,((x,R.TEMPLATE),(ad,R.STRUCTURAL),(e,R.CATALYST)),(y,),geometry=reported_adapter_end('PrototypeDSadapter',kind='cohesive',polarity='3-prime',terminal_sequence='T',statement='SourceFig1C/Methods explicitly3dTadapter; exactwholeprototype andsourceDNA A-tailedpolaritynotfullyassembled.'))
 x=cleanup(g,'lig_cleanup',y,'0.8xAMPure,23ulwaterelution')
 x=pcr(g,'pcr1',x,'KAPAHiFiHotStartRealTime50ul,2uM MWS13/MWS20,stopatFluorescentStandard3;typically6–8cycles',(( 'SourceMWS13',oligo_scene('MWS13',MWS13)),('SourceMWS20',oligo_scene('MWS20',MWS20))))
 x=cleanup(g,'pcr1_cleanup',x,'0.8xAMPure,40ulwaterelution')
 blockers=[g.input('blocker_'+str(i),name+' actualsource12inosines;not12bpUMIclaim',K.OLIGO,chemical_state={'as_written':seq,'inosines':12,'sugar':'DNA','rendering_limit':'NativeinosinenotrepresentedbySceneexactDNAalphabet'}) for i,(name,seq) in enumerate((('MWS60',MWS60),('MSW61',MSW61)))]
 probe=g.input('capture_probes','Source21IDTTP53xGenLockdownprobes;7fragmentpools,exon7at40percent/exons8–9at90percent,final0.75pmol/ul',K.OLIGO);conditions=g.input('capture_conditions','SourcehistoricalIDTcaptureprotocol;75ulM270streptavidinbeads;fullworkupunreviewed',K.SELECTION);hyb=g.product('hybridized','Sourcehybridizationmixture andallunbound/probe/blocker species',K.MOLECULAR);g.react('hybrid_capture','Source one-round hybridizationcapture',S.HYBRIDIZATION,((x,R.TEMPLATE),(probe,R.STRUCTURAL),*((b,R.STRUCTURAL) for b in blockers),(conditions,R.SELECTION)),(hyb,));x=cleanup(g,'capture_workup',hyb,'Sourcehybridcapture retained/unboundphysicalseparation;fullwash/elutiondetailsunreviewed')
 x=pcr(g,'pcr2',x,'KAPAHiFiHotStart,0.8uM MWS13/indexedMWS21;20cycles98C30s/60C45s/72C45s',(('SourceMWS13',oligo_scene('MWS13',MWS13)),('SourceindexedMWS21;indexsequenceunassigned',None)))
 x=cleanup(g,'pcr2_cleanup',x,'Source0.8xAMPurepurification');x=assay(g,'sample_qubit',x,'SourceQubitHSassay aliquot')
 other=g.input('other_libraries','Independentsource sequencing-pool libraries',K.LIBRARY,primary=True);pooled=g.product('pool','Source diluted pooledlibraries',K.MOLECULAR);g.react('pooling','Source physicalsequencingpool',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pooled,));x=assay(g,'pool_tapestation',pooled,'SourceAgilent4200TapeStationpoolQC');x=assay(g,'pool_kapa',x,'SourceKAPAlibraryquantificationassay')
 phix=g.input('phix','Independentsourceabout1percentPhiXcontrolDNA',K.MOLECULAR,primary=True);load=g.product('loading_pool','SourcePhiXspikedloadingpool',K.MOLECULAR);g.react('phix_pool','ActualsourcePhiXspike',S.POOL,((x,R.TEMPLATE),(phix,R.TEMPLATE)),(load,))
 return finish(g,load,'MiSeqv3600cyclekit;actualrun/indexmaskunassembled. Prototype10bpUMI differsold12bp/conserved5bp;sampleindexfromMWS21unassigned. Duplexconsensus/clippingcomputational,notDNAcleanups. Mainwording5prime-dA-tailed versus3prime-dTadapterretained asunresolvedsourcepolaritywording;noguessedcorrection.')
def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph();c('Negativebindingselectionexplicitlykeepssupernatant','KEEP SUPERNATANT' in g.steps['negative_selection'].label,True);c('OnephysicalcapturebeforepostPCR',g.steps['pcr2'].inputs[0].material,'capture_workup_kept');c('ActualPhiXindependentinput',g.materials['phix'].primary,True);c('PCRgeometrynotclaimedcomplete',bool(g.publication_problems()),True)
 src=Source(sources_root()/'crispr-ds__10.1101+gr.235291.118'/'gr.235291.118_PMC6169890.xml.txt','PublishedoriginalMethodsPCR/capture')
 if have(c,src,label='Publishedmain'):
  t=''.join(src.path.read_text().split());c('WholeactualPCRprimers',all(s in t for s in (MWS13,MWS20)),True);c('WholeinosineblockerstranscribedwithoutUMInormalization',MWS60 in t and MSW61 in t,True)
 c.report()
if __name__=='__main__':selftest()
