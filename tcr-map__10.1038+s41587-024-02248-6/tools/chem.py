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

SOURCE='PublishedPMC11994455 main74/104/110/252–280;MOESM1reportingsummaryandMOESM2UBMTAheadersadministrativeonly.3originalcapreached;ngsprep2019TScan/2023TScanIIreferencesnotwholerecipeinheritance.'
def graph_for_selection(method):
 if method not in {'facs','magnetic'}:raise ValueError('ActualsourceFACSormagneticbranchrequired')
 g=ProtocolGraph('TCR-MAP — '+method+'sourcebranch',authority=SOURCE,representation='published_instance')
 targets=g.input('target_cells','SourcepremadeG5-targetAPClibrary;antigen/MHC/circuitstockconstructionunassembled',K.BIOLOGICAL,primary=True);jurkats=g.input('srt_jurkats','IndependentpremadeSrtA-Jurkat/TCRstock;upstreamengineeringunassembled',K.BIOLOGICAL,primary=True);mix=g.product('cell_mixture','SourceAPC/Jurkatcellmixture;nativeunresponsive/reactivecellsallretained',K.MOLECULAR);g.react('coculture','SourceTCR-MAPcellcoculture;actualstockidentitiesunassigned',S.POOL,((targets,R.TEMPLATE),(jurkats,R.TEMPLATE)),(mix,))
 substrate=g.input('biotin_peptide','SourceLPETG-biotinsubstrate;biotinisnotDNAindex',K.MOLECULAR);labeled=g.product('labeled_cells','Sourcecovalentcell-surfacebiotintransfermixtureplusunlabeledcells/freepeptide/reagents;noefficiencyorspecificityprediction',K.MOLECULAR);g.react('sortase_transfer','Sourcecell-associatedSrtA-mediatedLPETGbiotintransfer',S.MAGIC,((mix,R.TEMPLATE),(substrate,R.STRUCTURAL)),(labeled,),reason='Sourcefunctionaltransferknown;completecell-circuit/protein-junctionchemistryunassembled. SrtAispartofinputcellstock,notinventednewfreeenzymereagent.')
 if method=='facs':
  x=cleanup(g,'facs',labeled,'SourceantiCD40/BV421streptavidinstain30min4C,PBEwash,SonyMA900sorttop1percentbiotinylatedtargetcells;carryoverunquantified')
 else:
  x=cleanup(g,'wash_before_magnet',labeled,'Sourcemagneticbranch:enzyme-freedissociation/PBS2mMEDTA,threePBEwashes')
  x=cleanup(g,'magnet1',x,'SourceMiltenyistreptavidinmicrobeads/magneticcolumn;boundtargetcellelution')
  sampled,retained=split(g,'round1_save',x,('Sourcefirst-roundcellaliquotsavedforgDNA/NGS;exactfractionunassigned','Sourceboundcelleluateforfurthergrowth/re-enrichment'))
  g.mark_terminal(Terminal(sampled.id,Disposition.UNRESOLVED,'Actualsavedfirst-roundgDNA/NGSbranch;itsseparatepreparationnotassembledhere.'))
  x=op(g,'regrow',retained,'Sourceeluateplatedandgrown1week;subsequentcoculture/enrichmentcontinues;division/carryovergeometryunassembled')
  others=g.input('round2_jurkats','IndependentsourceSrtA-Jurkat/TCRre-enrichmentstock',K.BIOLOGICAL,primary=True);r2=g.product('round2_mix','SourcegrownAPC/secondJurkatstockmixture',K.MOLECULAR);g.react('round2_coculture','Sourcefurthercoculturebeforesecondcapture',S.POOL,((x,R.TEMPLATE),(others,R.TEMPLATE)),(r2,));st2=g.input('round2_biotin','SourcefreshLPETG-biotinre-enrichmentstock;actualamountunassigned',K.MOLECULAR);z=g.product('round2_labeled','Sourcefurtherbiotin-transfermixture;fullgeometryunassembled',K.MOLECULAR);g.react('round2_transfer','Sourcefurthercoculture/biotinylation',S.MAGIC,((r2,R.TEMPLATE),(st2,R.STRUCTURAL)),(z,),reason='Knownsourcefurthercocultureandtwo-consecutive-enrichmentbranch;exactmolecularinventoryunassembled.');x=cleanup(g,'magnet2',z,'Sourcesecondconsecutivestreptavidinmagneticcolumnisolation;fullinterveningwashvolumesunassembled')
 kit=g.input('genejet_stock','SourceGeneJETGenomicDNAPurificationkit;historicalinternalchemistryunreviewed',K.CONDITION);dna=g.product('gdna','SourcepurifiedgDNAofisolatedcells;retainedcarrier/carryoverandmassunassigned',K.MOLECULAR);w=g.product('genejet_waste','Actualkitremovedphysicalfractions;compositionsunassembled',K.WASTE);g.react('genejet','SourcegDNAextraction',S.MAGIC,((x,R.TEMPLATE),(kit,R.CONDITION)),(dna,w),reason='SourcegDNApurificationestablished;exacthistoricalkitlysis/washpartitionsunassembled.');g.discard(w,'ActualGeneJETphysicalseparation;noidealpurityprediction.')
 x=op(g,'cited_ngs_preparation',dna,'SourceantigenlibraryIlluminapreparationfollowingTScan2019ref13andTScanII2023ref15;actualselectedPCRs/oligos/cleanupnotassembled')
 return finish(g,x,'SourceMiSeqorNextSeq,STANDARDIlluminasequencingprimers;noTScancustomprimerimport. Sourceinputlibrary≥40xrepresentationisindependentcontrolpreparation,notpooledmolecularnormalizer. Barcodedpeptideidentity/readenrichmentnotmoleculeUMI. Exactindexvalues/cycles/librarymoleculargeometryunassembled.')
def protocol_graph():return graph_for_selection('facs')
def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();f=protocol_graph();m=graph_for_selection('magnetic');c('FACSdefaultdoesnotinventgrowth','regrow' not in f.steps,True);c('Magneticfirst-roundsaveisphysical','round1_save_0' in m.terminals,True);c('CitedNGSboundaryconsumespurifiedDNA',f.steps['cited_ngs_preparation'].inputs[0].material,'gdna');c('NoexactsampleindexorPCRrecipeinherited',bool(f.publication_problems()),True)
 try:graph_for_selection('singlecell');bad=False
 except ValueError:bad=True
 c('Unspecifiedselectionbranchrejected',bad,True)
 src=Source(sources_root()/'tcr-map__10.1038+s41587-024-02248-6'/'s41587-024-02248-6_PMC11994455.xml.txt','PublishedmainNGSMethods280;transferandselectionboundaries')
 if have(c,src,label='PublishedMain'):
  t=''.join(src.path.read_text().split());c('Knownstandardsequencingprimerboundaryverified','usingthestandardIlluminasequencingprimers' in t,True);c('ActualcitedTScanIIversionverified','10.1016/j.cell.2023.10.024' in t,True)
 c.report()
if __name__=='__main__':selftest()
