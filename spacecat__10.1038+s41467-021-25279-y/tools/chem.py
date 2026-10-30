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

SOURCE='PublishedPMC8371137.1 mainstaining258–262/tissue298–308/analysis312;mandatoryformalcorrectionPMC8417131readfirst:Fig3Cfurthestrightpaneltwoaccidentalwhitedotsremoved,noassertedchemistrycorrection. ActuallargeSITable1 notacquireddueto120MBdownloadcap,notpublisherinaccessibility.'
CORRECTION='10.1038/s41467-021-25672-7 corrects two white dots accidentally introduced in Fig3C furthest-rightpanel; original correction gives no recipe change.'
def read1_interval():
 s=Scene();s.strand('Source R1 identifier interval only',[Segment('12ntcellbarcode','N'*12,placeholder=True,feature=feature('spacecat_cell','cell_barcode','unknown')),Segment('8ntmolecularUMI','N'*8,placeholder=True,feature=feature('spacecat_umi','umi','random'))]);return s

def protocol_graph():
 g=ProtocolGraph('SPACECAT — selected calceinNVOC tissue-to-SeqWell route',authority=SOURCE,representation='published_instance')
 x=g.input('live_section','Source~1mm lungtumor tissue slice; livecell phenotypicROI,actualspecimenidentityunassigned; cellculture/organoid routesseparate',K.BIOLOGICAL,primary=True)
 x=op(g,'stain',x,'SourcecalceinNVOC staining≥30min37C5%CO2;sourceprinted6M DMSOstock diluted1:100PEG200 thenfinal1:5000fortissue;actualstockunitnotnormalized ornewmolarspeciespredicted')
 x=op(g,'photoactivate',x,'Sourcestainingmediumretainedfortissue;405nm12s pulses4x ROIillumination;fluorescenceencodesselectionhistory,notDNAbarcode')
 x=op(g,'dissociate',x,'SourcegentleMACS mechanical/enzymaticdissociation afterphotoactivation;exactkitinternalchemistryunassembled')
 x=op(g,'viability_stain',x,'Sourcecalceinviolet1:1000PBS0 indicatespostdissociationviability;notanothersequencingindex')
 x=cleanup(g,'facs',x,'SourceSonySH800ZselectcalceinNVOC+ANDviabilitycalceinviolet+;unstimulatednegativecontrolsets gate;collect500ulRPMI10%FBS/onice;purity/recoveryunquantified')
 x=op(g,'seqwell_load',x,'SourceChemGenesbarcodedcapturebeads/preparedPDMSSeqWellarray;15000cellsorallFACSisolatedsettlebygravity. Completehistoricalbead/TSO/adapteroligosunassembled')
 x=cleanup(g,'array_wash',x,'Source4xPBS0thenRPMIwashbeforesealing;unloadedcells/fractionsremainphysicalwaste')
 x=op(g,'seal_lysis_capture',x,'Sourceplasmafunctionalized0.01umporepolycarbonatemembrane,37Cdry40min;lysis5Mguanidiniumthiocyanate/1mMEDTA/1%betaME/0.05%sarkosyl20minRT;hybridize2MNaCl8%PEG40minRT')
 x=cleanup(g,'collect_beads',x,'SourcecollectmRNAcapturebeadsin2MNaCl/3mMMgCl2/20mMTris/8%PEGwashbuffer,pellet1000g3min;actualarray/solution fractionsseparated')
 tso=g.input('tso','Source5ul100uM5primeSeqWellTSO Table1; exactmodifiedoligounassembled',K.OLIGO,state=unknown('ActualTable1SeqWellTSO'))
 rt=g.input('rt','SourceMaximaHMinus10ul in200ulRTmix40ulH2O/40ul5xbuffer/80ul30%PEG/20uldNTP/5ulTSO/5ulRNaseinhibitor',K.CATALYST);xrt=g.product('rt_product','Sourcebead-boundcDNA/RNAcomplexeswithoriginalRNA,freeTSO/reagents andconversion/carryoverclasses; fullfirstcopy/TSgeometryunassembled',K.MOLECULAR);g.react('rt_ts','SourceRT/TS30minRT thenovernight52C endoverend;beadprimerisalreadyininputpreparation',S.MAGIC,((x,R.TEMPLATE),(tso,R.STRUCTURAL),(rt,R.CATALYST)),(xrt,),reason='KnowncompoundsourceRT/TSwithalready-carriedbeadprimer;fullTable1oligo/endsunassembled,noinventednewprimeraddition.');x=cleanup(g,'rt_wash',xrt,'Source1000g1minbeadwashes,1xTE0.5%SDS/2xTE0.01%Tween')
 x=op(g,'exoi',x,'Source170ulH2O/20ul10xExoIbuffer/10ulExoI37C1h endoverend;hydrolysisproductsretaineduntilwash',S.DIGESTION,'SourceExonucleaseI')
 x=cleanup(g,'exoi_wash',x,'SourceTE-SDSonce/TE-TWtwice physicalbeadwashes')
 x=op(g,'alkali',x,'Source500ul0.1MNaOH5minRT denaturation;bothstrandclasses/reagentsremainuntilwash')
 x=cleanup(g,'alkali_wash',x,'Source500ulTE-TWonce/TEonce;retainedbeads,otherphysicalfractions')
 primer=g.input('second_primer','Source2uldN-SMRToligoTable1;wholebody/randomprimingwindowunassembled',K.OLIGO,state=unknown('ActualdN-SMRTsourceprimer'))
 enzyme=g.input('klenow','Source5ulKlenow/200ulsecondstrandmix53ulH2O/40ulMaxima5xbuffer/80ul30%PEG/20uldNTP/2ulprimer',K.CATALYST);y=g.product('second_product','Sourcesecond-strand/cDNAcopiesplusbeads/originaltemplates/freeprimer/reagents;exactproductendsunassembled',K.MOLECULAR);g.react('second_strand','SourceKlenow/randompriming37C1h',S.EXTENSION,((x,R.TEMPLATE),(primer,R.STRUCTURAL),(enzyme,R.CATALYST)),(y,));x=cleanup(g,'second_wash',y,'Source2xTE-TW/1xTE/1xH2Obeadwashes then240ulresuspension')
 aliquots=g.product('wta_aliquots','Source24parallel10ulbeadaliquotsfrom240ul;allstockallocated,not24independentbiologicalcells',K.MOLECULAR);g.react('allocate_wta','Source240uldistributed24x10ul',S.GENERIC,((x,R.TEMPLATE),),(aliquots,),changes_molecule=False)
 x=pcr(g,'wta',aliquots,'Source24parallelwells,each40ulpremix25ul2xKAPA+14.6ulH2O+0.4ul100uMISPCR plus10ulbeads;95C3min;4x98C20s/65C45s/72C3min;12x98C20s/67C20s/72C3min;72C5min',(('SourceISPCRsingleprimeridentityTable1unassembled',unknown('ISPCRbodyunassembled;sourceonespeciesPCRgeometryunavailable')),))
 x=op(g,'wta_pool',x,'Source post-WTA workflow follows earlierSMARTseq2libraryworkup;parallelPCRproductcollection retained,exactpool-transferdetailsunassembled')
 x=cleanup(g,'wta_clean',x,'SourceAMPure0.6xfollowed0.8x postWTAcleanup')
 x=assay(g,'tapestation_qc',x,'Separate sourceTapestationhsD5000peak~1000bp/primerpeaksQC')
 x=assay(g,'qubit_qc',x,'Separate sourceQubitHighSensitivityDNAquantification')
 used,rest=split(g,'tag_allocation',x,('Source900pgcDNAforNexteraXT','UntransferredsourceWTAstock;mass/fateunassigned'));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'NotallamplifiedcDNAistagmented.'))
 x=op(g,'nextera_xt',used,'SourceNexteraXTDNASamplePreparationKitfrom900pg;exacthistoricalloadedends/indexprimers/tagmentationPCRunassembled')
 x=cleanup(g,'final_clean',x,'SourcefinalamplifiedlibraryAMPure0.6xfollowed0.8x purification')
 other=g.input('other_arrays','TwootherindependentlypreparedsourceSeqWellarraylibrariesforactual3arraypool;identities/massesunassigned',K.LIBRARY,primary=True);pool=g.product('array_pool','Source3arrayphysicallibrarypool',K.LIBRARY);g.react('pool_arrays','Source3SeqWellarraylibrariespooledforsequencing',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,))
 return finish(g,pool,'SourceNextSeq500/550HighOutputv2kit75cycles,PEreadstructureR1=20/R2=50/I1=8;customSeqWellCR1P Table1bodyunassembled. Main312R1=12ntcellbarcode+8ntUMI,notcalcein-spatialDNAbarcode. Exactmolecularprimerlanding/indexmembersunassembled. CorrectedFig3Cwhitedotsimage-onlynotice: '+CORRECTION+' Co-cultureSMARTseq2/organoidcollectionchoicesnotmergedintoKProute. Computation/QCfiltersnotphysicalDNAremoval.')

def selftest():
 from checks import Check,Source,have
 from storage import sources_root
 c=Check();g=protocol_graph();c('DyesarenotDNAspatialbarcode',tuple(seg.feature.role for seg in read1_interval().strands['Source R1 identifier interval only'].segs),('cell_barcode','umi'));c('SecondstrandafterphysicalNaOHwash',g.steps['second_strand'].inputs[0].material,'alkali_wash_kept');c('NexteraaliquotretainsWTAstock','tag_allocation_1' in g.terminals,True);c('WTAusesactual24aliquotcollection',g.steps['wta'].inputs[0].material,'wta_aliquots');c('NoexactTable1orfullreadgeometryclaim',bool(g.publication_problems()),True)
 corr=Source(sources_root()/'spacecat__10.1038+s41467-021-25279-y'/'s41467-021-25672-7_PMC8417131.xml.txt','Mandatoryformalcorrection')
 if have(c,corr,label='Formalcorrection'):
  c('ActualcorrectionisFig3Ctwoaccidentaldots','two white dots were inadvertently introduced' in corr.path.read_text(),True)
 c.report()
if __name__=='__main__':selftest()
