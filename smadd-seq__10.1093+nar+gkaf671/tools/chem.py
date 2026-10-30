"""Private SMAdd-seq selected spheroplast source component, not ONT recipe reconstruction."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
class B:
    def __init__(self,title,authority,start):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,K.MOLECULAR,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source copresent reaction mixture; incomplete/unreacted classes and carryover unassigned'):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']='Positive source operation; historical component formulation/whole molecular products unassembled, not source absence.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def split(self,key,label,branches):
        yy=tuple(self.g.product(key+'_'+name,txt,K.MOLECULAR)for name,txt,disp in branches);self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),yy)
        for y,(_,txt,disp)in zip(yy,branches):
            if disp is not None:self.g.mark_terminal(Terminal(y.id,disp,txt))
        self.x=yy[0];return yy
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Source retained physical fraction; recovery/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated wash/workup fraction',K.WASTE);self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Source separated workup; no perfect specificity/yield.');self.x=y

ANGELICIN_STOCK_UL=50
ANGELICIN_STOCK_UM=10000
DILUENT_UL=490
SPHEROPLAST_UL=500
TOTAL_UL=SPHEROPLAST_UL+ANGELICIN_STOCK_UL+DILUENT_UL
LITERAL_CALCULATED_UM=ANGELICIN_STOCK_UL*ANGELICIN_STOCK_UM/TOTAL_UL
DECLARED_UM=500
WELL_UL=250
UNTRANSFERRED_UL=TOTAL_UL-3*WELL_UL
def protocol_graph():
    b=B('SMAdd-seq selected spheroplast preparation','PublishedPMC12276009main84–98/SI FigS1/S2. Upstreamref30preparationoutsideboundary; historicalONTkitinternalsunassembled.','Sourcepreparedspheroplastsfrom250mllogculture; exactsourcepredecessorchemistryunassembled');g=b.g
    b.clean('spheroplast_wash','Source2xTMSorbitolwash/resuspend500ul1.2Msorbitol20mMTris8/1mMMgCl2')
    b.op('add_angelicin','Source50ul10mMangelicinin490ulbufferaddedto500ulsuspension;5minRT',S.GENERIC,(b.stock('angelicin','SourcephysicalangelicinSigmaA0956 stockin100percentethanol; actualtagprecursor',K.MOLECULAR,R.STRUCTURAL),b.cond('mix_conditions','Source5minRT; declared500uM vs literal1040ul/~480.8uM retained as sourcearithmeticdifference')),product='Sourceangelicin/spheroplastmixture, pre-UVnoncovalent/intercalated/freeclassesunassigned; notyetassumedcovalentadduct')
    branches=[('w'+str(i),'Source250ulwell'+str(i),None)for i in range(1,4)]+[('remainder','Nominal290uluntransferredsourcepreirradiationmixture; actualfateunassigned',Disposition.UNRESOLVED)]
    wells=b.split('plate_allocation','Literal1040ulsourceinput→3x250ulwells+nominal290uluntransferred',branches);recovered=[]
    for j,well in enumerate(wells[:3],1):
        b.x=well
        for i in range(1,8):
            b.op('uv'+str(j)+'_'+str(i),'Source365nmUVA5min onice/water, round'+str(i)+'of7',S.CHEMICAL_CONVERSION,(b.cond('uva'+str(j)+'_'+str(i),'SourceStratalinker2400power5.0/4EIKOF15T8BLbulbs; photochemicalconditions'),),product='SourceDNAangelicinadduct/unreacted/crosslinkedsideclasses; position/occupancy/yieldunassigned, notpuremonoadductorperfectaccessibility')
            b.op('dark'+str(j)+'_'+str(i),'Sourcealternating5mindarkincubationonice/water',S.THERMAL_TREATMENT,(b.cond('dark_conditions'+str(j)+'_'+str(i),'Sourcedark/icewater5min'),))
        b.op('well_recovery'+str(j),'SourcesamewellrecoverywithicecoldTMSorbitol100ulwashdescriptor; perwellallocationunassembled',S.GENERIC,(b.cond('recovery_buffer'+str(j),'Sourcewashrecoverybuffer; exactdistributedvolumeunassembled'),));recovered.append(b.x)
    pooled=g.product('irradiated_pool','Actualpooled3wellcontents+sourcewashrecovery',K.MOLECULAR);g.react('pool_wells','Sourcepoolallirradiatedwellsinto1.5mllowadhesiontube',S.POOL,tuple((x,R.TEMPLATE)for x in recovered),(pooled,));b.x=pooled
    b.clean('hmw_extract','SourceMonarchHMWTissueT3060Lextraction/wideboretips; retainHMWDNA, exactinternalstocks/physicalsubfractionsunassembled')
    b.split('dna_qc','ActualQubitBroadRange/TapeStationQC aliquotsseparatefromDNAproductionstock',[('production','RemainingHMWDNAstock',None),('qubit','SourceQubitBroadRangeQCaliquot; volumeunassigned',Disposition.QC),('tapestation','SourcegenomicDNAQCaliquot; volumeunassigned',Disposition.QC)])
    b.split('library_input','Source3–4ugHMWDNAselectedforONTkit; otherstockretained',[('selected','Source3–4ugselectedDNA',None),('remainder','UntransferredHMWDNAstock; actualmass/fateunassigned',Disposition.UNRESOLVED)])
    b.op('ont_kit','SourcehistoricalSQKLSK110librarypreparation; internalrepair/ligation/purificationunassembled',S.MAGIC,(b.stock('ont_preparation','ActualhistoricalONTkitpreparationmolecularstock; exactcomponents/endsunassembled',K.MOLECULAR,R.STRUCTURAL),),product='SourceONTproductionlibrary; adapterends/fullDNAproductsandresidualclassesunassembled')
    b.split('loading','Sourceapproximately1.5ugfinishedlibraryselectedforflowcell',[('selected','Sourceapproximately1.5ugloadingaliquot',None),('remainder','Untransferredfinishedlibrarystock; actualmass/fateunassigned',Disposition.UNRESOLVED)])
    lib=g.product('library','SelectedangelicinDNAONTlibraryaliquot',K.LIBRARY);g.react('load','SourceR9.4.1FLOMIN106MinION24h/MUXscan6h',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceMinIONR9.4.1/LSK11024h6hMUX; exacthistoricaladapter/motorgeometryunassembled. Nuclei/controltrainingbranchesnotmerged; nanoporeadductcallnotchemicalbaseconversion.'))
    g.add_source_item(SourceItem('unassembled','Sourcevolume/concentrationarithmeticdifference500uMvs1040ulliteralmixture, interpretedassourcevariantsnotnormalized.','Wholelibrarygeometry/upstreamref30kitdetailsunassembled; noabsenceaudit.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('Literalnominaluntransferredmixture',UNTRANSFERRED_UL,290);c('Declaredconcentrationnotrewritten',DECLARED_UM,500);c('Literalrecipeconcentrationdifferent',round(LITERAL_CALCULATED_UM,3),480.769);c('Noangelicinreadditionduringphotocycles',sum(any(p.material=='angelicin'for p in s.inputs)for s in g.steps.values()),1);c('Threewellphotocyclecount',sum(k.startswith('uv')for k in g.steps),21);c('Eachwellactual7darkintervals',sum(k.startswith('dark')for k in g.steps),21);c('Pooltakes3producedwellfractions',len(g.steps['pool_wells'].inputs),3);c('ProductionDNAexcludesQC',g.steps['library_input'].inputs[0].material,'dna_qc_production');c('KitconsumesselectedDNAaliquot',g.steps['ont_kit'].inputs[0].material,'library_input_selected');c('Loadnotwholekitstock',g.steps['load'].inputs[0].material,'loading_selected');c('Remainderunresolved',g.terminals['plate_allocation_remainder'].disposition,Disposition.UNRESOLVED)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'nar_gkaf671_PMC12276009.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1093/nar/gkaf671')):
        t=p.read_text();c('Originalliteral490ulbuffer', '490' in t and '500' in t and 'seven'in t,True)
    c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
