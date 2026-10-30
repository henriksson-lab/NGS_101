from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
class B:
    def __init__(self,title,authority,start,kind=K.BIOLOGICAL):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,kind,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source mixture including incomplete/unreacted/off-target classes; occupancy/recovery unassigned',reason=None):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']=reason or 'Positive source operation; complete occupied molecular inputs/products and historical formulation unassembled, not a missing-information conclusion.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label,extra=()):
        y=self.g.product(key+'_kept','Actual source retained physical fraction; yield/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated source workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label))+tuple(extra),(y,w));self.g.discard(w,'Actual workup fraction excluded; no perfect specificity or molecular yield claim.');self.x=y
    def split(self,key,label,selected,remainder,disposition=Disposition.UNRESOLVED):
        y=self.g.product(key+'_selected',selected,K.MOLECULAR);r=self.g.product(key+'_remainder',remainder,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r));self.g.mark_terminal(Terminal(r.id,disposition,remainder));self.x=y;return y

def unknown_scene():
    s=Scene();s.strand('unassembled',[Segment('source products unassembled','[unassembled]',placeholder=True)]);return s

ORDERS={'FWD': 'CAGATGTTTTCAGCTGTGACTGAT',
 'REV': 'CTTGGACAGTGGTACTGCAATAATT',
 'gRNA_A_template': 'AAGCTAATACGACTCACTATAGGTTTTCTTTTCACCTTCCCTGGTTTTAGAGCTAGAAATAGCAAG',
 'gRNA_B_template': 'AAGCTAATACGACTCACTATAGGGCAGAGGCTCGTTAATGTGTGTTTTAGAGCTAGAAATAGCAAG',
 'gRNA_VP_template': 'AAGCTAATACGACTCACTATAGGGTGGAGAACTGACTCCGCCTGTTTTAGAGCTAGAAATAGCAAG'}
AUTHORITY='Publishedauthor-institute2018mainMethods1602–1653/1670–1682;publishedSI1S1/S5p12;actualTable2MOESM4ZIP/XMLHS2C4:G4.3originals, no regular4C/2019companionrecipe inherited.'
def scene(label,seq=None):
    s=Scene();s.strand('source',[Segment(label,seq or'[whole source oligo unassembled]',placeholder=seq is None)]);return s

def components(b,key,label,specs):
    stock=b.g.input(key+'_stock',label,K.MOLECULAR);out=[];ins=[]
    for suffix,kind,role,desc,state in specs:
        x=b.g.product(key+'_'+suffix,desc,kind,state=state);out.append(x);ins.append((x,role))
    b.g.react(key+'_components','CopresentcomponentsofONEpreparedstock, notfreshindependentadditions',S.GENERIC,((stock,R.TEMPLATE),),tuple(out),changes_molecule=False);return tuple(ins)
def clean(b,key,label,beads='ActualsourceAMPureXPallocation'):
    b.clean(key,label,extra=(b.stock(key+'_selection',beads,K.SELECTION,R.SELECTION),))
def protocol_graph():
    b=B('MC-4C selectedfetal-liverHS2longcontactpreparation',AUTHORITY,'Actualsourcepreparedcrosslinkedfetal-liverchromatin; initialfixation/regular4Cbaseconditionsoutsidepartialgraph');g=b.g
    b.op('first_digest','SourceDpnIIdigestion500ul',S.DIGESTION,(b.stock('dpnii','ActualDpnII'),b.cond('first_digest_conditions','Completefirstdigestionconditionsunassembled; MboIHAP1alternativeoutsidegraph')))
    b.op('first_ligation','Sourcefirstligation2ml',S.MAGIC,(b.stock('first_ligation_preparation','Actualsourcefirstligationpreparation; catalyst/wholeends/conditionsunassembled',K.MOLECULAR,R.STRUCTURAL),))
    b.op('reverse_crosslink','Sourcecrosslinkreversal beforeDNAprecipitation',S.MAGIC,(b.cond('reverse_crosslink_conditions','Delegatedregular4Cconditionsunassembled'),))
    clean(b,'first_precip','Actual20ulPbeads/2mlisopropanol/2x80percentEtOH/resuspendsecondaryrestrictionbuffer','Actual20ulNucleoMagPbeads')
    b.op('hindiii','SourceHS2HindIIIsecondarydigestovernight/heatinactivate',S.DIGESTION,(b.stock('hindiii_stock','ActualHindIII'),b.cond('hindiii_conditions','Sourceappropriatebuffer/overnight; Man1ASacIalternativeoutsidegraph')))
    b.op('circularize','SourcedilutedligationDNA5ngul',S.MAGIC,(b.stock('second_ligation_preparation','Actualsourcecircularizationpreparation, catalyst/ends/wholeproductsunassembled',K.MOLECULAR,R.STRUCTURAL),))
    clean(b,'circle_precip','Pbeads10ulperml/2x80percentEtOH','ActualNucleoMagPbeadallocation');clean(b,'column','QiagenPCRpurificationremovesresidualPbeads','ActualsourceQiagenPCRpurificationcolumn/reagents')
    b.split('cas9_template','Actual20ugtemplateallocatedfrompreparedMC4Cstock','Actual20ugCas9input','ActualuntransferredMC4Cstock');template=b.x;rnps=[]
    for name in ('gRNA_A_template','gRNA_VP_template','gRNA_B_template'):
        key=name.lower();stock=g.input(key+'_pcr_source','Actualoverlap-PCRpreparedstock; Table2DNAforwardtemplate andcommonreversebodyunassembled',K.MOLECULAR)
        f=g.product(key+'_forward','ActualTable2DNAforwardprimer, NOTfunctionalRNAguide',K.OLIGO,state=scene(name,ORDERS[name]));r=g.product(key+'_reverse','Actualsourceoverlapreverseprimer; citedNakayamaoriginalbodyunassembled',K.OLIGO,state=scene('wholeoverlapreverseunassembled'));mix=g.product(key+'_reaction','Actualsame-stockPCRcatalyst/NT/buffer, wholeformulationunassembled',K.MOLECULAR)
        g.react(key+'_parts','Copresentcomponentsofoneactualoverlap-PCRstock',S.GENERIC,((stock,R.TEMPLATE),),(f,r,mix),changes_molecule=False);b.x=f
        b.op(key+'_template','SourceoverlapPCRmakesT7transcriptiontemplate',S.MAGIC,((r,R.STRUCTURAL),(mix,R.STRUCTURAL)),reason='Actual65ntDNAforwardorderretained; wholeoverlapreverse/template/scaffoldunassembled, notpretendedfullgRNA.')
        b.op(key+'_ivt','MegashortscriptT7RNAtranscription',S.MAGIC,(b.stock(key+'_ivt_stock','ActualsourceT7kitreaction, catalyst/rNTP/bufferinventoryunassembled',K.MOLECULAR,R.STRUCTURAL),))
        clean(b,key+'_rna_clean','4xAMPureRNApurification/DEPCwater')
        guide=g.product(key+'_rna','ActualpurifiedRNAguidepreparation; fullRNAsequence/endsunassembled',K.OLIGO,state=scene('RNAguidewholebodyunassembled'),chemical_state={'polymer':'RNA','source_DNA_forward':ORDERS[name],'full_RNA_body':'NotinferredfrompartialDNAforwardprimer'})
        g.react(key+'_guide_stock','SourceRNAguidepreparation usedforRNP, purity/transcriptclassesunassembled',S.GENERIC,((b.x,R.TEMPLATE),),(guide,),changes_molecule=False);b.x=guide
        b.split(key+'_600ng','Actual600ngguideallocation','Actual600ngRNAguideinput','ActualuntransferredRNAguidepreparation')
        b.op(key+'_rnp','Source600ngguide+15pmolCas9in300ulRT30min',S.MAGIC,(b.stock(key+'_cas9','Actualallocated15pmolCas9protein',K.CATALYST,R.CATALYST),b.cond(key+'_rnp_condition','SourceRNPpreincubation; completeoccupiedcomplexunassembled')));rnps.append(b.x)
    b.x=template;b.op('cas9_deplete','AddallthreepreincubatedRNPs to20ugDNA37C3–6h',S.MAGIC,tuple((x,R.STRUCTURAL)for x in rnps),reason='Sourcephysicalcleavage/depletionofuninformativeMC4Ctemplateclasses; wholefunctionalguideRNAs/genomictarget/ends/occupancyunassembled, notperfectspecificity.',product='Sourcecleaved/uncleavedmulti-contactDNAclasses; fullproducts/fatesunassembled')
    b.op('cas9_inactivate','Add1/25vol10percentSDS70C5min',S.QUENCHING,(b.cond('cas9_stop','ActualsourceSDSprogramme'),));clean(b,'cas9_clean','.6xAMPure')
    aliquots=tuple(g.product('pcr_template_'+str(i),'Actual100ngtemplateallocationfor25ulPCR'+str(i),K.MOLECULAR)for i in range(96));rest=g.product('pcr_template_rest','ActualuntransferredCas9-treatedtemplatestock',K.MOLECULAR)
    g.react('allocate96','Source96separate100ngtemplatesfromonepreparedstock',S.SPLIT,((b.x,R.TEMPLATE),),aliquots+(rest,));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Actualstockremainder, quantity/yieldunassembled.'));products=[]
    for i,start in enumerate(aliquots):
        key='pcr'+str(i);b.x=start
        ins=components(b,key,'Actualone25ulsourcePhireIIPCRpreparation; enzyme/NT/bufferdosesunassembled',[('f',K.OLIGO,R.STRUCTURAL,'ActualHS2forwardprimer',scene('HS2forward',ORDERS['FWD'])),('r',K.OLIGO,R.STRUCTURAL,'ActualHS2reverseprimer',scene('HS2reverse',ORDERS['REV'])),('enzyme',K.CATALYST,R.CATALYST,'ActualPhireII',None),('nt',K.MOLECULAR,R.STRUCTURAL,'Same-PCRdNTPcomponents, sourceamountsunassembled',None),('buffer',K.CONDITION,R.CONDITION,'Same-PCRbuffer, fullformulationunassembled',None)])
        b.op(key,'98C30s/31cycles98C10s57C20s72C90s/72C5min',S.PCR,ins,geometry=PCRGeometry.unavailable(scene('WholeMC4Ctarget/ampliconunassembled'),reason='SourceHS2wholeprimersknown, actualconcatemertemplates/occupiedlanding/productsunassembled.'));products.append(b.x)
    pool=g.product('pcr_pool','Actualpoolofall96PCRreactions, nominal96x25ul',K.MOLECULAR);g.react('pool96','SourcepoolALL96PCRproducts',S.POOL,tuple((x,R.TEMPLATE)for x in products),(pool,));b.x=pool
    b.split('select_1ml','Actual1mLselectedfrompooledPCRstock','Actual1mLcleanup/size-selectioninput','ActualuntransferredPCRpoolremainder, quantity/fateunassembled');clean(b,'pcr_cleanup_size','.6xAMPuresimultaneouscleanupandsizeselection')
    b.op('pippin','SourcePippinHT1.5–8kbselection',S.MAGIC,(b.stock('pippin_selection','ActualsourcePippinHTselectionpreparation',K.SELECTION,R.SELECTION),),reason='Physicalsize-selectionbeforesequencing; completefragmentlengthdistribution/recoveryunassembled.')
    b.op('ont_finish','HistoricalONTkitfinish:SQK-NSK007/LSK108/LSK208, actualper-runchoiceunassembled',S.MAGIC,(b.stock('ont_kit','ActualhistoricalONTpreparedlibraryreagents; adapter/motor/NT/catalystformulationunassembled',K.MOLECULAR,R.STRUCTURAL),),reason='Sourcehistoricalkitoptions, no exactkit-flowcellpairing or borrowedcurrentONTconstructor.')
    b.split('loading','Actualsequencingallocation, amount/runpoolingunassembled','ActualNGScontribution','Untransferredlibraryclasses/stock, fateunassembled')
    lib=g.product('library','ActualMC4CMinIONlong-DNAlibraryclass, wholeends/moleculesunassembled',K.LIBRARY);g.react('load','SourceMinIONR9/R9.4/R9.5alternatives, actualassignmentunassembled',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceMinIONlongreadendpoint. Exacthistoricalkit-flowcellmapping/adapters/endproducts/readlengthunassembled; PCR/read-ligation/duplicatefilteringcomputational. Far-cis/transfragmentfingerprintsNOTsyntheticUMIs.'))
    g.add_source_item(SourceItem('scope','PublishedTable2HS2DNAguide-template-forwardoligos retainedasDNA, notinferredfullfunctionalRNA. Threeactual600ngRNA/15pmolCas9/300ulpreincubationsaddedto20ugMC4Ctemplate.','96actual100ngtemplatealiquots/96PCRpool/selected1mL andphysicalPippinrange; no copiedwholetemplateor computationalfilterdrawnaswetstep.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import zipfile,xml.etree.ElementTree as E
    c=Check();g=protocol_graph();c('96realPCRs',sum(x.kind==S.PCR for x in g.steps.values()),96);c('OneDNAstockallocated96plusremainder',len(g.steps['allocate96'].outputs),97);c('All96pooled',len(g.steps['pool96'].inputs),96);c('Sizecleanupusesselected1mL',g.steps['pcr_cleanup_size'].inputs[0].material,'select_1ml_selected');c('PCR0usesallocated100ng',g.steps['pcr0'].inputs[0].material,'pcr_template_0');c('PCR95usesdifferentallocated100ng',g.steps['pcr95'].inputs[0].material,'pcr_template_95');c('ThreepreincubatedRNPs',len(g.steps['cas9_deplete'].inputs),4);c('SDSafterphysicalCas9',g.steps['cas9_inactivate'].inputs[0].material,'cas9_deplete_product');c('PrimerDNAforwardnotcalledRNA',g.materials['grna_a_template_forward'].kind,K.OLIGO)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41588-018-0161-5_41588_2018_161_MOESM4_ESM.xlsx'
    if have(c,Source(p,'https://doi.org/10.1038/s41588-018-0161-5')):
        z=zipfile.ZipFile(p);ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};ss=[''.join(x.itertext())for x in E.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',ns)];cells={}
        for x in E.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//s:c',ns):
            v=x.find('s:v',ns)
            if v is not None and x.attrib.get('t')=='s':cells[x.attrib['r']]=ss[int(v.text)]
        for field,col in zip(('FWD','REV','gRNA_A_template','gRNA_VP_template','gRNA_B_template'),'CDEFG'):c('OriginalHS2wholeDNA '+field,cells[col+'4'],ORDERS[field])
    c('Privatepartialunready',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
