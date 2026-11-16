from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import (Scene,Segment,Construct,revcomp,feature,complement_segments,
                     MolecularState,Workflow,workflow_panel)
from protocol_model import role,closed_circle
from restriction import DPNII
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SourceCoverage,SequencingPlan,Terminal,Disposition
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
        ivt=components(b,key+'_ivt','Actual selected T7kit reaction inventory; packaging/stockcounts/dosesunassigned',[('enzyme',K.CATALYST,R.CATALYST,'Known source T7 RNA polymerase functional contribution; no dose inferred',None),('nt',K.MOLECULAR,R.STRUCTURAL,'Required rNTP contribution from same selected IVT kit supplies; formulation/dosesunassigned',None),('buffer',K.CONDITION,R.CONDITION,'Buffer conditions from same selected IVT kit supplies; formulation/dosesunassigned',None)])
        b.op(key+'_ivt','MegashortscriptT7RNAtranscription',S.MAGIC,ivt)
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
    g.add_source_item(SourceItem('scope','PublishedTable2HS2DNAguide-template-forwardoligos retainedasDNA, notinferredfullfunctionalRNA. Threeactual600ngRNA/15pmolCas9/300ulpreincubationsaddedto20ugMC4Ctemplate.','96actual100ngtemplatealiquots/96PCRpool/selected1mL andphysicalPippinrange; no copiedwholetemplateor computationalfilterdrawnaswetstep.'));return _publication_scope(g)

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
    c('Source-qualified page passes normal gate',g.publishable,True);c.report()



PAGE_FILE='mc-4c.html'
REVIEWED_BASE=(
 'source','dpnii','first_digest_conditions','first_digest','first_ligation_preparation','first_ligation','reverse_crosslink_conditions','reverse_crosslink','first_precip_condition','first_precip_selection','first_precip','hindiii_stock','hindiii_conditions','hindiii','second_ligation_preparation','circularize','circle_precip_condition','circle_precip_selection','circle_precip','column_condition','column_selection','column','cas9_template',
 'cas9_deplete','cas9_stop','cas9_inactivate','cas9_clean_condition','cas9_clean_selection','cas9_clean','allocate96','pool96','select_1ml','pcr_cleanup_size_condition','pcr_cleanup_size_selection','pcr_cleanup_size','pippin_selection','pippin','ont_kit','ont_finish','loading','load','library')
REVIEWED_GUIDE=('pcr_source','parts','template','ivt_stock','ivt_components','ivt','rna_clean_condition','rna_clean_selection','rna_clean','guide_stock','600ng','cas9','rnp_condition','rnp')
def _publication_scope(g):
    g.cover_source_item(SourceCoverage('scope','modelled','library'))
    g.declare_reviewed_scope(REVIEWED_BASE,AUTHORITY)
    # Source explicitly requires three RNP preparations and96separatePCRs.
    # These literal templates expand that physical source inventory, never inspect graph nodes.
    for prefix in ('grna_a_template','grna_vp_template','grna_b_template'):
        g.declare_reviewed_scope(tuple(prefix+'_'+suffix for suffix in REVIEWED_GUIDE),AUTHORITY)
    for number in range(96):
        g.declare_reviewed_scope(('pcr'+str(number)+'_stock','pcr'+str(number)+'_components','pcr'+str(number)),AUTHORITY)
    return g.require_publishable()

def contact_parts():
    return [role('viewpoint end toward contacts',inferred=True),
            Segment('first DpnII ligation junction',DPNII.site,inferred=True),
            role('contact fragment A',inferred=True),
            Segment('next DpnII ligation junction',DPNII.site,inferred=True),
            role('other contacts; order unknown',inferred=True),
            role('other viewpoint end',inferred=True)]


def class_scene(kind):
    sc=Scene()
    if kind in {'chromatin','fragments'}:
        for j,label in enumerate(('viewpoint fragment','contact fragment A','contact fragment B')):
            top=[role(label if kind=='fragments'else label.replace('fragment','genomic interval'),inferred=True)];bottom=complement_segments(top)
            sc.strand('f'+str(j),top,label='INFERRED — '+label)
            sc.anneal('b'+str(j),bottom,to='f'+str(j),pair=(bottom[0].name,top[0].name),
                      unpaired=(bottom[0].name,),label='Complement class')
        sc.note('b2','Separate restriction fragments; occupied contacts / ligation probability unassembled.'if kind=='fragments'else'Uncut genomic interval references in crosslinked chromatin, not pre-created restriction fragments; exact loci / contacts unassembled.')
        return sc
    if kind=='circle':
        return closed_circle('INFERRED — one reference strand of a duplex multi-contact circle',
            tuple(contact_parts()+[role('viewpoint interval between inverse primers; VP guide target',inferred=True)]),
            closure_label='secondary ligation closure; exact HindIII junction unassembled')
    if kind=='amplicon':
        parts=[Segment('HS2 FWD primer-derived end',ORDERS['FWD'],'r1')]+contact_parts()+[
            Segment('HS2 REV primer-derived end',revcomp(ORDERS['REV']),'r2')]
    elif kind=='ont':
        parts=[role('historical ONT end; unresolved',inferred=True),
            Segment('HS2 FWD boundary',ORDERS['FWD'],'r1')]+contact_parts()+[
            Segment('HS2 REV boundary',revcomp(ORDERS['REV']),'r2'),
            role('other ONT end; unresolved',inferred=True)]
    elif kind=='linearized':
        parts=[role('VP-cut to FWD interval',inferred=True)]+contact_parts()+[
            role('REV to VP-cut interval',inferred=True)]
    else:parts=contact_parts()
    bottom=complement_segments(parts);sc.strand('top',parts,label='INFERRED contact DNA')
    sc.anneal('bottom',bottom,to='top',pair=(bottom[-1].name,parts[0].name),
              unpaired=tuple(p.name for p in bottom if p.is_role_token()),label='DNA complement')
    if kind=='linearized':sc.note('bottom','VP-guide cleavage opens the circle between the inverse-primer sites; contact-spanning primer directions now converge on a linear template.')
    return sc


def molecular_walkthrough():
    w=Workflow(MolecularState('INFERRED — crosslinked chromatin before DpnII digestion',
                              tuple(class_scene('chromatin').rows())),track_graph=False)
    def step(title,kind,note):w.react(title,class_scene(kind).rows(),name=title,note=note)
    step('1a · DpnII restriction digestion in crosslinked chromatin','fragments',
         'Digest source crosslinked chromatin to create restriction fragments. The unknown interval references are not declared exact genomic cuts, positions or contacts. Subsequent joining is a separate operation below.')
    step('1b · Proximity ligation retains multi-contact concatemers','concatemer',
         'Published Results and Figure1 describe 3C-derived DNA concatemers containing fragments captured in spatial proximity. The fetal-liver route uses DpnII in 500 µl and a 2 ml first ligation. Canonical DpnII junctions are shown, but actual genomic fragment identities/orientations, ligation ends and regular4C conditions remain unassembled. Brackets are variable classes, not an invented genomic sequence.')
    step('2 · Reverse crosslinks and purify; digest with HindIII overnight','concatemer',
         'Release ligation products from the chromatin complexes, precipitate on NucleoMag P-beads and perform the secondary digest. The source uses HindIII for HS2; SacI for Man1A is a different branch. This trims the original long concatemers into PCR-compatible multi-contact templates; exact HindIII coordinates remain unassigned.')
    step('3 · Heat-inactivate the enzyme; dilute ligation at DNA 5 ng/µl','circle',
         'The secondary restriction fragments are circularized by diluted ligation. A circle can contain the viewpoint and several contact fragments. The drawn reference-strand order is illustrative; source-specific fragment lengths and exact closure bases are not established by this region-level map. Purify on P-beads and a Qiagen PCR column.')
    step('4 · Three Cas9/guide preparations cut viewpoint and flanking targets','linearized',
         'The Results distinguish two purposes: the viewpoint guide cuts between the outward-facing inverse-PCR primers to reduce rolling-circle amplification, while flanking-fragment guides suppress abundant uninformative templates. Thus Cas9 is not simply a perfect removal step. Each source guide allocation is 600 ng with 15 pmol Cas9 in 300 µl for 30 min at room temperature; add all three to 20 µg template for 3–6 h at 37°C. Actual occupied targets and cleavage products remain unassembled.')
    step('5 · SDS inactivate Cas9, purify and allocate 96 separate PCRs','linearized',
         'Add 1/25 volume of 10% SDS, 70°C 5 min, followed by 0.6× AMPure. Each of 96 physical PCRs receives 100 ng treated template in 25 µl. These are distinct template allocations, not 96 copies of one molecule. Cut/unreacted/off-target classes remain in the source mixture until the stated workup; neither guide sequence templates nor analysis filtering establish perfect depletion.')
    step('6 · HS2 FWD / REV viewpoint PCR, 31 cycles','amplicon',
         'Use the literal HS2 primer orders shown below with PhireII: 98°C 30 s; 31 cycles of 98°C 10 s / 57°C 20 s / 72°C 90 s; 72°C 5 min. The expected bounded product contains a primer-derived FWD end, captured contact fragments and the reverse complement of REV at the opposite reference-strand boundary. The outer VP-cut intervals remain on parental templates; PCR bounds the new product rather than degrading those intervals. No complete genomic primer-site alignment or exact amplicon is claimed.')
    step('7 · Pool all 96 products; clean 1 ml and select 1.5–8 kb','amplicon',
         'The source pools after PCR, selects 1 ml for 0.6× AMPure cleanup/size selection, then uses Pippin HT to retain 1.5–8 kb products. Source fragment filters affect which long molecules reach sequencing; they do not introduce molecular barcodes.')
    step('8 · Historical ONT preparation and MinION sequencing','ont',
         'The publication lists SQK-NSK007, SQK-LSK108 or SQK-LSK208 and the contemporary R9/R9.4/R9.5 flow cells. It does not establish their assignment for this displayed molecule. Bracketed adapter/motor/end classes mark that unresolved boundary; no current ONT adapter recipe is substituted. Long reads preserve several contact junctions; computational far-cis/trans-fragment fingerprints are not synthetic UMIs.')
    return w


def render_page():
    from page import head,caveat,table
    from chemdraw import panel
    from protocol_views import assemble_page,WorkflowSection,PlanLimitationSection
    g=protocol_graph()
    return assemble_page(g,[head('MC-4C'),'<main class="wrap" data-protocol-complete="true"><h1>MC-4C</h1><p>Long-read multi-contact chromatin libraries with viewpoint-specific amplification and Cas9 depletion of uninformative template classes.</p><p><a href="https://doi.org/10.1038/s41588-018-0161-5">Defining published article</a> · <a href="01_mc-4c.html">Research note</a></p>',
      caveat('Selected published fetal-liver HS2 preparation, beginning with prepared crosslinked chromatin. Regular4C base conditions, occupied genomic contact molecules, complete guideRNA bodies, ligation/PCR products and historical ONT adapters/kit-flowcell assignment remain unassembled. Partial DNA guide-template primers below are not reconstructed functional guideRNAs. The workflow records source operations and actual stocks, not complete reaction validation or perfect Cas9 depletion.'),
      '<h2>Step by step: contacts → circles → viewpoint-spanning reads</h2>',
      '<p>The visible workflow distinguishes a multi-contact concatemer, a secondary-digest circle, Cas9-opened template and the bounded PCR product. Its bracketed regions are inferred molecular classes, not measured genomic coordinates. Exact source primer boundaries are shown only where amplification contributes those primers.</p>',
      workflow_panel(molecular_walkthrough(),caption='MC-4C selected HS2 route: proximity ligation, circularization, Cas9 opening, inverse PCR and long-read library boundary'),
      '<h2>Source HS2 primer inventory</h2>',table(('DNA oligo','Whole source ordering'),tuple((name,ORDERS[name])for name in ('FWD','REV','gRNA_A_template','gRNA_VP_template','gRNA_B_template'))),
      '<h2>Read setup and scope</h2><p>The source lists SQK-NSK007, SQK-LSK108 or SQK-LSK208 and R9/R9.4/R9.5 flow cells according to the run. Exact kit-to-flowcell assignment and primer landing remain unassembled; no current ONT constructor is substituted.</p><p>Primary/secondary contact classification, validity filtering and far-cis/trans-fragment fingerprints are computational analyses. They do not physically remove molecules and are not synthetic UMIs. Supplementary Figure5 compares ligation-derived contact classes; it is not an extra wet cleanup.</p>',
      '<details><summary>Complete source preparation, all96PCRs and unresolved transitions</summary>',WorkflowSection('Actual physical preparation and fractions'),PlanLimitationSection('library'),'</details></main>'])

if __name__=='__main__':selftest()
