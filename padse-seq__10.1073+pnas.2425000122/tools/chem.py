"""Private PADSE-seq fixed-cell~30nmHER2/proximityprotein reporter preparation."""
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

from chemdraw import feature
from pcr import PCRGeometry
SP='CTT'+'A'*9+'GCTGAGGCTGAGGCTGAGGCACATCTTC'
BP={'her1':'CTTCAGTTAACCAAATAAACAACAGTCTTTAAACGTACGAGTCTTTATGTGCCTCAGCCAAAA','epcam':'ATATGTATAACCAAATAAACAACAGTCTTTAACTTCGTGAGTCTTTATGTGCCTCAGCCAAAA','pdl1':'TAAATTGCAACCAAATAAACAACAGTCTTTAAAACAGTGAGTCTTTATGTGCCTCAGCCAAAA'}
ASSORI='ATGTGCCTCAGCCTCAGCCTCAGCAAAA'
TP_HEAD='CTTAGCGACGGATCCTCAGCACATTCAGTTTCATCTCAAG'
TP_TAIL='GGCACATAAAGACTCAAAA'
TP=TP_HEAD+'N'*30+TP_TAIL
FP='AACCAAATAAACAACAGTCTTTAA'
RP='TCAGCACATTCAGTTTCATCTCAAG'
SOURCE_ORDERINGS={'sp':'N3-'+SP,**{k:'N3-'+v for k,v in BP.items()},'assori':ASSORI,'tp':TP,'fp':FP,'rp':RP}
def literal(seq,azide=False):
    identity=next((k for k,v in BP.items()if v==seq),None)
    ident=feature('padse.protein.'+identity,'feature_barcode','fixed',note='Sourceassignedprotein-IDprobebodycontext; minimalbarcode-onlyintervalunassembled')if identity else None
    s=Scene();s.strand('source',[Segment('source_whole',seq,feature=ident)],mod5='SourceazideN3(non-nucleotide),notNmer'if azide else'Unassigned');return s
def template():
    s=Scene();s.strand('template',[Segment('head',TP_HEAD),Segment('molecular_id','N'*30,placeholder=True,feature=feature('padse.tp.mid','inline_barcode','random',note='SourceTPmolecularID; ordinaryUMIusage/minimalcontextnotassumed')),Segment('tail',TP_TAIL)]);return s
def local_candidate(a,b,inferred=True):
    """Longest exact antiparallel span, conditional context; nopreferredbinding/end-product claim."""
    candidates=[];brc=revcomp(b)
    for i in range(len(a)):
        for j in range(len(brc)):
            n=0
            while i+n<len(a)and j+n<len(brc)and a[i+n]==brc[j+n]:n+=1
            if n:candidates.append((n,i,j))
    n,i,j=max(candidates);bs=len(b)-(j+n);s=Scene();ap=[];bp=[]
    if i:ap.append(Segment('a_left',a[:i],inferred=inferred))
    ap.append(Segment('a_pair',a[i:i+n],inferred=inferred))
    if i+n<len(a):ap.append(Segment('a_right',a[i+n:],inferred=inferred))
    if bs:bp.append(Segment('b_left',b[:bs],inferred=inferred))
    bp.append(Segment('b_pair',b[bs:bs+n],inferred=inferred))
    if bs+n<len(b):bp.append(Segment('b_right',b[bs+n:],inferred=inferred))
    s.strand('a',ap,label='INFERRED candidate placement, completeinputbasesretained; reactionoccupancyunassigned');s.anneal('b',bp,to='a',pair=('b_pair','a_pair'),unpaired=tuple(p.name for p in bp if p.name!='b_pair'));return s

def activated_bp_end():
    s=Scene();s.strand('activated',[Segment('unassembled_body','[source nicked BP body]',placeholder=True)],mod5='SourceAb/linker-associatedcontextunassembled',mod3='OH — sourceactivatedBPprimerend');return s

def protocol_graph():
    b=B('PADSE-seq selected fixedcellHER2/threeproximityprotein branch','PublishedPMC12012555main125–135/165 andSI133–187/124–131/TableS1originalp28/p29. NotSanger/origami/fluorescencevalidation. WholeSDAproducts/outsourcedNGSlibraryunassembled.','Source10,000cells oncollagenAcoveredcoverslip/PDMSchamberovernight; selectedsourcecelllineunassigned');g=b.g
    b.op('fix','Source4percentformaldehyde15minRT',S.CROSS_LINKING,(b.cond('formaldehyde','Source4percentformaldehyde'),));cells=b.x
    # Four independently sourcedAb-DNAfactory aliquots. Wholemanufacturerantibodystructure/dosingintoassayunassembled.
    probes={}
    for name,seq in [('her2',SP)]+list(BP.items()):
        b.x=g.input(name+'_antibody','Source10ul1mgmlanti'+name+'Abfactoryaliquot; manufacturerbodyunassembled',K.MOLECULAR)
        b.op(name+'_dbco','Source.5mMDBCOsulfoNHS1xPBS1hRTrotation',S.CHEMICAL_CONVERSION,(b.stock(name+'_dbco_stock','SourcephysicalDBCOsulfoNHSreactant',K.MOLECULAR,R.STRUCTURAL),b.cond(name+'_dbco_conditions','SourcePBS1hRT')),product='SourceDBCOmodified/unmodifiedantibodyclasses/freecrosslinker; substitution/yieldunassigned')
        b.clean(name+'_dbco_clean','SourceAmicon30kpurification; exactfreevsboundrecoveryunassigned')
        dna=g.input(name+'_azide_dna','Source15ul15uM5azideDNA '+('N3-'+seq),K.OLIGO,state=literal(seq,azide=True),chemical_state={'source_ordering':'N3-'+seq,'5prime':'azidechemicalgroup,N3isNOTthreerandomnucleotides','source_protein_identity':name})
        b.op(name+'_click','SourceazideDNA+DBCOAb4C1weekgentlerotation',S.CHEMICAL_CONVERSION,((dna,R.STRUCTURAL),b.cond(name+'_click_conditions','Source4C1week')),product='SourceAb-DNAconjugateclassesandremainingunconjugatedmaterial; exactlinkerstructure/attachmentstoichiometryunassembled')
        b.clean(name+'_conjugate_clean','SourceAmicon30kpurify/store4C')
        used,rest=b.split(name+'_assay_aliquot','SourceactualassayAdaaliquotselected; volumeunassigned',[('selected','Selected'+name+'Adaassayaliquot',None),('rest','Untransferred'+name+'Adafactorystock',Disposition.UNRESOLVED)]);probes[name]=used
    b.x=cells
    b.op('block','Source20ulperwell/QuickBlock1hRT',S.GENERIC,(b.cond('blocker','SourceQuickBlockbuffer'),))
    b.op('her2_recognition','SourceHER2Ada highsaltdilutionstock/overnight4C; finalNaClconcentrationunassigned',S.BINDING,((probes['her2'],R.STRUCTURAL),b.cond('her2_binding_conditions','SourceQuickBlockdilutionbuffer literal5MNaCldescriptor; actualfinaldilutionunassigned')))
    b.clean('her2_wash','Source3xnextstepbufferwash; retainfixedcells/Ada-boundclasses')
    ass=g.input('assori','Source500nMAssoriblockingoligo',K.OLIGO,state=literal(ASSORI))
    b.op('sp_block','SourceAssori500nM/1ulspermDNA.5mgml/4ulformamide/2ul20xSSC5h37C',S.HYBRIDIZATION,((ass,R.STRUCTURAL),b.stock('sperm_block','Sourcephysical1ulspermDNA.5mgmlblockingstock',K.MOLECULAR,R.STRUCTURAL),b.cond('sp_block_conditions','Source20ulassay/SSC/formamide5h37C; completebindingoccupancyunassigned')))
    b.clean('sp_block_wash','Source3xnextstepbufferwash')
    b.op('bp_recognition','SourceHER1/EpCAM/PDL1AdarecognitionfollowingHER2conditions; actualadditiontiming/dilutionsunassigned',S.BINDING,tuple((probes[k],R.STRUCTURAL)for k in BP)+(b.cond('bp_binding_conditions','SourcecorrespondingAdaantibody-bindingconditions'),))
    b.clean('bp_wash','Source3xnextstepbufferwash')
    b.op('activate','Source5UNtBbvCI1xCutSmart3h37C',S.DIGESTION,(b.stock('nickase_initial','Source5UNtBbvCIinitialactivation'),b.cond('activation_buffer','Source1xCutSmart3h37C')),product='Copresentactivated/inactiveantibody-associatedBP/freeproducts/reusedSPclasses; exactnickcoordinates/efficiencynotinferred')
    b.clean('activation_wash','Source3xnextstepbufferwash; releasedDNA/reagentsphysicalwashfractions')
    tp=g.input('tp','Source1uMTPcel randommolecularIDtemplate',K.OLIGO,state=template(),chemical_state={'source_ordering':TP,'source_molecular_ID_length':30,'ordinary_UMI_semantics':'unassigned'})
    b.op('tp_hybrid','SourceTPcel1uM/spermDNA/formamide2xSSC5h37C',S.HYBRIDIZATION,((tp,R.STRUCTURAL),b.stock('sperm_template','SourceSECONDfresh1ulspermDNA.5mgmlstock',K.MOLECULAR,R.STRUCTURAL),b.cond('tp_conditions','Source4ul100percentformamide/2xSSC5h37C')))
    b.clean('template_wash','Source3xnextstepbufferwash; occupancy/recoveryunassigned')
    templ=g.product('retained_tp','Samewellco-retainedTPtemplatecomponent',K.MOLECULAR,state=template());primer=g.product('retained_bp_primers','Samewellpreviouslyconjugated/activatedBP3OHprimerclasses; exactnick/wholeproductsunassembled',K.OLIGO,state=activated_bp_end());other=g.product('well_copresent','Othercopresentcells/Ada/SP/uncut/freeclasses; notphysicallyfractionated',K.MOLECULAR)
    g.react('sda_components','Trackcopresenttemplate/previousBPprimer/othercellcompartmentcomponents, NOnewprimerstockorphysicalseparation',S.GENERIC,((b.x,R.TEMPLATE),),(templ,primer,other),changes_molecule=False);b.x=templ
    b.op('sda','Source10UNtBbvCI/5UKlenow3to5exo−/500nMdNTP1xCutSmart3h37C',S.EXTENSION,((primer,R.STRUCTURAL),(other,R.STRUCTURAL),b.stock('nickase_sda','SourceactualSECOND10UNtBbvCIaddition'),b.stock('klenow','Source5UKlenowFragment3to5exo−'),b.stock('sda_nt','Source500nMdNTPmolecularsubstrate; perbase/totalconventionunassigned',K.MOLECULAR,R.STRUCTURAL),b.cond('sda_buffer','Source1xCutSmart3h37C')),product='Copresentcells/Ada/SPandSDA/bidirectionallyextendedreporters/proteinIDs+molecularIDs/offproductclasses; wholejunctionsandcopyyieldunassembled')
    sup=g.product('sda_supernatant','SourcecollectedSDAreportercontainingsupernatant',K.MOLECULAR);surface=g.product('surface_residual','Sourcecell/surface-associatedresidualfraction',K.MOLECULAR);g.react('collect_supernatant','SourcecollectSDA supernatant, notentirecellcompartment',S.EXTRACTION,((b.x,R.TEMPLATE),b.cond('collection','Sourcephysicalsupernatantcollection')),(sup,surface));g.mark_terminal(Terminal(surface.id,Disposition.UNRESOLVED,'Residualcells/Ada/DNAfateunassigned.'));b.x=sup
    b.op('heat_stop','Source80C20mininactivation',S.THERMAL_TREATMENT,(b.cond('heat_conditions','Source80C20min'),))
    b.op('pcr','Source1UHotStartTaq/100nMeachRPcelFP142s/500nMdNTP30cycles',S.PCR,(b.stock('pcr_rp','SourceRPcel100nM',K.OLIGO,R.STRUCTURAL,state=literal(RP)),b.stock('pcr_fp','SourceFP142s100nM; notfluorescentprimer/Sangerextensions',K.OLIGO,R.STRUCTURAL,state=literal(FP)),b.stock('taq','Source1UHotStartTaqDNApolymerase'),b.stock('pcr_nt','Source500nMdNTPactualPCRaddition',K.MOLECULAR,R.STRUCTURAL),b.cond('pcr_conditions','95C10min30x95C15s/60C30s/72C20s;72C5min')),geometry=PCRGeometry.unavailable(unknown_scene(),reason='WholesourceSDA/PCRproducts/primerlandingunassembled; sourceexactprimersretained.'))
    b.split('outsource_aliquot','Source5ngconvertedPCRDNAselectedforSangon',[('selected','Source5ngPCRDNAforsequencingprovider',None),('rest','UntransferredPCRstock; mass/fateunassigned',Disposition.UNRESOLVED)])
    b.op('provider_library','SourceSangonBiotechlibrarypreparation/NGS; exactplatform/recipe/rununassigned',S.MAGIC,(b.stock('provider_preparation','Actualoutsourcedsequencinglibraryconstructionmolecularpreparation; exactstock/adapter/indexendsunassembled',K.MOLECULAR,R.STRUCTURAL),),product='Sourceproductionsequencinglibrary; wholeadapted/indexedproductsunassembled')
    lib=g.product('library','SourceproviderproductionNGSlibraryaliquot; assayproteinproximitywithDNAreporters',K.LIBRARY);g.react('load','SourceNGSproviderboundary, no inventedIlluminaplatform/endgeometry',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceNGSafter5ngPCRDNA→Sangon; platform/runmask/libraryadapter/index/readlandingunassigned. Barcode/molecularIDcountsnotphysicalcleanup. Fixedcell30nmbranchonly; noorigami/Sangerrecipeinheritance.'))
    g.add_source_item(SourceItem('unassembled','Knownsource5azidewholeprobes/TPN30/proteinlabelsretained; fullAbDNAbackbone/nick/SDAproducts/providerlibraryunassembled.','Localcomputedpairingsareconditionalcandidateviews, notwholechemicalvalidation orsourceabsence.'));return g

def unknown_scene():
    s=Scene();s.strand('unknown',[Segment('unassembled','[whole source SDA/PCR product unassembled]',placeholder=True)]);return s
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import re
    c=Check();g=protocol_graph();c('AzideisnotthreeNnucleotides',sum(len(x.top)for x in literal(SP,True).strands['source'].segs),len(SP));c('Fourconjugatefactoryremainders',sum(k.endswith('_assay_aliquot_rest')for k in g.terminals),4);c('UseactualconjugatedHER2aliquot',any(p.material=='her2_assay_aliquot_selected'for p in g.steps['her2_recognition'].inputs),True);c('TPnt30retained',TP.count('N'),30);c('Source3OHBPprimersco-retainednotnewstock',any(p.material=='retained_bp_primers'and p.role==R.STRUCTURAL for p in g.steps['sda'].inputs),True);c('Initialnick5UandSDA10Useparatestocks',g.producer_of('sda_product').id,'sda');c('SDAactualnucleotides',any(p.material=='sda_nt'and p.role==R.STRUCTURAL for p in g.steps['sda'].inputs),True);c('Collectedsupernatantnotcellstock',g.steps['heat_stop'].inputs[0].material,'sda_supernatant');c('PCRactualNTcontribution',any(p.material=='pcr_nt'and p.role==R.STRUCTURAL for p in g.steps['pcr'].inputs),True);c('Provideronly5ngaliquot',g.steps['provider_library'].inputs[0].material,'outsource_aliquot_selected');c('Surfacefateunresolved',g.terminals['surface_residual'].disposition,Disposition.UNRESOLVED)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'pnas.2425000122_supp_pnas.2425000122.sapp.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1073/pnas.2425000122')):
        t=re.sub(r'\s+','',p.read_text())
        for k,v in SOURCE_ORDERINGS.items():c('Originalwholemodifiedordering '+k,v in t,True)
    c('SPAssoricandidateexactpaircomputable',len(local_candidate(SP,ASSORI).strands),2);c('Privatepartial',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
