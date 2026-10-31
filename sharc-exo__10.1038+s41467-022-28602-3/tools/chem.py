"""Protocol-owned source-qualified incomplete molecular graph."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, revcomp, feature
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan
from pcr import PCRGeometry

def clean(g,key,x,label):
    condition=g.input(key+'_condition',label,K.CONDITION)
    kept=g.product(key+'_kept','Source recovered preparation; exact carryover unquantified',K.LIBRARY if x.kind is K.LIBRARY else K.MOLECULAR)
    rest=g.product(key+'_rest','Actual physically separated workup fractions; no perfect molecular partition inferred',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(kept,rest));g.discard(rest,'Actual source workup separation.');return kept

def op(g,key,x,label,kind=S.MAGIC,*,reagent='source reagents; precise geometry unassembled',catalyst=False,product=None,reason='Source operation known; internal molecular geometry and remnant fates unassembled.',molecular_reagent=False,additional=(),product_kind=K.MOLECULAR):
    stock=g.input(key+'_reagent',reagent,K.CATALYST if catalyst else (K.MOLECULAR if molecular_reagent else K.CONDITION))
    out=g.product(key+'_product',product or label+'; accompanying physical reaction remnants retained',product_kind)
    g.react(key,label,kind,((x,R.TEMPLATE),(stock,R.CATALYST if catalyst else (R.STRUCTURAL if molecular_reagent else R.CONDITION)))+tuple(additional),(out,),reason=reason);return out

def oligo_scene(name,sequence,*,mod5='',mod3=''):
    sc=Scene();sc.strand(name,[Segment(name,sequence)],mod5=mod5,mod3=mod3);return sc


from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry
from protocol_graph import Terminal,Disposition
SOURCE='VanDamme etal.2022 DOI10.1038/s41467-022-28602-3 publishedPMC8854666 Methods175–223/Fig2; publishedSI chemicalcharacterizationandSupplementaryFig1; reporting-summaryidentityonly.'
ADAPTER_BODY='AGATCGGAAGAGCGGTTCAG'
@dataclass(frozen=True)
class CrosslinkedPair:
    chains:tuple[str,str]
    sites:tuple[int,int]
    five_ends:tuple[str,str]
    three_ends:tuple[str,str]
    def __post_init__(self):
        if len(self.chains)!=2 or len(self.sites)!=2 or any(not s or set(s)-set('ACGU') for s in self.chains) or any(type(i)is not int or not 0<=i<len(s) for s,i in zip(self.chains,self.sites)):raise ValueError('selected pair needs two actual RNA chains and explicit2OHcrosslink sites')
        if len(self.five_ends)!=2 or len(self.three_ends)!=2 or any(e not in {'OH','phosphate','unassigned'} for e in self.five_ends+self.three_ends):raise ValueError('explicit per-chain terminal chemistry required')
    def trim_selected(self,ends,retained_three_ends):
        # Caller-selected exonuclease retention bounds, not predicted RNaseR endpoints.
        if len(ends)!=2 or any(type(e)is not int or not i<e<=len(s) for s,i,e in zip(self.chains,self.sites,ends)):raise ValueError('retained3primebounds must preserve selectedcrosslink sites')
        kept=tuple(s[:e] for s,e in zip(self.chains,ends));released=tuple(s[e:] for s,e in zip(self.chains,ends))
        return CrosslinkedPair(kept,self.sites,self.five_ends,retained_three_ends),released
    def ligate_selected(self):
        # This selected source-compatible3OH/5P class joins suppliedchain1tochain2;
        # population polarity/alternative ligations/efficiency remainunassembled.
        if self.three_ends[0]!='OH' or self.five_ends[1]!='phosphate':raise ValueError('selected Rnl1 class requires explicit3OH and5phosphate partners, not inferred unknown ends')
        return JoinedRNA(self.chains[0]+self.chains[1],len(self.chains[0]),(self.sites[0],len(self.chains[0])+self.sites[1]))
@dataclass(frozen=True)
class JoinedRNA:
    sequence:str
    junction:int
    bridge_sites:tuple|None
    def __post_init__(self):
        if not self.sequence or set(self.sequence)-set('ACGU') or not 0<self.junction<len(self.sequence):raise ValueError('selected ligatedRNA needs actual sequence andjunction')
        if self.bridge_sites is not None and (len(self.bridge_sites)!=2 or any(type(i)is not int or not 0<=i<len(self.sequence) for i in self.bridge_sites)):raise ValueError('actualcrosslinkpositionsrequired')
    def reverse_bridge_selected(self):
        if self.bridge_sites is None:raise ValueError('selected complete reversal requires an existingbridge')
        return replace(self,bridge_sites=None),'unassembled hydrolyzed crosslink-reagent residue products'
    def complementary_insert_reference(self):
        # A reference relation only, not an RT-primer-derived fullfirststrand or predictedRToutcome.
        return ''.join({'A':'T','C':'G','G':'C','U':'A'}[b] for b in self.sequence[::-1])

def adapter():
    sc=Scene();sc.strand('adapter',[Segment('body',ADAPTER_BODY),Segment('terminal_ddc','C',note='Actualsource3ddCcytidinenucleotide')],mod5='rApp',mod3='ddC:no3OH');sc.ligation_end('adapter','body',through='terminal_ddc',kind='single-strand',owner='adapter',polarity='5-prime');return sc

def unknown_ends():
    sc=Scene();sc.strand('unassembled',[Segment('ends','[exactRNAend/junctiongeometryunassembled]',placeholder=True)]);return sc

def graph_for_branch(branch='vivo',trim=True):
    if branch not in {'vivo','p4p6'} or type(trim)is not bool:raise ValueError('actualsourcevivoorP4P6branch andtrimmingcontrolrequired')
    g=ProtocolGraph('SHARC-exo — '+branch+('/exo' if trim else '/noexo'),authority=SOURCE,representation='published_instance')
    x=g.input('sample','Source'+('HeLa/HEK293Tcellpopulation' if branch=='vivo' else 'P4P6T7transcribed/denaturingPAGEpurifiedRNA; exactsourcebody notassigned'),K.MOLECULAR,primary=True)
    if branch=='p4p6':x=op(g,'refold',x,'SourceP4P6RNA90°C5min/ice/Mgaddition thenwarmingtoRT',S.THERMAL_TREATMENT,reagent='sourceP4P6refoldingbufferandMg2+')
    crosslink_environment=g.input('crosslink_environment','Actual source30minRT exposure; distinct concentration/control experiments, not one mixed dose',K.CONDITION)
    x=op(g,'crosslink',x,'SourceDPIcrosslink30minRT'+(';source0/5/12.5/25mMconditions' if branch=='vivo' else ';5µl50mMDPIinto45µlRNA;matchedDMSOcontrol'),S.CROSS_LINKING,reagent='Actual prepared CDI-activated dipicolinicacid(DPI)/DMSO molecular reagent; concentration/control choice unassigned in generic trace',molecular_reagent=True,additional=((crosslink_environment,R.CONDITION),),product='SelectedRNA2OHcrosslinkclasses plusmonoadducts/unmodifiedRNA/cellorbufferremnants; noexactdistance/site/yieldprediction')
    if branch=='vivo':
        x=clean(g,'cellwash',x,'Source crosslinksolutionremoval andtwicePBSwash')
        x=op(g,'lysis',x,'Source6MGuSCNlysis/manualshake/EDTA/PBS/needlepassage',S.EXTRACTION,reagent='sourceTNAextractionbuffers')
        x=op(g,'pk',x,'SourceProteinaseK1mg/ml37°C1h',S.DEPROTEINIZATION,catalyst=True,reagent='sourceProteinaseK',product='CrosslinkedRNAplusDNA/proteindigestionproducts/buffers; PKisnotphysicalremoval')
        x=clean(g,'tna_precipitate',x,'SourceNaOAc/phenol/isopropanolTNAprecipitation/two70%ethanolwashes')
        x=op(g,'dnase',x,'SourceTURBODNase37°C20min',S.DIGESTION,catalyst=True,reagent='sourceTURBODNase',product='RNAplusDNAdigestionproducts/reagentsuntilsubsequentphysicalrecovery')
        x=clean(g,'rna_recover',x,'SourceNaOAc/phenol/isopropanolRNArecovery')
    else:x=clean(g,'invitro_recover',x,'SourceP4P6RNAethanolprecipitation')
    mn=g.input('rnaseiii_mn','Actual5mM MnCl2 contribution',K.MOLECULAR)
    shortcut=g.input('rnaseiii_buffer','Actual1x supplied Shortcut buffer',K.CONDITION)
    x=op(g,'rnaseiii',x,'SourceRNaseIII/Mn2+37°C5minfragmentation',S.DIGESTION,catalyst=True,reagent='Actual10ul source RNaseIII catalyst',additional=((mn,R.STRUCTURAL),(shortcut,R.CONDITION)),product='CrosslinkedRNAfragmentcomplexes,unlinkedRNAandcleavageproductsretaineduntilworkup; exactendsunassembled')
    x=clean(g,'fragment_recover',x,'SourcephenolstopandNaOAc/GlycoBlue/ethanolrecovery')
    x=clean(g,'dd2d',x,'Source8%first/16%seconddenaturingureaMOPS2Dgels;exciseabovediagonal/gelcrushRNArecovery')
    if branch=='p4p6':
        selected=g.product('p4p6_selected_fraction','Source DD2D RNA aliquot selected for'+('RNaseR' if trim else 'noRNaseR')+' branch',K.MOLECULAR)
        companion=g.product('p4p6_companion_fraction','Paired source DD2D aliquot for'+('noRNaseR' if trim else 'RNaseR')+' comparison; separate physical material',K.MOLECULAR)
        g.react('p4p6_split','Published main221 divides DD2D RNA into two actual fractions',S.SPLIT,((x,R.TEMPLATE),),(selected,companion));g.mark_terminal(Terminal(companion.id,Disposition.UNRESOLVED,'Paired control library follows its separate graph_for_branch invocation; no aliquot duplicated into selected route.'));x=selected
    if trim:
        trim_atp=g.input('rnaser_atp','Actual5mM ATP molecular contribution to RNaseR',K.MOLECULAR)
        trim_buffer=g.input('rnaser_buffer','Actual1xRNaseR digestion buffer',K.CONDITION)
        x=op(g,'rnaser',x,'SourceRNaseR'+('45°C2/12/24h' if branch=='vivo' else '37°C2hP4P6fraction')+'; noexactcrosslink-to-stopdistanceimposed',S.DIGESTION,catalyst=True,reagent='Actual20U source RNaseR catalyst',additional=((trim_atp,R.STRUCTURAL),(trim_buffer,R.CONDITION)),product='SelectedtrimmedcrosslinkedRNAs plusreleasedresidueproducts/untrimmedRNAs/enzymes; noinstantphysicalpurification')
        x=clean(g,'trim_recover',x,'SourceNaOAc/GlycoBlue/ethanolrecoveryaftertrimming')
    lig=g.input('proximity_rnl1','Actual5ul sourceT4RNA ligase1 catalyst',K.CATALYST);y=g.product('proximity_mix','NewphosphodiesterjoinedRNAclassesretainchemicalcrosslinks,plusalternative/unligatedRNA/reagents',K.MOLECULAR)
    lig_atp=g.input('proximity_atp','Actual1ul0.1mM ATP contribution',K.MOLECULAR)
    lig_inhibitor=g.input('proximity_inhibitor','Actual1ulSuperaseIn contribution',K.MOLECULAR)
    lig_buffer=g.input('proximity_buffer','Actual2ul10xligation buffer',K.CONDITION)
    g.react('proximity_ligate','SourceRNAproximityligation BEFORE decrosslinking',S.PROXIMITY_LIGATION,((x,R.TEMPLATE),(lig,R.CATALYST),(lig_atp,R.STRUCTURAL),(lig_inhibitor,R.STRUCTURAL),(lig_buffer,R.CONDITION)),(y,),geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown_ends(),),'SourceRNAproximityligationknown; noexactterminalgeometryoruniqueproductpopulationassigned.'))
    x=clean(g,'lig_recover',y,'Sourceboil2min,centrifugeremovalofprecipitate,thenethanolrecovery')
    x=op(g,'decrosslink',x,'Source1xboratedecrosslinkbufferpH11,45°C2h',S.CHEMICAL_CONVERSION,reagent='sourceborate/heat',product='Selectedreversed2OHbridges retainnewRNAphosphodiesterjunctions; reagenthydrolysisproducts/remnantsremainuntilrecovery; quantitativecompletenessnotassumed')
    x=clean(g,'debridge_recover',x,'Sourceethanol/GlycoBluerecovery')
    a=g.input('ddc_adapter','Source/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/',K.ADAPTER,state=adapter());lig=g.input('adapter_rnl1','SourceT4RNAligase1RT3h',K.CATALYST);y=g.product('adapter_mix','3adapterligatedRNAplusfree/unligatedadapters/reagents; terminalddC remainsanactualCbaseplusblock',K.MOLECULAR)
    adapter_supplies=g.input('adapter_supplies','Actual2ulDMSO/5ulPEG8000/1ul0.1MDTT/1ulSuperaseIn; molecularinventory grouping, not premix',K.MOLECULAR)
    adapter_buffer=g.input('adapter_buffer','Actual2ul10xT4RNA ligase buffer',K.CONDITION)
    g.react('adapter_ligate','SourceadenylatedddCadapterligation',S.LIGATION,((x,R.TEMPLATE),(a,R.STRUCTURAL),(lig,R.CATALYST),(adapter_supplies,R.STRUCTURAL),(adapter_buffer,R.CONDITION)),(y,),geometry=AdapterLigationGeometry.known('single-strand',adapter(),statement='Publishedmain209rAppadapter/3ddC; populationRNAend/junctiongeometrynotfullyassembled.'))
    dead_buffer=g.input('adapter_digest_buffer','Actual3ul10xRecJf buffer contribution',K.CONDITION)
    dead_inhibitor=g.input('adapter_digest_inhibitor','Actual1ulSuperaseIn contribution',K.MOLECULAR)
    x=op(g,'free_adapter_digest',y,'SourceRecJfand5primeDeadenylase37°C1h',S.DIGESTION,catalyst=True,reagent='Actual2ulRecJf and1ul5primeDeadenylase catalyst contributions; not a documented premix',additional=((dead_buffer,R.CONDITION),(dead_inhibitor,R.STRUCTURAL)),product='RNAplusadapterdigestionproducts/reagents; noidealperfectligationoradapterclearance')
    x=clean(g,'adapter_clean',x,'SourceZymoRNAcleanConcentrator5')
    p=g.input('barcoded_rt','SourcebarcodedRTprimer; originalsequence/identifiergeometryunassembled',K.OLIGO);rt=g.input('ssiv','Actual1ul sourceSSIV catalyst; source temperature program stated in step',K.CATALYST);y=g.product('rt_mix','SourcecDNAproductswithoriginalRNAtemplates/unextendedprimer/enzymeremnants; exactprimerbackboneunassembled',K.MOLECULAR)
    nts=g.input('rt_nt','Actual1ul10mMdNTP contribution supplied before primer annealing',K.MOLECULAR)
    primed=g.product('primed','Actual RNA/1pmolbarcodedprimer/dNTP mixture after65C5minthenice; whole binding geometry unassembled',K.MOLECULAR)
    prime_condition=g.input('prime_condition','Actual65C5minthenrapidice',K.CONDITION)
    g.react('primer_anneal','Actual source RNA/primer/dNTP annealing preparation',S.HYBRIDIZATION,((x,R.TEMPLATE),(p,R.STRUCTURAL),(nts,R.STRUCTURAL),(prime_condition,R.CONDITION)),(primed,))
    buffer=g.input('rt_buffer','Actual4ul5xMn2+buffer:Tris-HCl8.3/CH3COOK/MnCl2 source1x50mM/75mM/1.5mM',K.CONDITION)
    dtt=g.input('rt_dtt','Actual2ul0.1MDTT contribution',K.MOLECULAR)
    inhibitor=g.input('rt_inhibitor','Actual1ulSuperaseIn contribution',K.MOLECULAR)
    g.react('rt','SourceSSIVreverse transcription AFTERchemicalreversal',S.REVERSE_TRANSCRIPTION,((primed,R.TEMPLATE),(rt,R.CATALYST),(buffer,R.CONDITION),(dtt,R.STRUCTURAL),(inhibitor,R.STRUCTURAL)),(y,))
    x=op(g,'rna_digest',y,'SourceRNaseHplusRNaseA/T1mix37°C30min',S.DIGESTION,catalyst=True,reagent='sourceRNaseH/RNaseA/T1mix',product='cDNAplusRNAdigestionproducts,reagentsandremainingtemplatesuntilcleanup')
    x=clean(g,'dna_cleanup',x,'SourceZymoDNAcleanConcentrator5')
    cat=g.input('circligaseii','Actual1ul sourceCircLigaseII catalyst; Mn2+ and buffer supplied separately',K.CATALYST);y=g.product('circles','SourcecircularcDNApopulationplusunclosed/reagentmaterial; originalprimerend/spacergeometryunassembled',K.MOLECULAR)
    circle_mn=g.input('circle_mn','Actual1ul50mMMnCl2 contribution',K.MOLECULAR)
    circle_buffer=g.input('circle_buffer','Actual10xCircLigaseIIbuffer contribution; volume unassigned',K.CONDITION)
    g.react('circularize','SourcecDNAcircularization',S.CIRCULARIZATION,((x,R.TEMPLATE),(cat,R.CATALYST),(circle_mn,R.STRUCTURAL),(circle_buffer,R.CONDITION)),(y,),geometry=CircularizationGeometry.unavailable('Sourceexplicitcircularization; barcodedprimer/end/backbonegeometryunreviewed, noordinaryspacer-freecircleinvented.'))
    supplies=g.input('pcr_supplies','Actual source PCR supply inventory; stockcount,packaging,precursorrecipe unreviewed',K.MOLECULAR)
    pf=g.product('pcr_forward','Actual source PCR forward functional primer class; body/landing unreviewed',K.OLIGO)
    pr=g.product('pcr_reverse','Actual source PCR reverse functional primer class; body/landing unreviewed',K.OLIGO)
    pe=g.product('pcr_enzyme','Actual source PCR catalyst class; identity unreviewed',K.CATALYST)
    pn=g.product('pcr_nt','Actual PCR nucleotide class; composition/dose unreviewed',K.MOLECULAR)
    pb=g.product('pcr_buffer','Actual PCR buffer class; composition/dose unreviewed',K.CONDITION)
    g.react('pcr_components','Expose actual functional PCR supplies; not a documented master stock',S.GENERIC,((supplies,R.STRUCTURAL),),(pf,pr,pe,pn,pb),changes_molecule=False)
    library=g.product('pcr_product','Actual source amplified DNA library and reaction companions; complete assembly unreviewed',K.LIBRARY)
    g.react('pcr_boundary','Source libraryPCR explicitly delegated toreference3; no exact recipe inherited',S.PCR,((y,R.TEMPLATE),(pf,R.STRUCTURAL),(pr,R.STRUCTURAL),(pe,R.CATALYST),(pn,R.STRUCTURAL),(pb,R.CONDITION)),(library,),geometry=PCRGeometry.unavailable(Scene(),reason='Actual source PCR operation known; directly cited precursor and exact primers/landing/cycles unreviewed.'))
    x=library
    x=clean(g,'size',x,'Source6%nativeTBEgel175bpandaboveclass/crushextraction')
    g.sequence(x,SequencingPlan(acquisition='unspecified',limitation=('Publishedmain221P4P6MiSeqexplicitlyreported; sourcekit/cycles/readprimersunassembled.' if branch=='p4p6' else 'Publishedmainin-vivoNGSexplicit; MiSeqreportedforseparateP4P6branch,notassignedtothisbranch. Exactcycles/platform/readlayoutunassembled.')))
    g.add_source_item(SourceItem('pending','Actualspecificextraction/ligation/librarypredecessors andbarcodedRTprimeroriginal','PublishedMethods181/201/217references1/2/3notfollowedinthree-originalcap; nosourceabsenceclaim.'))
    return g

def protocol_graph():return graph_for_branch()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();r=CrosslinkedPair(('ACGUAC','GGAUU'),(1,2),('phosphate','phosphate'),('OH','OH'));kept,released=r.trim_selected((4,4),('OH','OH'));joined=kept.ligate_selected();reversed_rna,products=joined.reverse_bridge_selected()
    c('exotrimretainsselectedanchor-bearingRNA',kept.chains,('ACGU','GGAU'))
    c('physicalreleasedresidueidentitiesremainseparate',released,('AC','U'))
    c('proximityligationcreatesRNAjunctionwithoutremovingcrosslink',joined.bridge_sites,(1,6))
    c('debridgepreservesligatedRNAsequenceandjunction',(reversed_rna.sequence,reversed_rna.junction),('ACGUGGAU',4))
    c('bridgeproductsnotmistakenforphysicalRNAdeletion',(reversed_rna.bridge_sites,bool(products)),(None,True))
    c('complementaryinsertreferenceisnotRNAdeletion',reversed_rna.complementary_insert_reference(),'ATCCACGT')
    try:r.trim_selected((1,4),('OH','OH'));bad=False
    except ValueError:bad=True
    c('selectedtrimmingcannotdeleteitsretainedcrosslinkanchor',bad,True)
    c('ddCisarealterminalcytidine',adapter().strands['adapter'].segs[-1].top,'C')
    g=protocol_graph()
    c('DPI is an actual molecular contributor',g.materials['crosslink_reagent'].kind,K.MOLECULAR)
    c('RNaseR ATP is separate from enzyme',g.materials['rnaser_atp'].kind,K.MOLECULAR)
    c('RT consumes already prepared RNA/primer/dNTP mixture',g.steps['rt'].inputs[0].material,'primed')
    c('Circle Mn2+ contribution is not an enzyme',g.materials['circle_mn'].kind,K.MOLECULAR)
    c('PCR functional components derive from source supply inventory',g.producer_of('pcr_forward').id,'pcr_components')
    c('RNaseRnegativecontrolnotgivenexotrim','rnaser' not in graph_for_branch(trim=False).steps,True)
    c('all source branches preserve geometry/source limits',all(graph_for_branch(b,t).publication_problems() for b in ('vivo','p4p6') for t in (True,False)),True)
    c('P4P6pairedfraction is physically retained separately',graph_for_branch('p4p6').terminals['p4p6_companion_fraction'].disposition,Disposition.UNRESOLVED)
    p=Source(sources_root()/'sharc-exo__10.1038+s41467-022-28602-3'/'s41467-022-28602-3_PMC8854666.xml.txt','get_sources --slugs sharc-exo__10.1038+s41467-022-28602-3;publishedMethods209/213')
    if have(c,p,label='originalchemistrytranscription'):
        t=p.path.read_text();c('sourceprintsactualrAppddCadapter','/5rApp/'+ADAPTER_BODY+'/3ddC/' in t,True);c('sourceprints10hRTafter45Cdecrosslink','42°Cfor10h' in ''.join(t.split()) and '2hat45°C' in ''.join(t.split()),True)
    c.report()
if __name__=='__main__':selftest()
