"""Slic-seq: source-qualified SI-geometry partial model, published conflicts retained.

No reaction recipe is imported from its cited library predecessor (reference29).
"""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from chemdraw import Scene,Segment,revcomp,feature
from illumina import TRUSEQ_SMALL_RNA_RA3, INDEX2_PRIMER_RC, INDEX1_PRIMER
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,
    Terminal,Disposition,SequencingPlan,SourceItem)

AUTHORITY='PublishedPMC10484733Methods119/Data127/Results+actualPDFp2–4; publishedSIoriginalp13TableS1. Closedtwo-originalpass.'
SI_RNA_ARM=INDEX2_PRIMER_RC[:19]
MAIN_RNA_ARM=TRUSEQ_SMALL_RNA_RA3
DNA_ARM='A'+INDEX1_PRIMER[:21]
RT_PRIMER='A'+revcomp(SI_RNA_ARM)

@dataclass(frozen=True)
class EndClass:
    origin: str
    three: str
    def __post_init__(self):
        if (self.origin,self.three) not in {('old','vicinal-OH'),('old','oxidized-aldehyde'),('EndoV','OH')}:
            raise ValueError('Only source-selected qualitative end classes are represented')
    def oxidize_old(self):
        if (self.origin,self.three)!=('old','vicinal-OH'):
            raise ValueError('This local relation describes periodate acting on the selected original vicinal-OH end')
        return EndClass('old','oxidized-aldehyde')
    def require_three_oh(self):
        if self.three not in ('OH','vicinal-OH'):
            raise ValueError('Oxidized old end is not the selected 3OH ligation acceptor')
        return self

def selected_end_classes():
    return {'old_end':EndClass('old','vicinal-OH').oxidize_old(),
            'new_upstream_end':EndClass('EndoV','OH'),
            'new_downstream_five':'P',
            'scope':'Qualitative source EndoV cleavage class; other/unreacted ends and exact RNA substrates unassembled'}

def random_n10():
    return Segment('source DNA5 randomN10','NNNNNNNNNN',placeholder=True,
        feature=feature('slicseq.fragment_barcode','inline_barcode','random',
            note='Source describes duplicate discrimination; random-hexamer(N10) wording and exact processing window mapping unresolved. No cell/sample or RNA3N6 role inherited.'))

def adapter_scene(version,which):
    if version not in ('si','main') or which not in ('RNA3','DNA5'):
        raise ValueError('Choose explicit source version and RNA3/DNA5 component')
    if which=='RNA3':
        parts=[Segment('source RNA3 diversityN6','NNNNNN',note='Source says sequencing quality/diversity, not assigned a UMI role'),
               Segment('RNA3 fixed arm',SI_RNA_ARM if version=='si' else MAIN_RNA_ARM)]
        mod5='rApp activation, not backbone nucleotide'
    else:
        parts=[random_n10(),Segment('DNA5 fixed arm',DNA_ARM)];mod5='P'
    if version=='main':
        parts.append(Segment('appended terminalddC','C',note='Actual appended dideoxycytidine nucleotide with no3OH'))
    sc=Scene();sc.strand('adapter',parts,mod5=mod5,
        mod3='SpC3 non-nucleotide blocker' if version=='si' else 'ddC, no3OH')
    return sc

def adapter_length(version,which):
    return sum(len(s.top) for s in adapter_scene(version,which).strands['adapter'].segs)

def source_ordering(version,which):
    if which=='RNA3':bases='NNNNNN'+(SI_RNA_ARM if version=='si' else MAIN_RNA_ARM)
    elif which=='DNA5':bases='NNNNNNNNNN'+DNA_ARM
    else:raise ValueError('Unknown adapter')
    if version=='main':return ('/5rApp/' if which=='RNA3' else '/5Phos/')+bases+'/3ddC/'
    if version=='si':return ('5’rAPP-' if which=='RNA3' else '5’Phos-')+bases+'-3SpC3'
    raise ValueError('Unknown source version')

def rt_primer_scene():
    sc=Scene();sc.strand('RT',[Segment('unpaired5A','A'),Segment('RT complement',revcomp(SI_RNA_ARM))]);return sc

def si_rt_binding_scene():
    """Explicitly SI-qualified intended product, default inferred RNA/junction."""
    sc=Scene();sc.strand('template',[
        Segment('upstream RNA','[unassembled RNA]',placeholder=True,inferred=True),
        Segment('source RNA3 diversityN6','NNNNNN'),Segment('RNA3 fixed arm',SI_RNA_ARM)],
        mod3='SpC3 non-nucleotide blocker')
    sc.anneal('RT',[Segment('unpaired5A','A'),Segment('RT complement',revcomp(SI_RNA_ARM))],
        to='template',pair=('RT complement','RNA3 fixed arm'),unpaired=('unpaired5A',))
    return sc

def intended_acceptor_scene(kind):
    sc=Scene();sc.strand('acceptor',[Segment('whole '+kind,'[unassembled '+kind+']',placeholder=True,inferred=True)],mod3='INFERRED intended3OH')
    sc.ligation_end('acceptor','whole '+kind,kind='single-strand',polarity='3-prime',owner='substrate')
    return sc

def si_geometry(which):
    sc=adapter_scene('si',which)
    first=sc.strands['adapter'].segs[0].name
    sc.ligation_end('adapter',first,kind='single-strand',polarity='5-prime',owner='adapter')
    return AdapterLigationGeometry.known('single-strand',sc,intended_acceptor_scene('RNA' if which=='RNA3' else 'cDNA'),
        statement='Explicit SI-qualified donor geometry only. Main prints different RNA3 arm and appendedddC blockers. Actual-usedvariant/whole substrate/junction unresolved; intendedacceptor inferred by construction.')

def unavailable_scene():
    sc=Scene();sc.strand('unknown',[Segment('whole products','[unassembled products]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph('Slic-seq treatedHEK293T, explicitSI-qualifiedpartialgeometry',authority=AUTHORITY,representation='published_instance')
    x=g.input('cells','SourceHEK293Tcells; actualtranscripts/inosineoccupancy unassigned',K.BIOLOGICAL,primary=True)
    def reagent(key,label,kind=K.CONDITION,role=R.CONDITION,**kw):return g.input(key,label,kind,**kw),role
    def op(key,label,kind,inputs=(),product='Source mixture including unreacted/other/off-product classes; whole molecules/occupancy unassembled',chemical_state=None,**kw):
        nonlocal x
        y=g.product(key+'_product',product,K.MOLECULAR,chemical_state=chemical_state)
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        y=g.product(key+'_kept','Actual retainedfraction, imperfect molecular recovery/carryover unassigned',K.MOLECULAR)
        w=g.product(key+'_waste','Actual discardedworkup fractions; notideal specificity',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),reagent(key+'_conditions',label)),(y,w));g.discard(w,'Physical workup');x=y
    op('extract','SourceTRIzol totalRNAisolation',S.MAGIC,
       (reagent('trizol','SourcepreparedTRIzolreagent',K.MOLECULAR,R.STRUCTURAL),),reason='Sourceextractionpositive; detailedhistoricalkit/phases/recoveryunassembled.')
    clean('rna_isolation','SourceRNA recovery fromTRIzol extraction; exactworkup delegated')
    op('polya','SourceDynabeadsOligodT25 polyA enrichment',S.MAGIC,
       (reagent('polya_beads','SourcepreparedoligodT25beadstock',K.SELECTION,R.SELECTION),),reason='SourcepolyAcapturepositive; exactkitoperations/sourceamountsunassembled.')
    clean('polya_workup','Source physically recover polyA RNA; otherfractions removed withimperfect recovery')
    op('dnase','SourceTURBODNase treatment followingpolyAselection',S.MAGIC,
       (reagent('dnase_enzyme','SourceTURBODNasecatalyst',K.CATALYST,R.CATALYST),reagent('dnase_conditions','SourcedelegatedDNasetreatmentconditions')),
       reason='SourceDNAdigestionpositive; allproduct/off-productclasses andstockcleanupunassembled.')
    op('fragment','SourceAmbionRNAfragmentation70C10min',S.MAGIC,
       (reagent('fragment_reagent','SourcepreparedAmbionfragmentationreagent',K.MOLECULAR,R.STRUCTURAL),reagent('fragment_conditions','Source70C10minconditions')),
       reason='Sourcefragmentationpositive; exactRNAcleavages/endclassesunassembled.')
    op('pnk','SourceT4PNK endrepairbeforeperiodate; fullmixunassigned',S.MAGIC,
       (reagent('pnk_enzyme','SourceT4PNKcatalyst',K.CATALYST,R.CATALYST),reagent('pnk_conditions','SourcePNKconditionsunassembled; noATPadditionassumed')),
       reason='SourcePNKpositive; no5Pcreation/endclaim or reference29volumesinherited.')
    op('periodate','Source50mMNaIO4/25C30min/dark blocksoldRNA3ends',S.CHEMICAL_CONVERSION,
       (reagent('periodate_reactant','Preparedperiodatereactant, notcatalyst',K.MOLECULAR,R.STRUCTURAL),reagent('periodate_conditions','Source25C30mindark')),
       product='SourceoxidizedoldterminalaldehydeRNAclasses plusunreacted/backgroundRNA; no perfectoxidationassumed',chemical_state={'selected_old_end':'oxidized-aldehyde','other_classes':'unreacted/otherendsunassigned'})
    op('endov','SourceEndoV37C60min;10mMTris7.5/.5mMMn/50mMK/1mMDTT/5percentglycerol',S.DIGESTION,
       (reagent('endov_enzyme','SourceEndonucleaseV NEBM0305S',K.CATALYST,R.CATALYST),reagent('endov_conditions','SourceMn-containingbuffer; noMg/Ca substitution')),
       product='SelectedupstreamRNA3OH/downstreamRNA5P cleavageclasses atsecondbond3toI; oldoxidizedends/uncut/offtargetRNA retained asotherclasses, noexclusiveenrichmentclaim',
       chemical_state={'selected_upstream_three':'OH','selected_downstream_five':'P','position':'secondphosphodiester3toI','old_ends':'oxidized-class retained','scope':'qualitative classes, wholeRNA/atommap/yieldsunassembled'})
    a=g.input('si_rna3','ExplicitSI-qualifiedRNA3adapterN6+19ntarm/5rApp/3SpC3; actualusedvariantunresolved',K.ADAPTER,state=adapter_scene('si','RNA3'),chemical_state={'five_prime':'rApp','three_prime':'SpC3','backbone_length':adapter_length('si','RNA3')})
    op('rna3_ligate','SourceRNA3adapterligation, SI-qualifiedlocalgeometry; ligaseidentity/dose/conditionsunassigned',S.LIGATION,
       ((a,R.STRUCTURAL),reagent('rna3_ligase','SourceRNA3ligationcatalystunassembled; identitynotborrowedfromref29',K.CATALYST,R.CATALYST),reagent('rna3_conditions','SourceRNA3ligationconditionsunassembled')),
       geometry=si_geometry('RNA3'),product='OneRNA3ligationmixture: intendednewRNA3OH-to-activatedadapterproducts pluscopresentreleasedAMP/freeRNA/adapter/offproducts; wholejunctionunassembled',
       chemical_state={'copresent':'ActivationAMP released frombackbone butnotphysicallyseparated; fullbyproductfatesunassigned'})
    op('deadenylase','Source5primeDeadenylase beforeRecJf',S.MAGIC,
       (reagent('deadenylase_enzyme','Source5primeDeadenylase',K.CATALYST,R.CATALYST),),reason='Sourceadapterdeactivationpositive; no reactionvolumes/timesinherited.',product='Mixtureafterdeadenylase, copresentAMP/freeadapter/RNA/offproducts, actualmolecularfatesunassigned')
    op('recjf','SourceRecJf adapterprocessing afterdeadenylase',S.MAGIC,
       (reagent('recjf_enzyme','SourceRecJf',K.CATALYST,R.CATALYST),),reason='Sourcesaysremoveexcessadapter; fullsubstrates/specificity/dosesunassembled.',product='MixtureafterRecJf, copresentbyproducts/AMPandRNAproducts; noimaginedphysicalseparation')
    clean('rna3_clean','Sourcepurification afterRNA3linkerprocessing; reagent/amount/elutionunassigned')
    rp=g.input('rt_primer','Source20ntSSIII RTprimer',K.OLIGO,state=rt_primer_scene())
    op('rt','SourceSuperScriptIIIcDNAsynthesis, SI-qualified19bpprimerviewonly',S.MAGIC,
       ((rp,R.STRUCTURAL),reagent('ssiii','SourceSSIIIcatalyst',K.CATALYST,R.CATALYST),reagent('rt_nts','RTnucleotidecomponentrequiredbycopying; physicalstock/doseunassembled, notfreshinheritedref29addition',K.MOLECULAR,R.STRUCTURAL),reagent('rt_conditions','SourceRTconditionsunassembled')),
       reason='SourceRTpositive; completeoccupiedRNA/cDNA/endsandstockcompositionunassembled. No ExoI/RNAhydrolysis imported.')
    b=g.input('si_dna5','ExplicitSI-qualifiedDNA5adapterN10+22ntarm/5P/3SpC3; actualusedvariantunresolved',K.ADAPTER,state=adapter_scene('si','DNA5'),chemical_state={'five_prime':'P','three_prime':'SpC3','backbone_length':adapter_length('si','DNA5')})
    op('dna5_ligate','SourceDNA5adapterligation aftercDNAsynthesis, SI-qualifiedlocalgeometry',S.LIGATION,
       ((b,R.STRUCTURAL),reagent('dna5_ligase','SourceDNA5ligationcatalystunassembled; identitynotinherited',K.CATALYST,R.CATALYST),reagent('dna5_conditions','SourceDNA5ligationconditionsunassembled; noATPvolumeinvented')),
       geometry=si_geometry('DNA5'),product='IntendedcDNA3OH-to-adapter5Pclasses plusnativeRNA/cDNA/freeadaptor/offproducts; wholeproductsunassembled')
    clean('dna5_clean','Sourcepurification afterDNA5ligation; exacthistoricalworkupunassembled')
    prod=g.product('post_qpcr_stock','Actualremainingproductionstock after qPCR sampling, amountunassigned',K.MOLECULAR)
    qa=g.product('qpcr_aliquot','ActualsourceqPCRaliquot, amount/reagentmixunassigned',K.MOLECULAR)
    g.react('qpcr_split','SourceqPCRtoassesscycles takesactualQCaliquot',S.SPLIT,((x,R.TEMPLATE),),(prod,qa))
    g.mark_terminal(Terminal(qa.id,Disposition.QC,'SourceqPCRsample; wholeassay/reagentcomponentsunassembled, notproductionlibrary'));x=prod
    pf=g.input('neb_forward','SourceNEBprimercomponent, wholebody/indexidentityunassembled',K.OLIGO)
    pr=g.input('neb_reverse','OtherNEBprimercomponent, wholebody/indexidentityunassembled',K.OLIGO)
    stock=g.input('pcr_supplies','Aggregate productionPCR reagent supplies; physical packaging/numberofstocks/composition/volumes unassigned',K.MOLECULAR)
    pe=g.product('pcr_enzyme','Polymerase supply component; stock packaging/identity/dose unassigned',K.CATALYST)
    pn=g.product('pcr_nts','Nucleotide supply component; stock packaging/dose unassigned',K.MOLECULAR)
    pc=g.product('pcr_buffer','Condition supply component; stock packaging/composition unassigned',K.CONDITION)
    g.react('pcr_components','Represent copresent PCR supply components without asserting a premix or separate additions',S.GENERIC,((stock,R.TEMPLATE),),(pe,pn,pc),changes_molecule=False)
    op('pcr','SourceNEBprimerPCR, cycleschosenbyqPCR; conditionsunassigned',S.PCR,
       ((pf,R.STRUCTURAL),(pr,R.STRUCTURAL),(pe,R.CATALYST),(pn,R.STRUCTURAL),(pc,R.CONDITION)),
       geometry=PCRGeometry.unavailable(unavailable_scene(),reason='FullNEBprimerbodies/index/occupiedtemplate/wholecopiedproductsunassembled; nativeblockersnotassignedtonewDNAcopies.'))
    clean('library_clean','SelectedsourceNEBNextSamplePurificationBeads; lowmeltingagarosealternativeexcluded')
    loaded=g.product('loading_fraction','Actualsequencerloadingfraction; amountunassigned',K.MOLECULAR)
    rest=g.product('retained_stock','Untransferredlibrarystock; fateunassigned',K.MOLECULAR)
    g.react('loading_split','Actualloadingallocationfrompurifiedlibrarystock',S.SPLIT,((x,R.TEMPLATE),),(loaded,rest))
    g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Retainedstocknotautomaticallydiscarded'))
    lib=g.product('library','SourceSlicseqNGSlibrary, wholemolecules/readlandingunassembled',K.LIBRARY)
    g.react('load','SourceHiSeqXTenPE150',S.GENERIC,((loaded,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='SourcePE150 and randombarcodeprocessing; exactNEBprimers/index/readlanding/cyclewindowmappingunassembled. N10versusrandomhexamerwordingandRNA3N6diversityroleretained; computationalfiltersnotphysicalcleanup.'))
    g.add_source_item(SourceItem('adapter_discrepancy','Main119/PDFp2 RNA3N6+smallRNAarm/3ddC andDNA5N10+22ntarm/3ddC; SIoriginalp13TableS1 insteadRNA3N6+19ntarm/3SpC3 andDNA5samebody/3SpC3. BothgiveRT20.','LocalgeometryexplicitlySI-qualified; noactual-usedvariantclaim, sourceerrorassignment or main/SI merging. Mainalternativespreservedaswholemodifiedcomponents.'))
    g.add_source_item(SourceItem('barcode_discrepancy','Main127/PDFp3 says randomhexamer(N10), describesduplicate-discriminatingbarcodes; ResultsPDFp4 separatelysaysRNA3N6 improvessequencequality.','RetainphysicalN10andN6; N10inlinefeaturehasduplicate-discriminationpurposewithprocessing-mappinglimit, noRNA3N6UMI/cell/sampleassignment or target-ampliconinheritance.'))
    g.add_source_item(SourceItem('predecessor_limit','Main119citesreference29asmodifiedlibrarypredecessor.','No reference29reactionvolumes/ligaseidentities/PNKATP/RTprogram/hydrolysis/PCRrecipeinherited. Onlylocallyreadmain/SIfactsused; closedtwo-originalpass.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();root=sources_root()/HERE.parent.name
    for version,name in [('main','nar_gkad604_PMC10484733.xml.txt'),('si','nar_gkad604_supp_gkad604_supplemental_files_revised_NAR-SI_4.pdf.txt')]:
        p=root/name
        if have(c,Source(p,'https://doi.org/10.1093/nar/gkad604')):
            t=''.join(p.read_text().split()).replace('′','’').upper()
            for which in ('RNA3','DNA5'):c('Wholemodifiedsource '+version+which,source_ordering(version,which).upper() in t,True)
            c('WholeRTsource '+version,RT_PRIMER in t,True)
    blocked=selected_end_classes()['old_end']
    try:blocked.require_three_oh();bad=False
    except ValueError:bad=True
    c('Oldoxidizedendcannotbechosen3OHacceptor',bad,True)
    c('FreshEndoVendcanbechosen3OHacceptor',selected_end_classes()['new_upstream_end'].require_three_oh().origin,'EndoV')
    c('UnoxidizedcontrolendstillhasOH',EndClass('old','vicinal-OH').require_three_oh().origin,'old')
    c('MainDNAphysicalddCaddsC',adapter_length('main','DNA5')-adapter_length('si','DNA5'),1)
    c('MainRNAarmnotreplacedbySI',MAIN_RNA_ARM==SI_RNA_ARM,False)
    c('SIboundtemplatehasnoactivationAMP',si_rt_binding_scene().strands['template'].mod5,'')
    c('SIupstreamRNAinferreddefault',si_rt_binding_scene().strands['template'].segs[0].inferred,True)
    c('RNA3N6hasnosequencingidentifierrole',adapter_scene('si','RNA3').strands['adapter'].segs[0].feature,None)
    c('DNA5N10structuredinlinefeature',adapter_scene('si','DNA5').strands['adapter'].segs[0].feature.role,'inline_barcode')
    c('PeriodateafterPNK',g.steps['periodate'].inputs[0].material,'pnk_product')
    c('EndoVafterblocking',g.steps['endov'].inputs[0].material,'periodate_product')
    c('RNA3ligationafternewendcreation',g.steps['rna3_ligate'].inputs[0].material,'endov_product')
    c('RecJfafterdeadenylase',g.steps['recjf'].inputs[0].material,'deadenylase_product')
    c('ProductionretainsactualpostQCstock',g.steps['pcr'].inputs[0].material,'post_qpcr_stock')
    c('qPCRfractionisQC',g.terminals['qpcr_aliquot'].disposition,Disposition.QC)
    c('NoNEBprimerbodyborrowed',g.materials['neb_forward'].state,None)
    c('Privatepartial',bool(g.publication_problems()),True)
    c.report()

if __name__=='__main__':selftest()
