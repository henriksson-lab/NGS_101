"""Protocol-owned source-limited model; no finished-publication marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Scene, Segment, Construct, feature, revcomp, oligo, panel, workflow_from_protocol_graph, workflow_panel
from protocol_model import role, molecule
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, Terminal, Disposition, SourceItem, SequencingPlan
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from pcr import PCRGeometry
from page import head, info

def mixture(*species):
    sc = Scene()
    for n, (label, parts) in enumerate(species):
        sc.strand('s' + str(n), list(parts), label=label)
    return sc

def cleanup(g, ident, source, template, state):
    reagent = g.input(ident + '_condition', source, K.CONDITION)
    out = g.product(ident, 'recovered source preparation', K.MOLECULAR, state=state)
    waste = g.product(ident + '_waste', 'physically separated workup fractions; exact molecular partition unquantified', K.WASTE)
    g.react(ident + '_step', source, S.CLEANUP, ((template, R.TEMPLATE), (reagent, R.CONDITION)), (out, waste))
    g.discard(waste, 'Actual source physical separation; no perfect recovery/depletion asserted.')
    return out

def end_boundary(label):
    sc = molecule('unassembled ligation-ready geometry', (role(label),))
    return AdapterLigationGeometry('unavailable', 'unavailable', (sc,), label + '; no source-absence claim.')

def render_page():
    import html
    g = protocol_graph()
    parts = [head(TITLE), '<main class="wrap"><h1>' + html.escape(TITLE) + '</h1>',
             '<p><a href="' + NOTES + '">Research notes</a></p>', info(SOURCE),
             info('Private incomplete model. ' + CAVEAT)]
    if callable(globals().get('oligos')):
        parts.extend(['<h2>Source ordering oligos</h2>', *oligos()])
    for label, terminal in PREVIEWS:
        parts.extend(['<h2>' + html.escape(label) + '</h2>', workflow_panel(workflow_from_protocol_graph(g, terminal))])
    parts.extend([info(READOUT), '</main>'])
    return '\n'.join(parts)

def pop(label):
    return mixture(('source DNA',(role(label),)),('other physical material',(role('other original DNA/RNA/protein or reaction remnants; exact partition and recovery unquantified'),)))
def op(g,key,label,kind,x,out_state,reagent_label,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
    reagent=g.input(key+'_stock',reagent_label,reagent_kind)
    product=g.product(key+'_product',label,K.MOLECULAR,state=out_state)
    kw={'reason':'Source operation known; exact molecular substrate/product relation and mixture fractions unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    ports=[(x,R.TEMPLATE),(reagent,reagent_role)]
    if reagent_role==R.STRUCTURAL and kind==S.CHEMICAL_CONVERSION:
        cond=g.input(key+'_conditions',label+'; sourceconditions, notextramoleculardose',K.CONDITION);ports.append((cond,R.CONDITION))
    g.react(key,label,kind,ports,(product,),**kw)
    return product


def native(g,x,limitation):
    # Source acquisition may be stated without fabricated sequencing primers/cycle counts.
    g.sequence(x,SequencingPlan(limitation=limitation))

def pending(g,label):
    g.add_source_item(SourceItem('unreviewed_assembly',label,'Named source operations are represented; historical end/kit/primer assembly and complete original-version applicability unreviewed.'))



PAGE_FILE='fc-cet.html'
TITLE='fC-CET'
NOTES='01_fc-cet.html'
SOURCE='AssociatedmanuscriptNIHMS714853 mainMethods93–131/Fig1e; publishedNature21September2015updateonlyauthoraffiliations; DOI10.1038/nmeth.3569.'
CAVEAT='Physicalcyclized5fC-AI/tagremnant isdistinctfromPCRderivedTcall. MainprintsDBCO400mM; exactfinalVORunitapplicabilityunreviewed, notsilentlychangedtoµM. OriginalTELP/MightyAmpadaption internalgeometryunreviewed. Sourceupdatecorrectsaffiliations, notchemistry.'
READOUT='SourceHiSeq2500, suggested single-end100bporlonger; suggestionnotfixedsourceacquiredcycleclaim. SourcepreAIinput,AIinput,pulldownoutput separatelysequencedwithsamerecipe; exactAIinputclicksamplingpointunassembled. ValidationMiSeq150amplicons/fCAB/MABseparatefromgenome-wideTELProute.'
PREVIEWS=(('Originalgenome-wideAIcaptureandlibrary','library_output'),)
@dataclass(frozen=True)
class CyclizedDNA:
    reference:str
    sites:frozenset[int]
    stage:str='native'
    def __post_init__(self):
        if not self.reference or set(self.reference)-set('ACGT') or any(type(i)is not int or not 0<=i<len(self.reference) or self.reference[i]!='C' for i in self.sites) or self.stage not in {'native','ai','biotin','released'}:raise ValueError('actualC-origin5fCpositions and sourcechemicalstage required')
    def advance(self,stage):
        if (self.stage,stage) not in {('native','ai'),('ai','biotin'),('biotin','released')}:raise ValueError('selectedsourcechemicalorderrequired')
        return CyclizedDNA(self.reference,self.sites,stage)
    @property
    def modifications(self):return {i:{'native':'5fC','ai':'cyclized5fC-AI','biotin':'cyclized5fC-AI-clickeddisulfidebiotin','released':'cyclized5fC-AI-linkerremnant afterS-Sreduction'}[self.stage] for i in self.sites}
    def selected_first_copy(self):
        if self.stage=='native':raise ValueError('selectedmismatchcopyrequiresexplicitcyclizedsubstrate')
        return ''.join('A' if i in self.sites else {'A':'T','C':'G','G':'C','T':'A'}[b] for i,b in reversed(list(enumerate(self.reference))))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    dna=g.input('dna','Source mESC DNA-containing biological preparations; genotypes and biologicalreplicates separate, notalready-purifiedDNA',K.BIOLOGICAL,state=pop('sourcegenomicDNAwithnative5fC/othermodifiedcytosines'),primary=True)
    pk=g.input('extract_pk','ActualsourceproteinaseK',K.CATALYST);sds=g.input('extract_sds','ActualsourceSDS extractionreagent',K.MOLECULAR)
    x=g.product('extract_product','Source extractedDNA plusRNA/proteolytic/cellularremnants',K.MOLECULAR,state=pop('DNA,RNA/proteolytic/cellularremnants'))
    g.react('extract','SourceSDS/proteinaseKgenomicDNApreparation',S.MAGIC,((dna,R.TEMPLATE),(pk,R.CATALYST),(sds,R.STRUCTURAL)),(x,),reason='Positive source extraction; biologicalgenotype/replicatepreparation notonepurifiedDNAstock, fullformulation/productpartition unassembled.')
    x=cleanup(g,'organic','Sourcephenol/chloroformextractionandethanolprecipitation',x,x.state)
    x=op(g,'fragment','SourceNEBdsDNAFragmentaseapproximately100–400bp',S.FRAGMENTATION,x,x.state,'sourceFragmentase',K.CATALYST,R.CATALYST)
    preai=g.product('preai','SourcepreAIinputcomparisonDNA',K.MOLECULAR,state=x.state);target=g.product('label_input','SourceDNAallocatedtolabel/capture',K.MOLECULAR,state=x.state)
    g.react('preai_split','SourceparallelpreAI/input-labelling allocation; physicalstock/replicatemappingunassembled',S.MAGIC,((x,R.TEMPLATE),),(preai,target),reason='Sourcecomparisons and2biologicalreplicatesknown; notproof ofa singlecommon-tubephysicalsplit oramounts.')
    target=op(g,'ai','SourceAI suspension100mMMESpH6,37C24h850rpm',S.CHEMICAL_CONVERSION,target,pop('cyclized5fC-AIadduct-bearingDNA andunlabelled/backgroundmolecules'),'ActualAI5-(2-azidoethyl)-1,3-indandione molecularcontributor; templateadductnotDNA-T',K.MOLECULAR,R.STRUCTURAL)
    target=cleanup(g,'ai_cleanup','SourcegenomicDNAQIAquickPCRpurification',target,target.state)
    ai_input=g.product('ai_input','SourceAIinputcomparison; exactsamplingpointrelativeclickunreviewed',K.MOLECULAR,state=pop('AI-inputsourceDNA; exactbiotinconjugationstateunassigned'))
    capture_input=g.product('click_input','SourceAI-labeledDNAforbiotinconjugation',K.MOLECULAR,state=target.state)
    cond=g.input('ai_input_allocation','SourceAIinputversuspulldownpreparations; exactsamplingtimingunreviewed',K.CONDITION)
    g.react('ai_control','SourceAIinputcomparisonallocation',S.MAGIC,((target,R.TEMPLATE),(cond,R.CONDITION)),(ai_input,capture_input),reason='SourceAIcontrolidentityknown; precisephysicalclicksamplingtimingunassembled, no same-stateassertion.')
    x=op(g,'click','SourceDBCO-S-S-PEG3-biotinA112-10,mainprints400mM37C2h; unitnotnormalized',S.CHEMICAL_CONVERSION,capture_input,pop('cyclized5fC-AI-clickedS-SbiotinDNA plusunreactedspecies/reagents'),'ActualDBCO-S-S-PEG3-biotin molecularcontributor; main400mM finalversionapplicabilityunreviewed',K.MOLECULAR,R.STRUCTURAL)
    x=cleanup(g,'click_cleanup','SourceQIAquickPCRpurificationafterclick',x,x.state)
    beads=g.input('beads','SourceMyOneStreptavidinC1',K.SELECTION)
    bound=g.product('bound','SourcecapturedDNAduplex/strandcomplexes',K.MOLECULAR,state=x.state);other=g.product('unbound','Companionunboundcapturematerial',K.MOLECULAR,state=pop('unboundDNA/reagents; exactdisposalunreviewed'))
    g.react('capture','Sourcebiotin-DNAcaptureB&WpH7.5plus0.1%Tween20',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Noquantitativecapture/disposalinferred.'))
    x=cleanup(g,'capture_wash','Sourcecanonicalwash5times then1xSSCwash',bound,bound.state)
    x=op(g,'alkali','Sourcefresh0.15MNaOHRT10min;duplexdenaturationnotcovalentbiotinrelease',S.MAGIC,x,pop('bead-associatedbiotinadductDNAstrands anddissociatednonbiotinDNA/othermaterial'),'Actualsource0.15MNaOH',K.MOLECULAR,R.STRUCTURAL)
    x=cleanup(g,'strand_wash','Sourcewash0.1MNaOH thenB&W thenTE; removesdissociatedcompanionstrandswithoutperfectspecificity',x,pop('retainedbiotin-taggedAI-DNAstrandclasswithpossiblebackground; exactendsunassigned'))
    x=op(g,'dtt','Source50mMDTT37C2h,S-Slinkercleavage;5fC-AIadductnoterased',S.CHEMICAL_CONVERSION,x,pop('releasedAI-DNAwithsourcecleavedlinkerremnant plusbeads/tagfragments/reagents'),'ActualsourceDTT molecularreductant',K.MOLECULAR,R.STRUCTURAL)
    condition=g.input('release_condition','SourcesupernatantcollectionafterDTT',K.CONDITION)
    released=g.product('released','SourcecollectedDNAcontainingsupernatant; publishedSI names5fC-AI-SHclass, notunmodifiedC',K.MOLECULAR,state=pop('AI-DNAwithcleavedlinkerremnant andpossiblecarryover'))
    residual=g.product('bead_residual','Sourcepostreleaseresidualbeadmaterial',K.MOLECULAR,state=pop('beads/biotinlinkerfragments/possibleunreleasedDNA'))
    g.react('collect','SourcephysicallycollectDTTsupernatant',S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(released,residual));g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Sourceultimatebeaddisposalunassigned.'))
    output=cleanup(g,'desalt','SourceMicroBioSpinP6removesDTTfromcollectedsupernatant',released,released.state)
    for branch,current in (('preai',preai),('ai',ai_input),('output',output)):
        stock=g.input('telp_supplies_'+branch,'AggregateactualsourceTELP/MightyAmp-adaptation supplies; stockcount/packaging/internalorder unassembled',K.MOLECULAR)
        mighty=g.product('mighty_'+branch,'ActualcopresentMightyAmp primer-extensioncatalyst',K.CATALYST);oligos=g.product('telp_oligos_'+branch,'ActualsourceTELP primer/adaptor functionalinventory; exactwholebodies/species unreviewed',K.OLIGO);other=g.product('telp_other_'+branch,'OtheractualsourceTELP supplies includingextensionnucleotidefunctions; identities/doses unassembled',K.MOLECULAR)
        g.react('telp_components_'+branch,'ActualTELPsupplycomponents, notnew independentdoses orpremixclaim',S.GENERIC,((stock,R.TEMPLATE),),(mighty,oligos,other),changes_molecule=False)
        installed=g.product('telp_'+branch+'_product','Sourceadaptor-installed/copy-containingpreparation; occupiedstrands/ligatedjunctions/internalorderunassembled',K.MOLECULAR,state=pop('SourceTELP/MightyAmpadaptation withoriginalDNA/remnants retained'))
        g.react('telp_'+branch,'SourceTELPcontinuationwithoneMightyAmpon-beadextensionbeforePCR; internalgeometryunassembled',S.MAGIC,((current,R.TEMPLATE),(mighty,R.CATALYST),(oligos,R.STRUCTURAL),(other,R.STRUCTURAL)),(installed,),reason='Actualsourceadaptation; no predecessorrecipe, stockpackaging or wholeoligobody borrowed.');current=installed
        f=g.input('forward_'+branch,'SourceNEBindexedPCRprimer,wholepolymer/memberunreviewed',K.OLIGO);r=g.input('reverse_'+branch,'SourceNEBpairedPCRprimer,wholepolymerunreviewed',K.OLIGO)
        mix=g.input('pcr_mix_'+branch,'ActualsourceNEBNext2xPCRMasterMix',K.MOLECULAR)
        cat=g.product('enzyme_'+branch,'CopresentactualNEBNextmixpolymerase',K.CATALYST);nt=g.product('nt_'+branch,'CopresentactualNEBNextmixnucleotidecomponent',K.MOLECULAR);buf=g.product('pcr_buffer_'+branch,'CopresentactualNEBNextmixbuffer',K.CONDITION)
        g.react('pcr_components_'+branch,'OneactualsourcePCRmixcomponents, noextrafreshstockadditions',S.GENERIC,((mix,R.TEMPLATE),),(cat,nt,buf),changes_molecule=False)
        lib=g.product('amp_'+branch,'Source'+branch+'PCRDNApreparation;originaltemplates/adductDNA/remnantsnotassumedabsent',K.MOLECULAR,state=pop('PCRcopiedDNA plusoriginaltemplates; modification-dependentcopyclasses unquantified'))
        g.react('pcr_'+branch,'SourceNEBNextindexedPCRafterTELPadaptation',S.PCR,((current,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION)),(lib,),geometry=PCRGeometry.unavailable(Scene(),reason='WholeTELP/adapter/primerassemblyunreviewed; no perfectcopyoractualindexvalueinferred.'))
        qa=g.product('bioanalyzer_'+branch,'Actual Agilent2100Bioanalyzer analyticalallocation; amountunassigned',K.QC);retained=g.product('checked_'+branch,'Source libraryremaining afteranalyticalallocation',K.MOLECULAR)
        g.react('qa_'+branch,'Source Agilent2100Bioanalyzer libraryQC allocation',S.SPLIT,((lib,R.TEMPLATE),),(qa,retained));g.mark_terminal(Terminal(qa.id,Disposition.QC,'Actual sourceanalyticalQA, notwholelibrarytransformation.'))
        loaded=g.product('library_'+branch,'Actual selectedsourcePCR-DNA loadingcontribution; amountunassigned',K.LIBRARY);remaining=g.product('loading_remainder_'+branch,'Untransferred library stock; fateunassigned',K.MOLECULAR)
        g.react('load_'+branch,'Actual selectedHiSeq2500contribution, noentirestockconsumptionassumption',S.SPLIT,((retained,R.TEMPLATE),),(loaded,remaining));g.mark_terminal(Terminal(remaining.id,Disposition.UNRESOLVED,'UntransferredsourceDNAstock.'))
        native(g,loaded,'SourceHiSeq2500;single-end100bporlongerissuggested,exactacquiredcycles/primergeometryunreviewed.')
    g.declare_reviewed_scope(('dna','extract_pk','extract_sds','extract','organic_condition','organic_step','fragment_stock','fragment','preai_split','ai_stock','ai_conditions','ai','ai_cleanup_condition','ai_cleanup_step','ai_input_allocation','ai_control','click_stock','click_conditions','click','click_cleanup_condition','click_cleanup_step','beads','capture','capture_wash_condition','capture_wash_step','alkali_stock','alkali','strand_wash_condition','strand_wash_step','dtt_stock','dtt_conditions','dtt','release_condition','collect','desalt_condition','desalt_step'),SOURCE+'; literal selectedextraction/adduct/capture/release targets only. TELP/fullsourceobligation notcovered.')
    g.add_source_item(SourceItem('original_telp','SourceTELPprotocol23andoriginalkitgeometry;finalVORmethodapplicability','ActualTELPpredecessor/kitversionnotreadincap;MightyAmpadaptationretainedlocally.'))
    return g

def selftest():
    from checks import Check
    c=Check();x=CyclizedDNA('ACGT',frozenset({1}));r=x.advance('ai').advance('biotin').advance('released');g=protocol_graph()
    c('DTTreleasepreservesAIadductandlinkerremnant','5fC-AI-linkerremnant' in r.modifications[1],True)
    c('PhysicalmodifiedtemplateisnotconvertedtoT',r.reference,'ACGT')
    c('SelectedfirstcopyopposesadductwithA',r.selected_first_copy(),'ACAT')
    c('Sourcecomparisonsremainseparateacquisitions',set(g.sequencing),{'library_preai','library_ai','library_output'})
    c('Physicalcollectionfollowssulfidereduction',g.steps['collect'].inputs[0].material,'dtt_product')
    c('NaOHdoesnotitselfdiscardorcleavebiotin',g.steps['alkali'].outputs,('alkali_product',))
    c('SourceAI is physicalmolecularcontributor',g.materials['ai_stock'].kind,K.MOLECULAR)
    c('PCR consumes sameactualmixnucleotides',g.steps['pcr_output'].inputs[4].material,'nt_output')
    c('ActualBioanalyzer analyticalfractionretained',g.terminals['bioanalyzer_output'].disposition,Disposition.QC)
    c('FinalDNAstockremainder notdiscarded',g.terminals['loading_remainder_output'].disposition,Disposition.UNRESOLVED)
    c('Geometrylimitsremainexplicit',bool(g.publication_problems()))
    c('Privatepreviewhasnofinishedmarker','data-protocol-complete=' not in render_page())
    c.report()
