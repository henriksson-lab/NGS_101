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
    g.react(key,label,kind,((x,R.TEMPLATE),(reagent,reagent_role)),(product,),**kw)
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
    dna=g.input('dna','OriginalmESCgenomicDNA, genotypespreparedseparately',K.MOLECULAR,state=pop('sourcegenomicDNAwithnative5fC/othermodifiedcytosines'),primary=True)
    x=op(g,'extract','SourceSDS/proteinaseKgenomicDNApreparation',S.MAGIC,dna,pop('DNA,RNA/proteolytic/cellularremnants'),'sourceSDS/proteinaseK')
    x=cleanup(g,'organic','Sourcephenol/chloroformextractionandethanolprecipitation',x,x.state)
    x=op(g,'fragment','SourceNEBdsDNAFragmentaseapproximately100–400bp',S.FRAGMENTATION,x,x.state,'sourceFragmentase',K.CATALYST,R.CATALYST)
    preai=g.product('preai','SourcepreAIinputcomparisonDNA',K.MOLECULAR,state=x.state);target=g.product('label_input','SourceDNAallocatedtolabel/capture',K.MOLECULAR,state=x.state)
    g.react('preai_split','SourcepreAIinputseparatefromlabelingpreparation',S.SPLIT,((x,R.TEMPLATE),),(preai,target))
    target=op(g,'ai','SourceAI suspension100mMMESpH6,37C24h850rpm',S.CHEMICAL_CONVERSION,target,pop('cyclized5fC-AIadduct-bearingDNA andunlabelled/backgroundmolecules'),'sourceAI5-(2-azidoethyl)-1,3-indandione; templateadductnotDNA-T')
    target=cleanup(g,'ai_cleanup','SourcegenomicDNAQIAquickPCRpurification',target,target.state)
    ai_input=g.product('ai_input','SourceAIinputcomparison; exactsamplingpointrelativeclickunreviewed',K.MOLECULAR,state=pop('AI-inputsourceDNA; exactbiotinconjugationstateunassigned'))
    capture_input=g.product('click_input','SourceAI-labeledDNAforbiotinconjugation',K.MOLECULAR,state=target.state)
    cond=g.input('ai_input_allocation','SourceAIinputversuspulldownpreparations; exactsamplingtimingunreviewed',K.CONDITION)
    g.react('ai_control','SourceAIinputcomparisonallocation',S.MAGIC,((target,R.TEMPLATE),(cond,R.CONDITION)),(ai_input,capture_input),reason='SourceAIcontrolidentityknown; precisephysicalclicksamplingtimingunassembled, no same-stateassertion.')
    x=op(g,'click','SourceDBCO-S-S-PEG3-biotinA112-10,mainprints400mM37C2h; unitnotnormalized',S.CHEMICAL_CONVERSION,capture_input,pop('cyclized5fC-AI-clickedS-SbiotinDNA plusunreactedspecies/reagents'),'sourceDBCOreagent; mainprintedconcentrationfinalversionunreviewed')
    x=cleanup(g,'click_cleanup','SourceQIAquickPCRpurificationafterclick',x,x.state)
    beads=g.input('beads','SourceMyOneStreptavidinC1',K.SELECTION)
    bound=g.product('bound','SourcecapturedDNAduplex/strandcomplexes',K.MOLECULAR,state=x.state);other=g.product('unbound','Companionunboundcapturematerial',K.MOLECULAR,state=pop('unboundDNA/reagents; exactdisposalunreviewed'))
    g.react('capture','Sourcebiotin-DNAcaptureB&WpH7.5plus0.1%Tween20',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),(beads,R.SELECTION)),(bound,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Noquantitativecapture/disposalinferred.'))
    x=cleanup(g,'capture_wash','Sourcecanonicalwash5times then1xSSCwash',bound,bound.state)
    x=op(g,'alkali','Sourcefresh0.15MNaOHRT10min;duplexdenaturationnotcovalentbiotinrelease',S.MAGIC,x,pop('bead-associatedbiotinadductDNAstrands anddissociatednonbiotinDNA/othermaterial'),'sourcealkali')
    x=cleanup(g,'strand_wash','Sourcewash0.1MNaOH thenB&W thenTE; removesdissociatedcompanionstrandswithoutperfectspecificity',x,pop('retainedbiotin-taggedAI-DNAstrandclasswithpossiblebackground; exactendsunassigned'))
    x=op(g,'dtt','Source50mMDTT37C2h,S-Slinkercleavage;5fC-AIadductnoterased',S.CHEMICAL_CONVERSION,x,pop('releasedAI-DNAwithsourcecleavedlinkerremnant plusbeads/tagfragments/reagents'),'sourceDTT')
    condition=g.input('release_condition','SourcesupernatantcollectionafterDTT',K.CONDITION)
    released=g.product('released','SourcecollectedDNAcontainingsupernatant',K.MOLECULAR,state=pop('AI-DNAwithcleavedlinkerremnant andpossiblecarryover'))
    residual=g.product('bead_residual','Sourcepostreleaseresidualbeadmaterial',K.MOLECULAR,state=pop('beads/biotinlinkerfragments/possibleunreleasedDNA'))
    g.react('collect','SourcephysicallycollectDTTsupernatant',S.CLEANUP,((x,R.TEMPLATE),(condition,R.CONDITION)),(released,residual));g.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,'Sourceultimatebeaddisposalunassigned.'))
    output=cleanup(g,'desalt','SourceMicroBioSpinP6removesDTTfromcollectedsupernatant',released,released.state)
    for branch,current in (('preai',preai),('ai',ai_input),('output',output)):
        current=op(g,'telp_'+branch,'SourceTELPcontinuationwithoneMightyAmpon-beadprimerextensionbeforePCR;adaptorligationknown,internalorderunassembled',S.MAGIC,current,pop('sourceadaptor-installed/copy-containinglibrarypreparation; originaltemplates/adductDNA/remnantsretained'),'sourcecitedTELPwithlocalMightyAmpadaptation; no borrowedTELPoligos')
        f=g.input('forward_'+branch,'SourceNEBindexedPCRprimer,wholepolymer/memberunreviewed',K.OLIGO);r=g.input('reverse_'+branch,'SourceNEBpairedPCRprimer,wholepolymerunreviewed',K.OLIGO);cat=g.input('enzyme_'+branch,'SourceNEBNext2xPCRMasterMix',K.CATALYST)
        lib=g.product('library_'+branch,'Source'+branch+'PCRlibrary;originalmolecules/remnantsnotassumedabsent',K.LIBRARY,state=pop('PCRcopiedDNA; modification-dependentcopyclasses unquantified, plus'+('originalpreAIgenomicDNA' if branch=='preai' else 'originalAIinputDNAwithunassignedtagstate' if branch=='ai' else 'originalreleasedAIadductDNA')+'/templates/remnants'))
        g.react('pcr_'+branch,'SourceNEBNextindexedPCRafterTELPadaptation',S.PCR,((current,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(lib,),geometry=PCRGeometry.unavailable(Scene(),reason='WholeTELP/adapter/primerassemblyunreviewed; no perfectcopyoractualindexvalueinferred.'))
        native(g,lib,'SourceHiSeq2500;single-end100bporlongerissuggested,exactacquiredcycles/primergeometryunreviewed.')
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
    c('Geometrylimitsremainexplicit',bool(g.publication_problems()))
    c('Privatepreviewhasnofinishedmarker','data-protocol-complete=' not in render_page())
    c.report()
