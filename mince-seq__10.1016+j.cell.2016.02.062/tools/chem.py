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

def cleanup(g, ident, source, template, state, supplies=()):
    reagent = g.input(ident + '_condition', source, K.CONDITION)
    out = g.product(ident, 'recovered source preparation', K.MOLECULAR, state=state)
    waste = g.product(ident + '_waste', 'physically separated workup fractions; exact molecular partition unquantified', K.WASTE)
    g.react(ident + '_step', source, S.CLEANUP, ((template, R.TEMPLATE), (reagent, R.CONDITION))+tuple(supplies), (out, waste))
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
def op(g,key,label,kind,x,out_state,reagent_label,reagent_kind=K.CONDITION,reagent_role=R.CONDITION,extra_inputs=()):
    reagent=g.input(key+'_stock',reagent_label,reagent_kind)
    product=g.product(key+'_product',label,K.MOLECULAR,state=out_state)
    kw={'reason':'Source operation known; exact molecular substrate/product relation and mixture fractions unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(key,label,kind,((x,R.TEMPLATE),(reagent,reagent_role))+tuple(extra_inputs),(product,),**kw)
    return product


def native(g,x,limitation):
    # Source acquisition may be stated without fabricated sequencing primers/cycle counts.
    g.sequence(x,SequencingPlan(limitation=limitation))

def pending(g,label):
    g.add_source_item(SourceItem('unreviewed_assembly',label,'Named source operations are represented; historical end/kit/primer assembly and complete original-version applicability unreviewed.'))



PAGE_FILE='mince-seq.html'
TITLE='MINCE-seq'
NOTES='01_mince-seq.html'
SOURCE='Associated manuscriptNIHMS772603 maininitialResults/Fig1B/Methods1194–1198; DOI10.1016/j.cell.2016.02.062. ChemistrySIgated.'
CAVEAT='SelectedS2 EdUpulse route;1hthymidinechase is separate sourcecondition. CovalentEdU-clickDNA labeling differs from selectedproteinoccupancy/inferred nucleosome maps. Exact publishedSI/linked librarygeometry unreviewed; no universal protection/capture/recovery assumed.'
READOUT='Source input and streptavidin-enrichedDNA sequenced with fragment-capturing library described as25–500bp. Illumina acknowledgment establishes platform, but exact run/cycles/primers unreviewed. TandemH3-MINCEChIP validation not merged into baseline.'
PREVIEWS=(('EdUpulse capture and parallelinput','library_output'),)
def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    x=g.input('cells','SourceS2 cellpreparation; no pooledcelltypes',K.BIOLOGICAL,state=pop('chromatinDNA/RNA/proteins and othercellmaterial'),primary=True)
    x=op(g,'pulse','Source10µMEdU10min; rapid harvest/icePBS',S.MAGIC,x,pop('nascentEdU-bearingDNA and unlabelledDNA/cellmaterial; exactsites and fractions unassigned'),'Actual EdU contribution10uM;1h20uMthymidinechase is a separate condition, not an extra stock here',K.MOLECULAR,R.STRUCTURAL)
    fixation=g.input('fix_programme','Source fixation15min',K.CONDITION)
    x=op(g,'fix','Source1%formaldehyde15min',S.CROSS_LINKING,x,pop('sourcecrosslinkedchromatin withEdU-bearingDNA and othermaterial'),'Actual source formaldehyde contribution1percent',K.MOLECULAR,R.STRUCTURAL,extra_inputs=((fixation,R.CONDITION),))
    quenching=g.input('quench_programme','Source quench5min',K.CONDITION)
    x=op(g,'quench','Source250mMglycine5min',S.QUENCHING,x,x.state,'Actual source glycine contribution250mM',K.MOLECULAR,R.STRUCTURAL,extra_inputs=((quenching,R.CONDITION),))
    x=cleanup(g,'fixed_wash','Source3PBSwashes',x,x.state)
    x=op(g,'permeabilize','Source0.25%Triton/PBS30minRT',S.MAGIC,x,x.state,'Actual source TritonX100 contribution0.25percent',K.MOLECULAR,R.STRUCTURAL)
    bsa=g.input('bsa','Actual BSA protein contribution in source0.5percentBSA/PBS wash; residual carryover unassigned',K.MOLECULAR)
    x=cleanup(g,'bsa_wash','Source0.5%BSA/PBSwash/PBSresuspension',x,x.state,supplies=((bsa,R.STRUCTURAL),))
    tag=g.input('click_tag','Actual source biotinTEGazide first addition100uM; chemical tag, not sequencing identifier',K.MOLECULAR)
    reductant=g.input('click_ascorbate','Actual source sodiumascorbate second addition10mM',K.MOLECULAR)
    copper=g.input('click_copper','Actual source CuSO4 third addition2mM',K.MOLECULAR)
    clicktime=g.input('click_programme','Source1hroomtemperature incubation after ordered additions',K.CONDITION)
    clicked=g.product('click_product','Source clicked preparation and reaction remnants',K.MOLECULAR,state=pop('covalentlyclickedbiotin-EdUDNA and unclicked/unlabelledDNA/proteins/RNA/reagents'))
    g.react('click','Source ordered biotinTEGazide→ascorbate→CuSO4 additions,1hRT; no premix assumed',S.CHEMICAL_CONVERSION,((x,R.TEMPLATE),(tag,R.STRUCTURAL),(reductant,R.STRUCTURAL),(copper,R.STRUCTURAL),(clicktime,R.CONDITION)),(clicked,));x=clicked
    x=cleanup(g,'click_wash','Source3PBSwashes afterclick',x,x.state)
    x=op(g,'lysis','SourceSDS/EDTA/Trisice10min,10fold dilutionbuffer andBransonsonication40stotal30%power2.5son/5soff',S.MAGIC,x,x.state,'sourcecelllysis/dilution/sonication; no invariantfragmentends assumed')
    x=op(g,'mnase','Source37C5minprewarm thenMNase; dose/time titratedto di-/trinucleosomes',S.DIGESTION,x,pop('MNasefragmented protein-associatedbiotinDNA plusunlabelledDNA,cleavageremnants/RNA/proteins'),'sourceMNase',K.CATALYST,R.CATALYST)
    pk=g.input('stop_pk','Actual source proteinaseK20ug/ml contribution before65C12h',K.CATALYST)
    x=op(g,'stop','Source2mMEGTA/0.9%SDS quench plusPK20µg/ml,65C12h crosslinkreversal/proteolysis',S.MAGIC,x,pop('DNA withcovalentclickedEdUbiotin pluspeptides/RNA/remnants; no physicaldebrisdisposal'),'Actual source EGTA/SDS contributions;65C12h programme; packaging unassigned',K.MOLECULAR,R.STRUCTURAL,extra_inputs=((pk,R.CATALYST),))
    x=cleanup(g,'organic_one','FirstsourcePCIextraction/aqueousrecovery',x,x.state)
    x=op(g,'rnase','SourceaqueousphaseRNaseA37C30min;RNAfragmentspersist',S.DIGESTION,x,pop('biotinDNA/unlabelledDNA,RNAcleavagefragments/otherremnants'),'sourceRNaseA',K.CATALYST,R.CATALYST)
    x=cleanup(g,'organic_two','SecondsourcePCIextraction',x,x.state)
    x=cleanup(g,'precipitate','Source70%isopropanol/30mMNaOAc/glycogen−80Catleast30min;spin4C16000rcf30min/75%ethanolwash/0.1xTE',x,x.state)
    inp=g.product('input_aliquot','Source savedinputDNAaliquot beforestreptavidin',K.MOLECULAR,state=x.state);capture=g.product('capture_input','RemainingDNAforcapture',K.MOLECULAR,state=x.state)
    g.react('input_split','Source reservesinputaliquotbeforebiotinpull-down',S.SPLIT,((x,R.TEMPLATE),),(inp,capture))
    beads=g.input('beads','SourceM-280streptavidinDynabeads perhistoricalmanufacturerprotocol',K.SELECTION)
    bound=g.product('bound','SourcecapturedbiotinDNAcomplex pluspossiblebackground',K.MOLECULAR,state=capture.state);other=g.product('capture_other','Nonselectedcapturefraction',K.MOLECULAR,state=pop('unboundDNA/reagents; exactdispositionunreviewed'))
    g.react('capture','SourcebiotinDNApull-down',S.AFFINITY_CAPTURE,((capture,R.TEMPLATE),(beads,R.SELECTION)),(bound,other));g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Source exact flowthroughdispositionunreviewed.'))
    x=cleanup(g,'bound_wash','Sourcestreptavidinbeadwashes',bound,bound.state)
    x=op(g,'release','Sourcebeadsin0.1xTE/1%SDS/PK20µg/ml65C15min; no covalentbiotin-cleavingchemistryinferred',S.DEPROTEINIZATION,x,pop('sourceworkupDNA and bead/proteinpeptide/remnantmixture; exactaffinity dissociationunassembled'),'sourceproteinaseK',K.CATALYST,R.CATALYST)
    x=cleanup(g,'release_organic','SourcePCIextractionafterPK/SDSbeadincubation',x,x.state)
    x=cleanup(g,'release_precip','SourceDNAprecipitation asabove',x,x.state)
    for branch,current in (('input',inp),('output',x)):
        library=g.product('library_'+branch,'Source'+branch+'sequencinglibrary',K.LIBRARY,state=current.state)
        stock=g.input('kit_'+branch,'Source fragment-capturinglibrary25–500bp; exactSI/linkedrecipeunreviewed',K.CONDITION)
        g.react('library_'+branch+'_prep','Source '+branch+' librarypreparationcontinuation',S.MAGIC,((current,R.TEMPLATE),(stock,R.CONDITION)),(library,),reason='Source endpointknown; preciseSI/linkedlibrary chemistry and wholeconstruct unreviewed.')
        native(g,library,'SourceIlluminaNGS ofinputandpulldown; exactrun/primer/cyclegeometryunreviewed.')
    g.add_source_item(SourceItem('gated_si','OriginalSupp_Info_and_tables.docx','https://pmc.ncbi.nlm.nih.gov/articles/instance/4855302/bin/NIHMS772603-supplement-Supp_Info_and_tables.docx; PMCdownloadgate.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('BSA wash consumes actual protein contribution',g.steps['bsa_wash_step'].inputs[2].material,'bsa')
    c('Click consumes three actual ordered additions, not a fictitious premix',tuple(i.material for i in g.steps['click'].inputs[1:4]),('click_tag','click_ascorbate','click_copper'))
    c('Pre-extraction proteolysis consumes actual PK contribution',g.steps['stop'].inputs[2].material,'stop_pk')
    c('Original detailed sequencing source stays uncovered','gated_si' in g.source_coverage,False)
    c('Sourceinputreservedbeforecapture',g.steps['capture'].inputs[0].material,'capture_input')
    c('Inputandpulldownareseparatesequencinglibraries',set(g.sequencing),{'library_input','library_output'})
    c('RNAcleavageprecedessecondphysicalextraction',g.steps['organic_two_step'].inputs[0].material,'rnase_product')
    c('PKdoesnotitselfdisposeofcapturemixture',g.steps['release'].outputs,('release_product',))
    c('Bothlibrarybranchesremaingeometrylimited',bool(g.publication_problems()))
    c('Privatepreviewhasnofinishedmarker','data-protocol-complete=' not in render_page())
    c.report()
