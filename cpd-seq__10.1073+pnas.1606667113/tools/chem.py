"""Protocol-owned source-limited model; no finished-publication marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Scene, Segment, Construct, feature, revcomp, oligo, panel, workflow_from_protocol_graph, workflow_panel
from protocol_model import role, molecule
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, Terminal, Disposition, SourceItem
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

PAGE_FILE = 'cpd-seq.html'
TITLE = 'CPD-seq — lesion-derived 3-prime ends with nonbiotin strand release'
NOTES = '01_cpd-seq.html'
SOURCE = '<a href="https://doi.org/10.1073/pnas.1606667113">Mao et al. (2016)</a>; published PMC4987812 main and inline SI Methods/Figures1,S1.'
CAVEAT = 'A1 reagent branch represented; the source also uses A2/A3 sample codes. Source oligos are independent local declarations. Complete final-library polymer, barcode cycle boundaries, intermediate endoV chemistry and physical Ion run configuration remain unassembled. Background unblocked ends and APE1-related off-target/offset signals are not assumed absent. Lesion-bearing/template strands and debris remain until source-described separation; alkali releases the nonbiotin strand.'
READOUT = 'IonProton source readout; exact read length/run primer and barcode spans unassembled. Infer CPD dinucleotide adjacent to mapped read5-prime end on the opposite strand. Demultiplexing and nondipyrimidine/pUC19 filtering are computational, not physical DNA loss.'
PREVIEWS = (('Selected A1 library preparation', 'library'),)
TRP1_TOP = 'CCTCTCTATGGGCAGTCGGTGATT'
TRP1_BOTTOM = 'ATCACCGACTGCCCATAGAGAGGC'
A1_TOP = 'ATCCTCTTCTGAGTCGGAGACACGCAGGGATGAGATGGC'
A1_BOTTOM = 'CCATCTCATCCCTGCGTGTCTCCGACTCAGAAGAGGATNNNNNN'
PRIMER_A = 'CCATCTCATCCCTGCGTGTCTCCGAC'
PRIMER_TRP1 = 'CCTCTCTATGGGCAGTCGGTGATT'
SAMPLE = feature('cpd_a1_sample', 'sample_index', 'unknown')
SOURCE_OLIGO_EXCEPTIONS = {'trp1_top': 'Source inline SI1289 prints the final T after the English word phosphorothioate, which the general parser does not rejoin. Local source test checks the complete original wording and terminal T; source oligo is not changed.'}
def oligos():
    return [oligo('trP1 top', [Segment('trP1 top',TRP1_TOP)],three='; final linkage phosphorothioate-3′'),
            oligo('trP1 bottom',[Segment('trP1 bottom',TRP1_BOTTOM)],five='5′-phosphate-',three='-dideoxy-3′'),
            oligo('A1 top',[Segment('A1 top',A1_TOP)],five='5′-phosphate-',three='-dideoxy-3′'),
            oligo('A1 bottom',[Segment('A1 bottom',A1_BOTTOM)],five='5′-biotin-',three='-C3 phosphoramidite-3′'),
            oligo('Primer A',[Segment('Primer A',PRIMER_A)]),oligo('Primer trP1',[Segment('Primer trP1',PRIMER_TRP1)])]
def reagent(label,seq,*,mod5='',mod3=''):
    sc=Scene();sc.strand('oligo',[Segment(label,seq)],label=label,mod5=mod5,mod3=mod3);return sc
@dataclass(frozen=True)
class LesionState:
    stage: str
    def __post_init__(self):
        if self.stage not in {'input','a_tail','blocked','endoV','ape1','sap','denatured','ligated','captured','released','copied','library'}:
            raise ValueError('unknown CPD preparation stage')
    def rows(self):
        sc=Scene()
        if self.stage in {'released','copied','library'}:
            if self.stage=='released':
                sc.strand('template',[role('NONbiotin genomic strand with covalently ligated A1-top; exact end assembly unassigned')],label='NaOH-released DNA template')
                sc.strand('other_released',[role('possible other nonbiotin DNA released with the selected class; chemistry and recovery unassigned')],label='other physical released DNA, not assumed absent')
            else:
                sc.strand('library',[role('source adapter/PCR end'),role('genomic copy'),role('A1 sample-code region; exact boundaries unassigned',molecular_feature=SAMPLE),role('opposite source end')],label='primer-derived DNA/library roles')
                sc.strand('templates',[role('residual released DNA templates and other molecules; quantitative removal unassigned')],label='possible template carryover')
            return sc.rows()
        descriptions={'input':'source UV-damaged genomic DNA; actual CPD coordinates/sequence unknown','a_tail':'source end-repaired 3-prime-A-tailed genomic DNA','blocked':'trP1-bearing DNA; pre-existing free3-prime ends TdT/ddA-blocked, incomplete blocking possible','endoV':'T4endoV lesion-cleavage products; intermediate terminal chemical group unassembled','ape1':'source APE1-generated free3-prime-OH immediately upstream of lesion','sap':'new source3-prime-OH with genomic5-prime phosphate removed by SAP','denatured':'separated DNA strands; source new3-prime-OH retained','ligated':'nonbiotin genomic strand covalently joined at new3-prime-OH to A1-top5-prime phosphate','captured':'A1-top/genomic strand associated with biotin-A1-bottom on streptavidin'}
        sc.strand('genomic',[role(descriptions[self.stage])],label='source genomic DNA class')
        if self.stage in {'endoV','ape1','sap','denatured','ligated','captured'}:
            sc.strand('other_dna',[role('other original DNA strands/lesion-derived products; chemical details and recovery unassigned')],label='other physical DNA classes remain')
        if self.stage in {'ligated','captured'}:
            sc.strand('biotin',[role('independent biotin-A1-bottom; noncovalent association to genomic/A1-top')],label='biotin capture strand, not a genomic covalent appendage')
        sc.note('genomic','Stage history only; no complete conversion, global cut prediction, product yield or precise genomic coordinates asserted.')
        return sc.rows()

def protocol_graph():
    g=ProtocolGraph(TITLE,authority='Published PMC4987812 SI Methods1289–1297 and main FigureS1; selected A1 branch')
    dna=g.input('dna','source extracted yeast genomic DNA after UV treatment; repair-time branch upstream',K.MOLECULAR,state=LesionState('input'),primary=True)
    sonic=g.input('sonication','Bioruptor300,15cycles30sON/OFF; near400bp',K.CONDITION)
    cut=g.product('sonicated','fragmented genomic DNA',K.MOLECULAR,state=dna.state)
    g.react('shear','Sonicate genomic DNA',S.FRAGMENTATION,((dna,R.TEMPLATE),(sonic,R.CONDITION)),(cut,))
    clean=cleanup(g,'fragment_clean','Ethanol precipitation and1.2xAMPure selection',cut,cut.state)
    repair=g.input('end_modules','NEBNext End Repair and dA-Tailing Modules; exact internal recipe deferred',K.CATALYST)
    tailed=g.product('tailed','source end-repaired and A-tailed DNA',K.MOLECULAR,state=LesionState('a_tail'))
    g.react('end_processing','Source sequential end repair and dA tailing',S.END_REPAIR,((clean,R.TEMPLATE),(repair,R.CATALYST)),(tailed,))
    top=g.input('trp1_top','source trP1-top; terminal phosphorothioate-T',K.OLIGO,state=reagent('trP1 top',TRP1_TOP,mod3='terminal phosphorothioate linkage'))
    bottom=g.input('trp1_bottom','source5-phosphate/trP1-bottom3-dideoxy oligo',K.OLIGO,state=reagent('trP1 bottom',TRP1_BOTTOM,mod5='phosphate',mod3='dideoxy'))
    lig1=g.input('lig1','NEBNext Quick Ligation,16C16h',K.CATALYST)
    adapted=g.product('trp1_ligated','trP1-bearing source DNA',K.MOLECULAR,state=molecule('adapter-bearing source DNA',(role('source trP1 duplex-bearing genomic DNA; full junction unassembled'),)))
    geometry=reported_adapter_end('source trP1 T-overhang end',kind='cohesive',polarity='3-prime',terminal_sequence='T',body_sequence=TRP1_TOP[:-1],statement='Source explicitly states T overhang and opposing5-prime phosphate. Full duplex/junction assembly is not inferred from this one strand.')
    g.react('ligate_trp1','Ligate source trP1 adapter',S.LIGATION,((tailed,R.TEMPLATE),(top,R.STRUCTURAL),(bottom,R.STRUCTURAL),(lig1,R.CATALYST)),(adapted,),geometry=geometry)
    ready=cleanup(g,'trp1_clean','AMPure cleanup after trP1 ligation',adapted,adapted.state)
    tdt=g.input('tdt','Source TdT enzyme,37C2h; supply packaging unassigned',K.CATALYST)
    ddatp=g.input('ddatp','Source dideoxyATP molecular substrate for end blocking; dose/stock packaging unassigned',K.MOLECULAR)
    blocked=g.product('blocked','pre-existing free3-prime ends blocked with ddA; physical DNA retained',K.MOLECULAR,state=LesionState('blocked'))
    g.react('block_old_ends','TdT/ddATP block pre-existing free3-prime ends',S.TAILING,((ready,R.TEMPLATE),(tdt,R.CATALYST),(ddatp,R.STRUCTURAL)),(blocked,))
    clean=cleanup(g,'block_clean','AMPure cleanup after TdT/ddATP',blocked,blocked.state)
    endov=g.input('endov','30U T4endoV,37C2h',K.CATALYST)
    cleaved=g.product('endov_cut','lesion-associated breaks and other DNA retained',K.MOLECULAR,state=LesionState('endoV'))
    g.react('endov_digest','T4endoV creates lesion-associated single-strand breaks',S.DIGESTION,((clean,R.TEMPLATE),(endov,R.CATALYST)),(cleaved,))
    clean=cleanup(g,'endov_clean','Source purification after T4endoV',cleaved,cleaved.state)
    ape=g.input('ape1','15U APE1,37C1.5h',K.CATALYST)
    hydroxyl=g.product('hydroxyl','new source ligatable3-prime-OH upstream of lesion',K.MOLECULAR,state=LesionState('ape1'))
    g.react('make_3oh','APE1 generates ligatable lesion-associated3-prime-OH',S.DIGESTION,((clean,R.TEMPLATE),(ape,R.CATALYST)),(hydroxyl,))
    clean=cleanup(g,'ape_clean','Source purification after APE1',hydroxyl,hydroxyl.state)
    sap=g.input('sap','5U shrimp alkaline phosphatase,37C1h',K.CATALYST)
    dephos=g.product('dephosphorylated','source genomic5-prime phosphate removed; new3-prime-OH retained',K.MOLECULAR,state=LesionState('sap'))
    g.react('dephosphorylate','Remove genomic5-prime phosphates',S.DEPHOSPHORYLATION,((clean,R.TEMPLATE),(sap,R.CATALYST)),(dephos,))
    stop=g.input('sap_heat','65C15min phosphatase inactivation',K.CONDITION)
    stopped=g.product('sap_stopped','DNA plus inactivated SAP until purification',K.MOLECULAR,state=dephos.state)
    g.react('stop_sap','Heat-inactivate SAP',S.THERMAL_TREATMENT,((dephos,R.TEMPLATE),(stop,R.CONDITION)),(stopped,))
    clean=cleanup(g,'sap_clean','Source purification after SAP heat treatment',stopped,stopped.state)
    heat=g.input('denaturation','95C5min; immediately ice',K.CONDITION)
    separated=g.product('separated','denatured DNA strands retained',K.MOLECULAR,state=LesionState('denatured'))
    g.react('denature','Denature before source A-adapter ligation',S.THERMAL_TREATMENT,((clean,R.TEMPLATE),(heat,R.CONDITION)),(separated,))
    atop=g.input('a1_top','source A1-top5-phosphate/3-dideoxy',K.OLIGO,state=reagent('A1 top',A1_TOP,mod5='phosphate',mod3='dideoxy'))
    abottom=g.input('a1_bottom','source A1-bottom5-biotin/3-C3; N6 is not an asserted UMI',K.OLIGO,state=reagent('A1 bottom',A1_BOTTOM,mod5='biotin',mod3='C3 spacer'))
    lig2=g.input('lig2','240pmol source A duplex; NEBNext Quick Ligation,16C16h',K.CATALYST)
    aligated=g.product('a_ligated','nonbiotin covalent genomic/A1-top plus independent biotin strand',K.MOLECULAR,state=LesionState('ligated'))
    g.react('ligate_a1','Ligate source A1 adapter at new genomic3-prime-OH',S.LIGATION,((separated,R.TEMPLATE),(atop,R.STRUCTURAL),(abottom,R.STRUCTURAL),(lig2,R.CATALYST)),(aligated,),geometry=end_boundary('Source A1 oligos and new3OH/adapter5P known; generic N6 splint and complete genomic junction not aligned to an assigned genome'))
    beads=g.input('beads','source streptavidin beads; detailed workup not assigned',K.SELECTION)
    retained=g.product('retained','biotin-associated DNA population on beads',K.MOLECULAR,state=LesionState('captured'))
    flow=g.product('capture_other','nonretained capture material; exact physical fate unassembled',K.MOLECULAR,state=molecule('other capture fraction',(role('other DNA/adapters and reaction products; partition unquantified'),)))
    g.react('capture','Capture through independent biotin-A1-bottom strand',S.AFFINITY_CAPTURE,((aligated,R.TEMPLATE),(beads,R.SELECTION)),(retained,flow))
    g.mark_terminal(Terminal(flow.id,Disposition.UNRESOLVED,'Not carried into represented release; exact wash/disposal conditions are not reconstructed.'))
    naoh=g.input('naoh','0.15M NaOH release',K.CONDITION)
    released=g.product('released','NONbiotin DNA released for copying',K.MOLECULAR,state=LesionState('released'))
    bead_strand=g.product('bead_strand','biotin adapter strand remains bead-associated',K.MOLECULAR,state=reagent('biotin A1-bottom remaining on beads',A1_BOTTOM,mod5='biotin',mod3='C3 spacer'))
    g.react('release','Alkali separates nonbiotin DNA from biotin capture strand',S.STRAND_SEPARATION,((retained,R.TEMPLATE),(naoh,R.CONDITION)),(released,bead_strand))
    g.mark_terminal(Terminal(bead_strand.id,Disposition.UNRESOLVED,'Bead-associated biotin strand is not the released genomic template; final disposal not reconstructed.'))
    primer_a=g.input('primer_a','source primerA for copying and PCR',K.OLIGO,state=reagent('Primer A',PRIMER_A))
    polymerase=g.input('copy_polymerase','source dsDNA synthesis polymerase; exact formulation deferred',K.CATALYST)
    copied=g.product('copied','source dsDNA synthesis product and residual template molecules',K.MOLECULAR,state=LesionState('copied'))
    g.react('copy','Copy released nonbiotin template using primerA',S.EXTENSION,((released,R.TEMPLATE),(primer_a,R.STRUCTURAL),(polymerase,R.CATALYST)),(copied,))
    clean=cleanup(g,'copy_clean','AMPure purification of dsDNA; source diagnosticPCR is QC',copied,copied.state)
    primer_trp1=g.input('primer_trp1','source primertrP1',K.OLIGO,state=reagent('Primer trP1',PRIMER_TRP1))
    pcr=g.input('pcr_condition','six-cycle standard source PCR; complete formulation deferred',K.CONDITION)
    amplified=g.product('amplified','source6-cycle library product; exact assembly unassigned',K.MOLECULAR,state=LesionState('library'))
    g.react('pcr','Six-cycle libraryPCR',S.PCR,((clean,R.TEMPLATE),(primer_a,R.STRUCTURAL),(primer_trp1,R.STRUCTURAL),(pcr,R.CONDITION)),(amplified,),geometry=PCRGeometry.unavailable(molecule('source dsDNA library',(role('full adapter/strand assembly not yet validated'),)),reason='Printed primers verified locally; full source adapter/template topology and run architecture remain unassembled.'))
    library=cleanup(g,'library','Final AMPure cleanup; selected A1 library (other source sample codes outside this route)',amplified,amplified.state)
    g.mark_terminal(Terminal(library.id,Disposition.UNRESOLVED,'Source IonProton endpoint; final assembly, source sample-code boundaries and run geometry unassembled.'))
    g.declare_reviewed_scope(('dna','sonication','shear','fragment_clean_condition','fragment_clean_step','end_modules','end_processing','trp1_top','trp1_bottom','lig1','ligate_trp1','trp1_clean_condition','trp1_clean_step','tdt','ddatp','block_old_ends','block_clean_condition','block_clean_step','endov','endov_digest','endov_clean_condition','endov_clean_step','ape1','make_3oh','ape_clean_condition','ape_clean_step','sap','dephosphorylate','sap_heat','stop_sap','sap_clean_condition','sap_clean_step','denaturation','denature','a1_top','a1_bottom','lig2','ligate_a1','beads','capture','naoh','release','primer_a','copy_polymerase','copy','copy_clean_condition','copy_clean_step','primer_trp1','pcr_condition','pcr','library_condition','library_step'),g.authority)
    g.add_source_item(SourceItem('companion_protocol','Direct emRiboSeq/EndoSeq protocol version and remaining assembly boundaries','Child Methods1268/ref19 explicitly cites DOI10.1038/nprot.2015.099; original not read within capped pass.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();source=Source(sources_root()/HERE.parent.name/'pnas.1606667113_PMC4987812.html.txt','Published PMC4987812 inlineSI1289')
    if have(c,source):
        text=source.path.read_text()
        c('trP1-top complete source spelling includes final phosphorothioate-T','CCTCTCTATGGGCAGTCGGTGAT-phosphorothioate-T' in text)
        for label,seq in [('trP1 bottom',TRP1_BOTTOM),('A1 top',A1_TOP),('A1 bottom',A1_BOTTOM),('primer A',PRIMER_A),('primer trP1',PRIMER_TRP1)]:
            # The final trP1 primer occurs continuously in the same source paragraph.
            c(label+' direct source transcription',seq in text)
    g=protocol_graph();c('Private source boundary blocks publication',bool(g.publication_problems()))
    c('No finished page marker','data-protocol-complete=' not in render_page())
    c('Alkali release outputs distinct nonbiotin and bead-retained strands',g.steps['release'].outputs,('released','bead_strand'))
    c.report()
