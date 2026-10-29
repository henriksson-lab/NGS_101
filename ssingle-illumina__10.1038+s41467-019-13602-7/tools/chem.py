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

def selftest():
    from checks import Check
    c = Check(); g = protocol_graph()
    c('Pending source/assembly boundary prevents finished publication', bool(g.publication_problems()))
    c('Private preview has no finished marker', 'data-protocol-complete=' not in render_page())
    c.report()

import illumina as il
PAGE_FILE='ssingle-illumina.html'
TITLE='SSiNGLe-ILM — native3-prime-OH tails and chimeric-primer copying'
NOTES='01_ssingle-illumina.html'
SOURCE='<a href="https://doi.org/10.1038/s41467-019-13602-7">Cao et al. (2019)</a>; published PMC6925131 principle77/rationale93 and preparationMethods203–223.'
CAVEAT='Selected K562 native-break Illumina route, including source cell-line SDS permeabilization. Optional SAP/nicking controls, PBMC handling and Helicos/SMS route are excluded. Native accessible3OH can include SSBs,DSBs and repair products; MNase-created3P is a distinct non-tailing class. Free dT50-rT3 primer has terminal RNA residues, whereas an extended product has a DNA3OH. Tail lengths, population yields, full RNA/DNA copying topology and final run/index assembly remain unassigned.'
READOUT='Illumina paired-end150bp. Source Read1 beginsG10 and Read2 beginsT12; call the break at first base after Read2T12. Internal-priming/mapping/untailed-control filtering is computational. The Helicos-specific fill-and-lock offset is not applied here.'
PREVIEWS=(('Selected Illumina preparation','library'),)
P5G10=il.NEBNEXT_UNIVERSAL_PRIMER+'G'*10+'HN'
P7T12=il.P7+'CGTGAT'+il.TRUSEQ_READ2+'T'*12+'VN'
def chimeric_primer():
    sc=Scene();sc.strand('primer',[Segment('source50deoxythymidines','T'*50),role('source three terminal ribothymidines(rT); not DNA-T or RNA-U',length_bp=3)],label='dT50-rT3 primer',mod3='RNA-ended3-primeOH; source says free primer resists TdT tailing');return sc

def primer_parts(name,seq):
    return [Segment(name+' fixed body',seq[:-2]),Segment(name+' degenerate anchor',seq[-2:],placeholder=True)]

def oligos():
    return [oligo('P5G10',primer_parts('P5G10',P5G10)),oligo('P7T12',primer_parts('P7T12',P7T12)),'<p>Chimeric linear primer as written: oligo-d(T)50-r(T)3; three terminal ribothymidines are preserved as RNA chemistry.</p>']
@dataclass(frozen=True)
class EndHistory:
    stage:str
    def __post_init__(self):
        if self.stage not in {'nuclei','fragmented','denatured','polya','blocked_A','linear','polyc','blocked_C','library'}:
            raise ValueError('unsupported SSiNGLe stage')
    def rows(self):
        sc=Scene()
        if self.stage in {'nuclei','fragmented','denatured','polya','blocked_A'}:
            parts=[role('native break-bearing genomic DNA')]
            if self.stage in {'polya','blocked_A'}:parts.append(role('TdT polya tail; length distribution unassigned'))
            sc.strand('native',parts,label='native eligible3OH class',mod3='source ddN cap; identity/efficiency unassigned' if self.stage=='blocked_A' else '3-primeOH')
            if self.stage!='nuclei':sc.strand('mnase',[role('MNase fragment end; exact sequence unassigned')],label='artificial fragmentation-end class',mod3='phosphate; not a TdT3OH substrate in selected source route')
        else:
            if self.stage=='library':
                sc.strand('library',[role('source P5G10-derived end'),role('genomic copied sequence'),role('P7T12-derived end; complete assembly and index meaning unassigned')],label='source PCR library roles')
            else:
                parts=[role('source dT50-rT3 primer-derived5-prime region, including three internal RNA residues'),role('newly copied genomic DNA')]
                if self.stage in {'polyc','blocked_C'}:parts.append(role('source TdT polyc tail; length distribution unassigned'))
                sc.strand('copy',parts,label='extended DNA/RNA-chimeric product',mod3='ddC' if self.stage=='blocked_C' else 'DNA3-primeOH')
            sc.strand('originals',[role('original genomic strands retain polya/ddN caps or MNase3P; no new polyc assigned')],label='original DNA remains a distinct physical class')
            sc.strand('free_primer',[role('possible residual free dT50-rT3 primer; terminal RNA, no polyc assigned')],label='RNA-ended free-primer class; recovery unquantified')
        return sc.rows()

def step(g,i,label,kind,x,reagent,typ,port,state):
    r=g.input(i+'_reagent',reagent,typ);y=g.product(i+'_product',label,K.MOLECULAR,state=state);g.react(i,label,kind,((x,R.TEMPLATE),(r,port)),(y,));return y

def protocol_graph():
    g=ProtocolGraph(TITLE,authority='Published PMC6925131 native-break K562 Illumina branch; preparation203–223,rationale77/93')
    cells=g.input('cells','source K562 cells; native-break branch without optional nuclease/SAP controls',K.MOLECULAR,state=EndHistory('nuclei'),primary=True)
    fixed=step(g,'fix','1%formaldehyde10min',S.CROSS_LINKING,cells,'source formaldehyde fixation',K.CONDITION,R.CONDITION,cells.state)
    stopped=step(g,'quench','Source prints1.375Mglycine5min; no concentration correction invented',S.QUENCHING,fixed,'source glycine quench',K.CONDITION,R.CONDITION,fixed.state)
    nuclei=cleanup(g,'nuclei','Cell collection/PBSwash;PIPES/KCl/NP40lysis;MNase-buffer nuclear washes',stopped,stopped.state)
    permeabilized=step(g,'sds','Cell-line0.3%SDS37C1h then1.8%Triton5min',S.THERMAL_TREATMENT,nuclei,'source cell-line permeabilization conditions',K.CONDITION,R.CONDITION,nuclei.state)
    nuclei=cleanup(g,'sds_clean','Collect nuclei and twoMNase-buffer washes',permeabilized,permeabilized.state)
    digested=step(g,'mnase','MNase/RNAIf30min on ice; source150–500bp preparation',S.DIGESTION,nuclei,'source1200–3000U MNase and200U RNAIf',K.CATALYST,R.CATALYST,EndHistory('fragmented'))
    stopped=step(g,'stop_mnase','EDTA/SDS stop; products remain until DNA recovery',S.QUENCHING,digested,'source EDTA and nuclei lysis buffer',K.CONDITION,R.CONDITION,digested.state)
    protein=step(g,'proteinase','ProteinaseK55C45min then65C45min; reverse crosslinking',S.DEPROTEINIZATION,stopped,'source proteinaseK',K.CATALYST,R.CATALYST,stopped.state)
    dna=cleanup(g,'dna_clean','KCl/isopropanol precipitation;70%ethanol;2xVAHTS recovery',protein,EndHistory('fragmented'))
    denatured=step(g,'denature_a','95C5min then ice before TdT',S.THERMAL_TREATMENT,dna,'source TdTbuffer/CoCl2',K.CONDITION,R.CONDITION,EndHistory('denatured'))
    tailed=step(g,'polya','TdT/dATP37C30min; only eligible native3OH class shown',S.TAILING,denatured,'5U TdT/source dATP for ILM',K.CATALYST,R.CATALYST,EndHistory('polya'))
    capped=step(g,'cap_a','ddNTP37C30min; source ddN member unassigned',S.TAILING,tailed,'source TdT/ddNTP blocking',K.CATALYST,R.CATALYST,EndHistory('blocked_A'))
    stopped=step(g,'stop_a','70C10min TdT inactivation',S.THERMAL_TREATMENT,capped,'70C10min',K.CONDITION,R.CONDITION,capped.state)
    dna=cleanup(g,'tail_clean','2xVAHTS cleanup before linear copying',stopped,stopped.state)
    primer=g.input('linear_primer','source dT50-rT3 chimeric primer',K.OLIGO,state=chimeric_primer())
    taq=g.input('taq_linear','Tiangen Taq;10source primer-extension cycles94/55/72C',K.CATALYST)
    linear=g.product('linear','extended products, old genomic strands and possible residual free primer',K.MOLECULAR,state=EndHistory('linear'))
    g.react('linear_copy','Ten-cycle linear amplification using source chimeric primer',S.EXTENSION,((dna,R.TEMPLATE),(primer,R.STRUCTURAL),(taq,R.CATALYST)),(linear,))
    dna=cleanup(g,'linear_clean','2xVAHTS cleanup after linear copying',linear,linear.state)
    denatured=step(g,'denature_c','95C5min then ice before polyc tailing',S.THERMAL_TREATMENT,dna,'source TdTbuffer/CoCl2',K.CONDITION,R.CONDITION,dna.state)
    ctailed=step(g,'polyc','TdT/dCTP37C30min on new DNA3OH; no polyc on old ddN/3P or RNA-ended free primer',S.TAILING,denatured,'source TdT/dCTP',K.CATALYST,R.CATALYST,EndHistory('polyc'))
    capped=step(g,'cap_c','TdT/ddCTP37C30min',S.TAILING,ctailed,'source TdT/ddCTP',K.CATALYST,R.CATALYST,EndHistory('blocked_C'))
    stopped=step(g,'stop_c','70C10min TdT inactivation',S.THERMAL_TREATMENT,capped,'70C10min',K.CONDITION,R.CONDITION,capped.state)
    dna=cleanup(g,'c_tail_clean','2xVAHTS before source PCR',stopped,stopped.state)
    forward=g.input('p5g10','source P5G10 primer; HN not an assigned UMI',K.OLIGO,state=molecule('P5G10',tuple(primer_parts('P5G10',P5G10))))
    reverse=g.input('p7t12','source P7T12 primer; VN anchoring, no UMI assignment',K.OLIGO,state=molecule('P7T12',tuple(primer_parts('P7T12',P7T12))))
    taq=g.input('taq_pcr','Tiangen Taq; source initial3min94C,55Ccycle,37Cslow-ramp cycle,then30cycles94/72C',K.CATALYST)
    amplified=g.product('amplified','source libraryPCR products plus unassigned residual molecular classes',K.MOLECULAR,state=EndHistory('library'))
    g.react('pcr','Source anchored-primer libraryPCR; full product assembly unresolved',S.PCR,((dna,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(taq,R.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(molecule('source tailed products',(role('full RNA/DNA-template and variable-tail PCR topology unassembled'),)),reason='Source primers are verified; complete primer/template and final product assembly remains unreviewed.'))
    library=cleanup(g,'library','Final2xVAHTS cleanup; source Qubit assessment',amplified,amplified.state)
    g.mark_terminal(Terminal(library.id,Disposition.UNRESOLVED,'Source IlluminaPE150; exact library/run/index topology unassembled.'))
    g.declare_reviewed_scope(('cells','fix_reagent','fix','quench_reagent','quench','nuclei_condition','nuclei_step','sds_reagent','sds','sds_clean_condition','sds_clean_step','mnase_reagent','mnase','stop_mnase_reagent','stop_mnase','proteinase_reagent','proteinase','dna_clean_condition','dna_clean_step','denature_a_reagent','denature_a','polya_reagent','polya','cap_a_reagent','cap_a','stop_a_reagent','stop_a','tail_clean_condition','tail_clean_step','linear_primer','taq_linear','linear_copy','linear_clean_condition','linear_clean_step','denature_c_reagent','denature_c','polyc_reagent','polyc','cap_c_reagent','cap_c','stop_c_reagent','stop_c','c_tail_clean_condition','c_tail_clean_step','p5g10','p7t12','taq_pcr','pcr','library_condition','library_step'),g.authority)
    g.add_source_item(SourceItem('remaining_source_and_assembly','Remaining source SI and complete RNA/DNA/primer/run assembly','Published source SI and main Figure1/Methods; full original coverage and template/run reconstruction not completed.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();source=Source(sources_root()/HERE.parent.name/'s41467-019-13602-7_PMC6925131.xml.txt','Published SSiNGLe Methods219–223')
    if have(c,source):
        text=source.path.read_text();upper=text.upper()
        c('P5G10 transcribed from defining original',P5G10 in upper)
        c('P7T12 transcribed from defining original',P7T12 in upper)
        c('Terminal ribothymidine chemistry retained from source','oligo-d(T)50-r(T)3' in text)
    g=protocol_graph();c('Incomplete assembly blocks finished publication',bool(g.publication_problems()))
    c('Private preview has no finished marker','data-protocol-complete=' not in render_page())
    c.report()
