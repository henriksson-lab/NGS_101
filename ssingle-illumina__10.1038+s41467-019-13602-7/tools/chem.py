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
PREVIEWS=(('Selected Illumina preparation','library_loading'),)
P5G10=il.NEBNEXT_UNIVERSAL_PRIMER+'G'*10+'HN'
P7T12=il.P7+'CGTGAT'+il.TRUSEQ_READ2+'T'*12+'VN'
SOURCE_OLIGO_EXCEPTIONS={'linear_primer':'Defining source specifies oligo-d(T)50-r(T)3 in prose, not a literal53T sequence. Model retains53physical thymine nucleotides with50deoxy/3ribo sugars; exactsourceprose, backbonewidth and sugarblocks checked locally. Narrow quantified mixed-sugar transcription exception, not replacement ofrTbyDNA-TorRNA-U.'}
@dataclass(frozen=True)
class ChimericPrimerState:
    deoxy_count:int=50
    ribo_count:int=3
    base:str='T'
    def __post_init__(self):
        if (self.deoxy_count,self.ribo_count,self.base)!=(50,3,'T'):raise ValueError('sourcechimericprimer requires50deoxyT plus3riboT')
    @property
    def backbone(self):return self.base*(self.deoxy_count+self.ribo_count)
    @property
    def sugar_blocks(self):return((0,50,'2-deoxyribose'),(50,53,'ribose'))

def chimeric_primer():
    p=ChimericPrimerState();sc=Scene();sc.strand('primer',[Segment('source50deoxythymidines',p.backbone[:50]),Segment('source3RIBOthymidines',p.backbone[50:],note='Thymine bases on RIBOSE, notDNA-T orRNA-U')],label='Source dT50-rT3: two distinct sugar tracts',mod3='RNA3OH; source free-primer TdT resistance is substrate-class behavior, not a missing3OH');return sc

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

def step(g,i,label,kind,x,reagent,typ,port,state,extra=()):
    r=g.input(i+'_reagent',reagent,typ);y=g.product(i+'_product',label,K.MOLECULAR,state=state);g.react(i,label,kind,((x,R.TEMPLATE),(r,port))+tuple(extra),(y,));return y

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
    dna_qa=g.product('dna_qa','Actual source fragmentedDNA agarose/Qubit analytical allocations; amountsunassigned',K.QC)
    selected=g.product('tailing_input','Actual100ng selected fragmentedDNA for tailing',K.MOLECULAR,state=dna.state)
    source_stock=g.product('dna_stock','Untransferred source DNA stock; amount/fateunassigned',K.MOLECULAR,state=dna.state)
    g.react('dna_allocate','ActualDNAQA/100ngtailing/remainingstock allocations',S.SPLIT,((dna,R.TEMPLATE),),(dna_qa,selected,source_stock));g.mark_terminal(Terminal(dna_qa.id,Disposition.QC,'Source actual1%gel/Qubit assessments; no productionreturn'));g.mark_terminal(Terminal(source_stock.id,Disposition.UNRESOLVED,'Untransferred DNA stock'));dna=selected
    denatured=step(g,'denature_a','95C5min then ice before TdT',S.THERMAL_TREATMENT,dna,'source TdTbuffer/CoCl2',K.CONDITION,R.CONDITION,EndHistory('denatured'))
    tailed=step(g,'polya','TdT/dATP37C30min; only eligible native3OH class shown',S.TAILING,denatured,'Actual5U TdT forA-tailing',K.CATALYST,R.CATALYST,EndHistory('polya'),extra=((g.input('datp','Actual2µl10mMdATP forILM',K.MOLECULAR),R.STRUCTURAL),))
    ddnt=g.input('ddntp','Actual2µl10mMddNTP; sourcechoiceamongfourmembersunassigned',K.MOLECULAR)
    capped=g.product('cap_a_product','ActualA-tail blocking with carriedTdT; residualuncappedfractionunassigned',K.MOLECULAR,state=EndHistory('blocked_A'))
    g.react('cap_a','ActualddNTPaddition37°C30min; carriedA-tailTdT, nofreshenzyme',S.TAILING,((tailed,R.TEMPLATE),(g.materials['polya_reagent'],R.CATALYST),(ddnt,R.STRUCTURAL)),(capped,))
    stopped=step(g,'stop_a','70C10min TdT inactivation',S.THERMAL_TREATMENT,capped,'70C10min',K.CONDITION,R.CONDITION,capped.state)
    dna=cleanup(g,'tail_clean','2xVAHTS cleanup; elute16.2µl ALLusedforlinearcopy',stopped,stopped.state)
    primer=g.input('linear_primer','source dT50-rT3 chimeric primer',K.OLIGO,state=chimeric_primer(),chemical_state=ChimericPrimerState())
    taq=g.input('taq_linear','Tiangen Taq;10source primer-extension cycles94/55/72C',K.CATALYST)
    linear=g.product('linear','extended products, old genomic strands and possible residual free primer',K.MOLECULAR,state=EndHistory('linear'))
    g.react('linear_copy','Ten-cycle linear amplification using source chimeric primer',S.EXTENSION,((dna,R.TEMPLATE),(primer,R.STRUCTURAL),(taq,R.CATALYST),(g.input('linear_nts','Actual0.4µl10mMdNTP',K.MOLECULAR),R.STRUCTURAL),(g.input('linear_buffer','Actual2µl10xPCRbuffer',K.CONDITION),R.CONDITION)),(linear,))
    dna=cleanup(g,'linear_clean','2xVAHTS cleanup afterlinearcopy; actual15µleluate',linear,linear.state)
    denatured=step(g,'denature_c','95C5min then ice before polyc tailing',S.THERMAL_TREATMENT,dna,'source TdTbuffer/CoCl2',K.CONDITION,R.CONDITION,dna.state)
    ctailed=step(g,'polyc','TdT/dCTP37C30min on new DNA3OH; no polyc on old ddN/3P or RNA-ended free primer',S.TAILING,denatured,'Actualfresh5UTdT forC-tailing',K.CATALYST,R.CATALYST,EndHistory('polyc'),extra=((g.input('dctp','Actual2µl10mMdCTP',K.MOLECULAR),R.STRUCTURAL),))
    ddc=g.input('ddctp','Actual2µl10mMddCTP capaddition',K.MOLECULAR)
    capped=g.product('cap_c_product','ActualC-tail cap with carriedC-tailTdT, nofreshenzyme',K.MOLECULAR,state=EndHistory('blocked_C'))
    g.react('cap_c','ActualddCTPaddition37°C30min; carriedC-tailTdT',S.TAILING,((ctailed,R.TEMPLATE),(g.materials['polyc_reagent'],R.CATALYST),(ddc,R.STRUCTURAL)),(capped,))
    stopped=step(g,'stop_c','70C10min TdT inactivation',S.THERMAL_TREATMENT,capped,'70C10min',K.CONDITION,R.CONDITION,capped.state)
    dna=cleanup(g,'c_tail_clean','2xVAHTS; actual15.2µleluate ALLusedin20µlPCR',stopped,stopped.state)
    forward=g.input('p5g10','source P5G10 primer; HN not an assigned UMI',K.OLIGO,state=molecule('P5G10',tuple(primer_parts('P5G10',P5G10))))
    reverse=g.input('p7t12','source P7T12 primer; VN anchoring, no UMI assignment',K.OLIGO,state=molecule('P7T12',tuple(primer_parts('P7T12',P7T12))))
    taq=g.input('taq_pcr','Tiangen Taq; source initial3min94C,55Ccycle,37Cslow-ramp cycle,then30cycles94/72C',K.CATALYST)
    amplified=g.product('amplified','source libraryPCR products plus unassigned residual molecular classes',K.MOLECULAR,state=EndHistory('library'))
    g.react('pcr','Source anchored-primer libraryPCR; full product assembly unresolved',S.PCR,((dna,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(taq,R.CATALYST),(g.input('pcr_nts','Actual0.4µl10mMdNTP',K.MOLECULAR),R.STRUCTURAL),(g.input('pcr_buffer','Actual1xPCRbuffer',K.CONDITION),R.CONDITION)),(amplified,),geometry=PCRGeometry.unavailable(molecule('source tailed products',(role('full RNA/DNA-template and variable-tail PCR topology unassembled'),)),reason='Source primers are verified; complete primer/template and final product assembly remains unreviewed.'))
    library=cleanup(g,'library','Final2xVAHTS cleanup; source Qubit assessment',amplified,amplified.state)
    qa=g.product('library_qa','Actual QubitdsDNAHS analytical allocation from20µlfinaleluate; amountunassigned',K.QC)
    stock=g.product('library_stock','Actual source library stock afterQAallocation',K.LIBRARY)
    g.react('library_qa_allocate','Actualfinal20µl libraryeluate analytical/stock allocation',S.SPLIT,((library,R.TEMPLATE),),(qa,stock));g.mark_terminal(Terminal(qa.id,Disposition.QC,'Actual source Qubit assessment, not an extra inventedlibraryaliquot'))
    loading=g.product('library_loading','Actual IlluminaPE150 loadingcontribution; amountunassigned',K.LIBRARY)
    rest=g.product('untransferred_library','Untransferred source librarystock; amount/fateunassigned',K.LIBRARY)
    g.react('library_load_allocate','Actualsource loading allocation',S.SPLIT,((stock,R.TEMPLATE),),(loading,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Actualsource library remainder'))
    g.mark_terminal(Terminal(loading.id,Disposition.UNRESOLVED,'Actualsource DNAloading; broadremainingSI/mixedRNA-DNA/primer/runassemblySourceItem remains unmet, notwaived.'))
    g.declare_reviewed_scope(('cells','fix_reagent','fix','quench_reagent','quench','nuclei_condition','nuclei_step','sds_reagent','sds','sds_clean_condition','sds_clean_step','mnase_reagent','mnase','stop_mnase_reagent','stop_mnase','proteinase_reagent','proteinase','dna_clean_condition','dna_clean_step','denature_a_reagent','denature_a','polya_reagent','polya','cap_a','stop_a_reagent','stop_a','tail_clean_condition','tail_clean_step','linear_primer','taq_linear','linear_copy','linear_clean_condition','linear_clean_step','denature_c_reagent','denature_c','polyc_reagent','polyc','cap_c','stop_c_reagent','stop_c','c_tail_clean_condition','c_tail_clean_step','p5g10','p7t12','taq_pcr','pcr','library_condition','library_step'),g.authority)
    g.declare_reviewed_scope(('dna_allocate','dna_qa','tailing_input','dna_stock','datp','ddntp','linear_nts','linear_buffer','dctp','ddctp','pcr_nts','pcr_buffer','library_qa_allocate','library_qa','library_stock','library_load_allocate','library_loading','untransferred_library'),g.authority)
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
    mixed=g.materials['linear_primer'].chemical_state
    c('Source physicalT backbone and distinct sugar tracts',(len(mixed.backbone),mixed.sugar_blocks),(53,((0,50,'2-deoxyribose'),(50,53,'ribose'))))
    c('A-cap carries the actual originalTdT',g.steps['cap_a'].inputs[1].material,'polya_reagent')
    c('C-cap carries the actual C-tailTdT',g.steps['cap_c'].inputs[1].material,'polyc_reagent')
    c('Actual100ngsourceallocation precedes denaturation',g.steps['denature_a'].inputs[0].material,'tailing_input')
    c('Actual libraryQAfraction separate',g.terminals['library_qa'].disposition,Disposition.QC)
    c.report(exit_on_fail=True)
