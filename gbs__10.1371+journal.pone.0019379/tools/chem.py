"""Original GBS: unphosphorylated adapters and nicked precursor cohorts."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
import illumina as il
from chemdraw import Scene, Segment, Construct, revcomp, feature, panel, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, SourceItem, Terminal, Disposition
from protocol_model import molecule, role
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from page import head, info
PAGE_FILE = 'gbs.html'
SOURCE = 'Elshire et al.2011 DOI10.1371/journal.pone.0019379;published PMC3087801 Methods83/96–109 and Discussion195–203;published TextS1/TableS1.'
BARCODE = 'CTCC'  # Published TableS1 first member;no experimental sample/well assignment.
COMMON = 'C' + il.PE_READ2_ARM
FWD = il.NEBNEXT_UNIVERSAL_PRIMER
REV = il.P7 + 'CGGT' + COMMON
INDEX = feature('gbs.sample', 'sample_index', 'fixed', note='Published TableS1 CTCC selected as a component example;no study well assignment.')
_BARCODE_SOURCE_EQUATION = 'Published main83 defines canonical Read1 body+barcode and CWG+reverse-complement barcode+reverse-complement body;TableS1 explicitly contains selected CTCC member. Full assembled example is an equation across originals,not a contiguous source oligo;generic strict matcher does not instantiate lowercase xxxx/yyyy. Canonical body/complement,4nt selected code and5primeOH remain explicit. Local original-dependent checks and independent whole-component/equation review completed;no synthetic source twin or absence claim.'
SOURCE_OLIGO_EXCEPTIONS = {'barcode_top': _BARCODE_SOURCE_EQUATION, 'barcode_bottom': _BARCODE_SOURCE_EQUATION, 'common_top': 'Published main83 prints CWG+the complete fixed common-adapter body. The drawn CAG strand is one explicitly labelled member of source W=A/T,not a literal fully resolved source string or unique study molecule. Strict sourceguard correctly does not infer that resolution;component equation/allowed member independently checked against the originals. Common bottom and both PCRprimers remain directly checked;no inferred ratio,synthetic source text or absence claim.'}

@dataclass(frozen=True)
class NickedFragment:
    """Source-selected ligation class,not a prediction of all molecules or yield."""
    ends: tuple[str,str]
    def __post_init__(self):
        if len(self.ends)!=2 or any(x not in {'barcode','common'} for x in self.ends):
            raise ValueError('Two source adapter identities must be explicit')
    @property
    def unsealed_adapter_5oh_junctions(self):return tuple(range(2))
    def selected_copy(self):return CopiedFragment(self.ends,self)

@dataclass(frozen=True)
class CopiedFragment:
    """Selected new continuous PCR descendant;original nicked template persists."""
    ends: tuple[str,str]
    parent: NickedFragment
    def __post_init__(self):
        if self.ends!=self.parent.ends:raise ValueError('Copying cannot replace source adapter identities')
    @property
    def source_readout_class(self):
        return 'heterotypic source sequencing class' if self.ends[0]!=self.ends[1] else 'homotypic source limited-readout class'

def cohort_classes():
    return tuple(NickedFragment(x) for x in (('barcode','common'),('common','barcode'),('barcode','barcode'),('common','common')))

def barcode_adapter(overhang='CAG'):
    if overhang not in {'CAG','CTG'}:raise ValueError('Choose an explicit member of source CWG mixture')
    sc=Scene();sc.strand('top',[Segment('source Read1 body',il.TRUSEQ_READ1),Segment('source TableS1 sample index',BARCODE,feature=INDEX)],mod5='OH',label='Barcode adapter;unmodified5primeOH')
    sc.anneal('bottom',[Segment('ApeKI5prime cohesive end',overhang),Segment('source barcode complement',revcomp(BARCODE),feature=INDEX),Segment('Read1 complementary body',revcomp(il.TRUSEQ_READ1))],to='top',pair=('Read1 complementary body','source Read1 body'),mod5='OH',label='Unmodified bottom;selected '+overhang+' member of source CWG mixture')
    sc.ligation_end('bottom','ApeKI5prime cohesive end',kind='cohesive',polarity='5-prime',owner='adapter');return sc

def common_adapter(overhang='CAG'):
    if overhang not in {'CAG','CTG'}:raise ValueError('Choose an explicit member of source CWG mixture')
    sc=Scene();sc.strand('top',[Segment('ApeKI5prime cohesive end',overhang),Segment('common adapter complementary body',revcomp(COMMON))],mod5='OH',label='Common adapter;unmodified5primeOH')
    sc.anneal('bottom',[Segment('source common adapter body',COMMON)],to='top',pair=('source common adapter body','common adapter complementary body'),mod5='OH',label='Unmodified complementary strand')
    sc.ligation_end('top','ApeKI5prime cohesive end',kind='cohesive',polarity='5-prime',owner='adapter');return sc

def ordering(parts):
    sc=Scene();sc.strand('oligo',parts,mod5='OH',label='Source unmodified ordering oligo');return sc

def protocol_graph():
    g=ProtocolGraph('GBS — original ApeKI single-well route',authority=SOURCE,representation='published_instance')
    def op(key,label,kind,material,reagent_label,condition=False,additional=()):
        reagent=g.input(key+'_reagent',reagent_label,K.CONDITION if condition else K.CATALYST)
        out=g.product(key,label+';source mixture,carryover and yield unquantified',K.MOLECULAR)
        g.react(key+'_step',label,kind,((material,R.TEMPLATE),(reagent,R.CONDITION if condition else R.CATALYST))+tuple(additional),(out,));return out
    def recover(key,material,label):
        reagent=g.input(key+'_procedure',label,K.SELECTION);out=g.product(key,'Source recovered preparation;purity unquantified',K.MOLECULAR);waste=g.product(key+'_waste','Actual excluded workup fractions;composition unassigned',K.WASTE)
        g.react(key+'_step',label,S.CLEANUP,((material,R.TEMPLATE),(reagent,R.SELECTION)),(out,waste));g.discard(waste,'Source physical separation,not perfect depletion.');return out
    def allocate(key,material,label,disposition=Disposition.UNRESOLVED,kind=K.MOLECULAR):
        kept=g.product(key,'Selected source aliquot;exact polymer composition unassembled',K.MOLECULAR);rest=g.product(key+'_remainder',label,kind)
        g.react(key+'_split','Source physical aliquot and independent companion',S.SPLIT,((material,R.TEMPLATE),),(kept,rest));g.mark_terminal(Terminal(rest.id,disposition,label));return kept
    def qc(key,material,label):
        retained=g.product(key,'Retained29ul source stock after1ul assay from30ul eluate; exact loading fraction unreviewed',K.MOLECULAR,state=material.state);assay=g.product(key+'_assay',label,K.QC)
        g.react(key+'_allocate','Analytical aliquot separate from production',S.SPLIT,((material,R.TEMPLATE),),(retained,assay));g.mark_terminal(Terminal(assay.id,Disposition.QC,label));return retained
    dna=g.input('dna','Source genomic DNA100ng;upstream CTAB extraction outside selected library graph',K.MOLECULAR,primary=True)
    # The two stocks are annealed and quantified independently,then mixed.
    btop=g.input('barcode_top','Unmodified barcode top;TableS1 CTCC component',K.OLIGO,state=ordering([Segment('source Read1 body',il.TRUSEQ_READ1),Segment('source sample index',BARCODE,feature=INDEX)]))
    bbottom=g.input('barcode_bottom','Unmodified barcode bottom;selected CAG member of CWG mixture',K.OLIGO,state=ordering([Segment('source cohesive end','CAG'),Segment('sample index complement',revcomp(BARCODE),feature=INDEX),Segment('source Read1 complementary body',revcomp(il.TRUSEQ_READ1))]))
    ctop=g.input('common_top','Unmodified common top;selected CAG member of CWG mixture',K.OLIGO,state=ordering([Segment('source cohesive end','CAG'),Segment('source common complementary body',revcomp(COMMON))]))
    cbottom=g.input('common_bottom','Unmodified common bottom',K.OLIGO,state=ordering([Segment('source common body',COMMON)]))
    condition=g.input('anneal_condition','95°C2min,0.1°C/s cooling,25°C30min,4°C hold',K.CONDITION)
    b=g.product('barcode_duplex','Source barcode-adapter stock;shown component does not assign pool stoichiometry',K.MOLECULAR,state=barcode_adapter())
    c=g.product('common_duplex','Source common-adapter stock;shown component does not assign pool stoichiometry',K.MOLECULAR,state=common_adapter())
    g.react('anneal_barcode','Anneal barcode top/bottom stock separately',S.HYBRIDIZATION,((btop,R.TEMPLATE),(bbottom,R.STRUCTURAL),(condition,R.CONDITION)),(b,))
    g.react('anneal_common','Anneal common top/bottom stock separately',S.HYBRIDIZATION,((ctop,R.TEMPLATE),(cbottom,R.STRUCTURAL),(condition,R.CONDITION)),(c,))
    b=qc('barcode_quant',b,'Source PicoGreen assay aliquot;not returned to barcode stock')
    c=qc('common_quant',c,'Source PicoGreen assay aliquot;not returned to common stock')
    stock=g.product('adapter_stock','Source1:1mixed barcode/common stocks;not a predicted ligation product ratio',K.MOLECULAR)
    g.react('adapter_mix','Source mixes separately quantified adapter stocks1:1',S.POOL,((b,R.TEMPLATE),(c,R.TEMPLATE)),(stock,))
    pre=g.product('preplate','Source DNA plus1:1barcode/common adapters dried in the same well before digestion',K.MOLECULAR)
    g.react('preplate_step','Source adapter/DNA plate setup;no invented simultaneous digestion/ligation',S.GENERIC,((dna,R.TEMPLATE),(stock,R.STRUCTURAL)),(pre,),changes_molecule=False)
    digestion_buffer=g.input('digest_buffer','Actual source1xNEBBuffer3 in20ul digestion',K.MOLECULAR)
    cut=op('digest','ApeKI75°C2h with both unmodified adapters already present',S.DIGESTION,pre,'Actual source3.6U ApeKI',additional=((digestion_buffer,R.STRUCTURAL),))
    ligation_supply=g.input('ligation_supply','Actual30ul source solution containing1.66xligasebufferwithATP and640cohesiveendU T4ligase',K.MOLECULAR)
    ligase=g.product('ligase','T4ligase component of sameactual30ul solution',K.CATALYST)
    atp=g.product('ligation_atp','ATP component of same source ligase buffer; exact concentration unassigned',K.MOLECULAR)
    ligation_buffer=g.product('ligation_buffer','Other copresent ligation-buffer constituents; exact formulation unreviewed',K.MOLECULAR)
    g.react('ligation_components','Expose copresent components of actual added ligation solution; no extra doses',S.GENERIC,((ligation_supply,R.STRUCTURAL),),(ligase,atp,ligation_buffer),changes_molecule=False)
    nicked=g.product('nicked','Source mixture:heterotypic and both homotypic partially ligated classes,free adapters/ligase/unligated DNA retained',K.MOLECULAR)
    geometry=AdapterLigationGeometry.known('cohesive',barcode_adapter(),common_adapter(),molecule('Genomic fragment boundary',(role('Source ApeKI target sequence/end/junction geometry unassembled'),)),statement='Source5primeCWG cohesive ends;all adapter5prime ends areOH. Selected ligation joins genomic5P to adapter3OH and leaves opposite adapter5OH junction unsealed;this is not a fully sealed duplex product.')
    # One physical mixture supplies substrate and already-present structural adapters.
    # Reusing its identifier in two roles does not introduce a second copy or aliquot.
    g.react('ligate','Source single-well ligation leaves adapter5OH nicks;no perfect product distribution',S.LIGATION,((cut,R.TEMPLATE),(cut,R.STRUCTURAL),(ligase,R.CATALYST),(atp,R.STRUCTURAL),(ligation_buffer,R.STRUCTURAL)),(nicked,),geometry=geometry)
    nicked=op('inactivate','Source65°C30min ligase inactivation;not physical removal',S.THERMAL_TREATMENT,nicked,'Source heat treatment',condition=True)
    selected=allocate('pool_aliquot',nicked,'Source remainder of each50µl well after5µl pooled aliquot;later fate unassigned')
    other=g.input('other_samples','Other independently prepared sample aliquots;48or96plex source sets',K.MOLECULAR)
    pooled=g.product('pool','Source pooled selected aliquots;all adapter-end classes remain',K.MOLECULAR)
    g.react('pool_samples','Pool source5µl aliquots after ligase inactivation',S.POOL,((selected,R.TEMPLATE),(other,R.TEMPLATE)),(pooled,))
    pooled=recover('pool_cleanup',pooled,'Source main105QIAquickPCRpurification,50µl elution;Fig2 says size-exclusion column,scope unharmonized')
    selected=allocate('pcr_aliquot',pooled,'Source remainder of50µl eluate after2µl PCR allocation;later fate unassigned')
    f=g.input('forward','Actual25pmol source PCRprimerA',K.OLIGO,state=Construct([Segment('source PCRA',FWD)]));r=g.input('reverse','Actual25pmol source PCRprimerB',K.OLIGO,state=Construct([Segment('source PCRB',REV)]))
    mix=g.input('taq_mix','Actual1x Taq MasterMix in source50ulPCR; stockvolume/formulation unreviewed',K.MOLECULAR)
    taq=g.product('taq','Polymerase component of sameactualsource master mix',K.CATALYST)
    nts=g.product('pcr_nts','Copresent source master-mix dNTPs; precise recipe/dose unassembled',K.MOLECULAR)
    other=g.product('pcr_other','Other source master-mix components; formulation unassigned',K.MOLECULAR)
    g.react('pcr_components','Expose copresent functions of oneactualsource master mix',S.GENERIC,((mix,R.STRUCTURAL),),(taq,nts,other),changes_molecule=False)
    copied=g.product('pcr','Selected new PCR descendants plus original nicked DNA,displaced/unligated adapter strands and reaction remnants;all end classes retained',K.MOLECULAR)
    g.react('amplify','Source initial72°C5min fill/copy then18PCRcycles;polymerase copying does not ligase-seal the original5OH nick',S.PCR,((selected,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(taq,R.CATALYST),(nts,R.STRUCTURAL),(other,R.STRUCTURAL)),(copied,),geometry=PCRGeometry.unavailable(molecule('PCR boundary',(role('Full nicked genomic templates,displacement and descendant duplex geometry unassembled'),)),reason='Known selected fill/copy relation does not validate full heterogeneous molecular products or yields.'))
    copied=recover('pcr_cleanup',copied,'Source PCRcleanup,30µl final elution;no fragment-size excision')
    copied=qc('experion',copied,'Source1µl Experion assay;dimer-size assessment is not molecular size selection')
    loaded=g.product('loading','Source selected library contribution from post-QA stock; loadingvolume unassigned',K.MOLECULAR)
    rest=g.product('loading_remainder','Untransferred library stock; presence/volume/fateunassigned',K.MOLECULAR)
    g.react('loading_allocate','Select actualsequencing contribution; no whole29ul loading assertion',S.SPLIT,((copied,R.TEMPLATE),),(loaded,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'No manufacturer loading volume or residual amount inferred.'))
    lib=g.product('library','Source accepted GBS library;heterotypic/homotypic and residual species not perfectly partitioned',K.LIBRARY)
    g.react('assessment','Source acceptance/reconstruction decision if excessive dimers;not physical fragment removal',S.GENERIC,((loaded,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(acquired_roles=('Read 1',),read_lengths={'Read 1':86},limitation='Source GenomeAnalyzerII single-end86nt,inline4–8nt sample tags,no separate index read. Homotypic limited-readout classes are not transformed into heterotypic molecules;actual flowcell behavior and full primer landing unassembled. Computational64nt truncation/polyA padding is not physical DNA processing.'))
    g.declare_reviewed_scope(('digest_reagent','digest_buffer','digest_step','ligation_supply','ligation_components','ligate','inactivate_reagent','inactivate_step','pool_aliquot_split','other_samples','pool_samples','pool_cleanup_procedure','pool_cleanup_step','pcr_aliquot_split','forward','reverse','taq_mix','pcr_components','amplify','pcr_cleanup_procedure','pcr_cleanup_step','experion_allocate','loading_allocate','assessment'),SOURCE+' MainMethods105–109/Discussion197; literal inspected operational subset, not complete source geometry or readout coverage.')
    g.add_source_item(SourceItem('readout','Historical flowcell/Solexa and complete nicked-library/read-primer geometry','Source discussion197–203 distinguishes original nicked inputs,selected PCR descendants and homotypic flowcell limitations;exact molecular recovery/cluster growth not predicted.'))
    return g

def render_page():
    return '\n'.join([head('GBS'),'<main class="wrap"><h1>GBS</h1>','<p><a href="01_gbs.html">Research notes</a></p>',info(SOURCE),info('Private incomplete model. Both adapter types are unphosphorylated;ligation leaves one unsealed junction per end. Selected PCR copies are distinct from original nicked templates. All three end classes persist;no inferred fragment-size selection or conversion to a Y adapter.'),panel(barcode_adapter().rows(),caption='Published CTCC barcode-adapter component;no assigned study sample'),panel(common_adapter().rows(),caption='Source common-adapter component'),workflow_panel(workflow_from_protocol_graph(protocol_graph(),'library')),'</main>'])

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();classes=cohort_classes();copy=classes[0].selected_copy();g=protocol_graph()
    c('All source adapter-end combinations retained',{x.ends for x in classes},{('barcode','common'),('common','barcode'),('barcode','barcode'),('common','common')})
    c('Selected copying preserves adapter identities',copy.ends,copy.parent.ends)
    c('Original template retains both unsealed5OH junctions',copy.parent.unsealed_adapter_5oh_junctions,(0,1))
    c('Homotypic copies remain limited-readout classes',NickedFragment(('barcode','barcode')).selected_copy().source_readout_class,'homotypic source limited-readout class')
    c('Source adapters have no5prime phosphate',tuple(s.mod5 for sc in (barcode_adapter(),common_adapter()) for s in sc.strands.values()),('OH',)*4)
    c('Ligation consumes the same already-adapter-bearing mixture',{x.role for x in g.steps['ligate'].inputs if x.material=='digest'},{R.TEMPLATE,R.STRUCTURAL})
    c('Experion aliquot is QC rather than production',g.terminals['experion_assay'].disposition,Disposition.QC)
    c('Ligation ATP and enzyme derive from same actual added solution',tuple(g.producer_of(k).id for k in ('ligase','ligation_atp')),('ligation_components','ligation_components'))
    c('Source loading leaves QA material outside acquisition',g.steps['loading_allocate'].inputs[0].material,'experion')
    c('Original and copied whole molecular geometry incomplete',bool(g.publication_problems()))
    original=Source(sources_root()/'gbs__10.1371+journal.pone.0019379'/'journal.pone.0019379_PMC3087801.xml.txt','Published PMC3087801 Methods83/105/Discussion197')
    table=Source(sources_root()/'gbs__10.1371+journal.pone.0019379'/'journal.pone.0019379_supp_pone.0019379.s004.docx.txt','Published TableS1 via standard source fetch')
    if have(c,original,table,label='Original component equations and barcode member'):
        text=original.path.read_text();codes=table.path.read_text().split()
        c('Selected barcode is in the published set',BARCODE in codes)
        c('Whole source adapter equation fixed flanks',('5′-'+il.TRUSEQ_READ1+'xxxx' in text and '5′-CWGyyyy'+revcomp(il.TRUSEQ_READ1) in text))
        c('Whole common strands and PCRprimers printed in original',all(x in text for x in ('CWG'+revcomp(COMMON),COMMON,FWD,REV)))
    c.report()
