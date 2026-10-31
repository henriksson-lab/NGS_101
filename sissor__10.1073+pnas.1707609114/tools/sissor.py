"""Published SISSOR strand separation/partitioning and two source library routes."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
from chemdraw import Scene,Segment,revcomp,feature,workflow_from_protocol_graph,workflow_panel
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SequencingPlan,SourceItem
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from page import head,info
PAGE_FILE='sissor.html'
SOURCE='Chu et al.2017 DOI10.1073/pnas.1707609114,publishedPMC5703283 mainResults64–74/Methods164,SIAppendixI.1–I.5'
CHAMBERS=24
PRIMER_NOTE='Main Methods164/originalPDFp6 prints 84µM “3′ phosphorylated random hexamers”. That wording is retained;full source priming/initiation chemistry is unassembled. No ordinary3′OH,phosphorothioate protection or deblocking reaction is substituted.'
CHAMBER_INDEX=feature('sissor.chamber_library','sample_index','unknown',note='Identifies a separately recovered chamber library;not a unique genomic molecule or cell barcode.')

@dataclass(frozen=True)
class Strand:
    ident:str
    dna5to3:str
    inferred:bool=True
    def __post_init__(self):
        if not self.ident or not self.dna5to3 or set(self.dna5to3)-set('ACGT'):raise ValueError('Explicit finite DNA strand and stable identity required')

def separate_duplex(ident,top,*,inferred=True):
    """Caller-supplied finite duplex example;not an inferred human genome sequence."""
    first=Strand(ident+'.top',top,inferred)
    second=Strand(ident+'.bottom',revcomp(top),inferred)
    # Verify the relation by placement before returning either physical strand.
    sc=Scene();sc.strand('top',[Segment('top',first.dna5to3,inferred=inferred)])
    sc.anneal('bottom',[Segment('bottom',second.dna5to3,inferred=inferred)],to='top',pair=('bottom','top'))
    return first,second

@dataclass(frozen=True)
class Partition:
    chambers:tuple[tuple[Strand,...],...]
    unrecovered:tuple[Strand,...]

def partition_selected(strands,assignment,*,n_chambers=CHAMBERS):
    """Explicit caller-selected destinations,including None;not a stochastic predictor."""
    strands=tuple(strands);ids=[x.ident for x in strands]
    if not isinstance(n_chambers,int) or isinstance(n_chambers,bool) or n_chambers<1:raise ValueError('Positive chamber count required')
    if len(set(ids))!=len(ids) or set(assignment)!=set(ids):raise ValueError('Account for every supplied strand exactly once')
    bins=[[] for _ in range(n_chambers)];lost=[]
    for strand in strands:
        dest=assignment[strand.ident]
        if dest is None:lost.append(strand)
        elif not isinstance(dest,int) or isinstance(dest,bool) or not 0<=dest<n_chambers:raise ValueError('Destination outside explicitly supplied chambers')
        else:bins[dest].append(strand)
    return Partition(tuple(tuple(x) for x in bins),tuple(lost))

def unknown(label):
    sc=Scene();sc.strand('source',[role(label)],label='Source molecular assembly boundary;not absent-disclosure evidence');return sc

def index_inventory(label):
    sc=Scene();sc.strand('index',[role(label,molecular_feature=CHAMBER_INDEX)],label='Known chamber-index meaning;sequence,width,member and strand order unassembled');return sc

def graph_for_route(route='core'):
    if route not in {'core','transposon'}:raise ValueError('Choose a source-described library route')
    g=ProtocolGraph('SISSOR / '+route,authority=SOURCE,representation='published_instance')
    def inp(key,label,kind=K.CATALYST,state=None,**kw):return g.input(key,label,kind,state=state,**kw)
    def op(key,x,label,kind,reagents=(),**kw):
        y=g.product(key+'_mixture',label+' preparation;original and reaction species retained until actual workup',K.MOLECULAR)
        templates=tuple((part,R.TEMPLATE) for part in x) if isinstance(x,tuple) else ((x,R.TEMPLATE),)
        g.react(key,label,kind,(*templates,*reagents),(y,),**kw);return y
    def separate(key,x,label,reagents):
        strands=tuple(g.product(key+'_'+side,'Separated '+side+' DNA strand population in the SAME source mixture;residual duplex/lysate fate unquantified',K.MOLECULAR) for side in ('first','complementary'))
        g.react(key,label,S.STRAND_SEPARATION,((x,R.TEMPLATE),*reagents),strands);return strands
    def clean(key,x,label):
        condition=inp(key+'_procedure',label,K.SELECTION);y=g.product(key,'Recovered source DNA;carryover/recovery unquantified',K.MOLECULAR)
        waste=g.product(key+'_reject','Physically removed workup fractions;exact composition unassigned',K.WASTE)
        g.react(key+'_step',label,S.CLEANUP,((x,R.TEMPLATE),(condition,R.SELECTION)),(y,waste));g.discard(waste,'Actual source physical separation;no perfect recovery or chemical-class purity predicted.');return y
    def pcr(key,x,f,r,poly,label):
        return op(key,x,label,S.PCR,((f,R.STRUCTURAL),(r,R.STRUCTURAL),(poly,R.CATALYST)),geometry=PCRGeometry.unavailable(Scene(),reason='Source reaction/primer roles known;full source primers,genomic strands and historical kit end geometry unassembled.'))
    cell=inp('cell',('SourcePGP1Cell3 CoREroute' if route=='core' else 'SourcePGP1Cell2 transposonroute')+';culture/harvest outside graph',K.BIOLOGICAL,primary=True)
    device=inp('device','Pre-coated PDMS device:sourceBSA/randomhexamer/dNTP coating,purge/dry,UV15min BEFOREcellloading;full hardwarepreparation outside graph',K.CONDITION)
    x=op('capture',cell,'Source-imaged cell capture in the prepared device',S.GENERIC,((device,R.CONDITION),),changes_molecule=False)
    als=inp('als','SourceALS400mMKOH/10mMDTT/1%Tween20;ringmixer5Hz10min20°C',K.CONDITION)
    x=separate('lyse_separate',x,'Source alkaline cell lysis and chromosomal strand separation',((als,R.CONDITION),))
    representative=g.product('representative','One representative source partition;actual strand/homolog mixture unknown',K.MOLECULAR,state=unknown('Separated source ssDNA fragments in one chamber'))
    remaining=g.product('remaining_partitions','Other23source chamber preparations;actual strand/homolog mixtures unknown',K.MOLECULAR)
    unrecovered=g.product('unrecovered','Possible DNA/material not recovered to anychamber;amount unknown,may be empty',K.MOLECULAR)
    g.react('partition','Randomly distribute source fragments into24chambers;one representative class expanded',S.PARTITION,tuple((part,R.TEMPLATE) for part in x),(representative,remaining,unrecovered))
    g.mark_terminal(Terminal(unrecovered.id,Disposition.UNRESOLVED,'Source discusses possible material loss;not assigned a quantity,genotype or disposal fate.'))
    ns=inp('ns','Source on-chipNS400mMHCl/600mMTrisHCl/1%Tween20',K.CONDITION)
    x=op('neutralize',representative,'Neutralize source partition and push into MDA chamber',S.QUENCHING,((ns,R.CONDITION),))
    mix=inp('mda_mix','Actual source NexteraMDA mastermix allocation; source1xbuffer/84uM3prime-phosphorylatedhexamers/2mMdNTPs/150uMdUTP/0.84ugPhi29; stockvolume/per-chamberdose and initiation unassembled',K.MOLECULAR)
    primer=g.product('mda_hexamers','Primer component of same actual MDA stock. '+PRIMER_NOTE,K.OLIGO,state=unknown('Source randomhexamer/MDA end notation unresolved'))
    phi=g.product('mda_phi29','Phi29 component of same actual MDA stock; not an additional enzyme dose',K.CATALYST)
    nt=g.product('mda_nucleotides','Copresent dNTP/dUTP molecular substrates in same actual MDA stock',K.MOLECULAR)
    buf=g.product('mda_buffer','Buffer component of same actual MDA stock; historical formulation unassembled',K.CONDITION)
    g.react('mda_components','Functional component views of one actual source MDA stock',S.MAGIC,((mix,R.TEMPLATE),),(primer,phi,nt,buf),changes_molecule=False,reason='Source main164 explicitly lists these copresent MasterMix components. This is inventory, not physical separation, independent fresh doses or a verified molecular initiation mechanism.')
    x=op('mda',x,'Source MDA30°C15h in separate chambers;original templates and dU-containing copies/carryover',S.MAGIC,((primer,R.STRUCTURAL),(phi,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION)),reason=PRIMER_NOTE+' Observed amplification is source-reported;this is not a modeled extension from a guessed primer terminus.')
    retrieve=inp('retrieval','Flushfeedlinesair/TE;push eachchamber product to its own tip withTE;collect5µL',K.CONDITION)
    x=op('retrieve',x,'Separately recover one source chamber product;no chamber pooling',S.GENERIC,((retrieve,R.CONDITION),),changes_molecule=False)
    heat=inp('inactivate_phi29','65°C15min',K.CONDITION);x=op('mda_heat',x,'Heat-inactivate Phi29;not physical purification',S.THERMAL_TREATMENT,((heat,R.CONDITION),))
    if route=='core':
        q=inp('core_enzymes','SourceUDG+EndonucleaseIV/exo-Klenowbuffer37°C2h,65°C15min')
        x=op('core_incision',x,'CoRE uracil-base excision/incision in copied DNA;not entire-strand disposal',S.DIGESTION,((q,R.CATALYST),))
        pol=inp('core_klenow','Source exo-Klenow+dNTPs37°C1h,75°C15min')
        x=op('core_fill',x,'Source Klenow generation of double-stranded fragmentation products from nick/gap positions',S.MAGIC,((pol,R.CATALYST),),reason='SI operational CoRE description known;exact AP/end/flap/strand-displacement products not assembled. Released DNA and nucleobase species persist until actual recovery.')
        adapter=inp('neb_loop','Historical sourceNEBNext loop adapter containing uracil;full polymer/ends unreviewed',K.ADAPTER,state=unknown('Source NEBNext loop adapter'))
        ligase=inp('core_ligation','SourceT4DNAligase/KAPArapidbuffer25°C10min')
        geometry=AdapterLigationGeometry('unavailable','unavailable',(unknown('Source CoRE genomic/adapter end geometry not assembled'),),'Source names loop-adapter ligation;historical adapter/genomic junction geometry unassembled,not claimed unreported.')
        x=op('core_ligate',x,'Source NEBNext loop-adapter ligation',S.LIGATION,((adapter,R.STRUCTURAL),(ligase,R.CATALYST)),geometry=geometry)
        user=inp('loop_user','SourceUSER37°C15min')
        x=op('loop_open',x,'SourceUSER processes uracil-containing loop;other susceptible DNA/products not assigned vanished',S.DIGESTION,((user,R.CATALYST),))
        x=clean('loop_cleanup',x,'Source50µLAMPureXP purification beforePCR')
        f=inp('core_f','SourceforwardPCRprimer;historical exactpolymer unreviewed',K.OLIGO)
        r=inp('core_index','SourceNEBNext sequencing index primer;actualchambermember/width unassigned',K.OLIGO,state=index_inventory('Chamber sequencing index'))
        poly=inp('core_poly','SourceKAPASYBRmastermix;12cycles98/65/72')
        x=pcr('core_pcr',x,f,r,poly,'Source12cycle indexed PCR')
        out=clean('core_library',x,'PAGEgel recovery300–700bp')
    else:
        als2=inp('library_als','SourceSIlibraryALS400mMKOH/10mMDTT/1%Tween20,roomtemperature3minthenice',K.CONDITION)
        x=separate('mda_denature',x,'Denature source MDA amplicon aliquot for second-strand processing',((als2,R.CONDITION),))
        ns2=inp('library_ns','SourceSI NS400mMHCl/600mMTrisHCl;do not importonchipTweenaddition',K.CONDITION)
        x=op('library_neutralize',x,'Neutralize source library aliquot',S.QUENCHING,((ns2,R.CONDITION),))
        hexa=inp('unmodified_hexamers','SI250ng unmodifiedrandomhexamers;distinct from mainMDAendnotation',K.OLIGO)
        pol1=inp('pol1','DNApolymeraseI/Ampligasebuffer/NEBbuffer2/dNTPs37°C1h,65°C10min')
        x=op('second_strand',x,'Source PolI treatment creates transposase-accessible duplex DNA',S.EXTENSION,((hexa,R.STRUCTURAL),(pol1,R.CATALYST)))
        seal=inp('ampligase','Ampligase37°C10minthen65°C10min')
        x=op('nick_seal',x,'Source nick sealing;chemical debris remain',S.NICK_SEALING,((seal,R.CATALYST),))
        x=clean('ethanol_recovery',x,'Sourceethanol precipitation/elution')
        tag=inp('nextera','Historical Nexteratransposase/HMWbuffer55°C5min')
        tags=inp('tag_adapters','Unassembled source transposase-loaded adapter components',K.ADAPTER)
        x=op('tagment',x,'Source transposon tagmentation',S.TAGMENTATION,((tag,R.CATALYST),(tags,R.STRUCTURAL)))
        protease=inp('transposase_protease','Sourceprotease50°C10min/70°C20min')
        x=op('inactivate_transposase',x,'Protease inactivate transposase;no physical protein removal asserted',S.DEPROTEINIZATION,((protease,R.CATALYST),))
        fill=inp('tag_klenow','Exo-Klenow+dNTPs37°C15min/75°C20min')
        x=op('tag_fill',x,'Source single-stranded end filling;full product geometry unassembled',S.END_REPAIR,((fill,R.CATALYST),))
        a=inp('adapter1','Source adapter1PCRprimer;historical source polymerunreviewed',K.OLIGO)
        b=inp('barcode_adapter','Sourcechamber barcodeadapterPCRprimer;member/widthunassigned',K.OLIGO,state=index_inventory('Chamber barcodeadapter index'))
        kapa=inp('tag_poly','SourceKAPASYBRfast mastermix')
        x=pcr('tag_pcr1',x,a,b,kapa,'First source7cyclePCR adds adapter/index sequences')
        p1=inp('illumina1','SourceIlluminaprimer1;fullsourcepolymerunreviewed',K.OLIGO)
        p2=inp('illumina2','SourceIlluminaprimer2;fullsourcepolymerunreviewed',K.OLIGO)
        x=pcr('tag_pcr2',x,p1,p2,kapa,'Second source7cyclePCR;no unprinted interveningcleanup added')
        out=clean('tag_library',x,'SourcePAGEQC/AMPurepurification then300–700bpgel recovery')
    parallel=g.product('other_chamber_libraries','Other23individually amplified/barcoded same-cell chamber libraries;matched source route,full molecular geometry unassembled',K.LIBRARY)
    g.react('parallel_preparations','Summarize known parallel neutralization/MDA/recovery/'+route+' library processing forother23chambers',S.MAGIC,((remaining,R.TEMPLATE),),(parallel,),reason='Same source operation program as representative chamber;separate molecular assignments,barcodes and exact products not individually expanded. No unrecovered DNA is assigned to these libraries.')
    pooled=g.product('cell_pool','Physical pool of separately barcoded chamber libraries from the SAME source cell',K.LIBRARY)
    g.react('pool_libraries','Source combine separately indexed same-cell chamber libraries before sequencing',S.POOL,((out,R.TEMPLATE),(parallel,R.TEMPLATE)),(pooled,))
    g.sequence(pooled,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':100,'Read 2':100},limitation='Source main reports deepPE100. SIpilotSE36/50 andadditionalSE100/PE100 acquisitions are separate options. Full historical primers/adapter/index landing and exact per-run assignments remain unassembled.'))
    g.declare_reviewed_scope(COMMON_SCOPE+(CORE_SCOPE if route=='core' else TRANSPOSON_SCOPE),SOURCE)
    g.add_source_item(SourceItem('historical_geometry','Historical MDA/end chemistry,NEBNext/tagmentation primers and product/read geometry','Source main164 andSI I.3–I.5;citedPeters2012/Adey2010 originals remain unreviewed.'))
    return g


# Literal source inventory from main Fig1/Results74/Methods164 and SI I.3-I.5.
# Source-dependent unresolved polymers stay behind historical_geometry.
COMMON_SCOPE=tuple('cell device capture als lyse_separate partition ns neutralize mda_mix mda_components mda retrieval retrieve inactivate_phi29 mda_heat parallel_preparations pool_libraries cell_pool'.split())
CORE_SCOPE=tuple('core_enzymes core_incision core_klenow core_fill neb_loop core_ligation core_ligate loop_user loop_open loop_cleanup_procedure loop_cleanup_step core_f core_index core_poly core_pcr core_library_procedure core_library_step'.split())
TRANSPOSON_SCOPE=tuple('library_als mda_denature library_ns library_neutralize unmodified_hexamers pol1 second_strand ampligase nick_seal ethanol_recovery_procedure ethanol_recovery_step nextera tag_adapters tagment transposase_protease inactivate_transposase tag_klenow tag_fill adapter1 barcode_adapter tag_poly tag_pcr1 illumina1 illumina2 tag_pcr2 tag_library_procedure tag_library_step'.split())

def protocol_graph():return graph_for_route('core')
def render_page():
    return head('SISSOR — source-limited strand/compartment model')+'<main class="wrap"><h1>SISSOR</h1><p><a href="01_sissor.html">Research notes</a></p>'+info('Private incomplete published-source model. Two source library routes are independent options for chamber products. No loss-free strand separation,one-haplotype chamber,full primer geometry or error-free amplification is claimed.')+info(PRIMER_NOTE)+''.join('<h2>'+r+'</h2>'+workflow_panel(workflow_from_protocol_graph(graph_for_route(r),'cell_pool')) for r in ('core','transposon'))+'</main>'

def selftest():
    from checks import Check
    from storage import sources_root
    c=Check();a,b=separate_duplex('example','ACGTCAGGTTACGATC');p=partition_selected((a,b),{a.ident:0,b.ident:1})
    c('Complementary fragments can occupy different chambers',(p.chambers[0],p.chambers[1]),((a,),(b,)))
    same=partition_selected((a,b),{a.ident:0,b.ident:0})
    c('A chamber may contain multiple template strands',len(same.chambers[0]),2)
    missing=partition_selected((a,b),{a.ident:0,b.ident:None})
    c('Unrecovered DNA remains accounted separately',missing.unrecovered,(b,))
    try:partition_selected((a,b),{a.ident:0});bad=False
    except ValueError:bad=True
    c('No implicit strand disappearance with incomplete assignments',bad,True)
    try:partition_selected((a,b),{a.ident:0,b.ident:24});bad=False
    except ValueError:bad=True
    c('Destinations must be within the explicit chamber set',bad,True)
    c('Caller-supplied unverified sequences are inferred by default',a.inferred,True)
    g=protocol_graph();t=graph_for_route('transposon')
    c('Capture/retrieval preserve molecular identity',(g.steps['capture'].changes_molecule,g.steps['retrieve'].changes_molecule),(False,False))
    c('One stock supplies primer and nucleotide components',tuple(g.steps['mda_components'].outputs),('mda_hexamers','mda_phi29','mda_nucleotides','mda_buffer'))
    c('MDA nucleotide contribution is molecular substrate',next(i.role for i in g.steps['mda'].inputs if i.material=='mda_nucleotides'),R.STRUCTURAL)
    c('Required historical geometry remains uncovered',(g.publication_problems(),t.publication_problems()),(["authoritative source item 'historical_geometry' is not accounted for"],)*2)
    c('Both source routes remain incomplete',(g.publishable,t.publishable),(False,False))
    c('No finished-page marker','data-protocol-complete=' in render_page(),False)
    c('Separate bulk library routes do not physically combine','core_pcr' not in t.steps and 'tag_pcr1' not in g.steps,True)
    c('Protease does not erase transposase digestion material',t.steps['tag_fill'].inputs[0].material,'inactivate_transposase_mixture')
    main=sources_root()/HERE.parent.name/'pnas.1707609114_PMC5703283.xml.txt'
    if main.exists():
        text=main.read_text();c('Original end wording retained without replacing it','3′ phosphorylated random hexamers' in text,True)
    else:c.skip('Published end-notation transcription','External source unavailable')
    c.report()
