"""Bounded, source-local molecular workflow; not a complete protocol review."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, panel, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, Terminal, Disposition, SourceItem
from pcr import PCRGeometry
from page import head, info
PAGE_FILE='chem.html'

def unknown_template(label):
    return Construct([Segment(label,'X'*30,placeholder=True)])

def clean(g,key,material,label,selection,kind=K.MOLECULAR):
    sel=g.input(key+'_selection',selection,K.SELECTION)
    kept=g.product(key+'_kept',label,kind)
    waste=g.product(key+'_waste','Separated material from '+key,K.WASTE)
    g.react(key,label,S.CLEANUP,((material,R.TEMPLATE),(sel,R.SELECTION)),(kept,waste))
    g.discard(waste,'Removed only at this source-described physical separation; recovery is not predicted')
    return kept

def pcr(g,key,template,forward,reverse,label,conditions):
    f=g.input(key+'_forward',forward,K.OLIGO)
    r=g.input(key+'_reverse',reverse,K.OLIGO)
    cat=g.input(key+'_enzyme',conditions,K.CATALYST)
    out=g.product(key+'_product',label,K.MOLECULAR)
    g.react(key,label,S.PCR,((template,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST)),(out,),
        geometry=PCRGeometry.unavailable(unknown_template('unassigned target sequence'),
            reason='This capped model does not assign the study-specific target substrate and complete primer-pool landing geometry.'))
    return out

def aliquot(g,key,material,label):
    selected=g.product(key+'_selected',label,K.MOLECULAR)
    remainder=g.product(key+'_remainder','Untransferred '+key+' material',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((material,R.TEMPLATE),),(selected,remainder))
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'The source specifies an aliquot; remaining reaction material is not claimed discarded.'))
    return selected
from adapter_ligation import AdapterLigationGeometry
SOURCE='https://doi.org/10.1038/s41592-018-0259-9; Online Methods, published DART-seq protocol and SupplementaryTable1'
TITLE='DART-seq — targeted capture, published protocol/IgG interface'
LIMIT='This is the targeted-capture DART-seq, distinct from the APOBEC m6A method. Only the terminal12T bead window and IgG extension interface are base-resolved; full bead architecture, template switching, adapter geometry and read coordinates remain deferred.'
IGG='GGGGAAGACCGATGGGCCCTTGGTGG'
SPLINT='CGGTCTTCCCCAAAAAAAAAAAA'

def ligation_scene():
    sc=Scene()
    sc.strand('bead',[Segment('terminal dT window','T'*12)],label='Terminal12T pairing window of bead oligo; upstream bead sequence omitted')
    sc.anneal('splint',[Segment('primer arm',SPLINT[:-12]),Segment('A arm',SPLINT[-12:])],to='bead',pair=('A arm','terminal dT window'))
    sc.anneal('custom',[Segment('5prime binding arm',IGG[:11]),Segment('target extension',IGG[11:])],to='splint',pair=('5prime binding arm','primer arm'),mod5='p')
    sc.ligation_end('custom','5prime binding arm',kind='splinted',polarity='5-prime',owner='adapter')
    sc.ligation_end('bead','terminal dT window',kind='splinted',polarity='3-prime',owner='substrate')
    return sc

def primer_splint_scene():
    sc=Scene()
    sc.strand('custom',[Segment('5prime binding arm',IGG[:11]),Segment('target extension',IGG[11:])],mod5='p')
    sc.anneal('splint',[Segment('primer arm',SPLINT[:-12]),Segment('A arm',SPLINT[-12:])],to='custom',pair=('primer arm','5prime binding arm'))
    return sc

def sealed_window_scene():
    sc=Scene()
    sc.strand('joined',[Segment('terminal dT window','T'*12),Segment('custom5prime arm',IGG[:11]),Segment('target extension',IGG[11:])],label='Ligated terminal bead window; upstream bead architecture omitted')
    sc.anneal('splint',[Segment('primer arm',SPLINT[:-12]),Segment('A arm',SPLINT[-12:])],to='joined',pair=('A arm','terminal dT window'))
    sc.junction('joined','terminal dT window','custom5prime arm','sealed bead3′–primer5′ phosphodiester junction')
    sc.note('splint','Noncovalent splint remains in reaction until source bead washes; shown interface is not a complete bead oligo')
    return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    bead=g.input('bead','Original Drop-seq oligo(dT) bead; terminal pairing window only represented',K.MOLECULAR,primary=True)
    custom=g.input('custom','IgG_CustomPrimer /5Phos/'+IGG,K.OLIGO,state=Construct([Segment('IgG_CustomPrimer',IGG)]))
    splint=g.input('splint','IgG_Splint',K.OLIGO,state=Construct([Segment('IgG_Splint',SPLINT)]))
    ligase=g.input('ligase','T4DNA ligase with PEG4000;37°C1h',K.CATALYST)
    ready=g.product('ready','IgG primer/splint duplex annealed95°C→14°C',K.MOLECULAR,state=primer_splint_scene())
    g.react('anneal','Anneal phosphorylated custom primer to splint',S.HYBRIDIZATION,
        ((custom,R.TEMPLATE),(splint,R.STRUCTURAL)),(ready,))
    joined=g.product('joined','Ligated bead3prime IgG extension plus excess/noncovalent oligos',K.MOLECULAR,state=sealed_window_scene())
    g.react('ligation','Join bead3primeOH to custom primer5primephosphate',S.LIGATION,
        ((bead,R.TEMPLATE),(ready,R.STRUCTURAL),(ligase,R.CATALYST)),(joined,),
        geometry=AdapterLigationGeometry.known('splinted',ligation_scene(),statement='Source IgG primer11-base interface and terminal12T pairing window only; no full bead-length or upstream architecture inference.'))
    heated=g.product('heated','Heat-treated ligation mixture, including free splint/oligos',K.MOLECULAR)
    heat=g.input('heat','65°C3min then ice',K.CONDITION)
    g.react('inactivate','Thermally inactivate ligase',S.THERMAL_TREATMENT,((joined,R.TEMPLATE),(heat,R.CONDITION)),(heated,))
    washed=clean(g,'bead_washes',heated,'Washed DART-seq beads','TE-SDS once and TE-TW twice; remove wash supernatants')
    cells=g.input('cells','Single B cells and transcripts in source droplet workflow',K.BIOLOGICAL,primary=True)
    tso=g.input('tso','Source template-switch oligo AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG',K.OLIGO)
    rt=g.input('rt_enzyme','Maxima reverse transcriptase; source droplet lysis and RT buffers',K.CATALYST)
    cdna=g.product('cdna','Target-enriched and ordinary polyadenylated cDNA with bead-derived cell identity/UMI',K.MOLECULAR)
    g.react('rt','Capture transcripts and reverse transcribe using modified beads',S.REVERSE_TRANSCRIPTION,
        ((cells,R.TEMPLATE),(washed,R.STRUCTURAL),(tso,R.STRUCTURAL),(rt,R.CATALYST)),(cdna,))
    washed_cdna=clean(g,'post_rt_washes',cdna,'Washed cDNA-bearing beads','TE-SDS,twiceTE-TW,Tris washes')
    exo=g.input('exo','ExonucleaseI;37°C45min',K.CATALYST)
    digested=g.product('digested','Exonuclease-treated beads plus digestion products',K.MOLECULAR)
    g.react('exonuclease','Digest unextended primers afterRT',S.DIGESTION,((washed_cdna,R.TEMPLATE),(exo,R.CATALYST)),(digested,))
    clean_cdna=clean(g,'post_exo_washes',digested,'Washed cDNA-bearing beads afterExoI','TE-SDS/TE-TW/water bead washes')
    smart=g.input('smart','One SMART PCR primer species',K.OLIGO,state=Construct([Segment('SMART PCR primer','AAGCAGTGGTATCAACGCAGAGT')]))
    polymerase=g.input('polymerase','KAPAHiFi;4cycles65°C then11cycles67°C annealing',K.CATALYST)
    amplicons=g.product('cdna_amplicons','Amplified source DART-seq cDNA',K.MOLECULAR)
    g.react('cdna_pcr','Single-species SMART-handle cDNA amplification',S.PCR,
        ((clean_cdna,R.TEMPLATE),(smart,R.STRUCTURAL),(polymerase,R.CATALYST)),(amplicons,),
        geometry=PCRGeometry.unavailable(unknown_template('full cDNA handle architecture'),symmetric_species=True,
            reason='Published protocol specifies one SMART primer species; full source cDNA inward-site architecture is deferred, not assumed from matching platform.'))
    amplicons=clean(g,'cdna_cleanup',amplicons,'Purified amplified cDNA','0.6×SPRI AMPureXP')
    tn=g.input('transposase','NexteraXT FC-131-1096 tagmentation enzyme',K.CATALYST)
    adapters=g.input('loaded_adapters','Historical source Nextera loaded adapters; ends unassigned',K.ADAPTER)
    tagged=g.product('tagged','Tagmented amplified cDNA',K.MOLECULAR)
    g.react('tagment','Nextera cDNA tagmentation',S.TAGMENTATION,((amplicons,R.TEMPLATE),(tn,R.CATALYST),(adapters,R.STRUCTURAL)),(tagged,))
    indexed=pcr(g,'index_pcr',tagged,'Source indexPCR primer1','Source indexPCR primer2','Indexed DART-seq products','Source Nextera libraryPCR')
    lib=clean(g,'library_cleanup',indexed,'Targeted-capture DART-seq library','Source library purification',K.LIBRARY)
    g.sequence(lib,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':20,'Read 2':130},limitation='Main Methods and published detailed protocol specify NextSeq500,Read1=20,Read2=130; detailed protocol also specifies customRead1 primer. Primer landing and historical full-bead barcode/UMI widths remain deferred.'))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Full source primer-pool, intermediate and final-strand assembly',LIMIT+' Source material is available for follow-up; not an absence finding.'))
    return g

def detail():return panel(ligation_scene().rows(),caption='Source IgG splinted-ligation interface (partial bead window)')+panel(sealed_window_scene().rows(),caption='Separate sealed product, with noncovalent splint still present')

REVIEWED=('bead', 'custom', 'splint', 'ligase', 'ready', 'joined', 'heat', 'heated', 'bead_washes_selection', 'bead_washes_kept', 'bead_washes_waste', 'cells', 'tso', 'rt_enzyme', 'cdna', 'post_rt_washes_selection', 'post_rt_washes_kept', 'post_rt_washes_waste', 'exo', 'digested', 'post_exo_washes_selection', 'post_exo_washes_kept', 'post_exo_washes_waste', 'smart', 'polymerase', 'cdna_amplicons', 'cdna_cleanup_selection', 'cdna_cleanup_kept', 'cdna_cleanup_waste', 'transposase', 'loaded_adapters', 'tagged', 'index_pcr_forward', 'index_pcr_reverse', 'index_pcr_enzyme', 'index_pcr_product', 'library_cleanup_selection', 'library_cleanup_kept', 'library_cleanup_waste', 'anneal', 'ligation', 'inactivate', 'bead_washes', 'rt', 'post_rt_washes', 'exonuclease', 'post_exo_washes', 'cdna_pcr', 'cdna_cleanup', 'tagment', 'index_pcr', 'library_cleanup')
def extra_checks(c,g):
    from storage import sources_root
    source=sources_root()/Path(__file__).resolve().parents[1].name
    original=source/SOURCE_FILE
    if original.exists():
        text=''.join(original.read_text().split())
        for oligo in TRANSCRIBED:
            c('Verified original transcription '+oligo,oligo in text)
    else:c.skip('Original transcription checks','External original unavailable')
    source_checks(c,g)


SOURCE_FILE='s41592-018-0259-9_supp_NIHMS1512056-supplement-Sup_Table_1.xlsx.txt'
TRANSCRIBED=('/5Phos/'+IGG,SPLINT)
def source_checks(c,g):
    from chemdraw import revcomp
    c('IgG splint arm is reverse complement of primer5prime arm',revcomp(IGG[:11]),SPLINT[:-12])
    c('Sealed strand appends custom primer once', ''.join(x.top for x in sealed_window_scene().strands['joined'].segs),'T'*12+IGG)
    c('Geometry is a source-known partial splinted interface',g.steps['ligation'].geometry.kind,'splinted')
    sc=ligation_scene()
    c('Ligation termini are adjacent in computed scene',sc.strands['bead'].span('terminal dT window')[1],sc.strands['custom'].span('5prime binding arm')[0])

def render_page():
    return '\n'.join([head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_dart_targeted.html">Research notes</a></p>',
        info(LIMIT),detail(),''.join(workflow_panel(workflow_from_protocol_graph(protocol_graph(),endpoint)) for endpoint in ['library_cleanup_kept', 'bead_washes_kept', 'ready']),'</div>'])

def selftest():
    from checks import Check
    c=Check(); g=protocol_graph()
    issues=g.publication_problems()
    c('Only explicit unassembled-source coverage blocks publication',issues,["authoritative source item 'unassembled_source' is not accounted for"])
    c('Private preview does not claim protocol completion','data-protocol-complete="true"' in render_page(),False)
    extra_checks(c,g)
    c.report(exit_on_fail=True)
