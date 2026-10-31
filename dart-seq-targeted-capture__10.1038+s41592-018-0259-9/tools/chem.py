"""Bounded, source-local molecular workflow; not a complete protocol review."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, panel, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, Terminal, Disposition, SourceItem
from pcr import PCRGeometry
from protocol_model import role
import illumina as il
from page import head, info
PAGE_FILE='chem.html'

def unknown_template(label):
    return Construct([role(label+'; occupied length and sequence unassembled')])

def clean(g,key,material,label,selection,kind=K.MOLECULAR):
    sel=g.input(key+'_selection',selection,K.SELECTION)
    kept=g.product(key+'_kept',label,kind,chemical_state=material.chemical_state)
    waste=g.product(key+'_waste','Separated material from '+key,K.WASTE)
    g.react(key,label,S.CLEANUP,((material,R.TEMPLATE),(sel,R.SELECTION)),(kept,waste))
    g.discard(waste,'Removed only at this source-described physical separation; recovery is not predicted')
    return kept

def pcr(g,key,template,forward,reverse,label,conditions):
    f=g.input(key+'_forward',forward,K.OLIGO,state=Construct([Segment('Source P5hybrid',P5_HYBRID)]) if key=='index_pcr' else None,chemical_state={'source_terminal_bonds':'T*A*C; two phosphorothioate bonds'} if key=='index_pcr' else None)
    r=g.input(key+'_reverse',reverse,K.OLIGO)
    cat=g.input(key+'_enzyme',conditions,K.CATALYST)
    nt=g.input(key+'_nucleotides','Aggregate source index-amplification nucleotide supplies; formula/dose/packaging unassembled',K.MOLECULAR)
    out=g.product(key+'_product',label,K.MOLECULAR)
    g.react(key,label,S.PCR,((template,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST),(nt,R.STRUCTURAL)),(out,),
        geometry=PCRGeometry.unavailable(unknown_template('unassigned target sequence'),
            reason='This capped model does not assign the study-specific target substrate and complete primer-pool landing geometry.'))
    return out

def aliquot(g,key,material,label):
    selected=g.product(key+'_selected',label,material.kind,state=material.state,chemical_state=material.chemical_state)
    remainder=g.product(key+'_remainder','Untransferred '+key+' material; further actual reactions/amounts unassembled',material.kind,state=material.state,chemical_state=material.chemical_state)
    g.react(key,label,S.SPLIT,((material,R.TEMPLATE),),(selected,remainder))
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'The source specifies an aliquot; remaining reaction material is not claimed discarded.'))
    return selected
from adapter_ligation import AdapterLigationGeometry
SOURCE='https://doi.org/10.1038/s41592-018-0259-9; exactpublishedMOESM3 DART-seqprotocol andTable1; associatedAMPMC6378878 context, no whole finalVOR/AMequivalence'
TITLE='DART-seq — targeted capture, published protocol/IgG interface'
LIMIT='This is the targeted-capture DART-seq, distinct from the APOBEC m6A method. Published detailedprotocol supplies fullbead baseclasses and knownprimers; theterminal12T IgG Scene is a partial display only. Fullsourcepool/intermediate/finalstrandassembly remains REQUIRED and unimplemented; actualtagmentedproducts/readlandings unresolved. No sourceabsence inferred.'
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
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    bead=g.input('bead','Original Drop-seq oligo(dT) bead; terminal pairing window only represented',K.MOLECULAR,primary=True,state=source_bead(),chemical_state={'source_written':BEAD_WRITTEN,'attachment':'Bead linkage chemistry unassembled; source J12cellbarcode/N8UMI/T30 are exact declared classes'})
    bead_first=g.product('bead_first','Source12000beads for representative first ligation',K.MOLECULAR,state=bead.state)
    bead_others=g.product('bead_others','Actual remaining5–9 ligation-reaction bead contributions;12000perreaction; count within source6–10total range unassigned',K.MOLECULAR,state=bead.state)
    bead_rest=g.product('bead_rest','Unused initial bead stock, if any; amount/fate unassigned',K.MOLECULAR)
    g.react('bead_array_allocation','Allocate source6–10actualreaction contributions without choosing an invented exactcount',S.SPLIT,((bead,R.TEMPLATE),),(bead_first,bead_others,bead_rest));g.mark_terminal(Terminal(bead_rest.id,Disposition.UNRESOLVED,'Source unused stock outside selected6–10reactionarray unquantified.'))
    custom=g.input('custom','IgG_CustomPrimer /5Phos/'+IGG,K.OLIGO,state=Construct([Segment('IgG_CustomPrimer',IGG)]),chemical_state={'5prime':'source5Phos','source_ordering':'/5Phos/'+IGG})
    splint=g.input('splint','IgG_Splint',K.OLIGO,state=Construct([Segment('IgG_Splint',SPLINT)]))
    ligase=g.input('ligase','Source2µlT4DNA ligase per representative reaction;37C1h',K.CATALYST)
    ready=g.product('ready','IgG primer/splint duplex annealed95°C→14°C',K.MOLECULAR,state=primer_splint_scene())
    g.react('anneal','Anneal phosphorylated custom primer to splint',S.HYBRIDIZATION,
        ((custom,R.TEMPLATE),(splint,R.STRUCTURAL)),(ready,))
    others=g.input('other_toeholds','Other actual source annealed Bcell-panel toehold duplexes; whole pool assemblies not implemented',K.ADAPTER)
    pool=g.product('toehold_pool','Actual equimolar source toehold pool; selected IgG component base-resolved, other members unexpanded',K.MOLECULAR)
    g.react('toehold_pooling','Pool actual source annealed probe species; no synthetic IgG-only assay invented',S.POOL,((ready,R.TEMPLATE),(others,R.TEMPLATE)),(pool,))
    first_probe=g.product('first_probe','Source16µl pooled probe contribution for representative reaction',K.MOLECULAR)
    other_probe=g.product('other_probe','Source16µl per actual remaining5–9 reactions; cardinality unassigned',K.MOLECULAR)
    probe_rest=g.product('probe_rest','Untransferred prepared probe stock, if any; amount/fate unassigned',K.MOLECULAR)
    g.react('probe_array_allocation','Select actual per-reaction probe contributions',S.SPLIT,((pool,R.TEMPLATE),),(first_probe,other_probe,probe_rest));g.mark_terminal(Terminal(probe_rest.id,Disposition.UNRESOLVED,'Untransferred probe stock outside actual selected reactions.'))
    ligation_conditions=g.input('ligation_conditions','Source perreaction:16µlprobe+40µl50percentPEG4000+40µlT4buffer+72µlwater+2µlligase;12000beads',K.CONDITION)
    joined=g.product('joined','Ligated bead3prime IgG extension plus excess/noncovalent oligos',K.MOLECULAR,state=sealed_window_scene(),chemical_state={'selected_full_source_bead':source_bead(extended=True),'other_source_pool_members':'copresent, unimplemented'})
    g.react('ligation','Join bead3primeOH to custom primer5primephosphate',S.LIGATION,
        ((bead_first,R.TEMPLATE),(first_probe,R.STRUCTURAL),(ligase,R.CATALYST),(ligation_conditions,R.CONDITION)),(joined,),
        geometry=AdapterLigationGeometry.known('splinted',ligation_scene(),statement='Source IgG primer11-base interface and terminal12T pairing window only; no full bead-length or upstream architecture inference.'))
    heated=g.product('heated','Heat-treated ligation mixture, including free splint/oligos',K.MOLECULAR,chemical_state=joined.chemical_state)
    heat=g.input('heat','65°C3min then ice',K.CONDITION)
    g.react('inactivate','Thermally inactivate ligase',S.THERMAL_TREATMENT,((joined,R.TEMPLATE),(heat,R.CONDITION)),(heated,))
    other_ligase=g.input('other_ligation_supplies','Actual repeated source2µlT4ligase per each remaining5–9reaction; stock counts/formulation unassigned',K.CATALYST)
    other_products=g.product('other_ligation_products','Actual remaining ligation products after37C1h/65C3min/ice; array not resolved into invented exact tube count',K.MOLECULAR)
    g.react('other_ligations','Source remaining parallel ligation/cooling array',S.MAGIC,((bead_others,R.TEMPLATE),(other_probe,R.STRUCTURAL),(other_ligase,R.CATALYST)),(other_products,),reason='Source specifies6–10parallelreactions total. One representative reaction plus actual remaining5–9reactionarray, full pool products unimplemented; no fixed cardinality/recovery inferred.')
    pooled_beads=g.product('pooled_beads','Actual source6–10cooledligation preparations pooled',K.MOLECULAR,chemical_state={'selected_full_source_bead':source_bead(extended=True),'other_source_pool_members':'actual remaining source array and other members unimplemented'})
    g.react('bead_product_pool','Pool actual distinct representative and remaining reaction products',S.POOL,((heated,R.TEMPLATE),(other_products,R.TEMPLATE)),(pooled_beads,))
    washed=clean(g,'bead_washes',pooled_beads,'Washed DART-seq beads','TE-SDS once and TE-TW twice; remove wash supernatants')
    cells=g.input('cells','Single B cells and transcripts in source droplet workflow',K.BIOLOGICAL,primary=True)
    capture_cond=g.input('droplet_conditions','Source microfluidic cell/bead co-encapsulation/lysis and RNAhybridization; occupancy/yield unassigned',K.CONDITION)
    emulsion=g.product('captured_emulsion','Source droplet mixture with bead-captured RNA and other cellular material; not yet reverse-transcribed',K.MOLECULAR)
    g.react('droplet_capture','Celllysis/RNAcapture on modified and unmodified bead oligos',S.BINDING,((cells,R.TEMPLATE),(washed,R.STRUCTURAL),(capture_cond,R.CONDITION)),(emulsion,))
    breaker=g.input('emulsion_break_conditions','Source SSC/perfluorooctanol emulsionbreak and physical bead recovery',K.CONDITION)
    recovered=g.product('recovered_bead_rna','Actual RNA-bearing recovered bead stock; background/cellmaterial carryover unassigned',K.MOLECULAR,chemical_state=washed.chemical_state)
    oil=g.product('oil_reject','Source physically removed excess oil/phase fractions',K.WASTE)
    g.react('break_emulsion','Break emulsion and collect aqueous RNA-bearing beads BEFORE RT',S.EXTRACTION,((emulsion,R.TEMPLATE),(breaker,R.CONDITION)),(recovered,oil));g.discard(oil,'Source physically removed oil/phase fractions; not computational filtering.')
    recovered=clean(g,'pre_rt_washes',recovered,'Recovered source RNA-bearing beads ready for bulkRT','Source6xSSC and5xMaxima-buffer washes; no complete-purity inference')
    template=g.product('rt_template','Same-tube bead-associated capturedRNA template component; whole RNA unassembled',K.MOLECULAR)
    bead_primers=g.product('rt_bead_primers','Same recovered-bead-stock primer components, not fresh free oligos; selectedIgG fullbead component shown',K.OLIGO,state=source_bead(extended=True))
    g.react('rt_bead_components','Logical inventory of recovered bead/RNA stock; components are not physical fractions',S.MAGIC,((recovered,R.TEMPLATE),),(template,bead_primers),changes_molecule=False,reason='Source RTuses primers already on captured beads, not an additional free beadprimer dose. Other actual sourceprimer classes and occupied RNA remain unassembled.')
    mix=g.input('rt_mix','Source fresh200µlpreparedRTmix perbeadaliquot:80water/40Maxima5xbuffer/40Ficoll/20dNTP/5RNaseinhibitor/10Maxima/5TSO',K.MOLECULAR)
    tso=g.product('tso','Copresent source5µl100µMTSO, template-switchacceptor; not a new RNApriming stock',K.OLIGO,state=Construct([Segment('TSO DNA prefix',TSO_DNA),Segment('terminal RNA guanosines','GGG',note='Source rGrGrG')]),chemical_state={'sugar':'DNA prefix plus source terminalrGrGrG'})
    rt=g.product('rt_enzyme','Copresent source10µlMaximaH-RT catalyst component',K.CATALYST)
    nt=g.product('rt_nucleotides','Copresent source20µl10mMdNTP component',K.MOLECULAR)
    rtbuf=g.product('rt_buffer','Copresent sourceMaxima/Ficoll/buffer conditions',K.CONDITION)
    g.react('rt_mix_components','Inventory of one actual source freshly prepared200µlRTmix',S.MAGIC,((mix,R.TEMPLATE),),(tso,rt,nt,rtbuf),changes_molecule=False,reason='Explicit sourcepreparedRTmix; same-stock components, not added fresh separately.')
    cdna=g.product('cdna','Target-enriched and ordinarypolyA cDNA from actualpost-break bulkRT/template-switch preparation; full occupied products unassembled',K.MOLECULAR)
    g.react('rt','Source bulkRT after emulsionbreak/washes:30minroomtemperature then45C2h',S.REVERSE_TRANSCRIPTION,((template,R.TEMPLATE),(bead_primers,R.STRUCTURAL),(tso,R.STRUCTURAL),(nt,R.STRUCTURAL),(rt,R.CATALYST),(rtbuf,R.CONDITION)),(cdna,))
    washed_cdna=clean(g,'post_rt_washes',cdna,'Washed cDNA-bearing beads','TE-SDS,twiceTE-TW,Tris washes')
    exo=g.input('exo','ExonucleaseI;37°C45min',K.CATALYST)
    digested=g.product('digested','Exonuclease-treated beads plus digestion products',K.MOLECULAR)
    g.react('exonuclease','Digest unextended primers afterRT',S.DIGESTION,((washed_cdna,R.TEMPLATE),(exo,R.CATALYST)),(digested,))
    clean_cdna=clean(g,'post_exo_washes',digested,'Washed cDNA-bearing beads afterExoI','TE-SDS/TE-TW/water bead washes')
    clean_cdna=aliquot(g,'pcr_beads',clean_cdna,'Selected actual2000bead PCRcontribution; other sourcePCRtubes/stock unexpanded')
    smart=g.input('smart','One SMART PCR primer species',K.OLIGO,state=Construct([Segment('SMART PCR primer','AAGCAGTGGTATCAACGCAGAGT')]))
    pcr_mix=g.input('pcr_mix','Actual source25µl2xKAPAHiFiReadyMix for50µl/2000beadPCR; otherrecipe24.6µlwater+0.4µl100µMSMARTprimer',K.MOLECULAR)
    polymerase=g.product('polymerase','Copresent KAPAHiFi polymerase; source4cycles65C then11cycles67C',K.CATALYST)
    pcr_nt=g.product('pcr_nucleotides','Copresent sourceReadyMix nucleotide class; no separatefreshdose',K.MOLECULAR)
    pcr_buffer=g.product('pcr_buffer','Copresent sourceReadyMix buffer class',K.CONDITION)
    g.react('pcr_mix_components','Components of one actual source25µlReadyMix addition',S.MAGIC,((pcr_mix,R.TEMPLATE),),(polymerase,pcr_nt,pcr_buffer),changes_molecule=False,reason='Source explicitly names2xKAPAHiFiHotstartReadyMix; one actual stock addition decomposed into functional classes.')
    amplicons=g.product('cdna_amplicons','Amplified source DART-seq cDNA',K.MOLECULAR)
    g.react('cdna_pcr','Single-species SMART-handle cDNA amplification',S.PCR,
        ((clean_cdna,R.TEMPLATE),(smart,R.STRUCTURAL),(polymerase,R.CATALYST),(pcr_nt,R.STRUCTURAL),(pcr_buffer,R.CONDITION)),(amplicons,),
        geometry=PCRGeometry.unavailable(unknown_template('full cDNA handle architecture'),symmetric_species=True,
            reason='Published protocol specifies one SMART primer species; full source cDNA inward-site architecture is deferred, not assumed from matching platform.'))
    amplicons=clean(g,'cdna_cleanup',amplicons,'Purified amplified cDNA','0.6×SPRI AMPureXP')
    amplicons=qc(g,'cdna_qubit',amplicons,'SourceQubitDNAHS cDNA measurement; actualaliquot amount unassigned')
    amplicons=qc(g,'cdna_bioanalyzer',amplicons,'SourceBioAnalyzer cDNA QC; actualaliquot amount unassigned')
    amplicons=aliquot(g,'nextera_input',amplicons,'Actual source600pg cDNA forNexteraXT; transfer volume unassigned')
    tagstock=g.input('tagmentation_inventory','Aggregate actual sourceNexteraXTFC-131-1096 tagging reagent inventory; packaging/stockcounts unassigned',K.MOLECULAR)
    tn=g.product('transposase','Source tagmentation catalyst component; exact loadedcomplex/formulation unassembled',K.CATALYST)
    adapters=g.product('loaded_adapters','Same-source loaded adapter components; whole original ends unassembled, not freshly added free adapter stock',K.ADAPTER)
    g.react('tagmentation_components','Source tagging inventory components, not invented freshfreeadapter addition',S.MAGIC,((tagstock,R.TEMPLATE),),(tn,adapters),changes_molecule=False,reason='Source specifiesNexteraXTfragment/tag; exacthistorical kit stockformulation and original loadedstrands unimplemented.')
    tagged=g.product('tagged','Tagmented amplified cDNA',K.MOLECULAR)
    g.react('tagment','Nextera cDNA tagmentation',S.TAGMENTATION,((amplicons,R.TEMPLATE),(tn,R.CATALYST),(adapters,R.STRUCTURAL)),(tagged,))
    indexed=pcr(g,'index_pcr',tagged,'Source custom P5-PCRhybrid primer with terminalT*A*C','Source indexPCR primer2','Indexed DART-seq products','Source Nextera libraryPCR')
    lib=clean(g,'library_cleanup',indexed,'Targeted-capture DART-seq library','Source library purification',K.LIBRARY)
    lib=qc(g,'final_qubit',lib,'SourceQubitDNAHS final library measurement; aliquot unassigned')
    lib=qc(g,'final_bioanalyzer',lib,'SourceBioAnalyzer final library QC; aliquot unassigned')
    lib=aliquot(g,'loading',lib,'Selected actualNextSeq500 library contribution; loading/pooling amounts unassigned')
    g.sequence(lib,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':20,'Read 2':130},limitation='Main Methods and published detailed protocol specify NextSeq500,Read1=20,Read2=130; detailed protocol also specifies customRead1 primer. Primer landing and historical full-bead barcode/UMI widths remain deferred.'))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    g.declare_reviewed_scope(ADDITIONAL_REVIEWED,SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Full source primer-pool, intermediate and final-strand assembly',LIMIT+' Source material is available for follow-up; not an absence finding.'))
    return g

def detail():
    sc=Scene();sc.strand('source_bead',list(source_bead()),label='PublishedMOESM3 bead classes; sourceJ12/N8 actual members unassigned')
    return panel(sc.rows(),caption='Full source-provided bead baseclasses; T7,handle,J12cellbarcode,N8UMI,T30. Exact beadcoupling chemistry not assigned.')+panel(ligation_scene().rows(),caption='Source IgG splinted-ligation interface (partial terminal window)')+panel(sealed_window_scene().rows(),caption='Selected sealed IgG component; other actual pool members not yet assembled')+info('Source customRead1: '+CUSTOM_READ1+'; source P5hybrid: '+P5_HYBRID[:-2]+'*A*C. TSO DNAprefix '+TSO_DNA+' followed by rGrGrG. Fullfinalstrand and actual readlanding remain unimplemented; no unseen historicalNextera index/adapter assumed.')

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
        info('PRIVATE MODEL-SCOPE HOLD — required fullsourcepool/intermediate/finalstrandassembly remains unimplemented. '+SOURCE),info(LIMIT),detail(),''.join(workflow_panel(workflow_from_protocol_graph(protocol_graph(),endpoint)) for endpoint in ['loading_selected']),'</div>'])

def selftest():
    from checks import Check
    c=Check(); g=protocol_graph()
    issues=g.publication_problems()
    c('Only explicit unassembled-source coverage blocks publication',issues,["authoritative source item 'unassembled_source' is not accounted for"])
    c('Private preview does not claim protocol completion','data-protocol-complete="true"' in render_page(),False)
    extra_checks(c,g)
    ports=lambda key:tuple(p.material for p in g.steps[key].inputs)
    c('actual arraypool combines distinct source reactions',ports('bead_product_pool'),('heated','other_ligation_products'))
    c('other reactions use allocated other bead contribution','bead_others' in ports('other_ligations'),True)
    c('other reactions use allocated probe contribution','other_probe' in ports('other_ligations'),True)
    c('no RTcatalyst in droplet capture',any(p.role is R.CATALYST for p in g.steps['droplet_capture'].inputs),False)
    c('RT uses inventory from actual postbreak wash',ports('rt_bead_components'),('pre_rt_washes_kept',))
    c('TSO and enzyme derive from one actual200ulmix',ports('rt_mix_components'),('rt_mix',))
    c('no fresh loaded adapter supplied as independent dose',ports('tagmentation_components'),('tagmentation_inventory',))
    c('known custom P5 hybrid actually used',g.materials['index_pcr_forward'].state.top(),P5_HYBRID)
    c('known RTprimers preserve cell identity',source_bead().get('sourceJ12 cellbarcode').feature.role,'cell_barcode')
    c('full source cell feature survives bead ligation',g.materials['joined'].chemical_state['selected_full_source_bead'].get('sourceJ12 cellbarcode').feature,source_bead().get('sourceJ12 cellbarcode').feature)
    c('known N8 source semanticUMI',source_bead().get('sourceN8 UMI').feature.role,'umi')
    c('source selected600pg remainder retained',g.terminals['nextera_input_remainder'].disposition,Disposition.UNRESOLVED)
    c('source actual finalQC typed',g.materials['final_bioanalyzer_qc'].kind,K.QC)
    c('required fullassembly not excluded','unassembled_source' in g.source_coverage,False)
    from checks import Source,have
    from storage import sources_root
    source=Source(sources_root()/Path(__file__).resolve().parents[1].name/'s41592-018-0259-9_supp_NIHMS1512056-supplement-2.pdf.txt','ExactpublishedMOESM3protocolp3')
    if have(c,source,label='DART exactpublished protocol orderings'):
        text=''.join(source.path.read_text().split())
        for name,body in [('bead',BEAD_WRITTEN),('TSO',TSO_DNA+'rGrGrG'),('customP5PS',P5_HYBRID[:-2]+'*A*C'),('customRead1',CUSTOM_READ1)]:
            c('whole original protocol '+name,body in text,True)
    c.report(exit_on_fail=True)

BEAD_HANDLE='AAGCAGTGGTATCAACGCAGAGTAC'
BEAD_WRITTEN='TTTTTTT'+BEAD_HANDLE+'J'*12+'N'*8+'T(30)'
TSO_DNA='AAGCAGTGGTATCAACGCAGAGTGAAT'
CUSTOM_READ1='GCCTGTCCGCGG'+BEAD_HANDLE
P5_HYBRID=il.P5+CUSTOM_READ1
def source_bead(*,extended=False):
    parts=[Segment('sourceT7','T'*7),Segment('source PCRhandle',BEAD_HANDLE),Segment('sourceJ12 cellbarcode','N'*12,placeholder=True,feature=feature('dart.cell','cell_barcode','unknown',note='Source J12split-pool class; actual member unassigned')),Segment('sourceN8 UMI','N'*8,placeholder=True,feature=feature('dart.umi','umi','random')),Segment('sourceT30','T'*30)]
    if extended:parts.append(Segment('selected source IgGextension',IGG))
    return Construct(parts,name='Publishedprotocol sourcebead class'+(' withselectedIgGcomponent' if extended else ''))
def qc(g,key,stock,label):
    measured=g.product(key+'_qc',label,K.QC)
    remaining=g.product(key+'_remaining','Remaining source stock after actual QC aliquot',stock.kind,chemical_state=stock.chemical_state)
    g.react(key,label,S.SPLIT,((stock,R.TEMPLATE),),(measured,remaining))
    g.mark_terminal(Terminal(measured.id,Disposition.QC,'Actualsource QC contribution, not sequencing material.'))
    return remaining
ADDITIONAL_REVIEWED=('bead_array_allocation','other_toeholds','toehold_pooling','probe_array_allocation','ligation_conditions','other_ligation_supplies','other_ligations','bead_product_pool','droplet_conditions','droplet_capture','emulsion_break_conditions','break_emulsion','pre_rt_washes_selection','pre_rt_washes','rt_bead_components','rt_mix','rt_mix_components','pcr_beads','pcr_nucleotides','pcr_mix','pcr_mix_components','cdna_qubit','cdna_bioanalyzer','nextera_input','tagmentation_inventory','tagmentation_components','index_pcr_nucleotides','final_qubit','final_bioanalyzer','loading','loading_selected')
