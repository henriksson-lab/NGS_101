"""Original RAD-seq source-owned components and operational graph."""
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

PAGE_FILE = 'rad-seq.html'
SOURCE = 'Baird et al.2008 DOI10.1371/journal.pone.0003376;published PMC2557064 Results72/Fig1 and Methods141/145.'
P1_HEAD = il.NEBNEXT_UNIVERSAL_PRIMER
P2_A = 'CTCAGGCATCACTCGATTCCTCCGAGAACAA'
P2_B = il.P7[:21] + 'CGGAGGAATCGAGTGATGCCTGAGT'
FORWARD = il.P5[:20]
REVERSE = il.P7[:21]
INDEX = feature('rad.sample', 'sample_index', 'unknown', note='Source designed4/5nt sample barcode;member not assigned,not a random UMI.')
SOURCE_OLIGO_EXCEPTIONS = {'p1': 'Published PMC2557064 Methods141 prints lowercase xxxx for the SbfI4nt barcode (and xxxxx for EcoRI5nt); complete flanking oligos and5prime phosphate read in the original. Model imports exact canonical58nt body and its computed reverse complement, retaining the source barcode width as a structured unassigned sample index. Existing strict text matcher does not join lowercase x runs; no sequence-absence claim or synthetic source twin. Independent reviewer compared all four original P1 strands and the P2/PCR oligos to model components;whole-source transcription matches. Not an absence audit.'}


def p1_scene(enzyme='SbfI'):
    if enzyme not in {'SbfI', 'EcoRI'}:
        raise ValueError('Only the two source adapter variants are declared')
    width = 4 if enzyme == 'SbfI' else 5
    index = lambda: Segment('sample barcode', 'N' * width, placeholder=True, feature=INDEX)
    top = [Segment('source P1 body', P1_HEAD), index()]
    bottom = [index(), Segment('P1 complementary body', revcomp(P1_HEAD))]
    if enzyme == 'SbfI': top.append(Segment('source3prime cohesive end', 'TGCA'))
    else: bottom.insert(0, Segment('source5prime cohesive end', 'AATT'))
    sc = Scene(); sc.strand('top', top, label=enzyme+' P1 top')
    sc.anneal('bottom', bottom, to='top', pair=('P1 complementary body','source P1 body'), mod5='P', label='Source phosphorylated bottom')
    strand,end,polarity=('top','source3prime cohesive end','3-prime') if enzyme=='SbfI' else ('bottom','source5prime cohesive end','5-prime')
    sc.ligation_end(strand,end,kind='cohesive',polarity=polarity,owner='adapter')
    return sc


def p2_ordering():
    # Separate source ordering strands, not a falsely paired full Y structure.
    sc = Scene()
    sc.strand('a', [Segment('source P2 top', P2_A)], mod5='P', label='P2 top5prime phosphorylated')
    sc.strand('b', [Segment('source P2 bottom body', P2_B[:-1]),Segment('terminal3prime T',P2_B[-1])], label='P2 bottom;source divergent arm')
    sc.ligation_end('b','terminal3prime T',kind='cohesive',polarity='3-prime',owner='adapter')
    return sc


def p2_binding_reference():
    # Derived complementarity window only; no full Y folding or ligation claim.
    matches = [P2_A[:n] for n in range(1,len(P2_A)+1) if revcomp(P2_A[:n]) in P2_B]
    stem = max(matches, key=len)
    sc = Scene();sc.strand('a',[Segment('derived common window',stem)],label='P2 top component window')
    sc.anneal('b',[Segment('complementary window',revcomp(stem))],to='a',pair=('complementary window','derived common window'),label='P2 bottom component window')
    return sc


def protocol_graph():
    g = ProtocolGraph('RAD-seq — original SbfI library route', authority=SOURCE, representation='published_instance')
    def op(key,label,kind,material,reagent_label,*,molecular=False,additional=()):
        is_condition=kind in {S.THERMAL_TREATMENT,S.FRAGMENTATION}
        reagent=g.input(key+'_reagent',reagent_label,K.MOLECULAR if molecular else K.CONDITION if is_condition else K.CATALYST)
        product=g.product(key+'_mixture',label+';original and newly generated material/reaction remnants not quantitatively resolved',K.MOLECULAR)
        kwargs={'reason':'Source named preparation; exact historical kit components and molecular end products unreviewed.'} if kind is S.MAGIC else {}
        g.react(key,label,kind,((material,R.TEMPLATE),(reagent,R.STRUCTURAL if molecular else R.CONDITION if is_condition else R.CATALYST))+tuple(additional),(product,),**kwargs);return product
    def cleanup(key,material,label):
        procedure=g.input(key+'_procedure',label,K.SELECTION)
        keep=g.product(key,'Source recovered preparation;purity/recovery unquantified',K.MOLECULAR)
        reject=g.product(key+'_reject','Physically excluded fractions;composition unassigned',K.WASTE)
        g.react(key+'_step',label,S.CLEANUP,((material,R.TEMPLATE),(procedure,R.SELECTION)),(keep,reject))
        g.discard(reject,'Source actual physical separation,not perfect depletion.');return keep
    def ligate(key,material,adapter,label):
        published=key=='p1_ligate'
        enzyme=g.input(key+'_ligase','Actual source1000U T4DNA ligase' if published else 'INFERRED — functional ligase for reported P2 joining; actual identity/dose unreviewed',K.CATALYST,inferred=not published)
        additions=()
        if published:
            atp=g.input('p1_atp','Actual source1ul100mM rATP contribution',K.MOLECULAR)
            buffer=g.input('p1_supplement','Actual1ul10xEcoRIbuffer and5ulwater; separate source supply allocations, packagingunassigned',K.MOLECULAR)
            additions=((atp,R.STRUCTURAL),(buffer,R.STRUCTURAL))
        product=g.product(key+'_mixture','Selected source ligated DNA plus free adapter/ligase/possible unligated DNA',K.MOLECULAR)
        geometry=AdapterLigationGeometry.known('cohesive',adapter.state,molecule('Source substrate boundary',(role('Full genomic fragment ends and covalent product junction unassembled'),)),statement='Source SbfI P1 has3prime TGCA;source P2 bottom has terminal3prime T beyond the derived complementary stem. Complete genomic substrate and sealed junction geometry remain unassembled. Separate P2 ordering strands do not assert a complete paired Y fold.')
        g.react(key,label if published else 'INFERRED — functional ligase inventory for source P2 joining atRT',S.LIGATION,((material,R.TEMPLATE),(adapter,R.STRUCTURAL),(enzyme,R.CATALYST))+additions,(product,),geometry=geometry);return product
    dna=g.input('dna','Source genomic DNA;sample origins/quantity not treated as separate protocol entries',K.MOLECULAR,primary=True)
    dna=op('digest','Source SbfI digest37°C15min;EcoRI is a separate source adapter branch',S.DIGESTION,dna,'Source SbfI;actual genomic fragment sequences not assigned')
    dna=op('digest_inactivation','Source65°C20min heat inactivation',S.THERMAL_TREATMENT,dna,'Source heat treatment;not physical enzyme removal')
    p1=g.input('p1','Source SbfI P1 adapter,structured4nt sample index and phosphorylated bottom',K.ADAPTER,state=p1_scene())
    dna=ligate('p1_ligate',dna,p1,'Source P1 ligation20min before pooling/shearing')
    dna=op('ligation_inactivation','Source65°C20min heat inactivation;free adapters/remnants persist',S.THERMAL_TREATMENT,dna,'Source heat treatment')
    additional=g.input('other_samples','Other independently digested/barcoded P1-ligated samples;exact members unassigned',K.MOLECULAR)
    pooled=g.product('pool','Source multiplexed preparations;sample identity remains in P1 barcode',K.MOLECULAR)
    g.react('sample_pool','Source physical pooling after P1 tagging',S.POOL,((dna,R.TEMPLATE),(additional,R.TEMPLATE)),(pooled,))
    dna=op('shear','Random sonication average500bp;includes fragments with and without P1',S.FRAGMENTATION,pooled,'Source Bioruptor or Branson sonicator;cleavage positions/end chemistry unassembled')
    dna=cleanup('first_gel',dna,'Source1%agarose300–700bp excision/MinElute extraction')
    dna=op('blunt','Source Quick Blunting Kit end polishing',S.MAGIC,dna,'Actual historical NEB Quick Blunting molecular supplies; internal components and full ends unassembled',molecular=True)
    dna=cleanup('blunt_cleanup',dna,'Source Quick Spin column')
    adenine=g.input('adenine_supply','Source-written adenine(Fermentas) molecular contribution; exact nucleotide form/formulation unreviewed',K.MOLECULAR)
    dna=op('a_tail','Source Klenow exo-minus/adenine37°C adenylation',S.TAILING,dna,'Source Klenow exo-minus15U;duration unassigned',additional=((adenine,R.STRUCTURAL),))
    dna=cleanup('tail_cleanup',dna,'Source purification after adenylation;precise cartridge not reassigned')
    p2=g.input('p2','Source divergent P2 ordering strands;complete Y geometry not assembled',K.ADAPTER,state=p2_ordering())
    dna=ligate('p2_ligate',dna,p2,'Source P2 ligation after random shearing/end polishing/A-tailing')
    dna=cleanup('p2_cleanup',dna,'Source purification and50µl elution')
    aliquot=g.product('pcr_input','Source5µl aliquot for preparative PCR',K.MOLECULAR)
    remainder=g.product('pcr_remainder','Independent remainder of source50µl elution;later fate unassigned',K.MOLECULAR)
    g.react('pcr_split','Source5µl PCR allocation with unallocated physical remainder',S.SPLIT,((dna,R.TEMPLATE),),(aliquot,remainder))
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'Source describes PCR aliquot;remaining preparation fate not assigned in model.'))
    primers=g.input('pcr_primer_mix','Actual source5ul10uM modifiedSolexa amplification-primer mix',K.MOLECULAR)
    f=g.product('forward','P1-forward20nt component of sameactualprimer mix',K.OLIGO,state=Construct([Segment('P1-forward',FORWARD)]))
    r=g.product('reverse','P2-reverse21nt component of sameactualprimer mix',K.OLIGO,state=Construct([Segment('P2-reverse',REVERSE)]))
    g.react('primer_components','Expose actual copresent PCR-primer mixture components; not extra doses',S.GENERIC,((primers,R.STRUCTURAL),),(f,r),changes_molecule=False)
    mix=g.input('pcr_mastermix','Actual source50ul Phusion Master Mix; exact historical formulation unreviewed',K.MOLECULAR)
    poly=g.product('phusion','Polymerase functional component of same actual master mix',K.CATALYST)
    nts=g.product('pcr_nts','Functional nucleotide component of same actual master mix; precise recipe unreviewed',K.MOLECULAR)
    other=g.product('pcr_other','Other copresent master-mix components; exact composition unassembled',K.MOLECULAR)
    g.react('pcr_components','Expose actual source master-mix components; no fresh independent doses',S.GENERIC,((mix,R.STRUCTURAL),),(poly,nts,other),changes_molecule=False)
    water=g.input('pcr_water','Actual source40ulwater added to PCR',K.MOLECULAR)
    copied=g.product('pcr','New P1-initiated PCR copies plus original templates,unligated/non-P1 background and reaction remnants',K.MOLECULAR)
    g.react('amplify','P1-initiated copying creates complementary P2 reverse-primer site before selective amplification;not guaranteed background elimination',S.PCR,((aliquot,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(poly,R.CATALYST),(nts,R.STRUCTURAL),(other,R.STRUCTURAL),(water,R.STRUCTURAL)),(copied,),geometry=PCRGeometry.unavailable(molecule('Source PCR boundary',(role('Full sheared-template/end/first-copy products unassembled'),)),reason='Known intended Y-adapter gate is distinct from a checked whole molecular PCR relation;source background/yield not predicted.'))
    dna=cleanup('final_gel',copied,'Source final300–700bp gel purification and10nM dilution')
    lib=g.product('library','Source original RAD library;full acquired primer/end/run geometry unassembled',K.LIBRARY)
    g.react('dilute','Source final dilution for acquisition',S.GENERIC,((dna,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Source Illumina Genome Analyzer/Solexa protocols;actual acquired cycles and complete primer geometry not yet assigned. Barcode+26nt computational trimming is not a physical run length.'))
    g.declare_reviewed_scope(('dna','digest_reagent','digest','digest_inactivation_reagent','digest_inactivation','p1','p1_atp','p1_supplement','p1_ligate_ligase','p1_ligate','ligation_inactivation_reagent','ligation_inactivation','other_samples','sample_pool','shear_reagent','shear','first_gel_procedure','first_gel_step','blunt_reagent','blunt','blunt_cleanup_procedure','blunt_cleanup_step','adenine_supply','a_tail_reagent','a_tail','tail_cleanup_procedure','tail_cleanup_step','p2','p2_ligate','p2_cleanup_procedure','p2_cleanup_step','pcr_split','pcr_primer_mix','primer_components','pcr_mastermix','pcr_components','pcr_water','amplify','final_gel_procedure','final_gel_step','dilute'),SOURCE+' Literal main141 operational scope; whole Y-product/kit/readout obligations unfulfilled.')
    g.declare_reviewed_scope(('p2_ligate_ligase',),SOURCE+' Inferred functional ligase requirement for reported P2 joining; no historical enzyme identity/dose asserted.')
    g.add_source_item(SourceItem('historical_readout','Historical Solexa/Phusion protocols and complete source library geometry','Main141 refers to historical product guidelines;not replaced by a modern kit. Source Fig1 references Coyne2004 for Y-adapter design;direct predecessor original unreviewed.'))
    return g


def render_page():
    return '\n'.join([head('RAD-seq'),'<main class="wrap"><h1>RAD-seq</h1>','<p><a href="01_rad_seq.html">Research notes</a></p>',info(SOURCE),info('Private incomplete source model. INFERRED — P2 functional ligase is required for reported joining; its actual identity/dose remains unreviewed. Original restriction tagging then random shearing differs from ddRAD. Exact source adapter sequences are separate from unassembled full genomic junctions,selected copying and historical readout. No absence or full-validation claim.'),panel(p1_scene().rows(),caption='Source SbfI P1 adapter component'),panel(p2_ordering().rows(),caption='Separate P2 ordering strands'),panel(p2_binding_reference().rows(),caption='Derived P2 complementarity window only'),workflow_panel(workflow_from_protocol_graph(protocol_graph(),'library')),'</main>'])


def selftest():
    from checks import Check, Source, have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('P1 bottom phosphate retained',p1_scene().strands['bottom'].mod5,'P')
    c('Sample identifier is not UMI',p1_scene().strands['top'].segs[1].feature.role,'sample_index')
    c('Exact source P2 reverse primer is truncated canonical P7',P2_B.startswith(REVERSE),True)
    c('Random fragmentation follows physical sample pooling',g.steps['shear'].inputs[0].material,'pool')
    c('Preparative PCR uses source aliquot not remainder',g.steps['amplify'].inputs[0].material,'pcr_input')
    c('Physical PCR remainder has explicit unresolved fate',g.terminals['pcr_remainder'].disposition,Disposition.UNRESOLVED)
    c('Both PCR primers derive from same actual source mixture',tuple(g.producer_of(k).id for k in ('forward','reverse')),('primer_components','primer_components'))
    c('Full reaction geometry remains incomplete',bool(g.publication_problems()))
    c('Private preview has no finished marker','data-protocol-complete=' not in render_page())
    original=Source(sources_root() / 'rad-seq__10.1371+journal.pone.0003376' / 'journal.pone.0003376_PMC2557064.xml.txt','Published PMC2557064 Methods141;standard get_sources.py')
    if have(c,original,label='Whole source P1 transcription with lowercase barcode placeholders'):
        text=original.path.read_text()
        for enzyme in ('SbfI','EcoRI'):
            written=[]
            for strand in p1_scene(enzyme).strands.values():
                body=''.join(seg.top.replace('N','x') if seg.placeholder else seg.top for seg in strand.segs)
                written.append('5′-'+('Phos-' if strand.mod5=='P' else '')+body+'-3′')
            c('Whole printed '+enzyme+' P1 top/bottom and phosphate',tuple(value in text for value in written),(True,True))
        for key,label in (('forward','P1-forward primer'),('reverse','P2-reverse primer')):
            body=g.materials[key].state.top()
            c('Produced '+label+' matches original whole ordering',label+': 5′-'+body+'-3′' in text,True)
    c.report()
