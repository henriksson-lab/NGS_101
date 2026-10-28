"""DeKosky 2013 linkage chemistry, with source-limited sequencing handoff."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, complement_segments, revcomp, panel, protocol_graph_panels
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S,
    SequencingPlan, Terminal, Disposition)
from pcr import PCRGeometry
from page import head, info

PAGE_FILE='paired_vh_vl.html'
SOURCE='https://doi.org/10.1038/nbt.2492; Methods memory B-cell pairing; published Supplementary Fig.1/Tables5-6'
H_TAIL='TATTCCCATGGCGCGCC'
L_TAIL='GGCGCGCCATGGGAATA'
PRIMER_EXAMPLES={
    'CHrev-AHX89':'CGCAGTAGCGGTAAACGGC',
    'CLrev-BRH06':'GCGGATAACAATTTCACACAGG',
    'hVH4-DP63-fwd-OE':'TATTCCCATGGCGCGCCCAGGTGCAGCTACAGCAGTGGG',
    'hVK5-fwd-OE':'GGCGCGCCATGGGAATAGCCGAAACGACACTCACGCAGTCTC',
}
NESTED={
    'hIgM-rev-OEnested':'ATGGGCCCTGGGTTGGGGCGGATGCACTCC',
    'hIgKC-rev-OEnested':'GTGCGGCCGCAGATGGTGCAGCCACAGTTC',
}

@dataclass(frozen=True)
class LinkageReference:
    """Checks the published primer-tail relation, not transcript or pairing fidelity."""
    heavy_tail: str = H_TAIL
    light_tail: str = L_TAIL

    def __post_init__(self):
        if revcomp(self.heavy_tail)!=self.light_tail:
            raise ValueError('Heavy/light primer tails must be reverse complements')

    def heavy(self):
        return Construct([
            Segment('N-terminal CH1','X'*12,placeholder=True),
            Segment('VH reverse-strand region','X'*28,placeholder=True),
            Segment('heavy tail complement',revcomp(self.heavy_tail))],name='Heavy active strand')

    def light(self):
        return Construct([
            Segment('light shared tail',self.light_tail),
            Segment('VL region','X'*28,placeholder=True),
            Segment('N-terminal light constant region','X'*12,placeholder=True)],name='Light reference strand')

    def product(self):
        left,right=self.heavy(),self.light()
        if left.segments[-1].top!=right.segments[0].top:
            raise ValueError('The declared fragments must share the same overlap once')
        return Construct([*left.segments,*right.segments[1:]],
            name='VH:VL linked reference (variable regions unassigned)')

    def annealed(self):
        left,right=self.heavy(),self.light()
        sc=Scene();sc.strand('heavy',list(left),label='heavy active strand')
        sc.anneal('light',complement_segments(list(right)),to='heavy',
            pair=("light shared tail'",'heavy tail complement'),label='light opposite strand')
        sc.note('heavy','Only the primer-derived overlap is exact; gene lengths, transcripts and amplification specificity are unassigned')
        return sc

LINK=LinkageReference()

def primer_pool():
    sc=Scene()
    for i,(name,seq) in enumerate(PRIMER_EXAMPLES.items()):
        sc.strand('primer'+str(i),[Segment(name,seq)],label=name)
    return sc

def protocol_graph():
    g=ProtocolGraph('Paired VH:VL — linked-repertoire branch',authority=SOURCE)
    cells=g.input('cells','B cells already deposited into 125-pL microwells',K.BIOLOGICAL,primary=True)
    lysis=g.input('lysis_buffer','LiCl/Tris +1% lithium dodecyl sulfate, EDTA and DTT; dialysis-covered wells',K.CONDITION)
    lysate=g.product('lysate','Lysed cells within sealed wells',K.MOLECULAR)
    g.react('lysis','Lyse cells under dialysis membrane',S.PERMEABILIZATION,
        ((cells,R.TEMPLATE),(lysis,R.CONDITION)),(lysate,))
    beads=g.input('poly_dt','Poly(dT) magnetic beads; capture in wells, recover and wash before RT-PCR',K.SELECTION)
    captured=g.product('captured','Bead-bound heavy/light mRNAs',K.MOLECULAR)
    unbound=g.product('unbound','Unretained lysate and wash material',K.WASTE)
    g.react('capture','Capture polyadenylated mRNA and magnetically recover/ wash beads',S.AFFINITY_CAPTURE,
        ((lysate,R.TEMPLATE),(beads,R.SELECTION)),(captured,unbound))
    g.discard(unbound,'Removed during bead recovery and washes before emulsification')
    primers=g.input('linkage_primers','Complete Table5 constant/variable primer pools; displayed exact representatives are not the whole pool',K.OLIGO,state=primer_pool())
    mix=g.input('one_step_mix','Quanta OneStep Fast RT-PCR mix; bead suspension emulsified in mineral oil/surfactants',K.CATALYST)
    first=g.product('first_strands','Heavy/light first-strand cDNA in bead-containing emulsion compartments',K.MOLECULAR)
    g.react('rt','Reverse transcribe at 55°C for 30 min, then 94°C for 2 min',S.REVERSE_TRANSCRIPTION,
        ((captured,R.TEMPLATE),(primers,R.STRUCTURAL),(mix,R.CATALYST)),(first,))
    ch=g.input('ch_amplifier','CHrev-AHX89 amplification primer',K.OLIGO,state=Construct([Segment('CHrev-AHX89',PRIMER_EXAMPLES['CHrev-AHX89'])]))
    cl=g.input('cl_amplifier','CLrev-BRH06 amplification primer',K.OLIGO,state=Construct([Segment('CLrev-BRH06',PRIMER_EXAMPLES['CLrev-BRH06'])]))
    linked=g.product('linked','Desired linked VH:VL duplex (~850 bp population)',K.MOLECULAR,state=LINK.product(),chemical_state=LINK)
    byproduct=g.product('byproduct','Complementary-strand by-product drawn in published Supplementary Fig.1c',K.MOLECULAR)
    g.react('linkage_pcr','Incorporate V-primer tails, anneal and extend through shared overlap in emulsion PCR',S.PCR,
        ((first,R.TEMPLATE),(primers,R.STRUCTURAL),(ch,R.STRUCTURAL),(cl,R.STRUCTURAL),(mix,R.CATALYST)),(linked,byproduct),
        geometry=PCRGeometry.unavailable(LINK.annealed(),reason='Primer tails are explicit, but individual V(D)J templates and full pool landing coordinates are unassigned. The checked tail complement alone does not validate whole amplicons.'))
    g.mark_terminal(Terminal(byproduct.id,Disposition.UNRESOLVED,
        'Supplementary Fig.1c explicitly depicts this product class. Its later molecular fate was not resolved in the capped reading; no absent-information claim is made.'))
    cleanup=g.input('cleanup','Ether emulsion extraction and silica spin-column recovery',K.SELECTION)
    recovered=g.product('recovered','Recovered linked cDNA',K.MOLECULAR)
    waste=g.product('emulsion_waste','Discarded oil/ether phases and unretained column material',K.WASTE)
    g.react('recover','Break emulsion, remove oil and recover cDNA',S.CLEANUP,
        ((linked,R.TEMPLATE),(cleanup,R.SELECTION)),(recovered,waste))
    g.discard(waste,'Oil/ether upper phases and column wash material discarded during recovery')
    heavy_nested=g.input('heavy_nested','Table6 nested heavy constant-region primer pool',K.OLIGO)
    light_nested=g.input('light_nested','Table6 nested light constant-region primer pool',K.OLIGO)
    taq=g.input('taq','Taq polymerase in ThermoPol buffer',K.CATALYST)
    amplified=g.product('amplified','Nested amplified VH:VL products',K.MOLECULAR)
    g.react('nested_pcr','Nested PCR (39 cycles in the memory B-cell branch)',S.PCR,
        ((recovered,R.TEMPLATE),(heavy_nested,R.STRUCTURAL),(light_nested,R.STRUCTURAL),(taq,R.CATALYST)),(amplified,),
        geometry=PCRGeometry.unavailable(LINK.product(),reason='Nested primer sequences are published, but exact chain templates and primer-to-template coordinates are not assigned in this representation.'))
    gel=g.input('gel','Agarose gel electrophoresis and band extraction',K.SELECTION)
    selected=g.product('selected','Gel-extracted approximately 850-bp linked product',K.MOLECULAR)
    rejected=g.product('gel_rejected','Other DNA bands and unrecovered gel material',K.WASTE)
    g.react('gel_select','Gel-purify linked-size product',S.CLEANUP,
        ((amplified,R.TEMPLATE),(gel,R.SELECTION)),(selected,rejected))
    g.discard(rejected,'Excluded by band selection; no assertion that the Fig.1c class is quantitatively removed')
    lib=g.product('library','Linked product submitted for MiSeq sequencing',K.LIBRARY)
    g.react('sequencing_handoff','Source-limited library preparation handoff',S.MAGIC,
        ((selected,R.TEMPLATE),),(lib,),reason='The inspected workflow ends at gel extraction and submission for MiSeq. Exact sequencing-adapter construction/primer geometry remains unreviewed; no ligation or kit recipe is invented.')
    g.sequence(lib,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':250,'Read 2':250},
        limitation='2×250 paired-end MiSeq interrogates CDR-H3/CDR-L3 at opposite outer portions. Sequencing adapters/priming sites and exact cycle coordinates remain unassigned in this pass.'))
    g.declare_reviewed_scope(('cells','lysis_buffer','lysate','lysis','poly_dt','captured','unbound','capture',
        'ch_amplifier','cl_amplifier','linkage_primers','one_step_mix','first_strands','rt','linked','byproduct','linkage_pcr',
        'cleanup','recovered','emulsion_waste','recover','heavy_nested','light_nested','taq','amplified',
        'nested_pcr','gel','selected','gel_rejected','gel_select','library','sequencing_handoff'),SOURCE)
    return g

def render_page():
    return '\n'.join([head('Paired VH:VL — bounded linkage model'),'<div class="wrap"><h1>Paired VH:VL</h1>',
        '<p class="research-notes"><a href="01_paired_vh_vl.html">Research notes</a></p>',
        info('Physical heavy/light-chain linkage preserves a receptor pair in one amplicon. Only the primer-derived overlap is exact here; chain templates, amplification specificity, by-product fate and sequencing-library handoff remain unresolved.'),
        protocol_graph_panels(protocol_graph(),cls='long'),
        panel(LINK.annealed().rows(),caption='Source-derived complementary primer tails; unknown receptor bodies are schematic and not invented sequence.'),'</div>'])

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Source overlap tails are complementary',revcomp(H_TAIL),L_TAIL)
    c('Linked reference retains one overlap',LINK.product().top().count(L_TAIL),1)
    c('Source-drawn by-product remains explicitly unresolved',g.terminals['byproduct'].disposition,Disposition.UNRESOLVED)
    c('No arbitrary library ligation inferred',g.steps['sequencing_handoff'].kind,S.MAGIC)
    from storage import sources_root
    source=sources_root()/Path(__file__).resolve().parents[1].name/'nbt.2492_41587_2013_BFnbt2492_MOESM35_ESM.pdf.txt'
    if source.exists():
        text=''.join(source.read_text().split())
        for name,seq in {**PRIMER_EXAMPLES,**NESTED}.items():c('Published primer transcription '+name,seq in text)
    else:c.skip('Primer transcription','External published supplement unavailable')
    c.report()
