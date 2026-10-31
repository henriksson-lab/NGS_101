"""DeKosky 2013 linkage chemistry, with source-limited sequencing handoff."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, complement_segments, revcomp, panel, protocol_graph_panels
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S,
    SequencingPlan, Terminal, Disposition)
from pcr import PCRGeometry
from page import head, info, caveat, table
from protocol_views import assemble_page, WorkflowSection

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
            Segment('N-terminal CH1','[constant context]',placeholder=True,inferred=True),
            Segment('VH reverse-strand region','[variable context]',placeholder=True,inferred=True),
            Segment('heavy tail complement',revcomp(self.heavy_tail),inferred=True)],name='Heavy active strand')

    def light(self):
        return Construct([
            Segment('light shared tail',self.light_tail,inferred=True),
            Segment('VL region','[variable context]',placeholder=True,inferred=True),
            Segment('N-terminal light constant region','[constant context]',placeholder=True,inferred=True)],name='Light reference strand')

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
    g=ProtocolGraph('Paired VH:VL — linked-repertoire branch',authority=SOURCE,representation='published_instance')
    cells=g.input('cells','Actual source-prepared memory-B/IM9-control cohort already deposited into 125-pL microwells; no perfect occupancy/pairing claim',K.BIOLOGICAL,primary=True)
    lysis=g.input('lysis_buffer','LiCl/Tris +1% lithium dodecyl sulfate, EDTA and DTT; dialysis-covered wells',K.CONDITION)
    lysate=g.product('lysate','Lysed cells within sealed wells',K.MOLECULAR)
    g.react('lysis','Lyse cells under dialysis membrane',S.PERMEABILIZATION,
        ((cells,R.TEMPLATE),(lysis,R.CONDITION),(g.input('blocked_membrane','Actual prepared BSA-blocked dialysis membrane; cellulose carrier and protein blocker, not a DNA barcode',K.MOLECULAR),R.STRUCTURAL)),(lysate,))
    beads=g.input('poly_dt','Poly(dT) magnetic beads; capture in wells, recover and wash before RT-PCR',K.SELECTION)
    captured=g.product('captured','Bead-bound heavy/light mRNAs',K.MOLECULAR)
    unbound=g.product('unbound','Removed dialysis membrane and unretained lysate/wash material; RNA recovery is not perfect',K.WASTE)
    g.react('capture','Capture polyadenylated mRNA and magnetically recover/ wash beads',S.AFFINITY_CAPTURE,
        ((lysate,R.TEMPLATE),(beads,R.SELECTION)),(captured,unbound))
    g.discard(unbound,'Removed during bead recovery and washes before emulsification')
    mix=g.input('one_step_mix','Actual 2.85 ml prepared Quanta OneStep Fast RT-PCR reaction mixture with 0.05 wt% BSA and Table 5 primers; commercial-tube packaging unassigned',K.MOLECULAR)
    rt_enzyme=g.product('rt_enzyme','RT enzyme class from actual prepared reaction stock; identity/dose unassigned',K.CATALYST)
    dna_pol=g.product('linkage_polymerase','DNA polymerase class from the same actual prepared reaction stock',K.CATALYST)
    nts=g.product('linkage_nt','Nucleotides copresent in the actual prepared RT-PCR reaction stock',K.MOLECULAR)
    buffer=g.product('linkage_buffer','Actual copresent reaction-buffer class; full formulation unassembled',K.CONDITION)
    bsa=g.product('linkage_bsa','Actual 0.05 wt% BSA molecular carrier in the prepared reaction mixture',K.MOLECULAR)
    subset=Scene()
    for name in ('hVH4-DP63-fwd-OE','hVK5-fwd-OE'):
        subset.strand(name,[Segment(name,PRIMER_EXAMPLES[name])])
    primers=g.product('linkage_primers','Remaining actual Table 5 pool beyond CHrev/CLrev components; two exact representatives shown, other member classes unassembled',K.OLIGO,state=subset)
    ch=g.product('ch_amplifier','Same-pool CHrev-AHX89 component, not an extra fresh primer addition',K.OLIGO,state=Construct([Segment('CHrev-AHX89',PRIMER_EXAMPLES['CHrev-AHX89'])]))
    cl=g.product('cl_amplifier','Same-pool CLrev-BRH06 component, not an extra fresh primer addition',K.OLIGO,state=Construct([Segment('CLrev-BRH06',PRIMER_EXAMPLES['CLrev-BRH06'])]))
    g.react('one_step_inventory','Copresent components of the actual prepared RT-PCR mixture; no physical separation or duplicate primer doses',S.GENERIC,((mix,R.TEMPLATE),),(rt_enzyme,dna_pol,nts,buffer,bsa,primers,ch,cl),changes_molecule=False)
    first=g.product('first_strands','Heavy/light first-strand cDNA in bead-containing emulsion compartments',K.MOLECULAR)
    g.react('rt','Emulsify source bead-reaction stock; RT at55°C for30min then94°C for2min, actual well-array multiplicity unassigned',S.REVERSE_TRANSCRIPTION,
        ((captured,R.TEMPLATE),(primers,R.STRUCTURAL),(ch,R.STRUCTURAL),(cl,R.STRUCTURAL),(rt_enzyme,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION),(bsa,R.STRUCTURAL),(g.input('emulsion_oil','Actual 9 ml mineral-oil phase with4.5% Span80/0.4% Tween80/0.05% TritonX100',K.MOLECULAR),R.STRUCTURAL),(g.input('emulsion_array','Actual 5 min low-speed agitation;100 µl emulsion per PCR well, number of wells unassigned',K.CONDITION),R.CONDITION)),(first,))
    linked=g.product('linked','INFERRED — desired linked VH:VL context within actual emulsion products (~850bp class); not a pure preparation',K.MOLECULAR,state=LINK.product(),chemical_state=LINK)
    byproduct=g.product('byproduct','Complementary-strand by-product drawn in published Supplementary Fig.1c',K.MOLECULAR)
    g.react('linkage_pcr','INFERRED occupied chain contexts — source tail incorporation/annealing/extension in emulsion PCR',S.PCR,
        ((first,R.TEMPLATE),(primers,R.STRUCTURAL),(ch,R.STRUCTURAL),(cl,R.STRUCTURAL),(dna_pol,R.CATALYST),(nts,R.STRUCTURAL),(buffer,R.CONDITION)),(linked,byproduct),
        geometry=PCRGeometry.unavailable(LINK.annealed(),reason='Primer tails are explicit, but individual V(D)J templates and full pool landing coordinates are unassigned. The checked tail complement alone does not validate whole amplicons.'))
    cleanup=g.input('cleanup','Ether emulsion extraction and silica spin-column recovery',K.SELECTION)
    recovered=g.product('recovered','Recovered cDNA mixture including linked and unresolved coproduct classes; not selectively purified by identity',K.MOLECULAR,chemical_state={'copresent_classes':'Desired linked material and source Fig1c coproduct; exact identity/amount/fates unassembled'})
    waste=g.product('emulsion_waste','Discarded oil/ether phases and unretained column material',K.WASTE)
    g.react('recover','Break emulsion, remove oil and recover cDNA',S.CLEANUP,
        ((linked,R.TEMPLATE),(byproduct,R.TEMPLATE),(cleanup,R.SELECTION)),(recovered,waste))
    g.discard(waste,'Oil/ether upper phases and column wash material discarded during recovery')
    heavy_nested=g.input('heavy_nested','Table 6 nested heavy constant-region primer pool',K.OLIGO)
    light_nested=g.input('light_nested','Table 6 nested light constant-region primer pool',K.OLIGO)
    chosen=g.product('nested_template','Actual 4 µl selected from 40 µl recovered cDNA for nested PCR',K.MOLECULAR,chemical_state=recovered.chemical_state)
    rest=g.product('recovered_stock','Nominal 36 µl remaining recovered cDNA stock',K.MOLECULAR,chemical_state=recovered.chemical_state)
    g.react('nested_input','Actual 4-of-40 µl nested-PCR input allocation',S.SPLIT,((recovered,R.TEMPLATE),),(chosen,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferred recovered source cDNA stock'))
    supply=g.input('nested_supplies','Actual ThermoPol/Taq nested-PCR aggregate supplies; stock count/packaging unassigned',K.MOLECULAR)
    taq=g.product('taq','Actual Taq component from nested-PCR supplies',K.CATALYST)
    nested_nt=g.product('nested_nt','Actual nucleotide component from nested-PCR supplies; no invented dose',K.MOLECULAR)
    nested_buffer=g.product('nested_buffer','Actual ThermoPol buffer class and source 200 µl/400 nM reaction context',K.CONDITION)
    g.react('nested_inventory','Actual nested-PCR supply classes; not a master-mix packaging claim',S.GENERIC,((supply,R.TEMPLATE),),(taq,nested_nt,nested_buffer),changes_molecule=False)
    amplified=g.product('amplified','Nested amplified product population; coproduct classes/fates remain unassigned',K.MOLECULAR,chemical_state={'source Fig1c_coproduct':'May persist; no sequence/yield or perfect removal assigned'})
    g.react('nested_pcr','Nested PCR (39 cycles in the memory B-cell branch)',S.PCR,
        ((chosen,R.TEMPLATE),(heavy_nested,R.STRUCTURAL),(light_nested,R.STRUCTURAL),(taq,R.CATALYST),(nested_nt,R.STRUCTURAL),(nested_buffer,R.CONDITION)),(amplified,),
        geometry=PCRGeometry.unavailable(LINK.product(),reason='Nested primer sequences are published, but exact chain templates and primer-to-template coordinates are not assigned in this representation.'))
    gel=g.input('gel','Agarose gel electrophoresis and band extraction',K.SELECTION)
    selected=g.product('selected','Gel-selected approximately 850 bp product population; exact coproduct partition unassigned',K.MOLECULAR,chemical_state=amplified.chemical_state)
    rejected=g.product('gel_rejected','Other DNA bands and unrecovered gel material',K.WASTE)
    g.react('gel_select','Gel-purify linked-size product',S.CLEANUP,
        ((amplified,R.TEMPLATE),(gel,R.SELECTION)),(selected,rejected))
    g.discard(rejected,'Excluded by band selection; no assertion that the Fig.1c class is quantitatively removed')
    lib=g.product('library','Actual linked-repertoire preparation submitted for MiSeq; whole finishing and residual coproduct classes unassembled',K.LIBRARY,chemical_state=amplified.chemical_state)
    g.react('sequencing_handoff','Source-limited library preparation handoff',S.MAGIC,
        ((selected,R.TEMPLATE),),(lib,),reason='The inspected workflow ends at gel extraction and submission for MiSeq. Exact sequencing-adapter construction/primer geometry remains unreviewed; no ligation or kit recipe is invented.')
    g.sequence(lib,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':250,'Read 2':250},
        limitation='2×250 paired-end MiSeq interrogates CDR-H3/CDR-L3 at opposite outer portions. Sequencing adapters/priming sites and exact cycle coordinates remain unassigned in this pass.'))
    g.declare_reviewed_scope(('cells','lysis_buffer','lysate','lysis','poly_dt','captured','unbound','capture',
        'ch_amplifier','cl_amplifier','linkage_primers','one_step_mix','first_strands','rt','linked','byproduct','linkage_pcr',
        'cleanup','recovered','emulsion_waste','recover','heavy_nested','light_nested','taq','amplified',
        'nested_pcr','gel','selected','gel_rejected','gel_select','library','sequencing_handoff'),SOURCE)
    g.declare_reviewed_scope(('emulsion_oil','emulsion_array','blocked_membrane', 'one_step_inventory', 'nested_input', 'nested_supplies', 'nested_inventory'),'Same source Methods103–107 and publishedFig1/Table 5–6: prepared blocked membrane/2.85ml RT-PCR reaction mixture withBSA andcopresent primer components,4-of40ul cDNA nested template, actual source supply classes. No extra primer doses or pure-product separation implied.')
    return g.require_publishable()

def render_page():
    g=protocol_graph()
    return assemble_page(g,(head('Paired VH:VL — physical chain linkage'),'<main class="wrap" data-protocol-complete="true"><h1>Paired VH:VL</h1><p>mRNA captured from microwell-isolated cells is converted to physically linked heavy/light-chain amplicons in an emulsion. Native pairing is reported by chain linkage, not by an invented cell barcode.</p><p><a href="https://doi.org/10.1038/nbt.2492">DeKosky et al. (2013)</a> · <a href="01_paired_vh_vl.html">Research note and source scope</a></p>',
        caveat('<strong>INFERRED — occupied chain contexts.</strong> Source primer tails are known, but receptor sequences, gene lengths, occupied primer sites and complete library finishing remain unassembled. The underlying article is a final-version author manuscript, compared with the published supplement; full final/AM equivalence is not asserted. Complementary tails do not prove amplification specificity or correct single-cell pairing.'),
        '<h2>Preparation</h2><ol><li>Capture mRNA on poly(dT) beads in 125-pL microwells sealed by a BSA-blocked dialysis membrane. Recover and wash beads before emulsification.</li><li>Use the actual prepared Quanta RT-PCR reaction mixture with BSA and Table 5 primers. RT is followed by tail incorporation, annealing and overlap-extension PCR.</li><li>Break the emulsion and recover cDNA, carrying source coproduct classes with unresolved identity/fates.</li><li>Select 4 µl of 40 µl recovered cDNA for nested PCR; preserve the nominal 36 µl remaining stock. Gel-select the approximately 850 bp product population.</li><li>Retain a source-limited handoff to the actual MiSeq 2×250 endpoint. No unreviewed adapter kit or whole read-primer geometry is supplied.</li></ol>',
        '<h2>Known primer relation</h2>',panel(LINK.annealed().rows(),caption='INFERRED — occupied heavy/light contexts with the source-known complementary 17-base tails; symbolic chain bodies have no assigned gene length'),panel(primer_pool().rows(),caption='Four exact source Table 5 representatives; remaining degenerate pool members and transcript coverage are not reconstructed'),
        '<h2>Read setup and source limits</h2><p>The source approximately 850 bp product places CDR-H3 and CDR-L3 near opposite outer portions for paired-end 250-base reads. Actual per-cycle coordinates, sequencing adapters/indexing and complete receptor molecules remain unassigned. Supplementary Fig1c’s coproduct remains copresent through recovery; the gel is not claimed to remove it quantitatively. Optional separate-chain amplification and antibody-expression validation are different branches and are not substituted here.</p>',
        '<details><summary>Full molecular workflow, physical fractions and unresolved chemistry</summary>',WorkflowSection('INFERRED occupied chain views — source-qualified linked-repertoire preparation'),'</details></main>'))

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('Source overlap tails are complementary',revcomp(H_TAIL),L_TAIL)
    c('Linked reference retains one overlap',LINK.product().top().count(L_TAIL),1)
    c('Source coproduct remains copresent into recovery',any(i.material=='byproduct'for i in g.steps['recover'].inputs),True)
    c('No imagined coproduct separation', 'byproduct'in g.terminals,False)
    c('No arbitrary library ligation inferred',g.steps['sequencing_handoff'].kind,S.MAGIC)
    from storage import sources_root
    source=sources_root()/Path(__file__).resolve().parents[1].name/'nbt.2492_41587_2013_BFnbt2492_MOESM35_ESM.pdf.txt'
    if source.exists():
        text=''.join(source.read_text().split())
        for name,seq in {**PRIMER_EXAMPLES,**NESTED}.items():c('Published primer transcription '+name,seq in text)
    else:c.skip('Primer transcription','External published supplement unavailable')
    c('CHrev is from actual prepared reaction inventory',g.producer_of('ch_amplifier').id,'one_step_inventory')
    c('CLrev is not an extra external dose',g.producer_of('cl_amplifier').id,'one_step_inventory')
    c('RT uses same-stock nucleotide class',any(i.material=='linkage_nt'for i in g.steps['rt'].inputs),True)
    c('BSA is a molecular reaction component',g.materials['linkage_bsa'].kind,K.MOLECULAR)
    c('Nested PCR uses actual4ul selected stock',g.steps['nested_pcr'].inputs[0].material,'nested_template')
    c('Remaining recovered stock not sequenced',g.terminals['recovered_stock'].disposition,Disposition.UNRESOLVED)
    c('Occupied receptor context defaults inferred',all(seg.inferred for seg in LINK.product()),True)
    c('Normal publication gate',g.publication_problems(),[])
    c.report()
