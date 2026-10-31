"""Source-scoped original Shahi Abseq, not a BD kit or CITE-style RT recipe."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, panel, protocol_graph_panels
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S,
    SequencingPlan, Terminal, Disposition, SourceDiscrepancy, SourceClaim, AnnotationSlot)
from protocol_views import render_annotations, assemble_page, WorkflowSection
from page import head, info, caveat, table

PAGE_FILE='abseq.html'
SOURCE='https://doi.org/10.1038/srep44447; Results Microfluidic workflow/Bioinformatic workflow; Methods'
P5_PRIMER=il.P5[:20]
P7_PRIMER=il.P7
PROBE_3TAQ1=il.INDEX2_PRIMER_RC[:23]
FILTER_HOMOLOGY='GCAGTGGTATCAACGCAGAG'
ANTIBODY_IDENTITIES={'CD3':'TTATAAC','CD19':'TTAATTG'}
UMI_NT=8

# Only the source-declared identity-bearing components are represented here.
# Their polymer order, flanks and full strand boundaries are NOT reconstructed.
def identity_components():
    sc=Scene()
    sc.strand('antibody identity',[Segment('CD3 identity',ANTIBODY_IDENTITIES['CD3'],
        feature=feature('abseq.antibody','feature_barcode','fixed'))],label='CD3 tag component')
    sc.strand('CD19 identity',[Segment('CD19 identity',ANTIBODY_IDENTITIES['CD19'],feature=feature('abseq.antibody','feature_barcode','fixed'))],label='Alternative CD19 tag component; not joined to CD3')
    sc.strand('umi',[Segment('UMI','N'*UMI_NT,'umi',placeholder=True,
        feature=feature('abseq.umi','umi','random'))],label='UMI component')
    sc.note('umi','Separate components: positions and complete tag strand are unresolved in this pass')
    return sc

def primer_scene():
    sc=Scene()
    sc.strand('p5',[Segment('source P5 primer',P5_PRIMER)],label='P5 primer')
    sc.strand('p7',[Segment('source P7 primer',P7_PRIMER)],label='P7 primer')
    return sc


def supplies(g,key,label):
    stock=g.input(key,label,K.MOLECULAR)
    enzyme=g.product(key+'_enzyme','Actual polymerase class from '+label,K.CATALYST)
    nts=g.product(key+'_nts','Actual nucleotide class from the same source supplies; no fresh extra dose',K.MOLECULAR)
    buffer=g.product(key+'_buffer','Actual copresent reaction-buffer class; formulation unassembled',K.CONDITION)
    g.react(key+'_inventory','Copresent source supply classes; not physical separation or independently added doses',S.GENERIC,((stock,R.TEMPLATE),),(enzyme,nts,buffer),changes_molecule=False)
    return enzyme,nts,buffer

def unknown_product():
    sc=Scene();sc.strand('unassembled',[Segment('Actual linked product context unassembled','[unassembled product]',placeholder=True,inferred=True)]);return sc

def protocol_graph():
    g=ProtocolGraph('Abseq — original droplet SOE-PCR branch',authority=SOURCE,representation='published_instance')
    cells=g.input('cells','Actual Jurkat/Raji cell population for surface-protein measurement; mixing proportions unassigned',K.BIOLOGICAL,primary=True)
    conjugates=g.input('conjugates','Purified antibody–DNA conjugates: 5′ thiol oligos linked with SM(PEG)6; protein/DNA complex, not pure DNA',K.MOLECULAR,chemical_state={'source_components':'CD3/CD19 antibody proteins and thiol-modified DNA reporters with SM(PEG)6 linkage; whole tag bodies/stoichiometry unassembled'})
    stained=g.product('stained','Antibody-DNA bound cells plus excess conjugates',K.MOLECULAR)
    g.react('stain','Block/Fc-block and stain cells with conjugated antibodies on ice',S.BINDING,
        ((cells,R.TEMPLATE),(conjugates,R.STRUCTURAL),(g.input('blocking_mix','Actual prepared blocking buffer contributions: FBS/goat/mouse serum proteins and salmon-sperm dsDNA; stock allocations unassigned',K.MOLECULAR),R.STRUCTURAL),(g.input('fc_block','Actual 2.5 µg Fc Block contribution; source prepared protein reagent',K.MOLECULAR),R.STRUCTURAL)),(stained,))
    wash=g.input('facs_buffer','Actual FACS buffer PBS/2% FBS protein-carrier contributions; five 1 ml washes',K.MOLECULAR)
    washed=g.product('washed','Washed DNA-labeled cells',K.MOLECULAR)
    unbound=g.product('unbound','Unbound conjugates in washes',K.WASTE)
    g.react('wash','Remove unbound conjugates by cell washing',S.CLEANUP,
        ((stained,R.TEMPLATE),(wash,R.STRUCTURAL)),(washed,unbound))
    g.discard(unbound,'Removed with cell-staining washes')
    pk=g.input('proteinase','Actual proteinase K contribution at 1.5 µg/µl in source lysis stream',K.CATALYST)
    released=g.product('released','Antibody DNA tags released in single-cell lysate droplets',K.MOLECULAR)
    g.react('proteolyse','Protease lysis/tag release at 55°C for 30 min',S.DEPROTEINIZATION,
        ((washed,R.TEMPLATE),(pk,R.CATALYST),(g.input('lysis_environment','Source 100 mM Tris/2% Tween20; cell stream resuspended in 100 µl FACS buffer with 17% Optiprep; coflow droplet interface unassembled',K.CONDITION),R.CONDITION),(g.input('droplet_oil','Actual HFE7500/PEG-PFPE oil/surfactant contributions for cell droplets; no genome-derived material implied',K.MOLECULAR),R.STRUCTURAL)),(released,))
    heat=g.input('heat','95°C for 15 min',K.CONDITION)
    inactive=g.product('inactive','Cell-lysate droplets with protease heat-inactivated',K.MOLECULAR)
    g.react('heat_inactivate','Inactivate proteinase K before PCR',S.THERMAL_TREATMENT,
        ((released,R.TEMPLATE),(heat,R.CONDITION)),(inactive,))
    randomer=g.input('barcode_template','Actual dilute random cell-barcode templates with flanks; ~1 template per10 droplets, blank/multiple-template classes not quantified here',K.OLIGO,chemical_state={'unpositioned_identifier':feature('abseq.cell','cell_barcode','random',note='Source10mer versus8N notation unresolved; no exact width or full strand assigned')})
    barcode_forward=g.input('barcode_forward','Forward primer for clonal barcode amplification; exact sequence unassigned',K.OLIGO)
    barcode_reverse=g.input('barcode_reverse','Reverse primer for clonal barcode amplification; exact sequence unassigned',K.OLIGO)
    barcode_polymerase,barcode_nt,barcode_buffer=supplies(g,'barcode_polymerase','Actual barcode-PCR aggregate supplies; stock packaging/formulation unassigned')
    barcode=g.product('clonal_barcode','Clonal barcode-containing droplet class within source emulsion; blank/other compartments remain unassembled',K.MOLECULAR)
    g.react('barcode_pcr','Amplify dilute single barcode molecules in droplets',S.PCR,
        ((randomer,R.TEMPLATE),(barcode_forward,R.STRUCTURAL),(barcode_reverse,R.STRUCTURAL),(barcode_polymerase,R.CATALYST),(barcode_nt,R.STRUCTURAL),(barcode_buffer,R.CONDITION)),(barcode,),
        geometry=PCRGeometry.unavailable(identity_components(),reason='The barcode template/flanking-primer full strands are not established by this capped source pass.'))
    pool=g.product('merged','Mother droplets after source electrocoalescence and mixing; productive lineage shown, occupancy/incorrect mergers unassembled',K.MOLECULAR)
    reagent,soe_nt,soe_buffer=supplies(g,'soe_reagents','Actual Platinum Multiplex PCR Master Mix allocation in source reagent droplets')
    tween=g.input('soe_tween','Actual 2% w/v Tween20 additive, separate from prepared master mix',K.MOLECULAR)
    peg=g.input('soe_peg','Actual 2% w/v PEG6000 additive, separate from prepared master mix',K.MOLECULAR)
    g.react('droplet_merge','Electrocoalesce one lysate and one clonal barcode droplet with PCR reagents',S.POOL,
        ((inactive,R.TEMPLATE),(barcode,R.TEMPLATE),(reagent,R.CATALYST),(soe_nt,R.STRUCTURAL),(soe_buffer,R.CONDITION),(tween,R.STRUCTURAL),(peg,R.STRUCTURAL)),(pool,))
    daughters=tuple(g.product('daughter_'+str(i),'Actual equal daughter portion '+str(i)+' of mixed mother droplets; not a cloned PCR product',K.MOLECULAR)for i in range(1,5))
    g.react('split_four','Source physically splits mixed mother droplets into four equal PCR-stable portions BEFORE cycling',S.SPLIT,((pool,R.TEMPLATE),),daughters)
    p5=g.input('p5','Source-written 20-base P5 primer',K.OLIGO,state=Construct([Segment('P5',P5_PRIMER)]))
    p7=g.input('p7','Source-written P7 primer',K.OLIGO,state=Construct([Segment('P7',P7_PRIMER)]))
    linked=g.product('linked','Actual collection of SOE products from four separately cycled daughter portions; background/coproduct classes unassembled',K.MOLECULAR)
    g.react('soe','Actual four-portion parallel SOE-PCR interface: 11 cycles, no pooling before cycling or fresh reagent-drop addition',S.PCR,
        tuple((x,R.TEMPLATE)for x in daughters)+((p5,R.STRUCTURAL),(p7,R.STRUCTURAL),(reagent,R.CATALYST),(soe_nt,R.STRUCTURAL),(soe_buffer,R.CONDITION)),(linked,),
        geometry=PCRGeometry.unavailable(identity_components(),reason='Source states shared terminal homology and a linked product, but complete input strand order is unresolved here. No exact overlap/product geometry is claimed.'))
    column=g.input('column','Perfluorooctanol emulsion break and Zymo DNA Clean & Concentrator',K.SELECTION)
    recovered=g.product('recovered','Column-recovered DNA',K.MOLECULAR)
    column_waste=g.product('column_waste','Oil and unretained recovery material',K.WASTE)
    g.react('recover','Break emulsion and recover DNA on column',S.CLEANUP,
        ((linked,R.TEMPLATE),(column,R.SELECTION)),(recovered,column_waste))
    g.discard(column_waste,'Discarded during emulsion break/column purification')
    exo=g.input('exo','Exonuclease I',K.CATALYST)
    treated=g.product('treated','Linked DNA after unlinked ssDNA tag digestion',K.MOLECULAR)
    g.react('exo_digest','Digest remaining unlinked ssDNA antibody tags',S.DIGESTION,
        ((recovered,R.TEMPLATE),(exo,R.CATALYST)),(treated,))
    bulk_pol,bulk_nt,bulk_buf=supplies(g,'bulk_supplies','Actual bulk-PCR aggregate supplies; exact formulation/stock counts not assigned from SOE conditions')
    bulk_f=g.input('bulk_forward','Actual bulk forward primer class; whole order unassigned for this stage',K.OLIGO)
    bulk_r=g.input('bulk_reverse','Actual bulk reverse primer class; whole order unassigned for this stage',K.OLIGO)
    library=g.product('library','Bulk-amplified linked Abseq library stock; background/coproduct classes unassembled',K.LIBRARY)
    g.react('bulk_pcr','Second bulk PCR supplies sequencing input',S.PCR,
        ((treated,R.TEMPLATE),(bulk_f,R.STRUCTURAL),(bulk_r,R.STRUCTURAL),(bulk_pol,R.CATALYST),(bulk_nt,R.STRUCTURAL),(bulk_buf,R.CONDITION)),(library,),
        geometry=PCRGeometry.unavailable(unknown_product(),reason='Bulk-PCR strand architecture and conditions remain unassigned in this bounded pass.'))
    qc=g.product('final_qc','Actual Qubit/Bioanalyzer/qPCR analytical allocations; amounts unassigned',K.QC)
    stock=g.product('final_stock','Library stock after actual analytical allocations',K.LIBRARY)
    g.react('final_qc_allocate','Actual source library QC allocations, not production DNA',S.SPLIT,((library,R.TEMPLATE),),(qc,stock));g.mark_terminal(Terminal(qc.id,Disposition.QC,'Source analytical allocations'))
    load=g.product('sequenced','Actual selected 15 pM MiSeq loading allocation',K.LIBRARY);rest=g.product('library_stock','Untransferred prepared library stock; amount unassigned',K.LIBRARY)
    g.react('load_allocate','Actual sequencing loading contribution',S.SPLIT,((stock,R.TEMPLATE),),(load,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Untransferred source stock'))
    g.sequence(load,SequencingPlan(acquired_roles=('Read 1',),read_lengths={'Read 1':50},
        limitation='Article specifies MiSeq v2 single-read 50-cycle sequencing, but exact primer site and per-cycle antibody/UMI/cell-barcode spans are not reconstructed.'))
    g.declare_reviewed_scope(('cells','conjugates','stain','stained','facs_buffer','washed','unbound','wash',
        'proteinase','released','proteolyse','heat','inactive','heat_inactivate','barcode_template',
        'barcode_forward','barcode_reverse','barcode_polymerase','clonal_barcode','barcode_pcr','merged','soe_reagents',
        'droplet_merge','p5','p7','linked','soe','column','recovered','column_waste','recover','exo',
        'treated','exo_digest','library','bulk_pcr'),SOURCE)
    g.annotate(SourceDiscrepancy(target='barcode_template',title='Cell-barcode width differs within the article',
        claims=(SourceClaim('Published Results, barcode-emulsion generation and bioinformatics','Random 10-mers; cell barcodes described as 10 bp.'),
                SourceClaim('Published Methods, Bioinformatics','Cell Barcodes printed as NNNN NNNN (eight N characters).')),
        model_choice='No exact cell-barcode length or full molecule is assigned in this model.',
        evidence='Both descriptions are in the defining article; supplementary/reference checks remain incomplete.',
        uncertainty='Whether the Methods notation is abbreviated or represents a different barcode construct remains unresolved.',
        slot=AnnotationSlot.BEFORE_WORKFLOW))
    g.declare_reviewed_scope(('blocking_mix', 'fc_block', 'lysis_environment', 'droplet_oil', 'barcode_polymerase_inventory', 'soe_reagents_inventory', 'soe_tween', 'soe_peg', 'split_four', 'bulk_supplies', 'bulk_supplies_inventory', 'bulk_forward', 'bulk_reverse', 'final_qc_allocate', 'load_allocate', 'sequenced'),'Same publishedPMC5349531 Methods132/140/144 and Results68–82: source molecular carriers/PCR supply classes/four actual pre-PCR daughter portions/analytical and15pM loading allocations. Array inventory views are not fresh duplicate doses; bulk formulation/primers unassembled.')
    return g.require_publishable()

def render_page():
    g=protocol_graph()
    return assemble_page(g,(head('Abseq — original droplet protein reporters'),'<main class="wrap" data-protocol-complete="true"><h1>Abseq</h1><p>DNA-tagged antibodies report cell-surface proteins. Droplet SOE-PCR physically links reporter identity and UMI material to a clonal cell barcode.</p><p><a href="https://doi.org/10.1038/srep44447">Shahi et al. (2017)</a> · <a href="01_abseq.html">Research note</a></p>',
        caveat('This source-qualified frontend keeps full conjugate/tag strands, occupied SOE overlap geometry and bulk-PCR/finishing details unassembled. It is the original Shahi method, not a current BD kit or antibody RT recipe. Productive barcode-containing lineages are represented; source dilution (~1template per10drops) and stochastic mergers do not guarantee every compartment is occupied or correctly paired.'),
        '<h2>Preparation</h2><ol><li>Block and stain Jurkat/Raji cells with purified antibody–DNA complexes, then wash away unbound material.</li><li>Lyse in proteinase-K droplets, heat-inactivate the protease and separately prepare a dilute clonal-barcode emulsion.</li><li>Merge lysate, barcode and reagent droplets; mix and split each mother droplet into four actual portions before 11 cycles of SOE-PCR.</li><li>Break the emulsion, recover DNA, digest unlinked ssDNA with ExoI and carry out source-positive bulk PCR with its detailed supplies left unassigned.</li><li>Allocate Qubit/Bioanalyzer/qPCR material and sequence a 15 pM library contribution on MiSeq v2: single-read 50 cycles.</li></ol>',
        '<h2>Known components</h2>',panel(identity_components().rows(),caption='Source identity components displayed separately; CD3/CD19 alternatives and random eight-base UMI do not reconstruct full tag order'),panel(primer_scene().rows(),caption='Source-written20-base P5 and P7 primers for SOE analysis; no full library or bulk-PCR primer identity is inferred'),
        '<h2>Read and source limits</h2><p>The single 50-cycle read carries antibody/cell/UMI information, but exact cycle spans and full primer landing remain unresolved. The FAM 3Taq1 probe and bioinformatic filtering motif are different reported sequences; neither is silently used to reconstruct a shared overlap. PCR/NGS barcode clustering is computational filtering, not a physical molecular separation.</p>',
        AnnotationSlot.BEFORE_WORKFLOW,'<details><summary>Full source workflow, physical fractions and unresolved chemistry</summary>',WorkflowSection('Source-qualified productive-lineage preparation'),'</details></main>'))

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('No reverse transcription substituted for DNA SOE-PCR',any(x.kind is S.REVERSE_TRANSCRIPTION for x in g.steps.values()),False)
    c('Proteinase is inactivated upstream of droplet merger',g.steps['droplet_merge'].inputs[0].material,'inactive')
    c('Source disagreement keeps barcode width unassigned',g.materials['barcode_template'].state,None)
    from storage import sources_root
    source=sources_root()/Path(__file__).resolve().parents[1].name/'srep44447_PMC5349531.xml.txt'
    if source.exists():
        text=''.join(source.read_text().split())
        for seq in (P5_PRIMER,P7_PRIMER,PROBE_3TAQ1,FILTER_HOMOLOGY):c('Source transcription '+seq,seq in text)
    else:c.skip('Source transcription','External article unavailable')
    c('Conjugates are protein/DNA complexes',g.materials['conjugates'].kind,K.MOLECULAR)
    c('Four daughter portions exist before PCR',len(g.steps['split_four'].outputs),4)
    c('SOE consumes actual daughter templates',tuple(i.material for i in g.steps['soe'].inputs if i.role is R.TEMPLATE),tuple('daughter_'+str(i)for i in range(1,5)))
    c('Bulk reagents not inherited from SOE',any(i.material=='bulk_supplies_enzyme'for i in g.steps['bulk_pcr'].inputs),True)
    c('QC is separate from production',g.terminals['final_qc'].disposition,Disposition.QC)
    c('Known phenotype barcode meaning retained',g.materials['barcode_template'].chemical_state['unpositioned_identifier'].role,'cell_barcode')
    c('Normal publication gate',g.publication_problems(),[])
    c.report()
