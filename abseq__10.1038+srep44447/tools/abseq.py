"""Source-scoped original Shahi Abseq, not a BD kit or CITE-style RT recipe."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, panel, protocol_graph_panels
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S,
    SequencingPlan, SourceDiscrepancy, SourceClaim, AnnotationSlot)
from protocol_views import render_annotations
from page import head, info

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
    sc.strand('umi',[Segment('UMI','N'*UMI_NT,'umi',placeholder=True,
        feature=feature('abseq.umi','umi','random'))],label='UMI component')
    sc.note('umi','Separate components: positions and complete tag strand are unresolved in this pass')
    return sc

def primer_scene():
    sc=Scene()
    sc.strand('p5',[Segment('source P5 primer',P5_PRIMER)],label='P5 primer')
    sc.strand('p7',[Segment('source P7 primer',P7_PRIMER)],label='P7 primer')
    return sc

def protocol_graph():
    g=ProtocolGraph('Abseq — original droplet SOE-PCR branch',authority=SOURCE)
    cells=g.input('cells','Cells for surface-protein measurement',K.BIOLOGICAL,primary=True)
    conjugates=g.input('conjugates','Purified antibody–DNA conjugates: 5′ thiol oligos linked with SM(PEG)6',K.OLIGO)
    stained=g.product('stained','Antibody-DNA bound cells plus excess conjugates',K.MOLECULAR)
    g.react('stain','Block/Fc-block and stain cells with conjugated antibodies on ice',S.BINDING,
        ((cells,R.TEMPLATE),(conjugates,R.STRUCTURAL)),(stained,))
    wash=g.input('facs_buffer','FACS buffer: PBS + 2% FBS; five washes',K.SELECTION)
    washed=g.product('washed','Washed DNA-labeled cells',K.MOLECULAR)
    unbound=g.product('unbound','Unbound conjugates in washes',K.WASTE)
    g.react('wash','Remove unbound conjugates by cell washing',S.CLEANUP,
        ((stained,R.TEMPLATE),(wash,R.SELECTION)),(washed,unbound))
    g.discard(unbound,'Removed with cell-staining washes')
    pk=g.input('proteinase','Proteinase K lysis buffer in co-flow droplets',K.CATALYST)
    released=g.product('released','Antibody DNA tags released in single-cell lysate droplets',K.MOLECULAR)
    g.react('proteolyse','Protease lysis/tag release at 55°C for 30 min',S.DEPROTEINIZATION,
        ((washed,R.TEMPLATE),(pk,R.CATALYST)),(released,))
    heat=g.input('heat','95°C for 15 min',K.CONDITION)
    inactive=g.product('inactive','Cell-lysate droplets with protease heat-inactivated',K.MOLECULAR)
    g.react('heat_inactivate','Inactivate proteinase K before PCR',S.THERMAL_TREATMENT,
        ((released,R.TEMPLATE),(heat,R.CONDITION)),(inactive,))
    randomer=g.input('barcode_template','Single random cell-barcode template with flanking priming sequences; width conflicting in source',K.OLIGO)
    barcode_forward=g.input('barcode_forward','Forward primer for clonal barcode amplification; exact sequence unassigned',K.OLIGO)
    barcode_reverse=g.input('barcode_reverse','Reverse primer for clonal barcode amplification; exact sequence unassigned',K.OLIGO)
    barcode_polymerase=g.input('barcode_polymerase','PCR mixture for barcode emulsion',K.CATALYST)
    barcode=g.product('clonal_barcode','Clonal barcode amplicons within individual droplets',K.MOLECULAR)
    g.react('barcode_pcr','Amplify dilute single barcode molecules in droplets',S.PCR,
        ((randomer,R.TEMPLATE),(barcode_forward,R.STRUCTURAL),(barcode_reverse,R.STRUCTURAL),(barcode_polymerase,R.CATALYST)),(barcode,),
        geometry=PCRGeometry.unavailable(identity_components(),reason='The barcode template/flanking-primer full strands are not established by this capped source pass.'))
    pool=g.product('merged','Matched lysate/barcode/PCR droplets; mixed and divided into four PCR-stable portions',K.MOLECULAR)
    reagent=g.input('soe_reagents','Platinum Multiplex PCR mix, Tween20 and PEG6000',K.CATALYST)
    g.react('droplet_merge','Electrocoalesce one lysate and one clonal barcode droplet with PCR reagents',S.POOL,
        ((inactive,R.TEMPLATE),(barcode,R.TEMPLATE),(reagent,R.CATALYST)),(pool,))
    p5=g.input('p5','Source-written 20-base P5 primer',K.OLIGO,state=Construct([Segment('P5',P5_PRIMER)]))
    p7=g.input('p7','Source-written P7 primer',K.OLIGO,state=Construct([Segment('P7',P7_PRIMER)]))
    linked=g.product('linked','Cell barcode linked to antibody-tag molecules by SOE-PCR',K.MOLECULAR)
    g.react('soe','11 cycles of droplet SOE-PCR',S.PCR,
        ((pool,R.TEMPLATE),(p5,R.STRUCTURAL),(p7,R.STRUCTURAL),(reagent,R.CATALYST)),(linked,),
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
    library=g.product('library','Bulk-amplified linked Abseq sequencing library',K.LIBRARY)
    g.react('bulk_pcr','Second bulk PCR supplies sequencing input',S.PCR,
        ((treated,R.TEMPLATE),(p5,R.STRUCTURAL),(p7,R.STRUCTURAL),(reagent,R.CATALYST)),(library,),
        geometry=PCRGeometry.unavailable(identity_components(),reason='Bulk-PCR strand architecture and conditions remain unassigned in this bounded pass.'))
    g.sequence(library,SequencingPlan(acquired_roles=('Read 1',),read_lengths={'Read 1':50},
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
    return g

def render_page():
    g=protocol_graph()
    return '\n'.join([head('Abseq — bounded source model'),'<div class="wrap"><h1>Abseq</h1>',
        '<p class="research-notes"><a href="01_abseq.html">Research notes</a></p>',
        info('Original Shahi et al. Abseq: DNA reporter tags are linked to cell barcodes by droplet SOE-PCR. This capped model has unresolved full tag/primer geometry; the separate components below are not a reconstructed library.'),
        render_annotations(g,AnnotationSlot.BEFORE_WORKFLOW),protocol_graph_panels(g,cls='long'),
        panel(identity_components().rows(),caption='Source-declared identity components, displayed separately; no sequence order or cut geometry implied.'),'</div>'])

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
    c.report()
