"""Bounded, source-local molecular workflow; not a complete protocol review."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature, panel, workflow_from_protocol_graph, workflow_panel
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SequencingPlan, Terminal, Disposition, SourceItem
from pcr import PCRGeometry
from protocol_model import role
from page import head, info
PAGE_FILE='chem.html'

def unknown_template(label):
    return Construct([role(label+'; occupied length/sequence unassembled')])

def clean(g,key,material,label,selection,kind=K.MOLECULAR):
    sel=g.input(key+'_selection',selection,K.SELECTION)
    kept=g.product(key+'_kept',label,kind,chemical_state=material.chemical_state)
    waste=g.product(key+'_waste','Separated material from '+key,K.WASTE)
    g.react(key,label,S.CLEANUP,((material,R.TEMPLATE),(sel,R.SELECTION)),(kept,waste))
    g.discard(waste,'Removed only at this source-described physical separation; recovery is not predicted')
    return kept

def pcr(g,key,template,forward,reverse,label,conditions):
    f=g.input(key+'_forward',forward,K.OLIGO,state=primer_state(forward),chemical_state={'source_terminal_bonds':'C*C*G phosphorothioates'} if key=='adapters' else None)
    r=g.input(key+'_reverse',reverse,K.OLIGO,state=primer_state(reverse))
    stock=g.input(key+'_mix','PublishedSI actual50µlKAPAHiFi reagent; AM namesHotstartReadyMix; individual formulation unassembled',K.MOLECULAR)
    cat=g.product(key+'_enzyme','Copresent polymerizing catalyst class from actual source reagent',K.CATALYST)
    nt=g.product(key+'_nucleotides','Copresent nucleotide supply class; source formulation unassembled',K.MOLECULAR)
    buf=g.product(key+'_buffer','Source reagent buffer component',K.CONDITION)
    g.react(key+'_mix_inventory','Inventory components of one actual source50µl reagent addition',S.MAGIC,((stock,R.TEMPLATE),),(cat,nt,buf),changes_molecule=False,reason='Known source reagent addition decomposed into functional classes, not new free doses; full kit formulation not reconstructed.')
    cond=g.input(key+'_conditions',conditions,K.CONDITION)
    ports=[(template,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(cat,R.CATALYST),(nt,R.STRUCTURAL),(buf,R.CONDITION),(cond,R.CONDITION)]
    if key in ('outer','inner'):
        body=CROP_OUTER if key=='outer' else CROP_INNER
        crop=g.input(key+'_crop','SourceSI10µl12µM '+('CROPouter' if key=='outer' else 'CROPinner')+'; AM2µl differs; gene-relative excess not derived from unspecified panel count',K.OLIGO,state=Construct([Segment('CROP primer',body)]))
        ports.append((crop,R.STRUCTURAL))
    out=g.product(key+'_product',label,K.MOLECULAR,chemical_state={'inherited_cell_and_molecule_identifiers':template.chemical_state,'source_sample_index':r.state.get('sample index').feature if key=='adapters' else None,'occupied_product':'unassembled; no inferred whole-strand order'})
    g.react(key,label,S.PCR,ports,(out,),geometry=PCRGeometry.unavailable(unknown_template('unassigned target sequence'),reason='Study-specific target-pool substrates/full primer landing and final-strand assembly remain unimplemented; known universal/CROP polymers remain separate real inputs.'))
    return out

def aliquot(g,key,material,label):
    selected=g.product(key+'_selected',label,material.kind,chemical_state=material.chemical_state)
    remainder=g.product(key+'_remainder','Untransferred '+key+' material',material.kind,chemical_state=material.chemical_state)
    g.react(key,label,S.SPLIT,((material,R.TEMPLATE),),(selected,remainder))
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,'The source specifies an aliquot; remaining reaction material is not claimed discarded.'))
    return selected
SOURCE='https://doi.org/10.1038/s41592-020-0837-5; publishedSI MOESM1 steps16–41 andTableS6, associatedAMMethods189–207 comparison; finalVORprocedure equivalence unverified'
TITLE='TAP-seq — 10x3primev2/CROP-seq branch'
LIMIT='Model begins at purified GEM-RT cDNA, preserving source cell identity and UMI. Historical bead sequence, target-panel landing geometry and exact final-strand bases remain unassigned. CROP-seq and gene targets share the targeted library.'

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    cdna=g.input('cdna','35µl source10x3primev2 GEM-RT cDNA; RT template-switch oligo omitted upstream',K.MOLECULAR,primary=True,chemical_state={'cell_barcode':feature('tap.cell','cell_barcode','unknown',note='Source declares existing bead-derived cell identity; width/coordinate unassigned here'),'umi':feature('tap.umi','umi','unknown',note='Source declares existing UMI; no fresh UMI during targetedPCR, no16+10 assumption')})
    first=pcr(g,'outer',cdna,'Actual source gene-specific outer pool; whole members unassembled, CROPouter separate','PartialRead1 universal bead-handle primer','Outer targetedPCR products','PublishedSI representativePCR1x12 (recommended11–14); AMstandard11 differs. PrintedPCR1components101.5µl beforewaterad100µl unresolved; no negativewater or silent correction')
    first=clean(g,'outer_cleanup',first,'Purified outerPCR products','1.5×AMPureXP; publishedSIonce300µl+twice200µl80percentEtOH;30µl elution')
    first=qc(g,'outer_qubit',first,'PublishedSI1µl eluate forQubitDNAHS')
    second=pcr(g,'inner',aliquot(g,'outer_transfer',first,'10ng outerPCR product'),
        'Actual source gene-specific inner pool; whole members unassembled, CROPinner separate','PartialRead1 universal bead-handle primer','Semi-nested targetedPCR products','KAPA HiFi;8cycles')
    second=clean(g,'inner_cleanup',second,'Purified nestedPCR products','1.5×AMPureXP; publishedSIonce300µl+twice200µl80percentEtOH;30µl elution')
    second=qc(g,'inner_qubit',second,'PublishedSI1µl eluate forQubitDNAHS')
    third=pcr(g,'adapters',aliquot(g,'inner_transfer',second,'10ng nestedPCR product'),
        'Targeted10X primer','Selected Illumina N7XX reverse index primer','Illumina-adapted targetedPCR products','KAPA HiFi;8cycles')
    lib=clean(g,'adapter_cleanup',third,'TAP-seq gene/CROP-seq pooled targeted library','1.5×AMPureXP; publishedSIonce300µl+twice200µl80percentEtOH;30µl elution',K.LIBRARY)
    lib=qc(g,'final_qubit',lib,'PublishedSI1µl final eluate forQubitDNAHS')
    lib=qc(g,'final_bioanalyzer',lib,'PublishedSI1ng final eluate forBioAnalyzer; transfer volume unassigned')
    lib=aliquot(g,'loading',lib,'Selected source sequencing contribution; loading amounts/untransferred stock unassigned')
    g.sequence(lib,SequencingPlan(acquired_roles=('Read 1','Read 2','Index 1'),read_lengths={'Read 1':26,'Read 2':58,'Index 1':8},
        limitation='Main Methods: Read1 cell barcode+UMI26cycles; Read2 cDNA approximately58cycles; i7 sample index8cycles. Read2=58 is the representative approximate source schedule, not a fixed mandatory run length. Exact final strand geometry deferred.'))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    g.declare_reviewed_scope(ADDITIONAL_REVIEWED,SOURCE)
    g.add_source_item(SourceItem('unassembled_source','Full source primer-pool, intermediate and final-strand assembly',LIMIT+' Source material is available for follow-up; not an absence finding.'))
    return g

def qc(g,key,stock,label):
    measured=g.product(key+'_qc',label,K.QC)
    remaining=g.product(key+'_remaining','Remaining source stock after actualQC aliquot',stock.kind,chemical_state=stock.chemical_state)
    g.react(key,label,S.SPLIT,((stock,R.TEMPLATE),),(measured,remaining))
    g.mark_terminal(Terminal(measured.id,Disposition.QC,'Actual source QC aliquot; not sequenced or molecularly copied.'))
    return remaining

def detail():return ''

REVIEWED=('cdna', 'outer_forward', 'outer_reverse', 'outer_enzyme', 'outer_product', 'outer_cleanup_selection', 'outer_cleanup_kept', 'outer_cleanup_waste', 'outer_transfer_selected', 'outer_transfer_remainder', 'inner_forward', 'inner_reverse', 'inner_enzyme', 'inner_product', 'inner_cleanup_selection', 'inner_cleanup_kept', 'inner_cleanup_waste', 'inner_transfer_selected', 'inner_transfer_remainder', 'adapters_forward', 'adapters_reverse', 'adapters_enzyme', 'adapters_product', 'adapter_cleanup_selection', 'adapter_cleanup_kept', 'adapter_cleanup_waste', 'outer', 'outer_cleanup', 'outer_transfer', 'inner', 'inner_cleanup', 'inner_transfer', 'adapters', 'adapter_cleanup')
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


import illumina as il
import nextera as nx
PARTIAL_READ1='CTACACGACGCTCTTCCGATCT'
CROP_OUTER='AATGGACTATCATATGCTTACCGT'
CROP_INNER=nx.ADAPTOR_S7+'TGTGGAAAGGACGAAACACC'
TARGETED10X=il.P5+'TCTTTCCCTACACGACGCTCTTCCG'
SOURCE_FILE='s41592-020-0837-5_supp_EMS118361-supplement-Table_S6.xlsx.txt'
TRANSCRIBED=(PARTIAL_READ1,CROP_OUTER,CROP_INNER,TARGETED10X[:-3]+'C*C*G')
def primer_state(label):
    if label=='PartialRead1 universal bead-handle primer':return Construct([Segment('PartialRead1',PARTIAL_READ1)])
    if label=='Targeted10X primer':
        sc=Scene();sc.strand('primer',[Segment('Targeted10X',TARGETED10X)],label='Targeted10X; terminal C*C*G phosphorothioate linkages')
        return sc
    if label=='Selected Illumina N7XX reverse index primer':return Construct([
        Segment('P7',il.P7),Segment('sample index','N'*8,'idx',placeholder=True,feature=feature('tap.sample','sample_index','unknown')),
        Segment('N7 primer tail','GTCTCGTGGGCTCGG')])
    return None

def source_checks(c,g):
    c('Actual CROPouter stays in same targetedPCR reaction','outer_crop' in tuple(p.material for p in g.steps['outer'].inputs),True)
    c('Sample index role is structured',primer_state('Selected Illumina N7XX reverse index primer').get('sample index').feature.role,'sample_index')

def render_page():
    return '\n'.join([head(TITLE),'<div class="wrap"><h1>'+TITLE+'</h1>',
        '<p class="research-notes"><a href="01_tap_seq.html">Research notes</a></p>',
        info(LIMIT),info(SOURCE_NOTICES),detail(),''.join(workflow_panel(workflow_from_protocol_graph(protocol_graph(),endpoint)) for endpoint in ['loading_selected']),'</div>'])

def selftest():
    from checks import Check
    c=Check(); g=protocol_graph()
    issues=g.publication_problems()
    c('Only explicit unassembled-source coverage blocks publication',issues,["authoritative source item 'unassembled_source' is not accounted for"])
    c('Private preview does not claim protocol completion','data-protocol-complete="true"' in render_page(),False)
    extra_checks(c,g)
    ports=lambda key:tuple(p.material for p in g.steps[key].inputs)
    c('separate actual CROPinner component present','inner_crop' in ports('inner'),True)
    c('CROPouter exact polymer retained',g.materials['outer_crop'].state.top(),CROP_OUTER)
    c('CROPinner exact polymer retained',g.materials['inner_crop'].state.top(),CROP_INNER)
    c('one actual outer ReadyMix stock',ports('outer_mix_inventory'),('outer_mix',))
    c('actual nucleotide component consumed','outer_nucleotides' in ports('outer'),True)
    c('actual inner ReadyMix nucleotide consumed','inner_nucleotides' in ports('inner'),True)
    c('source finalQubit typedQC',g.materials['final_qubit_qc'].kind,K.QC)
    c('source finalBioAnalyzer typedQC',g.materials['final_bioanalyzer_qc'].kind,K.QC)
    c('source loading remainder retained',g.terminals['loading_remainder'].disposition,Disposition.UNRESOLVED)
    c('no invented30base target placeholder',list(unknown_template('test'))[0].is_role_token(),True)
    c('source existing cell identity structured',g.materials['cdna'].chemical_state['cell_barcode'].role,'cell_barcode')
    c('source UMI not redefined from26cycles',g.materials['cdna'].chemical_state['umi'].role,'umi')
    c('required fullassembly not excluded','unassembled_source' in g.source_coverage,False)
    c.report(exit_on_fail=True)

ADDITIONAL_REVIEWED=('outer_mix','outer_mix_inventory','outer_conditions','outer_crop','inner_mix','inner_mix_inventory','inner_conditions','inner_crop','adapters_mix','adapters_mix_inventory','adapters_conditions','outer_qubit','inner_qubit','final_qubit','final_bioanalyzer','loading','loading_selected')

SOURCE_NOTICES='PRIVATE MODEL-SCOPE HOLD: full target-pool/intermediate/final-strand assembly remains required and unimplemented. This model now follows publishedSI PCR settings: representative12-cyclePCR1 (11–14 advised),10µl12µMCROP primers and three80percentEtOHwashes. AssociatedAMMethods instead gives11cycles/2µlCROP/twowashes; these source-version/protocol differences are not automatically errors. INTRINSIC publishedSI PDFp17step16table prints35+4+2.5+10+50=101.5µl beforeH2O ad100µl. Values are retained; no negativewater or silent repair is inferred. Optional5µlGEMRT set-aside elsewhere does not silently change printed35µl tableinput. No full finalVOR/AMequivalence or newer10xv3/NextGEM recipe is assumed.'
