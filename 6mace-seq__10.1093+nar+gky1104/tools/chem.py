"""Protocol-owned source-limited model; no finished-publication marker."""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / 'lib')]
from chemdraw import Scene, Segment, Construct, feature, revcomp, oligo, panel, workflow_from_protocol_graph, workflow_panel
from protocol_model import role, molecule
from protocol_graph import ProtocolGraph, InputRole as R, MaterialKind as K, StepKind as S, Terminal, Disposition, SourceItem
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from pcr import PCRGeometry
from page import head, info

def mixture(*species):
    sc = Scene()
    for n, (label, parts) in enumerate(species):
        sc.strand('s' + str(n), list(parts), label=label)
    return sc

def cleanup(g, ident, source, template, state):
    reagent = g.input(ident + '_condition', source, K.CONDITION)
    out = g.product(ident, 'recovered source preparation', K.MOLECULAR, state=state)
    waste = g.product(ident + '_waste', 'physically separated workup fractions; exact molecular partition unquantified', K.WASTE)
    g.react(ident + '_step', source, S.CLEANUP, ((template, R.TEMPLATE), (reagent, R.CONDITION)), (out, waste))
    g.discard(waste, 'Actual source physical separation; no perfect recovery/depletion asserted.')
    return out

def end_boundary(label):
    sc = molecule('unassembled ligation-ready geometry', (role(label),))
    return AdapterLigationGeometry('unavailable', 'unavailable', (sc,), label + '; no source-absence claim.')

def render_page():
    import html
    g = protocol_graph()
    parts = [head(TITLE), '<main class="wrap"><h1>' + html.escape(TITLE) + '</h1>',
             '<p><a href="' + NOTES + '">Research notes</a></p>', info(SOURCE),
             info('Private incomplete model. ' + CAVEAT)]
    if callable(globals().get('oligos')):
        parts.extend(['<h2>Source ordering oligos</h2>', *oligos()])
    for label, terminal in PREVIEWS:
        parts.extend(['<h2>' + html.escape(label) + '</h2>', workflow_panel(workflow_from_protocol_graph(g, terminal))])
    parts.extend([info(READOUT), '</main>'])
    return '\n'.join(parts)

def selftest():
    from checks import Check
    c = Check(); g = protocol_graph()
    c('Pending source/assembly boundary prevents finished publication', bool(g.publication_problems()))
    from source_oligo_check import candidates
    c('Known PG polymers reach the actual source transcription guard',{x.material for x in candidates(g)},{'reverse_adapter','pg162','forward_ligation_adapter'})
    c('Full end and TruSeq geometry obligation stays uncovered','pending_assembly' in g.source_coverage,False)
    c('Private preview has no finished marker', 'data-protocol-complete=' not in render_page())
    c('Every physical DNA output has a declared downstream operation or terminal', not any(m.kind in {K.MOLECULAR,K.LIBRARY} and m.id not in g.terminals and not any(port.material==m.id for step in g.steps.values() for port in step.inputs) for m in g.materials.values()))
    for step_id in ('lambda','recj'):
        c(step_id+' digestion does not physically discard products',all(g.materials[mid].kind is not K.WASTE for mid in g.steps[step_id].outputs))
    from checks import Source,have
    from storage import sources_root
    src=Source(sources_root()/HERE.parent.name/'nar_gky1104_supp_gky1104_supplemental_files_Table_S1.xlsx.txt','Published gky1104 TableS1, standard get_sources bundle')
    if have(c,src,label='published adapter transcription'):
        table=src.path.read_text()
        for name,(seq,mod5,mod3) in SOURCE_OLIGOS.items():
            c(name+' exact published table row including terminal phosphate', any([col.strip() for col in line.split('\t')[:2]]==[name,('/5Phos/' if mod5 else '')+seq] for line in (v.strip() for v in table.splitlines())))
    c.report()

PAGE_FILE='6mace-seq.html'
TITLE='6mACE-seq (DNA)'
NOTES='01_6mace-seq.html'
SOURCE='<a href="https://doi.org/10.1093/nar/gky1104">Defining article</a>; Published main Methods82–84; published TableS1 rows3–7; Supplementary Methods/legends relevant headings and FigS1–S2'
CAVEAT='DNA-specific protection branch. Antibody protection is incomplete/selective in practice; no perfect specificity, complete digestion or complete protein removal is asserted. Exact assembled end transitions and full TruSeq PCR/read geometry are deferred. Lambda and RecJf are separate reactions with washes; RNA XRN1 is not substituted.'
READOUT='Source MiSeq/NextSeq single-end; read-start pileups report protected6mA positions computationally, not a chemical conversion of adenine.'
PREVIEWS=(("Source-qualified partial molecular route","library"),)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    current=g.input('dna','source genomic DNA/chromatin; genomic bases unassigned',K.MOLECULAR,state=molecule('source DNA',(role('genomic DNA; coordinates and sequence unassigned'),)),primary=True)
    reagent=g.input('shear_reagent','Covaris near150bp; RNaseA/T1 and AMPure source workup',K.CONDITION)
    product=g.product('shear_product','Covaris near150bp; RNaseA/T1 and AMPure source workup',K.MOLECULAR,state=molecule('source product classes',(role('fragmented RNA-depleted genomic DNA; exact residual RNA unquantified'),)))
    g.react('shear','Covaris near150bp; RNaseA/T1 and AMPure source workup',S.MAGIC,((current,R.TEMPLATE),(reagent,R.CONDITION)),(product,),reason='Source operation known; exact molecular assembly deliberately unreviewed.')
    current=product
    reagent=g.input('crosslink_reagent','Actual anti6mA antibody supply; overnight4C then254nmUV6x0.15J/cm2 source programme, no assumed complete binding/crosslinking',K.MOLECULAR)
    product=g.product('crosslink_product','Anti6mA incubation then254nmUV6x0.15J/cm2',K.MOLECULAR,state=molecule('source product classes',(role('genomic DNA population with UV antibody adducts at selected sites and uncrosslinked background'),)))
    crosslink_program=g.input('crosslink_program','Source antibody incubation overnight4C, split50ul UV aliquots and254nm6x0.15J/cm2',K.CONDITION)
    g.react('crosslink','Anti6mA incubation then254nmUV6x0.15J/cm2',S.CHEMICAL_CONVERSION,((current,R.TEMPLATE),(reagent,R.STRUCTURAL),(crosslink_program,R.CONDITION)),(product,))
    current=product
    reagent=g.input('capture_reagent','ProteinA bead binding and source sequential buffer washes',K.SELECTION)
    product=g.product('capture_product','ProteinA bead binding and source sequential buffer washes',K.MOLECULAR,state=molecule('source product classes',(role('bead-associated DNA-antibody complexes and possible nonspecific DNA'),)))
    other=g.product('capture_other','physically separated source flowthrough/washes; partition unquantified',K.WASTE)
    g.react('capture','ProteinA bead binding and source sequential buffer washes',S.AFFINITY_CAPTURE,((current,R.TEMPLATE),(reagent,R.SELECTION)),(product,other))
    g.discard(other,'Actual source isolation/washes; no perfect specificity asserted.')
    current=product
    reagent=g.input('endrepair_reagent','T4DNA polymerase12C30min then source bead washes',K.CATALYST)
    product=g.product('endrepair_product','T4DNA polymerase12C30min then source bead washes',K.MOLECULAR,state=molecule('source product classes',(role('end-repaired captured DNA; exact strand architecture unassembled'),)))
    g.react('endrepair','T4DNA polymerase12C30min then source bead washes',S.MAGIC,((current,R.TEMPLATE),(reagent,R.CATALYST)),(product,),reason='Source operation known; exact molecular assembly deliberately unreviewed.')
    current=product
    current=cleanup(g,'endrepair_wash','Actual source sequential bead washes after endrepair',current,current.state)
    reagent=g.input('pnk_reagent','T4PNK37C30min then source bead washes',K.CATALYST)
    product=g.product('pnk_product','T4PNK37C30min then source bead washes',K.MOLECULAR,state=molecule('source product classes',(role('phosphorylated captured DNA class'),)))
    g.react('pnk','T4PNK37C30min then source bead washes',S.MAGIC,((current,R.TEMPLATE),(reagent,R.CATALYST)),(product,),reason='Source operation known; exact molecular assembly deliberately unreviewed.')
    current=product
    current=cleanup(g,'pnk_wash','Actual source sequential bead washes after pnk',current,current.state)
    adapter=Scene()
    for name in ('PG160','PG161'):
        seq,mod5,mod3=SOURCE_OLIGOS[name]
        adapter.strand(name.lower(),[Segment(name,seq)],label=name+' reverse-adapter strand',mod5=mod5,mod3=mod3)
    current=source_ligation(g,'reverse','PG160/PG161 T4DNA ligation25C2h',current,adapter)
    current=cleanup(g,'reverse_wash','Actual source sequential bead washes after reverse',current,current.state)
    reagent=g.input('fill_reagent','Phi29 fill and source bead washes',K.CATALYST)
    product=g.product('fill_product','Phi29 fill and source bead washes',K.MOLECULAR,state=molecule('source product classes',(role('Phi29 fill products on captured DNA; no exact end product inferred'),)))
    g.react('fill','Phi29 fill and source bead washes',S.MAGIC,((current,R.TEMPLATE),(reagent,R.CATALYST)),(product,),reason='Source operation known; exact molecular assembly deliberately unreviewed.')
    current=product
    current=cleanup(g,'fill_wash','Actual source sequential bead washes after fill',current,current.state)
    reagent=g.input('pnk2_reagent','Second T4PNK then source bead washes',K.CATALYST)
    product=g.product('pnk2_product','Second T4PNK then source bead washes',K.MOLECULAR,state=molecule('source product classes',(role('rephosphorylated captured DNA population'),)))
    g.react('pnk2','Second T4PNK then source bead washes',S.MAGIC,((current,R.TEMPLATE),(reagent,R.CATALYST)),(product,),reason='Source operation known; exact molecular assembly deliberately unreviewed.')
    current=product
    current=cleanup(g,'pnk2_wash','Actual source sequential bead washes after pnk2',current,current.state)
    reagent=g.input('lambda_reagent','Lambda37C30min',K.CATALYST)
    product=g.product('lambda_product','Lambda37C30min',K.MOLECULAR,state=mixture(('protected DNA',(role('antibody-protected DNA class; no perfect protection asserted'),)),('Lambda digestion debris',(role('degraded DNA and other reaction fragments retained until source washes'),))))
    g.react('lambda','Lambda37C30min',S.DIGESTION,((current,R.TEMPLATE),(reagent,R.CATALYST)),(product,))
    current=product
    current=cleanup(g,'lambda_wash_product','Actual source bead washes after Lambda',current,molecule('recovered DNA',(role('washed protected DNA; retained/carryover debris unquantified'),)))
    reagent=g.input('recj_reagent','RecJf37C30min',K.CATALYST)
    product=g.product('recj_product','RecJf37C30min',K.MOLECULAR,state=mixture(('protected DNA',(role('protected DNA fraction; exact ends unassembled'),)),('RecJf digestion debris',(role('degraded DNA and other reaction products retained until source washes'),))))
    g.react('recj','RecJf37C30min',S.DIGESTION,((current,R.TEMPLATE),(reagent,R.CATALYST)),(product,))
    current=product
    current=cleanup(g,'recj_wash_product','Actual source bead washes after RecJf',current,molecule('recovered DNA',(role('washed protected DNA fraction; no perfect depletion assumed'),)))
    reagent=g.input('release_reagent','SDS elution and proteinaseK overnight65C followed ethanol precipitation',K.CONDITION)
    product=g.product('recovered','SDS elution and proteinaseK overnight65C followed ethanol precipitation',K.MOLECULAR,state=molecule('source product classes',(role('recovered protected ssDNA with possible residual peptide/adduct material; no perfect chemical reversal asserted'),)))
    g.react('release','SDS elution and proteinaseK overnight65C followed ethanol precipitation',S.MAGIC,((current,R.TEMPLATE),(reagent,R.CONDITION)),(product,),reason='Source operation known; exact molecular assembly deliberately unreviewed.')
    current=product
    primer_scene=Scene();primer_scene.strand('pg162',[Segment('PG162',SOURCE_OLIGOS['PG162'][0])],label='source PG162')
    primer=g.input('pg162','PG162 annealed to reverse adapter',K.OLIGO,state=primer_scene)
    phi=g.input('phi29','Phi29 DNApolymerase30C30min',K.CATALYST)
    extended=g.product('extended','PG162/Phi29 extension product',K.MOLECULAR,state=downstream_state('protected DNA-derived duplex extension product; exact terminus geometry unassembled'))
    g.react('extend','PG162 annealing and Phi29 extension30C30min',S.EXTENSION,((current,R.TEMPLATE),(primer,R.STRUCTURAL),(phi,R.CATALYST)),(extended,))
    adapter=Scene()
    for name in ('PG163','PG164'):
        adapter.strand(name.lower(),[Segment(name,SOURCE_OLIGOS[name][0])],label=name+' source forward-adapter strand')
    ligated=source_ligation(g,'forward_ligation','Forward adapter T4DNA ligase25C2h',extended,adapter)
    purified=cleanup(g,'forward_cleanup','Source1.7xAMPure purification',ligated,ligated.state)
    amplified=source_pcr(g,'pcr','Phusion/TruSeqPCR human18cycles; bacterial13cycle branch',purified)
    current=cleanup(g,'library','Source0.8xAMPure primer/adapter-dimer removal',amplified,amplified.state)
    g.mark_terminal(Terminal(current.id,Disposition.UNRESOLVED,'Full molecular end-transition assembly and source TruSeq PCR primer/run geometry remain unreviewed. No sequence-absence claim.'))
    # Independently reviewed published preparation inventory; literal inputs and
    # operations, never derived by enumerating graph construction or graph nodes.
    g.declare_reviewed_scope(('dna','shear_reagent','shear','crosslink_reagent','crosslink_program','crosslink','capture_reagent','capture','endrepair_reagent','endrepair','endrepair_wash_condition','endrepair_wash_step','pnk_reagent','pnk','pnk_wash_condition','pnk_wash_step','reverse_adapter','reverse_ligase','reverse','reverse_wash_condition','reverse_wash_step','fill_reagent','fill','fill_wash_condition','fill_wash_step','pnk2_reagent','pnk2','pnk2_wash_condition','pnk2_wash_step','lambda_reagent','lambda','lambda_wash_product_condition','lambda_wash_product_step','recj_reagent','recj','recj_wash_product_condition','recj_wash_product_step','release_reagent','release','pg162','phi29','extend','forward_ligation_adapter','forward_ligation_ligase','forward_ligation','forward_cleanup_condition','forward_cleanup_step','pcr_forward','pcr_reverse','pcr_enzyme','pcr','library_condition','library_step'),g.authority)
    g.add_source_item(SourceItem('pending_assembly','Unreviewed source/assembly boundary','Full molecular end-transition assembly and source TruSeq PCR primer/run geometry remain unreviewed. No sequence-absence claim.'))
    return g

# Published TableS1, source strings retain the terminal phosphate.
SOURCE_OLIGOS = {
    'PG160': ('TGGAATTCTCGGGTGCC','',''),
    'PG161': ('CCTTGGCACCCGAGAATTCCA','5′ phosphate',''),
    'PG162': ('GCCTTGGCACCCGAGAATTCCA','',''),
    'PG163': ('GTTCAGAGTTCTACAGTCCGACGATC','',''),
    'PG164': ('GATCGTCGGACTGTAGAACTCTGAACCACCA','',''),
}
def oligos():
    return [oligo(name,[Segment(name,seq)],five=(mod5+' 5′-' if mod5 else '5′-'),three='-3′'+mod3) for name,(seq,mod5,mod3) in SOURCE_OLIGOS.items()]

def downstream_state(label):
    return mixture(('source-derived DNA',(role(label),)),('possible carryover',(role('background DNA, residual templates and reaction components; exact recovery unquantified'),)))

def source_ligation(g,ident,label,template,adapter):
    donor=g.input(ident+'_adapter','source adapter strands, assembly not yet reconstructed',K.ADAPTER,state=adapter)
    ligase=g.input(ident+'_ligase',label,K.CATALYST)
    out=g.product(ident+'_product',label,K.MOLECULAR,state=downstream_state('source adapter-ligated DNA roles; exact junction assembly unreviewed'))
    geometry=AdapterLigationGeometry('unavailable','unavailable',(adapter,), 'Source-known adapter identity; substrate end and sealed junction assembly unreviewed, no source-absence claim.')
    g.react(ident,label,S.LIGATION,((template,R.TEMPLATE),(donor,R.STRUCTURAL),(ligase,R.CATALYST)),(out,),geometry=geometry)
    return out

def source_pcr(g,ident,label,template):
    forward=g.input(ident+'_forward','source first PCR primer; precise indexed order variant/run assignment unreviewed',K.OLIGO,state=molecule('first PCR primer',(role('source primer; exact study-specific identity unassigned'),)))
    reverse=g.input(ident+'_reverse','source second PCR primer; precise indexed order variant/run assignment unreviewed',K.OLIGO,state=molecule('second PCR primer',(role('source primer; exact study-specific identity unassigned'),)))
    polymerase=g.input(ident+'_enzyme',label,K.CATALYST)
    out=g.product(ident+'_product',label,K.LIBRARY,state=downstream_state('source amplicon library roles; primer landing and exact end geometry unassembled'))
    g.react(ident,label,S.PCR,((template,R.TEMPLATE),(forward,R.STRUCTURAL),(reverse,R.STRUCTURAL),(polymerase,R.CATALYST)),(out,),geometry=PCRGeometry.unavailable(template.state,reason='Source PCR known; exact selected primer identities and full library geometry unreviewed.'))
    return out
