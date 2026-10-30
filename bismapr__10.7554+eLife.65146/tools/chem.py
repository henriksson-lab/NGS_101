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
    c('Private preview has no finished marker', 'data-protocol-complete=' not in render_page())
    c('Every DNA molecular output continues to an operation or explicit terminal',not any(m.kind in {K.MOLECULAR,K.LIBRARY} and m.id not in g.terminals and not any(port.material==m.id for step in g.steps.values() for port in step.inputs) for m in g.materials.values()))
    c('UDG products remain physical DNA/material rather than discarded output',all(g.materials[mid].kind is not K.WASTE for mid in g.steps['udg'].outputs))
    c('PCR directly follows UDG in source reaction envelope, before physical cleanup',g.steps['pcr'].inputs[0].material=='udg_product')
    c('Selected BisMapR frontend has no RNaseA digestion operation',not any('RNaseA' in step.label and 'WITHOUT' not in step.label and step.kind is S.DIGESTION for step in g.steps.values()))
    c.report()

def unassembled(label):
    return mixture(('source-derived DNA',(role(label),)),('possible carryover',(role('other original DNA/RNA, reaction fragments or protein carryover; quantities and exact partition unassigned'),)))

def operation(g,ident,label,kind,x,state,reagent,*,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
    stock=g.input(ident+'_stock',reagent,reagent_kind)
    out=g.product(ident+'_product',label,K.MOLECULAR,state=state)
    kw={'reason':'Source operation established; exact molecular geometry/mixture relation not computed.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(ident,label,kind,((x,R.TEMPLATE),(stock,reagent_role)),(out,),**kw)
    return out

def ligate(g,ident,label,x,adapter_label):
    adapter_scene=molecule(adapter_label,(role('source adapter; exact historical sequence/modifications and geometry unassembled'),))
    adapter=g.input(ident+'_adapter',adapter_label,K.ADAPTER,state=adapter_scene)
    enzyme=g.input(ident+'_enzyme','source T4DNA ligase',K.CATALYST)
    out=g.product(ident+'_product',label,K.MOLECULAR,state=unassembled('source adapter-ligated product roles; no exact sealed junction claimed'))
    geometry=AdapterLigationGeometry('unavailable','unavailable',(adapter_scene,), 'Historical exact end/junction assembly remains unreviewed; no source-absence claim.')
    g.react(ident,label,S.LIGATION,((x,R.TEMPLATE),(adapter,R.STRUCTURAL),(enzyme,R.CATALYST)),(out,),geometry=geometry)
    return out

def amplify(g,ident,label,x,primer_label):
    f=g.input(ident+'_forward',primer_label+' forward member; exact member unassigned',K.OLIGO,state=molecule('forward primer',(role('source primer identity unassembled'),)))
    r=g.input(ident+'_reverse',primer_label+' reverse member; exact member unassigned',K.OLIGO,state=molecule('reverse primer',(role('source primer identity unassembled'),)))
    enzyme=g.input(ident+'_enzyme','source Q5DNApolymerase',K.CATALYST)
    out=g.product(ident+'_product',label,K.LIBRARY,state=unassembled('source amplified library and residual original templates/cleaved species; no molecular deletion assigned to PCR'))
    g.react(ident,label,S.PCR,((x,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),(out,),geometry=PCRGeometry.unavailable(x.state,reason='Exact original adapters/primer/member/strand geometry unassembled.'))
    return out

PAGE_FILE='bismapr.html'
TITLE='BisMapR'
NOTES='01_bismapr.html'
SOURCE='<a href="https://doi.org/10.7554/eLife.65146">Defining article</a>; Published PMC7901872 Fig1/Results52–76 and Methods142–156; source-cited2019 MapR frontend inspected separately'
CAVEAT='Selected2021 BisMapR route omitsRNaseA from source-cited MapR stop buffer. Do not use original MapR RNaseA digestion or infer protocol-wide equality between2019 and2020 MapR versions. C-to-U on exposed DNA and dUTP second-strand marking create distinct physical species; UDG does not physically remove DNA from the tube. Exact enzymes/end cleavage/adapter geometry and complete uracil depletion remain unassembled.'
READOUT='Source NextSeq500PE38; Read1 mapping strand reports hybrid-template DNA. Computational strand separation is not physical molecular purification.'
PREVIEWS=(("Source-qualified operation graph","library"),)

@dataclass(frozen=True)
class HybridPopulation:
    stage: str
    def __post_init__(self):
        if self.stage not in {'native','cut','rna_processed','bisulfite','second_strand','udg'}:raise ValueError('unsupported selected source population')
    def rows(self):
        sc=Scene()
        sc.strand('template',[role('original hybrid-template genomic DNA; source coordinates and sequence unassigned')],label='original template DNA')
        if self.stage in {'native','cut','rna_processed','bisulfite','second_strand','udg'}:
            sc.strand('displaced',[role('original displaced genomic strand; exposed cytosines selected as C-to-U class' if self.stage in {'bisulfite','second_strand'} else 'UDG-treated displaced-DNA/AP-product class; fragmented state not computed' if self.stage=='udg' else 'original displaced genomic DNA; no U label assigned')],label='displaced DNA, independent strand')
        sc.strand('rna',[role('RNA-fragment population and possible residual hybrid RNA; no physical discard' if self.stage in {'rna_processed','second_strand','udg'} else 'RNA associated noncovalently with template DNA; not a covalent DNA extension')],label='RNA or RNA-derived products')
        if self.stage in {'second_strand','udg'}:
            sc.strand('new_dna',[role('new dUTP-containing DNA second strand' if self.stage=='second_strand' else 'UDG-treated dU-DNA/AP-product class; exact secondary fragmentation unassembled')],label='new DNA, independent polymer')
        if self.stage=='udg':sc.strand('free_bases',[role('excised uracil and other reaction products remain until actual cleanup')],label='reaction products, not sequenced base calls')
        sc.strand('other',[role('other genomic DNA/protein and possible background or digestion fragments; recovery unquantified')],label='other physical material')
        return sc.rows()

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE,representation='published_instance')
    current=g.input('cells','source mESCs,5million; native R-loop population',K.MOLECULAR,state=HybridPopulation('native'),primary=True)
    current=cleanup(g,'cell_wash','Source-cited MapR cell wash/centrifugation',current,current.state)
    current=operation(g,'immobilize','Source-cited washed ConAbead immobilization',S.BINDING,current,current.state,'source-cited ConAbeads',reagent_kind=K.SELECTION,reagent_role=R.SELECTION)
    current=cleanup(g,'magnetic_transfer','Source-cited immobilized-cell magnetic transfer',current,current.state)
    current=operation(g,'permeabilize','Source-cited digitonin permeabilization',S.GENERIC,current,current.state,'source-cited Dig-wash buffer')
    current=operation(g,'target','Source-cited RHdelta-MNase binding',S.MAGIC,current,current.state,'source RNaseH-inactive/MNase fusion; precise2019/2020 condition equivalence unreviewed',reagent_kind=K.CATALYST,reagent_role=R.CATALYST)
    current=cleanup(g,'target_wash','Source-cited target wash before activation',current,current.state)
    calcium=g.input('calcium','source calcium activation of already introduced fusion',K.CONDITION)
    fragmented=g.product('mnase_product','source fusion-targeted cleavage population before physical release',K.MOLECULAR,state=HybridPopulation('cut'))
    g.react('mnase','Activate already bound MNase with source calcium; do not add another enzyme aliquot',S.DIGESTION,((current,R.TEMPLATE),(g.materials['target_stock'],R.CATALYST),(calcium,R.CONDITION)),(fragmented,))
    current=fragmented
    current=operation(g,'stop_release','Source-cited STOP/release WITHOUT RNaseA; retained hybrid RNA required',S.MAGIC,current,current.state,'source STOP mixture with RNaseA explicitly omitted; spike-in details/source revision unassembled')
    current=cleanup(g,'supernatant','Source-cited MapR release centrifugation/supernatant collection',current,current.state)
    current=operation(g,'proteolyze','Source-cited MapR SDS/proteinaseK handling; no RNaseA substitution',S.DEPROTEINIZATION,current,current.state,'source-cited ProteinaseK; exact2019/2020 quantities unresolved',reagent_kind=K.CATALYST,reagent_role=R.CATALYST)
    current=cleanup(g,'organic','Source-cited MapR DNA extraction, retaining RNA needed for hybrid protection; no perfect purity asserted',current,current.state)
    current=cleanup(g,'ethanol','Source-cited MapR DNA recovery; full version-specific conditions unassembled',current,current.state)
    current=operation(g,'bisulfite','Non-denaturing EZDNA Methylation-GoldCTreagent25C3h; selected exposed C-to-U',S.CHEMICAL_CONVERSION,current,HybridPopulation('bisulfite'),'source10µlsample+10µlwater+130µlCTreagent')
    current=operation(g,'desulfonate','Source kit chemical desulfonation; full reactive end intermediates unassembled',S.CHEMICAL_CONVERSION,current,current.state,'source EZDNA Methylation-Gold desulfonation reagents/conditions; exact kit protocol unreviewed')
    current=cleanup(g,'desulfonation','Source kit desulfonation and column purification;20µlM-Elution',current,current.state)
    current=operation(g,'second_strand','RNaseH/dUTP second-strand synthesis16C1h; RNA-derived fragments retained',S.MAGIC,current,HybridPopulation('second_strand'),'source NEBNextUltraIIDirectional E7760 second-strand buffer/enzyme mix; exact kit polymerase roles unassembled')
    current=cleanup(g,'second_cleanup','Source1.8xAMPure cleanup after second-strand synthesis',current,current.state)
    current=operation(g,'end_repair','Source EnzymaticsEnd-RepairMix',S.MAGIC,current,current.state,'source End-RepairMix; exact end transitions unassembled')
    current=operation(g,'a_tail','Source Klenowexo-minus A-tailing',S.MAGIC,current,current.state,'source Klenowexo-minus and nucleotides')
    current=cleanup(g,'preligation','Source MinElute purification before ligation',current,current.state)
    current=ligate(g,'adapter_ligation','Source NEBE7600 Illumina adapters/T4DNA ligation',current,'source NEBE7600 adapter; historical hairpin/ends unassembled')
    current=cleanup(g,'size_select','SourceAMPure fragments greater150bp',current,current.state)
    current=operation(g,'udg','Source UDG included in PCR mix; uracil excision is not physical disposal',S.DIGESTION,current,HybridPopulation('udg'),'source Enzymatics uracilDNA glycosylase; adapter hairpin processing also source-reported, exact products unassembled',reagent_kind=K.CATALYST,reagent_role=R.CATALYST)
    current=amplify(g,'pcr','Source Q5 dual-index PCR; original/UDG-treated DNA and reaction fragments remain',current,'source NEBE7600 dual-index primers')
    current=cleanup(g,'library','Source MinElute post-PCR purification',current,current.state)
    g.mark_terminal(Terminal(current.id,Disposition.UNRESOLVED,'Source NextSeq500PE38; exact molecular/read-primer/index geometry not assembled.'))
    g.add_source_item(SourceItem('pending_originals','MapR2020, named vendor documents and exact UDG/adapter/primer assembly','Methods144 cites Yan2020 DOI10.1002/cpmb.113 and Yan2019; ZymoD5005/NEBE7760/E7600 vendor originals unreviewed.'))
    return g
