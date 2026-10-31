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
    c('Every DNA molecular output continues or has explicit terminal fate',not any(m.kind in {K.MOLECULAR,K.LIBRARY} and m.id not in g.terminals and not any(port.material==m.id for step in g.steps.values() for port in step.inputs) for m in g.materials.values()))
    selected=FootprintDNA('ACCGC',frozenset({1,4}))
    c('Selected bisulfite chemistry is DNA C-to-U, preserving all positions',selected.chemical_dna,'AUCGU')
    c('Complementary firstDNA and reference-oriented descendant are distinct',(selected.first_copy,selected.reference_oriented_copy),('ACGAT','ATCGT'))
    c('PCR report does not alter original DNA chemical class',(selected.original_dna,selected.chemical_dna),('ACCGC','AUCGU'))
    c.raises('Cannot bisulfite-convert a supplied guanine as cytosine',lambda:FootprintDNA('ACCGC',frozenset({3})))
    c.raises('Cannot silently combine directlysis with S9.6 enrichment',lambda:graph_for_branch(True,True))
    c('Genomic bisulfite receives the actual spiked preparation',g.steps['bisulfite'].inputs[0].material,'spiked_genomic')
    c('Control physical class is retained through downstream population states',all(g.materials[k].state.control for k in ('spiked_genomic','bisulfite_product','amplified','pooled','library')),True)
    c('Control spike requires linearization then denaturation',g.steps['control_denaturation'].inputs[0].material,'linearized_control')
    c('Control is copresent rather than a separate terminal',('denatured_control' in g.terminals,tuple(p.material for p in g.steps['control_spike'].inputs)),(False,('restriction_cleanup','denatured_control')))
    c('Optional S9.6 route does not assert a source-specified post-IP spike schedule','control_spike' in graph_for_branch(False,True).steps,False)
    c('Separate direct-lysis paragraph does not acquire an invented control spike','control_spike' in graph_for_branch(True,False).steps,False)
    c('Source direct and optionalIP branches retain publication boundaries',all(graph_for_branch(a,b).publication_problems() for a,b in ((True,False),(False,True))),True)
    c.report()

def population(label):
    return mixture(('source DNA',(role(label),)),('other physical material',(role('other template strands, RNA/protein or reaction fragments; actual molecular recovery unquantified'),)))

def op(g,ident,label,kind,current,new_state,reagent_label,reagent_kind=K.CONDITION,reagent_role=R.CONDITION):
    reagent=g.input(ident+'_reagent',reagent_label,reagent_kind)
    product=g.product(ident+'_product',label,K.MOLECULAR,state=new_state)
    kw={'reason':'Source operation established; exact molecular reaction geometry and population fractions unassembled.'} if kind in {S.MAGIC,S.GENERIC} else {}
    if kind is S.GENERIC:kw['changes_molecule']=False
    g.react(ident,label,kind,((current,R.TEMPLATE),(reagent,reagent_role)),(product,),**kw)
    return product

PAGE_FILE='smrf-seq.html'
TITLE='SMRF-seq'
NOTES='01_smrf-seq.html'
SOURCE='<a href="https://doi.org/10.1016/j.jmb.2020.02.014">Defining article</a>; Associated accepted-manuscript PMC7669280 Methods122–148 and Results52–69; linked published supporting originals gated'
CAVEAT='Nondenaturing bisulfite selects exposed DNA C-to-U chemical footprints; PCR-derived reference-orientedC-to-T is not physical RNA mutation or loss. Original genomic DNA-extraction and optionalDRIP recipes are source-cited, not imported by family name. Historical PacBio001-143-835-08 and original primerTableS2 remain unreviewed/gated. Direct-lysis conversion is a distinct source branch withD5020, not interchangeable withD5030extracted-DNA branch.'
READOUT='Source PacBioRSII6h movies, selected long amplicons and SMRTbell templates; no syntheticUMI, samplebarcode or standard short-read geometry is inferred.'
PREVIEWS=(("Source-qualified molecular operation route","library"),)

@dataclass(frozen=True)
class FootprintDNA:
    original_dna:str
    converted:frozenset[int]
    def __post_init__(self):
        if not self.original_dna or set(self.original_dna)-set('ACGT') or not isinstance(self.converted,frozenset) or any(type(i)is not int or not 0<=i<len(self.original_dna) or self.original_dna[i]!='C' for i in self.converted):raise ValueError('explicit original DNA and selected actual cytosine positions required')
    @property
    def chemical_dna(self):return ''.join('U' if i in self.converted else b for i,b in enumerate(self.original_dna))
    @property
    def first_copy(self):return self.chemical_dna.translate(str.maketrans('ACGTU','TGCAA'))[::-1]
    @property
    def reference_oriented_copy(self):return revcomp(self.first_copy)

@dataclass(frozen=True)
class SMRFPopulation:
    stage:str
    control:bool=False
    def rows(self):
        if self.stage not in {'native','restricted','bisulfite','amplified','library'}:raise ValueError('unsupported SMRF-seq population')
        sc=Scene()
        sc.strand('dna',[role('original looped-out genomic strand with selected exposed C-to-U class; full desulfonation intermediates unassembled' if self.stage in {'bisulfite','amplified','library'} else 'original genomic DNA strands with native R-loop class; no coordinates assigned')],label='original genomic DNA')
        sc.strand('paired_dna',[role('other original genomic DNA including hybrid-paired strand; no ideal protection or conversion rate asserted')],label='other original DNA')
        sc.strand('rna',[role('original hybrid-associated RNA and possible residual/degraded RNA; no physical RNA deletion by DNA bisulfite/PCR assumed')],label='RNA-associated material')
        if self.stage in {'amplified','library'}:sc.strand('daughter',[role('selected PCR descendant duplexes; DNA-U copying yields reference-orientedT calls, no remaining uracil required in descendants')],label='daughter DNA, independent molecules')
        if self.control:sc.strand('plasmid_control',[role('source linearized/denatured plasmid spike and possible descendants; conversion, recovery and any copying unassembled')],label='source plasmid-control material')
        sc.strand('other',[role('other nucleic acids, protein/reagent fragments and possible background; recovery unquantified')],label='other physical material')
        return sc.rows()

def graph_for_branch(direct=False,enriched=False):
    if direct and enriched:raise ValueError('source direct-conversion branch is not the extracted-DNA S9.6 enrichment branch')
    g=ProtocolGraph(TITLE+' — '+('directlysis' if direct else 'extractedDNA')+(' S9.6-enriched' if enriched else ''),authority=SOURCE,representation='published_instance')
    current=g.input('cells','source NTERA-2 cells; extracted and direct-conversion branches are separate preparations',K.MOLECULAR,state=SMRFPopulation('native'),primary=True)
    if direct:
        current=cleanup(g,'pellet','Source approximately30000cells pelleted and resuspended in1mlDPBS',current,current.state)
        current=op(g,'direct_lysis','Source D5020ProcedureB M-Digestion buffer/proteinaseK25min37C',S.MAGIC,current,current.state,'source M-Digestion buffer and kit proteinaseK; exact lysis/peptide products unassembled')
        reagent='source D5020direct sodiumbisulfite37C3.5h, gentle rotation; downstream manufacturer procedure unreviewed'
    else:
        current=op(g,'dna_extraction','Source DRIP-protocol DNA extraction with proteinaseK shortened to2h',S.MAGIC,current,current.state,'mainMethods130cited standardDRIP original; exact lysis/extraction/cleanup version unreviewed')
        current=op(g,'restriction','Source target-flanking REs80–120U overnight37C; fragments no more15kb',S.DIGESTION,current,SMRFPopulation('restricted'),'source study-specific restriction enzymes; exact locus sets in gated original tables',K.CATALYST,R.CATALYST)
        current=cleanup(g,'restriction_cleanup','Source1xAMPure after restriction digestion',current,current.state)
        if enriched:
            antibody=g.input('s96','source antiRNA:DNAhybridS9.6 immunoprecipitation per directly citedDRIP original',K.SELECTION)
            selected=g.product('selected','source S9.6-enriched hybrid-containing DNA; no perfect specificity',K.MOLECULAR,state=current.state)
            other=g.product('other_ip_fraction','unselected IP fractions, exact further disposition unassembled',K.MOLECULAR,state=population('other DNA/RNA; no physical disappearance asserted'))
            g.react('s96_enrich','Source optionalS9.6 enrichment before bisulfite, not mandatory',S.AFFINITY_CAPTURE,((current,R.TEMPLATE),(antibody,R.SELECTION)),(selected,other))
            g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Detailed source-cited DRIP wash/elution/disposal remains unreviewed.'));current=selected
        reagent='source ZymoD5030conversion reagent37C2h, gentle rotation, NO prior DNA denaturation'
    if not direct and not enriched:
        control=g.input('control_plasmid','Source positive-control plasmid stock; actual identity, whole sequence and amount unassembled',K.MOLECULAR)
        linear=g.product('linearized_control','Source linearized control plasmid; actual restriction enzyme, junctions and end geometry unassembled',K.MOLECULAR)
        g.react('control_linearization','Source positive-control plasmid linearization; no enzyme or end recipe invented',S.MAGIC,((control,R.TEMPLATE),),(linear,),reason='Associated Methods138 explicitly specifies a linearized control; exact production chemistry unreviewed.')
        denatured=g.product('denatured_control','Source denatured linearized plasmid control; whole strand distribution/conditions unassembled',K.MOLECULAR)
        g.react('control_denaturation','Source denatured positive control; temperature, mechanism and duration unassembled',S.MAGIC,((linear,R.TEMPLATE),),(denatured,),reason='Associated Methods138 explicitly specifies denatured plasmid; no thermal or chemical denaturation recipe supplied.')
        spiked=g.product('spiked_genomic','Actual genomic preparation copresent with source linearized/denatured plasmid spike; amounts and sequence geometry unassembled',K.MOLECULAR,state=SMRFPopulation('restricted',control=True))
        g.react('control_spike','Source Methods138 positive-control plasmid spiked into genomic samples and treated accordingly',S.POOL,((current,R.TEMPLATE),(denatured,R.TEMPLATE)),(spiked,))
        current=spiked
    current=op(g,'bisulfite','Source non-denaturing DNA bisulfite; selected exposed C-to-U chemical footprints',S.CHEMICAL_CONVERSION,current,SMRFPopulation('bisulfite',control=not direct and not enriched),reagent)
    current=op(g,'desulfonate','Source kit desulfonation; exact historical conditions unreviewed',S.CHEMICAL_CONVERSION,current,current.state,'source D5020orD5030manufacturer desulfonation instruction; selected branch retained')
    current=cleanup(g,'bisulfite_recovery','Source kit sample recovery after desulfonation; detailed vendor operation unassembled',current,current.state)
    f=g.input('forward','source locus-specificTableS2forward primer, exact gated sequence unreviewed',K.OLIGO)
    r=g.input('reverse','source locus-specificTableS2reverse primer, exact gated sequence unreviewed',K.OLIGO)
    enzyme=g.input('phusion_u','source PhusionU F555S;0.2mMdNTPs,0.8µMprimers,30–35cycles; optional1Mbetaine',K.CATALYST)
    amplified=g.product('amplified','source2–4kb locus amplicons plus original/carryover material',K.MOLECULAR,state=SMRFPopulation('amplified',control=not direct and not enriched))
    g.react('pcr','Source uracil-tolerant PhusionU locusPCR; original DNA/RNA not physically erased',S.PCR,((current,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(enzyme,R.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(Scene(),reason='Exact gated locus-specific primers and template/strand assembly unreviewed.'))
    current=cleanup(g,'pcr_cleanup','Source1xAMPure after PCR',amplified,amplified.state)
    other=g.input('other_amplicons','source other nonoverlapping amplicons, less20products per run; equal amounts',K.MOLECULAR,state=population('independently prepared source target amplicons'))
    pooled=g.product('pooled','source pooled nonoverlapping amplicons; no synthetic sample index assumed',K.MOLECULAR,state=SMRFPopulation('amplified',control=not direct and not enriched))
    g.react('pool','Source equal-amount amplicon pool',S.POOL,((current,R.TEMPLATE),(other,R.TEMPLATE)),(pooled,))
    current=cleanup(g,'concentrate','Source1–2µg pooled PCRproducts concentrated with1xAMPure',pooled,pooled.state)
    kit=g.input('pacbio_kit','source PacBio2kbTemplatePreparation PN001-143-835-08; NO prior damage repair, all kit washes0.8xAMPure, ligation25C1h',K.CONDITION)
    library=g.product('library','source SMRTbell library; exact ends/adapters/ligation junctions and remaining kit operations unassembled',K.LIBRARY,state=SMRFPopulation('library',control=not direct and not enriched))
    g.react('smrtbell_boundary','Source historical SMRTbell construction with explicit reported modifications; no modern kit substitute',S.MAGIC,((current,R.TEMPLATE),(kit,R.CONDITION)),(library,),reason='Main reports historical protocol/version, omission and modified ligation/washes; exact ordered internal kit geometry unavailable without its original.')
    g.mark_terminal(Terminal(library.id,Disposition.UNRESOLVED,'Source PacBioRSII6h movies with QC sizing/quantitation; exact SMRTbell/run-primer geometry unassembled.'))
    g.add_source_item(SourceItem('gated_si','Original supplementalMethods/figures/tables, including target primerTableS2','PMC7669280 NIHMS1568681supplements1–13 canonical requests gated; original contents uninspected.'))
    g.add_source_item(SourceItem('pending_vendor','Exact PacBio001-143-835-08, ZymoD5030/D5020 and DRIP predecessor originals','MainMethods130/138/140/148; exact cited source editions/chemistry remain unreviewed.'))
    return g

def protocol_graph():return graph_for_branch()
