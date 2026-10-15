from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from batch_ngs import seg, truseq_library
from base_conversion import BasePath
from hybrid_base_protocols import make_spec
from hybrid_base_protocols import Definition, _conversion_scene, _part, _role_library_scene


def tab_path(base: str) -> BasePath:
    """TAB-seq protection/oxidation, bisulfite conversion and PCR readout."""
    if base == "5hmC":
        return BasePath(("5hmC", "5gmC", "5gmC", "C"))
    if base == "5mC":
        return BasePath(("5mC", "5caC", "U", "T"))
    if base == "C":
        return BasePath(("C", "C", "U", "T"))
    raise ValueError("TAB-seq input must be C, 5mC or 5hmC")

DEFINITION = Definition("TAB-seq — TET-assisted bisulfite sequencing",
        "10.1016/j.cell.2012.04.027", "Yu et al., Cell (2012)",
        "β-GT protects 5hmC as 5gmC, then mTet1 oxidizes 5mC to 5caC; after bisulfite, only original 5hmC reads as C.",
        ("TAB-treated insert",), "TruSeq reads interrogate the converted insert; an i7 read identifies the library.", "01_tab_seq.html", truseq=True)

FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(
    [seg(DEFINITION.payload[0], "X" * 34, placeholder=True)],
    DEFINITION.title + " final library", dual_index=False)
SPEC = make_spec(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS, source_detail=" Detailed chemistry was checked against the authors’ later protocol.")

def _conversion_graph(definition, final: Construct, primers: tuple) -> ProtocolGraph:
    """Typed source-provenance graph for the three cytosine-conversion methods."""
    source = f'{definition.source}; doi:{definition.doi}'
    graph = ProtocolGraph(definition.title, authority=source)
    initial_state = _conversion_scene((('original C', 'C'), ('original 5mC', '5mC'), ('original 5hmC', '5hmC')), 'genomic cytosines')
    dna = graph.input('genomic_dna', 'genomic DNA', MaterialKind.BIOLOGICAL, state=initial_state, primary=True)
    shear = graph.input('shearing_condition', 'mechanical shearing condition', MaterialKind.CONDITION)
    fragmented = graph.product('fragmented_dna', 'fragmented genomic DNA', MaterialKind.MOLECULAR, state=initial_state)
    graph.react('fragment', 'Shear genomic DNA', StepKind.FRAGMENTATION, ((dna, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragmented,))
    controls = graph.input('conversion_controls', 'C/5mC/5hmC control molecules', MaterialKind.MOLECULAR)
    beta_gt = graph.input('beta_gt', 'T4 beta-glucosyltransferase', MaterialKind.CATALYST)
    udp_glucose = graph.input('udp_glucose', 'UDP-glucose', MaterialKind.CONDITION)
    protected = graph.product('glucosylated_dna', '5hmC-glucosylated genomic DNA', MaterialKind.MOLECULAR, state=_conversion_scene((('original C', 'C'), ('original 5mC', '5mC'), ('original 5hmC', '5gmC')), 'after beta-GT'))
    graph.react('glucosylate_5hmc', 'Glucosylate 5hmC', StepKind.CHEMICAL_CONVERSION, ((fragmented, InputRole.TEMPLATE), (controls, InputRole.SPIKE_IN), (beta_gt, InputRole.CATALYST), (udp_glucose, InputRole.CONDITION)), (protected,))
    mtet = graph.input('mtet1', 'mTet1 dioxygenase', MaterialKind.CATALYST)
    tet_conditions = graph.input('tet_cofactors', 'mTet1 oxidation cofactors', MaterialKind.CONDITION)
    oxidized = graph.product('tet_oxidized_dna', 'mTet1-oxidized protected DNA', MaterialKind.MOLECULAR, state=_conversion_scene((('original C', 'C'), ('original 5mC', '5caC'), ('original 5hmC', '5gmC')), 'after mTet1'))
    graph.react('oxidize_5mc', 'Oxidize 5mC to 5caC', StepKind.CHEMICAL_CONVERSION, ((protected, InputRole.TEMPLATE), (mtet, InputRole.CATALYST), (tet_conditions, InputRole.CONDITION)), (oxidized,))
    conversion_input = oxidized
    reviewed = ['genomic_dna', 'shearing_condition', 'fragment', 'conversion_controls', 'beta_gt', 'udp_glucose', 'glucosylate_5hmc', 'mtet1', 'tet_cofactors', 'oxidize_5mc']
    repair = graph.input('end_repair_enzymes', 'WGBS end-repair enzyme mix', MaterialKind.CATALYST)
    repaired = graph.product('end_repaired_dna', 'end-repaired TAB-seq DNA', MaterialKind.MOLECULAR, state=conversion_input.state)
    graph.react('end_repair', 'End-repair treated DNA', StepKind.END_REPAIR, ((conversion_input, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input('a_tailing_enzyme', 'WGBS dA-tailing enzyme', MaterialKind.CATALYST)
    tailed = graph.product('da_tailed_dna', 'dA-tailed TAB-seq DNA', MaterialKind.MOLECULAR, state=repaired.state)
    graph.react('da_tail', 'Add one 3-prime dA to each end', StepKind.TAILING, ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input('wgbs_adapters', 'WGBS library adapters', MaterialKind.ADAPTER, inferred=True)
    ligase = graph.input('dna_ligase', 'DNA ligase', MaterialKind.CATALYST)
    ligated = graph.product('ligated_precursor', 'adapter-ligated TAB-seq DNA', MaterialKind.MOLECULAR, state=tailed.state, inferred=True)
    graph.react('adapter_ligation', 'Ligate WGBS adapters', StepKind.LIGATION, ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (ligated,), geometry=AdapterLigationGeometry.unavailable('The defining source names its WGBS preparation but does not print the adapter-end oligos.'))
    conversion_input = ligated
    reviewed.extend(('end_repair_enzymes', 'end_repair', 'a_tailing_enzyme', 'da_tail', 'wgbs_adapters', 'dna_ligase', 'adapter_ligation'))
    bisulfite = graph.input('bisulfite_reagent', 'bisulfite conversion reagent', MaterialKind.CONDITION)
    converted = graph.product('bisulfite_converted_dna', 'bisulfite-converted library template', MaterialKind.MOLECULAR, state=_conversion_scene(tuple(((base, tab_path(base).states[-1]) for base in ('C', '5mC', '5hmC'))), 'converted readout'))
    graph.react('bisulfite_conversion', 'Bisulfite-convert the treated DNA', StepKind.CHEMICAL_CONVERSION, ((conversion_input, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)), (converted,))
    conversion_input = converted
    reviewed.extend(('bisulfite_reagent', 'bisulfite_conversion'))
    forward = graph.input('library_forward_primer', 'library PCR forward primer', MaterialKind.OLIGO, inferred=not bool(primers))
    reverse = graph.input('library_reverse_primer', 'indexed library PCR reverse primer', MaterialKind.OLIGO, inferred=not bool(primers))
    polymerase = graph.input('library_polymerase', 'library amplification polymerase', MaterialKind.CATALYST)
    library = graph.product('sequencing_library', definition.title + ' sequencing library', MaterialKind.LIBRARY, state=final, inferred=not bool(primers))
    graph.react('library_pcr', 'Index and amplify the sequencing library', StepKind.PCR, ((conversion_input, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,), geometry=PCRGeometry.unavailable(conversion_input.state if isinstance(conversion_input.state, Construct) else final, reason='The source identifies the library PCR but does not print one immutable selected primer pair.'))
    plan = SequencingPlan(final, tuple(primers), tuple((p.role for p in primers))) if primers else SequencingPlan(limitation=definition.caveat or definition.sequencing)
    graph.sequence(library, plan)
    reviewed.extend(('library_forward_primer', 'library_reverse_primer', 'library_polymerase', 'library_pcr', 'sequencing_library'))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    return _conversion_graph(d, final, primers)
SPEC.protocol_graph = protocol_graph
