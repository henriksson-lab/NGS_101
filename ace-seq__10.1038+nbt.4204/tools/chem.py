from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from chemdraw import Segment, oligo
from base_conversion import BasePath
from hybrid_base_protocols import make_spec
from hybrid_base_protocols import Definition, _conversion_scene, _part, _role_library, _role_library_scene


def ace_path(base: str) -> BasePath:
    """ACE-seq glucosylation, APOBEC3A deamination and PCR readout."""
    if base == "5hmC":
        return BasePath(("5hmC", "5gmC", "5gmC", "C"))
    if base == "5mC":
        return BasePath(("5mC", "5mC", "U", "T"))
    if base == "C":
        return BasePath(("C", "C", "U", "T"))
    raise ValueError("ACE-seq input must be C, 5mC or 5hmC")

ACE_QC_PRIMERS = (
    ("lambda control forward", "AAAGGCCACTATCAGGCAGCTTTGT"),
    ("lambda control reverse", "GAAGAAATGTTTCCTGTAATGGAAGATG"),
    ("T4 control forward", "AGCTTATTTGAAATTAACTGGGTTG"),
    ("T4 control reverse", "TCTTCTCCAGAATAAATGAATTTACCT"),
)


def _ace_oligos():
    return tuple(oligo(name, (Segment(name, seq),)) for name, seq in ACE_QC_PRIMERS)



DEFINITION = Definition("ACE-seq — APOBEC-coupled epigenetic sequencing",
        "10.1038/nbt.4204", "Schutsky et al., Nature Biotechnology (2018)",
        "β-GT protects 5hmC before APOBEC3A deaminates C and 5mC, giving a bisulfite-free readout in which only original 5hmC remains C.",
        ("ACE-treated insert",), "The defining experiment used single-end NextSeq sequencing of Accel Methyl-NGS libraries.", "01_ace_seq.html",
        caveat="The supplement publishes four spike-in QC primers, shown below. They amplify the lambda and T4 conversion controls and are not library adapters. Accel Methyl-NGS adapter and sequencing-primer oligos remain proprietary and are shown by role.")

FINAL_LIBRARY, SEQ_PRIMERS = _role_library(
    DEFINITION.title + " final library", DEFINITION.payload), ()
SPEC = make_spec(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS, source_detail="", oligos=_ace_oligos)

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
    denaturing = graph.input('denaturing_condition', 'heat-denaturing condition', MaterialKind.CONDITION)
    apobec = graph.input('apobec3a', 'APOBEC3A', MaterialKind.CATALYST)
    treated = graph.product('ace_treated_dna', 'APOBEC3A-deaminated protected DNA', MaterialKind.MOLECULAR, state=_conversion_scene(tuple(((base, ace_path(base).states[-1]) for base in ('C', '5mC', '5hmC'))), 'ACE readout'))
    graph.react('apobec_deamination', 'Denature and deaminate C and 5mC', StepKind.CHEMICAL_CONVERSION, ((protected, InputRole.TEMPLATE), (denaturing, InputRole.CONDITION), (apobec, InputRole.CATALYST)), (treated,))
    adapters = graph.input('accel_adapters', 'Accel Methyl-NGS adapter reagents', MaterialKind.ADAPTER, inferred=True)
    ligase = graph.input('adapter_installation_enzymes', 'Accel adapter-installation enzyme mix', MaterialKind.CATALYST, inferred=True)
    ligated = graph.product('adapter_tagged_dna', 'Accel-adapter-tagged ACE-treated DNA', MaterialKind.MOLECULAR, state=_role_library_scene('ACE-seq adapter-tagged template', ('ACE-treated insert',), ligated=False), inferred=True)
    unavailable = AdapterLigationGeometry.unavailable('Accel Methyl-NGS adapter-end geometry is proprietary and not disclosed by the defining source.')
    graph.react('adapter_installation', 'Install Accel Methyl-NGS adapters', StepKind.LIGATION, ((treated, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (ligated,), geometry=unavailable)
    conversion_input = ligated
    reviewed = ['genomic_dna', 'shearing_condition', 'fragment', 'conversion_controls', 'beta_gt', 'udp_glucose', 'glucosylate_5hmc', 'denaturing_condition', 'apobec3a', 'apobec_deamination', 'accel_adapters', 'adapter_installation_enzymes', 'adapter_installation']
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
