from types import SimpleNamespace
from chemdraw import Construct, Scene
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind
from chemdraw import Segment, oligo
from base_conversion import BasePath
from protocol_model import role


def _conversion_scene(states, label):
    """Representative base sites, not an asserted genomic sequence or duplex."""
    scene = Scene()
    for origin, state in states:
        scene.strand(origin, [Segment(state, state, placeholder=True, length_bp=1)],
                     label=f"{label}: {origin} → {state}", at=0)
    return scene


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



DEFINITION = SimpleNamespace(
        title="ACE-seq — APOBEC-coupled epigenetic sequencing",
        doi="10.1038/nbt.4204", source="Schutsky et al., Nature Biotechnology (2018)",
        summary="β-GT protects 5hmC before APOBEC3A deaminates C and 5mC, giving a bisulfite-free readout in which only original 5hmC remains C.",
        sequencing="The defining experiment used single-end NextSeq sequencing of Accel Methyl-NGS libraries.", notes="01_ace_seq.html",
        caveat="INFERRED — exact Accel-NGS Methyl-Seq adapter and sequencing-primer assignments remain unresolved in this model; vendor originals are deferred. The supplement publishes four spike-in QC primers, shown below; these assay conversion controls and do not construct the library. The main Methods specify single-end NextSeq 500 using a 150-cycle kit; the actual acquisition mask is unassigned. Deamination leaves U at original C/5mC positions; T is the PCR readout.")

FINAL_LIBRARY = Construct([
    role("Accel left adapter, assignment unresolved", inferred=True),
    role("PCR copy of ACE-treated insert", inferred=True),
    role("Accel right adapter, assignment unresolved", inferred=True),
], name="ACE-seq final library (role-level model)")
SEQ_PRIMERS = ()
SPEC = SimpleNamespace(
    TITLE=DEFINITION.title, NOTES=DEFINITION.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.source}</a>. Publisher Supplementary Protocol 2 (MOESM6).',
    SUMMARY=DEFINITION.summary, CAVEAT=DEFINITION.caveat,
    FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=SEQ_PRIMERS,
    FINAL_CAPTION="INFERRED — role-level Accel library; exact adapters and insert length are unassigned.",
    SEQUENCING_INTRO=DEFINITION.sequencing, READ_LENGTHS=None,
    SEQUENCING_ENDING=DEFINITION.sequencing,
    SEQUENCING_UNAVAILABLE="Exact kit oligos, sequencing-primer binding and acquisition mask remain unassigned pending vendor-original review.",
    oligos=_ace_oligos,
)

def _conversion_graph(definition, final: Construct, primers: tuple) -> ProtocolGraph:
    """ACE-seq declarations; delegated kit steps do not imply known end geometry."""
    source = f'{definition.source}; doi:{definition.doi}; Supplementary Protocol 2'
    graph = ProtocolGraph(definition.title, authority=source)
    initial_state = _conversion_scene(tuple((b, b) for b in ('C', '5mC', '5hmC')), 'genomic sites')
    dna = graph.input('genomic_dna', 'genomic DNA', MaterialKind.BIOLOGICAL, state=initial_state, primary=True)
    controls = graph.input('conversion_controls', 'methylated lambda and T4-hmC controls', MaterialKind.MOLECULAR)
    graph.annotate(SourceDiscrepancy(
        'conversion_controls', 'Spike-in quantity differs between sources',
        (SourceClaim('Main Methods, PMC6453757 JATS line118', '0.5% total by mass'),
         SourceClaim('Publisher MOESM6 Supplementary Protocol2 Materials, line1307', '0.5% (w/w) of each control')),
        model_choice='The graph records the presence of methylated lambda and T4-hmC controls without selecting a quantity.',
        evidence='The defining article and its publisher supplement state different total versus per-control quantities.',
        uncertainty='The amount used for each experiment remains unresolved.'))
    shear = graph.input('shearing_condition', 'Covaris shearing to approximately 300 bp', MaterialKind.CONDITION)
    fragmented = graph.product('fragmented_dna', 'sheared sample and controls', MaterialKind.MOLECULAR, state=initial_state)
    graph.react('fragment', 'Shear genomic DNA and spike-in controls', StepKind.FRAGMENTATION,
                ((dna, InputRole.TEMPLATE), (controls, InputRole.SPIKE_IN), (shear, InputRole.CONDITION)), (fragmented,))
    beta_gt = graph.input('beta_gt', 'T4 beta-glucosyltransferase', MaterialKind.CATALYST)
    udp_glucose = graph.input('udp_glucose', 'UDP-glucose', MaterialKind.CONDITION)
    protected_state = _conversion_scene(tuple((b, ace_path(b).states[1]) for b in ('C', '5mC', '5hmC')), 'after beta-GT')
    protected = graph.product('glucosylated_dna', '5hmC-glucosylated DNA', MaterialKind.MOLECULAR, state=protected_state)
    graph.react('glucosylate_5hmc', 'Glucosylate 5hmC', StepKind.CHEMICAL_CONVERSION,
                ((fragmented, InputRole.TEMPLATE), (beta_gt, InputRole.CATALYST), (udp_glucose, InputRole.CONDITION)), (protected,))
    denaturing = graph.input('denaturing_condition', 'DMSO, heat denaturation and snap freezing', MaterialKind.CONDITION)
    denatured = graph.product('denatured_dna', 'denatured protected DNA (representative base sites)', MaterialKind.MOLECULAR, state=protected_state)
    graph.react('denature', 'Denature and snap freeze before A3A addition', StepKind.THERMAL_TREATMENT,
                ((protected, InputRole.TEMPLATE), (denaturing, InputRole.CONDITION)), (denatured,))
    apobec = graph.input('apobec3a', 'APOBEC3A', MaterialKind.CATALYST)
    ramp = graph.input('a3a_conditions', 'MES/Tween buffer and 4–50 °C temperature ramp', MaterialKind.CONDITION)
    treated_state = _conversion_scene(tuple((b, ace_path(b).states[-2]) for b in ('C', '5mC', '5hmC')), 'deamination before library construction')
    treated = graph.product('ace_treated_dna', 'C/5mC-deaminated protected DNA', MaterialKind.MOLECULAR, state=treated_state)
    graph.react('apobec_deamination', 'Deaminate C and 5mC to U', StepKind.CHEMICAL_CONVERSION,
                ((denatured, InputRole.TEMPLATE), (apobec, InputRole.CATALYST), (ramp, InputRole.CONDITION)), (treated,))
    cleanup = graph.input('conversion_cleanup_kit', 'Zymo Oligo Clean and Concentrator', MaterialKind.SELECTION)
    recovered = graph.product('recovered_converted_dna', 'purified converted DNA', MaterialKind.MOLECULAR, state=treated_state)
    graph.react('conversion_cleanup', 'Purify converted DNA before library construction (recovery boundary)', StepKind.GENERIC,
                ((treated, InputRole.TEMPLATE), (cleanup, InputRole.SELECTION)), (recovered,), changes_molecule=False)
    kit = graph.input('accel_kit', 'Accel-NGS Methyl-Seq 30024 and indexing kit 36024', MaterialKind.CONDITION)
    library = graph.product('sequencing_library', 'purified ACE-seq sequencing library', MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react('accel_library_construction', 'INFERRED — delegated Accel library construction, qPCR-chosen amplification and final purification', StepKind.MAGIC,
                ((recovered, InputRole.TEMPLATE), (kit, InputRole.CONDITION)), (library,),
                reason='Supplementary Protocol 2 delegates library construction to the named Swift procedure. Exact intermediate strands, adapter ends and PCR primers await the applicable vendor originals; no ligation-ready or PCR-template geometry is assigned here.')
    graph.sequence(library, SequencingPlan(acquired_roles=('Read 1',), limitation=SPEC.SEQUENCING_UNAVAILABLE + ' Single-end NextSeq 500 with a 150-cycle kit does not establish a 150-cycle acquisition mask.'))
    graph.declare_reviewed_scope((
        'genomic_dna', 'conversion_controls', 'shearing_condition', 'fragment',
        'beta_gt', 'udp_glucose', 'glucosylate_5hmc', 'denaturing_condition',
        'denature', 'apobec3a', 'a3a_conditions', 'apobec_deamination',
        'conversion_cleanup_kit', 'conversion_cleanup', 'accel_kit',
        'accel_library_construction', 'sequencing_library',
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    return _conversion_graph(d, final, primers)
SPEC.protocol_graph = protocol_graph
