from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from batch_ngs import seg, truseq_library
import illumina as il
from endprep import ThreePrimeOverhangAdapter, cohesive_ligation_geometry, repair_and_dA_tail
from base_conversion import BasePath
from hybrid_base_protocols import make_spec
from hybrid_base_protocols import Definition, _conversion_scene, _part, _role_library_scene


def oxbs_path(base: str) -> BasePath:
    """Oxidative-bisulfite state and PCR readout for one genomic cytosine."""
    if base == "5hmC":
        return BasePath(("5hmC", "5fC", "U", "T"))
    if base == "5mC":
        return BasePath(("5mC", "5mC", "C"))
    if base == "C":
        return BasePath(("C", "U", "T"))
    raise ValueError("oxBS input must be C, 5mC or 5hmC")

TRUSEQ_Y_END = ThreePrimeOverhangAdapter(
    "TruSeq Y adapter", il.STEM + "T", il.STEM_COMPLEMENT)


DEFINITION = Definition("oxBS-seq — oxidative bisulfite sequencing",
        "10.1126/science.1220671", "Booth et al., Science (2012)",
        "Potassium perruthenate converts 5hmC to 5fC before bisulfite, so 5hmC and unmodified C read as T while 5mC alone reads as C.",
        ("oxidative-bisulfite-converted insert",), "TruSeq reads interrogate the converted insert; an i7 read identifies the library.", "01_oxbs_seq.html", truseq=True,
        caveat=("Library implementation: the sonication/TruSeq shell shown here follows the authors’ "
                "<a href='https://doi.org/10.1038/nprot.2013.115'>2013 detailed protocol</a>; "
                "the defining 2012 Science study demonstrates an oxRRBS branch. Its original supplementary "
                "library details remain under review. <b>Published caption discrepancy:</b> 2013 Supplementary "
                "Fig. 1b describes DNA made with d5fCTP and 5fC conversion, but one clause names 5hmC "
                "as that substrate. The model uses 5fC → U, supported by the surrounding caption and "
                "main protocol. The inconsistent clause is an apparent caption error; its intended wording is unresolved."))

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
    repair = graph.input('end_repair_enzymes', 'end-repair enzyme mix', MaterialKind.CATALYST)
    repaired = graph.product('end_repaired_dna', 'end-repaired genomic fragments', MaterialKind.MOLECULAR, state=_role_library_scene('end-repaired oxBS insert', ('unconverted genomic insert',), ligated=False))
    graph.react('end_repair', 'End-repair fragmented DNA', StepKind.END_REPAIR, ((fragmented, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input('a_tailing_enzyme', 'dA-tailing enzyme', MaterialKind.CATALYST)
    da_scene = _role_library_scene('dA-tailed oxBS insert', ('unconverted genomic insert',), ligated=False)
    tailed = graph.product('da_tailed_dna', 'dA-tailed oxBS fragments', MaterialKind.MOLECULAR, state=da_scene)
    graph.react('da_tail', 'Add one 3-prime dA to each end', StepKind.TAILING, ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input('truseq_adapters', 'TruSeq Y adapters', MaterialKind.ADAPTER)
    ligase = graph.input('dna_ligase', 'DNA ligase', MaterialKind.CATALYST)
    ligated = graph.product('ligated_precursor', 'TruSeq-adapter-ligated oxBS precursor', MaterialKind.MOLECULAR, state=_role_library_scene('adapter-ligated oxBS precursor', ('unconverted genomic insert',)))
    end_prep = repair_and_dA_tail(Construct([_part('unconverted genomic insert')]))
    geometry = cohesive_ligation_geometry(end_prep, TRUSEQ_Y_END, statement='The protocol specifies TruSeq adapters: their public 3′-T end pairs with the dA-tailed insert.')
    graph.react('adapter_ligation', 'Ligate TruSeq adapters', StepKind.LIGATION, ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (ligated,), geometry=geometry)
    oxidant = graph.input('perruthenate', 'potassium perruthenate', MaterialKind.CONDITION)
    oxidized = graph.product('oxidized_dna', 'selectively oxidized adapter-ligated DNA', MaterialKind.MOLECULAR, state=_conversion_scene((('original C', 'C'), ('original 5mC', '5mC'), ('original 5hmC', '5fC')), 'after oxidation'))
    graph.react('oxidize_5hmc', 'Oxidize 5hmC selectively', StepKind.CHEMICAL_CONVERSION, ((ligated, InputRole.TEMPLATE), (oxidant, InputRole.CONDITION)), (oxidized,))
    conversion_input = oxidized
    reviewed = ['genomic_dna', 'shearing_condition', 'fragment', 'end_repair_enzymes', 'end_repair', 'a_tailing_enzyme', 'da_tail', 'truseq_adapters', 'dna_ligase', 'adapter_ligation', 'perruthenate', 'oxidize_5hmc']
    bisulfite = graph.input('bisulfite_reagent', 'bisulfite conversion reagent', MaterialKind.CONDITION)
    converted = graph.product('bisulfite_converted_dna', 'bisulfite-converted library template', MaterialKind.MOLECULAR, state=_conversion_scene(tuple(((base, oxbs_path(base).states[-1]) for base in ('C', '5mC', '5hmC'))), 'converted readout'))
    graph.react('bisulfite_conversion', 'Bisulfite-convert the treated DNA', StepKind.CHEMICAL_CONVERSION, ((conversion_input, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)), (converted,))
    conversion_input = converted
    reviewed.extend(('bisulfite_reagent', 'bisulfite_conversion'))
    forward = graph.input('library_forward_primer', 'library PCR forward primer', MaterialKind.OLIGO, inferred=not bool(primers))
    reverse = graph.input('library_reverse_primer', 'indexed library PCR reverse primer', MaterialKind.OLIGO, inferred=not bool(primers))
    polymerase = graph.input('library_polymerase', 'library amplification polymerase', MaterialKind.CATALYST)
    library = graph.product('sequencing_library', definition.title + ' sequencing library', MaterialKind.LIBRARY, state=final, inferred=not bool(primers))
    graph.react('library_pcr', 'Index and amplify the sequencing library', StepKind.PCR, ((conversion_input, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL), (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,), geometry=PCRGeometry.unavailable(conversion_input.state if isinstance(conversion_input.state, Construct) else final, reason='The selected library-kit PCR primer pair and index member have not been declared from verified source evidence in this model. Published control-DNA primers are a separate reagent set.'))
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
