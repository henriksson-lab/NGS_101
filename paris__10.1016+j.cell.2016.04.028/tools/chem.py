"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from protocol_specs import Definition, get
from interaction_capture import ReversibleDuplexCrosslink
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from protocol_graph import SequencingPlan
from pcr import PCRGeometry
from chemdraw import Construct, Segment
import illumina as il
from adapter_ligation import AdapterLigationGeometry
from circular import CircularizationGeometry
from rna_library_operations import molecule

DEFINITION = Definition(
    "PARIS — psoralen analysis of RNA interactions and structures", "10.1016/j.cell.2016.04.028",
    "Lu et al., Cell (2016)",
    "AMT photocrosslinking, two-dimensional gel enrichment and proximity ligation preserve the two partners of an RNA duplex in one chimera.",
    "custom", ("RNA A-derived cDNA", "RNA B-derived cDNA"),
    "The original 2016 supplement distinguishes HeLa random-hexamer RT from HEK293/mES adapter-primed RT, which inherits Flynn 2016 icSHAPE library preparation. The circular-cDNA junction still needs molecular reconstruction before primer placement can be computed.",
    "RNA",
    "The original 2016 Extended Experimental Procedures specify two library routes. The graph follows the HEK293/mES adapter-ligation route inherited from Flynn 2016; the separate HeLa random-hexamer primer is shown below. The additionally labelled 2018 oligos document a later implementation with different RT-primer random/index tracts. Exact circular-cDNA geometry is not yet reconstructed."
)


from protocol_specs import psoralen_workflow
from paris_chem import oligos as paris_oligos


def _source_workflow():
    return psoralen_workflow('AMT', 'cross-linked duplexes selected above the diagonal on a 2D gel')

SPEC = get('paris', definition=DEFINITION,
           workflow_fn=_source_workflow,
           oligos_fn=lambda: (
               '<b>2016 HeLa random-hexamer RT primer</b> — 5′-/5phos/DDDNNXXXXNNNNNNTACCCTTCGCTTCACACACAAG/iSp18/GGATCC/iSp18/TACTGAACCGCNNNNNN-3′ (D = A/G/T; XXXX = four-base sample barcode; first N6 = duplicate tag; terminal N6 = RNA priming)',
               '<b>2016 HEK293/mES inherited RT-primer-1 example</b> — 5′-/5phos/DDDNNAACCNNNNAGATCGGAAGAGCGTCGTGGA/iSp18/GGATCC/iSp18/TACTGAACCGC-3′ (Flynn 2016 icSHAPE reagent list; AACC is its printed example index)',
               '<b>Later implementation: Lu 2018 oligos</b> — these are version-labelled comparison reagents, not an assertion that the 2016 RT primers were identical.',
               *paris_oligos()),
           sequencing_unavailable='The custom sequencing-primer bases are published; exact placement remains unresolved only because the circular-cDNA junction is not yet represented.')


def _local_graph(d, final, primers, *, reagent: str,
                   selection_label: str):
    """Build a psoralen contact graph from source-specific reagent choices."""
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    topology = ReversibleDuplexCrosslink(reagent, selection_label)
    rna = graph.input("rna", "cellular RNA containing paired regions",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("cellular RNA", ("RNA molecules",)), primary=True)
    crosslinker = graph.input("crosslinker", reagent + " 0.5 mg/ml in PBS; UV365 nm, 30 min on ice (mmc1 Step 1)",
                              MaterialKind.CONDITION)
    crosslinked = graph.product("crosslinked", reagent + "-crosslinked RNA duplexes",
                                MaterialKind.MOLECULAR,
                                state=topology.duplex())
    graph.react("crosslink", "Photocrosslink paired RNA with " + reagent,
                StepKind.CROSS_LINKING,
                ((rna, InputRole.TEMPLATE), (crosslinker, InputRole.CONDITION)),
                (crosslinked,))
    s1 = graph.input("s1_nuclease", "S1 nuclease in urea/SDS lysate, room temperature 10 min (mmc1 Step 1)",
                     MaterialKind.CATALYST)
    extracted = graph.product("extracted", "S1-digested, proteinase K/TRIzol-purified crosslinked RNA",
                              MaterialKind.MOLECULAR, state=crosslinked.state)
    graph.react("s1_extract", "S1 digestion, proteinase K (50 °C 30 min) and TRIzol extraction",
                StepKind.DIGESTION,
                ((crosslinked, InputRole.TEMPLATE), (s1, InputRole.CATALYST)), (extracted,))
    fragmentation = graph.input("fragmentation_condition",
                                "ShortCut RNase III with MnCl2, 37 °C 20 min (mmc1 Step 2)",
                                MaterialKind.CATALYST)
    fragments = graph.product("fragments", "fragmented crosslinked RNA duplexes",
                              MaterialKind.MOLECULAR, state=crosslinked.state)
    graph.react("fragment", "Fragment the crosslinked RNA with RNase III", StepKind.FRAGMENTATION,
                ((extracted, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CATALYST)), (fragments,))
    selector = graph.input("selector", selection_label, MaterialKind.SELECTION)
    selected = graph.product("selected", "enriched crosslinked RNA duplexes",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("rejected", "nonselected RNA fragments", MaterialKind.WASTE)
    graph.react("select_duplexes", "Enrich crosslinked RNA duplexes",
                StepKind.SELECTION,
                ((fragments, InputRole.TEMPLATE), (selector, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "does not retain the selected crosslinked duplex population")
    ligase = graph.input("proximity_ligase", "T4 RNA ligase 1 (NEB M0437M), room temperature 16–20 h; boiled 2 min to stop (mmc1 Step 3)",
                         MaterialKind.CATALYST)
    chimera = graph.product("chimera", "proximity-ligated RNA chimera",
                            MaterialKind.MOLECULAR, state=topology.chimera())
    graph.react("proximity_ligation", "Proximity-ligate the paired RNA fragments",
                StepKind.PROXIMITY_LIGATION,
                ((selected, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (chimera,))
    reversal = graph.input("photoreversal", "UV254 nm, 15 min on ice (mmc1 Step 3)",
                           MaterialKind.CONDITION)
    reversed_chimera = graph.product("reversed_chimera", "photoreversed RNA chimera",
                                     MaterialKind.MOLECULAR, state=chimera.state)
    graph.react("photoreverse", "Photoreverse the psoralen crosslink",
                StepKind.CHEMICAL_CONVERSION,
                ((chimera, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (reversed_chimera,))
    reviewed = ["rna", "crosslinker", "crosslink", "s1_nuclease", "s1_extract", "extracted",
                "fragmentation_condition", "fragment",
                "selector", "select_duplexes", "proximity_ligase", "proximity_ligation",
                "photoreversal", "photoreverse"]
    adapter = graph.input("rna_adapter", "Flynn 2016 preadenylylated 3′-blocked RNA linker /5rApp/AGATCGGAAGAGCGGTTCAG/ (3′-biotin or 3′-ddC variant not assigned for PARIS)", MaterialKind.ADAPTER)
    adapter_ligase = graph.input("adapter_ligase", "T4 RNA ligase 1 high concentration, 15% PEG8000, 25 °C 3 h (Flynn step 26)", MaterialKind.CATALYST)
    tagged = graph.product("tagged_rna", "adapter-tagged PARIS RNA chimera", MaterialKind.MOLECULAR,
                           state=molecule("adapter-tagged chimera", ("RNA chimera", "RNA library adapter")))
    graph.react("ligate_adapter", "Ligate the RNA-library adapter (2016 HEK293/mES route)", StepKind.LIGATION,
                ((reversed_chimera, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (tagged,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The original supplement inherits Flynn 2016 step 26 onward; adapter chemistry is source-located, but this graph has not reconstructed the ligation-end Scene."))
    rna_gel = graph.input("rna_gel", "denaturing urea PAGE; PARIS keeps >40-nt inserts (mmc1 Step 3; Flynn steps 27–39)",
                          MaterialKind.SELECTION)
    sized_rna = graph.product("sized_tagged_rna", "size-selected adapter-tagged chimera",
                              MaterialKind.MOLECULAR, state=tagged.state)
    free_linker = graph.product("free_linker", "unligated RNA linker and short products", MaterialKind.WASTE)
    graph.react("rna_size_select", "Gel-purify the adapter-ligated RNA", StepKind.SELECTION,
                ((tagged, InputRole.TEMPLATE), (rna_gel, InputRole.SELECTION)), (sized_rna, free_linker))
    graph.discard(free_linker, "left in the excluded gel region")
    rt_primer = graph.input("library_rt_primer", "Flynn 2016 adapter-annealing RT-primer-1", MaterialKind.OLIGO)
    rt = graph.input("library_rt_enzyme", "SuperScript III (Flynn step 42)", MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "PARIS first-strand cDNA", MaterialKind.MOLECULAR,
                         state=molecule("PARIS cDNA", ("RT primer and chimera cDNA",)))
    graph.react("library_rt", "Reverse-transcribe the adapter-tagged RNA", StepKind.REVERSE_TRANSCRIPTION,
                ((sized_rna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (cdna,))
    # Flynn steps 43-50 (streptavidin capture of the 3'-biotin linker, RNase/biotin
    # elution) are not drawn: the PARIS supplement does not say which blocked
    # linker variant it used, so whether this capture applies is unresolved.
    cdna_gel = graph.input("cdna_gel", "denaturing PAGE, cDNA above unextended RT primer (Flynn steps 52–57)",
                           MaterialKind.SELECTION)
    sized_cdna = graph.product("sized_cdna", "gel-purified PARIS cDNA", MaterialKind.MOLECULAR, state=cdna.state)
    free_primer = graph.product("free_rt_primer", "unextended RT primer", MaterialKind.WASTE)
    graph.react("cdna_size_select", "Gel-purify cDNA away from unextended RT primer", StepKind.SELECTION,
                ((cdna, InputRole.TEMPLATE), (cdna_gel, InputRole.SELECTION)), (sized_cdna, free_primer))
    graph.discard(free_primer, "left in the excluded gel region")
    circle_ligase = graph.input("circle_ligase", "CircLigase II with MnCl2, 60 °C 2 h (Flynn step 58)", MaterialKind.CATALYST)
    circle = graph.product("cdna_circle", "circularized PARIS cDNA", MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("circularize_cdna", "Circularize the cDNA before library PCR", StepKind.CIRCULARIZATION,
                ((sized_cdna, InputRole.TEMPLATE), (circle_ligase, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.known("single-strand intramolecular ligation",
                    statement="Inherited Flynn 2016 icSHAPE step 58 circularizes purified cDNA with CircLigase II; exact circular junction is not yet reconstructed in this model."))
    reviewed.extend(("rna_gel", "rna_size_select", "cdna_gel", "cdna_size_select",
                     "rna_adapter", "adapter_ligase", "ligate_adapter", "library_rt_primer",
                     "library_rt_enzyme", "library_rt", "circle_ligase", "circularize_cdna"))
    # Original PARIS explicitly inherits Flynn 2016 step 26 onward, including
    # its step 60 Solexa primer pair. This assignment belongs to this route.
    p5 = graph.input("library_forward", "Flynn 2016 P5-Solexa PCR primer",
                     MaterialKind.OLIGO,
                     state=Construct([Segment("P5-Solexa", il.NEBNEXT_UNIVERSAL_PRIMER)]),
                     source="Flynn 2016 doi:10.1038/nprot.2016.011 reagent list and step 60")
    p3 = graph.input("library_reverse", "Flynn 2016 P3-Solexa PCR primer",
                     MaterialKind.OLIGO,
                     state=Construct([Segment("P7", il.P7),
                         Segment("P3-Solexa arm", il.PE_READ2_PRIMER)]),
                     source="Flynn 2016 doi:10.1038/nprot.2016.011 reagent list and step 60")
    polymerase = graph.input("library_polymerase", "Phusion HF PCR master mix, qPCR-monitored cycles (Flynn steps 60–61)",
                             MaterialKind.CATALYST)
    library = graph.product("library_library", "PARIS sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Amplify the inherited adapter-route library", StepKind.PCR,
                ((circle, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p3, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(circle.state,
                    reason="The original inherited primer pair is declared; its circular-template binding sites and complete amplicon have not been assembled in this model."))
    plan = (SequencingPlan(library=final, primers=tuple(primers)) if primers
            else SequencingPlan(limitation=d.sequencing))
    graph.sequence(library, plan)
    reviewed.extend(("library_forward", "library_reverse", "library_polymerase",
                     "library_pcr", "library_library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS, reagent='AMT psoralen', selection_label='two-dimensional-gel duplex selection')

SPEC.protocol_graph = protocol_graph
