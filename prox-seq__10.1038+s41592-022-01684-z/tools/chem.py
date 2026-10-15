from adapter_ligation import AdapterLigationGeometry
from chemdraw import Scene, Segment, Workflow
from types import SimpleNamespace
from proximity_topology import ProximityLigationTopology

from rna_cell_assays import Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, branches, custom_library, single, token, unknown_region
from rna_cell_assays import feature

from rna_cell_assays import chain

def proximity_workflow() -> Workflow:
    a = unknown_region("antibody A barcode", tag="cbc",
              molecular_feature=feature("protein_barcode_a", "feature_barcode", "fixed"))
    b = unknown_region("antibody B barcode", tag="cbc",
              molecular_feature=feature("protein_barcode_b", "feature_barcode", "fixed"))
    umi = unknown_region("PLA UMI", tag="umi",
                         molecular_feature=feature("pla_umi", "umi", "random"))
    return chain("cell surface", branches("DNA-conjugated antibodies", (
        ("probe A", (a, token("A connector half"))),
        ("probe B", (token("B connector half"), b, umi, token("poly(A)"))))), (
        ("Bind paired antibody–oligo probes to extracellular proteins",
         branches("two bound targets", (("protein A", (a, token("A connector half"))),
                                         ("protein B", (token("B connector half"), b, umi,
                                                        token("poly(A)"))))),
         "Probe identity is encoded in one fixed barcode on each antibody conjugate."),
        ("Hybridize the universal connector and ligate only proximal probes",
         single("119-base PLA reporter", (a, token("A connector half"),
                                           token("B connector half"), b, umi, token("poly(A)")),
                junctions=(("A connector half", "B connector half",
                            "connector-splinted proximity ligation"),)),
         "The universal connector splints the A/B nick but is not inserted. Only a "
         "ligated A–B product contains both primer sites and is amplifiable."),
        ("Co-capture polyadenylated PLA reporters and cellular mRNA",
         branches("same cell barcode", (("PLA branch", (a, b, umi)),
                                         ("RNA branch", (token("transcript cDNA"),)))),
         "Drop-seq, Smart-seq2 or 10x can supply the single-cell partition and RNA branch."),
        ("PCR each branch and add sequencing indexes",
         single("Prox-seq PLA library", (token("sequencing arm"), a, b, umi,
                                          token("indexed sequencing arm"))),
         "The PLA read reports both protein barcodes and a molecule identifier."),
    ))

DEFINITION = Definition("Prox-seq — proximity sequencing", "10.1038/s41592-022-01684-z",
        "Vistain et al., <i>Nature Methods</i> (2022)", "01_prox-seq.html",
        "Pairs of antibody-conjugated DNA probes ligate only when their extracellular protein targets are close; two protein barcodes and a UMI are recovered with the same cell's mRNA.",
        "The complete 119-base PLA reporter is amplified as a distinct library; the exact read/index layout follows the chosen Drop-seq, Smart-seq2 or 10x single-cell implementation.",
        "Prox-seq intentionally supports several scRNA-seq back ends. The schematic therefore records reporter structure and identifier roles without claiming one universal sequencing-primer geometry.")

d = DEFINITION
final = custom_library('Prox-seq PLA library', (unknown_region('protein barcode A', tag='cbc', molecular_feature=feature('protein_barcode_a', 'feature_barcode', 'fixed')), unknown_region('protein barcode B', tag='cbc', molecular_feature=feature('protein_barcode_b', 'feature_barcode', 'fixed')), unknown_region('PLA UMI', tag='umi', molecular_feature=feature('pla_umi', 'umi', 'random'))))
primers = ()
reads = None
OLIGOS = ()
workflow_fn = proximity_workflow
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
    SEQ_PRIMERS=primers, FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=reads, SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE=("" if primers else "Exact sequencing-primer sequences are not fixed by the defining protocol or were not exposed in its accessible primary text; read entry and cycle roles are stated explicitly."),
    workflow=workflow_fn, oligos=lambda: OLIGOS,
)

SPEC.CAVEAT = ("The schematic combines the invariant PLA reporter across published Drop-seq, "
               "Smart-seq2 and 10x implementations. Supplementary Tables 7–10 print "
               "branch-specific oligos and custom sequencing primers; the schematic "
               "does not choose one platform-specific library end.")
SPEC.SEQUENCING_UNAVAILABLE = (
    "The paper publishes custom sequencing primers for its droplet and plate branches, "
    "but a single primer-binding map cannot describe all Prox-seq back ends at once.")
SPEC.oligos = lambda: (
    "<b>Universal ligation connector</b> — 5′-TTTCACGACACGACACGATTTAGGTC-3′",
    "<b>Droplet Read 1</b> — 5′-GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC-3′",
    "<b>Droplet Read 2</b> — 5′-AGATGCATCGTCTCGTGGGCTCGG-3′",
    "<b>Droplet i7</b> — 5′-CCGAGCCCACGAGACGATGCATCT-3′",
    "<b>Plate Read 1</b> — 5′-TATTTGGCATCGTCTCGTGGGCTCGG-3′",
    "<b>Plate i5</b> — 5′-CCGAGCCCACGAGACGATGCCAAATA-3′",
    "<b>Plate i7</b> — 5′-ACTCTGCGTTGATACCACTGCTTCGGC-3′",
)
SPEC.TOPOLOGY_PANELS = ProximityLigationTopology().panels()
_base_workflow = SPEC.workflow


def connector_geometry():
    sc = Scene(); sc.strand("connector", [Segment("connector A half", "X" * 10, placeholder=True),
                                           Segment("connector B half", "X" * 10, placeholder=True)],
                                      label="universal connector")
    sc.anneal("probe ends", [Segment("probe-B 5-prime ligating arm", "X" * 10, placeholder=True),
                              Segment("probe-A 3-prime ligating arm", "X" * 10, placeholder=True)],
              to="connector", pair=("probe-B 5-prime ligating arm", "connector B half"),
              label="adjacent antibody-probe ends", mod5="p",
              unpaired=("probe-B 5-prime ligating arm", "probe-A 3-prime ligating arm",
                        "connector A half", "connector B half"))
    sc.ligation_end("probe ends", "probe-B 5-prime ligating arm", kind="splinted",
                    polarity="5-prime", owner="adapter", text="probe 5′ nick partner")
    published = Scene()
    published.strand("published connector", [Segment(
        "full connector oligo", "TTTCACGACACGACACGATTTAGGTC")],
        label="source-written 26-base connector")
    return AdapterLigationGeometry.known(
        "splinted", sc, published,
        statement="The defining paper prints the 26-base connector and says 20 reporter bases "
                  "hybridize across the ligation junction. It does not assign the 20 paired "
                  "bases to exact positions on each antibody-probe end; the paired span "
                  "therefore remains role-level while the complete connector is shown separately.")


def workflow():
    old = _base_workflow()
    wf = Workflow(old.initial)
    for reaction in old.reactions:
        if "universal connector" in reaction.action:
            wf.ligate_adapter(reaction.action, connector_geometry(), reaction.after.rows,
                              name=reaction.after.name, note=reaction.note)
        else:
            wf.react(reaction.action, reaction.after.rows,
                     name=reaction.after.name, note=reaction.note)
    return wf


SPEC.workflow = workflow

def prox_seq_graph(definition, final: Construct) -> ProtocolGraph:
    """Prox-seq graph with covalent reporter and matched RNA-library branches."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    cell = graph.input("single_cell", "single cell with surface proteins and mRNA",
                       MaterialKind.BIOLOGICAL, state=proximity_workflow().initial, primary=True)
    antibodies = graph.input("antibody_probe_pairs", "antibody-conjugated A/B probe pairs",
                             MaterialKind.SELECTION)
    probe_oligos = graph.input("probe_oligos", "A/B protein-barcode probe oligos",
                               MaterialKind.OLIGO)
    bound = graph.product("probe_bound_cell", "cell with bound A/B antibody–oligo probes",
                          MaterialKind.MOLECULAR,
                          state=branches("two bound targets", (
                              ("protein A", (unknown_region("antibody A barcode", tag="cbc"),
                                             token("A connector half"))),
                              ("protein B", (token("B connector half"),
                                             unknown_region("antibody B barcode", tag="cbc"),
                                             unknown_region("PLA UMI", tag="umi"),
                                             token("poly(A)"))))))
    unbound = graph.product("unbound_probes", "unbound antibody–oligo probes",
                            MaterialKind.WASTE)
    graph.react("bind_probes", "Bind paired antibody–oligo probes to surface proteins",
                StepKind.BINDING,
                ((cell, InputRole.TEMPLATE), (antibodies, InputRole.SELECTION),
                 (probe_oligos, InputRole.STRUCTURAL)), (bound, unbound))
    graph.discard(unbound, "removed before proximity-reporter construction")
    connector = graph.input("universal_connector", "universal 20-base connector",
                            MaterialKind.OLIGO)
    hybridized = graph.product("connector_hybridized_probes",
                               "connector-hybridized proximal probe pair",
                               MaterialKind.MOLECULAR, state=bound.state)
    graph.react("hybridize_connector", "Hybridize the universal connector",
                StepKind.HYBRIDIZATION,
                ((bound, InputRole.TEMPLATE),
                 (connector, InputRole.STRUCTURAL)), (hybridized,))
    ligase = graph.input("dna_ligase", "DNA ligase", MaterialKind.CATALYST)
    reporter_state = single("119-base PLA reporter", (
        unknown_region("antibody A barcode", tag="cbc"),
        token("A connector half"), token("B connector half"),
        unknown_region("antibody B barcode", tag="cbc"),
        unknown_region("PLA UMI", tag="umi"), token("poly(A)")),
        junctions=(("A connector half", "B connector half",
                    "connector-splinted proximity ligation"),))
    ligated = graph.product("ligated_reporter_cell", "cell bearing ligated PLA reporters",
                            MaterialKind.MOLECULAR, state=reporter_state)
    graph.react("proximity_ligation", "Ligate only connector-bridged proximal probes",
                StepKind.LIGATION,
                ((hybridized, InputRole.TEMPLATE),
                 (connector, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The sources establish a 20-base connector-hybridized ligation region but no universal probe sequence across protein targets."))
    capture = graph.input("oligodt_capture_primer",
                          "platform-specific oligo(dT) cell-capture primer",
                          MaterialKind.OLIGO, inferred=True)
    rtase = graph.input("reverse_transcriptase", "reverse transcriptase",
                        MaterialKind.CATALYST)
    cell_cdna = graph.product("cell_cdna", "co-captured PLA-reporter and transcript cDNA",
                              MaterialKind.MOLECULAR,
                              state=branches("same cell identity", (
                                  ("PLA", (token("PLA reporter cDNA"),)),
                                  ("RNA", (token("transcript cDNA"),)))))
    graph.react("reverse_transcription", "Co-capture PLA reporters and cellular mRNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated, InputRole.TEMPLATE), (capture, InputRole.STRUCTURAL),
                 (rtase, InputRole.CATALYST)), (cell_cdna,))
    pla_cdna = graph.product("pla_cdna", "PLA-reporter cDNA branch",
                             MaterialKind.MOLECULAR, state=reporter_state)
    expression_cdna = graph.product("expression_cdna", "matched cellular cDNA branch",
                                    MaterialKind.MOLECULAR,
                                    state=single("RNA branch", (token("transcript cDNA"),)))
    graph.react("split_cdna", "Split PLA and transcriptome cDNA branches",
                StepKind.SPLIT, ((cell_cdna, InputRole.TEMPLATE),),
                (pla_cdna, expression_cdna))
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    pla_forward = graph.input("pla_forward_primer", "PLA reporter PCR forward primer",
                              MaterialKind.OLIGO)
    pla_reverse = graph.input("pla_reverse_primer", "indexed PLA reporter PCR reverse primer",
                              MaterialKind.OLIGO)
    pla_library = graph.product("pla_library", "Prox-seq PLA library",
                                MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("pla_pcr", "PCR and index the PLA reporter branch", StepKind.PCR,
                ((pla_cdna, InputRole.TEMPLATE), (pla_forward, InputRole.STRUCTURAL),
                 (pla_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (pla_library,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Read/index geometry depends on the selected single-cell back end."))
    graph.sequence(pla_library, SequencingPlan(limitation=definition.endpoint))
    expression_forward = graph.input("expression_forward_primer",
                                     "platform-specific expression PCR forward primer",
                                     MaterialKind.OLIGO, inferred=True)
    expression_reverse = graph.input("expression_reverse_primer",
                                     "platform-specific indexed expression PCR reverse primer",
                                     MaterialKind.OLIGO, inferred=True)
    expression_final = custom_library("Prox-seq matched RNA library", (
        token("cell identity"), token("transcript cDNA", length=44)))
    expression_library = graph.product("expression_library", "matched RNA library",
                                       MaterialKind.LIBRARY, state=expression_final,
                                       inferred=True)
    graph.react("expression_pcr", "PCR and index the matched RNA branch", StepKind.PCR,
                ((expression_cdna, InputRole.TEMPLATE),
                 (expression_forward, InputRole.STRUCTURAL),
                 (expression_reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (expression_library,),
                geometry=PCRGeometry.unavailable(
                    expression_final,
                    reason="The method deliberately supports Drop-seq, Smart-seq2 and 10x back ends."))
    graph.sequence(expression_library, SequencingPlan(
        limitation="The matched RNA endpoint follows the chosen Drop-seq, Smart-seq2 or 10x implementation."))
    graph.declare_reviewed_scope((
        "single_cell", "antibody_probe_pairs", "probe_oligos", "bind_probes",
        "universal_connector", "hybridize_connector", "dna_ligase",
        "proximity_ligation", "oligodt_capture_primer", "reverse_transcriptase",
        "reverse_transcription", "split_cdna", "pcr_polymerase", "pla_forward_primer",
        "pla_reverse_primer", "pla_pcr", "pla_library", "expression_forward_primer",
        "expression_reverse_primer", "expression_pcr", "expression_library",
    ), source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = 'prox-seq'
    return prox_seq_graph(d, final)
SPEC.protocol_graph = protocol_graph
