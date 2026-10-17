from adapter_ligation import AdapterLigationGeometry
from chemdraw import Scene, Segment, Workflow
from types import SimpleNamespace
from proximity_topology import ProximityLigationTopology

from rna_cell_assays import Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, branches, custom_library, single, token, unknown_region
from rna_cell_assays import feature

from rna_cell_assays import chain
from protocol_graph import Disposition, SourceClaim, SourceDiscrepancy, Terminal

# Supplementary Table 10 (MOESM3 XLSX) prints one probe-A and one probe-B oligo per
# antibody barcode; the fixed parts below are identical in every printed pair, and the
# 8-nt "Barcode" column is the target-specific window.  Written 5'->3' as printed:
#   A: /5AzideN/CGCATTGCATCGTCTCGTGGGCTCGGCHHHHACHHHHACHHHNGCAG[barcode]GATCGCTAAATCGTG
#   B: /5Phos/TCGTGTCGTGTCTAAAGTCC[barcode]ACAAAAAAAAAAAAAAAAAAAAAAAAAAAAA/3AzideN/
# The azides couple to DBCO-antibodies (published Methods, antibody-DNA conjugation).
PROBE_A_5_SITE = "CGCATTGCATCGTCTCGTGGGCTCGGC"   # contains Table 7 U_fwd GCATCGTCTCGTGGGCTC
PROBE_A_UMI = "HHHHACHHHHACHHHN"
PROBE_A_SPACER = "GCAG"
PROBE_A_LINKER = "GATCG"
PROBE_A_ARM = "CTAAATCGTG"                      # 3' arm paired to the connector
PROBE_B_ARM = "TCGTGTCGTG"                      # 5'-phosphorylated arm paired to the connector
PROBE_B_LINKER = "TCTAAAGTCC"
PROBE_B_SPACER = "AC"
PROBE_B_POLYA = "A" * 29
BARCODE_LEN = 8
CONNECTOR = "TTTCACGACACGACACGATTTAGGTC"        # plate-based Methods ligase solution; Suppl. Methods
PROBE_A_LEN = (len(PROBE_A_5_SITE) + len(PROBE_A_UMI) + len(PROBE_A_SPACER) + BARCODE_LEN
               + len(PROBE_A_LINKER) + len(PROBE_A_ARM))
PROBE_B_LEN = (len(PROBE_B_ARM) + len(PROBE_B_LINKER) + BARCODE_LEN + len(PROBE_B_SPACER)
               + len(PROBE_B_POLYA))
PRODUCT_LEN = PROBE_A_LEN + PROBE_B_LEN
PROBE_A_TABLE = ("/5AzideN/" + PROBE_A_5_SITE + PROBE_A_UMI + PROBE_A_SPACER + "[8-nt barcode]"
                 + PROBE_A_LINKER + PROBE_A_ARM)
PROBE_B_TABLE = ("/5Phos/" + PROBE_B_ARM + PROBE_B_LINKER + "[8-nt barcode]" + PROBE_B_SPACER
                 + PROBE_B_POLYA + "/3AzideN/")
BARCODE_A = feature("protein_barcode_a", "feature_barcode", "fixed")
BARCODE_B = feature("protein_barcode_b", "feature_barcode", "fixed")
PLA_UMI = feature("pla_umi", "umi", "random")


def probe_a_parts():
    return (Segment("probe-A U_fwd site", PROBE_A_5_SITE),
            Segment("PLA UMI", PROBE_A_UMI, "umi", placeholder=True, feature=PLA_UMI),
            Segment("probe-A spacer", PROBE_A_SPACER),
            Segment("antibody A barcode", "N" * BARCODE_LEN, "cbc", placeholder=True,
                    feature=BARCODE_A),
            Segment("probe-A linker", PROBE_A_LINKER),
            Segment("A connector arm", PROBE_A_ARM))


def probe_b_parts():
    return (Segment("B connector arm", PROBE_B_ARM),
            Segment("probe-B linker", PROBE_B_LINKER),
            Segment("antibody B barcode", "N" * BARCODE_LEN, "cbc", placeholder=True,
                    feature=BARCODE_B),
            Segment("probe-B spacer", PROBE_B_SPACER),
            Segment("poly(A)", PROBE_B_POLYA))


def reporter_state():
    return single(f"{PRODUCT_LEN}-nt ligated PLA product (Table 10 probes)",
                  probe_a_parts() + probe_b_parts(),
                  junctions=(("A connector arm", "B connector arm",
                              "connector-splinted proximity ligation"),))


def bound_state():
    return branches("two bound targets", (("protein A probe", probe_a_parts()),
                                           ("protein B probe", probe_b_parts())))


def proximity_workflow() -> Workflow:
    return chain("cell surface", branches("DNA-conjugated antibodies", (
        ("probe A (5′ azide)", probe_a_parts()),
        ("probe B (5′ phosphate, 3′ azide)", probe_b_parts()))), (
        ("Bind paired antibody–oligo probes to extracellular proteins", bound_state(),
         "Probe identity is encoded in the 8-nt Table 10 barcode on each antibody conjugate; "
         "probe A also carries the HHHHACHHHHACHHHN UMI window and the U_fwd site."),
        ("Hybridize the universal connector and ligate only proximal probes",
         reporter_state(),
         "T4 DNA ligase seals probe A's 3′ arm to probe B's 5′-phosphorylated arm on the "
         "26-nt connector, which is a splint and is not inserted. Only a ligated A–B product "
         "carries both the U_fwd site and the 3′ poly(A)."),
        ("Co-capture polyadenylated PLA reporters and cellular mRNA",
         branches("same cell barcode", (("PLA branch", probe_a_parts() + probe_b_parts()),
                                         ("RNA branch", (token("transcript cDNA"),)))),
         "Drop-seq, Smart-seq2 or 10x can supply the single-cell partition and RNA branch."),
        ("PCR each branch and add sequencing indexes",
         single("Prox-seq PLA library", (token("platform cell-barcode / capture arm"),
                                          *reversed_payload(),
                                          token("U_fwd-side indexed arm"))),
         "Shown as the strand complementary to the probes (poly(A)-capture side on the left). "
         "The PLA read reports the UMI and both protein barcodes."),
    ))


def reversed_payload():
    """Role-level library payload in capture-primer orientation (complement strand)."""
    return (token("oligo(dT) / poly(A)", length=12),
            unknown_region("protein barcode B", tag="cbc", molecular_feature=BARCODE_B),
            token("ligated connector region", length=20),
            unknown_region("protein barcode A", tag="cbc", molecular_feature=BARCODE_A),
            unknown_region("PLA UMI", tag="umi", molecular_feature=PLA_UMI))

DEFINITION = Definition("Prox-seq — proximity sequencing", "10.1038/s41592-022-01684-z",
        "Vistain et al., <i>Nature Methods</i> (2022)", "01_prox-seq.html",
        "Pairs of antibody-conjugated DNA probes ligate only when their extracellular protein targets are close; two protein barcodes and a UMI are recovered with the same cell's mRNA.",
        "The ligated PLA reporter is amplified as a distinct library; the exact read/index layout follows the chosen Drop-seq, Smart-seq2 or 10x single-cell implementation.",
        "Prox-seq intentionally supports several scRNA-seq back ends. The schematic therefore records reporter structure and identifier roles without claiming one universal sequencing-primer geometry.")

d = DEFINITION
final = custom_library('Prox-seq PLA library', reversed_payload())
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

SPEC.CAVEAT = ("The probe oligos and connector are drawn from Supplementary Table 10 and "
               "the Methods; the 8-nt target barcodes are shown as N. The schematic combines "
               "the PLA reporter across published Drop-seq, Smart-seq2 and 10x implementations: "
               "Supplementary Tables 7–9 print branch-specific PCR and custom sequencing "
               "primers, and the schematic does not choose one platform-specific library end. "
               "Treating HHHHACHHHHACHHHN as the probe-A UMI window is derived from the "
               "Supplementary Methods statement that the UMI is read from probe A.")
SPEC.SEQUENCING_UNAVAILABLE = (
    "The paper publishes custom primers and distinct run recipes: Drop-seq Read 1/i7/Read 2 "
    "20/8/85 cycles, 10x PLA 28/8/75, 10x cDNA 28/10/10/90 (Read 1/i7/i5/Read 2), "
    "and plate PLA Read 1/i5/i7 75/8/8. A single primer-binding map cannot describe all "
    "Prox-seq back ends at once.")
SPEC.oligos = lambda: (
    f"<b>Probe A oligo (Table 10, PLA.seq.Bead.N.A)</b> — 5′-{PROBE_A_TABLE}-3′",
    f"<b>Probe B oligo (Table 10, PLA.Seq.Bead.N.B)</b> — 5′-{PROBE_B_TABLE}-3′",
    f"<b>Universal ligation connector</b> — 5′-{CONNECTOR}-3′",
    "<b>U_fwd PLA PCR primer (Tables 7–8)</b> — 5′-GCATCGTCTCGTGGGCTC-3′",
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
    """Table 10 probe ends on the published connector; Scene checks all 20 pairs."""
    sc = Scene()
    sc.strand("probe A end", [Segment("probe-A linker", PROBE_A_LINKER),
                              Segment("A connector arm", PROBE_A_ARM)],
              label="probe A 3′ end (3′-OH)")
    sc.anneal("connector", [Segment("connector 5′ overhang", CONNECTOR[:3]),
                            Segment("connector B site", CONNECTOR[3:13]),
                            Segment("connector A site", CONNECTOR[13:23]),
                            Segment("connector 3′ overhang", CONNECTOR[23:])],
              to="probe A end", pair=("connector A site", "A connector arm"),
              label="universal connector (26 nt)",
              unpaired=("connector 5′ overhang", "connector 3′ overhang"))
    sc.anneal("probe B end", [Segment("B connector arm", PROBE_B_ARM),
                              Segment("probe-B linker", PROBE_B_LINKER)],
              to="connector", pair=("B connector arm", "connector B site"),
              label="probe B 5′ end", mod5="p", above=True,
              unpaired=("probe-B linker",))
    sc.ligation_end("probe B end", "B connector arm", kind="splinted",
                    polarity="5-prime", owner="adapter", text="5′-phosphate nick partner")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement="Supplementary Table 10 probe A ends in CTAAATCGTG-3′ and probe B begins "
                  "5′-phosphate-TCGTGTCGTG; the 26-nt connector pairs with both 10-nt arms "
                  "(20 bases, as the article states) and leaves 3-nt single-stranded ends. "
                  "The arm assignment is computed from the printed sequences; the paper does "
                  "not draw it.")


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
    probe_oligos = graph.input("probe_oligos", "Supplementary Table 10 azide-modified A/B probe oligos (DBCO click to antibodies)",
                               MaterialKind.OLIGO)
    bound = graph.product("probe_bound_cell", "cell with bound A/B antibody–oligo probes",
                          MaterialKind.MOLECULAR, state=bound_state())
    unbound = graph.product("unbound_probes", "unbound antibody–oligo probes",
                            MaterialKind.WASTE)
    graph.react("bind_probes", "Bind paired antibody–oligo probes to surface proteins",
                StepKind.BINDING,
                ((cell, InputRole.TEMPLATE), (antibodies, InputRole.SELECTION),
                 (probe_oligos, InputRole.STRUCTURAL)), (bound, unbound))
    graph.discard(unbound, "removed by two or three 1% BSA/PBS washes before ligation (published Methods, sample preparation)")
    connector = graph.input("universal_connector", "universal 26-nt connector oligo (TTTCACGACACGACACGATTTAGGTC)",
                            MaterialKind.OLIGO,
                            state=Construct([Segment("universal connector", CONNECTOR)],
                                            name="universal connector"))
    hybridized = graph.product("connector_hybridized_probes",
                               "connector-hybridized proximal probe pair",
                               MaterialKind.MOLECULAR, state=bound.state)
    graph.react("hybridize_connector", "Hybridize the universal connector",
                StepKind.HYBRIDIZATION,
                ((bound, InputRole.TEMPLATE),
                 (connector, InputRole.STRUCTURAL)), (hybridized,))
    ligase = graph.input("dna_ligase", "T4 DNA ligase (NEB), 130 U/ml with 1 mM rATP, 37 °C", MaterialKind.CATALYST)
    reporter = reporter_state()
    ligated = graph.product("ligated_reporter_cell", "cell bearing ligated PLA reporters",
                            MaterialKind.MOLECULAR, state=reporter)
    unligated = graph.product("unligated_probes", "unligated bound probe A and probe B oligos",
                              MaterialKind.MOLECULAR, state=bound_state(), inferred=True)
    graph.react("proximity_ligation", "Ligate only connector-bridged proximal probes",
                StepKind.LIGATION,
                ((hybridized, InputRole.TEMPLATE),
                 (connector, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated, unligated),
                geometry=connector_geometry())
    graph.mark_terminal(Terminal("unligated_probes", Disposition.UNRESOLVED,
        "The article says unligated probes are 'automatically discarded during the library "
        "preparation step' (JATS Results, Prox-seq design). Probe B alone still carries poly(A) "
        "and probe A alone the U_fwd site; lack of exponential PCR and later SPRI cleanups are "
        "not modelled as a physical separation."))
    graph.annotate(SourceDiscrepancy(
        "ligated_reporter_cell", "Length of the ligated PLA product",
        (SourceClaim("Published article (PMC11289786), Results paragraph on probe design",
                     "\"The complete PLA product spans 119 bases, with 20 of those bases "
                     "hybridized to a connector.\""),
         SourceClaim("Supplementary Table 10 (MOESM3 XLSX), every PLA.seq.Bead A/B pair",
                     f"Probe A is {PROBE_A_LEN} nt and probe B {PROBE_B_LEN} nt as printed, "
                     f"so the ligated product is {PRODUCT_LEN} nt.")),
        f"The schematic draws the printed Table 10 sequences, so its product is {PRODUCT_LEN} nt.",
        "The lengths are computed from the transcribed Table 10 oligos; the 20 connector-paired "
        "bases agree with the article.",
        "The article does not say how its 119 bases were counted (for example, a shorter "
        "poly(A) or excluded terminal region); neither value is used for any primer landing."))
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
                             MaterialKind.MOLECULAR, state=reporter)
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
        "proximity_ligation", "unligated_probes", "oligodt_capture_primer", "reverse_transcriptase",
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
