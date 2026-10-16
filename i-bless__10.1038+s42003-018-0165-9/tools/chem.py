"""Protocol-owned source, oligos, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
from topology import HairpinTopology

BLESS_PROXIMAL = "TACTACCTCGAGAGTTACGCTAGGGATAACAGGGTAATATAGTTTTTTTCTATATTACCCTGTTATCCCTAGCGTAACTCTCGAGGTAGTA"

BLESS_DISTAL = "CGTCGTCTCGAGAGTTACGCTAGGGATAACAGGGTAATATAGTTTTTTTCTATATTACCCTGTTATCCCTAGCGTAACTCTCGAGACGACG"

BLESS_PCR_PROXIMAL = "CCCTAGCGTAACTCTCGAGGTAGTA"

BLESS_PCR_DISTAL = "CTAGCGTAACTCTCGAGACGACG"

BLESS_PROXIMAL_READ_TAG = "TCGAGGTAGTA"

BLESS_DISTAL_READ_TAG = "TCGAGACGACG"

TITLE = 'i-BLESS — immobilized Breaks Labeling, Enrichment on Streptavidin and Sequencing'
NOTES = '01_i_bless.html'
SOURCE = 'Biernacka et al., Communications Biology (2018); <a href="https://doi.org/10.1038/s42003-018-0165-9">doi:10.1038/s42003-018-0165-9</a>.'
SUMMARY = 'Agarose beads immobilize yeast cells through lysis and direct break labelling; the otherwise BLESS-like proximal/distal hairpins retain which sequenced end is the native DSB.'
CAVEAT = "The published i-BLESS hairpins and insert-PCR primers are shown exactly. Final ThruPLEX geometry and study index assignment remain unassembled. The September 2016 guide is a conditional kit reference: it offers single- and dual-index products, so no single index configuration is assigned to every study library."
FINAL_CAPTION = 'Role-level final library; the exact outer kit assembly and orientation of each sequenced insert remain open model assignments.'
READ_LENGTHS = {'Read 1': 75, 'Read 2': 75}
ACQUIRED_ROLES = ('Read 1', 'Read 2')
SEQUENCING_INTRO = 'The defining study reports paired 75-base reads on HiSeq 2500/4000.'
SEQUENCING_ENDING = SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE = 'Exact study-specific primer-binding geometry, index configuration and technical index-cycle mask remain unassigned; computational tag removal does not establish these physical details.'
FINAL_LIBRARY = Construct([role('unassembled ThruPLEX end', inferred=True), role('BLESS-tagged genomic insert'), role('unassembled opposite ThruPLEX end', inferred=True)], name='role-level i-BLESS library')
SEQ_PRIMERS = ()

def oligos():
    return (f'proximal hairpin — 5′-/5Phos/{BLESS_PROXIMAL[:45]}/iBiodT/{BLESS_PROXIMAL[46:]}-3′', f'distal hairpin — 5′-/5Phos/{BLESS_DISTAL}-3′', f'proximal PCR primer — 5′-{BLESS_PCR_PROXIMAL}-3′', f'distal PCR primer — 5′-{BLESS_PCR_DISTAL}-3′')

def protocol_graph():
    source = SOURCE
    graph = ProtocolGraph(TITLE, authority=source)
    sample = graph.input("sample", "yeast cells",
                         MaterialKind.BIOLOGICAL,
                         state=molecule("DSB-bearing sample", (role("genomic DNA with native DSBs"),)),
                         primary=True)
    current = sample
    deproteinization = graph.input("deproteinization_condition",
                                   "agarose immobilization, spheroplasting and deproteinization condition",
                                   MaterialKind.CONDITION)
    protected = graph.product("protected_dna", "agarose-protected deproteinated yeast DNA",
                              MaterialKind.MOLECULAR, state=sample.state)
    graph.react("immobilize", "Immobilize and deproteinate cells in agarose beads",
                StepKind.DEPROTEINIZATION,
                ((sample, InputRole.TEMPLATE),
                 (deproteinization, InputRole.CONDITION)), (protected,))
    current = protected
    repair = graph.input("proximal_repair", "native-break blunting and phosphorylation enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired_breaks", "blunt 5-prime-phosphorylated native breaks",
                             MaterialKind.MOLECULAR, state=current.state)
    graph.react("repair_breaks", "Blunt and 5-prime-phosphorylate native breaks",
                StepKind.END_REPAIR,
                ((current, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    proximal = graph.input("proximal_adapter", "biotinylated proximal hairpin",
                           MaterialKind.ADAPTER, state=Construct([Segment("published proximal hairpin", BLESS_PROXIMAL)]), source="Supplementary Table 4; 5′ phosphate and internal biotin-dT shown in oligo list")
    proximal_ligase = graph.input("proximal_ligase", "proximal-hairpin DNA ligase",
                                  MaterialKind.CATALYST)
    labelled = graph.product("proximal_labelled", "proximal-hairpin-labelled native breaks",
                             MaterialKind.MOLECULAR,
                             state=molecule("proximal-labelled break", (
                                 role("biotinylated proximal hairpin"),
                                 role("break-adjacent genomic DNA")),
                                 junctions=(("biotinylated proximal hairpin",
                                             "break-adjacent genomic DNA",
                                             "proximal ligation"),)))
    graph.react("ligate_proximal", "Ligate the proximal hairpin in situ",
                StepKind.LIGATION,
                ((repaired, InputRole.TEMPLATE), (proximal, InputRole.STRUCTURAL),
                 (proximal_ligase, InputRole.CATALYST)), (labelled,),
                geometry=reported_adapter_end(
                    "BLESS proximal hairpin", kind="blunt",
                    statement="The repaired native DSB is blunt and 5′ phosphorylated before proximal-hairpin ligation."))
    sonication = graph.input("sonication", "Covaris S220 sonication", MaterialKind.CONDITION)
    first = graph.product("initial_sonicated", "initially sonicated encapsulated DNA", MaterialKind.MOLECULAR, state=labelled.state)
    graph.react("initial_sonication", "Initially sonicate encapsulated DNA", StepKind.FRAGMENTATION,
                ((labelled, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (first,))
    recovery = graph.input("agarose_recovery", "Zymoclean Large Fragment DNA Recovery Kit", MaterialKind.CONDITION)
    extracted = graph.product("recovered_dna", "DNA recovered from agarose", MaterialKind.MOLECULAR, state=first.state)
    graph.react("recover_dna", "Recover total DNA from agarose", StepKind.EXTRACTION,
                ((first, InputRole.TEMPLATE), (recovery, InputRole.CONDITION)), (extracted,))
    fragments = graph.product("fragments", "proximal-labelled fragments after second sonication (~400 bp)", MaterialKind.MOLECULAR, state=extracted.state)
    graph.react("fragment", "Sonicate recovered DNA again to about 400 bp", StepKind.FRAGMENTATION,
                ((extracted, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)), (fragments,))
    streptavidin = graph.input("streptavidin", "streptavidin capture beads",
                               MaterialKind.SELECTION)
    captured = graph.product("captured", "captured proximal-labelled break fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    rejected = graph.product("uncaptured", "unlabelled genomic fragments", MaterialKind.WASTE)
    graph.react("capture", "Capture proximal-labelled molecules on streptavidin",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "does not retain the biotinylated native-break adapter")
    distal_repair = graph.input("distal_repair", "distal-end blunting enzyme mix",
                                MaterialKind.CATALYST)
    distal_ready = graph.product("distal_ready", "blunt captured distal fragment ends",
                                 MaterialKind.MOLECULAR, state=captured.state)
    graph.react("repair_distal", "Blunt the fragmentation-created distal end",
                StepKind.END_REPAIR,
                ((captured, InputRole.TEMPLATE), (distal_repair, InputRole.CATALYST)),
                (distal_ready,))
    distal = graph.input("distal_adapter", "distinct distal hairpin",
                         MaterialKind.ADAPTER, state=Construct([Segment("published distal hairpin", BLESS_DISTAL)]), source="Supplementary Table 4; 5′ phosphate shown in oligo list")
    distal_ligase = graph.input("distal_ligase", "distal-hairpin DNA ligase",
                                MaterialKind.CATALYST)
    dual = graph.product("dual", "dual-hairpin break fragments", MaterialKind.MOLECULAR,
                         state=molecule("dual-hairpin break fragment", (
                             role("proximal hairpin"), role("break-to-random genomic DNA"),
                             role("distal hairpin")),
                             junctions=(("proximal hairpin", "break-to-random genomic DNA",
                                         "native-break ligation"),
                                        ("break-to-random genomic DNA", "distal hairpin",
                                         "distal ligation"))))
    graph.react("ligate_distal", "Ligate the distinct distal hairpin",
                StepKind.LIGATION,
                ((distal_ready, InputRole.TEMPLATE), (distal, InputRole.STRUCTURAL),
                 (distal_ligase, InputRole.CATALYST)), (dual,),
                geometry=reported_adapter_end(
                    "BLESS distal hairpin", kind="blunt",
                    statement="The fragmentation-created distal end is blunted before the distinct distal hairpin is ligated."))
    iscei = graph.input("iscei", "I-SceI hairpin-release nuclease", MaterialKind.CATALYST)
    released = graph.product("released", "linearized dual-tag break inserts",
                             MaterialKind.MOLECULAR, state=dual.state)
    graph.react("linearize", "Linearize the captured hairpin product with I-SceI",
                StepKind.DIGESTION,
                ((dual, InputRole.TEMPLATE), (iscei, InputRole.CATALYST)), (released,))
    insert_fwd = graph.input("insert_forward", "published PCR primer 1", MaterialKind.OLIGO,
                             state=Construct([Segment("PCR primer 1", BLESS_PCR_PROXIMAL)]))
    insert_rev = graph.input("insert_reverse", "published PCR primer 2", MaterialKind.OLIGO,
                             state=Construct([Segment("PCR primer 2", BLESS_PCR_DISTAL)]))
    insert_pol = graph.input("insert_polymerase", "insert-PCR polymerase (study-specific choice unassigned)", MaterialKind.CATALYST)
    amplicon = graph.product("insert_amplicon", "amplified BLESS-tagged insert", MaterialKind.MOLECULAR, state=released.state)
    graph.react("insert_pcr", "PCR-amplify I-SceI-linearized inserts before XhoI", StepKind.PCR,
                ((released, InputRole.TEMPLATE), (insert_fwd, InputRole.STRUCTURAL),
                 (insert_rev, InputRole.STRUCTURAL), (insert_pol, InputRole.CATALYST)), (amplicon,),
                geometry=PCRGeometry.unavailable(Construct([role("unassembled I-SceI-linearized template")]),
                    reason="Published whole PCR primers are positive inputs; full linearized-template assembly and exact landing sites remain unassembled locally."))
    cleanup = graph.input("insert_cleanup", "PCR-product purification (implementation unassigned)", MaterialKind.CONDITION)
    clean = graph.product("purified_insert", "purified insert-PCR product", MaterialKind.MOLECULAR, state=amplicon.state)
    cleanup_other = graph.product("insert_cleanup_other", "uncharacterized material separated during PCR-product purification", MaterialKind.MOLECULAR)
    graph.react("purify_insert", "Purify insert PCR products", StepKind.CLEANUP,
                ((amplicon, InputRole.TEMPLATE), (cleanup, InputRole.CONDITION)), (clean, cleanup_other))
    graph.mark_terminal(Terminal(cleanup_other.id, Disposition.UNRESOLVED, "Source names purified PCR product; removed fraction composition and handling remain unassigned."))
    xhoi = graph.input("xhoi", "XhoI", MaterialKind.CATALYST)
    trimmed = graph.product("trimmed", "XhoI-treated tagged inserts", MaterialKind.MOLECULAR,
                            state=Construct([role("XhoI-treated BLESS-tagged insert; endpoint assembly unassigned")]))
    graph.react("trim_tags", "XhoI-digest purified PCR products", StepKind.DIGESTION,
                ((clean, InputRole.TEMPLATE), (xhoi, InputRole.CATALYST)), (trimmed,))
    kit_repair = graph.input("kit_repair_reagent", "ThruPLEX Template Preparation reagent reference", MaterialKind.CATALYST, inferred=True)
    kit_ready = graph.product("kit_ready", "conditional repaired blunt DNA", MaterialKind.MOLECULAR, state=trimmed.state, inferred=True)
    graph.react("kit_repair", "INFERRED — dated guide reference: template repair", StepKind.END_REPAIR,
                ((trimmed, InputRole.TEMPLATE), (kit_repair, InputRole.CATALYST)), (kit_ready,))
    kit_adapter = graph.input("kit_adapter", "dated guide 5′-blocked stem-loop adapter; exact study mapping unassigned", MaterialKind.ADAPTER,
                              state=Construct([role("unassembled kit stem-loop adapter", inferred=True)]), inferred=True)
    kit_ligase = graph.input("kit_ligase", "ThruPLEX Library Synthesis enzyme reference", MaterialKind.CATALYST, inferred=True)
    kit_ligated = graph.product("kit_ligated", "conditional nicked adapter-ligated DNA", MaterialKind.MOLECULAR,
                               state=Construct([role("unassembled nicked ThruPLEX intermediate", inferred=True)]), inferred=True)
    graph.react("kit_ligate", "INFERRED — dated guide reference: stem-loop ligation", StepKind.LIGATION,
                ((kit_ready, InputRole.TEMPLATE), (kit_adapter, InputRole.STRUCTURAL), (kit_ligase, InputRole.CATALYST)), (kit_ligated,),
                geometry=AdapterLigationGeometry.unavailable("Guide describes blunt stem-loop adapter ligation with blocked 5′ ends and a nick at genomic 3′ ends. Study-specific reagent geometry is unassembled; no generic dA/Y-adapter is substituted."))
    kit_mix = graph.input("kit_amplification", "ThruPLEX extension/cleavage/index-amplification reagents reference", MaterialKind.CATALYST, inferred=True)
    kit_fwd = graph.input("kit_forward", "conditional kit forward primer reference", MaterialKind.OLIGO, state=Construct([role("unassigned kit forward primer", inferred=True)]), inferred=True)
    kit_rev = graph.input("kit_reverse", "conditional kit reverse/index primer reference", MaterialKind.OLIGO, state=Construct([role("unassigned kit reverse primer", inferred=True)]), inferred=True)
    library = graph.product("library", "i-BLESS sequencing library", MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    graph.react("kit_amplify", "INFERRED — dated guide reference: extension, cleavage and indexed amplification", StepKind.PCR,
                ((kit_ligated, InputRole.TEMPLATE), (kit_mix, InputRole.CATALYST), (kit_fwd, InputRole.STRUCTURAL), (kit_rev, InputRole.STRUCTURAL)), (library,),
                geometry=PCRGeometry.unavailable(Construct([role("unassembled kit amplification template", inferred=True)]),
                    reason="Conditional September 2016 guide reference; study index product/member, whole primer identities and binding geometry remain unassigned."))
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY, acquired_roles=ACQUIRED_ROLES,
                    read_lengths=READ_LENGTHS, limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope((
        "sample", "deproteinization_condition", "immobilize", "proximal_repair", "repair_breaks",
        "proximal_adapter", "proximal_ligase", "ligate_proximal", "sonication", "initial_sonication",
        "agarose_recovery", "recover_dna", "fragment", "streptavidin", "capture", "distal_repair",
        "repair_distal", "distal_adapter", "distal_ligase", "ligate_distal", "iscei", "linearize",
        "insert_forward", "insert_reverse", "insert_polymerase", "insert_pcr", "insert_cleanup",
        "purify_insert", "xhoi", "trim_tags"), "i-BLESS Methods labeling paragraph and Supplementary Table 4; unresolved implementation choices explicitly marked")
    graph.declare_reviewed_scope(("kit_repair_reagent", "kit_repair", "kit_adapter", "kit_ligase",
        "kit_ligate", "kit_amplification", "kit_forward", "kit_reverse", "kit_amplify", "library"),
        "QAM-108-003 September 2016 guide Principle/Workflow; conditional reference, not assigned historical study revision")
    return graph.require_publishable()


STEM_NT = 42
PROXIMAL_TOPOLOGY = HairpinTopology(
    "i-BLESS proximal hairpin", BLESS_PROXIMAL, 0,
    len(BLESS_PROXIMAL) - STEM_NT, STEM_NT,
    colour="var(--c-r1)")
DISTAL_TOPOLOGY = HairpinTopology(
    "i-BLESS distal hairpin", BLESS_DISTAL, 0,
    len(BLESS_DISTAL) - STEM_NT, STEM_NT,
    colour="var(--c-r2)")
