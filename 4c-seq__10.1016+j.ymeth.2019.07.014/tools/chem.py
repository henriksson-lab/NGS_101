
from conformation_replication_specs import _three_c
from conformation_replication_specs import AdapterLigationGeometry, Circular4CTemplate, CircularizationGeometry, Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, da_tailed, inverse_pcr_product, molecule, multiway_contact, role, sp

DEFINITION = Definition(
        "4C-seq — circular chromosome conformation capture sequencing",
        "10.1016/j.ymeth.2019.07.014", "Krijger et al., Methods (2020)",
        "Two restriction digests and two ligations convert a viewpoint-centred 3C contact into a circle that can be opened by outward-facing inverse PCR.",
        (role("viewpoint PCR / Illumina end"), role("viewpoint fragment"),
         role("captured reporter"), role("second viewpoint PCR / Illumina end")),
        "A viewpoint-specific reading primer enters the captured reporter beside the primary restriction junction; sample and flow-cell sequences are added by the two PCR stages.", False,
        "The genomic portions of viewpoint primers are locus-specific. Published universal tails and second-round primer sequences are listed below; the chosen genomic portions and index depend on the experiment.")


def _conformation_graph(key: str, definition, final: Construct,
                        primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """Typed graphs for 4C/Capture-C derivatives with distinct contact ligation."""
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    fixed = graph.input("fixed_chromatin", "formaldehyde-fixed chromatin",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("formaldehyde-fixed chromatin",
                                       ("viewpoint", "proximal chromatin")), primary=True)
    restriction_name = {"4c-seq": "DpnII", "ng-capture-c": "DpnII or NlaIII",
                        "umi-4c": "DpnII", "tri-c": "NlaIII"}[key]
    restriction = graph.input("primary_restriction_enzyme", restriction_name,
                              MaterialKind.CATALYST)
    digested = graph.product(
        "restriction_fragments", "restriction-digested fixed chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("restriction fragments", ("viewpoint fragment",
                                                   "proximal reporter fragments")))
    graph.react("primary_digest", f"Digest fixed chromatin with {restriction_name}",
                StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE),
                 (restriction, InputRole.CATALYST)), (digested,))
    ligase = graph.input("proximity_ligase", "proximity-ligation DNA ligase",
                         MaterialKind.CATALYST)
    if key == "tri-c":
        contact_model = multiway_contact(
            "NlaIII", "viewpoint fragment", ("reporter fragment 1", "reporter fragment 2"))
        contact_state = contact_model.scene("multiway 3C concatemer")
        contact_label = "multiway proximity-ligation concatemer"
    else:
        names = (("captured viewpoint fragment", "reporter fragment")
                 if key == "ng-capture-c" else
                 ("viewpoint fragment", "reporter fragment"))
        contact_model = _three_c(restriction_name, names)
        contact_state = contact_model.scene()
        contact_label = "proximity-ligated 3C contact"
    contact = graph.product("contact_product", contact_label,
                            MaterialKind.MOLECULAR, state=contact_state)
    graph.react("proximity_ligation", "Proximity-ligate restriction fragments",
                StepKind.PROXIMITY_LIGATION,
                ((digested, InputRole.TEMPLATE),
                 (ligase, InputRole.CATALYST)), (contact,))
    proteinase = graph.input("crosslink_reversal", "proteinase/cross-link reversal condition",
                             MaterialKind.CONDITION)
    deproteinized = graph.product("deproteinized_3c_dna", "deproteinized 3C DNA",
                                  MaterialKind.MOLECULAR, state=contact_state)
    graph.react("reverse_crosslinks", "Reverse cross-links and recover 3C DNA",
                StepKind.DEPROTEINIZATION,
                ((contact, InputRole.TEMPLATE),
                 (proteinase, InputRole.CONDITION)), (deproteinized,))
    reviewed = [
        "fixed_chromatin", "primary_restriction_enzyme", "primary_digest",
        "proximity_ligase", "proximity_ligation", "crosslink_reversal",
        "reverse_crosslinks",
    ]

    if key == "4c-seq":
        second = graph.input("secondary_restriction_enzyme", "Csp6I",
                             MaterialKind.CATALYST)
        secondary = graph.product("secondary_fragments", "Csp6I-digested 3C products",
                                  MaterialKind.MOLECULAR, state=contact_state)
        graph.react("secondary_digest", "Digest 3C products with Csp6I",
                    StepKind.DIGESTION,
                    ((deproteinized, InputRole.TEMPLATE),
                     (second, InputRole.CATALYST)), (secondary,))
        circle_model = Circular4CTemplate(contact_model, "Csp6I", "viewpoint fragment")
        circular = graph.product("circular_template", "recircularized 4C template",
                                 MaterialKind.MOLECULAR, state=circle_model.scene())
        graph.react("recircularize", "Ligate dilute DNA to form 4C circles",
                    StepKind.CIRCULARIZATION,
                    ((secondary, InputRole.TEMPLATE),
                     (ligase, InputRole.CATALYST)), (circular,),
                    geometry=CircularizationGeometry.known(
                        "double-strand intramolecular ligation",
                        statement="The 4C workflow dilute-ligates the secondary restriction fragments into circles."))
        inverse_forward = graph.input("inverse_forward_primer", "outward viewpoint primer 1",
                                      MaterialKind.OLIGO)
        inverse_reverse = graph.input("inverse_reverse_primer", "outward viewpoint primer 2",
                                      MaterialKind.OLIGO)
        polymerase = graph.input("pcr_polymerase", "PCR polymerase",
                                 MaterialKind.CATALYST)
        first_amplicon = graph.product("inverse_pcr_product", "viewpoint inverse-PCR product",
                                      MaterialKind.MOLECULAR,
                                      state=inverse_pcr_product(circle_model))
        graph.react("inverse_pcr", "Inverse-PCR outward from the viewpoint",
                    StepKind.PCR,
                    ((circular, InputRole.TEMPLATE),
                     (inverse_forward, InputRole.STRUCTURAL),
                     (inverse_reverse, InputRole.STRUCTURAL),
                     (polymerase, InputRole.CATALYST)), (first_amplicon,),
                    geometry=PCRGeometry.unavailable(
                        final, reason="Viewpoint primer sequences are locus-specific."))
        handle_forward = graph.input("handle_forward_primer", "Illumina-handle PCR forward primer",
                                     MaterialKind.OLIGO, inferred=True)
        handle_reverse = graph.input("handle_reverse_primer", "indexed Illumina-handle PCR reverse primer",
                                     MaterialKind.OLIGO, inferred=True)
        library = graph.product("sequencing_library", "4C-seq sequencing library",
                                MaterialKind.LIBRARY, state=final, inferred=True)
        graph.react("handle_pcr", "Add sample and Illumina sequences by PCR",
                    StepKind.PCR,
                    ((first_amplicon, InputRole.TEMPLATE),
                     (handle_forward, InputRole.STRUCTURAL),
                     (handle_reverse, InputRole.STRUCTURAL),
                     (polymerase, InputRole.CATALYST)), (library,),
                    geometry=PCRGeometry.unavailable(
                        final, reason="The exact viewpoint/handle combination varies by experiment."))
        graph.sequence(library, SequencingPlan(limitation=definition.sequencing))
        reviewed.extend((
            "secondary_restriction_enzyme", "secondary_digest", "recircularize",
            "inverse_forward_primer", "inverse_reverse_primer", "pcr_polymerase",
            "inverse_pcr", "handle_forward_primer", "handle_reverse_primer",
            "handle_pcr", "sequencing_library",
        ))
    else:
        sonication = graph.input("sonication", "DNA sonication condition",
                                 MaterialKind.CONDITION)
        sonicated = graph.product("sonicated_3c_dna", "sonicated 3C DNA",
                                  MaterialKind.MOLECULAR, state=deproteinized.state)
        graph.react("sonicate", "Sonicate recovered 3C DNA", StepKind.FRAGMENTATION,
                    ((deproteinized, InputRole.TEMPLATE),
                     (sonication, InputRole.CONDITION)), (sonicated,))
        repair = graph.input("end_repair_enzymes", "end-repair enzyme mix",
                             MaterialKind.CATALYST)
        repaired = graph.product("end_repaired_dna", "end-repaired sonicated 3C DNA",
                                 MaterialKind.MOLECULAR, state=sonicated.state)
        graph.react("end_repair", "End-repair sonicated 3C DNA", StepKind.END_REPAIR,
                    ((sonicated, InputRole.TEMPLATE),
                     (repair, InputRole.CATALYST)), (repaired,))
        tailing = graph.input("a_tailing_enzyme", "dA-tailing enzyme",
                              MaterialKind.CATALYST)
        tailed_state = da_tailed("dA-tailed 3C fragment", "viewpoint-contact fragment")
        tailed = graph.product("da_tailed_dna", "dA-tailed 3C DNA",
                               MaterialKind.MOLECULAR, state=tailed_state)
        graph.react("da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                    ((repaired, InputRole.TEMPLATE),
                     (tailing, InputRole.CATALYST)), (tailed,))
        adapters = graph.input("illumina_adapters", "Illumina library adapters",
                               MaterialKind.ADAPTER, inferred=True)
        adapter_ligase = graph.input("adapter_ligase", "adapter DNA ligase",
                                     MaterialKind.CATALYST)
        adapted = graph.product("adapter_ligated_dna", "adapter-ligated 3C DNA",
                                MaterialKind.MOLECULAR, state=final, inferred=True)
        graph.react("adapter_ligation", "Ligate Illumina adapters",
                    StepKind.LIGATION,
                    ((tailed, InputRole.TEMPLATE),
                     (adapters, InputRole.STRUCTURAL),
                     (adapter_ligase, InputRole.CATALYST)), (adapted,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The defining protocol identifies its library kit or adapter role but does not print one immutable adapter-end reagent."))
        reviewed.extend((
            "sonication", "sonicate", "end_repair_enzymes", "end_repair",
            "a_tailing_enzyme", "da_tail", "illumina_adapters", "adapter_ligase",
            "adapter_ligation",
        ))
        polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                                 MaterialKind.CATALYST)

        if key == "umi-4c":
            first_forward = graph.input("nested1_viewpoint_primer", "first nested viewpoint primer",
                                        MaterialKind.OLIGO)
            first_reverse = graph.input("nested1_adapter_primer", "first universal-adapter primer",
                                        MaterialKind.OLIGO)
            nested = graph.product("nested_amplicon", "first nested UMI-4C amplicon",
                                   MaterialKind.MOLECULAR, state=adapted.state)
            graph.react("nested_pcr1", "First nested viewpoint PCR", StepKind.PCR,
                        ((adapted, InputRole.TEMPLATE),
                         (first_forward, InputRole.STRUCTURAL),
                         (first_reverse, InputRole.STRUCTURAL),
                         (polymerase, InputRole.CATALYST)), (nested,),
                        geometry=PCRGeometry.unavailable(
                            final, reason="The viewpoint primer is locus-specific."))
            second_forward = graph.input("nested2_viewpoint_primer", "indexed nested viewpoint primer",
                                         MaterialKind.OLIGO)
            second_reverse = graph.input("nested2_adapter_primer", "indexed universal-adapter primer",
                                         MaterialKind.OLIGO)
            library = graph.product("sequencing_library", "UMI-4C sequencing library",
                                    MaterialKind.LIBRARY, state=final, inferred=True)
            graph.react("nested_pcr2", "Second nested viewpoint PCR", StepKind.PCR,
                        ((nested, InputRole.TEMPLATE),
                         (second_forward, InputRole.STRUCTURAL),
                         (second_reverse, InputRole.STRUCTURAL),
                         (polymerase, InputRole.CATALYST)), (library,),
                        geometry=PCRGeometry.unavailable(
                            final, reason="The indexed viewpoint primer remains experiment-specific."))
            graph.sequence(library, SequencingPlan(limitation=definition.sequencing))
            reviewed.extend((
                "pcr_polymerase", "nested1_viewpoint_primer", "nested1_adapter_primer",
                "nested_pcr1", "nested2_viewpoint_primer", "nested2_adapter_primer",
                "nested_pcr2", "sequencing_library",
            ))
        else:
            forward = graph.input("index_forward_primer", "indexed library PCR forward primer",
                                  MaterialKind.OLIGO)
            reverse = graph.input("index_reverse_primer", "indexed library PCR reverse primer",
                                  MaterialKind.OLIGO)
            indexed = graph.product("indexed_library", "indexed 3C library",
                                    MaterialKind.MOLECULAR, state=final)
            graph.react("index_pcr", "Index and amplify the 3C library", StepKind.PCR,
                        ((adapted, InputRole.TEMPLATE),
                         (forward, InputRole.STRUCTURAL),
                         (reverse, InputRole.STRUCTURAL),
                         (polymerase, InputRole.CATALYST)), (indexed,),
                        geometry=PCRGeometry.unavailable(
                            final, reason="Index choice varies while the final TruSeq endpoint is canonical."))
            capture1 = graph.input("capture_oligo_round1", "biotinylated viewpoint capture oligo round 1",
                                   MaterialKind.SELECTION)
            selected1 = graph.product("capture_round1", "first-round viewpoint-enriched library",
                                      MaterialKind.MOLECULAR, state=final)
            reject1 = graph.product("capture_round1_reject", "uncaptured round-1 library",
                                    MaterialKind.WASTE)
            graph.react("capture_round1_step", "First viewpoint capture",
                        StepKind.AFFINITY_CAPTURE,
                        ((indexed, InputRole.TEMPLATE),
                         (capture1, InputRole.SELECTION)), (selected1, reject1))
            graph.discard(reject1, "does not hybridize to the selected viewpoint probe")
            capture2 = graph.input("capture_oligo_round2", "biotinylated viewpoint capture oligo round 2",
                                   MaterialKind.SELECTION)
            library = graph.product("sequencing_library", definition.title + " captured library",
                                    MaterialKind.LIBRARY, state=final)
            reject2 = graph.product("capture_round2_reject", "uncaptured round-2 library",
                                    MaterialKind.WASTE)
            graph.react("capture_round2_step", "Second viewpoint capture",
                        StepKind.AFFINITY_CAPTURE,
                        ((selected1, InputRole.TEMPLATE),
                         (capture2, InputRole.SELECTION)), (library, reject2))
            graph.discard(reject2, "does not survive the second viewpoint capture")
            graph.sequence(library, SequencingPlan(
                final, tuple(primers), tuple(p.role for p in primers)))
            reviewed.extend((
                "pcr_polymerase", "index_forward_primer", "index_reverse_primer",
                "index_pcr", "capture_oligo_round1", "capture_round1_step",
                "capture_oligo_round2", "capture_round2_step", "sequencing_library",
            ))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = SPEC.FINAL_LIBRARY
    primers = SPEC.SEQ_PRIMERS
    key = '4c-seq'
    return _conformation_graph(key, d, final, primers)


from types import SimpleNamespace
from conformation_replication_specs import chain, _three_c
from chromosome_capture import ThreeCLibrary
from replication_labeling import CopyNumberReplication

def workflow():
    contact = _three_c("DpnII", ("viewpoint fragment", "contacted reporter fragment"))
    circle = Circular4CTemplate(contact, "Csp6I", "viewpoint fragment")
    return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "proximal chromatin")), (
        ("Digest chromatin with the primary restriction enzyme and proximity-ligate",
         contact.scene(), "The first junction records proximity in fixed nuclei."),
        ("Reverse crosslinks, digest with a second enzyme and ligate dilute DNA",
         circle.scene(), "The second ligation recircularizes shorter 4C templates."),
        ("Inverse-PCR outward from the chosen viewpoint", inverse_pcr_product(circle),
         "A two-stage PCR adds sample and Illumina sequences to viewpoint-specific products."),
    ))

d = DEFINITION
if d.standard_truseq:
    from batch_ngs import truseq_library
    FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(list(d.final_parts), d.title + " final library", dual_index=False)
else:
    FINAL_LIBRARY = Construct([role("Illumina-compatible left end"), *d.final_parts, role("Illumina-compatible right end")], name=d.title + " final library")
    SEQ_PRIMERS = ()
OLIGOS = (
    "<b>Reading inverse-PCR primer</b> — 5′-TACACGACGCTCTTCCGATCT-[0–5-nt spacer]-[viewpoint-specific sequence]-3′",
    "<b>Non-reading inverse-PCR primer</b> — 5′-ACTGGAGTTCAGACGTGTGCTCTTCCGATCT-[viewpoint-specific sequence]-3′",
    "<b>Second-round universal forward</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′",
)
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES="01_4c_seq.html",
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=SEQ_PRIMERS, FINAL_CAPTION=d.sequencing,
    SEQUENCING_INTRO=d.sequencing, READ_LENGTHS=None,
    SEQUENCING_ENDING=d.sequencing,
    SEQUENCING_UNAVAILABLE=("The defining protocol uses locus-specific or custom primer sequence that cannot be represented by one universal binding sequence."
                            if not SEQ_PRIMERS else ""),
    workflow=workflow, oligos=lambda: OLIGOS,
)
SPEC.protocol_graph = protocol_graph
