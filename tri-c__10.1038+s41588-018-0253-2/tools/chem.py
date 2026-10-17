
from chemdraw import Segment, feature
from conformation_replication_specs import AdapterLigationGeometry, Construct, Definition, InputRole, MaterialKind, PCRGeometry, ProtocolGraph, SequencingPlan, StepKind, da_tailed, molecule, multiway_contact, role, sp
from protocol_graph import SourceClaim, SourceDiscrepancy
ENDPOINT = Segment("random sonication endpoint identifier",
                   "[random shearing endpoint used as UMI]", "umi", placeholder=True,
                   feature=feature("sonication_endpoint", "umi", "random",
                                   note="physical shearing coordinate, not a synthesized barcode"))

DEFINITION = Definition(
        "Tri-C — targeted multiway chromatin contacts", "10.1038/s41588-018-0253-2",
        "Oudelaar et al., Nature Genetics (2018)",
        "Small NlaIII fragments, gentle sonication and viewpoint capture preserve multiple proximity-ligation junctions in one sequenced molecule.",
        (role("viewpoint"), role("reporter 1"), role("second proximity junction"),
         role("reporter 2")),
        "Paired-end Illumina reads (MiSeq v2 250 bp for the first erythroid experiment; NextSeq v2 150 bp for the second multiplexed experiment) identify the viewpoint and two or more reporters within one captured fragment from ~450-bp sonication. Index-read cycles are not stated.", True,
        "<b>Published source discrepancies:</b> the 2018 article Methods describes the biotinylated capture-probe pool at 2.9 nM total, whereas the directly cited <a href='https://doi.org/10.1038/protex.2018.113'>2018 Tri-C detailed protocol</a> (Capture enrichment, p. 7) specifies 2.9 nM <i>per oligonucleotide</i>. Its materials list (p. 2) names NEBNext E6040, but the library-preparation heading (p. 5) says E7370. The <a href='https://doi.org/10.21203/rs.2.1650/v2'>2019 revision</a> (pp. 3, 9) consistently names E6040 and presents E7645 Ultra II as a later low-input option. Both detailed versions say that 10 nmol in 3.46 µl gives 2.9 µM; arithmetic gives about 2.9 mM. The <a href='https://doi.org/10.1101/307405'>preprint</a> Methods instead specifies a 13-fmol pool without a matching volume. The schematic uses the E6040 reaction family and no probe concentration; the actual 2018 mixture and kit lot remain uncertain.")


ARTICLE = ("Oudelaar et al., Nature Genetics 50:1744 (2018), doi:10.1038/s41588-018-0253-2, "
           "Online Methods 'Tri-C - Experimental procedure'")
PROTEX = ("Oudelaar et al. 2018 Tri-C Protocol Exchange v1, doi:10.1038/protex.2018.113, "
          "'Sonication and addition of Illumina sequencing adaptors' and 'Capture enrichment' "
          "pp. 5-8, Tables 3 and 5 pp. 10-11 (cited as defining-article reference 49)")


def _tri_c_graph(definition, final: Construct,
                 primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """Tri-C graph from the defining Methods and the directly cited detailed method.

    The defining article says 3C libraries were "prepared as previously described"
    (its ref. 41, NG Capture-C) and capture used NimbleGen SeqCap EZ reagents "as
    previously described". Steps whose detail comes only from the cited 2018
    Protocol Exchange method carry that source explicitly; nothing is taken from
    the predecessor's own text.
    """
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    fixed = graph.input("fixed_chromatin", "cells fixed in 2% formaldehyde for 10 min, "
                        "lysed, Dounce-homogenized, SDS/Triton X-100 treated",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("formaldehyde-fixed chromatin",
                                       ("viewpoint", "proximal chromatin")),
                        primary=True, source=ARTICLE)
    restriction = graph.input("primary_restriction_enzyme",
                              "NlaIII, 3 x 600 U over 16-24 h at 37 °C", MaterialKind.CATALYST,
                              source=ARTICLE)
    digested = graph.product(
        "restriction_fragments", "NlaIII-digested fixed chromatin", MaterialKind.MOLECULAR,
        state=molecule("restriction fragments", ("viewpoint fragment",
                                                 "proximal reporter fragments")),
        source=ARTICLE)
    graph.react("primary_digest", "Digest fixed chromatin with NlaIII", StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)),
                (digested,), source=ARTICLE)
    heat = graph.input("digest_heat_inactivation", "65 °C for 20 min, then ice",
                       MaterialKind.CONDITION, source=ARTICLE)
    inactivated = graph.product("inactivated_digest", "heat-inactivated NlaIII digest",
                                MaterialKind.MOLECULAR, state=digested.state, source=ARTICLE)
    graph.react("inactivate_digest", "Heat-inactivate the digestion", StepKind.THERMAL_TREATMENT,
                ((digested, InputRole.TEMPLATE), (heat, InputRole.CONDITION)),
                (inactivated,), source=ARTICLE, changes_molecule=False)
    ligase = graph.input("proximity_ligase", "T4 DNA ligase, 720 U, overnight at 16 °C",
                         MaterialKind.CATALYST, source=ARTICLE)
    contact_model = multiway_contact(
        "NlaIII", "viewpoint fragment", ("reporter fragment 1", "reporter fragment 2"))
    contact_state = contact_model.scene("multiway 3C concatemer")
    contact = graph.product("contact_product", "multiway proximity-ligation concatemer",
                            MaterialKind.MOLECULAR, state=contact_state, source=ARTICLE)
    graph.react("proximity_ligation", "Proximity-ligate NlaIII fragments in nuclei",
                StepKind.PROXIMITY_LIGATION,
                ((inactivated, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
                (contact,), source=ARTICLE)
    proteinase = graph.input("crosslink_reversal",
                             "de-crosslinking with Proteinase K (3 U) overnight at 65 °C",
                             MaterialKind.CONDITION, source=ARTICLE)
    deproteinized = graph.product("deproteinized_3c_dna", "de-crosslinked 3C DNA",
                                  MaterialKind.MOLECULAR, state=contact_state, source=ARTICLE)
    graph.react("reverse_crosslinks", "Reverse cross-links and digest protein",
                StepKind.DEPROTEINIZATION,
                ((contact, InputRole.TEMPLATE), (proteinase, InputRole.CONDITION)),
                (deproteinized,), source=ARTICLE)
    extraction = graph.input("dna_extraction",
                             "RNase (7.5 mU), phenol-chloroform, sodium acetate/ethanol "
                             "precipitation at -80 °C, 70% ethanol wash",
                             MaterialKind.CONDITION, source=ARTICLE)
    library_3c = graph.product("purified_3c_library", "purified 3C library DNA",
                               MaterialKind.MOLECULAR, state=contact_state, source=ARTICLE)
    graph.react("extract_3c_dna", "RNase-treat and extract 3C DNA", StepKind.EXTRACTION,
                ((deproteinized, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)),
                (library_3c,), source=ARTICLE, changes_molecule=False)

    sonication = graph.input("sonication",
                             "Covaris S220, one 55-s cycle, 10% duty, intensity 4, "
                             "200 cycles/burst, target ~450 bp", MaterialKind.CONDITION,
                             source=ARTICLE)
    sonicated = graph.product("sonicated_3c_dna", "~450-bp sonicated 3C DNA",
                              MaterialKind.MOLECULAR, state=contact_state, source=ARTICLE)
    graph.react("sonicate", "Sonicate gently to retain multiple junctions",
                StepKind.FRAGMENTATION,
                ((library_3c, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (sonicated,), source=ARTICLE)
    beads_sonication = graph.input("sonication_ampure", "AMPure XP, 0.7:1 bead-to-sample",
                                   MaterialKind.SELECTION, source=ARTICLE)
    sized = graph.product("size_selected_3c_dna", "bead-recovered sonicated 3C DNA",
                          MaterialKind.MOLECULAR, state=contact_state, source=ARTICLE)
    sized_waste = graph.product("sonication_cleanup_supernatant",
                                "bead supernatant (fragments and solutes not bound at this ratio)",
                                MaterialKind.WASTE, source=ARTICLE)
    graph.react("sonication_cleanup", "Clean up and size-select after sonication",
                StepKind.CLEANUP,
                ((sonicated, InputRole.TEMPLATE), (beads_sonication, InputRole.SELECTION)),
                (sized, sized_waste), source=ARTICLE, changes_molecule=False)
    graph.discard(sized_waste, "supernatant removed from the magnetically held AMPure beads")
    graph.annotate(SourceDiscrepancy(
        "sonication_cleanup", "Post-sonication bead ratio",
        (SourceClaim("Defining article Online Methods, Tri-C experimental procedure",
                     "Clean-up and size selection after sonication used AMPure XP at 0.7:1."),
         SourceClaim("2018 Protocol Exchange method, Sonication step I.3 and note E (p. 5)",
                     "220 µl beads are added to the 120-µl sonicate (about 1.8x); "
                     "all clean-ups use 1.8 volumes and size selection is called unnecessary.")),
        "The graph uses the article's 0.7:1 size-selecting clean-up for the published experiment.",
        "The defining article describes the libraries actually sequenced; the detailed method "
        "is a general protocol.",
        "Which ratio each sequenced library received is not stated beyond the article's sentence."))

    repair = graph.input("end_repair_enzymes",
                         "NEBNext End Repair buffer and enzyme mix, 20 °C 30 min",
                         MaterialKind.CATALYST, source=PROTEX)
    repaired = graph.product("end_repaired_dna", "end-repaired sonicated 3C DNA",
                             MaterialKind.MOLECULAR, state=contact_state, source=PROTEX)
    graph.react("end_repair", "End-repair sonicated 3C DNA", StepKind.END_REPAIR,
                ((sized, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,), source=PROTEX)
    tailing = graph.input("a_tailing_enzyme", "NEBNext dA-tailing buffer and Klenow, 37 °C 30 min",
                          MaterialKind.CATALYST, source=PROTEX)
    tailed = graph.product("da_tailed_dna", "dA-tailed 3C DNA", MaterialKind.MOLECULAR,
                           state=da_tailed("dA-tailed 3C fragment", "viewpoint-contact fragment"),
                           source=PROTEX)
    graph.react("da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)),
                (tailed,), source=PROTEX)
    adapters = graph.input("illumina_adapters",
                           "Illumina TruSeq-compatible adaptor supplied with the NEBNext kit",
                           MaterialKind.ADAPTER, inferred=True, source=PROTEX)
    adapter_ligase = graph.input("adapter_ligase",
                                 "Quick Ligation buffer and T4 DNA ligase, 20 °C 15 min",
                                 MaterialKind.CATALYST, source=PROTEX)
    ligated = graph.product("adapter_ligated_dna", "adapter-ligated 3C DNA",
                            MaterialKind.MOLECULAR, state=final, inferred=True, source=PROTEX)
    graph.react("adapter_ligation", "Ligate Illumina adaptors", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The article names Illumina TruSeq adaptors added with NEBNext reagents and "
                    "the detailed method names the kit 'Adaptor'; neither prints the adaptor "
                    "strands, their uracil position or the ligation-ready end."),
                source=PROTEX)
    user = graph.input("user_enzyme", "USER Enzyme Mix (3 µl), 37 °C 15 min",
                       MaterialKind.CATALYST, source=PROTEX)
    user_treated = graph.product("user_treated_dna", "USER-treated adapter-ligated 3C DNA",
                                 MaterialKind.MOLECULAR, state=final, inferred=True,
                                 source=PROTEX)
    graph.react("user_treatment", "Treat the ligation with USER enzyme", StepKind.DIGESTION,
                ((ligated, InputRole.TEMPLATE), (user, InputRole.CATALYST)),
                (user_treated,), source=PROTEX,
                reason="The detailed method adds USER directly to the ligation. Its molecular "
                       "target (a uracil in the kit adaptor) is not printed in the inspected sources.")
    beads_ligation = graph.input("ligation_ampure", "AMPure XP, 1.8:1 bead-to-sample",
                                 MaterialKind.SELECTION, source=ARTICLE)
    cleaned = graph.product("cleaned_adapter_ligated_dna", "bead-recovered adapter-ligated DNA",
                            MaterialKind.MOLECULAR, state=final, inferred=True, source=ARTICLE)
    cleaned_waste = graph.product("ligation_cleanup_supernatant",
                                  "bead supernatant after ligation/USER",
                                  MaterialKind.WASTE, source=ARTICLE)
    graph.react("ligation_cleanup", "Bead clean-up after ligation and USER treatment",
                StepKind.CLEANUP,
                ((user_treated, InputRole.TEMPLATE), (beads_ligation, InputRole.SELECTION)),
                (cleaned, cleaned_waste), source=ARTICLE, changes_molecule=False)
    graph.discard(cleaned_waste, "supernatant removed from the magnetically held AMPure beads")

    polymerase = graph.input("pcr_polymerase", "Agilent Herculase II Fusion DNA polymerase",
                             MaterialKind.CATALYST, source=ARTICLE)
    forward = graph.input("index_forward_primer", "NEB Universal primer", MaterialKind.OLIGO,
                          source=PROTEX)
    reverse = graph.input("index_reverse_primer", "NEB Index primer (E7335/E7500 sets)",
                          MaterialKind.OLIGO, source=PROTEX)
    indexed = graph.product("indexed_library", "indexed adapter-ligated 3C library",
                            MaterialKind.MOLECULAR, state=final, source=ARTICLE)
    graph.react("index_pcr", "Seven-cycle indexing PCR", StepKind.PCR,
                ((cleaned, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (indexed,),
                geometry=PCRGeometry.unavailable(
                    final, reason="Article: 7-cycle Herculase II PCR; Protocol Exchange Table 3 "
                                  "names NEB Universal and Index primers without bases or the "
                                  "index chosen for each library."),
                source=ARTICLE)
    others = graph.input("other_indexed_libraries",
                         "other differentially indexed Tri-C libraries",
                         MaterialKind.MOLECULAR, source=ARTICLE)
    pooled = graph.product("multiplexed_library", "multiplexed indexed libraries, "
                           "1.5 µg per library", MaterialKind.MOLECULAR, state=final,
                           source=ARTICLE)
    graph.react("pool_libraries", "Pool libraries for one multiplexed capture", StepKind.POOL,
                ((indexed, InputRole.TEMPLATE), (others, InputRole.TEMPLATE)),
                (pooled,), source=ARTICLE, changes_molecule=False)

    blockers = graph.input("capture_blockers",
                           "COT DNA plus NimbleGen TS-HE Universal and TS-HE Index blocking oligos",
                           MaterialKind.CONDITION, source=PROTEX)
    beads = graph.input("streptavidin_beads", "M-270 streptavidin beads and SeqCap washes",
                        MaterialKind.SELECTION, source=PROTEX)
    capture_pcr = graph.input("post_capture_pcr_mix",
                              "KAPA HiFi HotStart ReadyMix", MaterialKind.CATALYST, source=PROTEX)
    post_lm_1 = graph.input("post_lm_pcr_oligo_1", "NimbleGen Post-LM-PCR Oligo 1",
                            MaterialKind.OLIGO, source=PROTEX)
    post_lm_2 = graph.input("post_lm_pcr_oligo_2", "NimbleGen Post-LM-PCR Oligo 2",
                            MaterialKind.OLIGO, source=PROTEX)
    current = pooled
    reviewed_rounds: list[str] = []
    for n, hours in ((1, "64-72 h"), (2, "24 h")):
        probes = graph.input(f"capture_oligo_round{n}",
                             "biotinylated 120-nt viewpoint capture oligos "
                             "(Supplementary Table 5 coordinates)",
                             MaterialKind.OLIGO, source=ARTICLE)
        hybrid = graph.product(f"capture_hybrids_round{n}",
                               f"round-{n} library/probe hybridization mix",
                               MaterialKind.MOLECULAR, state=final, source=PROTEX)
        graph.react(f"hybridize_round{n}",
                    f"Denature at 95 °C and hybridize with probes at 47 °C for {hours}",
                    StepKind.HYBRIDIZATION,
                    ((current, InputRole.TEMPLATE), (probes, InputRole.STRUCTURAL),
                     (blockers, InputRole.CONDITION)), (hybrid,), source=PROTEX,
                    changes_molecule=False)
        bound = graph.product(f"captured_round{n}", f"round-{n} bead-bound captured library",
                              MaterialKind.MOLECULAR, state=final, source=PROTEX)
        reject = graph.product(f"capture_round{n}_reject", f"round-{n} unbound/washed-off library",
                               MaterialKind.WASTE, source=PROTEX)
        graph.react(f"capture_round{n}_step", f"Round-{n} streptavidin capture and stringent washes",
                    StepKind.AFFINITY_CAPTURE,
                    ((hybrid, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                    (bound, reject), source=PROTEX, changes_molecule=False)
        graph.discard(reject, "removed in bead supernatants and 47 °C stringent washes")
        last = n == 2
        amplified = graph.product(
            "sequencing_library" if last else f"amplified_capture_round{n}",
            definition.title + " captured library" if last else
            f"round-{n} amplified captured library",
            MaterialKind.LIBRARY if last else MaterialKind.MOLECULAR, state=final,
            source=PROTEX)
        graph.react(f"post_capture_pcr_round{n}",
                    f"Amplify round-{n} capture off the beads (KAPA HiFi, Table 5)",
                    StepKind.PCR,
                    ((bound, InputRole.TEMPLATE), (post_lm_1, InputRole.STRUCTURAL),
                     (post_lm_2, InputRole.STRUCTURAL), (capture_pcr, InputRole.CATALYST)),
                    (amplified,),
                    geometry=PCRGeometry.unavailable(
                        final, reason="Protocol Exchange Table 5 names Post-LM-PCR Oligos 1 and 2 "
                                      "(no bases) and a 98/60/72 °C programme with 'go to step 2, "
                                      "14 cycles'; DNA is amplified on the beads, not eluted."),
                    source=PROTEX)
        reviewed_rounds += [f"capture_oligo_round{n}", f"hybridize_round{n}",
                            f"capture_round{n}_step", f"post_capture_pcr_round{n}"]
        current = amplified
    graph.sequence(current, SequencingPlan(final, tuple(primers), tuple(p.role for p in primers)))
    graph.declare_reviewed_scope([
        "fixed_chromatin", "primary_restriction_enzyme", "primary_digest",
        "digest_heat_inactivation", "inactivate_digest", "proximity_ligase",
        "proximity_ligation", "crosslink_reversal", "reverse_crosslinks", "dna_extraction",
        "extract_3c_dna", "sonication", "sonicate", "sonication_ampure", "sonication_cleanup",
        "end_repair_enzymes", "end_repair", "a_tailing_enzyme", "da_tail",
        "illumina_adapters", "adapter_ligase", "adapter_ligation", "user_enzyme",
        "user_treatment", "ligation_ampure", "ligation_cleanup", "pcr_polymerase",
        "index_forward_primer", "index_reverse_primer", "index_pcr",
        "other_indexed_libraries", "pool_libraries", "capture_blockers",
        "streptavidin_beads", "post_capture_pcr_mix", "post_lm_pcr_oligo_1",
        "post_lm_pcr_oligo_2", *reviewed_rounds, "sequencing_library",
    ], f"{ARTICLE}; {PROTEX}")
    return graph.require_publishable()


def protocol_graph():
    return _tri_c_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

from types import SimpleNamespace
from conformation_replication_specs import chain

def workflow():
    contact = multiway_contact("NlaIII", "viewpoint fragment",
                               ("reporter fragment 1", "reporter fragment 2"))
    sonicated = molecule("sonicated multiway 3C fragment",
                         ("random sonication endpoint", "viewpoint", "reporter 1",
                          "reporter 2", "random sonication endpoint"),
                         junctions=(("viewpoint", "reporter 1", "proximity ligation"),
                                    ("reporter 1", "reporter 2", "proximity ligation")))
    adapted = molecule("adapter-ligated multiway 3C fragment",
                       ("left TruSeq adapter", "random sonication endpoint", "viewpoint",
                        "reporter 1", "reporter 2", "right TruSeq adapter"),
                       junctions=(("left TruSeq adapter", "random sonication endpoint", "adapter ligation"),
                                  ("viewpoint", "reporter 1", "proximity ligation"),
                                  ("reporter 1", "reporter 2", "proximity ligation"),
                                  ("reporter 2", "right TruSeq adapter", "adapter ligation")))
    return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "two proximal regions")), (
        ("Digest with NlaIII and proximity-ligate in nuclei", contact.scene("multiway 3C concatemer"),
         "At least two ligation junctions in one molecule encode a concurrent contact."),
        ("Reverse crosslinks and sonicate gently to about 450 bp", sonicated,
         "Longer fragments preserve multiple short NlaIII products and their junctions."),
        ("Ligate adapters, treat with USER, index, then double-capture the viewpoint", adapted,
         "The detailed Tri-C method places USER treatment after adapter ligation. A biotinylated oligo targets the middle of the short viewpoint restriction fragment; each capture round is followed by PCR amplification off the streptavidin beads."),
    ))

d = DEFINITION
if d.standard_truseq:
    from batch_ngs import truseq_library
    FINAL_LIBRARY, SEQ_PRIMERS = truseq_library(list(d.final_parts), d.title + " final library", dual_index=False)
else:
    FINAL_LIBRARY = Construct([role("Illumina-compatible left end"), *d.final_parts, role("Illumina-compatible right end")], name=d.title + " final library")
    SEQ_PRIMERS = ()
OLIGOS = ()
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES="01_tri_c.html",
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
