import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, revcomp
from chromosome_capture import ThreeCLibrary
from conformation_replication_specs import Definition, chain, molecule, role
from endprep import EndPreparedDuplex, dA_tailed_scene
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind, Terminal)
SAMPLE_INDEX = Segment("adapter sample index", "I" * 6, "cbc", placeholder=True,
                       feature=feature("sample_index", "sample_index", "unknown"))
SUPPLEMENT = "https://doi.org/10.1038/nmeth.3922 (Supplementary Table 4)"
SUPP_FIG1 = "https://doi.org/10.1038/nmeth.3922 (Supplementary Figure 1c legend)"
SUPP_FIG1A = "https://doi.org/10.1038/nmeth.3922 (Supplementary Figure 1a legend: 3C DNA sonication)"
# GEO GSE77553/GSM2054403 belongs to the C-walks study (Olivares-Chauvet et al.,
# Nature 2016, PMID 27919068), not to the defining UMI-4C paper. Its "3way 4C"
# branch reuses UMI-4C; its reaction details are related-study context only.
RELATED = ("Related-study context only: GEO GSE77553/GSM2054403 ‘3way 4C’ extraction "
           "protocol (C-walks study, PMID 27919068); applicability to the defining "
           "UMI-4C experiments is not verified")

# Supplementary Table 4, transcribed verbatim (bases only; modifications kept below).
DS_PREFIX = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
UNIVERSAL_PRIMER_2 = "CAAGCAGAAGACGGCATACGA"
FORK_STRAND_1 = DS_PREFIX                      # printed identically, no modification
FORK_STRAND_2_5P = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"   # after (PO4)-
FORK_STRAND_2_3P = "ATCTCGTATGCCGTCTTCTGCTTG"            # after (NNNNNN)
import illumina as il
# The printed shell is the canonical TruSeq P5/Read 1, Index 1/Read 2 and P7 shell;
# P5 and the Read 1 primer site overlap by four bases (ACAC).
if not (DS_PREFIX == il.P5[:-4] + il.TRUSEQ_READ1 and FORK_STRAND_2_5P == il.INDEX1_PRIMER
        and FORK_STRAND_2_3P == il.P7_RC and il.P7.startswith(UNIVERSAL_PRIMER_2)):
    raise ValueError("UMI-4C Supplementary Table 4 transcription no longer matches the TruSeq shell")


def adapter_scene() -> Scene:
    """Annealed fork derived by computation from the two printed strands."""
    paired = FORK_STRAND_1[-13:-1]
    if revcomp(paired) != FORK_STRAND_2_5P[:12]:
        raise ValueError("printed fork strands no longer pair")
    sc = Scene()
    sc.strand("strand 1", [Segment("unpaired P5 / Read 1 arm", FORK_STRAND_1[:-13]),
                           Segment("paired fork core", paired),
                           Segment("3-prime T overhang", FORK_STRAND_1[-1])],
              label="forked adapter strand 1")
    sc.anneal("strand 2", [Segment("paired core complement", FORK_STRAND_2_5P[:12]),
                           Segment("unpaired Read 2 arm", FORK_STRAND_2_5P[12:]),
                           SAMPLE_INDEX,
                           Segment("unpaired P7 arm", FORK_STRAND_2_3P)],
              to="strand 1", pair=("paired core complement", "paired fork core"),
              unpaired=("unpaired Read 2 arm", SAMPLE_INDEX.name, "unpaired P7 arm"),
              mod5="phosphate", label="forked adapter strand 2")
    sc.ligation_end("strand 1", "3-prime T overhang", kind="cohesive", polarity="3-prime",
                    owner="adapter",
                    text="derived from the printed strands: single 3-prime T meets a dA-tailed end")
    return sc

ENDPOINT = Segment("random sonication endpoint identifier",
                   "[random shearing endpoint used as UMI]", "umi", placeholder=True,
                   feature=feature("sonication_endpoint", "umi", "random",
                                   note="physical shearing coordinate, not a synthesized barcode"))

DEFINITION = Definition(
        "UMI-4C — quantitative targeted chromosome contacts", "10.1038/nmeth.3922",
        "Schwartzman et al., Nature Methods (2016)",
        "Random shearing supplies molecule-specific endpoints before single-end adapter ligation and nested viewpoint PCR, enabling quantitative targeted 3C profiling.",
        (role("nested viewpoint / Illumina end"), role("viewpoint"),
         role("contact reporter"), role("universal adapter end")),
        "Paired-end reads (R1/R2 FASTQs in the authors' package): Read 1 begins with the DS bait primer and pad (Supplementary Table 1), Read 2 begins at the dA junction and random sonication endpoint. Read lengths and whether the six-base adapter index was read are not given in the inspected sources.", False,
        "The published ‘UMI’ is the genomic sonication endpoint, not a fixed-length random oligonucleotide. Supplementary Table 4 prints both forked-adapter strands, a six-base adapter index window, Universal Primer 2 and the US/DS bait primers; fork pairing is computed from those strands. The defining Online Methods were not available: dephosphorylation and heat release of the non-ligated adapter strand come from a related study's 3way 4C description (GEO GSE77553, C-walks paper), end repair and dA tailing are inferred (not stated in the inspected defining sources), and enzymes and cycle numbers remain unverified for the defining experiments. Supplementary Figure 1c reports linear Universal-Primer-2 products that are not sequenceable.")


# Substrate strands before adapter ligation. Both ends of a sheared fragment are random
# sonication ends; the reporter-side end is the one read after the adapter junction.
# The defining supplement does not state end repair or dA tailing: the single 3-prime dA
# is inferred from the unpaired 3-prime T computed on printed adapter strand 1.
INFERRED_TAILING = ("Not stated in the inspected defining supplement (Online Methods unread); "
                    "inferred from the single unpaired 3-prime T computed on printed adapter "
                    "strand 1 (Supplementary Table 4). The related-study 3way 4C description "
                    "also end-repairs and dA-tails.")


def contact_parts() -> list[Segment]:
    """One sonicated viewpoint-contact fragment, viewpoint side left, endpoint right."""
    return [Segment("viewpoint-side sonication end", "[viewpoint-side sheared end]",
                    "gene", placeholder=True),
            Segment("viewpoint fragment", "[viewpoint DNA]", "gene", placeholder=True),
            Segment("DpnII site", "GATC", "gene"),
            Segment("contact reporter", "[DpnII-ligated reporter DNA]", "insert",
                    placeholder=True),
            ENDPOINT]


def tailed_duplex(five_prime_phosphate: bool) -> EndPreparedDuplex:
    return EndPreparedDuplex(Construct(contact_parts(), name="end-repaired UMI-4C fragment"),
                             five_prime_phosphate=five_prime_phosphate,
                             three_prime_overhang="A", overhang_inferred=True)


def tailed_scene() -> Scene:
    """dA-tailed duplex; the 5-prime phosphate state after end repair is not reported."""
    return dA_tailed_scene(tailed_duplex(False),
                           "dA-tailed UMI-4C fragment (inferred 3-prime dA); "
                           "5-prime phosphate state unassigned")


def dephosphorylated_scene() -> Scene:
    return dA_tailed_scene(tailed_duplex(False),
                           "5-prime-dephosphorylated dA-tailed fragment "
                           "(related-study description; inferred)")


def protocol_graph() -> ProtocolGraph:
    """UMI-4C graph: supplement-printed oligos, related-study steps labelled as such."""
    definition, final, primers = DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    fixed = graph.input("fixed_chromatin", "formaldehyde-fixed chromatin",
                        MaterialKind.BIOLOGICAL,
                        state=molecule("formaldehyde-fixed chromatin",
                                       ("viewpoint", "proximal chromatin")), primary=True)
    restriction = graph.input("primary_restriction_enzyme", "DpnII",
                              MaterialKind.CATALYST)
    digested = graph.product(
        "restriction_fragments", "restriction-digested fixed chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("restriction fragments", ("viewpoint fragment",
                                                   "proximal reporter fragments")))
    graph.react("primary_digest", "Digest fixed chromatin with DpnII",
                StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE),
                 (restriction, InputRole.CATALYST)), (digested,))
    ligase = graph.input("proximity_ligase", "proximity-ligation DNA ligase",
                         MaterialKind.CATALYST)
    contact_state = ThreeCLibrary("DpnII", ("viewpoint fragment", "reporter fragment")).scene()
    contact = graph.product("contact_product", "proximity-ligated 3C contact",
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
    sonication = graph.input("sonication", "DNA sonication condition",
                             MaterialKind.CONDITION, source=SUPP_FIG1A)
    sonicated = graph.product("sonicated_3c_dna", "sonicated 3C DNA",
                              MaterialKind.MOLECULAR, state=deproteinized.state,
                              source=SUPP_FIG1A)
    graph.react("sonicate", "Sonicate recovered 3C DNA", StepKind.FRAGMENTATION,
                ((deproteinized, InputRole.TEMPLATE),
                 (sonication, InputRole.CONDITION)), (sonicated,), source=SUPP_FIG1A)
    repair = graph.input("end_repair_enzymes", "end-repair enzyme mix (unassigned)",
                         MaterialKind.CATALYST, inferred=True)
    repaired = graph.product("end_repaired_dna", "end-repaired sonicated 3C DNA",
                             MaterialKind.MOLECULAR, state=sonicated.state, inferred=True)
    graph.react("end_repair", "INFERRED — end-repair sonicated 3C DNA (not stated in inspected defining sources)", StepKind.END_REPAIR,
                ((sonicated, InputRole.TEMPLATE),
                 (repair, InputRole.CATALYST)), (repaired,), reason=INFERRED_TAILING)
    tailing = graph.input("a_tailing_enzyme", "dA-tailing enzyme (unassigned)",
                          MaterialKind.CATALYST, inferred=True)
    tailed = graph.product("da_tailed_dna", "dA-tailed 3C DNA",
                           MaterialKind.MOLECULAR, state=tailed_scene(), inferred=True)
    graph.react("da_tail", "INFERRED — add one 3-prime dA to each end (from the printed adapter 3-prime T)", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE),
                 (tailing, InputRole.CATALYST)), (tailed,), reason=INFERRED_TAILING)
    phosphatase = graph.input("alkaline_phosphatase", "alkaline phosphatase",
                              MaterialKind.CATALYST, source=RELATED, inferred=True)
    dephos = graph.product("dephosphorylated_dna", "5-prime-dephosphorylated dA-tailed 3C DNA",
                           MaterialKind.MOLECULAR, state=dephosphorylated_scene(),
                           source=RELATED, inferred=True)
    graph.react("dephosphorylate", "INFERRED — dephosphorylate 5-prime ends (related-study description)",
                StepKind.DEPHOSPHORYLATION,
                ((tailed, InputRole.TEMPLATE),
                 (phosphatase, InputRole.CATALYST)), (dephos,), source=RELATED,
                reason="Not stated in the inspected defining supplement; the C-walks study's 3way 4C description dephosphorylates so only the phosphorylated adapter strand can ligate.")
    adapters = graph.input("illumina_adapters", "printed forked adapter (strand 2 5-prime phosphate, six-base index window)",
                           MaterialKind.ADAPTER, state=adapter_scene(), source=SUPPLEMENT)
    adapter_ligase = graph.input("adapter_ligase", "adapter DNA ligase",
                                 MaterialKind.CATALYST)
    adapted = graph.product("adapter_ligated_dna", "adapter-ligated 3C DNA",
                            MaterialKind.MOLECULAR, state=final, inferred=True)
    graph.react("adapter_ligation", "Ligate Illumina adapters",
                StepKind.LIGATION,
                ((dephos, InputRole.TEMPLATE),
                 (adapters, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.known(
                    "cohesive", dephosphorylated_scene(), adapter_scene(),
                    statement="Adapter end source-known: fork pairing (12-bp core) and the single unpaired 3-prime T are computed from the two Supplementary Table 4 strands, strand 2 printed 5-prime (PO4). Substrate end inferred: the complementary single 3-prime dA is not stated in the inspected defining supplement, and the dephosphorylated substrate 5-prime ends come from the related-study 3way 4C description, under which only adapter strand 2 can be sealed (to the substrate 3-prime dA)."))
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    heat = graph.input("strand_release", "95 °C denaturation, ice, 1x SPRI",
                       MaterialKind.CONDITION, source=RELATED, inferred=True)
    released = graph.product("ligated_strand_template", "single-stranded template carrying ligated adapter strand 2",
                             MaterialKind.MOLECULAR, state=adapted.state,
                             source=RELATED, inferred=True)
    free_strand = graph.product("released_adapter_strand1", "released non-ligated adapter strand 1",
                                MaterialKind.WASTE, source=RELATED, inferred=True)
    graph.react("release_unligated_strand", "INFERRED — heat-release the non-ligated adapter strand (related-study description)",
                StepKind.STRAND_SEPARATION,
                ((adapted, InputRole.TEMPLATE),
                 (heat, InputRole.CONDITION)), (released, free_strand), source=RELATED,
                reason="Not stated in the inspected defining supplement; taken from the C-walks study's 3way 4C description.")
    graph.discard(free_strand, "related-study description: 1x SPRI after denaturation selects the UMI-4C template")
    universal = graph.input("universal_primer_2", "Universal Primer 2 (reverse primer in all reactions)",
                            MaterialKind.OLIGO,
                            state=Segment("Universal Primer 2", UNIVERSAL_PRIMER_2, "p7"),
                            source=SUPPLEMENT)
    first_forward = graph.input("nested1_viewpoint_primer", "upstream (US) bait primer, locus-specific",
                                MaterialKind.OLIGO, source=SUPPLEMENT)
    nested = graph.product("nested_amplicon", "first-round UMI-4C amplicon",
                           MaterialKind.MOLECULAR, state=adapted.state, inferred=True)
    graph.react("nested_pcr1", "First bait PCR: US primer + Universal Primer 2", StepKind.PCR,
                ((released, InputRole.TEMPLATE),
                 (first_forward, InputRole.STRUCTURAL),
                 (universal, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (nested,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The US bait primer is locus-specific; its use in the first round is inferred from the nested design (Supplementary Table 1 reads begin with the DS primer)."))
    second_forward = graph.input("nested2_viewpoint_primer", "downstream (DS) bait primer with printed P5/Read 1 prefix",
                                 MaterialKind.OLIGO, source=SUPPLEMENT)
    library = graph.product("sequencing_library", "UMI-4C sequencing library",
                            MaterialKind.LIBRARY, state=final, inferred=True)
    linear = graph.product("linear_universal_products", "candidate linear Universal-Primer-2 extension products",
                           MaterialKind.MOLECULAR, source=SUPP_FIG1, inferred=True)
    graph.react("nested_pcr2", "Second bait PCR: DS primer + Universal Primer 2", StepKind.PCR,
                ((nested, InputRole.TEMPLATE),
                 (second_forward, InputRole.STRUCTURAL),
                 (universal, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library, linear),
                geometry=PCRGeometry.unavailable(
                    final, reason="The DS bait sequence is locus-specific; only its Illumina prefix is universal."))
    graph.mark_terminal(Terminal(linear.id, Disposition.UNRESOLVED,
        "Supplementary Figure 1c quantifies linear (Universal-Primer-2-only) amplification in the final library and states such products are not compatible with Illumina sequencing; no physical separation is represented."))
    graph.sequence(library, SequencingPlan(
        final, tuple(primers), acquired_roles=tuple(p.role for p in primers
                                                    if p.role != "Index 1 (i7)"),
        limitation=definition.sequencing))
    graph.declare_reviewed_scope((
        "fixed_chromatin", "primary_restriction_enzyme", "primary_digest",
        "proximity_ligase", "proximity_ligation", "crosslink_reversal",
        "reverse_crosslinks", "sonication", "sonicate", "end_repair_enzymes",
        "end_repair", "a_tailing_enzyme", "da_tail", "alkaline_phosphatase",
        "dephosphorylate", "illumina_adapters", "adapter_ligase", "adapter_ligation",
        "strand_release", "release_unligated_strand", "universal_primer_2",
        "pcr_polymerase", "nested1_viewpoint_primer", "nested_pcr1",
        "nested2_viewpoint_primer", "nested_pcr2", "sequencing_library",
    ), source)
    return graph.require_publishable()


from types import SimpleNamespace

def workflow():
    contact = ThreeCLibrary("DpnII", ("viewpoint fragment", "reporter fragment"))
    one_end = molecule("single-end-adapter 3C fragment",
                       ("viewpoint", "contact junction", "reporter", "random sonication endpoint",
                        "single-end Illumina adapter"),
                       junctions=(("viewpoint", "contact junction", "proximity ligation"),
                                  ("random sonication endpoint", "single-end Illumina adapter", "adapter ligation")))
    return chain(molecule("formaldehyde-fixed chromatin", ("viewpoint", "proximal chromatin")), (
        ("Digest with DpnII, proximity-ligate and reverse crosslinks", contact.scene(),
         "This is a linear 3C library; UMI-4C does not perform 4C recircularization."),
        ("Sonicate, end-repair and add one 3-prime dA to each end (inferred)",
         tailed_scene(),
         "INFERRED — end repair and dA tailing are not stated in the inspected defining supplement; the single 3-prime dA is inferred from the printed adapter's unpaired 3-prime T. The random shearing coordinate is already fixed before adapter ligation."),
        ("Ligate the single-end adapter", one_end,
         "The random shearing coordinate acts as the molecular identifier; it is not a synthesized random barcode."),
        ("Amplify with two nested viewpoint PCRs", one_end,
         "The nested viewpoint primer and universal adapter primer create the sequencing library."),
    ))

d = DEFINITION
# Final library assembled from the printed oligos: the DS primer's P5/Read 1 prefix on
# the left; the ligated adapter strand 2 (Read 2 arm, six-base window, P7) on the right
# after the sonication endpoint. Reporter, pad and bait bodies are locus-specific roles.
FINAL_LIBRARY = Construct([
    Segment("P5 (DS prefix)", il.P5[:-4], "p5"),
    Segment("Read 1 primer site (DS prefix)", il.TRUSEQ_READ1, "r1"),
    Segment("DS bait primer", "[locus-specific DS bait sequence]", "gene", placeholder=True),
    Segment("bait pad", "[bait-specific pad]", "gene", placeholder=True),
    Segment("DpnII site", "GATC", "gene"),
    Segment("contact reporter", "[DpnII-ligated reporter DNA]", "insert", placeholder=True),
    ENDPOINT,
    Segment("dA junction", "A", inferred=True),
    Segment("Index 1 / Read 2 arm (adapter strand 2)", il.INDEX1_PRIMER, "r2"),
    SAMPLE_INDEX,
    Segment("P7 reverse complement (adapter strand 2)", il.P7_RC, "p7"),
], name=d.title + " final library")
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=tuple(p.role for p in SEQ_PRIMERS))
if _problems:
    raise ValueError("UMI-4C final library: " + "; ".join(_problems))
OLIGOS = (
    "<b>Forked adapter, strand 1 (Supplementary Table 4)</b> — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′",
    "<b>Forked adapter, strand 2 (Supplementary Table 4)</b> — 5′-(PO4)-GATCGGAAGAGCACACGTCTGAACTCCAGTCAC(NNNNNN)ATCTCGTATGCCGTCTTCTGCTTG-3′ (N6: adapter index window)",
    "<b>Universal Primer 2 (Supplementary Table 4)</b> — 5′-CAAGCAGAAGACGGCATACGA-3′ (reverse primer in all reactions)",
    "<b>G1p1 example viewpoint primers (Supplementary Table 4)</b> — upstream 5′-CATCCCCTCAACTCCACCTTAC-3′; downstream (DS) 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCTCAACTCCACCTTACCCACCAT-3′ (all DS primers carry this Illumina prefix); pad gatc",
)
SPEC = SimpleNamespace(
    TITLE=d.title, NOTES="01_umi_4c.html",
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
