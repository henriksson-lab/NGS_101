"""Targeted chromatin-cleavage and native occupancy library chemistries."""
from __future__ import annotations

from types import SimpleNamespace

from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import truseq_library
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, da_tailed, molecule, role


# Orsi et al., Current Protocols 2015, Basic Protocol 3.  The asterisk is an
# ordering phosphorothioate; the molecular library constructor below uses the
# corresponding unmodified bases and keeps the chemistry in the displayed oligos.
CHEC_UNIVERSAL_ADAPTER = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
CHEC_INDEXED_ADAPTER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
CHEC_P5_PCR = "AATGATACGGCGACCACCGAG"
CHEC_P7_PCR = "CAAGCAGAAGACGGCATACGAG"
HENIKOFF_DA_LIGATION = reported_adapter_end(
    "phosphorothioate-protected TruSeq adapter", kind="cohesive",
    polarity="3-prime", terminal_sequence="T",
    statement=("The cited Henikoff low-input procedure publishes the protected TruSeq "
               "adapter and its terminal dT/dA ligation compatibility."))


def henikoff_single_index_library(name: str, payload):
    """Published short-insert paired-end shell shared by ChEC and ORGANIC."""
    return truseq_library(
        list(payload), name, dual_index=False, i7_nt=6)


def chec_library():
    return henikoff_single_index_library(
        "ChEC-seq library", (role("DNA released by two nearby protein-tethered MNase cuts"),))


def chec_workflow():
    return chain(molecule("living budding yeast expressing a chromatin-protein–MNase fusion", (
        role("chromatin"), role("DNA-bound protein fused to 3xFLAG-MNase"))), (
        ("Wash cells and permeabilize for 5 minutes with 0.1% digitonin",
         molecule("permeabilized cells with inactive tethered MNase", (
             role("accessible intracellular chromatin"), role("protein-tethered MNase"))),
         "Low nuclear calcium keeps MNase largely inactive before the timed reaction. A matched soluble nuclear MNase strain measures accessibility and enzyme-sequence background."),
        ("Add 2 mM calcium and remove timed digestion aliquots",
         molecule("chromatin cut near occupied protein sites", (
             role("distal chromatin"), role("MNase cleavage end"), role("protein-protected binding region"),
             role("MNase cleavage end"), role("distal chromatin"))),
         "Calcium activates the fusion in situ. Closely spaced seconds-to-minutes samples distinguish rapid tether-proximal cleavage from slower mass-action background."),
        ("Quench each aliquot with EDTA/EGTA, SDS and high salt",
         molecule("stopped ChEC digest", (
             role("large background genomic DNA"), role("short two-cut protein-site fragment"))),
         "Chelators remove the essential calcium and SDS disrupts cells and protein complexes, freezing each kinetic endpoint."),
        ("Extract nucleic acid, digest RNA and retain the small-DNA supernatant",
         molecule("size-selected ChEC DNA", (role("short fragment bounded by two MNase cleavage ends"),)),
         "At a 2.5:1 AMPure-bead:sample ratio the large genomic DNA binds; the unbound fraction is recovered, enriching short protected fragments before library construction."),
        ("End-repair, phosphorylate and add one 3-prime dA",
         molecule("A-tailed ChEC fragment", (
             role("repaired MNase boundary"), role("protein-site insert"), role("repaired MNase boundary"), role("3-prime dA"))),
         "T4 DNA polymerase, Klenow and T4 PNK make ligatable blunt phosphorylated ends; Klenow exo-minus then creates the dA overhang."),
        ("Ligate the six-base indexed TruSeq adapter and PCR-enrich",
         molecule("single-index paired-end ChEC-seq library", (
             role("P5 and Read 1 arm"), role("two-cut ChEC insert"), role("Read 2 arm, six-base i7 and P7"))),
         "The low-input paired-end protocol preserves short inserts and uses KAPA polymerase with the published P5 and P7 enrichment primers."),
    ))


def chec_oligos():
    return (
        "TruSeq Universal Adapter — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′",
        "TruSeq Indexed Adapter — 5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT*G-3′",
        f"P5 library-enrichment primer — 5′-{CHEC_P5_PCR}-3′",
        f"P7 library-enrichment primer — 5′-{CHEC_P7_PCR}-3′",
    )


def nice_library():
    """Likely historical single-i7 E7370 product; inference is attached at source."""
    return truseq_library(
        [role("sonicated genomic DNA carrying an internal biotin-dA or biotin-dC label")],
        "NicE-seq library", dual_index=False, i7_nt=6, inferred_adapters=True)


def nice_workflow():
    return chain(molecule("native cells or cells fixed with 1% formaldehyde", (
        role("nuclei containing open and nucleosome-protected chromatin"),)), (
        ("Isolate nuclei with NP-40",
         molecule("intact nuclei with accessible chromatin exposed", (
             role("open chromatin"), role("protein-occluded chromatin"))),
         "Fixed cells are quenched with glycine first; native-cell NicE-seq omits fixation but otherwise uses the same labeling logic."),
        ("Nick exposed CCD sites and perform biotin nick translation",
         molecule("internally labelled open chromatin", (
             role("genomic DNA"), role("Nt.CviPII nick at accessible CCD"),
             role("DNA polymerase I tract containing biotin-dA and/or biotin-dC"), role("genomic DNA"))),
         "Nt.CviPII creates a single-strand 3-prime-OH at accessible CCD sites. DNA polymerase I extends from the nick while its 5-prime-to-3-prime exonuclease displaces the downstream strand, incorporating biotin-14-dATP and biotin-16-dCTP into an internal tract."),
        ("Stop with EDTA/RNase, reverse crosslinks and extract genomic DNA",
         molecule("purified biotin-marked genomic DNA", (
             role("long genomic molecule with sparse internal biotin labels at accessible sites"),)),
         "Proteinase K plus SDS at 65 degrees C removes chromatin and reverses formaldehyde links; phenol/chloroform extraction preserves the covalent nucleotide labels."),
        ("Sonicate the labelled genome to approximately 150 bp",
         molecule("randomly sheared labelled fragments", (
             role("random shear end"), role("biotin-labelled accessible-chromatin tract"), role("random shear end"))),
         "Unlike DNase-seq, the nick itself is not used as a fragment end: mechanical shearing creates both library boundaries, while internal biotin records accessibility."),
        ("End-repair and add one 3-prime dA to each end",
         da_tailed("dA-tailed labelled genomic fragment", "labelled or unlabelled genomic insert"),
         "The cohesive dA ends are shown before the adapter is introduced."),
        ("Ligate the NEBNext Ultra adapter",
         molecule("adapter-ligated fragments before selection", (
             role("NEBNext adapter"), role("labelled or unlabelled genomic insert"), role("NEBNext adapter"))),
         "The paper uses the discontinued E7370S kit before affinity capture, so every sheared molecule can receive adapters."),
        ("Capture internally biotinylated molecules on streptavidin beads",
         molecule("bead-bound accessible-chromatin library molecules", (
             role("streptavidin-bound internal biotin"), role("adapter-ligated genomic fragment"))),
         "Two-molar-NaCl B&W buffer and four stringent washes remove unlabelled adapter-ligated background; amplification is performed from the bead-bound fraction."),
        ("INFERRED — amplify ten cycles with a compatible historical NEBNext index set",
         molecule("single-index NicE-seq library", (
             role("P5 and Read 1 arm", inferred=True), role("biotin-enriched genomic insert"),
             role("Read 2 arm, sample index and P7", inferred=True))),
         "The paper gives E7370S and the PCR-cycle count but not the separately purchased NEBNext oligo-set catalogue number. The page shows the compatible six-base single-i7 geometry as an inference, not as a transcribed oligo choice."),
    ))


def organic_library():
    return henikoff_single_index_library(
        "ORGANIC profiling library",
        (role("native MNase footprint retained with the affinity-purified factor"),))


def organic_workflow():
    return chain(molecule("unfixed cells expressing an epitope-tagged factor or carrying a specific antibody epitope", (
        role("native chromatin with factor-bound and unbound DNA"),)), (
        ("Isolate nuclei without formaldehyde",
         molecule("native nuclei", (role("factor–DNA complexes preserved by low-salt conditions"),)),
         "ORGANIC deliberately omits crosslinking and sonication; cold buffers and the absence of ATP limit active redistribution during preparation."),
        ("Titrate MNase to a mostly mononucleosome ladder with shorter fragments",
         molecule("native chromatin trimmed at exposed DNA", (
             role("MNase boundary"), role("protein-protected DNA footprint"), role("MNase boundary"),
             role("nucleosome-protected DNA"))),
         "MNase acts endonucleolytically and exonucleolytically until bound protein blocks it. Under-digestion retains accessibility bias; over-digestion can disrupt complexes, so the ladder is a protocol checkpoint."),
        ("Solubilize digested chromatin in low salt",
         molecule("soluble native protein–DNA complexes", (
             role("epitope-tagged or antibody-recognized factor"), role("MNase-bounded protected fragment"))),
         "The factor remains noncovalently associated with its protected DNA while insoluble material is removed."),
        ("Affinity-purify the factor and its bound DNA",
         molecule("ORGANIC immunoprecipitate", (
             role("antibody or tag-affinity matrix"), role("factor"), role("short protected DNA fragment"))),
         "A matched total-digest input retains the MNase background. Stringent washes enrich stable direct protein–DNA complexes without a crosslink-reversal step."),
        ("Remove protein and recover 3–100 ng native footprint DNA",
         molecule("purified paired-boundary footprints", (
             role("MNase boundary"), role("20–40-bp single-factor footprint or longer co-occupied fragment"), role("MNase boundary"))),
         "Proteinase K and DNA purification leave the original paired MNase endpoints; fragment length is biological signal and must not be discarded."),
        ("End-repair and add one 3-prime dA to each end",
         da_tailed("dA-tailed native footprint", "native MNase footprint"),
         "The short footprint remains bounded by its two repaired MNase endpoints."),
        ("Ligate indexed adapters and PCR-enrich",
         molecule("single-index paired-end ORGANIC library", (
             role("P5 and Read 1 arm"), role("native MNase footprint"),
             role("Read 2 arm, six-base i7 and P7"))),
         "The authors' low-input workflow is tuned to retain very short inserts and uses the published phosphorothioate-protected TruSeq adapters and KAPA enrichment PCR."),
    ))


def _finish_henikoff_library(graph: ProtocolGraph, prefix: str, template, label: str,
                              final, primers, source: str, reviewed: list[str],
                              read_lengths) -> None:
    """Published low-input end repair through paired-end library construction."""
    repair = graph.input(prefix + "_repair_enzymes", label + " end-repair enzymes",
                         MaterialKind.CATALYST)
    repaired = graph.product(prefix + "_repaired", label + " end-repaired DNA",
                             MaterialKind.MOLECULAR, state=template.state)
    graph.react(prefix + "_repair", "End-repair " + label + " DNA",
                StepKind.END_REPAIR,
                ((template, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))
    kinase = graph.input(prefix + "_kinase", label + " T4 polynucleotide kinase",
                         MaterialKind.CATALYST)
    phosphorylated = graph.product(prefix + "_phosphorylated", label + " phosphorylated DNA",
                                   MaterialKind.MOLECULAR, state=repaired.state)
    graph.react(prefix + "_phosphorylate", "Phosphorylate repaired " + label + " ends",
                StepKind.PHOSPHORYLATION,
                ((repaired, InputRole.TEMPLATE), (kinase, InputRole.CATALYST)),
                (phosphorylated,))
    tailing = graph.input(prefix + "_tailing_enzyme", label + " Klenow exo-minus",
                          MaterialKind.CATALYST)
    tailed = graph.product(prefix + "_tailed", label + " dA-tailed DNA",
                           MaterialKind.MOLECULAR,
                           state=da_tailed(label + " dA-tailed DNA", label + " insert"))
    graph.react(prefix + "_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((phosphorylated, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)),
                (tailed,))
    adapters = graph.input(prefix + "_adapters", label + " published indexed TruSeq adapters",
                           MaterialKind.ADAPTER)
    ligase = graph.input(prefix + "_adapter_ligase", label + " adapter ligase",
                         MaterialKind.CATALYST)
    adapted = graph.product(prefix + "_adapted", label + " adapter-ligated DNA",
                            MaterialKind.MOLECULAR, state=final)
    graph.react(prefix + "_ligation", "Ligate indexed TruSeq adapters",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,), geometry=HENIKOFF_DA_LIGATION)
    p5 = graph.input(prefix + "_p5_primer", label + " published P5 enrichment primer",
                     MaterialKind.OLIGO)
    p7 = graph.input(prefix + "_p7_primer", label + " published P7 enrichment primer",
                     MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_pcr_polymerase", label + " KAPA PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    geometry = PCRGeometry.exact(
        final, PCRPrimer("P5 enrichment primer", CHEC_P5_PCR),
        PCRPrimer("P7 enrichment primer", CHEC_P7_PCR),
        statement="The cited Henikoff procedure publishes both enrichment-primer sequences.")
    graph.react(prefix + "_pcr", "PCR-enrich the " + label + " library", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=geometry)
    graph.sequence(library, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers), read_lengths))
    reviewed.extend((
        prefix + "_repair_enzymes", prefix + "_repair", prefix + "_kinase",
        prefix + "_phosphorylate", prefix + "_tailing_enzyme", prefix + "_tail",
        prefix + "_adapters", prefix + "_adapter_ligase", prefix + "_ligation",
        prefix + "_p5_primer", prefix + "_p7_primer", prefix + "_pcr_polymerase",
        prefix + "_pcr", prefix + "_library",
    ))


def _chec_graph(final, primers) -> ProtocolGraph:
    source = ("Zentner et al., Nature Communications (2015); doi:10.1038/ncomms9733; "
              "library chemistry: Orsi et al., Current Protocols (2015)")
    graph = ProtocolGraph("ChEC-seq", authority=source)
    cells = graph.input(
        "fusion_cells", "living yeast expressing the chromatin-protein–MNase fusion",
        MaterialKind.BIOLOGICAL,
        state=molecule("living budding yeast expressing a protein–MNase fusion", (
            role("chromatin"), role("DNA-bound protein fused to 3xFLAG-MNase"))),
        primary=True)
    digitonin = graph.input("digitonin", "0.1% digitonin permeabilization condition",
                            MaterialKind.CONDITION)
    chromatin = graph.product(
        "permeabilized_chromatin", "accessible intracellular chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("permeabilized cells", (role("accessible intracellular chromatin"),)))
    tethered_mnase = graph.product("tethered_mnase", "inactive protein-tethered MNase",
                                   MaterialKind.CATALYST)
    graph.react("permeabilize", "Permeabilize cells with digitonin",
                StepKind.PERMEABILIZATION,
                ((cells, InputRole.TEMPLATE), (digitonin, InputRole.CONDITION)),
                (chromatin, tethered_mnase))
    calcium = graph.input("calcium", "2 mM calcium activation condition",
                          MaterialKind.CONDITION)
    cut = graph.product(
        "cut_chromatin", "chromatin cut near occupied protein sites",
        MaterialKind.MOLECULAR,
        state=molecule("chromatin cut near occupied protein sites", (
            role("distal chromatin"), role("MNase cleavage end"),
            role("protein-protected binding region"), role("MNase cleavage end"),
            role("distal chromatin"))))
    graph.react("timed_digest", "Activate tethered MNase and remove timed aliquots",
                StepKind.DIGESTION,
                ((chromatin, InputRole.TEMPLATE), (tethered_mnase, InputRole.CATALYST),
                 (calcium, InputRole.CONDITION)), (cut,))
    quench = graph.input("quench_reagents", "EDTA/EGTA, SDS and high-salt quench",
                         MaterialKind.CONDITION)
    stopped = graph.product(
        "stopped_digest", "stopped ChEC digest", MaterialKind.MOLECULAR,
        state=molecule("stopped ChEC digest", (
            role("large background genomic DNA"), role("short two-cut protein-site fragment"))))
    graph.react("quench", "Quench each timed digestion aliquot", StepKind.QUENCHING,
                ((cut, InputRole.TEMPLATE), (quench, InputRole.CONDITION)), (stopped,))
    extraction = graph.input("extraction_condition", "nucleic-acid extraction and RNase condition",
                             MaterialKind.CONDITION)
    extracted = graph.product("extracted_dna", "extracted ChEC DNA",
                              MaterialKind.MOLECULAR, state=stopped.state)
    graph.react("extract", "Extract nucleic acid and digest RNA", StepKind.EXTRACTION,
                ((stopped, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)),
                (extracted,))
    beads = graph.input("ampure_selection", "2.5:1 AMPure size-selection rule",
                        MaterialKind.SELECTION)
    selected = graph.product(
        "selected_fragments", "short ChEC DNA fragments", MaterialKind.MOLECULAR,
        state=molecule("size-selected ChEC DNA", (
            role("short fragment bounded by two MNase cleavage ends"),)))
    large = graph.product("large_genomic_dna", "large bead-bound genomic DNA",
                          MaterialKind.WASTE)
    graph.react("size_select", "Retain the small-DNA supernatant", StepKind.SELECTION,
                ((extracted, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (selected, large))
    graph.discard(large, "binds the AMPure beads while the short ChEC fragments remain in supernatant")
    reviewed = [
        "fusion_cells", "digitonin", "permeabilize", "calcium", "timed_digest",
        "quench_reagents", "quench", "extraction_condition", "extract",
        "ampure_selection", "size_select",
    ]
    _finish_henikoff_library(
        graph, "chec", selected, "ChEC", final, primers, source, reviewed,
        {"Read 1": 25, "Index 1 (i7)": 6, "Read 2": 25})
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _organic_graph(final, primers) -> ProtocolGraph:
    source = ("Kasinathan et al., Nature Methods (2014); doi:10.1038/nmeth.2766; "
              "library chemistry: Orsi et al., Current Protocols (2015)")
    graph = ProtocolGraph("ORGANIC profiling", authority=source)
    cells = graph.input(
        "native_cells", "unfixed cells containing factor-bound chromatin",
        MaterialKind.BIOLOGICAL,
        state=molecule("unfixed cells", (
            role("native chromatin with factor-bound and unbound DNA"),)), primary=True)
    isolation = graph.input("nuclei_isolation", "native low-salt nuclei-isolation condition",
                            MaterialKind.CONDITION)
    nuclei = graph.product(
        "native_nuclei", "native nuclei preserving factor–DNA complexes",
        MaterialKind.MOLECULAR,
        state=molecule("native nuclei", (role("factor–DNA complexes"),)))
    graph.react("isolate_nuclei", "Isolate nuclei without formaldehyde",
                StepKind.EXTRACTION,
                ((cells, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)), (nuclei,))
    mnase = graph.input("mnase", "titrated micrococcal nuclease", MaterialKind.CATALYST)
    digested = graph.product(
        "digested_chromatin", "native chromatin trimmed at exposed DNA",
        MaterialKind.MOLECULAR,
        state=molecule("native chromatin trimmed at exposed DNA", (
            role("MNase boundary"), role("protein-protected DNA footprint"),
            role("MNase boundary"), role("nucleosome-protected DNA"))))
    graph.react("mnase_digest", "Titrate MNase to a mostly mononucleosome ladder",
                StepKind.DIGESTION,
                ((nuclei, InputRole.TEMPLATE), (mnase, InputRole.CATALYST)), (digested,))
    low_salt = graph.input("low_salt", "low-salt solubilization condition",
                           MaterialKind.CONDITION)
    soluble = graph.product(
        "soluble_complexes", "soluble native protein–DNA complexes",
        MaterialKind.MOLECULAR,
        state=molecule("soluble native protein–DNA complexes", (
            role("epitope-tagged or antibody-recognized factor"),
            role("MNase-bounded protected fragment"))))
    insoluble = graph.product("insoluble_material", "insoluble nuclear material",
                              MaterialKind.WASTE)
    graph.react("solubilize", "Solubilize digested chromatin in low salt",
                StepKind.SOLUBILIZATION,
                ((digested, InputRole.TEMPLATE), (low_salt, InputRole.CONDITION)),
                (soluble, insoluble))
    graph.discard(insoluble, "does not enter the soluble native protein–DNA fraction")
    total_input = graph.product("total_digest_input", "matched total-digest input aliquot",
                                MaterialKind.MOLECULAR, state=soluble.state)
    capture_input = graph.product("capture_input", "complexes assigned to affinity purification",
                                  MaterialKind.MOLECULAR, state=soluble.state)
    graph.react("split_input", "Retain matched total-digest input",
                StepKind.SPLIT, ((soluble, InputRole.TEMPLATE),),
                (total_input, capture_input), changes_molecule=False)
    affinity = graph.input("affinity_matrix", "factor-specific antibody or tag-affinity matrix",
                           MaterialKind.SELECTION)
    immunoprecipitate = graph.product(
        "immunoprecipitate", "factor-associated ORGANIC immunoprecipitate",
        MaterialKind.MOLECULAR,
        state=molecule("ORGANIC immunoprecipitate", (
            role("factor"), role("short protected DNA fragment"))))
    unbound = graph.product("unbound_complexes", "unbound native chromatin complexes",
                            MaterialKind.WASTE)
    graph.react("affinity_capture", "Affinity-purify factor-bound DNA",
                StepKind.AFFINITY_CAPTURE,
                ((capture_input, InputRole.TEMPLATE), (affinity, InputRole.SELECTION)),
                (immunoprecipitate, unbound))
    graph.discard(unbound, "is not stably associated with the selected factor")
    deproteinase = graph.input("deproteinization", "proteinase K and DNA-purification condition",
                              MaterialKind.CONDITION)
    footprints = graph.product(
        "footprints", "purified paired-boundary native footprints",
        MaterialKind.MOLECULAR,
        state=molecule("purified paired-boundary footprints", (
            role("MNase boundary"), role("20–40-bp or longer protected fragment"),
            role("MNase boundary"))))
    graph.react("remove_protein", "Remove protein and recover footprint DNA",
                StepKind.DEPROTEINIZATION,
                ((immunoprecipitate, InputRole.TEMPLATE),
                 (deproteinase, InputRole.CONDITION)), (footprints,))
    input_deproteinase = graph.input("input_deproteinization", "matched-input proteinase K condition",
                                    MaterialKind.CONDITION)
    input_dna = graph.product("input_dna", "purified total-digest input DNA",
                              MaterialKind.MOLECULAR, state=total_input.state)
    graph.react("remove_input_protein", "Remove protein from matched input",
                StepKind.DEPROTEINIZATION,
                ((total_input, InputRole.TEMPLATE),
                 (input_deproteinase, InputRole.CONDITION)), (input_dna,))
    reviewed = [
        "native_cells", "nuclei_isolation", "isolate_nuclei", "mnase", "mnase_digest",
        "low_salt", "solubilize", "split_input", "affinity_matrix", "affinity_capture",
        "deproteinization", "remove_protein", "input_deproteinization",
        "remove_input_protein",
    ]
    lengths = {"Read 1": 25, "Index 1 (i7)": 6, "Read 2": 25}
    _finish_henikoff_library(
        graph, "ip", footprints, "ORGANIC IP", final, primers, source, reviewed, lengths)
    _finish_henikoff_library(
        graph, "input", input_dna, "ORGANIC input", final, primers, source, reviewed, lengths)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def _nice_graph(final, primers) -> ProtocolGraph:
    source = ("Ponnaluri et al., Genome Biology (2017); doi:10.1186/s13059-017-1247-6; "
              "downstream library kit: NEBNext Ultra DNA E7370")
    graph = ProtocolGraph("NicE-seq", authority=source)
    cells = graph.input(
        "cells", "native or formaldehyde-fixed cells",
        MaterialKind.BIOLOGICAL,
        state=molecule("native or fixed cells", (
            role("nuclei containing open and nucleosome-protected chromatin"),)),
        primary=True)
    isolation = graph.input("nuclei_isolation", "NP-40 nuclei-isolation condition",
                            MaterialKind.CONDITION)
    nuclei = graph.product(
        "nuclei", "intact nuclei with accessible chromatin exposed",
        MaterialKind.MOLECULAR,
        state=molecule("intact nuclei", (role("open chromatin"),
                                          role("protein-occluded chromatin"))))
    graph.react("isolate_nuclei", "Isolate nuclei with NP-40", StepKind.EXTRACTION,
                ((cells, InputRole.TEMPLATE), (isolation, InputRole.CONDITION)), (nuclei,))
    nickase = graph.input("nickase", "Nt.CviPII nicking endonuclease",
                          MaterialKind.CATALYST)
    nicked = graph.product(
        "nicked_chromatin", "accessible CCD sites carrying single-strand nicks",
        MaterialKind.MOLECULAR,
        state=molecule("nicked accessible chromatin", (
            role("genomic DNA"), role("Nt.CviPII nick at accessible CCD"),
            role("genomic DNA"))))
    graph.react("nick", "Nick exposed CCD sites", StepKind.DIGESTION,
                ((nuclei, InputRole.TEMPLATE), (nickase, InputRole.CATALYST)), (nicked,))
    biotin_dntps = graph.input("biotin_dntps", "biotin-dATP and/or biotin-dCTP",
                               MaterialKind.MOLECULAR)
    polymerase_i = graph.input("polymerase_i", "DNA polymerase I",
                               MaterialKind.CATALYST)
    labelled = graph.product(
        "labelled_chromatin", "internally biotin-labelled open chromatin",
        MaterialKind.MOLECULAR,
        state=molecule("internally labelled open chromatin", (
            role("genomic DNA"), role("Nt.CviPII nick at accessible CCD"),
            role("DNA polymerase I tract containing biotin-dA and/or biotin-dC"),
            role("genomic DNA"))))
    graph.react("nick_translation", "Perform biotin nick translation",
                StepKind.EXTENSION,
                ((nicked, InputRole.TEMPLATE), (biotin_dntps, InputRole.STRUCTURAL),
                 (polymerase_i, InputRole.CATALYST)), (labelled,))
    reversal = graph.input("crosslink_reversal", "EDTA/RNase, proteinase K and extraction condition",
                           MaterialKind.CONDITION)
    purified = graph.product(
        "purified_dna", "purified biotin-marked genomic DNA",
        MaterialKind.MOLECULAR,
        state=molecule("purified biotin-marked genomic DNA", (
            role("long genomic molecule with sparse internal biotin labels"),)))
    graph.react("purify", "Stop labelling, reverse cross-links and extract DNA",
                StepKind.DEPROTEINIZATION,
                ((labelled, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)),
                (purified,))
    sonication = graph.input("sonication", "approximately 150-bp sonication condition",
                             MaterialKind.CONDITION)
    fragments = graph.product(
        "fragments", "randomly sheared internally labelled fragments",
        MaterialKind.MOLECULAR,
        state=molecule("randomly sheared labelled fragments", (
            role("random shear end"), role("biotin-labelled accessible-chromatin tract"),
            role("random shear end"))))
    graph.react("sonicate", "Sonicate labelled genomic DNA", StepKind.FRAGMENTATION,
                ((purified, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (fragments,))
    repair = graph.input("end_repair_mix", "NEBNext Ultra end-repair enzyme mix",
                         MaterialKind.CATALYST)
    repaired = graph.product("repaired_fragments", "end-repaired labelled fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    graph.react("end_repair", "End-repair labelled fragments", StepKind.END_REPAIR,
                ((fragments, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,))
    tailing = graph.input("tailing_enzyme", "NEBNext Ultra dA-tailing enzyme",
                          MaterialKind.CATALYST)
    tailed = graph.product(
        "tailed_fragments", "dA-tailed labelled genomic fragments",
        MaterialKind.MOLECULAR,
        state=da_tailed("dA-tailed labelled genomic fragment",
                        "labelled or unlabelled genomic insert"))
    graph.react("a_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapter = graph.input("nebnext_adapter", "NEBNext Ultra E7370 adapter",
                          MaterialKind.ADAPTER, inferred=True)
    ligase = graph.input("adapter_ligase", "NEBNext Ultra adapter ligase",
                         MaterialKind.CATALYST)
    adapted = graph.product("adapted_pool", "adapter-ligated labelled and unlabelled fragments",
                            MaterialKind.MOLECULAR, state=final, inferred=True)
    graph.react("adapter_ligation", "Ligate the NEBNext Ultra adapter",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The paper names E7370S but not the separately purchased NEBNext oligo set."))
    streptavidin = graph.input("streptavidin", "streptavidin affinity beads",
                              MaterialKind.SELECTION)
    bound = graph.product("bound_library", "bead-bound internally biotinylated molecules",
                          MaterialKind.MOLECULAR, state=final, inferred=True)
    unlabelled = graph.product("unlabelled_library", "unlabelled adapter-ligated background",
                               MaterialKind.WASTE)
    graph.react("capture", "Capture internally biotinylated molecules",
                StepKind.AFFINITY_CAPTURE,
                ((adapted, InputRole.TEMPLATE), (streptavidin, InputRole.SELECTION)),
                (bound, unlabelled))
    graph.discard(unlabelled, "is removed by stringent streptavidin-bead washes")
    forward = graph.input("pcr_forward", "compatible historical NEBNext forward primer",
                          MaterialKind.OLIGO, inferred=True)
    reverse = graph.input("pcr_reverse", "compatible historical NEBNext indexed reverse primer",
                          MaterialKind.OLIGO, inferred=True)
    pcr_polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                                 MaterialKind.CATALYST)
    library = graph.product("sequencing_library", "single-index NicE-seq library",
                            MaterialKind.LIBRARY, state=final, inferred=True)
    graph.react("library_pcr", "Amplify the bead-bound library for ten cycles",
                StepKind.PCR,
                ((bound, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr_polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The compatible E7370 endpoint is inferred because the oligo-set catalogue number is absent."))
    graph.sequence(library, SequencingPlan(
        final, tuple(primers), tuple(p.role for p in primers)))
    graph.declare_reviewed_scope((
        "cells", "nuclei_isolation", "isolate_nuclei", "nickase", "nick",
        "biotin_dntps", "polymerase_i", "nick_translation", "crosslink_reversal",
        "purify", "sonication", "sonicate", "end_repair_mix", "end_repair",
        "tailing_enzyme", "a_tail", "nebnext_adapter", "adapter_ligase",
        "adapter_ligation", "streptavidin", "capture", "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "sequencing_library",
    ), source)
    return graph.require_publishable()


def get(key: str):
    if key == "organic-profiling":
        final, primers = organic_library()
        spec = SimpleNamespace(
            TITLE="ORGANIC profiling — native affinity-purified chromatin footprints",
            NOTES="01_organic_profiling.html",
            SOURCE='<a href="https://doi.org/10.1038/nmeth.2766">Kasinathan et al., Nature Methods (2014)</a>; exact library chemistry from the authors’ <a href="https://doi.org/10.1002/0471142727.mb2131s110">detailed Current Protocols procedure</a>.',
            SUMMARY="MNase digests unfixed native chromatin until bound proteins block further cleavage. Low-salt solubilization and affinity purification isolate factor-associated footprints whose paired endpoints and insert lengths map direct occupancy.",
            CAVEAT="ORGANIC relies on stable native protein–DNA association: weak or rapidly exchanging complexes can be lost, while MNase sequence preference and digestion depth affect boundaries. A titrated MNase ladder, total-digest input and biological replicates are integral controls.",
            oligos=chec_oligos,
            workflow=organic_workflow,
            FINAL_LIBRARY=final,
            FINAL_CAPTION="Both insert ends are MNase boundaries around DNA that remained protected with the affinity-purified factor. Paired-end geometry is essential because a 20–40-bp single-factor footprint and longer co-occupied or nucleosomal fragments are distinguished by their full insert span.",
            SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": 25, "Index 1 (i7)": 6, "Read 2": 25},
            SEQUENCING_INTRO="The defining study used 25-cycle paired-end HiSeq 2500 sequencing. Standard TruSeq primers read inward from the two native MNase boundaries, while the six-cycle Index 1 read identifies the sample.",
            SEQUENCING_ENDING="",
            SEQUENCING_UNAVAILABLE="",
        )
        spec.protocol_graph = lambda: _organic_graph(final, primers)
        return spec
    if key == "nice-seq":
        final, primers = nice_library()
        spec = SimpleNamespace(
            TITLE="NicE-seq — nick-translation labelling of accessible chromatin",
            NOTES="01_nice_seq.html",
            SOURCE='<a href="https://doi.org/10.1186/s13059-017-1247-6">Ponnaluri et al., Genome Biology (2017)</a>; downstream library reagents from the named NEBNext Ultra DNA E7370 kit.',
            SUMMARY="Nt.CviPII nicks exposed CCD sites inside nuclei and DNA polymerase I converts each accessible nick into a tract containing biotin-dATP and/or biotin-dCTP. After random sonication and adapter ligation, streptavidin selects labelled fragments for PCR.",
            CAVEAT="The defining paper names the discontinued NEBNext Ultra DNA kit E7370S but not the separately purchased NEBNext oligo set. Adapter/index regions are therefore shown as an inferred compatible historical single-i7 construct; the labelling and selection chemistry itself is source-established.",
            oligos=lambda: (),
            workflow=nice_workflow,
            FINAL_LIBRARY=final,
            FINAL_CAPTION="The genomic insert contains an internal biotin-labelled nick-translation tract; its sequenced ends come from later random sonication. <inf>The six-base single-i7 NEBNext shell is a compatible historical E7370 implementation because the paper does not identify its separately purchased oligo set.</inf>",
            SEQ_PRIMERS=primers,
            READ_LENGTHS=None,
            SEQUENCING_INTRO="The deposited data are paired-end NextSeq 500 reads. Standard TruSeq-family Read 1, Index 1 and Read 2 primers bind the inferred E7370-compatible final construct; the paper does not make one genomic-read cycle count part of the method description.",
            SEQUENCING_ENDING="",
            SEQUENCING_UNAVAILABLE="",
        )
        spec.protocol_graph = lambda: _nice_graph(final, primers)
        return spec
    if key != "chec-seq":
        raise KeyError(key)
    final, primers = chec_library()
    spec = SimpleNamespace(
        TITLE="ChEC-seq — calcium-triggered chromatin endogenous cleavage sequencing",
        NOTES="01_chec_seq.html",
        SOURCE='<a href="https://doi.org/10.1038/ncomms9733">Zentner et al., Nature Communications (2015)</a>; library oligos from the cited Henikoff laboratory low-input paired-end protocol.',
        SUMMARY="A protein of interest is genetically fused to MNase. Digitonin admits calcium into living yeast, activating cleavage beside occupied binding sites; short two-cut fragments are selected and converted into a paired-end Illumina library.",
        CAVEAT="ChEC-seq reports calcium-triggered cleavage by a mobile protein–MNase fusion, not an antibody-captured footprint. A free nuclear MNase control and early kinetic samples are essential because slower cleavage contains accessibility, enzyme-bias and transient-collision signal.",
        oligos=chec_oligos,
        workflow=chec_workflow,
        FINAL_LIBRARY=final,
        FINAL_CAPTION="The selected insert is bounded by two MNase-created ends near a protein-binding event. End repair and dA-tailing erase the original overhang geometry but retain both genomic boundary coordinates; the historical construct has one six-base i7 index and no i5 index.",
        SEQ_PRIMERS=primers,
        READ_LENGTHS={"Read 1": 25, "Index 1 (i7)": 6, "Read 2": 25},
        SEQUENCING_INTRO="The defining data set used 25-base paired-end reads. Standard TruSeq Read 1 and Read 2 primers approach the short ChEC insert from opposite ends, and the six-cycle Index 1 read identifies the sample.",
        SEQUENCING_ENDING="",
        SEQUENCING_UNAVAILABLE="",
    )
    spec.protocol_graph = lambda: _chec_graph(final, primers)
    return spec
