"""Protocol-owned source, oligos, workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import truseq_library
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, da_tailed, molecule, role
from protein_dna_mapping import henikoff_single_index_library, _finish_henikoff_library

CHEC_UNIVERSAL_ADAPTER = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"

CHEC_INDEXED_ADAPTER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"

CHEC_P5_PCR = "AATGATACGGCGACCACCGAG"

CHEC_P7_PCR = "CAAGCAGAAGACGGCATACGAG"

HENIKOFF_DA_LIGATION = reported_adapter_end(
    "phosphorothioate-protected TruSeq adapter", kind="cohesive",
    polarity="3-prime", terminal_sequence="T",
    statement=("The cited Henikoff low-input procedure publishes the protected TruSeq "
               "adapter and its terminal dT/dA ligation compatibility."))


def organic_library():
    return henikoff_single_index_library(
        "ORGANIC profiling library",
        (role("native MNase footprint retained with the affinity-purified factor"),))


def workflow():
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


def oligos():
    return (
        "TruSeq Universal Adapter — 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′",
        "TruSeq Indexed Adapter — 5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT*G-3′",
        f"P5 library-enrichment primer — 5′-{CHEC_P5_PCR}-3′",
        f"P7 library-enrichment primer — 5′-{CHEC_P7_PCR}-3′",
    )


FINAL_LIBRARY, SEQ_PRIMERS = organic_library()
TITLE = 'ORGANIC profiling — native affinity-purified chromatin footprints'
NOTES = '01_organic_profiling.html'
SOURCE = '<a href="https://doi.org/10.1038/nmeth.2766">Kasinathan et al., Nature Methods (2014)</a>; exact library chemistry from the authors’ <a href="https://doi.org/10.1002/0471142727.mb2131s110">detailed Current Protocols procedure</a>.'
SUMMARY = 'MNase digests unfixed native chromatin until bound proteins block further cleavage. Low-salt solubilization and affinity purification isolate factor-associated footprints whose paired endpoints and insert lengths map direct occupancy.'
CAVEAT = 'ORGANIC relies on stable native protein–DNA association: weak or rapidly exchanging complexes can be lost, while MNase sequence preference and digestion depth affect boundaries. A titrated MNase ladder, total-digest input and biological replicates are integral controls.'
FINAL_CAPTION = 'Both insert ends are MNase boundaries around DNA that remained protected with the affinity-purified factor. Paired-end geometry is essential because a 20–40-bp single-factor footprint and longer co-occupied or nucleosomal fragments are distinguished by their full insert span.'
READ_LENGTHS = {'Read 1': 25, 'Index 1 (i7)': 6, 'Read 2': 25}
SEQUENCING_INTRO = 'The defining study used 25-cycle paired-end HiSeq 2500 sequencing. Standard TruSeq primers read inward from the two native MNase boundaries, while the six-cycle Index 1 read identifies the sample.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''

def protocol_graph() -> ProtocolGraph:
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
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
        graph, "ip", footprints, "ORGANIC IP", final, primers, source, reviewed, lengths,
        HENIKOFF_DA_LIGATION, CHEC_P5_PCR, CHEC_P7_PCR)
    _finish_henikoff_library(
        graph, "input", input_dna, "ORGANIC input", final, primers, source, reviewed, lengths,
        HENIKOFF_DA_LIGATION, CHEC_P5_PCR, CHEC_P7_PCR)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
