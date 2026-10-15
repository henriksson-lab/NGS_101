"""Protocol-owned source, oligos, workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import truseq_library
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import chain, da_tailed, molecule, role
from protein_dna_mapping import henikoff_single_index_library, _finish_henikoff_library


def nice_library():
    """Likely historical single-i7 E7370 product; inference is attached at source."""
    return truseq_library(
        [role("sonicated genomic DNA carrying an internal biotin-dA or biotin-dC label")],
        "NicE-seq library", dual_index=False, i7_nt=6, inferred_adapters=True)


def workflow():
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


FINAL_LIBRARY, SEQ_PRIMERS = nice_library()
TITLE = 'NicE-seq — nick-translation labelling of accessible chromatin'
NOTES = '01_nice_seq.html'
SOURCE = '<a href="https://doi.org/10.1186/s13059-017-1247-6">Ponnaluri et al., Genome Biology (2017)</a>; downstream library reagents from the named NEBNext Ultra DNA E7370 kit.'
SUMMARY = 'Nt.CviPII nicks exposed CCD sites inside nuclei and DNA polymerase I converts each accessible nick into a tract containing biotin-dATP and/or biotin-dCTP. After random sonication and adapter ligation, streptavidin selects labelled fragments for PCR.'
CAVEAT = ("<b>Published polymerase amount conflict.</b> Ponnaluri et al. 2017 Methods "
          "specify 10 U DNA polymerase I in the standard nick-translation reaction, "
          "while the authors’ GSM2579632 extraction metadata specify 50 U for that "
          "deposited sample. The model names DNA polymerase I without assigning a "
          "universal amount; which amount applies to other runs remains unresolved. "
          "The defining paper names the discontinued NEBNext Ultra DNA kit E7370S "
          "but not the separately purchased NEBNext oligo set. Adapter/index regions "
          "are shown as an inferred compatible historical single-i7 construct.")
FINAL_CAPTION = 'The genomic insert contains an internal biotin-labelled nick-translation tract; its sequenced ends come from later random sonication. <inf>The six-base single-i7 NEBNext shell is a compatible historical E7370 implementation because the paper does not identify its separately purchased oligo set.</inf>'
READ_LENGTHS = None
SEQUENCING_INTRO = 'The deposited data are paired-end NextSeq 500 reads. Standard TruSeq-family Read 1, Index 1 and Read 2 primers bind the inferred E7370-compatible final construct; the paper does not make one genomic-read cycle count part of the method description.'
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
def oligos():
    return ()

def protocol_graph() -> ProtocolGraph:
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
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
