"""Targeted chromatin-cleavage and native occupancy library chemistries."""
from __future__ import annotations

from types import SimpleNamespace

from batch_ngs import truseq_library
from protocol_model import chain, molecule, role


# Orsi et al., Current Protocols 2015, Basic Protocol 3.  The asterisk is an
# ordering phosphorothioate; the molecular library constructor below uses the
# corresponding unmodified bases and keeps the chemistry in the displayed oligos.
CHEC_UNIVERSAL_ADAPTER = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
CHEC_INDEXED_ADAPTER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
CHEC_P5_PCR = "AATGATACGGCGACCACCGAG"
CHEC_P7_PCR = "CAAGCAGAAGACGGCATACGAG"


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
        ("End-repair, dA-tail and ligate the NEBNext Ultra adapter",
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
        ("End-repair, dA-tail, ligate indexed adapters and PCR-enrich",
         molecule("single-index paired-end ORGANIC library", (
             role("P5 and Read 1 arm"), role("native MNase footprint"),
             role("Read 2 arm, six-base i7 and P7"))),
         "The authors' low-input workflow is tuned to retain very short inserts and uses the published phosphorothioate-protected TruSeq adapters and KAPA enrichment PCR."),
    ))


def get(key: str):
    if key == "organic-profiling":
        final, primers = organic_library()
        return SimpleNamespace(
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
    if key == "nice-seq":
        final, primers = nice_library()
        return SimpleNamespace(
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
    if key != "chec-seq":
        raise KeyError(key)
    final, primers = chec_library()
    return SimpleNamespace(
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
