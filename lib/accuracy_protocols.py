"""Source-backed models for amplification and consensus-accuracy protocols.

The helpers keep circular consensus, endogenous endpoint families, and explicit UMI
families distinct.  They are not interchangeable merely because all reduce errors.
"""
from __future__ import annotations

from types import SimpleNamespace

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp
from protocol_model import closed_circle, molecule, role


# Source-transcribed facts used by the two protocol models below.  Keeping these
# separate from drawing geometry lets the focused protocol suites check the paper,
# while construction functions continue to own molecular invariants.
TRUEPRIME_PRIMER_LENGTH_NT = (7, 9)
TRUEPRIME_AMPLIFICATION_C = 30
TRUEPRIME_INACTIVATION_C = 65
BOTSEQS_TYPICAL_DILUTION_FOLD = 100_000
BOTSEQS_INSERT_BP = 350
BOTSEQS_READ_BP = 100


def _workflow(initial: Scene, steps) -> Workflow:
    out = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        out.react(action, state.rows(), name=action, note=note)
    return out


def _trueprime_workflow() -> Workflow:
    template = molecule("denatured genomic DNA", (role("single-stranded genomic template"),))
    primed = molecule("de novo primed template", (
        role("single-stranded genomic template"), role("7–9-nt DNA primer")))
    product = molecule("branched multiple-displacement product", (
        role("long genomic copies"), role("displaced single strands"), role("new DNA primers")))
    input_dna = molecule("genomic DNA", (role("genomic DNA"),))
    return _workflow(input_dna, (
        ("Lyse the cell, denature DNA and neutralize", template,
         "The single-cell workflow lyses on ice before neutralization; purified DNA is alkaline-denatured."),
        ("TthPrimPol synthesizes short DNA primers directly from dNTPs", primed,
         "The principal primer products are 7–9 nt; no synthetic random primer is added."),
        ("phi29 DNA polymerase extends the primers and displaces downstream strands", product,
         "Displacement exposes new single-stranded template for another priming round."),
        ("Repeat priming and strand-displacement synthesis at 30 °C", product,
         "The paper used 3 h for WGA and 3 or 6 h for single cells, then 10 min at 65 °C."),
    ))


def _botseqs_workflow() -> Workflow:
    dna = molecule("sheared genomic DNA", (role("genomic fragments with native endpoints"),))
    ligated = molecule("TruSeq PCR-free molecules", (
        role("P5-side adapter"), role("genomic insert"), role("P7-side adapter")),
        junctions=(("P5-side adapter", "genomic insert", "adapter ligation"),
                   ("genomic insert", "P7-side adapter", "adapter ligation")))
    families = molecule("redundantly sequenced bottleneck", (
        role("Watson read family"), role("same two genomic endpoints"), role("Crick read family")))
    return _workflow(dna, (
        ("End-repair, dA-tail and ligate Illumina Y adapters", ligated,
         "Adapter orientation distinguishes Watson and Crick descendants."),
        ("Dilute the adapter-ligated library to about one haploid-genome equivalent", ligated,
         "For 500 ng input the typical dilution was 100,000-fold; it is imposed before PCR."),
        ("PCR-amplify and paired-end sequence each retained molecule redundantly", families,
         "The two genomic endpoints are the endogenous family key; no synthetic UMI is used."),
    ))


def _botseqs_final_rows(library: Construct):
    """Draw the actual adapter-ligation boundaries on the canonical TruSeq molecule."""
    sc = Scene.duplex(list(library), label=library.name,
                      unpaired=("genomic insert'",))
    sc.junction("top", "Read 1 arm", "genomic insert", "adapter ligation")
    sc.junction("top", "dA junction", "Index 1 / Read 2 arm", "adapter ligation")
    sc.labels("top")
    return sc.rows()


def _cirseq_workflow() -> Workflow:
    rna = molecule("purified viral RNA", (role("viral RNA"),))
    fragment = molecule("size-selected RNA", (role("85–100-nt RNA fragment", length_bp=92),))
    circ = closed_circle("circular RNA template", (role("RNA fragment"),),
                         closure_label="RNA intramolecular ligation")
    tandem = molecule("rolling-circle cDNA", (
        role("template copy 1"), role("template copy 2"), role("template copy 3")))
    ds = molecule("double-stranded tandem-repeat cDNA", (role("three or more template copies"),))
    w = Workflow(MolecularState("Starting material", tuple(rna.rows())))
    w.react("Fragment with Zn2+ and size-select", fragment.rows(), name="Fragment and size-select",
            note="Fragments are kept short enough for at least three copies in one read.")
    w.react("5′-phosphorylate and circularize RNA", circ.rows(), name="Circularize RNA",
            note="T4 RNA ligase 1 closes the native 3′→5′ junction.")
    w.react("Random-prime rolling-circle reverse transcription", tandem.rows(),
            name="Rolling-circle reverse transcription",
            note="Physically linked repeats derive independently from the same RNA template.")
    w.react("Synthesize second strand, blunt, dA-tail and ligate TruSeq adapters", ds.rows(),
            name="Complete the Illumina library",
            note="Size selection retains tandem inserts long enough for within-read consensus.")
    return w


def _duplex_repair_workflow() -> Workflow:
    damaged = molecule("damaged or nicked DNA duplex", (role("damaged duplex with ragged ends"),))
    excised = molecule("lesion-excised duplex", (role("short gaps"), role("trimmed overhangs")))
    blunt = molecule("blunt repaired duplex", (role("restricted fill-in"), role("remaining nicks")))
    sealed = molecule("nick-sealed duplex", (role("original duplex interior"),))
    ata = molecule("dA-tailed duplex", (role("original duplex interior"), role("3′ dA", length_bp=1)))
    adapted = molecule("duplex-UMI library", (
        role("custom duplex UMI adapter A", tag="umi",
             molecular_feature=feature("duplex_umi_a", "umi", "random")),
        role("repaired insert"),
        role("custom duplex UMI adapter B", tag="umi",
             molecular_feature=feature("duplex_umi_b", "umi", "random"))),
        junctions=(("custom duplex UMI adapter A", "repaired insert", "adapter ligation"),
                   ("repaired insert", "custom duplex UMI adapter B", "adapter ligation")))
    return _workflow(damaged, (
        ("Excise damaged bases and remove single-stranded overhangs", excised,
         "EndoIV, Fpg, UDG, T4 PDG and EndoVIII act with ExoVII."),
        ("Blunt ends and restrict fill-in", blunt,
         "T4 PNK and T4 DNA polymerase fill short gaps and residual short 5′ overhangs."),
        ("Seal nicks with HiFi Taq DNA ligase", sealed,
         "Nick sealing precedes dA-tailing so polymerase does not start from open nicks."),
        ("Restricted dA-tail with dATP as the only nucleotide", ata,
         "Klenow fragment exo− and Taq add the ligation-compatible 3′ dA."),
        ("Ligate custom dual-index duplex-UMI adapters and PCR-amplify", adapted,
         "Complementary UMI families recover both original strands."),
    ))


TIME_B_PREFIX = "TGGGTGGAGGGTGG"
TIME_BARCODE = "AAAGA"                 # TS301, one published representative
TIME_TN5_REVERSE = "/5Phos/CTGTCTCTTATACACATCT"
TIME_METHYLATED_A = "T/iMe-dC/GT/iMe-dC/GG/iMe-dC/AG/iMe-dC/GT/iMe-dC/AGATGTGTATAAGAGA/iMe-dC/AG"
TIME_B_TS301 = TIME_B_PREFIX + TIME_BARCODE + nx.ME
TIME_BLOCK_B = "CTGTCTCTTATACACATCTHHHHHCCACCCTCCACCCA"
TIME_BLOCK_A = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
TIME_N501 = "AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTC"
TIME_N701 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAAAGCAGTGGTATCAACGCAGATCTGGGTGGAGGGTGG"
TIME_I5 = "TAGATCGC"
TIME_I7 = "TCGCCTTA"
TIME_CUSTOM_R2 = "AAGCAGTGGTATCAACGCAGATC" + TIME_B_PREFIX
TIME_CUSTOM_I1 = revcomp(TIME_CUSTOM_R2)


def _time_library():
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5 (N501)", TIME_I5, "cbc",
            feature=feature("pool_index_i5", "sample_index", "whitelist",
                            whitelist="TIME-seq Supplementary Table 5 indexed PCR primers")),
        seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
        role("bisulfite-converted captured insert"),
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("Tn5 sample barcode reverse complement", revcomp(TIME_BARCODE), "cbc",
            feature=feature("tn5_sample_barcode", "inline_barcode", "whitelist",
                            whitelist="TIME-seq Supplementary Table 5 TS301–TS364")),
        seg("B-adaptor prefix reverse complement", revcomp(TIME_B_PREFIX)),
        seg("custom sequencing handle reverse complement",
            revcomp("AAGCAGTGGTATCAACGCAGATC"), "r2"),
        seg("i7 reverse complement (N701)", revcomp(TIME_I7), "cbc",
            feature=feature("pool_index_i7", "sample_index", "whitelist",
                            whitelist="TIME-seq Supplementary Table 5 indexed PCR primers")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="TIME-seq final library")
    primers = (sp.NEXTERA["R1"],
               sp.custom("Index 1 (i7)", "TIME-seq custom Index 1", TIME_CUSTOM_I1,
                         "Griffin et al. Supplementary Table 5"),
               sp.NEXTERA["I2"],
               sp.custom("Read 2", "TIME-seq custom Read 2", TIME_CUSTOM_R2,
                         "Griffin et al. Supplementary Table 5"))
    errors = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if errors:
        raise ValueError("TIME-seq final library: " + "; ".join(errors))
    return lib, primers


def _time_workflow() -> Workflow:
    dna = molecule("sample genomic DNA", (role("genomic DNA"),))
    tagged = molecule("sample-indexed tagmentation product", (
        role("methylated A adapter"), role("genomic fragment"),
        role("5-nt Tn5 sample barcode", tag="cbc",
             molecular_feature=feature("tn5_sample_barcode", "inline_barcode", "whitelist",
                                       whitelist="TIME-seq Supplementary Table 5 TS301–TS364")),
        role("B adapter")),
        junctions=(("methylated A adapter", "genomic fragment", "Tn5 strand transfer"),
                   ("genomic fragment", "5-nt Tn5 sample barcode", "Tn5 strand transfer")))
    pooled = molecule("pooled indexed DNA", (role("many Tn5-barcoded samples"),))
    captured = molecule("hybrid-captured pool", (role("biotin-RNA bait"), role("target loci")))
    converted = molecule("bisulfite-converted captured DNA", (role("converted target loci"),))
    return _workflow(dna, (
        ("Tagment each sample with methylated A and indexed B transposomes", tagged,
         "The 5-nt B-adapter barcode identifies the sample before pooling."),
        ("Pool samples and perform methylated end repair", pooled,
         "5-methyl-dCTP replaces dCTP so adapter cytosines survive bisulfite treatment."),
        ("Hybridize biotinylated RNA baits and enrich target loci", captured,
         "Capture occurs before bisulfite conversion."),
        ("Bisulfite-convert the captured DNA", converted,
         "Unmethylated cytosines become uracils while protected methylcytosines remain."),
        ("PCR with indexed P5 and P7 primers", molecule("final TIME-seq library", (
            role("P5/i5 end"), role("converted captured insert"), role("Tn5 sample barcode"),
            role("P7/i7 end"))),
         "Pool indices and the internal Tn5 barcode support two-level demultiplexing."),
    ))


def get(key: str):
    if key == "trueprime":
        final = Construct([role("TruePrime-amplified genomic DNA")], name="TruePrime WGA product")
        return SimpleNamespace(
            TITLE="TruePrime — primase-driven whole-genome amplification",
            NOTES="01_trueprime.html",
            SOURCE='<a href="https://doi.org/10.1038/ncomms13296">Picher et al., <i>Nature Communications</i> (2016)</a>.',
            SUMMARY="TthPrimPol makes 7–9-nt DNA primers directly on denatured genomic DNA; phi29 DNA polymerase extends them by multiple-displacement synthesis.",
            CAVEAT="TruePrime ends at amplified DNA. The paper used an unspecified NEBNext library preparation downstream, so no adapter sequence or sequencing-primer landing site is asserted here.",
            oligos=lambda: (), workflow=_trueprime_workflow, FINAL_LIBRARY=final,
            FINAL_CAPTION="The named protocol product is amplified genomic DNA, before sequencing adapters.",
            SEQ_PRIMERS=(), SEQUENCING_INTRO="", READ_LENGTHS=None,
            SEQUENCING_ENDING="The paper fragmented this product, prepared a NEBNext library and performed paired-end HiSeq 2500 sequencing.",
            SEQUENCING_UNAVAILABLE="The publication does not identify the NEBNext kit generation or publish its adapter and sequencing-primer sequences.")
    if key == "botseqs":
        # The 2016 paper establishes standard TruSeq Y adapters but does not establish a
        # dual-index configuration.  Model the period-appropriate single-index molecule.
        lib, primers = truseq_library([role("genomic insert")],
                                      "BotSeqS TruSeq PCR-free library",
                                      dual_index=False)
        return SimpleNamespace(
            TITLE="BotSeqS — bottleneck sequencing system", NOTES="01_botseqs.html",
            SOURCE='<a href="https://doi.org/10.1073/pnas.1607794113">Hoang et al., <i>PNAS</i> (2016)</a>.',
            SUMMARY="A TruSeq PCR-free library is diluted after adapter ligation and before PCR; native fragment endpoints group redundant Watson and Crick read families without a synthetic UMI.",
            CAVEAT="", oligos=lambda: (), workflow=_botseqs_workflow,
            FINAL_LIBRARY=lib, SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": BOTSEQS_READ_BP, "Read 2": BOTSEQS_READ_BP},
            final_rows=lambda: _botseqs_final_rows(lib),
            FINAL_CAPTION="Ordinary TruSeq paired-end molecule; error correction comes from the post-ligation bottleneck and shared genomic endpoints, not an added barcode.",
            SEQUENCING_INTRO="Paired reads establish both endpoints and P5/P7 orientation for Watson- and Crick-family consensus.")
    if key == "cirseq":
        lib, primers = truseq_library([role("tandem-repeat cDNA")], "CirSeq library", dual_index=False)
        return SimpleNamespace(
            TITLE="CirSeq — circular resequencing of RNA virus populations", NOTES="01_cirseq.html",
            SOURCE='<a href="https://doi.org/10.1038/nprot.2014.118">Acevedo and Andino, <i>Nature Protocols</i> (2014)</a>.',
            SUMMARY="Short viral RNA fragments are circularized and copied repeatedly by rolling-circle reverse transcription, placing at least three linked template copies inside one Illumina read.",
            CAVEAT="", oligos=lambda: (), workflow=_cirseq_workflow,
            FINAL_LIBRARY=lib, SEQ_PRIMERS=primers, READ_LENGTHS={"Read 1": 300, "Index 1 (i7)": 8, "Read 2": 300},
            FINAL_CAPTION="TruSeq adapters flank a double-stranded cDNA insert containing tandem copies of one circular RNA fragment.",
            SEQUENCING_INTRO="The protocol used single-end 300-cycle MiSeq reads; the binding view also shows the adapter's documented index and opposite-end primer sites.")
    if key == "duplex-repair":
        final = Construct([
            role("custom left duplex-UMI adapter", tag="umi",
                 molecular_feature=feature("duplex_umi_a", "umi", "random")),
            role("repaired insert"),
            role("custom right duplex-UMI adapter", tag="umi",
                 molecular_feature=feature("duplex_umi_b", "umi", "random")),
        ], name="Duplex-Repair sequencing library")
        return SimpleNamespace(
            TITLE="Duplex-Repair — damage-aware end repair for duplex sequencing", NOTES="01_duplex_repair.html",
            SOURCE='<a href="https://doi.org/10.1093/nar/gkab855">Xiong et al., <i>Nucleic Acids Research</i> (2022)</a>.',
            SUMMARY="Lesion excision, controlled blunting, nick sealing and dATP-only tailing are separated to minimize copying of one-strand damage into an apparent duplex mutation.",
            CAVEAT="The paper identifies custom dual-index duplex-UMI adapters but does not publish their sequences; adapter geometry and primer landing positions therefore remain unspecified.",
            oligos=lambda: (), workflow=_duplex_repair_workflow,
            FINAL_LIBRARY=final, FINAL_CAPTION="Published roles of the custom duplex-UMI library; exact adapter bases are unavailable.",
            SEQ_PRIMERS=(), READ_LENGTHS=None, SEQUENCING_INTRO="",
            SEQUENCING_ENDING="The authors sequenced custom duplex-UMI libraries on HiSeq 2500 or HiSeq X in paired-end mode.",
            SEQUENCING_UNAVAILABLE="The custom dual-index duplex-UMI adapter and primer oligos are not given in the accessible paper or supplement.")
    if key == "time-seq":
        lib, primers = _time_library()
        return SimpleNamespace(
            TITLE="TIME-seq — tagmentation-based indexing for methylation sequencing", NOTES="01_time_seq.html",
            SOURCE='<a href="https://doi.org/10.1038/s43587-023-00555-2">Griffin et al., <i>Nature Aging</i> (2024)</a>, Supplementary Information and Supplementary Table 5.',
            SUMMARY="Bisulfite-resistant Tn5 adapters install a 5-nt sample barcode before pooling; methylated end repair, hybrid capture, bisulfite conversion and indexed PCR produce the targeted methylation library.",
            CAVEAT="The displayed B-adapter uses published representative TS301; the paper publishes 64 alternative 5-nt sample indexes.",
            oligos=lambda: (
                f"<p><b>Tn5 reverse</b>: 5′-{TIME_TN5_REVERSE}-3′</p>",
                f"<p><b>Methylated adaptor A</b>: 5′-{TIME_METHYLATED_A}-3′</p>",
                f"<p><b>Indexed adaptor B, TS301</b>: 5′-{TIME_B_TS301}-3′</p>",
                f"<p><b>Adaptor B blocking oligo</b>: 5′-{TIME_BLOCK_B}-3′</p>",
                f"<p><b>Adaptor A blocking oligo</b>: 5′-{TIME_BLOCK_A}-3′</p>",
                f"<p><b>Indexed PCR primer N501</b>: 5′-{TIME_N501}-3′</p>",
                f"<p><b>Indexed PCR primer N701</b>: 5′-{TIME_N701}-3′</p>",
                f"<p><b>Custom Read 2 sequencing primer</b>: 5′-{TIME_CUSTOM_R2}-3′</p>",
                f"<p><b>Custom Index 1 sequencing primer</b>: 5′-{TIME_CUSTOM_I1}-3′</p>",
            ), workflow=_time_workflow,
            FINAL_LIBRARY=lib, SEQ_PRIMERS=primers,
            READ_LENGTHS={"Read 1": 145, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 5},
            FINAL_CAPTION="The internal five-base Tn5 barcode identifies the pre-pooling sample; i5 and i7 identify the PCR pool. TS301 is drawn as one published representative.",
            SEQUENCING_INTRO="Read 2 and Index 1 use the two published custom sequencing primers; Read 1 and Index 2 use the matching Nextera sites. The economical run described in the Supplementary Information uses 145–153 Read 1 cycles, 8 i7 cycles, optional 8 i5 cycles and 5 Read 2 cycles; this view uses the 145-cycle Read 1 endpoint.")
    raise KeyError(key)
