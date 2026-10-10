"""Shared molecular constructors for run-on, RNA probing and cell-state assays.

The public methods remain separate protocol pages.  These constructors centralize the
relationships that must not drift between them: a run-on label belongs at the nascent
RNA 3' end, an RT stop is upstream of a chemical adduct, genotype/lineage branches keep
their cellular identifiers, and a proximity report exists only after the two probe
oligos have been covalently joined.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import seqprimers as sp
from batch_ngs import nextera_library, seg
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp


def token(name: str, *, tag: str | None = None, length: int = 18,
          molecular_feature=None, inferred: bool = False) -> Segment:
    return Segment(name, "X" * length, tag, placeholder=True,
                   feature=molecular_feature, inferred=inferred)


def unknown_region(name: str, *, tag: str | None = None, molecular_feature=None) -> Segment:
    """A source-established role whose exact sequence or length is not established."""
    return Segment(name, f"[{name}]", tag, placeholder=True,
                   feature=molecular_feature)


def single(label: str, parts: tuple[Segment, ...], *, note: str = "",
           junctions: tuple[tuple[str, str, str], ...] = ()) -> Scene:
    """One selectable strand with checked, named covalent junctions."""
    sc = Scene()
    sc.strand("molecule", list(parts), label=label)
    for left, right, caption in junctions:
        sc.junction("molecule", left, right, caption)
    if note:
        sc.note("molecule", note)
    sc.labels("molecule")
    return sc


def branches(label: str, rows: tuple[tuple[str, tuple[Segment, ...]], ...]) -> Scene:
    sc = Scene()
    for row, parts in rows:
        sc.strand(row, list(parts), label=label)
        sc.labels(row)
    return sc


def chain(initial_name: str, initial: Scene,
          transitions: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState(initial_name, tuple(initial.rows())))
    for action, state, note in transitions:
        w.react(action, state.rows(), name=action, note=note)
    return w


def run_on(*, chromatin: bool) -> Workflow:
    """Nascent-RNA run-on; the label is structurally fixed at the extended 3' end."""
    substrate = "isolated chromatin" if chromatin else "isolated nuclei"
    labelled = "biotin-NTP run-on end" if chromatin else "BrU run-on end"
    selection = ("streptavidin selects biotin-labelled nascent RNA" if chromatin else
                 "anti-BrdU immunopurification is repeated three times")
    return chain(substrate, single(substrate, (token("engaged RNA polymerase"),
                                               token("nascent RNA"))), (
        ("Allow engaged polymerases to run on; prevent new initiation",
         single("extended nascent RNA", (token("pre-existing nascent RNA"),
                                          token(labelled, tag="umi", length=8))),
         "The nucleotide analogue is incorporated only at the RNA 3′ end extended in vitro."),
        ("Fragment RNA and select analogue-containing molecules",
         single("selected nascent RNA", (token("nascent RNA fragment"),
                                          token(labelled, tag="umi", length=8))), selection),
        ("Install the RNA-library adapter and reverse-transcribe",
         single("nascent-RNA cDNA", (token("library handle"),
                                      token("nascent-RNA-derived cDNA"))),
         "The genomic end of the read reports the position and orientation of engaged polymerase."),
        ("PCR to a sequencing-compatible library",
         single("run-on library", (token("left sequencing arm"),
                                    token("nascent-RNA-derived cDNA"),
                                    token("right sequencing arm"))),
         "The defining source does not expose one universal modern primer set for this library."),
    ))


def icshape() -> Workflow:
    modified = token("NAI-N3 2-prime-O adduct", tag="umi", length=6)
    return chain("RNA in living cells", single("folded cellular RNA", (token("RNA"),)), (
        ("Acylate flexible RNA 2′-hydroxyls with NAI-N3",
         single("chemically probed RNA", (token("RNA upstream"), modified,
                                           token("RNA downstream"))),
         "NAI-N3 supplies both the structure-dependent adduct and an azide enrichment handle."),
        ("Click DIBO-biotin to the azide, fragment, and repair RNA 3′ ends",
         single("biotinylated RNA fragment", (token("RNA fragment"), modified,
                                               token("DIBO-biotin", tag="cbc", length=6))),
         "Copper-free click chemistry makes only NAI-N3-modified fragments selectable."),
        ("Ligate the preadenylated 3′ RNA adapter and enrich on streptavidin",
         single("adapter-ligated selected RNA", (token("RNA fragment"),
                                                  token("3-prime RNA adapter")),
                junctions=(("RNA fragment", "3-prime RNA adapter", "RNA ligation"),)),
         "T4 RNL2 truncated K227Q joins the repaired RNA end to the blocked adapter."),
        ("Reverse-transcribe until the adduct stops synthesis",
         single("barcoded truncated cDNA", (token("13-nt experimental barcode + UMI", tag="cbc", length=13),
                                             token("cDNA ending one base before the adduct"))),
         "The RT primer carries 9 random bases and a 4-base experimental barcode."),
        ("Circularize cDNA with CircLigase II and PCR",
         single("icSHAPE library", (token("left Solexa arm"),
                                     token("13-nt experimental barcode + UMI", tag="cbc", length=13),
                                     token("RT-stop cDNA"), token("right Solexa arm"))),
         "A read's first transcript base marks the RT stop; analysis shifts one nucleotide to the adduct."),
    ))


def genotyping(kind: str) -> Workflow:
    if kind == "got":
        initial = single("10x first-strand cDNA", (
            token("target transcript cDNA"),
            token("10x cell barcode", tag="cbc", length=16,
                  molecular_feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                            whitelist="10x Chromium chemistry-specific whitelist")),
            unknown_region("chemistry-specific UMI", tag="umi",
                  molecular_feature=feature("umi", "umi", "random")),))
        return chain("barcoded 10x cDNA", initial, (
            ("Optionally spike a target primer into RT or cDNA amplification", initial,
             "The paper specifies 5′ RT/additive spike-ins and a 3′ cDNA-amplification spike-in when needed."),
            ("Split amplified cDNA into expression and genotyping branches",
             branches("same cell barcode and UMI", (
                 ("expression branch", (token("whole-transcriptome cDNA"),)),
                 ("genotype branch", (
                     token("target cDNA"),
                     token("cell barcode", tag="cbc", length=16),
                     unknown_region("chemistry-specific UMI", tag="umi",
                                    molecular_feature=feature("umi", "umi", "random")),
                 )),
             )),
             "Only a small cDNA aliquot enters targeted genotyping; the remainder follows standard 10x preparation."),
            ("Target-PCR from upstream of the variant toward the 10x end",
             single("GoT amplicon", (token("variant-containing target"),
                                      token("cell barcode", tag="cbc", length=16),
                                      unknown_region("chemistry-specific UMI", tag="umi",
                                          molecular_feature=feature("umi", "umi", "random")))),
             "The amplicon deliberately retains the cell barcode and UMI."),
            ("Add sequencing handles and sequence with the expression library",
             single("GoT library", (token("sequencing arm"), token("variant-containing target"),
                                     token("cell barcode", tag="cbc", length=16),
                                     unknown_region("chemistry-specific UMI", tag="umi",
                                         molecular_feature=feature("umi", "umi", "random")),
                                     token("indexed sequencing arm"))),
             "Read 1 reports cellular identifiers and Read 2 reaches the selected variant."),
        ))
    if kind == "target":
        start = branches("one FACS-indexed cell", (
            ("genomic molecules", (token("genomic DNA"),)),
            ("RNA molecules", (token("poly(A) RNA"),))))
        preamp = branches("one-tube targeted RT-PCR", (
            ("gDNA amplicons", (token("targeted genomic loci"),)),
            ("cDNA amplicons", (token("targeted coding loci"),)),
            ("transcriptome", (token("Smart-seq+ cDNA"),))))
        return chain("single cell", start, (
            ("Lyse mildly, reverse-transcribe, and multiplex-pre-amplify", preamp,
             "Genomic and transcript-specific primers amplify mutation loci alongside whole-transcriptome cDNA."),
            ("Split the common pre-amplification product", branches("matched cell", (
                ("whole-transcriptome branch", (token("Smart-seq+ cDNA"),)),
                ("targeted genotype branch", (token("gDNA and cDNA target amplicons"),)))),
             "Both branches originate from the same physical well and FACS record."),
            ("Tagment the transcriptome; perform two nested PCRs on genotype targets",
             branches("two sequencing libraries", (
                 ("Nextera transcriptome", (token("Nextera arm"), token("transcript cDNA"), token("Nextera arm"))),
                 ("CS1/CS2 genotype", (token("CS1 site"), token("mutation locus"), token("CS2 site"))))),
             "PCR1 adds CS1/CS2 tails; PCR2 adds the cell-specific Illumina-compatible index."),
        ))
    raise ValueError(kind)


def lineage() -> Workflow:
    return chain("integrated scGESTALT transgene", single("lineage cassette", (
        token("DsRed transcript"), token("CRISPR target array"), token("poly(A)"))), (
        ("Induce Cas9 and guide pools at successive developmental times",
         single("edited lineage cassette", (token("DsRed transcript"),
                                              token("cumulative target-array indels"), token("poly(A)"))),
         "Repair-derived insertions and deletions form a permanent, cumulative lineage record."),
        ("Capture cassette RNA and transcriptome by inDrops",
         single("cell-indexed cDNA", (token("transcript or lineage cDNA"),
                                      token("cell-barcode part 1", tag="cbc", length=8),
                                      token("cell-barcode part 2", tag="cbc", length=8),
                                      token("UMI", tag="umi", length=6))),
         "The inDrops cell identifier is the join key between expression and lineage libraries."),
        ("IVT-amplify, then split fragmented and unfragmented RNA",
         branches("same indexed amplified RNA", (
             ("transcriptome", (token("fragmented amplified RNA"),)),
             ("lineage", (token("unfragmented lineage RNA"),)))),
         "Fragmented RNA makes the expression library; unfragmented RNA preserves the full target array."),
        ("Nested-PCR the target array and add matched Illumina indexes",
         single("scGESTALT lineage library", (token("lineage barcode array"),
                                               token("cell-barcode part 1", tag="cbc", length=8),
                                               token("cell-barcode part 2", tag="cbc", length=8),
                                               token("UMI", tag="umi", length=6), token("library index", tag="cbc", length=8))),
         "The lineage and transcriptome libraries use the same multiplexing index so cells can be matched."),
    ))


def proximity() -> Workflow:
    a = unknown_region("antibody A barcode", tag="cbc",
              molecular_feature=feature("protein_barcode_a", "feature_barcode", "fixed"))
    b = unknown_region("antibody B barcode", tag="cbc",
              molecular_feature=feature("protein_barcode_b", "feature_barcode", "fixed"))
    umi = unknown_region("PLA UMI", tag="umi",
                         molecular_feature=feature("pla_umi", "umi", "random"))
    return chain("cell surface", branches("DNA-conjugated antibodies", (
        ("probe A", (a, token("A connector half"))),
        ("probe B", (token("B connector half"), b, umi, token("poly(A)"))))), (
        ("Bind paired antibody–oligo probes to extracellular proteins",
         branches("two bound targets", (("protein A", (a, token("A connector half"))),
                                         ("protein B", (token("B connector half"), b, umi,
                                                        token("poly(A)"))))),
         "Probe identity is encoded in one fixed barcode on each antibody conjugate."),
        ("Hybridize the universal connector and ligate only proximal probes",
         single("119-base PLA reporter", (a, token("universal connector"), b, umi, token("poly(A)")),
                junctions=(("antibody A barcode", "universal connector", "proximity ligation"),
                           ("universal connector", "antibody B barcode", "proximity ligation"))),
         "Only a ligated A–B product contains both primer sites and is amplifiable."),
        ("Co-capture polyadenylated PLA reporters and cellular mRNA",
         branches("same cell barcode", (("PLA branch", (a, b, umi)),
                                         ("RNA branch", (token("transcript cDNA"),)))),
         "Drop-seq, Smart-seq2 or 10x can supply the single-cell partition and RNA branch."),
        ("PCR each branch and add sequencing indexes",
         single("Prox-seq PLA library", (token("sequencing arm"), a, b, umi,
                                          token("indexed sequencing arm"))),
         "The PLA read reports both protein barcodes and a molecule identifier."),
    ))


def find_seq() -> Workflow:
    return chain("heterogeneous cell suspension", single("one cell", (token("genome"),
                                                                        token("poly(A) RNA"))), (
        ("Co-encapsulate one cell, molten allyl agarose, and acrydite oligo(dT)",
         single("hydrogel bead precursor", (token("cell genome"), token("cell RNA"),
                                             token("covalent oligo(dT) matrix"))),
         "APS/TEMED links acrydite oligo(dT) to allyl agarose."),
        ("Lyse strongly and reverse-transcribe mRNA on the hydrogel",
         single("genome-capture bead", (token("retained genome"), token("matrix-bound cDNA"))),
         "Proteinase K/LiDS exposes DNA while the covalent capture matrix retains the transcriptome."),
        ("Re-encapsulate beads, run marker-specific TaqMan PCR, and sort droplets",
         single("fluorescent selected bead", (token("retained genome"), token("matrix-bound cDNA"),
                                               token("TaqMan marker signal"))),
         "Marker PCR is a selection measurement; it is not the transcriptome library."),
        ("Recover selected cDNA and build Smart-seq2 or Smart-seq3 libraries",
         single("FIND-seq library", (token("Nextera left arm"), token("selected-cell cDNA"),
                                      token("Nextera right arm"))),
         "The 2024 protocol supports bulk bins or individually sorted beads and UMI-aware Smart-seq3."),
    ))


def custom_library(name: str, payload: tuple[Segment, ...]) -> Construct:
    return Construct([token("Illumina-compatible left arm"), *payload,
                      token("Illumina-compatible right arm")], name=name)


TARGET_R1 = "ACACTGACGACATGGTTCTACA"
TARGET_R2 = "TACGGTAGCAGAGACTTGGTCT"


def target_genotype_library() -> tuple[Construct, tuple[sp.SeqPrimer, ...]]:
    lib = Construct([token("P5 and cell index"), Segment("CS1", TARGET_R1, "r1"),
                     token("targeted gDNA or cDNA locus", length=40),
                     Segment("CS2 reverse complement", revcomp(TARGET_R2), "r2"),
                     token("P7")], name="TARGET-seq genotyping library")
    primers = (sp.custom("Read 1", "TARGET-seq CS1-seq", TARGET_R1,
                         "Rodriguez-Meira et al. 2019 Key Resources"),
               sp.custom("Read 2", "TARGET-seq CS2-seq", TARGET_R2,
                         "Rodriguez-Meira et al. 2019 Key Resources"))
    problems = sp.verify(lib, primers, required_roles=("Read 1", "Read 2"))
    if problems:
        raise ValueError("TARGET-seq primer geometry: " + "; ".join(problems))
    return lib, primers


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    notes: str
    summary: str
    endpoint: str
    caveat: str = ""


DEFS = {
    "gro-seq": Definition("GRO-seq — global run-on sequencing", "10.1126/science.1162228",
        "Core, Waterfall &amp; Lis, <i>Science</i> (2008)", "01_gro-seq.html",
        "A nuclear run-on incorporates BrUTP into short extensions made by already-engaged polymerases; repeated anti-BrdU purification selects those nascent RNAs for strand-specific sequencing.",
        "The original paper used 33-base single-end Illumina reads from the 5′ end of fragmented, BrU-labelled RNA; its historical library does not expose a current canonical primer set."),
    "chro-seq": Definition("ChRO-seq — chromatin run-on sequencing", "10.1038/s41588-018-0244-3",
        "Chu et al., <i>Nature Genetics</i> (2018)", "01_chro-seq.html",
        "Run-on transcription in isolated chromatin incorporates biotinylated nucleotides, allowing active polymerases to be mapped even in solid or degraded tissue.",
        "The final strand-specific Illumina read begins beside the run-on-labelled RNA 3′ end; exact sequencing oligos are implementation-specific in the defining paper."),
    "icshape": Definition("icSHAPE — in vivo click SHAPE", "10.1038/nprot.2016.011",
        "Flynn et al., <i>Nature Protocols</i> (2016)", "01_icshape.html",
        "NAI-N3 marks flexible RNA 2′-hydroxyls in cells, DIBO-biotin enables enrichment, and reverse-transcriptase stops encode nucleotide-resolution RNA structure.",
        "A 75-cycle single-end Illumina read begins with the 13-nt inline experimental barcode/UMI and then reports cDNA whose first mapped base locates the RT stop."),
    "got": Definition("GoT — Genotyping of Transcriptomes", "10.1038/s41586-019-1367-0",
        "Nam et al., <i>Nature</i> (2019)", "01_got.html",
        "A targeted amplicon branch is made from 10x cDNA while retaining its cell barcode and UMI, linking a selected transcribed variant to the cell's whole transcriptome.",
        "For 3′ v3, Read 1 uses 28 cycles for cell barcode and UMI, the 8-cycle index selects the library, and the 98- or 130-cycle Read 2 reaches the targeted variant.",
        "The defining paper supports 10x 3′ v2/v3 and 5′ inputs. The schematic shows the common identifier-preserving branch; exact UMI length follows the chosen 10x chemistry."),
    "target-seq": Definition("TARGET-seq — targeted genotyping and transcriptomics", "10.1016/j.molcel.2019.01.009",
        "Rodriguez-Meira et al., <i>Molecular Cell</i> (2019)", "01_target-seq.html",
        "Targeted gDNA and coding-DNA loci are co-amplified with Smart-seq+ whole-transcriptome cDNA from one FACS-indexed cell, then split into genotype and expression libraries.",
        "The genotype library uses published CS1/CS2 custom Read 1/Read 2 primers (151 cycles each) plus a 10-cycle cell index; its matched transcriptome branch is Nextera XT."),
    "scgestalt": Definition("scGESTALT — single-cell lineage recording", "10.1038/s41596-018-0058-x",
        "Raj &amp; Gagnon, <i>Nature Protocols</i> (2018)", "01_scgestalt.html",
        "Cumulative CRISPR indels in an expressed target array record ancestry; inDrops captures that lineage transcript alongside the cell transcriptome and supplies their shared cell identifier.",
        "The lineage run uses Read 1 = 250 cycles for the edited array, Index 1 = 8 cycles for cell-barcode part 1, Index 2 = 8 cycles for library index, and Read 2 = 14 cycles for cell-barcode part 2 plus 6-nt UMI."),
    "prox-seq": Definition("Prox-seq — proximity sequencing", "10.1038/s41592-022-01684-z",
        "Vistain et al., <i>Nature Methods</i> (2022)", "01_prox-seq.html",
        "Pairs of antibody-conjugated DNA probes ligate only when their extracellular protein targets are close; two protein barcodes and a UMI are recovered with the same cell's mRNA.",
        "The complete 119-base PLA reporter is amplified as a distinct library; the exact read/index layout follows the chosen Drop-seq, Smart-seq2 or 10x single-cell implementation.",
        "Prox-seq intentionally supports several scRNA-seq back ends. The schematic therefore records reporter structure and identifier roles without claiming one universal sequencing-primer geometry."),
    "find-seq": Definition("FIND-seq — nucleic-acid cytometry for rare cells", "10.1038/s41596-024-01021-y",
        "Shin et al., <i>Nature Protocols</i> (2024)", "01_find-seq.html",
        "A hydrogel bead retains one cell's genome and oligo(dT)-captured transcriptome while marker-specific TaqMan PCR drives fluorescent droplet sorting before whole-transcriptome library preparation.",
        "Recovered selected-cell cDNA is prepared by Smart-seq2 or UMI-aware Smart-seq3 and sequenced through its validated Nextera library ends."),
}


def get(key: str):
    d = DEFS[key]
    if key == "target-seq":
        final, primers = target_genotype_library()
        reads = {"Read 1": 151, "Read 2": 151}
    elif key == "find-seq":
        final, primers = nextera_library([token("selected-cell transcript cDNA", length=44)],
                                         "FIND-seq transcriptome library")
        reads = None
    elif key == "icshape":
        final = custom_library("icSHAPE library", (
            token("experimental barcode", tag="cbc", length=4,
                  molecular_feature=feature("experimental_barcode", "sample_index", "fixed")),
            token("UMI", tag="umi", length=9,
                  molecular_feature=feature("umi", "umi", "random")),
            token("RT-stop cDNA", length=40),))
        primers = (); reads = None
    elif key == "got":
        final = custom_library("GoT targeted library", (
            token("cell barcode", tag="cbc", length=16,
                  molecular_feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                            whitelist="10x Chromium chemistry-specific whitelist")),
            unknown_region("chemistry-specific UMI", tag="umi",
                  molecular_feature=feature("umi", "umi", "random")),
            token("variant-containing cDNA", length=44),))
        primers = (); reads = None
    elif key == "scgestalt":
        final = custom_library("scGESTALT lineage library", (
            token("lineage barcode array", length=44),
            token("cell-barcode part 1", tag="cbc", length=8,
                  molecular_feature=feature("cell_barcode_1", "cell_barcode", "combinatorial",
                                            group="cell_barcode", part="part 1")),
            token("cell-barcode part 2", tag="cbc", length=8,
                  molecular_feature=feature("cell_barcode_2", "cell_barcode", "combinatorial",
                                            group="cell_barcode", part="part 2")),
            token("UMI", tag="umi", length=6,
                  molecular_feature=feature("umi", "umi", "random")),
            token("library index", tag="cbc", length=8,
                  molecular_feature=feature("sample_index", "sample_index", "fixed")),))
        primers = (); reads = None
    elif key == "prox-seq":
        final = custom_library("Prox-seq PLA library", (
            unknown_region("protein barcode A", tag="cbc",
                  molecular_feature=feature("protein_barcode_a", "feature_barcode", "fixed")),
            unknown_region("protein barcode B", tag="cbc",
                  molecular_feature=feature("protein_barcode_b", "feature_barcode", "fixed")),
            unknown_region("PLA UMI", tag="umi",
                  molecular_feature=feature("pla_umi", "umi", "random")),))
        primers = (); reads = None
    else:
        final = custom_library(d.title + " library", (token("nascent-RNA-derived cDNA", length=44),))
        primers = (); reads = None

    workflow_fn = {"gro-seq": lambda: run_on(chromatin=False),
                   "chro-seq": lambda: run_on(chromatin=True),
                   "icshape": icshape,
                   "got": lambda: genotyping("got"),
                   "target-seq": lambda: genotyping("target"),
                   "scgestalt": lineage,
                   "prox-seq": proximity,
                   "find-seq": find_seq}[key]
    oligos = ()
    if key == "icshape":
        oligos = (
            "<b>3′ adapter</b> — 5′-/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/-3′",
            "<b>RT-primer-1 architecture</b> — 5′-/5Phos/DDDNN AACC NNNN AGATCGGAAGAGCGTCGTGGA /iSp18/ GGATCC /iSp18/ TACTGAACCGC-3′",
        )
    elif key == "target-seq":
        oligos = (f"<b>CS1-seq binding region</b> — 5′-{TARGET_R1}-3′",
                   f"<b>CS2-seq binding region</b> — 5′-{TARGET_R2}-3′")
    unavailable = "" if primers else "Exact sequencing-primer sequences are not fixed by the defining protocol or were not exposed in its accessible primary text; read entry and cycle roles are stated explicitly."
    return SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
        SEQ_PRIMERS=primers, FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint,
        READ_LENGTHS=reads, SEQUENCING_ENDING=d.endpoint,
        SEQUENCING_UNAVAILABLE=unavailable, workflow=workflow_fn,
        oligos=lambda: oligos,
    )
