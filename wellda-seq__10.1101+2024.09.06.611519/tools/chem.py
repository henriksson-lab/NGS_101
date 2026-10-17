"""wellDA-seq: protocol-owned oligos, final libraries and reaction graph.

Authority: Wang et al., bioRxiv 2024, doi:10.1101/2024.09.06.611519 v1 (posted
2024-09-10; no published version linked by the source fetcher on 2026-10-17).
Methods "Tn5 transposome assembly", "wellDA-seq procedure" Steps 1-3 and the
supplement's Extended Data Fig. 1 (ATAC and DNA library structures) and
Supplementary Tables 5-6 (cell-barcode primer plates).
"""
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import chain, molecule, role

CITATION = "Wang et al., bioRxiv (2024)"
DOI = "10.1101/2024.09.06.611519"

# --- Oligos printed in Methods (verbatim; XXXXXXXX = 8-nt well index) -------------
DEPI_ME_T1 = "GCCTCCCTCGCGCCATAGATGTGTATAAGAGACAG"   # Tn5-1 adapter 1, also ATAC Read 1 primer
DEPI_ME_T2 = "CTTGCCAGCCCGCTCAGAGATGTGTATAAGAGACAG"  # Tn5-1 adapter 2, also ATAC Read 2 primer
ME_REVERSE = "CTGTCTCTTATACACATCT"                   # 5'-/phos/, annealed to both
ATAC_IDX7 = "CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAG"   # custom Index 1 primer (ATAC)
ATAC_IDX5 = "CTGTCTCTTATACACATCTATGGCGCGAGGGAGGC"    # custom Index 2 primer (ATAC)
ATAC_S5_SPACER = "TGA"                               # wellDA_ATAC_S5XX: P5-TGA-index-handle
DNA_IN_LINK = "CTGAGTCGGAGACACGCA"                   # 5' tail of wellDA_DNA_IN_N7XX


def _handle(adapter: str) -> str:
    if not adapter.endswith(nx.ME):
        raise ValueError("DEpi adapter must end in the Tn5 mosaic end")
    return adapter[:-len(nx.ME)]


ATAC_HANDLE1 = _handle(DEPI_ME_T1)   # GCCTCCCTCGCGCCAT
ATAC_HANDLE2 = _handle(DEPI_ME_T2)   # CTTGCCAGCCCGCTCAG
if ME_REVERSE != nx.ME_RC:
    raise ValueError("ME reverse oligo must be the mosaic-end complement")
if (ATAC_IDX7, ATAC_IDX5) != (revcomp(DEPI_ME_T2), revcomp(DEPI_ME_T1)):
    raise ValueError("ATAC index primers must be the reverse complements of the DEpi adapters")


def _cb(name, ident, plate, *, rc=False):
    return Segment(name + (" (reverse complement)" if rc else ""), "X" * 8, "cbc",
                   placeholder=True,
                   feature=feature(ident, "cell_barcode", "combinatorial",
                                   group="well_id", part=name, whitelist=plate))


def _atac_insert():
    return Segment("accessible-chromatin insert (Tn5-1)", "N" * 12, placeholder=True)


def _dna_insert():
    return Segment("genomic insert (Tn5-2, after protease lysis)", "N" * 12, placeholder=True)


def atac_gap_filled() -> Construct:
    """Extended Data Fig. 1a 'Gap filling' state."""
    return Construct([
        Segment("DEpi handle 1", ATAC_HANDLE1, "s5"), Segment("mosaic end (t1)", nx.ME, "me"),
        _atac_insert(),
        Segment("mosaic end (t2, reverse complement)", nx.ME_RC, "me"),
        Segment("DEpi handle 2 (reverse complement)", revcomp(ATAC_HANDLE2), "s7"),
    ], name="gap-filled ATAC fragment")


def atac_library() -> Construct:
    """Extended Data Fig. 1a final ATAC library (on-chip co-enrichment = ATAC enrichment)."""
    return Construct([
        Segment("P5", il.P5, "p5"), Segment("spacer", ATAC_S5_SPACER),
        _cb("ATAC_CB1", "atac_cb1", "Supplementary Table 5 ATAC-CB1 (S5XX, 72)"),
        *atac_gap_filled(),
        _cb("ATAC_CB2", "atac_cb2", "Supplementary Table 5 ATAC-CB2 (N7XX, 72)", rc=True),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="wellDA-seq ATAC library")


def dna_gap_filled() -> Construct:
    """Extended Data Fig. 1b 'Gap filling' state (Illumina TDE1 adapters)."""
    return Construct([
        Segment("Nextera S5", nx.S5, "s5"), Segment("mosaic end (S5 side)", nx.ME, "me"),
        _dna_insert(),
        Segment("mosaic end (S7 side, reverse complement)", nx.ME_RC, "me"),
        Segment("Nextera S7 (reverse complement)", nx.S7_RC, "s7"),
    ], name="gap-filled DNA fragment")


def dna_coindexed() -> Construct:
    """Extended Data Fig. 1b 'DNA&ATAC co-enrichment PCR' product."""
    return Construct([
        Segment("P5", il.P5, "p5"),
        _cb("DNA_CB1", "dna_cb1", "Supplementary Table 6 DNA-CB1 (S5XX, 72)"),
        *dna_gap_filled(),
        _cb("DNA_CB2", "dna_cb2", "Supplementary Table 6 DNA-CB2 (IN_N7XX, 72)", rc=True),
        Segment("DNA_IN link", revcomp(DNA_IN_LINK)),
    ], name="co-indexed DNA product")


def dna_library() -> Construct:
    """Extended Data Fig. 1b final DNA library after the DNA enrichment PCR."""
    return Construct([
        *dna_coindexed(),
        Segment("DNA_CB3 (reverse complement)", "X" * 8, "cbc", placeholder=True,
                feature=feature("dna_cb3", "cell_barcode", "whitelist",
                                whitelist="Supplementary Table 6 enrichment primers (11 OUT_N7XX)",
                                note="Called a DNA cell barcode; its role relative to CB1/CB2 is not stated.")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="wellDA-seq DNA library")


ATAC_PRIMERS = (
    sp.custom("Read 1", "DEpi_ME_t1 (spiked-in custom Read 1)", DEPI_ME_T1, f"{CITATION}, Methods Step 3"),
    sp.custom("Index 1 (i7)", "Idx7 (spiked-in custom Index 1)", ATAC_IDX7, f"{CITATION}, Methods Step 3"),
    sp.custom("Index 2 (i5)", "Idx5 (spiked-in custom Index 2)", ATAC_IDX5, f"{CITATION}, Methods Step 3"),
    sp.custom("Read 2", "DEpi_ME_t2 (spiked-in custom Read 2)", DEPI_ME_T2, f"{CITATION}, Methods Step 3"),
)


def _verified_atac():
    lib = atac_library()
    problems = sp.verify(lib, ATAC_PRIMERS, required_roles=tuple(p.role for p in ATAC_PRIMERS))
    if problems:
        raise ValueError("wellDA-seq ATAC library: " + "; ".join(problems))
    return lib


TITLE = "wellDA-seq"
NOTES = "01_wellda_seq.html"
SOURCE = f'Defining source: <a href="https://doi.org/{DOI}">{CITATION}</a> (preprint v1, Methods and Extended Data Fig. 1).'
SUMMARY = ("Nuclei are bulk-tagmented with a custom Tn5 (Tn5-1, DEpi adapters) to mark open chromatin; single nuclei "
           "in ICELL8 nanowells are protease-lysed and re-tagmented with Illumina Tn5 (Tn5-2) across the whole genome. "
           "One on-chip PCR gives the ATAC and DNA fragments well-specific barcode pairs; pooled products are then "
           "enriched separately per modality.")
CAVEAT = ("Barcodes are drawn as 8-nt placeholders (Supplementary Tables 5–6). Read lengths are not stated in the preprint, "
          "and only the ATAC custom sequencing primers are printed; the DNA library's sequencing primers and index-read "
          "cycles are not modelled.")
FINAL_LIBRARY = _verified_atac()
SEQ_PRIMERS = ATAC_PRIMERS
DNA_LIBRARY = dna_library()
FINAL_CAPTION = ("ATAC library (Extended Data Fig. 1a). ATAC_CB1 and ATAC_CB2 together label one nanowell; the "
                 "DNA library from the same well carries a matched DNA_CB1/DNA_CB2 pair.")
EXTRA_MOLECULAR_STATES = (
    ("Gap-filled ATAC fragment", atac_gap_filled(),
     "Extended Data Fig. 1a: a Tn5-1 fragment after 72 °C gap filling, before the on-chip barcoding PCR."),
    ("Gap-filled DNA fragment", dna_gap_filled(),
     "Extended Data Fig. 1b: a Tn5-2 (Illumina TDE1) fragment after gap filling."),
    ("DNA library", DNA_LIBRARY,
     "DNA library (Extended Data Fig. 1b): Nextera S5/S7 ends from Tn5-2, DNA_CB1 and DNA_CB2 from the on-chip "
     "PCR, then DNA_CB3 and P7 from the DNA enrichment primer. Sequencing primers for this library are not printed."),
)
SEQUENCING_INTRO = ("Both libraries were sequenced on a NextSeq 2000. For the ATAC library the four custom primers "
                    "below were spiked in (Methods Step 3); the DNA library's primers and all read lengths are not stated.")
READ_LENGTHS = None
SEQUENCING_ENDING = ""
SEQUENCING_UNAVAILABLE = ""


def workflow():
    return chain(molecule("permeabilized nucleus", (role("genomic DNA"), role("accessible chromatin"))), (
        ("Bulk-tagment open chromatin with Tn5-1 (DEpi_ME_t1/t2)",
         molecule("ATAC-tagged nucleus", (role("DEpi-tagged accessible fragments"), role("closed chromatin"))), ""),
        ("Lyse single nuclei in nanowells and tagment again with Tn5-2",
         molecule("two tagged fragment classes", (role("DEpi ATAC fragments"), role("Nextera DNA fragments"))), ""),
        ("Gap-fill, four-primer on-chip PCR, pool and enrich per modality",
         molecule("ATAC and DNA libraries", (role("ATAC_CB1/CB2"), role("DNA_CB1/CB2/CB3"))), ""),
    ))


def oligos():
    x = "XXXXXXXX"
    return (
        f"DEpi_ME_t1 — 5′-{DEPI_ME_T1}-3′",
        f"DEpi_ME_t2 — 5′-{DEPI_ME_T2}-3′",
        f"ME reverse — 5′-/5Phos/{ME_REVERSE}-3′",
        f"wellDA_ATAC_S5XX — 5′-{il.P5}{ATAC_S5_SPACER}{x}{ATAC_HANDLE1}-3′",
        f"wellDA_ATAC_N7XX — 5′-{il.P7}{x}{ATAC_HANDLE2}-3′",
        f"wellDA_DNA_S5XX — 5′-{il.P5}{x}{nx.S5}-3′",
        f"wellDA_DNA_IN_N7XX — 5′-{DNA_IN_LINK}{x}{nx.S7}-3′",
        f"Bioo-PCR-F — 5′-{il.P5}-3′",
        f"Bioo-PCR-R — 5′-{il.P7}-3′",
        f"wellDA_DNA_OUT_N7XX — 5′-{il.P7}{x}{DNA_IN_LINK}-3′",
        f"Idx7 (ATAC Index 1) — 5′-{ATAC_IDX7}-3′",
        f"Idx5 (ATAC Index 2) — 5′-{ATAC_IDX5}-3′",
    )


def protocol_graph() -> ProtocolGraph:
    source = f"{CITATION}; doi:{DOI}"
    g = ProtocolGraph(TITLE, authority=source)
    reviewed = []

    def inp(i, label, kind, **kw):
        reviewed.append(i)
        return g.input(i, label, kind, **kw)

    def react(i, label, kind, inputs, outputs, **kw):
        reviewed.append(i)
        reviewed.extend(o.id for o in outputs)
        return g.react(i, label, kind, inputs, outputs, **kw)

    cells = inp("cells", "dissociated cells (cell line or breast tissue)", MaterialKind.BIOLOGICAL,
                state=molecule("cell", (role("genomic DNA"), role("chromatin"))), primary=True)
    lysis1 = inp("perm_buffer", "lysis buffer with NP40, Tween-20, 0.01% digitonin, 4–5 min", MaterialKind.CONDITION)
    nuclei = g.product("nuclei", "permeabilized nuclei", MaterialKind.MOLECULAR,
                       state=molecule("permeabilized nucleus", (role("genomic DNA"), role("accessible chromatin"))))
    react("permeabilize", "Permeabilize cells and wash", StepKind.PERMEABILIZATION,
          ((cells, InputRole.TEMPLATE), (lysis1, InputRole.CONDITION)), (nuclei,))

    tn5_1 = inp("tn5_1", "Tn5-1: DEpi_ME_t1 and DEpi_ME_t2 each annealed to 5′-phos ME reverse, "
                "loaded on Diagenode Tagmentase", MaterialKind.ADAPTER)
    tn5_1_enzyme = inp("tagmentase", "Diagenode Tagmentase (Tn5)", MaterialKind.CATALYST)
    atac_tagged = g.product("atac_tagged", "ATAC-tagged nuclei (DEpi adapters on open chromatin)",
                            MaterialKind.MOLECULAR,
                            state=molecule("ATAC-tagged nucleus", (role("DEpi-tagged accessible fragments"),
                                                                   role("closed chromatin"))))
    react("tagment_atac", "Bulk tagmentation with Tn5-1, 2× TD buffer, 37 °C 30 min", StepKind.TAGMENTATION,
          ((nuclei, InputRole.TEMPLATE), (tn5_1, InputRole.STRUCTURAL), (tn5_1_enzyme, InputRole.CATALYST)),
          (atac_tagged,))

    imaging = inp("cellselect", "ICELL8 dispensing, DAPI imaging and CellSelect single-nucleus calling",
                  MaterialKind.SELECTION)
    wells = g.product("single_nucleus_wells", "selected nanowells with one tagged nucleus",
                      MaterialKind.MOLECULAR, state=atac_tagged.state)
    rejected = g.product("rejected_wells", "empty, multiplet or degraded-nucleus wells", MaterialKind.WASTE)
    react("select_wells", "Dispense into 5,184-nanowell chip and select single-nucleus wells",
          StepKind.SELECTION, ((atac_tagged, InputRole.TEMPLATE), (imaging, InputRole.SELECTION)),
          (wells, rejected))
    g.discard(rejected, "not dispensed with downstream reagents")

    protease = inp("protease", "protease lysis buffer (Tween, Triton X-100), 54.5 °C 30 min then 75.3 °C 15 min",
                   MaterialKind.CATALYST)
    naked = g.product("lysed_genome", "deproteinized genome carrying the Tn5-1 ATAC fragments",
                      MaterialKind.MOLECULAR,
                      state=molecule("deproteinized genome", (role("DEpi-tagged ATAC fragments"),
                                                             role("previously closed genomic DNA"))))
    react("lyse", "Protease lysis removes chromatin proteins in each well", StepKind.DEPROTEINIZATION,
          ((wells, InputRole.TEMPLATE), (protease, InputRole.CATALYST)), (naked,))

    tde1 = inp("tn5_2", "Tn5-2: Illumina TDE1 (20034198) Nextera-adapter transposome", MaterialKind.ADAPTER)
    tde1_enzyme = inp("tde1_enzyme", "Illumina TDE1 Tn5", MaterialKind.CATALYST)
    dual = g.product("dual_tagged", "ATAC (DEpi) and DNA (Nextera) tagged fragments",
                     MaterialKind.MOLECULAR,
                     state=molecule("two tagged classes", (role("DEpi ATAC fragments"), role("Nextera DNA fragments"))))
    react("tagment_dna", "Second tagmentation with Tn5-2, 54.5 °C 10 min", StepKind.TAGMENTATION,
          ((naked, InputRole.TEMPLATE), (tde1, InputRole.STRUCTURAL), (tde1_enzyme, InputRole.CATALYST)),
          (dual,))

    neutral = inp("neutralization", "neutralization mix (EDTA, dNTP) with first-index primers", MaterialKind.CONDITION)
    gapfill_pol = inp("gapfill_polymerase", "KAPA HiFi HotStart polymerase (on-chip master mix, 72.1 °C 8 min)",
                      MaterialKind.CATALYST)
    gap_filled = g.product("gap_filled", "gap-filled ATAC and DNA fragments", MaterialKind.MOLECULAR,
                           state=molecule("gap-filled fragments", (
                               role("DEpi handle 1 · ME · ATAC insert · ME · DEpi handle 2"),
                               role("Nextera S5 · ME · genomic insert · ME · Nextera S7"))))
    react("gap_fill", "Neutralize, then fill the 9-nt Tn5 gaps at 72 °C", StepKind.EXTENSION,
          ((dual, InputRole.TEMPLATE), (neutral, InputRole.CONDITION), (gapfill_pol, InputRole.CATALYST)),
          (gap_filled,))

    s5 = inp("well_s5_primers", "wellDA_ATAC_S5XX + wellDA_DNA_S5XX (one well-specific pair)", MaterialKind.OLIGO)
    n7 = inp("well_n7_primers", "wellDA_ATAC_N7XX + wellDA_DNA_IN_N7XX (one well-specific pair)", MaterialKind.OLIGO)
    chip_pol = inp("chip_polymerase", "KAPA HiFi HotStart polymerase", MaterialKind.CATALYST)
    coindexed = g.product("coindexed", "well-barcoded ATAC and DNA products", MaterialKind.MOLECULAR,
                          state=molecule("co-indexed products", (role("ATAC_CB1 · ATAC fragment · ATAC_CB2"),
                                                                 role("DNA_CB1 · DNA fragment · DNA_CB2"))))
    react("chip_pcr", "On-chip co-enrichment PCR, 12 cycles, four primers per well", StepKind.PCR,
          ((gap_filled, InputRole.TEMPLATE),
           (s5, InputRole.STRUCTURAL), (n7, InputRole.STRUCTURAL), (chip_pol, InputRole.CATALYST)),
          (coindexed,), geometry=PCRGeometry.unavailable(
              atac_library(), reason="One reaction amplifies two adapter families with four primers; the exact "
              "products are drawn as the ATAC library and the co-indexed DNA state."))

    beads = inp("ampure_15", "AMPure XP 1.5×", MaterialKind.SELECTION)
    pooled = g.product("pooled", "pooled chip products", MaterialKind.MOLECULAR, state=coindexed.state)
    pool_waste = g.product("pool_waste", "primers and short products in bead supernatant", MaterialKind.WASTE)
    react("collect", "Collect all wells by centrifugation and purify (AMPure 1.5×)", StepKind.CLEANUP,
          ((coindexed, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (pooled, pool_waste))
    g.discard(pool_waste, "removed with the bead supernatant")

    atac_in = g.product("atac_aliquot", "30 ng pooled product for ATAC enrichment", MaterialKind.MOLECULAR,
                        state=atac_library())
    dna_in = g.product("dna_aliquot", "30 ng pooled product for DNA enrichment", MaterialKind.MOLECULAR,
                       state=dna_coindexed())
    react("split", "Take separate 30 ng aliquots for each modality", StepKind.SPLIT,
          ((pooled, InputRole.TEMPLATE),), (atac_in, dna_in), changes_molecule=False)

    def enrich(prefix, template_mat, template, fwd, rev, label, cycles, final):
        f = inp(prefix + "_fwd", fwd.name, MaterialKind.OLIGO)
        r = inp(prefix + "_rev", rev.name, MaterialKind.OLIGO)
        pol = inp(prefix + "_polymerase", "KAPA HiFi HotStart ReadyMix", MaterialKind.CATALYST)
        amplified = g.product(prefix + "_amplified", label + " enrichment product", MaterialKind.MOLECULAR,
                              state=final)
        react(prefix + "_enrich", f"{label} enrichment PCR ({cycles} cycles)", StepKind.PCR,
              ((template_mat, InputRole.TEMPLATE), (f, InputRole.STRUCTURAL), (r, InputRole.STRUCTURAL),
               (pol, InputRole.CATALYST)), (amplified,),
              geometry=PCRGeometry.exact(template, fwd, rev,
                                         statement=f"Methods Step 3 {label} enrichment primers"))
        b = inp(prefix + "_ampure", "AMPure XP 1.2×", MaterialKind.SELECTION)
        lib = g.product(prefix + "_library", f"wellDA-seq {label} library", MaterialKind.LIBRARY, state=final)
        waste = g.product(prefix + "_waste", "primers and short products in bead supernatant", MaterialKind.WASTE)
        react(prefix + "_cleanup", f"Purify the {label} library (AMPure 1.2×)", StepKind.CLEANUP,
              ((amplified, InputRole.TEMPLATE), (b, InputRole.SELECTION)), (lib, waste))
        g.discard(waste, "removed with the bead supernatant")
        return lib

    p5 = PCRPrimer("Bioo-PCR-F", il.P5)
    atac_lib = enrich("atac", atac_in, atac_library(), p5, PCRPrimer("Bioo-PCR-R", il.P7),
                      "ATAC", 12, FINAL_LIBRARY)
    dna_lib = enrich("dna", dna_in, dna_coindexed(), p5,
                     PCRPrimer("wellDA_DNA_OUT_N7XX", il.P7 + "N" * 8 + DNA_IN_LINK, anneal_nt=len(DNA_IN_LINK)),
                     "DNA", 7, DNA_LIBRARY)
    g.sequence(atac_lib, SequencingPlan(library=FINAL_LIBRARY, primers=ATAC_PRIMERS,
                                        limitation="Read lengths are not stated."))
    g.sequence(dna_lib, SequencingPlan(library=DNA_LIBRARY, limitation=(
        "The preprint prints no sequencing primers or read lengths for the DNA library; bcl2fastq "
        "demultiplexes DNA and ATAC by index without mismatches.")))

    g.annotate(SourceDiscrepancy(
        "well_s5_primers", "Barcode-table numbering",
        (SourceClaim("Preprint Methods Step 2 (4)–(5)",
                     "Index primers are cited as 'Supplementary Table 4-5'."),
         SourceClaim("Supplementary information table list and Extended Data Fig. 1",
                     "Supplementary Table 4 is tumour QC metrics; the ATAC barcodes are Table 5 and the DNA barcodes Table 6 (Extended Data Fig. 1 cites Tables S5 and S6).")),
        "The model cites Supplementary Tables 5 (ATAC) and 6 (DNA).",
        "The supplement's table titles and Extended Data Fig. 1 agree with each other.",
        "Only the cross-reference differs; no barcode sequence is affected."))
    g.annotate(SourceDiscrepancy(
        "atac_enrich", "Extra bases before P7 in the drawn ATAC library",
        (SourceClaim("Supplement Extended Data Fig. 1a, co-enrichment and ATAC enrichment products",
                     "The top strand reads …ATAC_CB2-ATCACGATCTCGTATGCCGTCTTCTGCTTG-3′ (six bases ATCACG between ATAC_CB2 and the P7 complement, with a matching bottom strand)."),
         SourceClaim("Methods Step 2 (5), Extended Data Fig. 1a primer list and Supplementary Table 5 ATAC-CB2 primers",
                     "wellDA_ATAC_N7XX is CAAGCAGAAGACGGCATACGAGAT-XXXXXXXX-CTTGCCAGCCCGCTCAG, so its copy places ATAC_CB2 directly against ATCTCGTATG…")),
        "The model copies the printed N7XX primer, with no ATCACG.",
        "All 72 ATAC-CB2 primers in Supplementary Table 5 and the Methods template lack these bases, and no other printed reagent supplies them.",
        "The six bases look like a drawing artefact, but the source does not say so."))
    g.declare_reviewed_scope(reviewed, source)
    return g.require_publishable()
