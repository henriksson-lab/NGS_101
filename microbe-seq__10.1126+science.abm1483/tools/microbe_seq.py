"""Microbe-seq molecular model.

Authority: Zheng et al., Science 2022 (doi:10.1126/science.abm1483) — not readable in
this audit (publisher paywall; the science.org supplement URL answers 403).  The green
open-access location listed by OpenAlex and Unpaywall, OSTI purl 1959222, was obtained on
2026-10-17 through the Internet Archive and is **the bioRxiv v1 preprint PDF itself**, not
an accepted manuscript, so it adds no published text.  Every reaction
condition and oligo below comes from the linked bioRxiv preprint 10.1101/2020.12.14.422699
v1 supplementary Methods (media-1, text lines 40–113) and Table S6 (media-3, lines
217–700).  They are therefore labelled as preprint-derived and drawn ``inferred``; the
published version may differ.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, feature, revcomp

from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "Microbe-seq — droplet single-microbe genome sequencing"
NOTES = "01_microbe-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1126/science.abm1483">Zheng et al., <i>Science</i> (2022)</a>.'
PREPRINT = ("bioRxiv 10.1101/2020.12.14.422699 v1 supplementary Methods and Table S6 "
            "(preprint-derived; the published Science version was not readable in this audit)")
SUMMARY = "Individual microbes are lysed and MDA-amplified in droplets, tagmented by reagent fusion, then fused with barcode-bead droplets for barcoding PCR."
CAVEAT = ("INFERRED — preprint-derived model. The published <i>Science</i> article and supplement were not readable "
          "(publisher 403; the only open-access copy indexed for this DOI, OSTI purl 1959222, turns out to be the bioRxiv v1 preprint PDF); "
          "conditions, oligos and read configuration come from the linked bioRxiv v1 supplement and may differ in the published version. "
          "Table S6 prints the barcode-1 (8–11 nt, 96 oligos) and barcode-2 (8 nt, 384 oligos) extension oligos; the initiating hydrogel oligo, "
          "the P5PE1 primer and the custom Read-1 primer are not printed. The final library shows one 8-nt barcode-1 representative.")

# Table S6 (preprint), constant arms of the split-and-pool extension oligos, 5′→3′.
BC1_OLIGO_5P = "AAGGCGTCACAAGCAATCACTC"           # Barcode 1 oligo: 5′ linker arm
BC1_OLIGO_3P = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAG"   # 3′ arm, anneals to the hydrogel oligo
BC2_OLIGO_5P = "TCTCTTATACACATCTGACGCTGCCGACGA"   # Barcode 2 oligo: 5′ Tn5-adapter arm
BC2_OLIGO_3P = BC1_OLIGO_5P                       # 3′ arm anneals to the barcode-1 linker copy
BC1_LENGTHS = (8, 9, 10, 11)                      # 24 oligos of each length
BC2_LENGTH = 8

if not nx.ADAPTOR_S5.startswith(revcomp(BC2_OLIGO_5P)):
    raise AssertionError("barcode-2 arm does not copy into the Nextera S5/ME adapter strand")
BEAD_TN5_ARM = revcomp(BC2_OLIGO_5P)               # bead primer 3′ end
TN5_COPIED = nx.ADAPTOR_S5[len(BEAD_TN5_ARM):]     # copied from the tagmented template


def _bead_barcode(part: str, n: int, label: str):
    return seg(label, "B" * n, "cbc", placeholder=True, inferred=True,
               feature=feature(f"microbeseq_bead_{part}", "cell_barcode", "whitelist",
                               whitelist=f"preprint Table S6 barcode-{part[-1]} oligos",
                               group="microbeseq_bead_barcode", part=part))


FINAL_LIBRARY = Construct([
    seg("P5PE1-derived end", "[P5PE1 primer: bases not printed]", "p5", placeholder=True),
    seg("hydrogel-oligo arm", revcomp(BC1_OLIGO_3P), "r1", inferred=True),
    _bead_barcode("bc1", BC1_LENGTHS[0], "barcode 1 (8–11 nt)"),
    seg("barcode linker", revcomp(BC1_OLIGO_5P), inferred=True),
    _bead_barcode("bc2", BC2_LENGTH, "barcode 2"),
    seg("bead Tn5-adapter arm", BEAD_TN5_ARM, "me", inferred=True),
    seg("copied mosaic-end bases", TN5_COPIED, "me", inferred=True),
    seg("single-microbe genomic DNA", "X" * 40, placeholder=True),
    seg("mosaic end (rc)", nx.ME_RC, "me", inferred=True),
    seg("Nextera i7 primer arm (rc)", nx.S7_RC, "r2", inferred=True),
    seg("i7 sample index", "I" * 8, "cbc", placeholder=True, inferred=True,
        feature=feature("sample_index_i7", "sample_index", "unknown")),
    seg("P7 (rc)", il.P7_RC, "p7", inferred=True),
], name="Microbe-seq library (preprint-derived)")
FINAL_CAPTION = ("INFERRED — preprint Table S6 arms. The bead primer carries barcode 1, a 22-nt linker and barcode 2 "
                 "upstream of the Tn5 S5/mosaic-end arm; the Nextera i7 primer adds the single 8-nt sample index.")
SEQ_PRIMERS = ()
SEQUENCING_LIMITATION = ("A custom Read-1 primer (IDT) is required, but its sequence is not printed in the preprint; "
                         "the Read-1 landing site is therefore not computed.")
SEQUENCING_ENDING = ("Preprint supplementary Methods: a custom Read-1 primer is required. In a 100-bp run, Read 1 is 45 bp of barcode only, "
                     "Index 1 is 8 bp and Read 2 contains microbial sequence. In a 300-bp run, Read 1 is 150 bp "
                     "(first 45 bp barcode, last 75 bp microbial, adapter between).")
READ_SPECS = (
    sp.read_spec("Read 1", 45, "custom Read-1 primer: bead barcode 1, linker, barcode 2 (100-bp run; 150 in a 300-bp run)",
                 "extract"),
    sp.read_spec("Index 1 (i7)", 8, "Nextera i7 sample index", "demultiplex"),
    sp.read_spec("Read 2", "remaining", "microbial genomic sequence", "include"),
)


def protocol_graph():
    graph = ProtocolGraph("Microbe-seq", authority=SOURCE)
    microbes = graph.input("microbes", "microbial suspension", MaterialKind.BIOLOGICAL, primary=True)
    loading = graph.input("droplet_loading", "single-microbe droplet loading condition", MaterialKind.CONDITION)
    droplets = graph.product("single_microbe_droplets", "droplets containing individual microbes", MaterialKind.BIOLOGICAL)
    excluded = graph.product("excluded_droplets", "empty or multiplet droplets", MaterialKind.WASTE)
    graph.react("partition", "Partition microbes into lysis droplets", StepKind.PARTITION,
                ((microbes, InputRole.TEMPLATE), (loading, InputRole.CONDITION)), (droplets, excluded))
    graph.discard(excluded, "empty and multiplet droplets do not yield accepted single-microbe libraries")
    lysis = graph.input("lysis_mix", "droplet lysis reagents", MaterialKind.CONDITION)
    dna = graph.product("released_genomes", "lysed single-microbe genomic DNA", MaterialKind.MOLECULAR)
    debris = graph.product("cell_debris", "insoluble microbial debris", MaterialKind.WASTE)
    graph.react("lysis", "Lyse microbes inside droplets", StepKind.SOLUBILIZATION,
                ((droplets, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (dna, debris))
    graph.discard(debris, "insoluble debris does not contribute to the sequencing library")
    random = graph.input("mda_primers", "MDA random primers", MaterialKind.OLIGO,
                         inferred=True, source=PREPRINT + " (MDA reagent primer bases not printed)")
    phi29 = graph.input("phi29", "phi29 DNA polymerase (MDA reagent droplet merge)", MaterialKind.CATALYST)
    mda = graph.product("mda_products", "single-microbe MDA products", MaterialKind.MOLECULAR)
    graph.react("mda", "Merge MDA-reagent droplets; amplify 30 °C 6–8 h, 65 °C 10 min",
                StepKind.EXTENSION,
                ((dna, InputRole.TEMPLATE), (random, InputRole.STRUCTURAL), (phi29, InputRole.CATALYST)),
                (mda,), source=PREPRINT)
    transposome = graph.input("transposome", "Nextera tagmentation reagent (commercial)", MaterialKind.ADAPTER,
                              source=PREPRINT)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagmented_dna", "tagmented MDA fragments with 9-nt gaps", MaterialKind.MOLECULAR)
    twaste = graph.product("tagmentation_waste", "fragments without an S5-adapter end", MaterialKind.WASTE)
    graph.react("tagment", "Merge Nextera droplets and tagment, 55 °C 10 min", StepKind.TAGMENTATION,
                ((mda, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)),
                (tagged, twaste), source=PREPRINT)
    graph.discard(twaste, "fragments lacking the S5-adapter end are not primed by the bead primer")
    bead = graph.input("barcode_bead_primer",
                       "photocleavable hydrogel-bead primer: hydrogel oligo + barcode 1 + linker + barcode 2 + Tn5 arm",
                       MaterialKind.OLIGO, inferred=True,
                       source=PREPRINT + " (Table S6 extension oligos; initiating hydrogel oligo not printed)")
    opposite = graph.input("opposite_primer", "opposite (S7-side) droplet primer: presence and identity not printed", MaterialKind.OLIGO,
                           inferred=True, source=PREPRINT + " (the PCR reagent primer content is not itemised)")
    polymerase = graph.input("droplet_polymerase", "Phusion PCR reagents", MaterialKind.CATALYST)
    barcoded = graph.product("barcoded_fragments", "bead-barcoded genomic amplicons", MaterialKind.MOLECULAR)
    graph.react("barcode_pcr", "Merge bead/PCR droplets: 72 °C 4 min gap fill, then 10 PCR cycles", StepKind.PCR,
                ((tagged, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (opposite, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (barcoded,),
                source=PREPRINT,
                geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                    reason="The hydrogel-oligo bases and the opposite droplet primer are not printed in the preprint."))
    size_filter = graph.input("droplet_size_filter", "droplet size-filter device", MaterialKind.SELECTION)
    kept = graph.product("filtered_droplets", "size-filtered barcoded droplets", MaterialKind.MOLECULAR)
    merged = graph.product("merged_droplets", "enlarged post-PCR merged droplets (possible multi-microbe)", MaterialKind.WASTE)
    graph.react("size_filter", "Remove droplets that merged during PCR", StepKind.SELECTION,
                ((barcoded, InputRole.TEMPLATE), (size_filter, InputRole.SELECTION)), (kept, merged),
                source=PREPRINT)
    graph.discard(merged, "the size-filter device removes most enlarged droplets")
    ampure = graph.input("ampure_1p1x", "PFO emulsion breaking, 1.1× AMPure", MaterialKind.SELECTION)
    pooled = graph.product("purified_amplicons", "purified pooled barcoded amplicons", MaterialKind.MOLECULAR)
    cwaste = graph.product("cleanup_waste", "oil phase and unretained small DNA", MaterialKind.WASTE)
    graph.react("break_purify", "Break the emulsion and purify the aqueous phase", StepKind.CLEANUP,
                ((kept, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (pooled, cwaste), source=PREPRINT)
    graph.discard(cwaste, "removed by phase separation and bead cleanup")
    p5pe1 = graph.input("p5pe1", "P5PE1 primer (IDT; bases not printed)", MaterialKind.OLIGO,
                        source=PREPRINT)
    i7 = graph.input("index_primer_i7", "Nextera i7 index primer (Illumina)", MaterialKind.OLIGO)
    phusion = graph.input("library_polymerase", "Phusion HF polymerase, 5–10 cycles", MaterialKind.CATALYST)
    amplified = graph.product("library_raw", "indexed library before cleanup", MaterialKind.MOLECULAR)
    graph.react("index_pcr", "Add P5 and the Nextera i7 sample index by PCR", StepKind.PCR,
                ((pooled, InputRole.TEMPLATE), (p5pe1, InputRole.STRUCTURAL),
                 (i7, InputRole.STRUCTURAL), (phusion, InputRole.CATALYST)), (amplified,),
                source=PREPRINT,
                geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                    reason="P5PE1 bases and its landing on the hydrogel-oligo arm are not printed."))
    ampure2 = graph.input("ampure_0p8x", "0.8× AMPure", MaterialKind.SELECTION)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY,
                            inferred=True, source=PREPRINT)
    small = graph.product("small_products", "primers and short products below the 0.8× cutoff", MaterialKind.WASTE)
    graph.react("final_cleanup", "Purify the library with 0.8× AMPure", StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (ampure2, InputRole.SELECTION)), (library, small),
                source=PREPRINT)
    graph.discard(small, "not retained by 0.8× AMPure beads")
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, limitation=SEQUENCING_LIMITATION,
        read_uses={"Read 1": "extract", "Index 1 (i7)": "demultiplex", "Read 2": "include"}))
    graph.declare_reviewed_scope((
        "microbes", "droplet_loading", "partition", "lysis_mix", "lysis", "mda_primers",
        "phi29", "mda", "transposome", "tn5", "tagment", "barcode_bead_primer",
        "opposite_primer", "droplet_polymerase", "barcode_pcr", "droplet_size_filter",
        "size_filter", "ampure_1p1x", "break_purify", "p5pe1", "index_primer_i7",
        "library_polymerase", "index_pcr", "ampure_0p8x", "final_cleanup", "library"), PREPRINT)
    return graph.require_publishable()
