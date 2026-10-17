"""Linder et al. m6A-miCLIP library and crosslink signatures.

Authority split (bounded pass 2026-10-17-breadth-first-12):
* Linder et al. 2015 Online Methods (PMC4487409 JATS lines 153–173) fix the
  reaction order and conditions: on-bead PNK 3′-dephosphorylation, T4 RNA ligase
  3′-adapter ligation, SDS-PAGE/nitrocellulose transfer with proteinase K release,
  SuperScript III RT, 6% TBE-urea cDNA size selection, iCLIP-style circularization
  and re-linearization, 18–21 PCR cycles, HiSeq 2500.  It prints no oligo.
* Grozhik et al. 2017 (same laboratory, PMC5562447 §2.3 lines 140–182) prints the
  L3 linker, cut oligo, RT1–RT16 and P5/P3 Solexa primers.  Their use in every 2015
  library is not stated, so every oligo-derived segment below is ``inferred``.
"""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Workflow, feature, revcomp
from circular import circularize_ssdna
from circular import CircularizationGeometry
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "miCLIP — single-nucleotide m6A/m6Am mapping"
NOTES = "01_miclip.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.3453">Linder et al. (2015)</a>; detailed open protocol: <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5562447/">Methods in Molecular Biology</a>.'
OLIGO_SOURCE = ('Grozhik et al. 2017 (same laboratory), section 2.3 “Oligonucleotides”; '
                'Linder et al. 2015 print no oligo and delegate the adapter, circularization and '
                're-linearization to the iCLIP protocol, so use of these exact oligos in each 2015 library is inferred.')
SUMMARY = "UV crosslinking fixes an anti-m6A antibody at a modified adenosine. Reverse transcriptase records that protein–RNA adduct as a characteristic substitution or a cDNA truncation before circularization and PCR."
CAVEAT = ("INFERRED oligo geometry: the 2015 article prints no oligonucleotide and refers to the iCLIP protocol. "
          "Adapter, RT-primer, cut-oligo and PCR-primer bases shown here are the same laboratory’s 2017 detailed miCLIP oligos "
          "(Grozhik et al., section 2.3); which RT barcode each 2015 library used is not stated. "
          "Antibody-dependent RT outcomes are alternative molecular products. Read lengths are not stated in the 2015 Methods; "
          "the 50-cycle value is the 2017 protocol’s recommendation (Note 25).")

# Grozhik 2017 §2.3, transcribed 5′→3′ (modifications recorded separately).
L3_LINKER = "AGATCGGAAGAGCGGTTCAG"            # /5rApp/ … /3ddC/, HPLC purified
CUT_OLIGO = "GTTCAGGATCCACGACGCTCTTC"         # … /3ddC/, PAGE purified
P5_SOLEXA = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
P3_SOLEXA = il.PE_PCR_PRIMER_2
RT_FIXED = "AGATCGGAAGAGCGTCGTGGATCCTGAACCGC"  # after /5phos/NN[BBBB]NNN; printed …GTG gatc CTG…
RT_BARCODES = ("aacc", "acaa", "attg", "aggt", "cgcc", "ccgg", "ctaa", "catt",
               "gcca", "gacc", "ggtt", "gtgg", "tccg", "tgcc", "tatt", "ttaa")  # RT1–RT16
INLINE_LAYOUT = (3, 4, 2)  # Read 1: random bases, sample barcode, random bases
BAMHI = "GGATCC"

if P5_SOLEXA != il.TRUSEQ_P5_FULL:
    raise AssertionError("P5 Solexa primer is not the P5 + TruSeq Read 1 arm")
if revcomp(CUT_OLIGO) not in RT_FIXED or BAMHI not in CUT_OLIGO:
    raise AssertionError("cut oligo does not reconstitute the RT-primer BamHI site")
if not RT_FIXED.endswith(revcomp(L3_LINKER)[:9]):
    raise AssertionError("RT primer 3′ end does not pair with the L3 linker")
_CUT = RT_FIXED.index(BAMHI) + 1                 # G^GATCC
RT_LEFT, RT_RIGHT = RT_FIXED[:_CUT], RT_FIXED[_CUT:]
L3_COPIED = revcomp(L3_LINKER)[len(RT_FIXED) - RT_FIXED.index(revcomp(L3_LINKER)[:9]):]
_P3_RC = revcomp(P3_SOLEXA)
if not (_P3_RC.startswith(L3_LINKER) and _P3_RC.endswith(il.P7_RC)):
    raise AssertionError("P3 Solexa primer does not span L3 and P7")
P3_SPACER = _P3_RC[len(L3_LINKER):-len(il.P7_RC)]


def read_identifiers():
    """Read-1 orientation NNN–BBBB–NN (reverse complement of the RT primer's NN-BBBB-NNN)."""
    return [
        seg("UMI bases 1–3", "N" * INLINE_LAYOUT[0], "umi", placeholder=True, inferred=True,
            feature=feature("miclip_umi_5p", "umi", "random", group="miclip_umi", part="bases 1–3")),
        seg("four-base sample barcode", "B" * INLINE_LAYOUT[1], "cbc", placeholder=True, inferred=True,
            feature=feature("inline_sample_barcode", "sample_index", "whitelist",
                            whitelist="RT1–RT16 barcodes (reverse complement in Read 1)")),
        seg("UMI bases 4–5", "N" * INLINE_LAYOUT[2], "umi", placeholder=True, inferred=True,
            feature=feature("miclip_umi_3p", "umi", "random", group="miclip_umi", part="bases 4–5")),
    ]


def rt_primer_identifiers():
    """RT-primer orientation /5phos/NN–BBBB–NNN; the same features as in Read 1."""
    a, b, c = read_identifiers()
    return [seg(c.name, c.top, "umi", placeholder=True, inferred=True, feature=c.feature),
            seg(b.name, b.top, "cbc", placeholder=True, inferred=True, feature=b.feature),
            seg(a.name, a.top, "umi", placeholder=True, inferred=True, feature=a.feature)]


def first_strand_cdna():
    """5′-phosphorylated RT product: RT primer, copied L3 5′ part, then cDNA of the RNA."""
    return Construct([
        *rt_primer_identifiers(),
        seg("RT-primer arm before BamHI cut", RT_LEFT, "r1", inferred=True),
        seg("RT-primer arm after BamHI cut", RT_RIGHT, "r2", inferred=True),
        seg("copy of L3 5′ part", L3_COPIED, "r2", inferred=True),
        seg("RNA-derived cDNA", "D" * 30, placeholder=True),
    ], name="miCLIP first-strand cDNA (5′-P … 3′-OH)")


def linearized_cdna():
    """Circle opened at G^GATCC: the RT-primer arm is split and the barcode moves inward."""
    s = {x.name: x for x in first_strand_cdna()}
    order = ("RT-primer arm after BamHI cut", "copy of L3 5′ part", "RNA-derived cDNA",
             *(x.name for x in rt_primer_identifiers()), "RT-primer arm before BamHI cut")
    return Construct([s[n] for n in order], name="BamHI-linearized cDNA")


_AFTER_INSERT = revcomp(L3_COPIED) + revcomp(RT_RIGHT)
if not _AFTER_INSERT.startswith(L3_LINKER):
    raise AssertionError("linearized cDNA does not carry the full L3 complement")


def pcr_template():
    """Duplex view whose bottom strand is the linearized cDNA (top = its complement)."""
    return Construct([
        seg("Read 1-side arm", revcomp(RT_LEFT), "r1", inferred=True),
        *read_identifiers(),
        seg("RNA insert", "R" * 30, placeholder=True),
        seg("L3 linker", L3_LINKER, "r2", inferred=True),
        seg("BamHI-overhang remnant", _AFTER_INSERT[len(L3_LINKER):], inferred=True),
    ], name="linearized cDNA as PCR template")


FINAL_LIBRARY = Construct([
    # The printed P5 Solexa primer is P5 and the Read 1 arm in one oligo (shared ACAC once).
    seg("P5", P5_SOLEXA[:-len(il.TRUSEQ_READ1)], "p5", inferred=True),
    seg("TruSeq Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=True),
    *read_identifiers(),
    seg("crosslinked RNA insert", "R" * 30, placeholder=True),
    seg("L3 linker", L3_LINKER, "r2", inferred=True),
    seg("P3 Solexa spacer", P3_SPACER, "r2", inferred=True),
    seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
], name="miCLIP library (2017 oligos, inferred for 2015)")
SEQ_PRIMERS = (sp.TRUSEQ["R1"],)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS, required_roles=("Read 1",)):
    raise ValueError("miCLIP final library: " + "; ".join(problems))
READ_LENGTHS = {"Read 1": 50}
SHOW_READ_CONTENT = True
FINAL_CAPTION = ("INFERRED — P5/P3 Solexa primers and L3 from the 2017 protocol. Read 1 begins with NNN–BBBB–NN: "
                 "five random bases form a split UMI and the central four bases are the reverse complement of the RT-primer barcode; "
                 "RNA sequence follows, starting at the base after an RT stop.")
SEQUENCING_INTRO = ("Read 1 lands on the TruSeq Read 1 arm carried by the P5 Solexa primer. Linder et al. used paired-end "
                    "sequencing for the short-insert (~30–40 nt) libraries, but print neither read lengths nor a Read 2 primer; "
                    "Read 2 is therefore not drawn. The 50-cycle length is the 2017 recommendation.")

_initial = Scene(); _initial.strand("RNA", [seg("fragmented poly(A) RNA", "R" * 42, placeholder=True)], label="fragmented poly(A) RNA")
_initial.mark("RNA", "fragmented poly(A) RNA", "30–130 nt fragments containing m6A or cap-proximal m6Am")
INITIAL_ROWS = tuple(_initial.rows()); INITIAL_NAME = "Fragmented poly(A) RNA"


def oligos():
    return [
        "5′-/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/-3′ — L3 linker (HPLC; 2017 §2.3)",
        "5′-GTTCAGGATCCACGACGCTCTTC/3ddC/-3′ — cut oligo (PAGE; 2017 §2.3)",
        "5′-/5phos/NN[BBBB]NNNAGATCGGAAGAGCGTCGTGgatcCTGAACCGC-3′ — RT1–RT16 (BBBB = " + ", ".join(RT_BARCODES) + "; RT1 printed NnaaccNNN)",
        "5′-" + P5_SOLEXA + "-3′ — P5 Solexa PCR primer",
        "5′-" + P3_SOLEXA + "-3′ — P3 Solexa PCR primer",
    ]


def _oligo(name: str, seq: str) -> Construct:
    """Ordering string as printed in Grozhik 2017 (modifications noted in the label)."""
    return Construct([seg(name, seq)], name=name)


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    rna = graph.input("fragmented_rna", "fragmented poly(A) RNA (30–130 nt)", MaterialKind.BIOLOGICAL,
                      state=INITIAL_ROWS, primary=True)
    antibody = graph.input("antibody", "anti-m6A antibody (1–5 µg)", MaterialKind.MOLECULAR)
    uv = graph.input("uv_254", "254-nm UV, 2 × 0.15 J cm⁻²", MaterialKind.CONDITION)
    crosslinked = graph.product("crosslinked_rna", "anti-m6A–RNA covalent adduct",
                                MaterialKind.MOLECULAR)
    graph.react("crosslink", "UV-crosslink anti-m6A to modified RNA", StepKind.CROSS_LINKING,
                ((rna, InputRole.TEMPLATE), (antibody, InputRole.STRUCTURAL),
                 (uv, InputRole.CONDITION)), (crosslinked,))
    beads = graph.input("protein_ag", "Protein A/G beads", MaterialKind.SELECTION)
    enriched = graph.product("enriched_rna", "bead-bound antibody–RNA complexes",
                             MaterialKind.MOLECULAR)
    unbound = graph.product("unbound_rna", "unbound RNA fragments", MaterialKind.WASTE)
    graph.react("immunopurify", "Capture on Protein A/G; high-salt, IP and PNK washes",
                StepKind.AFFINITY_CAPTURE,
                ((crosslinked, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (enriched, unbound))
    graph.discard(unbound, "RNA not retained on the beads is removed in the washes")
    pnk = graph.input("pnk", "T4 polynucleotide kinase, pH 6.5 buffer", MaterialKind.CATALYST)
    dephos = graph.product("dephosphorylated_rna", "on-bead RNA with 3′-OH ends",
                           MaterialKind.MOLECULAR)
    graph.react("dephosphorylate", "Dephosphorylate RNA 3′ ends on beads (30 min)",
                StepKind.DEPHOSPHORYLATION,
                ((enriched, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (dephos,))
    adapter = graph.input("adapter_3p", "L3 linker /5rApp/…/3ddC/ (pre-adenylated)",
                          MaterialKind.ADAPTER, inferred=True, source=OLIGO_SOURCE,
                          state=_oligo("L3 linker", L3_LINKER))
    ligase = graph.input("rna_ligase", "T4 RNA ligase (NEB), overnight", MaterialKind.CATALYST)
    ligated = graph.product("adapter_ligated", "3′-L3-ligated crosslinked RNA",
                            MaterialKind.MOLECULAR)
    graph.react("adapter_ligation", "Ligate the pre-adenylated 3′ linker on beads",
                StepKind.LIGATION,
                ((dephos, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "Single-strand RNA 3′-OH to adenylated linker 5′ end; the 2015 article names neither the linker sequence nor its ends, and the 2017 L3 (/5rApp/, /3ddC/) is only inferred to apply."))
    membrane = graph.input("page_membrane", "NuPAGE SDS-PAGE, nitrocellulose transfer, autoradiography",
                           MaterialKind.SELECTION)
    excised = graph.product("membrane_rna", "heavy/light-chain RNA complexes excised from membrane",
                            MaterialKind.MOLECULAR)
    free_rna = graph.product("membrane_waste", "free linker and non-crosslinked RNA outside the excised regions",
                             MaterialKind.WASTE)
    graph.react("membrane_select", "Elute, SDS-PAGE, transfer and excise antibody-chain regions",
                StepKind.SELECTION,
                ((ligated, InputRole.TEMPLATE), (membrane, InputRole.SELECTION)),
                (excised, free_rna))
    graph.discard(free_rna, "membrane regions outside RNA-crosslinked heavy and light chains are not excised")
    protk = graph.input("proteinase_k", "proteinase K", MaterialKind.CATALYST)
    released = graph.product("released_rna", "RNA with residual crosslinked peptide",
                             MaterialKind.MOLECULAR)
    graph.react("release", "Release RNA with proteinase K; phenol–chloroform extract and precipitate",
                StepKind.DEPROTEINIZATION,
                ((excised, InputRole.TEMPLATE), (protk, InputRole.CATALYST)), (released,))
    rt_primer = graph.input("rt_primer", "RT primer /5phos/NN[BBBB]NNN… (RT1–RT16)",
                            MaterialKind.OLIGO, inferred=True, source=OLIGO_SOURCE)
    rt = graph.input("reverse_transcriptase", "SuperScript III", MaterialKind.CATALYST)
    readthrough = graph.product("readthrough_cdna", "read-through mutation-bearing cDNA",
                                MaterialKind.MOLECULAR)
    stopped = graph.product("truncated_cdna", "crosslink-truncated cDNA", MaterialKind.MOLECULAR)
    graph.react("reverse_transcribe", "Copy across or stop at the residual peptide adduct",
                StepKind.REVERSE_TRANSCRIPTION,
                ((released, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (readthrough, stopped))
    pooled = graph.product("pooled_cdna", "pooled mutation and truncation cDNAs",
                           MaterialKind.MOLECULAR, state=first_strand_cdna(), inferred=True)
    graph.react("pool_rt_outcomes", "Pool both informative RT outcomes", StepKind.POOL,
                ((readthrough, InputRole.TEMPLATE), (stopped, InputRole.TEMPLATE)),
                (pooled,), changes_molecule=False)
    gel = graph.input("tbe_urea_gel", "6% TBE-urea gel", MaterialKind.SELECTION)
    sized = graph.product("sized_cdna", "cDNA 80–120 nt (total RNA) or 70–100 nt (mRNA)",
                          MaterialKind.MOLECULAR, state=first_strand_cdna(), inferred=True)
    off_size = graph.product("off_size_cdna", "cDNA outside the excised window, including unextended RT primer",
                             MaterialKind.WASTE)
    graph.react("size_select", "Size-select first-strand cDNA", StepKind.SELECTION,
                ((pooled, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (sized, off_size))
    graph.discard(off_size, "gel regions outside the excised window are left behind")
    circligase = graph.input("circligase", "CircLigase", MaterialKind.CATALYST)
    closed = circularize_ssdna(first_strand_cdna(), five_prime_phosphate=True)
    circle = graph.product("circular_cdna", "covalently closed first-strand cDNA",
                           MaterialKind.MOLECULAR, state=closed, inferred=True)
    graph.react("circularize", "Circularize first-strand cDNA", StepKind.CIRCULARIZATION,
                ((sized, InputRole.TEMPLATE), (circligase, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.known(
                    "single-strand intramolecular ligation",
                    statement="CircLigase joins the RT primer's 5′ phosphate to the cDNA 3′ hydroxyl (iCLIP procedure)."))
    cut_oligo = graph.input("cut_oligo", "cut oligo …/3ddC/", MaterialKind.OLIGO,
                            inferred=True, source=OLIGO_SOURCE,
                            state=_oligo("cut oligo", CUT_OLIGO))
    cutter = graph.input("restriction_enzyme", "BamHI", MaterialKind.CATALYST, inferred=True,
                         source=OLIGO_SOURCE)
    opened = graph.product("opened_cdna", "BamHI-linearized cDNA", MaterialKind.MOLECULAR,
                           state=linearized_cdna(), inferred=True)
    graph.react("linearize", "Anneal cut oligo and open the circle at the RT-primer BamHI site",
                StepKind.DIGESTION,
                ((circle, InputRole.TEMPLATE), (cut_oligo, InputRole.STRUCTURAL),
                 (cutter, InputRole.CATALYST)), (opened,),
                reason="A single cut opens the circle without excising cDNA; the barcode and UMI move from the 5′ end to an internal position adjacent to the Read 1 arm.")
    fwd = graph.input("pcr_forward", "P5 Solexa PCR primer", MaterialKind.OLIGO,
                      inferred=True, source=OLIGO_SOURCE,
                      state=_oligo("P5 Solexa PCR primer", P5_SOLEXA))
    rev = graph.input("pcr_reverse", "P3 Solexa PCR primer", MaterialKind.OLIGO,
                      inferred=True, source=OLIGO_SOURCE,
                      state=_oligo("P3 Solexa PCR primer", P3_SOLEXA))
    pol = graph.input("pcr_polymerase", "library PCR polymerase, 18–21 cycles", MaterialKind.CATALYST)
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY, inferred=True)
    template = pcr_template()
    graph.react("library_pcr", "PCR-amplify the linearized cDNA (18–21 cycles)", StepKind.PCR,
                ((opened, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
                 (rev, InputRole.STRUCTURAL), (pol, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.exact(
                    template, PCRPrimer("P5 Solexa", P5_SOLEXA, anneal_nt=len(revcomp(RT_LEFT))),
                    PCRPrimer("P3 Solexa", P3_SOLEXA),
                    statement="INFERRED — P5 Solexa 3′ end pairs with the RT-primer arm and P3 Solexa 3′ end with the L3 copy (2017 oligos)."))
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY,
                                           primers=tuple(SEQ_PRIMERS),
                                           read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope((
        "fragmented_rna", "antibody", "uv_254", "crosslink", "protein_ag",
        "immunopurify", "pnk", "dephosphorylate", "adapter_3p", "rna_ligase",
        "adapter_ligation", "page_membrane", "membrane_select", "proteinase_k", "release",
        "rt_primer", "reverse_transcriptase", "reverse_transcribe", "pool_rt_outcomes",
        "tbe_urea_gel", "size_select", "circligase", "circularize", "cut_oligo",
        "restriction_enzyme", "linearize", "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "library"), SOURCE)
    return graph.require_publishable()
