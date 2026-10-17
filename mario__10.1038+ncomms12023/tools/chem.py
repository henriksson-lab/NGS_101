"""Source-owned protocol definition and reaction graph.

Authority: Nguyen et al., Nat. Commun. 7, 12023 (2016), doi:10.1038/ncomms12023,
PMC4931010 JATS Methods (lines 116–261: cross-linking, lysis/RNase I, IPB biotinylation and
T1 immobilization, linker ligation, proximity ligation, T7-exonuclease terminal-linker
removal, RNase III shearing, RT-adapter ligation, RT, C1 pull-down/RNase H elution,
circularization, BamHI relinearization, DP5/DP3 pre-amplification, final PE PCR and the
oligo list) and Supplementary Fig. 1 (circularization strategy).  The ES-1 sample is drawn;
the ES-2/MEF and ES-indirect RT primers differ only in the inline identifiers and are listed
in the oligo table.  rRNA depletion branches (DSN for ES-1/ES-indirect, GeneRead for
ES-2/MEF) and the EGS/formaldehyde ES-indirect cross-link are not drawn.
"""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment, feature, revcomp
from circular import CircularizationGeometry
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)
from protocol_model import closed_circle, molecule, role
from protocol_specs import Definition, get

import mario_chem as mc

DEFINITION = Definition(
    "MARIO — mapping RNA interactome in vivo", "10.1038/ncomms12023",
    "Nguyen et al., Nature Communications (2016)",
    "UV-crosslinked RNPs are biotinylated on protein cysteines and immobilized on streptavidin; "
    "a biotinylated RNA linker is ligated to RNA 5′ ends and proximity-ligated to the 3′ end of a "
    "co-bound RNA, giving RNA1–linker–RNA2 chimeras that are selected through the linker biotin "
    "after reverse transcription.",
    "truseq", ("RNA1 fragment", "MARIO linker", "RNA2 fragment"),
    "Paired-end Illumina HiSeq reads (read length not stated). Read 1 starts with the inline "
    "random/sample barcode (NNNNXXXXNN in the generic Methods form) and then reads the RNA from "
    "the cDNA end; PE PCR Primer 2.0 adds no i7 index.",
    "RNA")

CAVEAT = (
    "<b>Drawn sample:</b> ES-1 (RT primer <code>/5Phos/NNAGGTNNN…</code>). Exact sequences come "
    "from the paper's oligo list; RNA fragments are role tokens. The Read 2 sequencing primer "
    "of the Illumina PE kit is not printed in the paper, so only Read 1 is placed. The "
    "linker-to-read transfer of the internal biotin-dT as T is derived. rRNA depletion "
    "(DSN for ES-1/ES-indirect; GeneRead for ES-2/MEF), the dual-cross-linked ES-indirect "
    "branch and the gel size selections are described but not drawn as separate molecules.")

LINKER_DNA = (mc.RNA_LINKER_LEFT + "T" + mc.RNA_LINKER_RIGHT).replace("U", "T")
ES1_INDEX = "AGGT"
RT_ARM5 = mc.RT_BODY[:19]       # AGATCGGAAGAGCGTCGTG, ends with the BamHI G
RT_ARM3 = mc.RT_BODY[19:]       # GATCCTGAACCGCTCTTCCGATCT
if RT_ARM3[4:] != revcomp(mc.RT_ADAPTER):
    raise AssertionError("RT primer 3′ end must anneal to the printed RT adapter")
if mc.PE1 != il.P5[:-4] + il.TRUSEQ_READ1 or not mc.PE1.endswith(mc.DP5):
    raise AssertionError("PE PCR Primer 1.0 must contain P5, TruSeq Read 1 and DP5")
if not mc.PE2.endswith(mc.DP3) or not mc.PE2.startswith(il.P7):
    raise AssertionError("PE PCR Primer 2.0 must contain P7 and DP3")

_UMI = feature("es1_random_barcode", "umi", "random", group="es1_random_barcode",
               note="Methods call the NN…NNNN random part a 6-nt duplicate-removal barcode")
_INDEX = feature("es1_sample_barcode", "sample_index", "whitelist",
                 whitelist="MARIO RT-primer sample barcodes AGGT/CGCC/CATT")


def _umi(part, n):
    return Segment(f"random barcode {part}", "N" * n, "umi",
                   feature=feature(_UMI.id + "." + part.split()[-1], "umi", "random",
                                   group=_UMI.group, part=part,
                                   note=_UMI.note))


def _linear_cdna() -> Construct:
    """BamHI-relinearized ES-1 cDNA strand, written 5′→3′."""
    return Construct([
        Segment("RT-primer 3′ part after BamHI", RT_ARM3),
        Segment("antisense chimera cDNA", "X" * 20, placeholder=True),
        _umi("part 1", 2),
        Segment("sample barcode", ES1_INDEX, "cbc", feature=_INDEX),
        _umi("part 2", 3),
        Segment("RT-primer 5′ arm", RT_ARM5),
    ], name="relinearized ES-1 cDNA")


def _final() -> tuple[Construct, tuple]:
    right = revcomp(mc.PE2)
    lib = Construct([
        Segment("P5 (shares ACAC with Read 1)", il.P5[:-4], "p5"),
        Segment("TruSeq Read 1 arm (PE 1.0)", il.TRUSEQ_READ1, "s5"),
        _umi("part 2", 3),
        Segment("sample barcode", revcomp(ES1_INDEX), "cbc", feature=_INDEX),
        _umi("part 1", 2),
        Segment("RNA1 fragment", "X" * 14, placeholder=True),
        Segment("MARIO linker copy", LINKER_DNA, "linker"),
        Segment("RNA2 fragment", "X" * 14, placeholder=True),
        Segment("RT adapter copy", right[:20]),
        Segment("PE 2.0 read-2 arm", right[20:-len(il.P7)], "s7"),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name=DEFINITION.title + " final library (ES-1)")
    primers = (sp.TRUSEQ["R1"],)
    problems = sp.verify(lib, primers, required_roles=("Read 1",))
    if problems:
        raise ValueError("; ".join(problems))
    return lib, primers


LINEAR_CDNA = _linear_cdna()
FINAL_LIBRARY, SEQ_PRIMERS = _final()
SPEC = get("mario", definition=DEFINITION, final_override=(FINAL_LIBRARY, SEQ_PRIMERS),
           oligos_fn=mc.oligos,
           sequencing_unavailable=(
               "The paper prints PE PCR Primers 1.0/2.0 but not the Illumina PE Read 2 "
               "sequencing primer; read lengths are not stated."))
SPEC.CAVEAT = CAVEAT


def _state(label, *names):
    return molecule(label, tuple(role(n) for n in names))


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    g = ProtocolGraph(d.title, authority=source)
    rnps = g.input("rnps", "E14 mES, MEF or S2 cells with protein-bound RNAs",
                   MaterialKind.BIOLOGICAL, state=_state("cellular RNPs", "protein-bound RNA"),
                   primary=True)
    uv = g.input("uv", "UV-C 254 nm, 400 mJ cm⁻² on ice", MaterialKind.CONDITION)
    fixed = g.product("fixed", "UV-cross-linked RNPs", MaterialKind.MOLECULAR, state=rnps.state)
    g.react("crosslink", "UV-cross-link RNA to bound proteins", StepKind.CROSS_LINKING,
            ((rnps, InputRole.TEMPLATE), (uv, InputRole.CONDITION)), (fixed,))
    rnase_i = g.input("rnase_i", "TURBO DNase, then RNase I 1:100, 37 °C 3 min",
                      MaterialKind.CATALYST)
    frag = g.product("fragmented", "lysate RNA fragments (5′-OH, 3′-P) on cross-linked proteins",
                     MaterialKind.MOLECULAR,
                     state=_state("fragmented RNP", "5′-HO-RNA fragment-3′-P", "protein"))
    g.react("fragment", "Lyse, DNase-treat and RNase I-fragment RNA", StepKind.FRAGMENTATION,
            ((fixed, InputRole.TEMPLATE), (rnase_i, InputRole.CATALYST)), (frag,))
    ipb = g.input("ipb", "EZ-link iodoacetyl-PEG2-biotin (cysteine biotinylation), after 500 mM "
                  "NaCl/0.3% SDS 65 °C; DTT quench; dialysis", MaterialKind.SELECTION)
    t1 = g.input("t1_beads", "MyOne Streptavidin T1 beads; denaturing/high-salt washes",
                 MaterialKind.SELECTION)
    bound = g.product("immobilized", "bead-immobilized cross-linked RNPs",
                      MaterialKind.MOLECULAR, state=frag.state)
    unbound = g.product("unbound", "supernatant and washes with non-cross-linked RNA",
                        MaterialKind.WASTE)
    g.react("immobilize", "Biotinylate proteins and immobilize RNPs on streptavidin",
            StepKind.AFFINITY_CAPTURE,
            ((frag, InputRole.TEMPLATE), (ipb, InputRole.SELECTION),
             (t1, InputRole.SELECTION)), (bound, unbound))
    g.discard(unbound, "aspirated supernatant and stringent washes")
    pnk = g.input("pnk_minus", "T4 PNK (3′-phosphatase minus) + ATP, 37 °C 1 h",
                  MaterialKind.CATALYST)
    p5rna = g.product("p5_rna", "RNA fragments with 5′-P, 3′-P retained", MaterialKind.MOLECULAR,
                      state=_state("5′-phosphorylated RNA", "5′-P-RNA fragment-3′-P"))
    g.react("phosphorylate_rna", "Phosphorylate RNA 5′ ends without removing 3′ phosphates",
            StepKind.PHOSPHORYLATION,
            ((bound, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (p5rna,))
    linker = g.input("linker", "biotinylated RNA linker (5′-OH), 30 µL 20 µM; exact state = 19-nt RNA arm after /iBiodT/",
                     MaterialKind.ADAPTER, state=mc.RNA_LINKER_RIGHT)
    rnl1 = g.input("rnl1", "T4 RNA ligase 1, PEG8000, 16 °C overnight", MaterialKind.CATALYST)
    linked = g.product("linker_rna", "5′-HO-linker–RNA-3′-P", MaterialKind.MOLECULAR,
                       state=molecule("linker-ligated RNA", (role("biotinylated linker"),
                                                             role("RNA2 fragment")),
                                      junctions=(("biotinylated linker", "RNA2 fragment",
                                                  "linker 3′-OH to RNA 5′-P"),)))
    g.react("ligate_linker", "Ligate the linker 3′ end to RNA 5′ ends", StepKind.LIGATION,
            ((p5rna, InputRole.TEMPLATE), (linker, InputRole.STRUCTURAL),
             (rnl1, InputRole.CATALYST)), (linked,),
            geometry=AdapterLigationGeometry.unavailable(
                "Single-stranded RNA ligation: the linker's 3′-OH joins variable RNA 5′-P ends; "
                "the 5′-OH linker cannot self-ligate. No duplex end geometry applies."))
    pnk_buf = g.input("pnk_ph65", "'T4 PNK (3′-phosphatase minus)' in pH 6.5 PNK buffer without "
                      "ATP, 37 °C 20 min (as printed)", MaterialKind.CATALYST)
    oh3 = g.product("oh3", "linker–RNA with 3′-OH (intended)", MaterialKind.MOLECULAR,
                    state=linked.state)
    g.react("dephosphorylate", "Remove RNA 3′ phosphates", StepKind.DEPHOSPHORYLATION,
            ((linked, InputRole.TEMPLATE), (pnk_buf, InputRole.CATALYST)), (oh3,))
    pnk2 = g.input("pnk_linker", "T4 PNK (3′-phosphatase minus) + ATP, 37 °C 1 h",
                   MaterialKind.CATALYST)
    ready = g.product("ready", "5′-P-linker–RNA-3′-OH", MaterialKind.MOLECULAR, state=linked.state)
    g.react("phosphorylate_linker", "Phosphorylate the linker 5′ end",
            StepKind.PHOSPHORYLATION, ((oh3, InputRole.TEMPLATE), (pnk2, InputRole.CATALYST)),
            (ready,))
    rnl1_dilute = g.input("rnl1_dilute", "T4 RNA ligase 1 in 15 mL with DMSO, 16 °C overnight",
                          MaterialKind.CATALYST)
    chimera = g.product("chimera", "5′-RNA1–linker–RNA2-3′ chimera", MaterialKind.MOLECULAR,
                        state=molecule("MARIO chimera", (role("RNA1 fragment"),
                                                         role("biotinylated linker"),
                                                         role("RNA2 fragment")),
                                       junctions=(("RNA1 fragment", "biotinylated linker",
                                                   "proximity ligation"),
                                                  ("biotinylated linker", "RNA2 fragment",
                                                   "linker ligation"))))
    g.react("proximity_ligate", "Proximity-ligate a co-bound RNA 3′-OH to the linker 5′-P",
            StepKind.PROXIMITY_LIGATION,
            ((ready, InputRole.TEMPLATE), (rnl1_dilute, InputRole.CATALYST)), (chimera,))
    elute = g.input("elution", "biotin/SDS elution at 95 °C, acid phenol–chloroform, "
                    "precipitation", MaterialKind.SELECTION)
    protk = g.input("proteinase", "Proteinase K, 55 °C 2 h", MaterialKind.CATALYST)
    rna_mix = g.product("rna_mix", "purified RNA: chimeras plus unligated linker–RNA2 and "
                        "linker-free RNA", MaterialKind.MOLECULAR, state=chimera.state)
    protein_waste = g.product("protein_waste", "beads, digested protein and organic phase",
                              MaterialKind.WASTE)
    g.react("extract", "Elute, deproteinize and purify RNA", StepKind.DEPROTEINIZATION,
            ((chimera, InputRole.TEMPLATE), (elute, InputRole.SELECTION),
             (protk, InputRole.CATALYST)), (rna_mix, protein_waste))
    g.discard(protein_waste, "beads, protein digest and organic phase discarded")
    guard = g.input("linker_cdna_oligo", "T*C*G*C*-protected linker-complement DNA oligo",
                    MaterialKind.OLIGO, state=mc.LINKER_CDNA)
    t7 = g.input("t7_exo", "T7 exonuclease, 37 °C 1 h; then rigorous TURBO DNase",
                 MaterialKind.CATALYST)
    trimmed = g.product("trimmed", "RNA with terminal (5′) linkers degraded; internal "
                        "linkers of chimeras retained", MaterialKind.MOLECULAR,
                        state=chimera.state)
    g.react("remove_terminal_linker", "Degrade biotin-bearing terminal linkers on "
            "linker–RNA2 with the hybridized guard oligo and T7 exonuclease",
            StepKind.DIGESTION,
            ((rna_mix, InputRole.TEMPLATE), (guard, InputRole.STRUCTURAL),
             (t7, InputRole.CATALYST)), (trimmed,), changes_molecule=False)
    rnase_iii = g.input("rnase_iii", "RNase III fragmentation kit (150–400 nt); 2.2× SPRISelect",
                        MaterialKind.CATALYST)
    sheared = g.product("sheared", "sheared chimera fragments", MaterialKind.MOLECULAR,
                        state=chimera.state)
    g.react("shear", "Shear RNA with RNase III", StepKind.FRAGMENTATION,
            ((trimmed, InputRole.TEMPLATE), (rnase_iii, InputRole.CATALYST)), (sheared,))
    adapter = g.input("rt_adapter", "/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/ RT adapter",
                      MaterialKind.ADAPTER, state=mc.RT_ADAPTER)
    rnl2 = g.input("rnl2", "T4 RNA ligase 2 truncated KQ, PEG8000, 16 °C overnight",
                   MaterialKind.CATALYST)
    adapted = g.product("adapted", "RNA with 3′ RT adapter", MaterialKind.MOLECULAR,
                        state=molecule("adapter-ligated chimera RNA",
                                       (role("RNA1 fragment"), role("linker"),
                                        role("RNA2 fragment"), role("RT adapter")),
                                       junctions=(("RNA2 fragment", "RT adapter",
                                                   "pre-adenylated 3′ ligation"),)))
    g.react("ligate_rt_adapter", "Ligate the pre-adenylated RT adapter to RNA 3′ ends",
            StepKind.LIGATION,
            ((sheared, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
             (rnl2, InputRole.CATALYST)), (adapted,),
            geometry=AdapterLigationGeometry.unavailable(
                "Single-stranded 3′ ligation of the 5′-adenylated, 3′-ddC adapter to variable "
                "RNA 3′-OH ends; no duplex end geometry applies."))
    rt_primer = g.input("rt_primer", "5′-phosphorylated ES-1 RT primer", MaterialKind.OLIGO,
                        state="NN" + ES1_INDEX + "NNN" + mc.RT_BODY)
    ssiii = g.input("ssiii", "SuperScript III with T4 gp32, 50 °C 45 min, 55 °C 20 min; "
                    "no heat inactivation", MaterialKind.CATALYST)
    hybrid = g.product("hybrid", "RNA–cDNA hybrid", MaterialKind.MOLECULAR,
                       state=_state("RNA–cDNA hybrid", "RNA (with biotin linker)",
                                    "cDNA from 5′-P RT primer"))
    g.react("reverse_transcribe", "Reverse-transcribe from the adapter-annealed RT primer",
            StepKind.REVERSE_TRANSCRIPTION,
            ((adapted, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (ssiii, InputRole.CATALYST)), (hybrid,))
    c1 = g.input("c1_beads", "MyOne C1 streptavidin; B&W washes", MaterialKind.SELECTION)
    captured = g.product("captured", "bead-bound biotin-linker RNA–cDNA hybrids",
                         MaterialKind.MOLECULAR, state=hybrid.state)
    no_linker = g.product("no_linker", "hybrids without a biotinylated linker",
                          MaterialKind.WASTE)
    g.react("capture", "Capture linker-containing hybrids on streptavidin",
            StepKind.AFFINITY_CAPTURE,
            ((hybrid, InputRole.TEMPLATE), (c1, InputRole.SELECTION)), (captured, no_linker))
    g.discard(no_linker, "unbound fraction and B&W washes discarded")
    rnase_h = g.input("rnase_h", "RNase H, 37 °C 1 h; 70 °C 20 min; 2.2× SPRISelect",
                      MaterialKind.CATALYST)
    cdna = g.product("cdna", "released single-stranded cDNA", MaterialKind.MOLECULAR,
                     state=_state("ss cDNA", "5′-P RT primer", "antisense chimera cDNA"))
    g.react("release", "Release cDNA by digesting the RNA strand", StepKind.DIGESTION,
            ((captured, InputRole.TEMPLATE), (rnase_h, InputRole.CATALYST)), (cdna,))
    circligase = g.input("circligase", "CircLigase II, betaine/MnCl₂, 60 °C 2 h",
                         MaterialKind.CATALYST)
    circle = g.product("circle", "circular ss cDNA", MaterialKind.MOLECULAR,
                       state=closed_circle("circular cDNA", (role("RT primer"),
                                                             role("antisense chimera cDNA")),
                                           closure_label="CircLigase cDNA 3′ to RT-primer 5′-P"))
    g.react("circularize", "Circularize cDNA (truncated cDNAs included)",
            StepKind.CIRCULARIZATION,
            ((cdna, InputRole.TEMPLATE), (circligase, InputRole.CATALYST)), (circle,),
            geometry=CircularizationGeometry.known(
                "single-strand intramolecular ligation",
                statement="Methods 'Circularization': CircLigase II joins the cDNA 3′ end to the "
                          "5′-phosphorylated RT-primer end (Supplementary Fig. 1)."))
    cut_oligo = g.input("cut_oligo", "Cut_oligo (3′ inverted dT)", MaterialKind.OLIGO,
                        state=mc.CUT_OLIGO)
    bamhi = g.input("bamhi", "FastDigest BamHI, 37 °C 30 min; 2× SPRISelect",
                    MaterialKind.CATALYST)
    linear = g.product("linear", "relinearized cDNA", MaterialKind.MOLECULAR,
                       state=molecule("relinearized ES-1 cDNA", tuple(LINEAR_CDNA)))
    g.react("relinearize", "Anneal Cut_oligo and cut the RT-primer BamHI site",
            StepKind.DIGESTION,
            ((circle, InputRole.TEMPLATE), (cut_oligo, InputRole.STRUCTURAL),
             (bamhi, InputRole.CATALYST)), (linear,))
    dp5 = g.input("dp5", "DP5 truncated forward primer", MaterialKind.OLIGO, state=mc.DP5)
    dp3 = g.input("dp3", "DP3 truncated reverse primer", MaterialKind.OLIGO, state=mc.DP3)
    hifi = g.input("nebnext_hifi", "NEBNext High-Fidelity 2× PCR Master Mix",
                   MaterialKind.CATALYST)
    pre_amp_state = Construct([
        Segment("DP3 site", mc.DP3),
        Segment("antisense chimera cDNA", "X" * 20, placeholder=True),
        *list(LINEAR_CDNA)[2:],
    ], name="DP5/DP3 pre-amplicon")
    pre = g.product("preamp", "DP5/DP3 pre-amplified cDNA (6 cycles; 150–350 bp gel cut)",
                    MaterialKind.MOLECULAR, state=molecule("DP pre-amplicon",
                                                           tuple(pre_amp_state)))
    g.react("preamplify", "Pre-amplify with DP5/DP3", StepKind.PCR,
            ((linear, InputRole.TEMPLATE), (dp5, InputRole.STRUCTURAL),
             (dp3, InputRole.STRUCTURAL), (hifi, InputRole.CATALYST)), (pre,),
            geometry=PCRGeometry.exact(
                LINEAR_CDNA, PCRPrimer("DP3", mc.DP3), PCRPrimer("DP5", mc.DP5),
                statement="DP3 equals the RT-primer 3′ part retained after BamHI; DP5 is "
                          "complementary to the RT-primer 5′ arm (Methods oligo list)."))
    pe1 = g.input("pe1", "Illumina PE PCR Forward Primer 1.0", MaterialKind.OLIGO,
                  state=mc.PE1)
    pe2 = g.input("pe2", "Illumina PE PCR Reverse Primer 2.0", MaterialKind.OLIGO,
                  state=mc.PE2)
    library = g.product("library", "MARIO sequencing library (250–550 bp)",
                        MaterialKind.LIBRARY, state=final)
    g.react("final_pcr", "Amplify with full PE 1.0/2.0 primers (titrated cycles)",
            StepKind.PCR,
            ((pre, InputRole.TEMPLATE), (pe1, InputRole.STRUCTURAL),
             (pe2, InputRole.STRUCTURAL), (hifi, InputRole.CATALYST)), (library,),
            geometry=PCRGeometry.exact(
                pre_amp_state, PCRPrimer("PE 2.0", mc.PE2, anneal_nt=len(mc.DP3)),
                PCRPrimer("PE 1.0", mc.PE1, anneal_nt=len(mc.DP5)),
                statement="PE 1.0 and PE 2.0 end in DP5 and DP3 (Methods oligo list)."))
    g.sequence(library, SequencingPlan(library=final, primers=tuple(primers)))
    g.annotate(SourceDiscrepancy(
        "rt_primer", "RT-primer random-base count differs",
        (SourceClaim("Methods 'Reverse transcription' (PMC4931010 JATS line 175)",
                     "generic /5Phos/NNXXXXNNNN… with a 6-nt random barcode"),
         SourceClaim("Supplementary Fig. 1 legend", "10-nt barcode"),
         SourceClaim("Methods 'Oligonucleotide sequences used in MARIO' (line 229)",
                     "ES-1 primer /5Phos/NNAGGTNNN… (5 random + 4 sample bases)")),
        "The model draws the explicitly printed ES-1 oligo (NN + AGGT + NNN).",
        "ES-2/MEF and ES-indirect primers in the same list have six random bases, matching "
        "the generic form.",
        "Whether the ES-1 primer was ordered with five or six random bases is not resolved."))
    g.annotate(SourceDiscrepancy(
        "dephosphorylate", "3′-dephosphorylation enzyme is a phosphatase-minus kinase",
        (SourceClaim("Methods 'Proximity ligation' (JATS line 149)",
                     "RNA 3′ ends dephosphorylated using the 3′-phosphatase activity of T4 PNK, "
                     "but the mix lists T4 PNK (3′-phosphatase minus)"),
         SourceClaim("Methods 'Ligation of a biotin-tagged RNA linker' (line 143)",
                     "the 3′-phosphatase-minus PNK was chosen because it does not convert "
                     "3′-P to 3′-OH")),
        "The graph keeps the stated 3′-OH outcome needed for proximity ligation and labels "
        "the enzyme as printed.",
        "Both statements are in the published Methods; the printed mix is a pH 6.5 PNK "
        "buffer without ATP (line 149).",
        "The enzyme actually used for 3′ dephosphorylation is unresolved."))
    g.declare_reviewed_scope(tuple(g.materials) + tuple(g.steps), source)
    return g.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)


SPEC.protocol_graph = protocol_graph
