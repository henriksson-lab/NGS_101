"""Source-owned protocol definition and reaction graph.

Authority: Nguyen et al., Nat. Commun. 7, 12023 (2016), doi:10.1038/ncomms12023,
PMC4931010 JATS Methods (lines 116–261: cross-linking, lysis/RNase I, IPB biotinylation and
T1 immobilization, linker ligation, proximity ligation, T7-exonuclease terminal-linker
removal, RNase III shearing, RT-adapter ligation, RT, C1 pull-down/RNase H elution,
circularization, BamHI relinearization, DP5/DP3 pre-amplification, final PE PCR and the
oligo list) and Supplementary Fig. 1 (circularization strategy).  The ES-1 sample is drawn;
the ES-2/MEF and ES-indirect RT primers are retained as source inventory. Their
GeneRead and dual-cross-link branches are separate; the drawn ES-1 route includes DSN.
"""
from __future__ import annotations

from types import SimpleNamespace
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, revcomp
from circular import CircularizationGeometry, circularize_ssdna
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind, Disposition, Terminal, MolecularRelation)
from protocol_model import molecule, role

import mario_chem as mc

TITLE = "MARIO — mapping RNA interactome in vivo"
NOTES = "01_mario.html"
SOURCE = '<a href="https://doi.org/10.1038/ncomms12023">Nguyen et al., Nature Communications (2016)</a>, Methods and Supplementary Figure1.'
SUMMARY = "Protein-bound RNAs receive an internally biotinylated RNA linker before dilute proximity ligation. Linker-bearing RNA/cDNA hybrids are captured after RT; cDNA circularization, BamHI reopening and PCR yield the library."
SEQUENCING_INTRO = "Published paired-end HiSeq sequencing; biological read lengths remain unassigned. The explicitly printed ES-1 RT primer supplies five random bases plus a four-base sample index, conflicting with the generic six-random-base description. The inline identifier is reverse-complemented in Read1."
SEQUENCING_ENDING = SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE = "Only the compatible Read1 primer is placed. The source-positive Read2 is unplaced; its exact primer and both read lengths remain unassigned."
READ_LENGTHS = {}
ACQUIRED_ROLES = ("Read 1", "Read 2")
READ_USES = {"Read 1": "extract_then_include", "Read 2": "include"}
SHOW_READ_CONTENT = False
READ_SPECS = (
    sp.read_spec("Read 1", "unassigned", "ES-1 copied inline identifier (3 random + ACCT + 2 random), then chimera-derived sequence", "extract_then_include"),
    sp.read_spec("Read 2", "unassigned", "Opposite chimera-derived end; source-positive paired read, primer geometry unassigned", "include"),
)
READ_SPECS_NOTE = "The paper declares paired biological reads; no external index acquisition recipe is inferred from the inline sample barcode."
FINAL_CAPTION = "Source-primer-derived ES-1 target amplicon; variable RNA-fragment lengths and actual DSN/PCR yield remain unassigned."

CAVEAT = (
    "<b>Drawn sample:</b> ES-1 (RT primer <code>/5Phos/NNAGGTNNN…</code>). Exact sequences come "
    "from the paper's oligo list; RNA fragments are role tokens. The exact Read2 sequencing primer "
    "is unassigned in this bounded model, so only the compatible Read1 site is placed. The "
    "linker-to-read transfer of the internal biotin-dT as T is derived. rRNA depletion "
    "(DSN selected here for ES-1; GeneRead applies to ES-2/MEF), the dual-cross-linked ES-indirect "
    "branch are separate source variants. The selected ES-1 pathway includes DSN rRNA-cDNA depletion. "
    "The two PCR steps use their actual ss or post-DSN precursors; exact first-cycle duplex geometry is not asserted.")

LINKER_DNA = (mc.RNA_LINKER_LEFT + "T" + mc.RNA_LINKER_RIGHT).replace("U", "T")
ES1_INDEX = "AGGT"
BAMHI_CUT = mc.RT_BODY.index("GGATCC") + 1  # G^GATCC; source RT-primer site
RT_ARM5 = mc.RT_BODY[:BAMHI_CUT]
RT_ARM3 = mc.RT_BODY[BAMHI_CUT:]
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
        role("antisense chimera cDNA; variable RT extent"),
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
        role("RNA1 fragment"),
        Segment("MARIO linker copy", LINKER_DNA, "linker"),
        role("RNA2 fragment"),
        Segment("RT adapter copy", right[:len(mc.RT_ADAPTER)]),
        Segment("PE 2.0 read-2 arm", right[len(mc.RT_ADAPTER):-len(il.P7)], "s7"),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name=TITLE + " final library (ES-1)")
    primers = (sp.TRUSEQ["R1"],)
    problems = sp.verify(lib, primers, required_roles=("Read 1",))
    if problems:
        raise ValueError("; ".join(problems))
    return lib, primers


LINEAR_CDNA = _linear_cdna()
FINAL_LIBRARY, SEQ_PRIMERS = _final()


def rt_primer_parts():
    return (_umi("part 1", 2), Segment("sample barcode", ES1_INDEX, "cbc", feature=_INDEX), _umi("part 2", 3), Segment("RT-primer body including BamHI", mc.RT_BODY))


def rnp_state(*, mod5="OH", mod3="P", linked=False):
    sc = Scene()
    for name in ("RNA1", "RNA2"):
        parts = list(linker_parts()) + [role(name + " fragment")] if linked else [role(name + " fragment")]
        sc.strand(name, parts, label="independent RNA in one retained protein complex", mod5=mod5, mod3=mod3)
        sc.note(name, "protein contact retains the RNAs; it is not a covalent RNA–protein–RNA strand")
    return sc


def linker_parts():
    return (Segment("RNA linker left arm", mc.RNA_LINKER_LEFT), Segment("internal biotin-dT", "T", note="source /iBiodT/"), Segment("RNA linker right arm", mc.RNA_LINKER_RIGHT))


def chimera_parts():
    return (role("RNA1 fragment; upstream terminal occupancy unassigned"), *linker_parts(), role("RNA2 fragment"))


def rt_adapter_parts():
    return (Segment("RT adapter core", mc.RT_ADAPTER, note="source /5rApp/ precursor"),
            Segment("terminal ddC", "C", note="source /3ddC/; blocked3′end"))


def hybrid_state():
    sc = Scene()
    sc.strand("RNA", (*chimera_parts(), *rt_adapter_parts()), label="RNA template with linker biotin and retained terminal ddC; paired extent varies", mod3="ddC")
    sc.strand("cDNA", (*rt_primer_parts(), role("antisense chimera cDNA; variable RT extent")), label="independent complementary cDNA; paired extent unassigned", mod5="P")
    return sc


def cdna_construct():
    return Construct((*rt_primer_parts(), role("antisense chimera cDNA; variable RT extent")), name="actual ES-1 first-strand cDNA")


def protocol_graph():
    source = SOURCE
    final, primers = FINAL_LIBRARY, SEQ_PRIMERS
    g = ProtocolGraph(TITLE, authority=source)
    rnps = g.input("rnps", "E14 mES ES-1 cells with protein-bound RNAs",
                   MaterialKind.BIOLOGICAL, state=rnp_state(mod5="unassigned", mod3="unassigned"),
                   primary=True)
    uv = g.input("uv", "UV-C 254 nm, 400 mJ cm⁻² on ice", MaterialKind.CONDITION)
    fixed = g.product("fixed", "UV-cross-linked RNPs", MaterialKind.MOLECULAR, state=rnps.state)
    g.react("crosslink", "UV-cross-link RNA to bound proteins", StepKind.CROSS_LINKING,
            ((rnps, InputRole.TEMPLATE), (uv, InputRole.CONDITION)), (fixed,))
    rnase_i = g.input("rnase_i", "TURBO DNase, then RNase I 1:100, 37 °C 3 min",
                      MaterialKind.CATALYST)
    frag = g.product("fragmented", "lysate RNA fragments (5′-OH, 3′-P) on cross-linked proteins",
                     MaterialKind.MOLECULAR,
                     state=rnp_state())
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
                      state=rnp_state(mod5="P"))
    g.react("phosphorylate_rna", "Phosphorylate RNA 5′ ends without removing 3′ phosphates",
            StepKind.PHOSPHORYLATION,
            ((bound, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (p5rna,))
    linker = g.input("linker", "whole internally biotinylated RNA linker (5′-OH), 30 µL 20 µM",
                     MaterialKind.ADAPTER, state=molecule("source RNA linker", linker_parts()))
    rnl1 = g.input("rnl1", "T4 RNA ligase 1, PEG8000, 16 °C overnight", MaterialKind.CATALYST)
    linked = g.product("linker_rna", "5′-HO-linker–RNA-3′-P", MaterialKind.MOLECULAR,
                       state=rnp_state(linked=True))
    g.react("ligate_linker", "Ligate the linker 3′ end to RNA 5′ ends", StepKind.LIGATION,
            ((p5rna, InputRole.TEMPLATE), (linker, InputRole.STRUCTURAL),
             (rnl1, InputRole.CATALYST)), (linked,),
            geometry=AdapterLigationGeometry.unavailable(
                "Single-stranded RNA ligation: the linker's 3′-OH joins variable RNA 5′-P ends; "
                "the 5′-OH linker cannot self-ligate. No duplex end geometry applies."))
    pnk_buf = g.input("pnk_ph65", "'T4 PNK (3′-phosphatase minus)' in pH 6.5 PNK buffer without "
                      "ATP, 37 °C 20 min (as printed)", MaterialKind.CATALYST)
    oh3 = g.product("oh3", "linker–RNA with 3′-OH (intended)", MaterialKind.MOLECULAR,
                    state=rnp_state(linked=True, mod3="OH"), inferred=True)
    g.react("dephosphorylate", "INFERRED intended outcome — remove RNA 3′ phosphates", StepKind.DEPHOSPHORYLATION,
            ((linked, InputRole.TEMPLATE), (pnk_buf, InputRole.CATALYST)), (oh3,))
    pnk2 = g.input("pnk_linker", "T4 PNK (3′-phosphatase minus) + ATP, 37 °C 1 h",
                   MaterialKind.CATALYST)
    ready = g.product("ready", "5′-P-linker–RNA-3′-OH", MaterialKind.MOLECULAR, state=rnp_state(linked=True, mod5="P", mod3="OH"), inferred=True)
    g.react("phosphorylate_linker", "Phosphorylate the linker 5′ end",
            StepKind.PHOSPHORYLATION, ((oh3, InputRole.TEMPLATE), (pnk2, InputRole.CATALYST)),
            (ready,))
    rnl1_dilute = g.input("rnl1_dilute", "T4 RNA ligase 1 in 15 mL with DMSO, 16 °C overnight",
                          MaterialKind.CATALYST)
    chimera = g.product("chimera", "5′-RNA1–linker–RNA2-3′ chimera", MaterialKind.MOLECULAR,
                        state=molecule("internal target chimera region; terminal classes vary", chimera_parts()))
    g.react("proximity_ligate", "Proximity-ligate a co-bound RNA 3′-OH to the linker 5′-P",
            StepKind.PROXIMITY_LIGATION,
            ((ready, InputRole.TEMPLATE), (rnl1_dilute, InputRole.CATALYST)), (chimera,))
    elute = g.input("elution", "biotin/SDS elution at 95 °C, acid phenol–chloroform, "
                    "precipitation", MaterialKind.SELECTION)
    protk = g.input("proteinase", "Proteinase K, 55 °C 2 h", MaterialKind.CATALYST)
    rna_mix = g.product("rna_mix", "purified RNA: chimeras plus unligated linker–RNA2 and "
                        "linker-free RNA", MaterialKind.MOLECULAR, state=chimera.state)
    protein_waste = g.product("protein_waste", "physically separated beads and non-retained extraction fractions",
                              MaterialKind.WASTE)
    protein_digest = g.product("protein_digest", "protein-derived digestion products; exact extraction partition unassigned", MaterialKind.MOLECULAR)
    g.react("extract", "Elute, deproteinize and purify RNA", StepKind.DEPROTEINIZATION,
            ((chimera, InputRole.TEMPLATE), (elute, InputRole.SELECTION),
             (protk, InputRole.CATALYST)), (rna_mix, protein_waste, protein_digest))
    g.discard(protein_waste, "physically separated beads and extraction fractions; this is not a quantitative assignment of every peptide")
    g.mark_terminal(Terminal(protein_digest.id, Disposition.UNRESOLVED,
        "Proteolysis does not establish physical disposal; exact peptide partition through RNA extraction/precipitation remains unassigned."))
    guard = g.input("linker_cdna_oligo", "T*C*G*C*-protected linker-complement DNA oligo",
                    MaterialKind.OLIGO, state=mc.LINKER_CDNA)
    t7 = g.input("t7_exo", "T7 exonuclease, 37 °C 1 h; then rigorous TURBO DNase",
                 MaterialKind.CATALYST)
    trimmed = g.product("trimmed", "RNA with terminal (5′) linkers degraded; internal "
                        "linkers of chimeras retained", MaterialKind.MOLECULAR,
                        state=chimera.state)
    degraded_linkers=g.product("terminal_digest_products", "terminal-linker RNA and guard-DNA degradation products", MaterialKind.MOLECULAR)
    g.react("remove_terminal_linker", "Degrade biotin-bearing terminal linkers on "
            "linker–RNA2 with the hybridized guard oligo and T7 exonuclease",
            StepKind.DIGESTION,
            ((rna_mix, InputRole.TEMPLATE), (guard, InputRole.STRUCTURAL),
             (t7, InputRole.CATALYST)), (trimmed, degraded_linkers))
    g.mark_terminal(Terminal(degraded_linkers.id, Disposition.UNRESOLVED, "Exonuclease/DNase chemistry does not establish physical disposal; partition through the following extraction/precipitation remains unassigned."))
    rnase_iii = g.input("rnase_iii", "RNase III fragmentation kit (150–400 nt); 2.2× SPRISelect",
                        MaterialKind.CATALYST)
    sheared = g.product("sheared", "sheared chimera fragments", MaterialKind.MOLECULAR,
                        state=molecule("target junction-bearing sheared RNA", chimera_parts()))
    g.react("shear", "Shear RNA with RNase III", StepKind.FRAGMENTATION,
            ((trimmed, InputRole.TEMPLATE), (rnase_iii, InputRole.CATALYST)), (sheared,))
    adapter = g.input("rt_adapter", "/5rApp/AGATCGGAAGAGCGGTTCAG/3ddC/ RT adapter",
                      MaterialKind.ADAPTER, state=Construct(rt_adapter_parts()))
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
                        state=molecule("5′phosphorylated ES-1 RT primer", rt_primer_parts()))
    ssiii = g.input("ssiii", "SuperScript III with T4 gp32, 50 °C 45 min, 55 °C 20 min; "
                    "no heat inactivation", MaterialKind.CATALYST)
    hybrid = g.product("hybrid", "RNA–cDNA hybrid", MaterialKind.MOLECULAR,
                       state=hybrid_state())
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
    rnase_h = g.input("rnase_h", "RNase H, 37 °C 1 h; collect supernatant,70 °C 20 min;2.2×SPRISelect",
                      MaterialKind.CATALYST)
    cdna_chemical = cdna_construct()
    cdna = g.product("cdna", "released single-stranded cDNA", MaterialKind.MOLECULAR,
                     state=molecule("5′phosphorylated first-strand cDNA", tuple(cdna_chemical)), chemical_state=cdna_chemical)
    rna_digest = g.product("rnase_h_products", "RNA hydrolysates after RNase H; biotin-bearing fraction remains with beads", MaterialKind.MOLECULAR)
    g.react("release", "Release cDNA by digesting the RNA strand", StepKind.DIGESTION,
            ((captured, InputRole.TEMPLATE), (rnase_h, InputRole.CATALYST)), (cdna, rna_digest))
    g.mark_terminal(Terminal(rna_digest.id, Disposition.UNRESOLVED, "Source collects cDNA-containing supernatant; exact hydrolysate partition between beads, supernatant and SPRI cleanup is unassigned."))
    circligase = g.input("circligase", "CircLigase II, betaine/MnCl₂, 60 °C 2 h",
                         MaterialKind.CATALYST)
    circle_chemical = circularize_ssdna(cdna_chemical, five_prime_phosphate=True)
    circle = g.product("circle", "circular ss cDNA", MaterialKind.MOLECULAR,
                       state=circle_chemical, chemical_state=circle_chemical)
    g.react("circularize", "Circularize cDNA (truncated cDNAs included)",
            StepKind.CIRCULARIZATION,
            ((cdna, InputRole.TEMPLATE), (circligase, InputRole.CATALYST)), (circle,),
            relations=(MolecularRelation.CIRCLE_CLOSURE,),
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
        role("antisense chimera cDNA; variable RT extent"),
        *list(LINEAR_CDNA)[2:-1],
        Segment("DP5-derived opposite arm", revcomp(mc.DP5)),
    ], name="DP5/DP3 pre-amplicon")
    pre = g.product("preamp", "ES-1 DP5/DP3 pre-amplified cDNA; cycles titrated to80–100ng, no pre-DSN gel",
                    MaterialKind.MOLECULAR, state=pre_amp_state)
    g.react("preamplify", "ES-1: pre-amplify with DP5/DP3 for DSN depletion", StepKind.PCR,
            ((linear, InputRole.TEMPLATE), (dp5, InputRole.STRUCTURAL),
             (dp3, InputRole.STRUCTURAL), (hifi, InputRole.CATALYST)), (pre,),
            geometry=PCRGeometry.unavailable(linear.state,
                reason="Actual BamHI-reopened first-strand cDNA is single stranded. Exact primer sequences are known, but a duplex with both initial primer sites is not supplied."))
    dsn_prep=g.input("dsn_preparation", "1.8×SPRI; ES-1 80–100ng, skip150–350gel;98°C2min then68°C5h reannealing", MaterialKind.CONDITION)
    dsn=g.input("dsn", "DSN1U/µl2µl;68°C25min, stopsolution and1.8×SPRI", MaterialKind.CATALYST)
    surviving=g.product("dsn_survivors", "ES-1 surviving cDNA population after DSN; native ss/ds ends unassigned", MaterialKind.MOLECULAR,
        state=molecule("actual post-DSN candidates", (role("remaining cDNA species; duplex extent and surviving PCR ends unassigned"),)))
    hydrolysed=g.product("dsn_hydrolysates", "DSN-derived DNA products; cleanup partition unassigned", MaterialKind.MOLECULAR)
    g.react("dsn_deplete", "ES-1: reanneal and digest rapidly forming abundant duplexes with DSN", StepKind.DIGESTION,
        ((pre,InputRole.TEMPLATE),(dsn_prep,InputRole.CONDITION),(dsn,InputRole.CATALYST)),(surviving,hydrolysed))
    g.mark_terminal(Terminal(hydrolysed.id, Disposition.UNRESOLVED, "DSN hydrolysis does not establish quantitative physical exclusion through the source1.8×SPRI cleanup."))
    pe1 = g.input("pe1", "Illumina PE PCR Forward Primer 1.0", MaterialKind.OLIGO,
                  state=mc.PE1)
    pe2 = g.input("pe2", "Illumina PE PCR Reverse Primer 2.0", MaterialKind.OLIGO,
                  state=mc.PE2)
    library = g.product("library", "MARIO PE-primer PCR products before size selection",
                        MaterialKind.MOLECULAR, state=final)
    g.react("final_pcr", "Amplify with full PE 1.0/2.0 primers (titrated cycles)",
            StepKind.PCR,
            ((surviving, InputRole.TEMPLATE), (pe1, InputRole.STRUCTURAL),
             (pe2, InputRole.STRUCTURAL), (hifi, InputRole.CATALYST)), (library,),
            geometry=PCRGeometry.unavailable(surviving.state,
                reason="Actual post-DSN DNA population; source primers are explicit but native surviving strand pairing/ends are unassigned, not a reconstructed pre-DSN duplex."))
    cleanup=g.input("final_selection", "1.8×SPRI and gel selection250–550bp", MaterialKind.SELECTION)
    selected=g.product("selected_library", "size-selected ES-1 MARIO library", MaterialKind.LIBRARY,state=final)
    excluded=g.product("final_excluded", "SPRI supernatants/washes and non-excised gel material", MaterialKind.WASTE)
    g.react("final_select", "Clean up and select250–550bp final products", StepKind.SELECTION,
        ((library,InputRole.TEMPLATE),(cleanup,InputRole.SELECTION)),(selected,excluded),changes_molecule=False)
    g.discard(excluded, "Source physical supernatant/wash removals and gel excision; no molecular hydrolysate is equated with this entire discarded fraction.")
    g.sequence(selected, SequencingPlan(library=final, primers=tuple(primers),acquired_roles=("Read 1",),
        unplaced_reads=(sp.read_spec("Read 2", "unassigned", "Source-positive paired biological read; primer unassigned", "include"),),limitation=SEQUENCING_UNAVAILABLE))
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
    g.declare_reviewed_scope(("rnps", "uv", "crosslink", "rnase_i", "fragment", "ipb", "t1_beads", "immobilize", "pnk_minus", "phosphorylate_rna", "linker", "rnl1", "ligate_linker", "pnk_ph65", "dephosphorylate", "pnk_linker", "phosphorylate_linker", "rnl1_dilute", "proximity_ligate", "elution", "proteinase", "extract", "linker_cdna_oligo", "t7_exo", "remove_terminal_linker", "rnase_iii", "shear", "rt_adapter", "rnl2", "ligate_rt_adapter", "rt_primer", "ssiii", "reverse_transcribe", "c1_beads", "capture", "rnase_h", "release", "circligase", "circularize", "cut_oligo", "bamhi", "relinearize", "dp5", "dp3", "nebnext_hifi", "preamplify", "dsn_preparation", "dsn", "dsn_deplete", "pe1", "pe2", "final_pcr", "final_selection", "final_select", "selected_library"), source)
    return g.require_publishable()


SPEC = SimpleNamespace(TITLE=TITLE, NOTES=NOTES, SOURCE=SOURCE, SUMMARY=SUMMARY, CAVEAT=CAVEAT,
    FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=SEQ_PRIMERS, READ_LENGTHS=READ_LENGTHS,
    ACQUIRED_ROLES=ACQUIRED_ROLES, READ_USES=READ_USES, SHOW_READ_CONTENT=SHOW_READ_CONTENT,
    READ_SPECS=READ_SPECS, READ_SPECS_NOTE=READ_SPECS_NOTE, FINAL_CAPTION=FINAL_CAPTION,
    SEQUENCING_INTRO=SEQUENCING_INTRO, SEQUENCING_ENDING=SEQUENCING_ENDING,
    SEQUENCING_UNAVAILABLE=SEQUENCING_UNAVAILABLE, oligos=mc.oligos, protocol_graph=protocol_graph)
