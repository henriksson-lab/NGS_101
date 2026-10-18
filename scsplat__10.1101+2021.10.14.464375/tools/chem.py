"""Source-owned scSPLAT protocol definition and reaction graph.

Authority: Raine et al., Scientific Reports 12:5772 (2022), doi:10.1038/s41598-022-09798-2
(PMC8986790), Methods "scWGBS library preparation and sequencing" and Supplementary Table S1.
The bioRxiv preprint 10.1101/2021.10.14.464375 names the directory but is superseded.
"""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, oligo
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind, Terminal,
                            Disposition)
from protocol_model import role

PUBLISHED = "Raine et al., Scientific Reports (2022), doi:10.1038/s41598-022-09798-2"
NEBNEXT = (PUBLISHED + " names NEBNext dual index primers; native full outer arms, "
           "selected members and physical widths remain unassigned in this pass")

# Supplementary Table S1 (MOESM1): SSO_H_1..SSO_H_32 share this 5′ tail, then a
# six-base cell barcode and eight H (A/C/T) random bases. SSO_H_1 is drawn as the
# representative family member; every member carries 5AmMC6.
SSO_TAIL = "TTCCCTACACGACGCTCTTCCGATCT"
SSO_H_1_BARCODE = "TAGCTT"
H8 = "HHHHHHHH"
# Published Methods / Table S1 SPLAT-adapter oligos.
SS1A = "GACGTGTGCTCTTCCGATCTNNNNNN"     # 5′AmMC6 … 3′AmMo
SS1B = "AGATCGGAAGAGCACACGTC"           # 5′Phos … 3′AmMo

if il.TRUSEQ_READ1[-len(SSO_TAIL):] != SSO_TAIL:
    raise AssertionError("SSO_H tail is not the 3′ 26 nt of TruSeq Read 1")
_READ2_TOP = "A" + il.INDEX1_PRIMER        # top-strand Read 2 / Index 1 arm
if not _READ2_TOP.startswith(SS1B):
    raise AssertionError("ss1b is not the 5′ start of the TruSeq Read 2 arm")

CELL_BARCODE = feature("scsplat_cell_barcode", "cell_barcode", "whitelist",
                       whitelist="32 six-base SSO_H_X barcodes, Supplementary Table S1",
                       note="Methods: first six bases of Read 1; demultiplexed with je")

TITLE = 'scSPLAT — scalable single-cell WGBS'
NOTES = '01_scsplat.html'
SOURCE = '<a href="https://doi.org/10.1038/s41598-022-09798-2">Raine et al., Scientific Reports2022</a>, Methods and SupplementaryTableS1. Directory retains superseded preprint identity.'
SUMMARY = 'Converted single-cell DNA receives a barcoded random-primer copy before pooling and bulk splinted ligation of the opposite adapter.'
CAVEAT = 'INFERRED — native NEBNext dual-index primer members, physical index widths and complete outer-arm architecture remain unassigned in this capped reconstruction. Exact source SSO_H and blocked ss1a/ss1b components are retained without a default completed TruSeq shell. Physical converted U strands are distinct from copied DNA. Nonamplification and ExoI hydrolysis do not establish physical removal.'
INDEX5 = feature('sample_index_i5','sample_index','unknown',note='Native member/width unassigned')
INDEX7 = feature('sample_index_i7','sample_index','unknown',note='Native member/width unassigned')
FINAL_LIBRARY = Construct([
    role('native NEBNext forward end; full arm/member/width unassigned',inferred=True,molecular_feature=INDEX5),
    Segment('SSO_H P5 tail',SSO_TAIL,'r1'),
    Segment('cell barcode','N'*6,'cbc',placeholder=True,feature=CELL_BARCODE),
    Segment('H8 random-primer copy',H8,placeholder=True),
    role('copied genomic DNA reflecting bisulfite conversion'),
    Segment('ss1b SPLAT adapter',SS1B,'r2'),
    role('native NEBNext opposite end; full arm/member/width unassigned',inferred=True,molecular_feature=INDEX7),
],name='source-qualified scSPLAT endpoint')
SEQ_PRIMERS = ()
FINAL_CAPTION = 'INFERRED native outer ends — source SSO_H/barcode/H8/copy/ss1b components; no default index widths or completed canonical shell.'
SEQUENCING_UNAVAILABLE = 'Native sequencing-primer sites and outer adapter/index architecture remain unassigned. Source PE150 and inline six-base Read1 barcode do not establish physical index widths or acquired index cycles.'
SEQUENCING_ENDING = 'Published PE150 on HiSeqX or NovaSeq6000; Read1 begins with the six-base cell barcode. Native index acquisition and primer-site geometry remain unassigned.'
SHOW_READ_CONTENT = False
READ_SPECS = (sp.read_spec('Read 1',150,'six-base inline cell barcode followed by random-primer/copied-DNA sequence','extract_then_include'),sp.read_spec('Read 2',150,'opposite copied-DNA end','include'),sp.read_spec('Index 1 (i7)','unassigned','Native pool index; acquisition/cycles unassigned','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','Native pool index; acquisition/cycles unassigned','undetermined'))


def sso_h_1() -> Scene:
    sc = Scene()
    sc.strand("SSO_H_1", [Segment("SSO_H P5 tail", SSO_TAIL, "r1"),
                          Segment("SSO_H_1 cell barcode", SSO_H_1_BARCODE, "cbc",
                                  feature=CELL_BARCODE),
                          Segment("H8 random primer", H8, placeholder=True)],
              label="SSO_H_1 (representative of SSO_H_1–32)", mod5="AmMC6")
    return sc


def splat_adapter() -> Scene:
    sc = Scene()
    sc.strand("ss1b", [Segment("ss1b ligating strand", SS1B)], label="SPLAT adapter",
              mod5="p", mod3="AmMo")
    sc.anneal("ss1a", [Segment("ss1a paired arm", SS1A[:-6]),
                       Segment("ss1a random splint", "N" * 6)],
              to="ss1b", pair=("ss1a paired arm", "ss1b ligating strand"),
              label="SPLAT adapter", mod5="AmMC6", mod3="AmMo",
              unpaired=("ss1a random splint",))
    sc.overhang("ss1a", "ss1a random splint", polarity="3-prime", owner="adapter",
                text="3′ random hexamer splint")
    sc.ligation_end("ss1b", "ss1b ligating strand", kind="splinted", polarity="5-prime",
                    owner="adapter", text="5′ phosphate joins the single-strand 3′-OH")
    return sc


def adapter_geometry() -> AdapterLigationGeometry:
    return AdapterLigationGeometry.known(
        "splinted", splat_adapter(),
        statement="Published Methods/Table S1 print ss1a (5′AmMC6/GACGTGTGCTCTTCCGATCTNNNNNN/3′AmMo) "
                  "and ss1b (5′Phos/AGATCGGAAGAGCACACGTC/3′AmMo); the protruding 3′ N6 of the bottom "
                  "strand hybridizes to the single-stranded 3′ end so T4 DNA ligase joins it to the "
                  "ss1b 5′ phosphate. The Scene verifies that the printed sequences pair.")


def oligos():
    return [
        oligo("SSO_H_1 (Table S1; SSO_H_2–32 differ only in the six-base barcode)",
              (Segment("tail", SSO_TAIL, "r1"), Segment("barcode", SSO_H_1_BARCODE, "cbc", feature=CELL_BARCODE),
               Segment("H8", H8, placeholder=True)), mods="/5AmMC6/"),
        oligo("ss1a (SPLAT adapter bottom strand)", (Segment("ss1a", SS1A),),
              mods="/5AmMC6/", three="/3AmMo/-3′"),
        oligo("ss1b (SPLAT adapter top strand; Methods also call the partner ss2b)",
              (Segment("ss1b", SS1B),), mods="/5Phos/", three="/3AmMo/-3′"),
    ]


def _barcoded_segments() -> list[Segment]:
    return [Segment("SSO_H P5 tail", SSO_TAIL, "r1"),
            Segment("cell barcode", "X" * 6, "cbc", placeholder=True, feature=CELL_BARCODE),
            Segment("H8 random-primer copy", H8, placeholder=True),
            Segment("Klenow copy of converted DNA", "X" * 28, placeholder=True)]


def _single(label: str, segs: list[Segment]) -> Scene:
    sc = Scene()
    sc.strand("strand", segs, label=label, mod5="AmMC6" if label in ("barcoded second strand", "SPLAT-ligated single strand") else "", mod3="AmMo" if label=="SPLAT-ligated single strand" else "")
    return sc


def _primed_duplex() -> Scene:
    """Barcoded copy on its converted template; the SSO tail and barcode stay unpaired."""
    sc = _single("barcoded second strand", _barcoded_segments())
    sc.anneal("template", [Segment("converted template", "x" * 28, placeholder=True),
                           Segment("H8 binding site", "D" * 8, placeholder=True)],
              to="strand", pair=("converted template", "Klenow copy of converted DNA"),
              label="original converted strand")
    return sc


def _ligated_strand() -> Scene:
    sc = _single("SPLAT-ligated single strand",
                 [*_barcoded_segments(), Segment("ss1b SPLAT adapter", SS1B, "r2")])
    sc.junction("strand", "Klenow copy of converted DNA", "ss1b SPLAT adapter",
                "splinted ligation")
    return sc


def _ss(label: str) -> Scene:
    return _single(label, [Segment(label, "X" * 28, placeholder=True)])


def protocol_graph():
    src = PUBLISHED + ", Methods scWGBS library preparation and sequencing"
    graph = ProtocolGraph(TITLE, authority=src)
    cell = graph.input("cell", "FACS-sorted single cell or nucleus in 2 µl Zymo M-digestion buffer with proteinase K",
                       MaterialKind.BIOLOGICAL, state=Construct([Segment("single-cell genomic DNA", "X" * 28, placeholder=True)], name="genomic DNA"), primary=True)
    bisulfite = graph.input("bisulfite", "15 µl EZ Methylation Direct conversion reagent (manufacturer's program)",
                            MaterialKind.CONDITION)
    converted_raw = graph.product("converted_raw", "bisulfite-treated DNA in conversion reagent",
                                  MaterialKind.MOLECULAR, state=_ss("bisulfite-converted DNA: physical U at converted C, retained methylcytosines"))
    graph.react("convert", "Lyse and bisulfite-convert DNA in each 384-plate well",
                StepKind.CHEMICAL_CONVERSION,
                ((cell, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)), (converted_raw,))
    column = graph.input("spin_plate", "Zymo 384-well spin plate: binding, M-wash, M-desulphonation 15 min RT, two washes",
                         MaterialKind.SELECTION)
    sso = graph.input("sso_h", "SSO_H_X barcoded strand-synthesis oligo, 500 nM in 7 µl elution buffer",
                      MaterialKind.OLIGO, state=sso_h_1(),
                      source=PUBLISHED + ", Supplementary Table S1 (SSO_H_1 drawn)")
    converted = graph.product("converted", "desulphonated single-stranded converted DNA eluted with SSO_H_X",
                              MaterialKind.MOLECULAR, state=_ss("single-stranded converted DNA with physical U and retained methylcytosines"))
    column_waste = graph.product("column_waste", "spin-plate flow-through, desulphonation buffer and washes",
                                 MaterialKind.WASTE)
    graph.react("desulphonate", "Desulphonate and wash on the 384 spin plate; elute with SSO_H_X",
                StepKind.CLEANUP,
                ((converted_raw, InputRole.TEMPLATE), (column, InputRole.SELECTION)),
                (converted, column_waste))
    graph.discard(column_waste, "removed by centrifugation through the spin plate")
    klenow = graph.input("klenow", "Klenow exo− (2.5 U), 2× Blue buffer, 1 mM dNTP; 95 °C 2 min then ice, "
                         "4 °C 5 min, 25 °C 5 min, 37 °C 60 min", MaterialKind.CATALYST)
    barcoded = graph.product("barcoded", "SSO_H-primed barcoded second strand on the converted template",
                             MaterialKind.MOLECULAR, state=_primed_duplex())
    graph.react("strand_synthesis", "Random-prime second-strand synthesis with SSO_H_X and Klenow exo−",
                StepKind.EXTENSION,
                ((converted, InputRole.TEMPLATE), (sso, InputRole.STRUCTURAL),
                 (klenow, InputRole.CATALYST)), (barcoded,))
    exo = graph.input("exonuclease_i", "Exonuclease I (0.5 µl), 37 °C 30 min", MaterialKind.CATALYST)
    digested = graph.product("exo_treated", "barcoded duplex after excess-oligo digestion",
                             MaterialKind.MOLECULAR, state=barcoded.state)
    spent = graph.product("digested_sso", "Exo I-degraded excess SSO_H_X", MaterialKind.MOLECULAR)
    graph.react("exo_digest", "Digest excess strand-synthesis oligos with Exonuclease I",
                StepKind.DIGESTION,
                ((barcoded, InputRole.TEMPLATE), (exo, InputRole.CATALYST)), (digested, spent))
    graph.mark_terminal(Terminal(spent.id, Disposition.UNRESOLVED,"ExoI hydrolysis is not separation; exact class recovery through subsequent beads is unassigned."))
    spri1 = graph.input("spri_pool", "AMPure XP 1.2×; elute 10 µl 10 mM Tris–HCl pH 8.5",
                        MaterialKind.SELECTION)
    pooled = graph.product("pooled", "pooled barcoded reactions from up to 32 wells",
                           MaterialKind.MOLECULAR, state=barcoded.state)
    pool_waste = graph.product("pool_spri_waste", "1.2× SPRI supernatant and washes", MaterialKind.WASTE)
    graph.react("pool_purify", "Pool up to 32 single-cell reactions and SPRI-purify",
                StepKind.CLEANUP,
                ((digested, InputRole.TEMPLATE), (spri1, InputRole.SELECTION)), (pooled, pool_waste))
    graph.discard(pool_waste, "removed with the bead supernatant")
    ssb = graph.input("et_ssb", "95 °C 2 min with 50 ng ET-SSB, then ice", MaterialKind.CONDITION)
    single = graph.product("single_strand", "released barcoded single strand", MaterialKind.MOLECULAR,
                           state=_single("barcoded second strand", _barcoded_segments()))
    original = graph.product("original_strand", "released original converted strand (no SSO_H tail)",
                             MaterialKind.MOLECULAR, state=_ss("original converted strand"))
    graph.react("denature", "Heat-denature to separate the second strand from the original strand",
                StepKind.STRAND_SEPARATION,
                ((pooled, InputRole.TEMPLATE), (ssb, InputRole.CONDITION)), (single, original))
    graph.mark_terminal(Terminal("original_strand", Disposition.UNRESOLVED,
                                 "Carried into the bulk SPLAT reaction without physical separation. Missing SSO_H tail limits the intended primer-defined route, but does not establish zero off-target amplification or removal; native ligation/recovery unresolved."))
    adapter = graph.input("splat_adapter", "SPLAT 3′ adapter ss1a/ss1b, 5 µM final",
                          MaterialKind.ADAPTER, state=splat_adapter(),
                          source=PUBLISHED + ", Methods and Supplementary Table S1")
    ligase = graph.input("t4_ligase", "T4 DNA ligase HC (15 U), T4 PNK (10 U), 1× ligase buffer, PEG4000; 20 °C 1 h",
                         MaterialKind.CATALYST)
    ligated = graph.product("ligated", "SPLAT-ligated barcoded single strand", MaterialKind.MOLECULAR,
                            state=_ligated_strand())
    graph.react("splat", "SPLAT: splinted ligation of ss1b to the 3′ end in bulk",
                StepKind.LIGATION,
                ((single, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=adapter_geometry())
    spri2 = graph.input("spri_ligation", "SPRI 1:1; elute 10 µl", MaterialKind.SELECTION)
    clean = graph.product("ligated_clean", "purified SPLAT-ligated strands", MaterialKind.MOLECULAR,
                          state=ligated.state)
    lig_waste = graph.product("ligation_spri_waste", "1:1 SPRI supernatant, free adapter fraction and washes",
                              MaterialKind.WASTE)
    graph.react("ligation_cleanup", "SPRI-purify the ligation (1:1)", StepKind.CLEANUP,
                ((ligated, InputRole.TEMPLATE), (spri2, InputRole.SELECTION)), (clean, lig_waste))
    graph.discard(lig_waste, "removed with the bead supernatant")
    forward = graph.input("pcr_i5", "NEBNext i5 dual-index primer (sequence not printed)",
                          MaterialKind.OLIGO, inferred=True, source=NEBNEXT)
    reverse = graph.input("pcr_i7", "NEBNext i7 dual-index primer (sequence not printed)",
                          MaterialKind.OLIGO, inferred=True, source=NEBNEXT)
    kapa = graph.input("kapa", "KAPA HiFi 2× PCR mastermix, 12–16 cycles", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "dual-indexed scSPLAT pool library", MaterialKind.LIBRARY,
                              state=FINAL_LIBRARY, inferred=True, source=NEBNEXT)
    pcr_residual=graph.product('pcr_residual','carried native blocked template strands and off-target PCR species; recovery unassigned',MaterialKind.MOLECULAR,state=clean.state)
    graph.react("index_pcr", "INFERRED — Amplify with NEBNext dual-index primers (arms not printed)",
                StepKind.PCR,
                ((clean, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (amplified,pcr_residual),
                geometry=PCRGeometry.unavailable(
                    clean.state, reason="Actual single-strand SPLAT-ligated precursor enters PCR; native primer/first-copy geometry unassigned. Final library is a separate endpoint."))
    graph.mark_terminal(Terminal(pcr_residual.id,Disposition.UNRESOLVED,'PCR copying does not chemically erase original5AmMC6/3AmMo strands or imply off-target removal. Class-specific recovery through two final bead selections remains unassigned.'))
    spri3 = graph.input("spri_pcr", "0.8× SPRI reagent for two successive source recoveries; final10µlTris pH8.0",
                        MaterialKind.SELECTION)
    first_recovered=graph.product('pcr_cleanup_first_retained','first0.8×SPRI-recovered PCR library',MaterialKind.LIBRARY,state=FINAL_LIBRARY,inferred=True,source=NEBNEXT)
    first_waste=graph.product('pcr_cleanup_first_waste','first0.8×SPRI supernatant/washes',MaterialKind.WASTE)
    graph.react('pcr_cleanup_first','First post-PCR0.8×SPRI recovery',StepKind.CLEANUP,((amplified,InputRole.TEMPLATE),(spri3,InputRole.SELECTION)),(first_recovered,first_waste),changes_molecule=False)
    graph.discard(first_waste,'Physically removed first recovery supernatant/washes; class depletion unassigned.')
    library = graph.product("library", "purified scSPLAT sequencing library", MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY, inferred=True, source=NEBNEXT)
    pcr_waste = graph.product("pcr_spri_waste", "0.8× SPRI supernatants (primers, short products)",
                              MaterialKind.WASTE)
    graph.react("pcr_cleanup", "Second post-PCR0.8×SPRI recovery; final10µl elution", StepKind.CLEANUP,
                ((first_recovered, InputRole.TEMPLATE), (spri3, InputRole.SELECTION)), (library, pcr_waste))
    graph.discard(pcr_waste, "removed with the bead supernatants")
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS, limitation=SEQUENCING_UNAVAILABLE, acquired_roles=("Read 1", "Read 2"),
        # PE150 is stated; index-read cycles are not. Pools carrying unique dual
        # indices were combined before sequencing, so the index reads demultiplex pools.
        read_lengths={"Read 1": 150, "Read 2": 150},
        read_uses={"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
                   "Index 2 (i5)": "demultiplex", "Read 2": "include"}))
    graph.declare_reviewed_scope(REVIEWED_TARGETS, src)
    graph.annotate(SourceDiscrepancy(
        "splat_adapter", "Published SPLAT-adapter partner name differs",
        (SourceClaim("Scientific Reports Methods, adapter preparation (PMC8986790 JATS line 131)",
                     "The SPLAT adapter was prepared from oligo ss1a and ss2b (Supplementary Table S1)."),
         SourceClaim("Same paragraph, printed oligo list; Supplementary Table S1",
                     "The second oligo is printed as ss1b: 5′Phos/AGATCGGAAGAGCACACGTC/3′AmMo; no ss2b entry exists.")),
        "The model draws the printed ss1b sequence as the partner of ss1a; the printed sequences pair over 20 nt.",
        "Both statements are in the published article; only the name differs.",
        "Whether ss2b is a typo for ss1b or a distinct reagent is not stated; the original SPLAT paper (Raine 2017, NAR) has not been checked."))
    return graph.require_publishable()


REVIEWED_TARGETS = ('cell','bisulfite','convert','spin_plate','sso_h','desulphonate','klenow','strand_synthesis','exonuclease_i','exo_digest','spri_pool','pool_purify','et_ssb','denature','original_strand','splat_adapter','t4_ligase','splat','spri_ligation','ligation_cleanup','pcr_i5','pcr_i7','kapa','index_pcr','spri_pcr','pcr_cleanup_first','pcr_cleanup_first_retained','pcr_cleanup','library')
