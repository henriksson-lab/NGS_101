"""CapTrap-seq dual full-length selection and long-read endpoints.

Every oligo below is transcribed from the published Nature Communications Methods
(Carbonell-Sala et al. 2024, PMC11211341 JATS, "CapTrap-seq full-length cDNA libraries").
The two linker strands are printed there without 5'/3' labels; only the reversal of
each printed lower strand pairs antiparallel with its upper strand, so the lower
strands are reversed here (a derived reading, checked by ``Scene.anneal``).
"""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

TITLE = "CapTrap-seq — cap- and poly(A)-selected full-length cDNA"
NOTES = "01_captrap-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41467-024-49523-3">Carbonell-Sala et al., '
          '<i>Nature Communications</i> (2024)</a>, Methods "CapTrap-seq full-length cDNA libraries" '
          '(full protocol cited as Nature Protocol Exchange).')
AUTHORITY = ("Carbonell-Sala et al. 2024 Nat Commun 15:5278 doi:10.1038/s41467-024-49523-3, Methods "
             "'CapTrap-seq full-length cDNA libraries' (PMC11211341)")
SUMMARY = ("Anchored oligo(dT) selects the RNA 3-prime end, cap trapping selects cDNA that "
           "reached the native 5-prime cap, and two splinted single-strand linker ligations "
           "add the handles used for UMI second-strand synthesis and LA-PCR before long-read "
           "library preparation.")
CAVEAT = ("The published Methods print both strands of each linker without 5′/3′ labels. "
          "The model reverses each printed lower strand because only that reading pairs "
          "antiparallel with its upper strand (derived). The cDNA-terminal bases opposite "
          "the 5′ linker's G-N5 tail and the 3′ linker's A5 tail are not drawn as sequence. "
          "The ONT SQK-LSK108/109 adapter-ligation step is shown as inferred: the article "
          "names the kits but not their adapter bases or end geometry. The same LA-PCR "
          "product was alternatively converted with the PacBio SMRTbell Express Template "
          "Prep Kit 2.0 (not drawn). The LA-PCR product is drawn as a full duplex; the "
          "article does not report its sequencing read layout beyond the platform kits.")

# Published Methods, verbatim bases (modifications as separate fields below).
OLIGO_DT_16VN = "TTTTTTTTTTTTTTTTVN"             # 5' [Phos]
LINKER5_UPPER = "GTGGTAUCAACGCAGAGUACGNNNNN"     # 3' -P
LINKER5_LOWER_PRINTED = "CACCATAGTTGCGTCTCATG"   # P- ... -P
LINKER3_UPPER = "AAAAAGCAUCGCUGTCTCUTAUACACAUCUCCGAGCCCACGAGAC"  # 3' -P
LINKER3_LOWER_PRINTED = "CGTAGCGACAGAGAATATGTGTAGAGGCTCGGGTGCTCTG"
UMI8_PRIMER = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNGTGGTATCAACGCAGAGTAC"
LA_FORWARD = "TCGTCGGCAGCGTC"
LA_REVERSE = "GTCTCGTGGGCTCGG"

# Derived 5'->3' orientation of the printed lower strands (see module docstring).
LINKER5_LOWER = LINKER5_LOWER_PRINTED[::-1]
LINKER3_LOWER = LINKER3_LOWER_PRINTED[::-1]

UMI = feature("umi", "umi", "random")


def _umi8_segments():
    handle = UMI8_PRIMER[:len(LA_FORWARD)]
    if handle != LA_FORWARD:
        raise ValueError("UMI8 primer no longer starts with the forward LA-PCR primer")
    i = UMI8_PRIMER.index("N")
    j = UMI8_PRIMER.rindex("N") + 1
    return [Segment("forward LA-PCR primer site", handle, "r1"),
            Segment("UMI8-primer spacer", UMI8_PRIMER[len(handle):i]),
            Segment("UMI (8 nt)", UMI8_PRIMER[i:j], "umi", feature=UMI),
            Segment("5′-linker-complementary 3′ end", UMI8_PRIMER[j:], "tso")]


def _cdna(n=30):
    return Segment("cap-selected cDNA", "X" * n, placeholder=True)


def _dt_cdna_segments():
    return [Segment("oligo(dT)16", OLIGO_DT_16VN[:16]),
            Segment("anchor VN", OLIGO_DT_16VN[16:], placeholder=True, bottom="BN"), _cdna()]


def five_linker_geometry() -> AdapterLigationGeometry:
    sc = Scene()
    sc.strand("5′ linker upper", [Segment("5′-linker paired handle", LINKER5_UPPER[:20], "tso"),
                                   Segment("G-N5 splint tail", LINKER5_UPPER[20:])],
              label="5′ linker upper (U-containing)", mod3="P")
    sc.anneal("5′ linker lower", [Segment("5′-linker lower strand", LINKER5_LOWER, "tso")],
              to="5′ linker upper", pair=("5′-linker lower strand", "5′-linker paired handle"),
              label="5′ linker lower", mod5="P", mod3="P")
    sc.strand("cDNA substrate", [_cdna(28)], label="released single-strand cDNA")
    sc.ligation_end("5′ linker lower", "5′-linker lower strand", kind="splinted",
                    polarity="5-prime", owner="adapter",
                    text="5′-phosphate beside the G-N5 tail")
    sc.ligation_end("cDNA substrate", "cap-selected cDNA", kind="splinted",
                    polarity="3-prime", owner="substrate", text="cDNA 3′ end")
    sc.labels("5′ linker upper")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement="Methods print GTGGTAUCAACGCAGAGUACGNNNNN-P over P-CACCATAGTTGCGTCTCATG-P; read "
                  "antiparallel (derived), the 3′-blocked G-N5 tail splints the lower strand's "
                  "5′ phosphate against the single-strand cDNA 3′ end.")


def three_linker_geometry() -> AdapterLigationGeometry:
    sc = Scene()
    sc.strand("3′ linker upper", [Segment("A5 splint tail", LINKER3_UPPER[:5]),
                                   Segment("3′-linker paired body", LINKER3_UPPER[5:], "p7")],
              label="3′ linker upper (U-containing)", mod3="P")
    sc.anneal("3′ linker lower", [Segment("3′-linker lower strand", LINKER3_LOWER, "p7")],
              to="3′ linker upper", pair=("3′-linker lower strand", "3′-linker paired body"),
              label="3′ linker lower")
    sc.strand("cDNA substrate", [Segment("oligo(dT)16", OLIGO_DT_16VN[:16]),
                                 Segment("anchor VN", OLIGO_DT_16VN[16:], placeholder=True), _cdna(12)],
              label="5′-linked cDNA, oligo(dT) end", mod5="P")
    sc.ligation_end("3′ linker lower", "3′-linker lower strand", kind="splinted",
                    polarity="3-prime", owner="adapter", text="lower-strand 3′-OH beside the A5 tail")
    sc.ligation_end("cDNA substrate", "oligo(dT)16", kind="splinted", polarity="5-prime",
                    owner="substrate", text="5′-phosphorylated oligo(dT) primer end")
    sc.labels("3′ linker upper")
    return AdapterLigationGeometry.known(
        "splinted", sc,
        statement="Methods print the 3′-blocked upper strand AAAAAGCAUCGCUG…GAC-P over a lower "
                  "strand without end modification; read antiparallel (derived), the 5′ A5 tail "
                  "splints the lower strand's 3′-OH against the 5′-phosphorylated OligodT 16VN end.")


def five_linked_scene() -> Scene:
    sc = Scene()
    sc.strand("cDNA", [*_dt_cdna_segments(), Segment("5′-linker lower strand", LINKER5_LOWER, "tso")],
              label="5′-linked cDNA", mod5="P", mod3="P")
    sc.anneal("5′ linker upper", [Segment("5′-linker paired handle", LINKER5_UPPER[:20], "tso"),
                                   Segment("G-N5 splint tail", LINKER5_UPPER[20:])],
              to="cDNA", pair=("5′-linker paired handle", "5′-linker lower strand"),
              label="upper strand (not ligated)", mod3="P", unpaired=("G-N5 splint tail",))
    sc.junction("cDNA", "cap-selected cDNA", "5′-linker lower strand", "splinted ligation")
    sc.labels("cDNA")
    return sc


def handled_cdna_segments():
    return [Segment("3′-linker lower strand", LINKER3_LOWER, "p7"), *_dt_cdna_segments(),
            Segment("5′-linker lower strand", LINKER5_LOWER, "tso")]


def doubly_linked_scene() -> Scene:
    sc = Scene()
    sc.strand("cDNA", handled_cdna_segments(), label="doubly linked cDNA", mod3="P")
    sc.junction("cDNA", "3′-linker lower strand", "oligo(dT)16", "3′ linker ligation")
    sc.junction("cDNA", "cap-selected cDNA", "5′-linker lower strand", "5′ linker ligation")
    sc.labels("cDNA")
    return sc


def second_strand_priming_scene() -> Scene:
    sc = Scene()
    sc.strand("cDNA", handled_cdna_segments(), label="SAP/USER-treated cDNA")
    sc.anneal("UMI8 primer", _umi8_segments(), to="cDNA",
              pair=("5′-linker-complementary 3′ end", "5′-linker lower strand"),
              label="second primer UMI8")
    sc.labels("cDNA")
    return sc


LA_PCR_PRODUCT = Construct(
    [*_umi8_segments(), Segment("cDNA (mRNA sense)", "X" * 30, placeholder=True),
     Segment("anchor complement", "NB", placeholder=True, bottom="NV"),
     Segment("poly(A) copy", revcomp(OLIGO_DT_16VN[:16])),
     Segment("3′-linker-derived end", revcomp(LINKER3_LOWER), "p7")],
    name="LA-PCR-amplified CapTrap cDNA")
FINAL_LIBRARY = Construct(
    [Segment("ONT adapter", "[SQK-LSK108/109 adapter]", "r2", placeholder=True, inferred=True),
     *LA_PCR_PRODUCT],
    name="ONT-adapted CapTrap-seq cDNA")
FINAL_CAPTION = ("INFERRED — LA-PCR product after ONT adapter ligation; the SQK-LSK108/109 "
                 "end-prep and adapter bases and the adapter positions are not given in the article.")
SEQ_PRIMERS = ()
SEQUENCING_ENDING = ("There is no synthesis sequencing primer on the ONT branch: the article states "
                     "that 500–1000 ng cDNA was prepared strictly following the SQK-LSK108 and "
                     "SQK-LSK109 adapter-ligation protocols and run on MinION R9.4 flow cells. "
                     "For PacBio Sequel I/II the same cDNA was converted with the SMRTbell Express "
                     "Template Prep Kit 2.0. No size selection was applied before loading.")


def oligos():
    return [
        "<b>OligodT 16VN</b> — 5′-[Phos]TTTTTTTTTTTTTTTTVN-3′",
        "<b>5′ linker, upper strand</b> — GTGGTAUCAACGCAGAGUACGNNNNN-P (as printed)",
        "<b>5′ linker, lower strand</b> — P-CACCATAGTTGCGTCTCATG-P (as printed; aligned 3′→5′ under the upper strand, derived)",
        "<b>3′ linker, upper strand</b> — AAAAAGCAUCGCUGTCTCUTAUACACAUCUCCGAGCCCACGAGAC-P (as printed)",
        "<b>3′ linker, lower strand</b> — CGTAGCGACAGAGAATATGTGTAGAGGCTCGGGTGCTCTG (as printed; aligned 3′→5′, derived)",
        "<b>second primer UMI8</b> — 5′-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNGTGGTATCAACGCAGAGTAC-3′",
        "<b>forward LA-PCR primer</b> — 5′-TCGTCGGCAGCGTC-3′",
        "<b>reverse LA-PCR primer</b> — 5′-GTCTCGTGGGCTCGG-3′",
    ]


def protocol_graph():
    source = AUTHORITY
    ont_source = source + " (ONT SQK-LSK108/109 adapter bases are not disclosed)"
    g = ProtocolGraph("CapTrap-seq", authority=source)
    rna = g.input("rna", "5 µg total RNA (capped, polyadenylated mRNA analysed)",
                  MaterialKind.BIOLOGICAL, primary=True)
    oligo_dt = g.input("anchored_oligo_dt", "OligodT 16VN (5′-phosphorylated)", MaterialKind.OLIGO,
                       state=OLIGO_DT_16VN)
    rt = g.input("reverse_transcriptase", "PrimeScript II reverse transcriptase", MaterialKind.CATALYST)
    hybrid = g.product("rna_cdna_hybrid", "RNA–cDNA hybrids (42 °C 60 min; RNAClean 1.8x)",
                       MaterialKind.MOLECULAR)
    g.react("reverse_transcription", "Prime at poly(A) and synthesize first-strand cDNA",
            StepKind.REVERSE_TRANSCRIPTION,
            ((rna, InputRole.TEMPLATE), (oligo_dt, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (hybrid,))
    periodate = g.input("periodate", "NaIO4 in NaOAc pH 4.5, 5 min on ice; Tris pH 8.5 stop",
                        MaterialKind.CONDITION)
    oxidized = g.product("oxidized_hybrids", "hybrids with RNA diols (m7G cap, RNA 3′ end) oxidized to aldehydes",
                         MaterialKind.MOLECULAR)
    g.react("oxidation", "Oxidize RNA diols with periodate", StepKind.CHEMICAL_CONVERSION,
            ((hybrid, InputRole.TEMPLATE), (periodate, InputRole.CONDITION)), (oxidized,))
    hydrazide = g.input("biotin_hydrazide", "Biotin (Long Arm) Hydrazide, NaOAc pH 6.0, 30 min 40 °C",
                        MaterialKind.CONDITION)
    biotinylated = g.product("biotinylated_hybrids", "biotinylated hybrids (cap and RNA 3′ end)",
                             MaterialKind.MOLECULAR)
    g.react("biotinylation", "Biotinylate the aldehydes", StepKind.CHEMICAL_CONVERSION,
            ((oxidized, InputRole.TEMPLATE), (hydrazide, InputRole.CONDITION)), (biotinylated,))
    rnase_one = g.input("rnase_one", "RNase ONE (single-strand RNA), 30 min 37 °C", MaterialKind.CATALYST)
    trimmed = g.product("trimmed_hybrids",
                        "hybrids keeping biotin only where cDNA protects the capped 5′ end",
                        MaterialKind.MOLECULAR)
    g.react("rnase_one_trim", "Digest unprotected single-strand RNA and its 3′-end biotin",
            StepKind.DIGESTION, ((biotinylated, InputRole.TEMPLATE), (rnase_one, InputRole.CATALYST)),
            (trimmed,))
    beads = g.input("streptavidin", "M-270 streptavidin beads in LiCl binding buffer",
                    MaterialKind.SELECTION)
    captured = g.product("captured_hybrids", "bead-bound cap-complete hybrids", MaterialKind.MOLECULAR)
    unbound = g.product("unbound_hybrids", "uncapped or truncated (non-biotinylated) hybrids",
                        MaterialKind.WASTE)
    g.react("cap_capture", "Capture biotinylated cap-complete hybrids and wash three times",
            StepKind.AFFINITY_CAPTURE, ((trimmed, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
            (captured, unbound))
    g.discard(unbound, "removed with the supernatant and three CapTrap washes")
    release = g.input("release_digest", "two 95 °C releases in RNase ONE buffer; RNase H + RNase ONE 30 min 37 °C",
                      MaterialKind.CATALYST)
    single = g.product("released_cdna", "single-strand cap-selected cDNA (AMPure 1.8x, dried)",
                       MaterialKind.MOLECULAR)
    g.react("release", "Heat-release cDNA and digest the RNA", StepKind.DIGESTION,
            ((captured, InputRole.TEMPLATE), (release, InputRole.CATALYST)), (single,))
    linker5 = g.input("linker5", "5′ linker duplex (denatured 5 min 55 °C)", MaterialKind.ADAPTER,
                      state=(LINKER5_UPPER, LINKER5_LOWER))
    ligase = g.input("ligase", "Mighty Mix DNA ligase (4 h 30 °C or 16 h 16 °C)", MaterialKind.CATALYST)
    linked5 = g.product("five_linked_cdna", "5′-linked cDNA", MaterialKind.MOLECULAR,
                        state=five_linked_scene())
    free5 = g.product("free_linker5", "non-incorporated 5′ linker", MaterialKind.WASTE)
    g.react("linker5_ligation", "Ligate the 5′ linker to the cDNA 3′ end", StepKind.LIGATION,
            ((single, InputRole.TEMPLATE), (linker5, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (linked5, free5), geometry=five_linker_geometry())
    g.discard(free5, "removed by two 1.8x AMPure XP purifications")
    linker3 = g.input("linker3", "3′ linker duplex (denatured 5 min 65 °C)", MaterialKind.ADAPTER,
                      state=(LINKER3_UPPER, LINKER3_LOWER))
    linked = g.product("doubly_linked_cdna", "doubly linked full-length cDNA", MaterialKind.MOLECULAR,
                       state=doubly_linked_scene())
    g.react("linker3_ligation", "Ligate the 3′ linker to the oligo(dT) 5′ end", StepKind.LIGATION,
            ((linked5, InputRole.TEMPLATE), (linker3, InputRole.STRUCTURAL),
             (ligase, InputRole.CATALYST)), (linked,), geometry=three_linker_geometry())
    sap_user = g.input("sap_user", "Shrimp Alkaline Phosphatase + USER, 30 min 37 °C then 5 min 95 °C",
                       MaterialKind.CATALYST)
    handled = g.product("handled_cdna", "single-strand cDNA carrying both lower-strand handles",
                        MaterialKind.MOLECULAR, state=second_strand_priming_scene())
    fragments = g.product("upper_strand_fragments", "U-excised upper linker strands and unligated linker",
                          MaterialKind.WASTE)
    g.react("sap_user_treatment", "Dephosphorylate and excise uracils to make single-strand handles",
            StepKind.DIGESTION, ((linked, InputRole.TEMPLATE), (sap_user, InputRole.CATALYST)),
            (handled, fragments))
    g.discard(fragments, "removed by 1.8x AMPure XP purification")
    umi_primer = g.input("umi8_primer", "second primer UMI8", MaterialKind.OLIGO, state=UMI8_PRIMER)
    kapa = g.input("kapa_hifi", "2x KAPA HiFi mix with DMSO (95 °C 5 min, 55 °C 5 min, 72 °C 30 min)",
                   MaterialKind.CATALYST)
    duplex = g.product("ds_cdna", "double-stranded UMI-tagged cDNA", MaterialKind.MOLECULAR)
    g.react("second_strand", "Prime on the 5′-linker handle and synthesize the second strand",
            StepKind.EXTENSION, ((handled, InputRole.TEMPLATE), (umi_primer, InputRole.STRUCTURAL),
                                 (kapa, InputRole.CATALYST)), (duplex,))
    exo = g.input("exonuclease_i", "Exonuclease I, 30 min 37 °C; AMPure 1.8x then 1.4x",
                  MaterialKind.CATALYST)
    cleaned = g.product("cleaned_second_strand", "purified double-stranded cDNA", MaterialKind.MOLECULAR)
    spent = g.product("excess_primer", "Exo I-degraded excess single-strand primer", MaterialKind.WASTE)
    g.react("exo_cleanup", "Degrade excess primer and purify", StepKind.DIGESTION,
            ((duplex, InputRole.TEMPLATE), (exo, InputRole.CATALYST)), (cleaned, spent))
    g.discard(spent, "nucleotides and short fragments removed by the two AMPure XP purifications")
    la_fwd = g.input("la_pcr_forward", "forward LA-PCR primer", MaterialKind.OLIGO, state=LA_FORWARD)
    la_rev = g.input("la_pcr_reverse", "reverse LA-PCR primer", MaterialKind.OLIGO, state=LA_REVERSE)
    la_taq = g.input("la_taq", "TaKaRa LA Taq, 16 cycles, two independent reactions merged",
                     MaterialKind.CATALYST)
    enriched = g.product("enriched_cdna", LA_PCR_PRODUCT.name, MaterialKind.MOLECULAR,
                         state=LA_PCR_PRODUCT)
    g.react("la_pcr", "Amplify between the two handles by LA-PCR (1x AMPure XP)", StepKind.PCR,
            ((cleaned, InputRole.TEMPLATE), (la_fwd, InputRole.STRUCTURAL), (la_rev, InputRole.STRUCTURAL),
             (la_taq, InputRole.CATALYST)), (enriched,),
            geometry=PCRGeometry.exact(
                LA_PCR_PRODUCT, PCRPrimer("forward LA-PCR primer", LA_FORWARD, anneal_nt=len(LA_FORWARD)),
                PCRPrimer("reverse LA-PCR primer", LA_REVERSE, anneal_nt=len(LA_REVERSE)),
                statement="Methods print TCGTCGGCAGCGTC (the UMI8-primer 5′ end) and GTCTCGTGGGCTCGG "
                          "(the 5′ end of the derived 3′-linker lower strand)."))
    ont_adapter = g.input("ont_adapter", "SQK-LSK108/109 end-prep and adapter-ligation reagents",
                          MaterialKind.ADAPTER, inferred=True, source=ont_source)
    ont_ligase = g.input("ont_ligase", "ONT kit ligation reagents", MaterialKind.CATALYST,
                         inferred=True, source=ont_source)
    library = g.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY,
                        inferred=True, source=ont_source)
    g.react("ont_ligation", "Prepare the ONT adapter-ligation library", StepKind.LIGATION,
            ((enriched, InputRole.TEMPLATE), (ont_adapter, InputRole.STRUCTURAL),
             (ont_ligase, InputRole.CATALYST)), (library,),
            geometry=AdapterLigationGeometry.unavailable(
                "The article names ONT SQK-LSK108 and SQK-LSK109 but not their end-prep or "
                "motor-adapter end sequence or geometry; the kit manuals were not read."),
            source=ont_source)
    g.sequence(library, SequencingPlan(limitation=SEQUENCING_ENDING))
    g.declare_reviewed_scope((
        "rna", "anchored_oligo_dt", "reverse_transcriptase", "reverse_transcription",
        "periodate", "oxidation", "biotin_hydrazide", "biotinylation", "rnase_one",
        "rnase_one_trim", "streptavidin", "cap_capture", "release_digest", "release",
        "linker5", "ligase", "linker5_ligation", "linker3", "linker3_ligation", "sap_user", "sap_user_treatment",
        "umi8_primer", "kapa_hifi", "second_strand", "exonuclease_i", "exo_cleanup",
        "la_pcr_forward", "la_pcr_reverse", "la_taq", "la_pcr", "ont_adapter", "ont_ligase", "ont_ligation",
        "library"), source)
    return g.require_publishable()
