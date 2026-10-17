"""SNARE-seq RNA and accessible-chromatin molecular states.

Authority: Chen, Lake & Zhang, Nat. Biotechnol. 37:1452 (2019), PMC6893138 JATS
Methods (tagmentation, nuclei barcoding, library preparation, NGS sequencing) and
Supplementary Table 4 ("Oligonucleotide sequences").  Drop-seq (Macosko 2015) is the
cited predecessor for droplets, bead primers and the Nextera XT cDNA library.
"""
import illumina as il
import nextera as nx
import seqprimers as sp

from chemdraw import Construct, Row, Scene, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, Input, InputRole, Material, MaterialKind,
                            ProtocolGraph, SequencingPlan, Step, StepKind, Terminal)

TITLE = "SNARE-seq — transcriptome and chromatin accessibility"
NOTES = "01_snare-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41587-019-0290-0">'
          'Chen et al., <i>Nature Biotechnology</i> (2019)</a>.')
SUMMARY = ("A phosphorylated poly(A) splint anneals across the complete Nextera S5+ME "
           "end of pre-tagmented chromatin. A Drop-seq bead oligo captures that tail and "
           "nuclear mRNA, transferring the same cell barcode to both branches.")
CAVEAT = ('<b>Published table cross-reference error:</b> The <a href="https://doi.org/10.1038/s41587-019-0290-0">article Methods</a> name “Table S3” for Nextera-R1-rc-polyA. '
          '<a href="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41587-019-0290-0/MediaObjects/41587_2019_290_MOESM5_ESM.xlsx">Table S3</a> instead contains differential-accessibility results. '
          'The actual /5Phos/ Nextera-R1-rc-polyA sequence is in the “Oligonucleotide sequences” sheet of '
          '<a href="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41587-019-0290-0/MediaObjects/41587_2019_290_MOESM6_ESM.xlsx">Table S4</a>; this model uses that printed strand. The source does not identify a particular library’s index-primer member. '
          '<br><b>Inferred RNA-library ends:</b> the cDNA library is made “with Nextera XT kit as described previously” (Drop-seq); its P5-side primer, '
          'the bases joining P5 to the Read1CustSeqB site, and the i7 adapter are not printed in the SNARE-seq sources, so those regions are drawn as inferred. '
          'HP10/HP11 are Illumina reagent primer mixes named, not printed, by the paper; their Nextera Read 1/Read 2 sites are inferred.')

BEAD_SMART = "AAGCAGTGGTATCAACGCAGAGTAC"
BEAD_SPACER = "T" * 8
CELL_NT, UMI_NT = 12, 8
# Supplementary Table 4 prints the 63-nt splint as ...CCGACG + A31; the first of those
# As completes the S5 complement, leaving a 30-nt poly(A) beyond the 33-nt Nextera arm.
BEAD_DT, SPLINT_POLYA = "T" * 30, "A" * 30
SPLINT = nx.ME_RC + nx.S5_RC + SPLINT_POLYA
NEXTERA_R1_BK = nx.S5 + nx.ME
NEXTERA_AD2_BK = nx.S7_RC
TSO = "AAGCAGTGGTATCAACGCAGAGTGAAT"  # followed by rGrGrG as ordered
TSO_PCR = "AAGCAGTGGTATCAACGCAGAGT"
NEXTERA_R2 = nx.S7 + nx.ME
# Supplementary Table 4, verbatim (first member of each index family; *: phosphorothioate).
P501_TSO = "AATGATACGGCGACCACCGAGATCTACACTAGATCGCGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC"
AD2_1 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGT"
READ1_CUSTOM = "GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC"
TSO_INDEX_SPACER = P501_TSO[len(il.P5) + 8:-len(BEAD_SMART)]
if (P501_TSO[:len(il.P5)] != il.P5 or TSO_INDEX_SPACER + BEAD_SMART != READ1_CUSTOM
        or revcomp(il.P7_RC) + "TCGCCTTA" + nx.S7 + nx.ME[:6] != AD2_1):
    raise AssertionError("Supplementary Table 4 index primers changed shape")


def _cell_segment():
    return Segment("12-nt cell barcode", "N" * CELL_NT, "cbc",
                   feature=feature("cell_barcode", "cell_barcode", "combinatorial",
                                   note="12 split-pool synthesis cycles (Drop-seq); "
                                   "ChemGenes Macosko201110 beads, no fixed list"))


def _umi_segment():
    return Segment("8-nt UMI", "N" * UMI_NT, "umi",
                   feature=feature("umi", "umi", "random"))


def _bead_segments():
    return [Segment("bead spacer", BEAD_SPACER),
            Segment("bead SMART handle", BEAD_SMART, "r1"),
            _cell_segment(), _umi_segment(), Segment("bead oligo-dT", BEAD_DT)]


def _tagged_segments():
    return [Segment("S5 handle", nx.S5, "s5"),
            Segment("left mosaic end", nx.ME, "me"),
            Segment("accessible genomic fragment", "X" * 42, placeholder=True),
            Segment("right mosaic-end complement", nx.ME_RC, "me"),
            Segment("S7-handle complement", nx.S7_RC, "s7")]


def tagmented_scene():
    sc = Scene.duplex(_tagged_segments(), label="productive TDE1-tagged chromatin")
    sc.junction("top", "left mosaic end", "accessible genomic fragment", "Tn5 transfer")
    sc.junction("top", "accessible genomic fragment", "right mosaic-end complement",
                "Tn5 transfer")
    sc.note("top", "Only an S5 end binds the capture splint; later PCR selects the asymmetric S5/S7 product.")
    sc.labels("top")
    return sc


def capture_scene():
    """Exact three-strand S5-end/splint/bead capture intermediate."""
    sc = Scene()
    sc.strand("tagged DNA", [Segment("S5 handle", nx.S5, "s5"),
                             Segment("mosaic end", nx.ME, "me"),
                             Segment("accessible genomic DNA", "X" * 36,
                                     placeholder=True)],
              label="tagmented accessible DNA")
    sc.anneal("phosphorylated capture splint", [
        Segment("mosaic-end reverse complement", nx.ME_RC, "me"),
        Segment("S5-handle reverse complement", nx.S5_RC, "s5"),
        Segment("30-nt poly(A) capture tail", SPLINT_POLYA),
    ], to="tagged DNA", pair=("mosaic-end reverse complement", "mosaic end"),
              label="Nextera-R1-rc-polyA", mod5="/5Phos/")
    sc.anneal("bead oligo", _bead_segments(), to="phosphorylated capture splint",
              pair=("bead oligo-dT", "30-nt poly(A) capture tail"),
              label="Drop-seq bead oligo")
    sc.mark("phosphorylated capture splint", "30-nt poly(A) capture tail",
            "hybridization transfers this chromatin end to the barcode bead")
    sc.labels("phosphorylated capture splint")
    sc.labels("bead oligo")
    return sc


def atac_template():
    return Construct([
        *_bead_segments(),
        Segment("copied S5 handle", nx.S5, "s5"),
        Segment("copied left mosaic end", nx.ME, "me"),
        Segment("accessible genomic fragment", "X" * 42, placeholder=True),
        Segment("right mosaic-end complement", nx.ME_RC, "me"),
        Segment("S7-handle complement", nx.S7_RC, "s7"),
    ], name="gap-filled and ligated SNARE-seq chromatin template")


def gap_filled_scene():
    sc = Scene.duplex(list(atac_template()), label="barcoded chromatin PCR template")
    sc.junction("bottom", "copied left mosaic end'", "accessible genomic fragment'",
                "T7 ligation after gap fill")
    sc.mark("top", "12-nt cell barcode", "same bead identity as the RNA branch")
    sc.labels("top")
    return sc


def preamp_product():
    # The PCR starts inside the SMART handle, so the bead's eight-T spacer is absent.
    return Construct(list(atac_template())[1:], name="SNARE-seq common ATAC preamp product")


def preamp_scene():
    sc = Scene.duplex(list(preamp_product()), label="Nextera-R2 / TSO-PCR amplicon")
    sc.mark("top", "12-nt cell barcode", "cell identity")
    sc.mark("top", "8-nt UMI", "used for RNA molecules")
    sc.labels("top")
    return sc


def atac_library():
    """P5xx-Tso / Ad2.x chromatin library (Supplementary Table 4 primers)."""
    return Construct([
        Segment("P5", il.P5, "p5"),
        Segment("i5 (P50x-Tso)", "J" * 8, "cbc", placeholder=True,
                feature=feature("sample_index_i5", "sample_index", "whitelist",
                                whitelist="Supplementary Table 4 P501–P506-Tso")),
        Segment("P5xx-Tso linker", TSO_INDEX_SPACER),
        *preamp_product(),
        Segment("i7 (Ad2.x)", "I" * 8, "cbc", placeholder=True,
                feature=feature("sample_index_i7", "sample_index", "whitelist",
                                whitelist="Supplementary Table 4 Ad2.1–Ad2.6")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name="SNARE-seq accessibility library")


def rna_library():
    """Nextera XT cDNA library; adapter ends beyond the bead handle are inferred."""
    return Construct([
        Segment("P5 (inferred)", il.P5, "p5", inferred=True),
        Segment("P5-to-SMART linker (inferred)", TSO_INDEX_SPACER, inferred=True),
        Segment("bead SMART handle", BEAD_SMART, "r1"),
        _cell_segment(), _umi_segment(), Segment("bead oligo-dT", BEAD_DT),
        Segment("cDNA (3′ end)", "X" * 42, placeholder=True),
        Segment("Nextera XT mosaic-end complement (inferred)", nx.ME_RC, "me", inferred=True),
        Segment("S7 complement (inferred)", nx.S7_RC, "s7", inferred=True),
        Segment("i7 (inferred, not read)", "I" * 8, "cbc", placeholder=True, inferred=True,
                feature=feature("sample_index_i7", "sample_index", "unknown")),
        Segment("P7 reverse complement (inferred)", il.P7_RC, "p7", inferred=True),
    ], name="SNARE-seq RNA library")


_SRC = "Chen et al. 2019 Methods (NGS sequencing) and Supplementary Table 4"
READ1_CUST = sp.custom("Read 1", "Read1CustSeqB", READ1_CUSTOM, _SRC)
HP10 = sp.custom("Index 1 (i7)", "HP10 in the Index 1 slot (Nextera Read 1 site, inferred)",
                 nx.READ1_PRIMER, "Chen et al. 2019 Methods names HP10; contents not printed",
                 note="75-cycle genomic read from the S5/ME end")
HP11 = sp.custom("Read 2", "HP11 (Nextera Read 2 site, inferred)", nx.READ2_PRIMER,
                 "Chen et al. 2019 Methods names HP11; contents not printed")
RNA = rna_library()
ATAC = atac_library()
RP = (READ1_CUST, HP11)
AP = (READ1_CUST, HP10, HP11)
FINAL_LIBRARIES = (
    ("RNA library", RNA, RP,
     "Read 1 (30 cycles, Read1CustSeqB) reads the 12-nt cell barcode, 8-nt UMI and oligo-dT; "
     "Read 2 (80 cycles, HP11) reads cDNA. Ends beyond the bead handle are inferred.",
     "HiSeq 2500: Read1CustSeqB for Read 1 (30 bp) and HP11 for Read 2 (80 bp); no index read is stated."),
    ("Accessibility library", ATAC, AP,
     "Read 1 (30 cycles) reads the same bead barcode and UMI; HP10 in the Index 1 slot and HP11 "
     "read 75 bp from each genomic end. The 8-bp Index 2 (i5) read is stated but its "
     "priming site is not printed, so it is not drawn.",
     "HiSeq 2500: Read1CustSeqB Read 1 (30 bp), HP10 Index 1 (75 bp), HP11 Read 2 (75 bp), "
     "plus an 8-bp Index 2 read for demultiplexing (not drawn; primer not printed)."),
)

# Methods "NGS Sequencing": the chromatin run also acquires an 8-bp Index 2 read whose
# priming site is not printed; it is declared as acquired without a represented site.
ATAC_INDEX2_READ = sp.read_spec("Index 2 (i5)", 8, "P50x-Tso i5 index", "demultiplex")


def oligos():
    return [f"<b>Nextera-R1-rc-polyA</b> — 5′-/5Phos/{SPLINT}-3′",
            f"<b>Nextera-R1-bk</b> — 5′-/5Phos/{NEXTERA_R1_BK}/3InvdT/-3′",
            f"<b>Nextera-Ad2-bk</b> — 5′-/5Phos/{NEXTERA_AD2_BK}/3InvdT/-3′",
            f"<b>TSO</b> — 5′-{TSO}rGrGrG-3′",
            f"<b>Tso-PCR</b> — 5′-{TSO_PCR}-3′",
            f"<b>Nextera-R2</b> — 5′-{NEXTERA_R2}-3′",
            f"<b>P501-Tso, example i5 primer</b> — 5′-{P501_TSO[:-2]}*A*C-3′",
            f"<b>Ad2.1_TAAGGCGA, example i7 primer</b> — 5′-{AD2_1}-3′",
            f"<b>Read1CustomSeqB</b> — 5′-{READ1_CUSTOM}-3′"]


def protocol_graph() -> ProtocolGraph:
    """RNA and chromatin captured on one bead, co-amplified, then split in two halves."""
    source = "Chen et al. 2019, doi:10.1038/s41587-019-0290-0"
    graph = ProtocolGraph(TITLE, authority=source)

    def add_input(key, label, kind, *, primary=False, state=None):
        graph.add_input(Material(key, label, kind, primary=primary, source=source, state=state))

    def add_step(key, label, kind, ports, products, *, geometry=None):
        products = tuple(products)
        graph.add_step(Step(key, label, kind, tuple(ports),
                            tuple(x.id for x in products), geometry=geometry,
                            source=source), products)

    def mol(key, label, state=None, kind=MaterialKind.MOLECULAR):
        return Material(key, label, kind, state=state, source=source)

    add_input("nucleus", "~200,000 permeabilized nuclei (poly(A) RNA + chromatin)",
              MaterialKind.BIOLOGICAL, primary=True)
    add_input("loaded_tn5_adapters", "TDE1-loaded Nextera adapters", MaterialKind.ADAPTER)
    add_input("tde1", "TDE1 Tagment DNA Enzyme (37 °C, 30 min)", MaterialKind.CATALYST)
    tagged = mol("tagmented_nucleus", "nucleus with productively tagmented chromatin",
                 tagmented_scene())
    tag_waste = mol("non_library_tn5_products",
                    "symmetric S5/S5 and S7/S7 Tn5 products", kind=MaterialKind.WASTE)
    add_step("tagment", "Tagment accessible chromatin in nuclei", StepKind.TAGMENTATION,
             (Input("nucleus", InputRole.TEMPLATE),
              Input("loaded_tn5_adapters", InputRole.STRUCTURAL),
              Input("tde1", InputRole.CATALYST)), (tagged, tag_waste))
    graph.mark_terminal(Terminal(tag_waste.id, Disposition.DISCARD,
                                 "candidate: not the asymmetric S5/S7 product; not selectively "
                                 "amplified by the S5-side capture and Nextera-R2/Tso-PCR pair"))

    add_input("capture_splint", "phosphorylated Nextera-R1-rc-polyA splint (in lysis buffer)",
              MaterialKind.OLIGO, state=SPLINT)
    add_input("barcode_bead_oligo", "Drop-seq barcoded bead primers (ChemGenes Macosko201110)",
              MaterialKind.OLIGO)
    captured = mol("captured_nucleus", "RNA and chromatin captured on one barcode bead",
                   capture_scene())
    add_step("capture", "Co-encapsulate nuclei and beads; 72 °C 5 min lysis releases Tn5; "
             "splint and poly(A) RNA hybridize to the bead", StepKind.HYBRIDIZATION,
             (Input(tagged.id, InputRole.TEMPLATE),
              Input("capture_splint", InputRole.STRUCTURAL),
              Input("barcode_bead_oligo", InputRole.STRUCTURAL)), (captured,))
    add_input("nextera_r1_bk", "Nextera-R1-bk blocking oligo (/5Phos/ … /3InvdT/)",
              MaterialKind.OLIGO, state=NEXTERA_R1_BK)
    add_input("nextera_ad2_bk", "Nextera-Ad2-bk blocking oligo (/5Phos/ … /3InvdT/)",
              MaterialKind.OLIGO, state=NEXTERA_AD2_BK)
    blocked = mol("washed_beads", "harvested beads washed with blocking oligos")
    add_step("block", "Break droplets; wash beads in 6X SSC with blocking oligos",
             StepKind.HYBRIDIZATION,
             (Input(captured.id, InputRole.TEMPLATE),
              Input("nextera_r1_bk", InputRole.STRUCTURAL),
              Input("nextera_ad2_bk", InputRole.STRUCTURAL)), (blocked,))

    add_input("tso", "Template Switch Oligo (rGrGrG 3′ end)", MaterialKind.OLIGO, state=TSO)
    add_input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    add_input("gap_fill_polymerase", "Hemo KlenTaq", MaterialKind.CATALYST)
    add_input("t7_ligase", "T7 DNA ligase", MaterialKind.CATALYST)
    rna_cdna = mol("barcoded_cdna", "cell-barcoded, UMI-bearing, template-switched cDNA",
                   (Row(chunks=[("SMART handle — cell barcode — UMI — dT — cDNA — TSO", None, False)]),))
    sealed = mol("sealed_atac", "gap-filled, splint-ligated, bead-extended chromatin template",
                 gap_filled_scene())
    add_step("rt_ligation", "One RT/ligation mix on the beads: RT with template switching, "
             "gap fill and splint-side ligation (RT 30 min, then 42 °C 90 min)",
             StepKind.REVERSE_TRANSCRIPTION,
             (Input(blocked.id, InputRole.TEMPLATE),
              Input("tso", InputRole.STRUCTURAL),
              Input("reverse_transcriptase", InputRole.CATALYST),
              Input("gap_fill_polymerase", InputRole.CATALYST),
              Input("t7_ligase", InputRole.CATALYST)), (rna_cdna, sealed))
    add_input("exonuclease_i", "Exonuclease I (37 °C, 45 min)", MaterialKind.CATALYST)
    exo_beads = mol("exo_treated_beads", "beads carrying barcoded cDNA and chromatin templates")
    exo_waste = mol("exo_digested", "single-stranded oligos digested by Exonuclease I",
                    kind=MaterialKind.WASTE)
    add_step("exo_i", "Exonuclease I treatment of the beads", StepKind.DIGESTION,
             (Input(rna_cdna.id, InputRole.TEMPLATE), Input(sealed.id, InputRole.TEMPLATE),
              Input("exonuclease_i", InputRole.CATALYST)), (exo_beads, exo_waste))
    graph.mark_terminal(Terminal(exo_waste.id, Disposition.DISCARD,
                                 "digested by Exonuclease I; which bead-bound species it removes "
                                 "is not stated (Drop-seq: unextended primers)"))

    add_input("tso_pcr", "Tso-PCR", MaterialKind.OLIGO, state=TSO_PCR)
    add_input("nextera_r2", "Nextera-R2", MaterialKind.OLIGO, state=NEXTERA_R2)
    add_input("pcr_polymerase", "PCR polymerase (not named)", MaterialKind.CATALYST)
    preamp = mol("common_preamp", "co-amplified cDNA and chromatin product (10,000 beads)",
                 preamp_product())
    add_step("preamp_pcr", "PCR 10,000 beads with Nextera-R2 / Tso-PCR, 16 cycles total; "
             "column purification", StepKind.PCR,
             (Input(exo_beads.id, InputRole.TEMPLATE),
              Input("tso_pcr", InputRole.STRUCTURAL), Input("nextera_r2", InputRole.STRUCTURAL),
              Input("pcr_polymerase", InputRole.CATALYST)), (preamp,),
             geometry=PCRGeometry.exact(
                 atac_template(), PCRPrimer("Tso-PCR", TSO_PCR),
                 PCRPrimer("Nextera-R2", NEXTERA_R2),
                 statement="Supplementary Table 4 prints both primers; drawn on the chromatin "
                 "product. cDNA amplifies through Tso-PCR sites at both ends."))
    rna_half = mol("rna_half", "half of the PCR product for the cDNA library")
    atac_half = mol("atac_half", "half of the PCR product for the chromatin library")
    add_step("split_halves", "Split the purified PCR product into two halves", StepKind.SPLIT,
             (Input(preamp.id, InputRole.TEMPLATE),), (rna_half, atac_half))

    add_input("spri_06x", "0.6X bead size selection", MaterialKind.SELECTION)
    rna_enriched = mol("rna_enriched", "Tso-PCR-enriched cDNA after 0.6X selections")
    rna_small = mol("rna_small_fragments", "fragments not retained by 0.6X beads",
                    kind=MaterialKind.WASTE)
    add_step("rna_enrich", "0.6X selection, Tso-PCR-alone enrichment PCR, 0.6X selection",
             StepKind.PCR,
             (Input(rna_half.id, InputRole.TEMPLATE), Input("tso_pcr", InputRole.STRUCTURAL),
              Input("pcr_polymerase", InputRole.CATALYST),
              Input("spri_06x", InputRole.SELECTION)), (rna_enriched, rna_small),
             geometry=PCRGeometry.unavailable(
                 RNA, reason="Single-primer Tso-PCR enrichment of full-length cDNA; the cDNA "
                 "state is not drawn base by base.", symmetric_species=True))
    graph.mark_terminal(Terminal(rna_small.id, Disposition.DISCARD,
                                 "removed by the two 0.6X bead size selections"))
    add_input("nextera_xt", "Nextera XT kit (tagmentation + index PCR; primers not printed)",
              MaterialKind.ADAPTER)
    add_input("rna_finish_primers", "Nextera XT PCR primers (P5-side primer and i7 not printed)",
              MaterialKind.OLIGO)
    rna_library = mol("rna_library", "SNARE-seq RNA library", RNA, MaterialKind.LIBRARY)
    add_step("finish_rna", "Nextera XT library “as described previously” (Drop-seq)",
             StepKind.PCR,
             (Input(rna_enriched.id, InputRole.TEMPLATE),
              Input("nextera_xt", InputRole.STRUCTURAL),
              Input("rna_finish_primers", InputRole.STRUCTURAL),
              Input("pcr_polymerase", InputRole.CATALYST)), (rna_library,),
             geometry=PCRGeometry.unavailable(
                 RNA, reason="The SNARE-seq sources do not print the Nextera XT P5-side primer "
                 "or i7 adapter used for cDNA."))

    add_input("p5xx_tso", "P50x-Tso index primer family", MaterialKind.OLIGO, state=P501_TSO)
    add_input("ad2x", "Ad2.x index primer family", MaterialKind.OLIGO, state=AD2_1)
    atac_pcr = mol("atac_indexed", "P5/P7-indexed chromatin amplicons", ATAC)
    add_step("finish_atac", "Index PCR with P5XX-Tso / Ad2.X", StepKind.PCR,
             (Input(atac_half.id, InputRole.TEMPLATE),
              Input("p5xx_tso", InputRole.STRUCTURAL), Input("ad2x", InputRole.STRUCTURAL),
              Input("pcr_polymerase", InputRole.CATALYST)), (atac_pcr,),
             geometry=PCRGeometry.exact(
                 ATAC, PCRPrimer("P501-Tso", P501_TSO), PCRPrimer("Ad2.1_TAAGGCGA", AD2_1),
                 statement="Supplementary Table 4; P501-Tso and Ad2.1 represent index families "
                 "whose per-library members are not stated."))
    add_input("page_gel", "PAGE gel excision of 225–1000 bp, column purification",
              MaterialKind.SELECTION)
    atac_library = mol("atac_library", "SNARE-seq accessibility library", ATAC,
                       MaterialKind.LIBRARY)
    atac_off = mol("atac_out_of_range", "chromatin amplicons outside 225–1000 bp",
                   kind=MaterialKind.WASTE)
    add_step("size_atac", "Excise 225–1000 bp from PAGE gel", StepKind.SELECTION,
             (Input(atac_pcr.id, InputRole.TEMPLATE), Input("page_gel", InputRole.SELECTION)),
             (atac_library, atac_off))
    graph.mark_terminal(Terminal(atac_off.id, Disposition.DISCARD,
                                 "left in the gel outside the excised 225–1000 bp band"))

    graph.mark_terminal(Terminal(rna_library.id, Disposition.SEQUENCING),
                        SequencingPlan(RNA, RP, tuple(p.role for p in RP),
                                       read_lengths={"Read 1": 30, "Read 2": 80},
                                       read_uses={"Read 1": "extract", "Read 2": "include"}))
    graph.mark_terminal(Terminal(atac_library.id, Disposition.SEQUENCING),
                        SequencingPlan(ATAC, AP, tuple(p.role for p in AP),
                                       read_lengths={"Read 1": 30, "Index 1 (i7)": 75,
                                                     "Read 2": 75},
                                       read_uses={"Read 1": "extract", "Index 1 (i7)": "include",
                                                  "Read 2": "include"},
                                       unplaced_reads=(ATAC_INDEX2_READ,)))

    reviewed_targets = (
        "nucleus", "loaded_tn5_adapters", "tde1", "tagment", "capture_splint",
        "barcode_bead_oligo", "capture", "nextera_r1_bk", "nextera_ad2_bk", "block",
        "tso", "reverse_transcriptase", "gap_fill_polymerase", "t7_ligase", "rt_ligation",
        "exonuclease_i", "exo_i", "tso_pcr", "nextera_r2", "pcr_polymerase", "preamp_pcr",
        "split_halves", "spri_06x", "rna_enrich", "nextera_xt", "rna_finish_primers",
        "finish_rna", "p5xx_tso", "ad2x", "finish_atac", "page_gel", "size_atac",
        "rna_library", "atac_library",
    )
    graph.declare_reviewed_scope(reviewed_targets, source)
    return graph.require_publishable()
