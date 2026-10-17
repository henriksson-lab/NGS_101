"""Protocol-owned source, molecular workflow and reaction graph.

Source scope (round 2026-10-17-breadth-first-11): published NAR article PMC11724362
Methods "Single-cell encapsulation and mRNA barcoding using inDrops platform",
"inDrops-2 (IVT) sequencing library preparation", "inDrops-2 (TS) sequencing library
preparation" and "Barcoded hydrogel bead synthesis" (JATS .txt lines 105-125), and the
corrected Supplementary Tables workbook, Table S1 (xlsx.txt lines 16-55), S2/S4 heads.
The two branches use different beads (v1 vs v2) and different library construction, so
they split before encapsulation rather than after a shared RT.
"""
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct
from pcr import PCRGeometry
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind, Terminal)
from protocol_model import chain, molecule, role
from extended_protocols import Definition, CELL, UMI, ident

DEFINITION = Definition(
    'inDrops-2', '10.1093/nar/gkae1312',
    'Juzėnas et al., Nucleic Acids Research (2025)',
    'Open-source droplet capture with two separately validated library routes: v1 T7-promoter '
    'beads followed by IVT linear amplification, or v2 beads with in-droplet template switching '
    'followed by cDNA PCR, enzymatic fragmentation and adapter ligation.',
    'Read 1 reads the bead cell barcode and UMI; Read 2 reports the transcript insert; '
    'i7 reads the PE2 sample index.',
    '01_indrops_2.html')

SAMPLE = ident("sample index", "sample_index", "sample_index", length=6)
TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = (f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>, '
          'Methods and corrected Supplementary Table S1.')
SUMMARY = DEFINITION.summary
CAVEAT = (
    'SOURCE DISCREPANCY — v1 (IVT) bead primer: the article Methods ("Synthesis of V1 beads") prints '
    'the full-length primer as /5Acryd/PC/CGATGACG…ATACCACCATGG…CCGATCT[11 positions]'
    'GAGTGATTGCTTGTGACGCCTT[8]NNNNNNNN(T)19, whereas corrected Table S1 prints '
    'CGATTGATCAACG…[10 positions]…[8]NNNNNN(T)19V; Table S2 bc2 primers carry an N8 window '
    'and (A)19, and Table S2 bc1 barcodes vary in length. The model draws only role-level '
    'barcode/UMI segments and does not choose a UMI length or 5′ stub. '
    'SOURCE DISCREPANCY — v2 (TS) bead stub: Methods prints /5Acryd/CGATGACG(PC)CTACACGACGCTCTTC-3′, '
    'while Table S1 lists a "V4 DNA stub" ending …CTCTTCCGATCT; the Methods full-length v2 primer '
    'ends in (T)18, Table S1 in (T)19. The model uses neither for drawn bases. '
    'The TS ligation adapter is shown as a reagent; its annealed end geometry is not drawn.')

FINAL_LIBRARY = Construct((
    role("PE1 / Illumina Read 1 end"), CELL, UMI, role("oligo(dT)"),
    role("3-prime transcript insert"), role("PE2 / Illumina Read 2 site"), SAMPLE,
    role("P7 end")), name=TITLE + " final library roles")
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = (
    'Reported acquisition: IVT branch R1 54, i7 8, R2 ≥35 cycles (NextSeq550/HiSeq2500); '
    'TS branch R1 28, i7 8, R2 35–92 cycles (MiSeq/HiSeq2500/NextSeq550/NovaSeq6000). '
    'The printed PE2 indexes are 6 nt, so the 8-cycle i7 read extends past the index.')
SEQUENCING_UNAVAILABLE = (
    "The article does not print the sequencing primers; this role-level model does not place them "
    "on the PE1/PE2 sites.")

_IVT_SEQ = ("IVT branch: R1 54 cycles (cell barcode + UMI), i7 8 cycles, R2 ≥35 cycles "
            "(insert), per Methods 'inDrops-2 (IVT) sequencing library preparation'.")
_TS_SEQ = ("TS branch: R1 28 cycles (8+4+8 cell barcode/linker + 8 UMI), i7 8 cycles, R2 35–92 "
           "cycles (insert), per Methods 'inDrops-2 (TS) sequencing library preparation'.")


def oligos():
    """Table S1 (corrected workbook) strings, verbatim apart from removed spaces."""
    return (
        "<b>IVT branch</b>",
        "V1 DNA stub primer 5′-/5Acryd/PC/CGATTGATCAACGTAATACGACTCACTATAGGGATACCATCTACACTCTTTCCCTACACGACGCTCTTCCG-3′",
        "PE2-N6 (2nd step RT primer) 5′-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNN-3′",
        "PE1 (PCR primer) 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGA-3′",
        "PE2-Ind1 (indexing primer) 5′-CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC-3′ (Ind2–6 differ in the 6-nt index)",
        "<b>Template-switching branch</b>",
        "Full length RT primer (v2 design) 5′-CTACACGACGCTCTTCCGATCT[12345678]CATG[12345678]NNNNNNNNTTTTTTTTTTTTTTTTTTT-3′",
        "TSO 5′-AAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3′",
        "cDNA FWD primer 5′-CTACACGACGCTCTTCCGATCT-3′",
        "cDNA REV primer 5′-AAGCAGTGGTATCAACGCAGAGT-3′",
        "Ligation FWD primer 5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCAC/3ddC/-3′",
        "Ligation REV primer 5′-/5AmMC6/GCTCTTCCGATCT-3′ (annealed with Ligation FWD)",
        "Forward PCR index primer 5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′",
        "PE2-ind1 5′-CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTG*T-3′ (* phosphorothioate; ind2–6 differ in the 6-nt index)",
    )


def workflow():
    return chain(molecule("cell and barcoded hydrogel bead in a droplet", (role("poly(A) RNA"), CELL, UMI)), (
        ("Photo-release the bead primer, lyse and reverse-transcribe (TSO present only in the TS branch)",
         molecule("barcoded first-strand cDNA", (CELL, UMI, role("oligo(dT)"), role("3-prime transcript cDNA"))),
         "Cell barcode and UMI are copied into every captured molecule."),
        ("IVT: second strand, T7 IVT, RNA fragmentation, PE2-N6 RT; TS: cDNA PCR, FS fragmentation, adapter ligation",
         molecule("adapter-bearing inDrops-2 fragment", (CELL, UMI, role("3-prime transcript insert"), role("PE2 site"))),
         "Two different routes reach the same read roles."),
        ("Index PCR",
         molecule("inDrops-2 library", (CELL, UMI, role("3-prime transcript insert"), SAMPLE)),
         "PE1/PE2-Ind (IVT) or Forward PCR index/PE2-ind (TS) primers."),
    ))


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}; Methods and corrected Table S1"
    graph = ProtocolGraph(d.title, authority=source)
    reviewed: list[str] = []

    def inp(ident_, label, kind):
        reviewed.append(ident_)
        return graph.input(ident_, label, kind)

    def step(ident_, label, kind, inputs, outputs, **kw):
        reviewed.append(ident_)
        graph.react(ident_, label, kind, inputs, outputs, **kw)

    def prod(ident_, label, parts, kind=MaterialKind.MOLECULAR):
        return graph.product(ident_, label, kind, state=molecule(label, parts))

    cells = graph.input(
        "cells", "single-cell suspension", MaterialKind.BIOLOGICAL,
        state=molecule("single-cell suspension", (role("poly(A) RNA"),)), primary=True)
    reviewed.append("cells")
    ivt_cells = graph.product("ivt_cells", "cells loaded with v1 beads (IVT branch)",
                              MaterialKind.BIOLOGICAL, state=cells.state)
    ts_cells = graph.product("ts_cells", "cells loaded with v2 beads (TS branch)",
                             MaterialKind.BIOLOGICAL, state=cells.state)
    step("choose_branch", "Choose the IVT (v1 beads) or TS (v2 beads + TSO) implementation",
         StepKind.SPLIT, ((cells, InputRole.TEMPLATE),), (ivt_cells, ts_cells),
         changes_molecule=False)

    first_parts = (CELL, UMI, role("oligo(dT)"), role("3-prime transcript cDNA"))

    # IVT branch -------------------------------------------------------------
    v1 = inp("v1_bead_primer", "photo-released v1 bead primer (T7 promoter, PE1 Read 1 site, "
             "two-part cell barcode, UMI, oligo(dT))", MaterialKind.OLIGO)
    maxima = inp("ivt_rt_enzyme", "Maxima H minus reverse transcriptase", MaterialKind.CATALYST)
    ivt_first = prod("ivt_first_strand", "v1-barcoded first-strand cDNA (T7 promoter 5′)",
                     (role("T7 promoter"), role("PE1 Read 1 site"), *first_parts))
    step("ivt_rt", "Encapsulate, photo-release (350 nm), lyse and reverse-transcribe in droplets",
         StepKind.REVERSE_TRANSCRIPTION,
         ((ivt_cells, InputRole.TEMPLATE), (v1, InputRole.STRUCTURAL),
          (maxima, InputRole.CATALYST)), (ivt_first,))
    sss = inp("ivt_second_strand_module", "NEBNext Ultra II non-directional second strand synthesis module",
              MaterialKind.CATALYST)
    ivt_ds = prod("ivt_ds_cdna", "double-stranded T7-promoter cDNA",
                  (role("T7 promoter"), role("PE1 Read 1 site"), *first_parts))
    step("ivt_second_strand", "Second-strand synthesis (2.5 h 16 °C, 20 min 65 °C)",
         StepKind.EXTENSION, ((ivt_first, InputRole.TEMPLATE), (sss, InputRole.CATALYST)),
         (ivt_ds,))
    t7 = inp("ivt_t7_polymerase", "HiScribe T7 High Yield RNA synthesis kit", MaterialKind.CATALYST)
    arna = prod("ivt_rna", "amplified antisense RNA",
                (role("PE1 Read 1 site"), CELL, UMI, role("oligo(dT)"), role("transcript RNA copy")))
    step("ivt", "Linearly amplify by T7 in vitro transcription (15 h, 37 °C)",
         StepKind.TRANSCRIPTION, ((ivt_ds, InputRole.TEMPLATE), (t7, InputRole.CATALYST)),
         (arna,))
    frag = inp("ivt_rna_fragmentation", "Ambion RNA fragmentation reagent (70 °C, 1.5 min)",
               MaterialKind.CONDITION)
    arna_frag = prod("ivt_rna_fragments", "barcode-bearing amplified RNA fragments",
                     (role("PE1 Read 1 site"), CELL, UMI, role("oligo(dT)"), role("3-prime transcript fragment")))
    other_frag = graph.product("ivt_internal_fragments", "internal RNA fragments lacking the barcode end",
                               MaterialKind.MOLECULAR,
                               state=molecule("internal RNA fragment", (role("transcript RNA fragment"),)))
    step("ivt_fragment", "Fragment the amplified RNA and stop with EDTA/SPRI",
         StepKind.FRAGMENTATION, ((arna, InputRole.TEMPLATE), (frag, InputRole.CONDITION)),
         (arna_frag, other_frag))
    graph.mark_terminal(Terminal(other_frag.id, Disposition.UNRESOLVED, (
        "Not separated in the source: internal fragments can be PE2-N6-primed but lack the PE1 "
        "end for indexing PCR, so they are not expected to amplify into the library.")))
    pe2n6 = inp("pe2_n6", "PE2-N6 random-hexamer RT primer", MaterialKind.OLIGO)
    maxima2 = inp("ivt_second_rt_enzyme", "Maxima H minus reverse transcriptase", MaterialKind.CATALYST)
    ivt_cdna = prod("ivt_cdna", "PE2-tagged second-step cDNA",
                    (role("PE1 Read 1 site"), CELL, UMI, role("oligo(dT)"),
                     role("3-prime transcript insert"), role("PE2 / Illumina Read 2 site")))
    step("ivt_second_rt", "Random-prime the RNA fragments with PE2-N6 and reverse-transcribe",
         StepKind.REVERSE_TRANSCRIPTION,
         ((arna_frag, InputRole.TEMPLATE), (pe2n6, InputRole.STRUCTURAL),
          (maxima2, InputRole.CATALYST)), (ivt_cdna,))
    pe1 = inp("ivt_pe1", "PE1 PCR primer", MaterialKind.OLIGO)
    pe2ind = inp("ivt_pe2_ind", "PE2-Ind1–6 indexing primer", MaterialKind.OLIGO)
    kapa = inp("ivt_kapa", "KAPA HiFi HotStart PCR mix", MaterialKind.CATALYST)
    ivt_lib = graph.product("ivt_library", "inDrops-2 (IVT) sequencing library",
                            MaterialKind.LIBRARY, state=final)
    step("ivt_pcr", "Index PCR (10 cycles), then 0.6–0.8× SPRI dual-size selection",
         StepKind.PCR,
         ((ivt_cdna, InputRole.TEMPLATE), (pe1, InputRole.STRUCTURAL),
          (pe2ind, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (ivt_lib,),
         geometry=PCRGeometry.unavailable(
             final, reason="Table S1 prints PE1 and PE2-Ind sequences; this role-level model "
             "does not place them on drawn bases."))
    graph.sequence(ivt_lib, SequencingPlan(limitation=_IVT_SEQ))
    reviewed.append("ivt_library")

    # Template-switching branch ---------------------------------------------
    v2 = inp("v2_bead_primer", "photo-released v2 bead primer (Read 1 site, 8+CATG+8 cell barcode, "
             "N8 UMI, oligo(dT))", MaterialKind.OLIGO)
    ts_rt_enzyme = inp("ts_rt_enzyme", "Maxima H minus reverse transcriptase", MaterialKind.CATALYST)
    ts_first = prod("ts_first_strand", "v2-barcoded first-strand cDNA",
                    (role("Read 1 site"), *first_parts))
    step("ts_rt", "Encapsulate, photo-release, lyse and reverse-transcribe in droplets",
         StepKind.REVERSE_TRANSCRIPTION,
         ((ts_cells, InputRole.TEMPLATE), (v2, InputRole.STRUCTURAL),
          (ts_rt_enzyme, InputRole.CATALYST)), (ts_first,))
    tso = inp("tso", "TSO AAGCAGTGGTATCAACGCAGAGTACATrGrGrG (50 µM in RT-lysis mix)",
              MaterialKind.OLIGO)
    switched = prod("ts_switched", "template-switched barcoded cDNA",
                    (role("Read 1 site"), *first_parts, role("TSO complement")))
    step("ts_template_switch", "Template switching in the same droplet RT reaction",
         StepKind.TEMPLATE_SWITCHING,
         ((ts_first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
          (ts_rt_enzyme, InputRole.CATALYST)), (switched,))
    fwd = inp("ts_cdna_fwd", "cDNA FWD primer", MaterialKind.OLIGO)
    rev = inp("ts_cdna_rev", "cDNA REV primer", MaterialKind.OLIGO)
    kapa2 = inp("ts_kapa", "2× KAPA HiFi HotStart ReadyMix", MaterialKind.CATALYST)
    amplified = prod("ts_amplified", "PCR-amplified full-length barcoded cDNA",
                     (role("Read 1 site"), *first_parts, role("TSO complement")))
    step("ts_cdna_pcr", "Amplify cDNA (14 cycles), 0.8× SPRI",
         StepKind.PCR,
         ((switched, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
          (rev, InputRole.STRUCTURAL), (kapa2, InputRole.CATALYST)), (amplified,),
         geometry=PCRGeometry.unavailable(
             amplified.state, reason="Table S1 prints cDNA FWD/REV; their landing is not drawn "
             "in this role-level model."))
    fs = inp("ts_fs_enzyme", "NEB FS DNA Library Prep Kit fragmentation/end-prep enzymes",
             MaterialKind.CATALYST)
    fragments = prod("ts_fragments", "end-repaired, dA-tailed barcode-end fragment",
                     (role("Read 1 site"), CELL, UMI, role("oligo(dT)"),
                      role("3-prime transcript insert"), role("3-prime dA")))
    internal = graph.product("ts_internal_fragments", "internal or TSO-end fragments",
                             MaterialKind.MOLECULAR,
                             state=molecule("internal fragment", (role("cDNA fragment"),)))
    step("ts_fragment", "Fragment 50 ng cDNA (8 min 37 °C), end-repair and dA-tail (30 min 65 °C)",
         StepKind.FRAGMENTATION, ((amplified, InputRole.TEMPLATE), (fs, InputRole.CATALYST)),
         (fragments, internal))
    adapter = inp("ts_ligation_adapter", "annealed Ligation FWD/REV adapter (55 nM)",
                  MaterialKind.ADAPTER)
    ligase = inp("ts_ligase", "NEB FS kit ligation reagents", MaterialKind.CATALYST)
    ligated = prod("ts_ligated", "adapter-ligated barcode-end fragment",
                   (role("Read 1 site"), CELL, UMI, role("oligo(dT)"),
                    role("3-prime transcript insert"), role("PE2 / Illumina Read 2 site")))
    step("ts_ligate", "Ligate the modified adapter (15 min, 20 °C), 0.8× SPRI",
         StepKind.LIGATION,
         ((fragments, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
          (ligase, InputRole.CATALYST)), (ligated,),
         geometry=AdapterLigationGeometry.unavailable(
             "Table S1 prints both adapter strands (5′-phosphate/3′-ddC FWD; 5′-amino REV) and says "
             "they are annealed; the resulting end geometry is not drawn in this role-level model."))
    graph.mark_terminal(Terminal(internal.id, Disposition.UNRESOLVED, (
        "Not separated in the source: fragments without the bead Read 1 end can ligate the adapter "
        "but lack the Forward PCR index primer site, so they are not expected to be amplified.")))
    fidx = inp("ts_forward_index", "Forward PCR index primer", MaterialKind.OLIGO)
    ridx = inp("ts_pe2_ind", "PE2-ind1–6 indexing primer", MaterialKind.OLIGO)
    q5 = inp("ts_q5", "NEBNext Ultra II Q5 Master Mix", MaterialKind.CATALYST)
    ts_lib = graph.product("ts_library", "inDrops-2 (TS) sequencing library",
                           MaterialKind.LIBRARY, state=final)
    step("ts_pcr", "Index PCR (14 cycles), then 0.6–0.8× SPRI dual-size selection",
         StepKind.PCR,
         ((ligated, InputRole.TEMPLATE), (fidx, InputRole.STRUCTURAL),
          (ridx, InputRole.STRUCTURAL), (q5, InputRole.CATALYST)), (ts_lib,),
         geometry=PCRGeometry.unavailable(
             final, reason="Table S1 prints both indexing primers; this role-level model does "
             "not place them on drawn bases."))
    graph.sequence(ts_lib, SequencingPlan(limitation=_TS_SEQ))
    reviewed.append("ts_library")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
