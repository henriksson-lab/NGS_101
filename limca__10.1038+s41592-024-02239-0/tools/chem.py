"""Protocol-owned source, molecular workflow and reaction graph."""
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from chemdraw import Construct
from protocol_model import chain, molecule, role, role_library
from extended_protocols import Definition
from chemdraw import feature

DEFINITION = Definition('LiMCA', '10.1038/s41592-024-02239-0', 'Wu et al., Nature Methods (2024)', 'Cytoplasm and nucleus of one picked cell are separated by gentle lysis and centrifugation: the supernatant goes through Smart-seq2, while the nuclear pellet is fixed, NlaIII-digested, proximity-ligated, lysed and amplified with META (20 META tags).', 'Two libraries per cell (Tn5-tagmented Smart-seq2 cDNA and META-amplified 3C DNA), each sequenced paired-end 2 × 150 bp (HiSeq X10 or NovaSeq); they are joined by experimental labels (tube of origin), not by a shared molecular barcode.', '01_limca.html')
# Cells are pooled after library PCR (Methods, Library construction), so each cell's library
# must carry a distinguishing index; its form, length and bases are not printed (inferred role).
CELL_INDEX = role("per-cell library index (not printed)", tag="cbc", inferred=True,
                  molecular_feature=feature("cell_library_index", "sample_index", "unknown"))
PAYLOAD = (role("Smart-seq2 insert or Hi-C junction"),)

TITLE = DEFINITION.title
NOTES = DEFINITION.notes
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.citation}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ("This pass read the Nature Methods article Methods (PMC11166570). They name the Smart-seq2 "
          "oligo-dT30VN, template-switch oligo and ISPCR primer but print no bases, and describe META only "
          "as 'as previously described' with 20 META tags. The META amplification and gDNA library steps, "
          "and the per-cell library indices, are therefore unresolved boundaries; the author protocol "
          "(Research Square pex-2423, ref. 44) was not retrievable.")
FINAL_LIBRARY = role_library(TITLE + " final library roles", PAYLOAD, index=CELL_INDEX)
SEQ_PRIMERS = ()
FINAL_CAPTION = DEFINITION.endpoint
SEQUENCING_INTRO = DEFINITION.endpoint
READ_LENGTHS = None
SEQUENCING_ENDING = DEFINITION.endpoint
SEQUENCING_UNAVAILABLE = "The article states paired-end 2 × 150 bp for both libraries but prints no adapter, index or sequencing-primer sequences (Vazyme TD502 for cDNA; META tags for gDNA), so no primer landing is asserted."

def workflow():
    return chain(molecule("one picked cell", (role("cytoplasmic mRNA"), role("nuclear chromatin"))), (
        ("Soft lysis (0.45% IGEPAL, 30 min on ice) and 500g spin: 5 µl supernatant vs nuclear pellet",
         molecule("matched compartments", (role("cytoplasm"), role("nucleus"))),
         "Physical nucleus–cytoplasm separation, as in Trio-seq; no molecular barcode links the halves."),
        ("Cytoplasm: Smart-seq2 (oligo-dT30VN, SSII + TSO, 21 cycles ISPCR, 0.7× AMPure)",
         molecule("full-length cDNA amplicons", (role("full-length cDNA"),)),
         "Primer bases are not printed in the article."),
        ("Nucleus: 2.5% PFA fix, NlaIII digest, T4 ligation, protease lysis",
         molecule("single-cell contact junctions", (role("restriction fragment A"), role("restriction fragment B")),
                  junctions=(("restriction fragment A", "restriction fragment B", "proximity ligation"),)),
         "No biotin fill-in or pull-down: a conventional 3C ligation inside the nucleus."),
        ("Nextera Tn5 (TD502) on 1–5 ng cDNA; META WGA (20 tags) on the 3C lysate",
         molecule("paired LiMCA libraries", (CELL_INDEX, role("Smart-seq2 insert or Hi-C junction"))),
         "INFERRED — META and gDNA library steps follow cited earlier work; index bases are not printed."),
    ))


def oligos():
    return ()


def _cleanup(graph, prefix, template, label, beads, out_label, state):
    selection = graph.input(prefix + "_beads", beads, MaterialKind.SELECTION)
    kept = graph.product(prefix + "_kept", out_label, MaterialKind.MOLECULAR, state=state)
    waste = graph.product(prefix + "_waste", label + " bead supernatant/unbound material",
                          MaterialKind.WASTE)
    graph.react(prefix, "Purify " + label + " with " + beads, StepKind.CLEANUP,
                ((template, InputRole.TEMPLATE), (selection, InputRole.SELECTION)),
                (kept, waste))
    graph.discard(waste, "not retained by the SPRI/AMPure bead step")
    return kept


def protocol_graph() -> ProtocolGraph:
    d = DEFINITION
    final = FINAL_LIBRARY
    source = f"{d.citation}; doi:{d.doi}; PMC11166570 Methods"
    graph = ProtocolGraph(d.title, authority=source)
    cell = graph.input(
        "cell", "one picked viable cell",
        MaterialKind.BIOLOGICAL,
        state=molecule("one picked cell", (
            role("cytoplasmic mRNA"), role("nuclear chromatin"))), primary=True)
    lysis = graph.input("soft_lysis", "7 µl soft lysis buffer (0.45% IGEPAL-CA630, SUPERaseIn), 30 min ice, vortex, 500g 5 min 4 °C",
                        MaterialKind.CONDITION)
    cytoplasm = graph.product("cytoplasm", "5 µl supernatant (cytoplasmic RNA)",
                              MaterialKind.MOLECULAR, state=cell.state)
    nucleus = graph.product("nucleus", "nuclear pellet in remaining 2 µl",
                            MaterialKind.MOLECULAR, state=cell.state)
    graph.react("separate_compartments", "Soft-lyse and centrifuge: move 5 µl supernatant to a new tube",
                StepKind.SPLIT, ((cell, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)),
                (cytoplasm, nucleus), changes_molecule=False)

    # RNA branch: Smart-seq2 as printed in the Methods (no bases printed).
    dt = graph.input("rna_rt_primer", "oligo-dT30VN (1 µM; full sequence not printed)", MaterialKind.OLIGO)
    rt = graph.input("rna_rt_enzyme", "SuperScript II reverse transcriptase (10 U/µl), betaine, MgCl2, GTP",
                     MaterialKind.CATALYST)
    first = graph.product("rna_first_strand", "oligo-dT-primed first-strand cDNA",
                          MaterialKind.MOLECULAR,
                          state=molecule("first-strand cDNA", (role("oligo-dT30VN"), role("cDNA"))))
    graph.react("rna_rt", "72 °C 5 min with oligo-dT/dNTP; RT 42 °C 90 min, 10 × (50 °C 2 min/42 °C 2 min)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cytoplasm, InputRole.TEMPLATE), (dt, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)),
                (first,))
    tso = graph.input("rna_tso", "template-switch oligonucleotide (1 µM; sequence and modifications not printed)",
                      MaterialKind.OLIGO)
    switched = graph.product("rna_cdna", "full-length cDNA carrying TSO complement",
                             MaterialKind.MOLECULAR,
                             state=molecule("template-switched cDNA", (role("TSO complement"), role("cDNA"), role("oligo-dT30VN"))))
    graph.react("rna_template_switch", "Template switching in the same RT reaction",
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)),
                (switched,))
    ispcr = graph.input("rna_ispcr", "ISPCR primer (0.28 µl of 10 µM; sequence not printed)", MaterialKind.OLIGO)
    kapa = graph.input("rna_kapa", "KAPA HiFi HotStart ReadyMix", MaterialKind.CATALYST)
    amplicons = graph.product("rna_amplicons", "Smart-seq2 cDNA amplicons", MaterialKind.MOLECULAR,
                              state=molecule("full-length cDNA amplicons", (role("full-length cDNA"),)))
    graph.react("rna_preamp", "ISPCR preamplification: 98 °C 3 min, 21 × (98 °C 20 s, 65 °C 30 s, 72 °C 4 min)",
                StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (ispcr, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)),
                (amplicons,), geometry=PCRGeometry.unavailable(
                    Construct([role("TSO-side handle"), role("full-length cDNA"), role("oligo-dT-side handle")],
                              name="Smart-seq2 cDNA roles"),
                    symmetric_species=True,
                    reason="The article names a single ISPCR primer but prints neither its sequence nor the TSO/oligo-dT handles it anneals to."))
    rna_clean = _cleanup(graph, "rna_ampure", amplicons, "cDNA amplicons", "0.7× AMPure XP",
                         "purified cDNA amplicons", amplicons.state)
    tn5 = graph.input("rna_tn5", "Nextera-type Tn5 (Vazyme TD502), 1–5 ng cDNA input", MaterialKind.CATALYST)
    tn5_adapters = graph.input("rna_tn5_adapters", "TD502 loaded Tn5 adapters (sequences not printed)",
                               MaterialKind.ADAPTER, inferred=True)
    tagmented = graph.product("rna_tagmented", "tagmented cDNA fragments", MaterialKind.MOLECULAR,
                              state=molecule("tagmented cDNA", (role("Tn5 adapter"), role("cDNA fragment"), role("Tn5 adapter"))))
    graph.react("rna_tagment", "Tagment 1–5 ng cDNA with Vazyme TD502 Tn5", StepKind.TAGMENTATION,
                ((rna_clean, InputRole.TEMPLATE), (tn5_adapters, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)),
                (tagmented,))
    rna_final = role_library("LiMCA cDNA library roles", (CELL_INDEX, role("cDNA fragment")))
    fwd = graph.input("rna_index_fwd", "TD502 kit index PCR primer (not printed)", MaterialKind.OLIGO, inferred=True)
    rev = graph.input("rna_index_rev", "TD502 kit index PCR primer, second (not printed)", MaterialKind.OLIGO, inferred=True)
    pol = graph.input("rna_index_pol", "TD502 kit PCR polymerase (not named)", MaterialKind.CATALYST, inferred=True)
    rna_lib = graph.product("rna_library", "indexed cDNA library (per cell)", MaterialKind.MOLECULAR,
                            state=rna_final, inferred=True)
    graph.react("rna_index_pcr", "TD502 'library preparation' PCR (conditions not given)", StepKind.PCR,
                ((tagmented, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL), (rev, InputRole.STRUCTURAL),
                 (pol, InputRole.CATALYST)), (rna_lib,),
                geometry=PCRGeometry.unavailable(rna_final, reason="Kit index primers and their sequences are not given in the article."))
    rna_pool = graph.input("rna_other_cells", "libraries from other cells", MaterialKind.MOLECULAR)
    pooled = graph.product("rna_pooled", "pooled cDNA libraries", MaterialKind.MOLECULAR, state=rna_final, inferred=True)
    graph.react("rna_pool", "Pool cells", StepKind.POOL,
                ((rna_lib, InputRole.TEMPLATE), (rna_pool, InputRole.TEMPLATE)), (pooled,), changes_molecule=False)
    first_sel = _cleanup(graph, "rna_spri_06", pooled, "pooled libraries", "0.6× SPRI", "0.6× SPRI-purified pool", rna_final)
    second_sel = graph.input("rna_spri_02_beads", "0.2× SPRI", MaterialKind.SELECTION)
    rna_seq = graph.product("rna_sequencing_library", "LiMCA cDNA sequencing library", MaterialKind.LIBRARY,
                            state=rna_final, inferred=True)
    rna_waste = graph.product("rna_spri_02_waste", "material not kept after 0.2× SPRI", MaterialKind.WASTE)
    graph.react("rna_spri_02", "Then purify with 0.2× SPRI (fraction kept not stated)", StepKind.CLEANUP,
                ((first_sel, InputRole.TEMPLATE), (second_sel, InputRole.SELECTION)), (rna_seq, rna_waste))
    graph.discard(rna_waste, "not retained by the second SPRI step")
    graph.sequence(rna_seq, SequencingPlan(limitation="cDNA library: paired-end 2 × 150 bp on HiSeq X10 or NovaSeq, 0.2–0.6 Gb per cell; adapter, index and primer sequences not printed."))

    # Nuclear branch: single-nucleus 3C.
    pfa = graph.input("formaldehyde", "8 µl 2.5% paraformaldehyde, RT 10 min; quench 0.25 M glycine (+ magnetic beads)",
                      MaterialKind.CONDITION)
    fixed = graph.product("fixed_nucleus", "fixed nucleus", MaterialKind.MOLECULAR, state=nucleus.state)
    graph.react("crosslink", "Fix the nucleus and quench", StepKind.CROSS_LINKING,
                ((nucleus, InputRole.TEMPLATE), (pfa, InputRole.CONDITION)), (fixed,))
    sds = graph.input("sds_triton", "Hi-C lysis buffer, 0.3% SDS 62 °C 10 min, Triton X-100 quench 37 °C 15 min",
                      MaterialKind.CONDITION)
    open_nucleus = graph.product("permeabilized_nucleus", "SDS-opened fixed nucleus", MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("permeabilize", "Lyse and open chromatin with SDS, quench with Triton", StepKind.PERMEABILIZATION,
                ((fixed, InputRole.TEMPLATE), (sds, InputRole.CONDITION)), (open_nucleus,))
    nla = graph.input("restriction_enzyme", "NlaIII (4 U/µl in rCutSmart; printed 'NIAIII')", MaterialKind.CATALYST)
    digested = graph.product("restriction_fragments", "NlaIII-digested fixed chromatin",
                             MaterialKind.MOLECULAR, state=fixed.state)
    graph.react("digest", "Digest with NlaIII, 37 °C 2 h rotating", StepKind.DIGESTION,
                ((open_nucleus, InputRole.TEMPLATE), (nla, InputRole.CATALYST)), (digested,))
    ligase = graph.input("proximity_ligase", "T4 DNA ligase (10 U/µl in 1× T4 buffer)", MaterialKind.CATALYST)
    contacts = graph.product(
        "contacts", "in-nucleus proximity-ligation junctions",
        MaterialKind.MOLECULAR,
        state=molecule("single-cell contact junctions", (
            role("restriction fragment A"), role("restriction fragment B")),
            junctions=(("restriction fragment A", "restriction fragment B",
                        "proximity ligation"),)))
    graph.react("proximity_ligation", "Buffer exchange, then T4 ligation RT 2.5 h rotating",
                StepKind.PROXIMITY_LIGATION,
                ((digested, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (contacts,))
    reversal = graph.input("crosslink_reversal", "2 µl lysis buffer with Qiagen protease; 50 °C 1 h, 65 °C 1 h, 70 °C 15 min",
                           MaterialKind.CONDITION)
    lysate = graph.product("contact_lysate", "protease-lysed single-cell 3C DNA (no purification described)",
                           MaterialKind.MOLECULAR, state=contacts.state)
    graph.react("reverse_crosslinks", "Protease lysis and heat; store at −80 °C",
                StepKind.DEPROTEINIZATION,
                ((contacts, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (lysate,))
    meta_tn5 = graph.input("meta_tn5", "Tn5 transposase (Vazyme S111-01)", MaterialKind.CATALYST)
    meta_tags = graph.input("meta_transposome", "20 META tags loaded on Tn5 (sequences not printed)",
                            MaterialKind.ADAPTER, inferred=True)
    tagged = graph.product("meta_tagged", "META-transposed 3C DNA", MaterialKind.MOLECULAR,
                           state=molecule("META-tagged 3C DNA", (role("META tag"), role("restriction fragment A"),
                                                                 role("restriction fragment B"), role("META tag"))))
    graph.react("meta_transposition", "META: Tn5 transposition with 20 META tags", StepKind.TAGMENTATION,
                ((lysate, InputRole.TEMPLATE), (meta_tags, InputRole.STRUCTURAL), (meta_tn5, InputRole.CATALYST)),
                (tagged,))
    hic_final = role_library("LiMCA 3C (gDNA) library roles", (CELL_INDEX, role("3C fragment or contact junction")))
    meta_rest = graph.input("meta_amplification", "META amplification and gDNA library construction 'as previously described' (refs 18, 26)",
                            MaterialKind.CONDITION)
    hic_seq = graph.product("hic_sequencing_library", "LiMCA 3C (gDNA) sequencing library", MaterialKind.LIBRARY,
                            state=hic_final, inferred=True)
    graph.react("meta_amplify_library", "META amplification and gDNA library construction", StepKind.MAGIC,
                ((tagged, InputRole.TEMPLATE), (meta_rest, InputRole.CONDITION)), (hic_seq,),
                reason="The LiMCA Methods defer META amplification and gDNA library preparation to refs 18/26 and the unretrieved author protocol; primers, cycles and indexing are not reconstructed.")
    graph.sequence(hic_seq, SequencingPlan(limitation="gDNA library: paired-end 2 × 150 bp on HiSeq X10 or NovaSeq, 3–6 Gb per cell; META tag, index and primer sequences not printed."))
    reviewed = [x for x in graph.materials] + [x for x in graph.steps]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
