"""Source-owned protocol definition and reaction graph.

Authority: Cho et al., Cell 184:3559 (2021), PMC8238917 JATS (published author
manuscript), STAR Methods lines 211-285 (seed library, 1st-Seq, HDMI-array
processing, tissue, RT, second strand, library PCR) and analysis lines 341-343.
"""
from __future__ import annotations

from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind, Terminal, Disposition)
from chemdraw import Construct, Scene, Segment, complement_segments, feature
from protocol_model import role
import seqscope_chem as sc
from seqscope_chem import oligos as seqscope_oligos

SOURCE_LABEL = "Cho et al., Cell (2021)"
DOI = "10.1016/j.cell.2021.05.010"
AUTHORITY = f"{SOURCE_LABEL}, PMC8238917 STAR Methods; doi:{DOI}"

SEQUENCING = (
    "1st-Seq: the HDMI-DraI seed is clustered on a MiSeq v3 flow cell and read single-end "
    "with the custom Read1-DraI primer (25 cycles for HDMI-DraI, 37 for HDMI32-DraI); the run "
    "stops after Read 1 and the flow cell is retrieved. 2nd-Seq: the PCR library is read "
    "paired-end (100–150 bp; Illumina and BGI platforms) with the TruSeq primers. Read 1 "
    "gives the HDMI (the first 20 bases for HDMI-DraI, 30 for HDMI32-DraI, matched to the "
    "reverse complement of the 1st-Seq read); the first 9 bases of Read 2 are the random-"
    "primer UMI, followed by transcript sequence. The i7 read length is not stated.")

FINAL_LIBRARY, SEQ_PRIMERS = sc.final_library()

TITLE = 'Seq-Scope — flow-cell spatial transcriptomics'
NOTES = '01_seq_scope.html'
SOURCE = AUTHORITY
SUMMARY = 'Selected published HDMI-DraI20 liver/colon route: decode DNA clusters as a spatial array, capture RNA, retain original first-strand DNA on the array and recover independently copied second-strand DNA for indexed library PCR.'
CAVEAT = 'On-instrument bridge amplification/linearization remains an explicit gap. Native RNA, anchored DNA, SBS daughters, copied second strands and original primers remain separate populations. The final library is derived from printed primers, not used as the incoming PCR template. HDMI32 and EcoRI/UMI-array alternatives are not expanded in this graph. RNA hydrolysis and DNA exonuclease digestion do not themselves constitute physical removal.'
FINAL_CAPTION = 'Printed-oligo-derived final library reference: HDMI20 is spatial identity, N9 is random-primer-derived UMI, and the 8-mer i7 is a sample index. Original primer phosphorothioate bonds stay on native primer lineages; copying their bases does not copy their chemical bonds.'
SEQUENCING_INTRO = SEQUENCING
SEQUENCING_ENDING = 'Bottom-tile first-run HDMI reads are reverse-complemented into the whitelist; destroyed top tiles are physically lost at disassembly. Copying Read2 N9 into Read1, error-correction/deduplication and polyA trimming are computational operations. Dataset-specific2nd-Seq cycles and i7 acquisition length remain unassigned;100–150bp PE is the source range.'
oligos = seqscope_oligos


def seed_parts():
    return [Segment('source seed P7',sc.SEED_5P),Segment('source seed TR1',sc.TR1),
            sc._hdmi('seqscope_hdmi20',sc.HDMI20),Segment('source HR1',sc.HR1),
            Segment('capture dT30',sc.POLYT+sc.DRAI[:3])]


def native(stage):
    scene=Scene()
    if stage=='seed':
        scene.strand('seed',[*seed_parts(),Segment('remaining DraI half',sc.DRAI[3:]),Segment('DraI adapter',sc.DRAI_ADAPTER),Segment('source P5rc',sc.SEED_3P)],label='native IDT ssDNA seed; not completed transcript library')
    elif stage in ('clusters','decoded','cut','trimmed','array','captured','hybrid','cleared','first','second'):
        parts=seed_parts()
        if stage in ('clusters','decoded'):parts += [Segment('remaining DraI half',sc.DRAI[3:]),Segment('DraI adapter',sc.DRAI_ADAPTER),Segment('source P5rc',sc.SEED_3P)]
        if stage in ('hybrid','cleared','first','second'):parts += [role('DNA extension copied from captured RNA; length/sequence variable')]
        scene.strand('anchored_DNA',parts,label='P7-bound cluster DNA; original capture lineage remains attached')
        if stage in ('clusters','decoded'):
            scene.strand('other_cluster_lineages',[role('possible native seed/other cluster strands; attachment, linearization and recovery unresolved')])
        if stage in ('cut','trimmed'):
            scene.strand('cleaved_ends',[role('released DraI-adapter/P5rc ends plus possible ExoI hydrolysates until actual wash')])
        if stage in ('decoded','cut','trimmed'):
            scene.strand('SBS_daughter',[role('Read1-DraI original primer plus25 copied DNA bases; retained duplex subset after DraI' if stage!='decoded' else 'Read1-DraI primer plus25 copied DNA bases; opposite to anchored seed')],label='independent SBS extension DNA; first-run HDMI read is reverse complement of spatial DNA')
        if stage in ('captured','hybrid','cleared'):
            scene.strand('RNA',[role('original tissue mRNA body; RNA nucleotides'),role('RNA polyA hybridized to capture dT')],label='independent RNA template; not covalently joined to array DNA')
        if stage=='second':
            scene.strand('second_DNA',[*complement_segments(list(sc.pcr1_product())),role("possible further copied seed 5-prime region; 3-prime extension extent unassigned")],label='PCR-eligible N9-primed DNA daughter reference; original surface DNA remains separate')
        scene.strand('residuals',[role('possible unextended DNA, unbound RNA/oligos or hydrolysis products; class-specific recovery unresolved')])
    elif stage=='rna':
        scene.strand('RNA',[role('tissue mRNA'),role('polyA RNA')],label='released tissue RNA plus possible other tissue polymers')
    elif stage=='eluted':
        scene.strand('eluted_DNA',[*complement_segments(list(sc.pcr1_product())),role("possible further copied seed 5-prime region; 3-prime extension extent unassigned")],label='eluted single-stranded daughter reference; further extension extent and recovery unresolved')
        scene.strand('residuals',[role('possible detached DNA/free oligos; class-specific recovery unresolved')])
    elif stage in ('pcr1','pcr2'):
        con=sc.pcr1_product() if stage=='pcr1' else FINAL_LIBRARY
        # Native product classes remain separate; final canonical read reference is not a PCR input.
        scene.strand('daughter_one',list(con),label='amplified DNA strand; copied bases do not inherit original primer PS bonds')
        scene.strand('daughter_two',complement_segments(list(con)),label='independent complementary amplified DNA')
        scene.strand('originals',[role('possible original incoming DNA, originalPS-primers and off-target products; class recoveries unresolved')])
        scene.strand('PS_incorporated_daughters',[role('DNA molecules directly incorporating original PCR primers retain their two3′ phosphorothioate bonds; copied lineages do not inherit those bonds')])
    return scene


def _strip(seq):
    return seq.replace("*", "")


def _cleanup(graph, reviewed, ident, label, template, selection_label, kept_label,
             waste_label, waste_reason, state=None, kind=StepKind.CLEANUP,
             material=MaterialKind.MOLECULAR):
    matrix = graph.input(ident + "_matrix", selection_label, MaterialKind.SELECTION)
    kept = graph.product(ident + "_kept", kept_label, material,
                         state=state if state is not None else template.state)
    waste = graph.product(ident + "_waste", waste_label, MaterialKind.WASTE)
    graph.react(ident, label, kind,
                ((template, InputRole.TEMPLATE), (matrix, InputRole.SELECTION)),
                (kept, waste))
    graph.discard(waste, waste_reason)
    reviewed.extend((ident + "_matrix", ident))
    return kept


def _local_graph():
    graph = ProtocolGraph(TITLE, authority=AUTHORITY)
    r = []
    exact = sc.graph_states()
    seed = graph.input("seed", "HDMI-DraI ssDNA seed library (IDT Ultramer)",
                       MaterialKind.OLIGO, state=native("seed"), primary=True)
    flow_cell = graph.input("flow_cell", "MiSeq v3 flow cell and on-instrument cluster reagents",
                            MaterialKind.CONDITION)
    clusters = graph.product("clusters", "HDMI clusters on the MiSeq flow cell",
                             MaterialKind.MOLECULAR,
                             state=native("clusters"))
    graph.react("cluster", "Load 100 pM HDMI-DraI and generate clusters on the MiSeq",
                StepKind.MAGIC, ((seed, InputRole.TEMPLATE), (flow_cell, InputRole.CONDITION)),
                (clusters,),
                reason="Native bridge amplification and strand linearisation remain unassigned within the inspected source scope; no platform recipe is supplied.")
    read1_drai = graph.input("read1_drai", "Read1-DraI custom Read 1 primer", MaterialKind.OLIGO,
                             state=exact["read1_drai"])
    sbs = graph.input("first_seq_sbs", "MiSeq v3 SBS polymerase and reagents (manual mode, 25 cycles single end, no re-synthesis)",
                      MaterialKind.CATALYST)
    decoded = graph.product("decoded", "coordinate-decoded clusters with Read1-DraI duplex",
                            MaterialKind.MOLECULAR,
                            state=native("decoded"))
    graph.react("first_seq", "1st-Seq: read the HDMI with Read1-DraI and retrieve the flow cell",
                StepKind.EXTENSION,
                ((clusters, InputRole.TEMPLATE), (read1_drai, InputRole.STRUCTURAL),
                 (sbs, InputRole.CATALYST)), (decoded,))
    drai = graph.input("drai", "DraI (NEB R0129), 1× CutSmart, 37 °C overnight", MaterialKind.CATALYST)
    cut = graph.product("drai_cut", "clusters with exposed 3′ oligo-dT", MaterialKind.MOLECULAR,
                        state=native("cut"))
    graph.react("drai_digest", "DraI cuts TTT|AAA to remove the P5 end and expose oligo-dT",
                StepKind.DIGESTION,
                ((decoded, InputRole.TEMPLATE), (drai, InputRole.CATALYST)), (cut,))
    exo1 = graph.input("exo1_array", "Exonuclease I (NEB M2903), 37 °C 45 min", MaterialKind.CATALYST)
    trimmed = graph.product("exo1_array_kept", "Read1-DraI-protected HDMI strands",
                            MaterialKind.MOLECULAR, state=native("trimmed"))
    graph.react("exo1_array_digest", "Exo I removes non-specific ssDNA", StepKind.DIGESTION,
                ((cut, InputRole.TEMPLATE), (exo1, InputRole.CATALYST)), (trimmed,))
    naoh = graph.input("naoh_array", "water, 3× 0.1 N NaOH 5 min, 0.1 M Tris pH 7.5, water washes",
                       MaterialKind.CONDITION)
    array = graph.product("hdmi_array", "single-stranded HDMI-array (oligo-dT 3′ ends)",
                          MaterialKind.MOLECULAR, state=native("array"))
    primer_waste = graph.product("read1_drai_removed", "wash fractions: SBS primer/daughters, cleaved ends and hydrolysates",
                                 MaterialKind.WASTE)
    graph.react("denature_array", "Denature and wash away the Read1-DraI strand",
                StepKind.STRAND_SEPARATION,
                ((trimmed, InputRole.TEMPLATE), (naoh, InputRole.CONDITION)), (array, primer_waste))
    graph.discard(primer_waste, "Physically removed by water/NaOH/Tris washes; digestion alone was not a discard.")
    disassembly=graph.input("disassembly_condition","source flow-cell disassembly and water removal of glass/debris",MaterialKind.CONDITION)
    bottom=graph.product("bottom_array","retained bottom HDMI tiles",MaterialKind.MOLECULAR,state=array.state)
    top=graph.product("top_array_lost","destroyed top tiles plus glass/debris",MaterialKind.WASTE)
    graph.react("disassemble","Expose bottom array; top tiles physically destroyed",StepKind.SPLIT,((array,InputRole.TEMPLATE),(disassembly,InputRole.CONDITION)),(bottom,top))
    graph.discard(top,"Source states top tiles destroyed at disassembly; first-run coordinate whitelist subsequently uses bottom tiles only.")
    array=bottom
    tissue = graph.input("tissue", "10 µm fresh-frozen section on the disassembled HDMI-array, "
                         "4% formaldehyde, H&E imaging", MaterialKind.BIOLOGICAL)
    release = graph.input("rna_release", "collagenase I 37 °C 20 min, then 1 mg/mL pepsin in 0.1 M HCl 37 °C 10 min",
                          MaterialKind.CONDITION)
    rna = graph.product("released_rna", "poly(A) RNA released from the fixed section",
                        MaterialKind.MOLECULAR, state=native("rna"))
    graph.react("release_rna", "Release RNA from the fixed tissue over the array",
                StepKind.PERMEABILIZATION,
                ((tissue, InputRole.TEMPLATE), (release, InputRole.CONDITION)), (rna,))
    captured = graph.product("captured", "tissue poly(A) RNA hybridized to oligo-dT",
                             MaterialKind.MOLECULAR,
                             state=native("captured"))
    graph.react("capture", "Capture released RNA on the HDMI-array oligo-dT",
                StepKind.HYBRIDIZATION,
                ((array, InputRole.TEMPLATE), (rna, InputRole.STRUCTURAL)), (captured,))
    rt = graph.input("maxima", "Maxima H− RT, Ficoll PM-400, actinomycin D, 42 °C overnight",
                     MaterialKind.CATALYST)
    hybrid = graph.product("rt_hybrid", "surface cDNA:RNA hybrid", MaterialKind.MOLECULAR,
                           state=native("hybrid"))
    graph.react("reverse_transcribe", "Extend the array oligo-dT on captured RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured, InputRole.TEMPLATE), (captured, InputRole.SELF_PRIMER),
                 (rt, InputRole.CATALYST)), (hybrid,))
    exo1b = graph.input("exo1_tissue", "Exonuclease I, 37 °C 45 min", MaterialKind.CATALYST)
    extended = graph.product("exo1_tissue_kept", "RNA-hybridized cDNA strands",
                             MaterialKind.MOLECULAR, state=native("hybrid"))
    graph.react("exo1_tissue_digest", "Exo I removes DNA not hybridized with mRNA",
                StepKind.DIGESTION,
                ((hybrid, InputRole.TEMPLATE), (exo1b, InputRole.CATALYST)), (extended,))
    pk = graph.input("proteinase_k", "tissue digestion buffer with proteinase K, 37 °C 40 min",
                     MaterialKind.CATALYST)
    cleared = graph.product("tissue_cleared", "array cDNA:RNA after tissue removal",
                            MaterialKind.MOLECULAR, state=native("cleared"))
    graph.react("digest_tissue", "Digest the tissue from the array", StepKind.DEPROTEINIZATION,
                ((extended, InputRole.TEMPLATE), (pk, InputRole.CATALYST)), (cleared,))
    naoh2 = graph.input("naoh_rna", "water, 3× 0.1 N NaOH 5 min, Tris, water washes",
                        MaterialKind.CONDITION)
    first = graph.product("first_strand", "single-stranded surface cDNA", MaterialKind.MOLECULAR,
                          state=native("first"))
    rna_waste = graph.product("rna_removed", "NaOH/water wash RNA products, digested tissue and DNA hydrolysates", MaterialKind.WASTE)
    graph.react("strip_rna", "Remove all mRNA with NaOH", StepKind.STRAND_SEPARATION,
                ((cleared, InputRole.TEMPLATE), (naoh2, InputRole.CONDITION)), (first, rna_waste))
    graph.discard(rna_waste, "removed by NaOH washes")
    random_primer = graph.input("random_primer", "TruSeq Read2-conjugated N9 random primer (IDT)",
                                MaterialKind.OLIGO,state=Construct([Segment("Read2 random-primer handle",sc.TR2_RANDOM[:-9]),Segment("nativeN9", "N"*9, feature=feature("seqscope_random_primer_umi","umi","random"))],"native randomprimer"))
    klenow = graph.input("klenow", "Klenow fragment (3′→5′ exo−), NEBuffer 2, 37 °C 2 h",
                         MaterialKind.CATALYST)
    duplex = graph.product("second_strand", "surface cDNA paired with random-primed second strand",
                           MaterialKind.MOLECULAR, state=native("second"))
    graph.react("second_strand_synthesis", "Random-primed second-strand synthesis on the array",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (random_primer, InputRole.STRUCTURAL),
                 (klenow, InputRole.CATALYST)), (duplex,))
    washed = _cleanup(graph, r, "wash_unbound", "Wash away DNA not bound to the array",
                      duplex, "three water washes", "array-bound second-strand duplexes",
                      "unbound primers and unanchored products", "washed off before elution")
    naoh3 = graph.input("naoh_elute", "2× 30 µL 0.1 N NaOH 5 min, neutralized with 3 M KOAc",
                        MaterialKind.CONDITION)
    eluted = graph.product("eluted", "eluted second strand", MaterialKind.MOLECULAR,
                           state=native("eluted"))
    surface = graph.product("surface_left", "array-anchored first-strand cDNA", MaterialKind.WASTE, state=native("first"))
    graph.react("elute_second_strand", "Elute the second strand with NaOH",
                StepKind.STRAND_SEPARATION,
                ((washed, InputRole.TEMPLATE), (naoh3, InputRole.CONDITION)), (eluted, surface))
    graph.discard(surface, "remains attached to the HDMI-array, which is not used further")
    purified = _cleanup(graph, r, "ampure_18", "AMPure XP 1.8× cleanup", eluted,
                        "AMPure XP 1.8×, 40 µL elution", "purified second strand",
                        "AMPure supernatant", "not bound by the beads")
    kapa = graph.input("kapa_pcr1", "KAPA HiFi HotStart ReadyMix, 13–15 cycles", MaterialKind.CATALYST)
    f1 = graph.input("pcr1_forward", "first-round forward primer (3′ phosphorothioates)",
                     MaterialKind.OLIGO, state=exact["pcr1_forward"])
    r1 = graph.input("pcr1_reverse", "first-round reverse primer (3′ phosphorothioates)",
                     MaterialKind.OLIGO, state=exact["pcr1_reverse"])
    pcr1 = graph.product("pcr1", "first-round PCR product", MaterialKind.MOLECULAR,
                         state=native("pcr1"))
    graph.react("pcr1_amplify", "First-round library PCR", StepKind.PCR,
                ((purified, InputRole.TEMPLATE), (f1, InputRole.STRUCTURAL),
                 (r1, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (pcr1,),
                geometry=PCRGeometry.unavailable(purified.state, reason="Actual incoming template is eluted ssDNA with N9 primer tail, not the surface duplex. First-cycle two-primer assembly remains unassigned; source primer sequences retained."))
    pcr1_clean = _cleanup(graph, r, "ampure_12", "AMPure XP 1.2× cleanup", pcr1,
                          "AMPure XP 1.2×", "purified first-round product",
                          "AMPure supernatant", "not bound by the beads")
    pcr2_template=graph.product("pcr2_template","10µl2nM PCR1 physical aliquot",MaterialKind.MOLECULAR,state=pcr1_clean.state)
    remainder=graph.product("pcr1_remainder","unsampled PCR1 preparation",MaterialKind.MOLECULAR,state=pcr1_clean.state)
    aliquot=graph.input("pcr2_aliquot_condition","source10µl of2nM first-round product",MaterialKind.CONDITION)
    graph.react("aliquot_pcr2","Allocate template for indexing PCR",StepKind.SPLIT,((pcr1_clean,InputRole.TEMPLATE),(aliquot,InputRole.CONDITION)),(pcr2_template,remainder))
    graph.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,"Physical unsampled preparation; downstream fate unspecified, not assumed discarded."))
    kapa2 = graph.input("kapa_pcr2", "KAPA HiFi HotStart ReadyMix, 8–9 cycles, 10 µL of 2 nM template",
                        MaterialKind.CATALYST)
    f2 = graph.input("pcr2_forward", "second-round P5 forward primer", MaterialKind.OLIGO,
                     state=exact["pcr2_forward"])
    r2 = graph.input("pcr2_reverse", "second-round P7 / 8-mer index reverse primer",
                     MaterialKind.OLIGO,state=Construct([Segment("sourceP7",sc.PCR2_REVERSE_5P),Segment("source8merindex","N"*8,feature=feature("seqscope_i7","sample_index","unknown")),Segment("sourceRead2andPS3end",sc.PCR2_REVERSE_3P,note="Two phosphorothioate bonds before terminal G,A on original primer; not copied as chemical bonds")],"native indexed reverse primer"))
    library = graph.product("pcr2", "indexed Seq-Scope library", MaterialKind.LIBRARY,
                            state=native("pcr2"))
    graph.react("pcr2_amplify", "Second-round indexing PCR", StepKind.PCR,
                ((pcr2_template, InputRole.TEMPLATE), (f2, InputRole.STRUCTURAL),
                 (r2, InputRole.STRUCTURAL), (kapa2, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(pcr2_template.state, reason="Actual cleaned PCR1 mixture enters indexing PCR. Exact primer sites on the reference can be computed separately; native mixture/PS-lineage assembly remains unassigned."))
    gel = _cleanup(graph, r, "gel_select", "Gel-excise 400–850 bp (Zymoclean)", library,
                   "agarose gel, Zymoclean Gel DNA Recovery", "400–850 bp library",
                   "out-of-range gel products", "outside the excised 400–850 bp window",
                   kind=StepKind.SELECTION, material=MaterialKind.LIBRARY)
    final = _cleanup(graph, r, "ampure_final", "AMPure XP 0.6–0.7× cleanup", gel,
                     "AMPure XP 0.6–0.7×", "final Seq-Scope library",
                     "AMPure supernatant", "not bound by the beads",
                     material=MaterialKind.LIBRARY)
    graph.sequence(final, SequencingPlan(library=FINAL_LIBRARY, primers=SEQ_PRIMERS, limitation=SEQUENCING_ENDING+" Actual pooled companion libraries and ratios remain unassigned; final molecular read-site reference does not establish class-pure recovery.",
                                         read_uses={"Read 1": "extract",
                                                    "Index 1 (i7)": "demultiplex",
                                                    "Read 2": "extract_then_include"}))
    graph.annotate(SourceDiscrepancy(
        "seed", "Seed end labels differ from the printed bases",
        (SourceClaim("STAR Methods, Generation of Seed HDMI-oligo Library (PMC8238917 JATS line 217)",
                     "backbone: (P5 sequence) (TR1) (HDMI) (HR1) (Oligo-dT) (DraI) (DraI-adapter) (P7 sequence)"),
         SourceClaim("Same section, HDMI-DraI and HDMI32-DraI sequences (JATS lines 221–225)",
                     "5′ CAAGCAGAAGACGGCATACGAGAT … GTGTAGATCTCGGTGGTCGCCGTATCATT 3′, i.e. the "
                     "Illumina P7 adapter and the reverse complement of P5")),
        "Segments are named for their bases (P7 adapter at the 5′ end, P5 complement at the 3′ "
        "end); the printed bases are used unchanged.",
        "The label/base identity is computed against the canonical Illumina sequences.",
        "The paper may label the ends by their side in the final library (TR1 side = P5 after "
        "PCR); it does not say so."))
    r.extend(("disassembly_condition","disassemble","pcr2_aliquot_condition","aliquot_pcr2","seed", "flow_cell", "cluster", "read1_drai", "first_seq_sbs", "first_seq",
              "drai", "drai_digest", "exo1_array", "exo1_array_digest", "naoh_array",
              "denature_array", "tissue", "rna_release", "release_rna", "capture", "maxima",
              "reverse_transcribe", "exo1_tissue", "exo1_tissue_digest", "proteinase_k",
              "digest_tissue", "naoh_rna", "strip_rna", "random_primer", "klenow",
              "second_strand_synthesis", "naoh_elute", "elute_second_strand", "kapa_pcr1",
              "pcr1_forward", "pcr1_reverse", "pcr1_amplify", "kapa_pcr2", "pcr2_forward",
              "pcr2_reverse", "pcr2_amplify", "ampure_final_kept"))
    graph.declare_reviewed_scope(r, AUTHORITY)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph()
