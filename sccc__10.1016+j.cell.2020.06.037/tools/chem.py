from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry
from types import SimpleNamespace
from protocol_model import chain, role_library
from multimodal_protocols import CELL, CircularizationGeometry, Definition, InputRole, MaterialKind, ProtocolGraph, SAMPLE, StepKind, UMI, _finish_pcr, molecule, role, role_library
DEFINITION = Definition("scCC — single-cell calling cards", "10.1016/j.cell.2020.06.037",
        "Moudgil et al., Cell (2020)", "Self-reporting transposons encode their genomic insertion flank in RNA; circularization brings droplet identifiers next to that junction so TF-binding history can be joined to gene expression.",
        "The calling-card read resolves the insertion flank and retains the same cell barcode and UMI as the parent 10x cDNA.", "01_sccc.html")


def _sccc_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cell = graph.input("cell", "cell carrying an integrated self-reporting transposon",
                       MaterialKind.BIOLOGICAL, state=_sccc().initial, primary=True)
    rna_polymerase = graph.input("rna_polymerase", "cellular RNA polymerase",
                                 MaterialKind.CATALYST)
    transcript = graph.product("srt_transcript", "self-reporting transcript containing the insertion flank",
                               MaterialKind.MOLECULAR,
                               state=molecule("SRT transcript", (role("reporter transcript"),
                                   role("terminal repeat"), role("genomic insertion flank"), role("poly(A)"))))
    graph.react("transcribe_srt", "Transcribe through the SRT into flanking genomic DNA",
                StepKind.TRANSCRIPTION,
                ((cell, InputRole.TEMPLATE), (rna_polymerase, InputRole.CATALYST)),
                (transcript,))
    bead_primer = graph.input("bead_primer", "droplet oligo(dT) cell-barcode/UMI primer",
                              MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("barcoded_cdna", "cell-barcoded SRT cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("barcoded SRT cDNA", (CELL, UMI,
                             role("insertion-flank cDNA"), role("terminal-repeat cDNA"))))
    graph.react("reverse_transcribe", "Reverse-transcribe and barcode the SRT transcript",
                StepKind.REVERSE_TRANSCRIPTION,
                ((transcript, InputRole.TEMPLATE), (bead_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    expression = graph.product("expression_branch", "expression-library aliquot",
                               MaterialKind.MOLECULAR, state=cdna.state)
    junction = graph.product("junction_branch", "SRT-junction aliquot",
                             MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("split_branches", "Split expression and SRT-junction branches",
                StepKind.SPLIT, ((cdna, InputRole.TEMPLATE),), (expression, junction),
                changes_molecule=False)
    reviewed = ["cell", "rna_polymerase", "transcribe_srt", "bead_primer",
                "reverse_transcriptase", "reverse_transcribe", "split_branches"]
    expression_final = role_library("scCC expression-library roles",
                                    (CELL, UMI, role("transcript cDNA")), index=SAMPLE)
    tso = graph.input("expression_tso", "10x_TSO added back only to expression aliquot", MaterialKind.OLIGO)
    expression_rt = graph.input("expression_rt", "expression-branch reverse transcriptase", MaterialKind.CATALYST)
    switched = graph.product("expression_switched", "TSO-completed expression cDNA", MaterialKind.MOLECULAR, state=expression.state)
    graph.react("restore_expression_tso", "Add back TSO and complete expression-branch RT", StepKind.TEMPLATE_SWITCHING,
                ((expression, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (expression_rt, InputRole.CATALYST)), (switched,))
    reviewed.extend(("expression_tso", "expression_rt", "restore_expression_tso"))
    _finish_pcr(graph, "expression", switched, expression_final,
                "Expression branch: droplet transcriptome reads linked by cell barcode.",
                reviewed, label="scCC expression")
    junction_forward = graph.input("srt_enrich_forward", "Bio_Illumina_Seq1_scCC_10X_3xPT", MaterialKind.OLIGO)
    junction_reverse = graph.input("srt_enrich_reverse", "Bio_Long_PB_LTR_3xPT", MaterialKind.OLIGO)
    junction_polymerase = graph.input("srt_enrich_polymerase", "SRT-enrichment PCR polymerase", MaterialKind.CATALYST)
    amplicon = graph.product("junction_amplicon", "phosphorylated internally biotinylated SRT amplicons", MaterialKind.MOLECULAR, state=junction.state)
    graph.react("amplify_junction", "Amplify the SRT junction with phosphorylated biotinylated primers", StepKind.PCR,
                ((junction, InputRole.TEMPLATE), (junction_forward, InputRole.STRUCTURAL),
                 (junction_reverse, InputRole.STRUCTURAL), (junction_polymerase, InputRole.CATALYST)), (amplicon,),
                geometry=PCRGeometry.unavailable(final, reason="Original Table S4 publishes both primers; their source-exact cDNA binding geometry has not yet been assembled in this graph."))
    reviewed.extend(("srt_enrich_forward", "srt_enrich_reverse", "srt_enrich_polymerase", "amplify_junction"))
    circle_ligase = graph.input("circle_ligase", "SRT amplicon circularization ligase",
                                MaterialKind.CATALYST)
    circle = graph.product("circle", "circular SRT amplicon", MaterialKind.MOLECULAR,
                           state=amplicon.state)
    graph.react("circularize", "Circularize the SRT amplicon", StepKind.CIRCULARIZATION,
                ((amplicon, InputRole.TEMPLATE), (circle_ligase, InputRole.CATALYST)), (circle,),
                geometry=CircularizationGeometry.known(
                    "double-strand intramolecular ligation",
                    statement="The SRT branch circularizes its amplicon before exonuclease selection."))
    exonuclease = graph.input("exonuclease", "linear-DNA exonuclease",
                              MaterialKind.CATALYST)
    cleaned = graph.product("clean_circle", "exonuclease-resistant SRT circles",
                            MaterialKind.MOLECULAR, state=circle.state)
    linear_waste = graph.product("linear_waste", "digested linear DNA", MaterialKind.WASTE)
    graph.react("remove_linear", "Remove residual linear DNA", StepKind.DIGESTION,
                ((circle, InputRole.TEMPLATE), (exonuclease, InputRole.CATALYST)),
                (cleaned, linear_waste))
    graph.discard(linear_waste, "does not contain retained circular SRT junction molecules")
    shear = graph.input("shear_condition", "DNA shearing condition", MaterialKind.CONDITION)
    fragments = graph.product("junction_fragments", "junction-proximal circle fragments",
                              MaterialKind.MOLECULAR, state=cleaned.state)
    graph.react("shear", "Shear the circularized product", StepKind.FRAGMENTATION,
                ((cleaned, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    capture = graph.input("capture_probe", "junction capture reagent", MaterialKind.SELECTION)
    captured = graph.product("captured_junction", "captured insertion-junction fragments",
                             MaterialKind.MOLECULAR, state=fragments.state)
    uncaptured = graph.product("uncaptured", "uncaptured circle fragments", MaterialKind.WASTE)
    graph.react("capture_junction", "Capture the junction-bearing branch",
                StepKind.AFFINITY_CAPTURE,
                ((fragments, InputRole.TEMPLATE), (capture, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "does not carry the selected SRT junction")
    reviewed.extend(("circle_ligase", "circularize", "exonuclease", "remove_linear",
                     "shear_condition", "shear", "capture_probe", "capture_junction"))
    endrepair = graph.input("end_repair_enzyme", "end-repair enzyme mixture", MaterialKind.CATALYST)
    repaired = graph.product("repaired_junction", "end-repaired captured junction fragments", MaterialKind.MOLECULAR, state=captured.state)
    graph.react("end_repair", "End-repair captured junction fragments on beads", StepKind.END_REPAIR,
                ((captured, InputRole.TEMPLATE), (endrepair, InputRole.CATALYST)), (repaired,))
    tailase = graph.input("tailing_enzyme", "A-tailing enzyme mixture", MaterialKind.CATALYST)
    tailed = graph.product("tailed_junction", "A-tailed captured junction fragments", MaterialKind.MOLECULAR, state=repaired.state)
    graph.react("a_tail", "A-tail the captured junction fragments", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailase, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("custom_adapters", "preannealed scCC P5/P7 adapters", MaterialKind.ADAPTER)
    ligase = graph.input("adapter_ligase", "adapter ligase", MaterialKind.CATALYST)
    tagged = graph.product("adapter_tagged_junction", "custom-adapter-tagged junction fragments", MaterialKind.MOLECULAR, state=tailed.state)
    graph.react("ligate_custom_adapters", "Ligate the custom scCC P5/P7 adapters", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL), (ligase, InputRole.CATALYST)), (tagged,),
                geometry=AdapterLigationGeometry.unavailable("Original Table S4 prints the custom adapters and the Methods specify annealing; this graph has not yet assembled the adapter-end Scene."))
    reviewed.extend(("end_repair_enzyme", "end_repair", "tailing_enzyme", "a_tail", "custom_adapters", "adapter_ligase", "ligate_custom_adapters"))
    _finish_pcr(graph, "junction", tagged, final, d.endpoint, reviewed,
                label="scCC insertion junction")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'sccc'
    return _sccc_graph(d, final)

def _sccc():
    return chain(molecule("cell carrying self-reporting transposons", (
        role("TF-directed piggyBac transposase"), role("integrated SRT"))), (
        ("Transcribe through the SRT terminal repeat into flanking genomic DNA",
         molecule("self-reporting RNA", (role("reporter transcript"),
                                          role("transposon terminal repeat"),
                                          role("genomic insertion flank"), role("poly(A)"))),
         "Removing the transposon's polyadenylation signal makes the insertion junction part of the RNA 3′ UTR."),
        ("Run droplet RT and split the barcoded first-strand product",
         molecule("SRT cDNA", (CELL, UMI, role("insertion-flank cDNA"),
                                role("terminal-repeat cDNA"))),
         "One aliquot becomes the expression library; the other is amplified specifically for SRT transcripts."),
        ("Circularize the SRT amplicon, remove linear DNA and shear",
         molecule("junction-proximal fragment", (CELL, UMI,
                                                  role("genomic insertion flank"),
                                                  role("terminal repeat"))),
         "Circularization moves the 3′ cell barcode and UMI next to the transposon–genome junction before short-read library construction."),
        ("Capture the junction-bearing branch and add custom Illumina adapters",
         molecule("scCC library", (role("custom P5 end"), CELL, UMI,
                                    role("insertion flank"), role("terminal repeat"),
                                    role("custom P7 end"))),
         "The shared droplet barcode joins insertion sites to the separately sequenced transcriptome."),
    ))

d = DEFINITION
FINAL_LIBRARY = role_library(d.title + " final library roles", (CELL, UMI, role("transposon–genome junction")), index=SAMPLE)
AUTHOR_PROTOCOL = "https://dx.doi.org/10.17504/protocols.io.xwifpce"
OLIGOS = (
    "<b>Biotinylated junction primer</b> — 5′-/5Phos/ACACTCTTTCCC/iBiodT/ACACGACGCTCTTCCGA*T*C*T-3′",
    "<b>Biotinylated piggyBac primer</b> — 5′-/5Phos/GCGTCAATTTTACGCAGAC/iBiodT/ATCTTTC*T*A*G-3′",
    "<b>scCC P5 adapter</b> — 5′-AATGATACGGCGACCACCGAGATCTTCACTCATTCCACACGACTCCTTGCCAGTCTC*T-3′",
    "<b>scCC P7 adapter</b> — 5′-/5Phos/GAGACTGGCAAGTACACGTCGCACTCACCATGA[index]ATCTCGTATGCCGTCTTCTGCTTG-3′",
    "<b>scCC P5 PCR primer</b> — 5′-AATGATACGGCGACCACCGAGATC-3′",
    "<b>scCC P7 PCR primer</b> — 5′-CAAGCAGAAGACGGCATACGAGAT-3′",
    "<b>scCC custom Read 2 primer</b> — 5′-CGTGTAGGGAAAGAGTGTGCGTCAATTTTACGCAGACTATCTTTCTAG-3′",
    "<b>scCC custom Index 1 primer</b> — 5′-GAGACTGGCAAGTACACGTCGCACTCACCATGA-3′",
)
SPEC = SimpleNamespace(
    TITLE=d.title,
    NOTES=d.notes,
    SOURCE=(f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>; '
            f'author <a href="{AUTHOR_PROTOCOL}">single-cell calling-cards library protocol</a>.'),
    SUMMARY=d.summary,
    CAVEAT="The final Cell Table S4 and author library protocol publish the junction primers and custom adapter oligos below. Their exact landing sites are not yet assembled in this role-level library model.",
    FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint,
    SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE="Final Cell Table S4 prints custom PCR and sequencing oligos, but their binding sites are not yet assembled in this role-level library model; index members vary by library.",
    workflow=_sccc,
    oligos=lambda: OLIGOS,
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
