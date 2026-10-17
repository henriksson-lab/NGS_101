"""Protocol-owned spatial assay and reaction graph."""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.spatial import (SPATIAL_BEAD, UMI, SAMPLE, _finish_pcr,
    slide_seq_library)
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from chemdraw import MolecularState, Scene, Workflow, duplex_rows, revcomp
from protocol_graph import Disposition, Terminal
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SourceClaim,
                            SourceDiscrepancy, StepKind)
from protocol_model import chain, molecule, role, role_library

from types import SimpleNamespace

TITLE = 'Slide-seq (v1)'
NOTES = '01_slide_seq.html'
SOURCE = 'Defining source: <a href="https://doi.org/10.1126/science.aaw1219">Rodriques et al., Science (2019)</a>.'
SUMMARY = 'A pre-decoded monolayer of 10-µm DNA-barcoded beads captures RNA released from a tissue section for spatially registered 3′ RNA sequencing.'
CAVEAT = ("Provenance: bead oligos, Table S1 primers, reaction conditions and the read structure "
          "are from the published Supplementary Materials of the author manuscript (PMC6927209, "
          "NIHMS1036108 supplement-1). The Science-hosted supplement PDF (aaw1219_Rodriques_SM.pdf) "
          "was not obtained, so its identity with this author-manuscript supplement is unchecked. "
          "Table S1 prints two apparent errors (TSO 5′ end; decoding primer Truseq -1); see the "
          "discrepancy notices.")

# Published Table S1 (supplement-1.docx text lines 310-325): SOLiD decoding primers, all
# /5Phos/. 'Truseq' anneals to the bead TruSeq handle and 'UP' to the bead linker; the -n
# primers are recessed by n 5' bases. Values as printed.
TRUSEQ_PCR_HANDLE = "CTACACGACGCTCTTCCGATCT"
BEAD_LINKER = "TCTTCAGCGTTCCCGAGA"
DECODING_PRIMERS = (
    ("Truseq", "AGATCGGAAGAGCGTCGTGTAG"), ("Truseq -1", "GATCGGAAGAGCGTCGTCTAG"),
    ("Truseq -2", "ATCGGAAGAGCGTCGTGTAG"), ("TruSeq-3", "TCGGAAGAGCGTCGTGTAG"),
    ("TruSeq-4", "CGGAAGAGCGTCGTGTAG"), ("UP", "TCTCGGGAACGCTGAAGA"),
    ("UP-1", "CTCGGGAACGCTGAAGA"), ("UP-2", "TCGGGAACGCTGAAGA"),
    ("UP-3", "CGGGAACGCTGAAGA"), ("UP-4", "GGGAACGCTGAAGA"))


def _recessed_series_mismatches():
    """Printed decoding primers that are not the 5'-recessed form of their parent."""
    parents = {"Truseq": revcomp(TRUSEQ_PCR_HANDLE), "UP": revcomp(BEAD_LINKER)}
    out = []
    for name, seq in DECODING_PRIMERS:
        family = "UP" if name.startswith("UP") else "Truseq"
        n = int(name.rsplit("-", 1)[1]) if "-" in name else 0
        expected = parents[family][n:]
        if seq != expected:
            out.append((name, seq, expected))
    return tuple(out)


DECODING_MISMATCHES = _recessed_series_mismatches()
FINAL_CAPTION = "The barcode read identifies a bead's decoded position and its UMI; the paired insert read identifies the captured transcript."
SEQUENCING_INTRO = "The barcode read identifies a bead's decoded position and its UMI; the paired insert read identifies the captured transcript."
READ_LENGTHS = {'Read 1': 42, 'Index 1 (i7)': 8, 'Read 2': 50}
READ_USES = {'Read 1': 'extract', 'Index 1 (i7)': 'demultiplex', 'Read 2': 'include'}
SEQUENCING_ENDING = ''
SEQUENCING_UNAVAILABLE = ''
oligos = lambda: ("<p>Published supplement bead oligo Sequence 2 (representative; 5′ linker): 5'- TTTTTTTTCTACACGACGCTCTTCCGATCT<cbc>JJJJJJJJ</cbc>TCTTCAGCGTTCCCGAGA<cbc>JJJJJJJ</cbc><umi>NNNNNNNN</umi>TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -3'</p>", "<p>SMART PCR primer: 5'- AAGCAGTGGTATCAACGCAGAGT -3'</p>", "<p>TruSeq PCR handle primer: 5'- CTACACGACGCTCTTCCGATCT -3'</p>", "<p>TruSeq5 tagmentation PCR primer: 5'- AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -3'</p>", "<p>Template-switch oligo, Table S1 as printed (r = RNA, + = LNA; see discrepancy notice): 5'- AAGCTGGTATCAACGCAGAGTGAATrG+GrG -3'</p>") + tuple(
    f"<p>SOLiD decoding primer {name} (Table S1 as printed): 5'- /5Phos/{seq} -3'</p>"
    for name, seq in DECODING_PRIMERS)
FINAL_LIBRARY, SEQ_PRIMERS = slide_seq_library()
_definition = SimpleNamespace(title=TITLE, doi='10.1126/science.aaw1219', citation='Rodriques et al., Science (2019)', endpoint="The barcode read identifies a bead's decoded position and its UMI; the paired insert read identifies the captured transcript.")


def _slide_seq():
    return chain(molecule("10-µm bead puck", (SPATIAL_BEAD, UMI, role("poly(dT) capture oligo"))), (
        ("Decode each bead position by SOLiD sequencing-by-ligation",
         molecule("mapped bead puck", (SPATIAL_BEAD, role("known x–y coordinate"))),
         "The decoded bead sequence is a spatial address before tissue is applied."),
        ("Transfer a fresh-frozen tissue section and hybridize released RNA (6× SSC, 15 min)",
         molecule("RNA captured on bead", (SPATIAL_BEAD, UMI, role("poly(dT)"),
                                             role("poly(A) RNA"))),
         "RNA inherits position by hybridizing to an oligo on a mapped bead."),
        ("Reverse-transcribe on the puck with Maxima H− and the template-switch oligo",
         molecule("Slide-seq cDNA", (SPATIAL_BEAD, UMI, role("3-prime transcript cDNA"),
                                    role("TSO complement"))),
         "Published supplement: 1 h at 42 °C; tissue is then digested with proteinase K/SDS."),
        ("Release and wash beads; amplify cDNA with TruSeq-handle and SMART PCR primers",
         molecule("amplified Slide-seq cDNA", (SPATIAL_BEAD, UMI,
             role("3-prime transcript cDNA"), role("SMART end"))),
         "4 + 9 cycles in KAPA HiFi HotStart ReadyMix, then 0.6× AMPure XP."),
        ("Tagment 600 pg with Nextera XT and amplify with Truseq5 + N7xx (12 cycles)",
         molecule("Slide-seq sequencing library", (SPATIAL_BEAD, UMI,
             role("3-prime transcript insert"), SAMPLE)),
         "Sequencing joins each transcript molecule to the bead position decoded on the puck."),
    ))

workflow = _slide_seq

def _slide_graph(d, final, primers):
    source = f"{d.citation}; doi:{d.doi}"
    supp = "Rodriques et al. 2019 Supplementary Materials (PMC6927209 author-manuscript supplement-1)"
    graph = ProtocolGraph(d.title, authority=source)
    puck = graph.input("puck", "DNA-barcoded 10-micrometre bead puck (ChemGenes beads, Drop-seq synthesis)",
                       MaterialKind.BIOLOGICAL, state=_slide_seq().initial, primary=True)
    decoder = graph.input("decoder", "ten 5′-phosphorylated decoding primers (Table S1: Truseq, "
                          "Truseq -1 to -4 on the bead TruSeq handle; UP, UP-1 to -4 on the bead linker), "
                          "5 µM in 4× SSC, and ligation mix (T4 DNA ligase rapid, 40× diluted SOLiD SR-75 "
                          "sequencing oligo); two ligations per primer, six constant bases, 14-base barcodes",
                          MaterialKind.OLIGO, source=supp + ", Puck Sequencing and Table S1")
    mapped = graph.product("mapped_puck", "position-decoded bead puck",
                           MaterialKind.MOLECULAR,
                           state=molecule("mapped bead puck", (SPATIAL_BEAD,
                               role("known x-y coordinate"))))
    graph.react("decode", "Decode each bead position by primer hybridization, two ligations "
                "(SOLiD buffer C/B cleavage between) and formamide stripping per primer",
                StepKind.HYBRIDIZATION,
                ((puck, InputRole.TEMPLATE), (decoder, InputRole.STRUCTURAL)), (mapped,),
                source=supp, changes_molecule=False)
    rna = graph.input("tissue_rna", "poly(A) RNA released from the 10-µm fresh-frozen tissue section",
                      MaterialKind.MOLECULAR)
    hyb = graph.input("hyb_buffer", "6× SSC with 2 U/µl RNase inhibitor, 15 min room temperature",
                      MaterialKind.CONDITION, source=supp)
    captured = graph.product("captured_rna", "RNA captured on mapped beads",
                             MaterialKind.MOLECULAR,
                             state=molecule("RNA captured on bead", (SPATIAL_BEAD, UMI,
                                 role("poly(A) RNA"))))
    graph.react("capture", "Capture tissue RNA on mapped bead oligo(dT)",
                StepKind.HYBRIDIZATION,
                ((mapped, InputRole.TEMPLATE), (rna, InputRole.STRUCTURAL),
                 (hyb, InputRole.CONDITION)), (captured,), source=supp)
    tso = graph.input("tso", "template-switch oligo AAGCTGGTATCAACGCAGAGTGAATrG+GrG (as printed)",
                      MaterialKind.OLIGO, source=supp + ", Table S1")
    rt_primer = graph.input("rt_primer", "bead-tethered oligo(dT) capture oligo (RT primer)",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "Maxima H− reverse transcriptase in Maxima RT buffer "
                     "with Ficoll PM-400, dNTPs and RNase inhibitor; 1 h 42 °C",
                     MaterialKind.CATALYST, source=supp)
    cdna = graph.product("cdna", "bead-barcoded 3-prime cDNA extended across the template-switch oligo",
                         MaterialKind.MOLECULAR,
                         state=molecule("Slide-seq cDNA", (SPATIAL_BEAD, UMI,
                             role("3-prime transcript cDNA"), role("TSO complement"))))
    graph.react("reverse_transcribe", "Reverse-transcribe captured RNA from the bead oligo(dT) "
                "with template switching",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (tso, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)), (cdna,), source=supp)
    digest = graph.input("tissue_digestion", "2× tissue digestion buffer (SDS, EDTA, proteinase K), "
                         "37 °C 40 min", MaterialKind.CATALYST, source=supp)
    digested = graph.product("digested", "bead cDNA in digested tissue/RT mixture",
                             MaterialKind.MOLECULAR, state=cdna.state)
    graph.react("digest_tissue", "Digest the tissue section in the RT mixture",
                StepKind.DEPROTEINIZATION,
                ((cdna, InputRole.TEMPLATE), (digest, InputRole.CATALYST)), (digested,),
                source=supp, changes_molecule=False)
    wash = graph.input("bead_washes", "pipette beads off the glass; wash buffer (Tris/EDTA/Tween-20), "
                       "3000 RCF spins, then water", MaterialKind.SELECTION, source=supp)
    washed = graph.product("washed_beads", "washed cDNA-bearing beads", MaterialKind.MOLECULAR,
                           state=cdna.state)
    wash_waste = graph.product("wash_waste", "glass substrate, digested tissue, RT reagents and "
                               "supernatants", MaterialKind.WASTE)
    graph.react("wash_beads", "Release beads from the glass and wash by centrifugation",
                StepKind.CLEANUP,
                ((digested, InputRole.TEMPLATE), (wash, InputRole.SELECTION)),
                (washed, wash_waste), source=supp)
    graph.discard(wash_waste, "glass discarded and supernatants removed after each spin")
    handle = graph.input("handle_primer", "Truseq_PCR_handle primer CTACACGACGCTCTTCCGATCT",
                         MaterialKind.OLIGO, source=supp + ", Table S1")
    smart = graph.input("smart_primer", "SMART PCR primer AAGCAGTGGTATCAACGCAGAGT",
                        MaterialKind.OLIGO, source=supp + ", Table S1")
    kapa = graph.input("kapa", "KAPA HiFi HotStart ReadyMix; 95 °C 3 min, 4 cycles (65 °C anneal), "
                       "9 cycles (67 °C anneal), 72 °C 5 min", MaterialKind.CATALYST, source=supp)
    amplified = graph.product("amplified_cdna", "amplified bead-barcoded cDNA",
                              MaterialKind.MOLECULAR,
                              state=molecule("amplified Slide-seq cDNA", (SPATIAL_BEAD, UMI,
                                  role("3-prime transcript cDNA"), role("SMART end"))))
    graph.react("cdna_pcr", "Amplify bead cDNA from the TruSeq handle and the SMART site",
                StepKind.PCR,
                ((washed, InputRole.TEMPLATE), (handle, InputRole.STRUCTURAL),
                 (smart, InputRole.STRUCTURAL), (kapa, InputRole.CATALYST)), (amplified,),
                source=supp, geometry=PCRGeometry.unavailable(amplified.state,
                    reason="Bead-oligo and TSO-end junction bases are not drawn; the printed TSO only "
                           "partly matches the SMART PCR primer (see discrepancy notice)."))
    ampure = graph.input("ampure_cdna", "AMPure XP 30 µl per 50 µl PCR (0.6×); elute 10 µl",
                         MaterialKind.SELECTION, source=supp)
    clean = graph.product("clean_cdna", "purified amplified cDNA (600 pg taken forward)",
                          MaterialKind.MOLECULAR, state=amplified.state)
    ampure_waste = graph.product("ampure_cdna_waste", "0.6× AMPure supernatant, primers and washes",
                                 MaterialKind.WASTE)
    graph.react("cdna_cleanup", "Purify amplified cDNA with 0.6× AMPure XP",
                StepKind.CLEANUP,
                ((amplified, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
                (clean, ampure_waste), source=supp)
    graph.discard(ampure_waste, "removed with the bead supernatant")
    transposome = graph.input("transposome", "Nextera XT tagmentation reagent (FC-131-1096)",
                              MaterialKind.ADAPTER, source=supp)
    tn5 = graph.input("tn5", "Tn5 transposase (in Nextera XT)", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cdna", "tagmented Slide-seq cDNA", MaterialKind.MOLECULAR,
                           state=cdna.state)
    waste = graph.product("tn5_waste", "internal and SMART-end tagmentation fragments without the "
                          "bead TruSeq handle", MaterialKind.MOLECULAR)
    graph.react("tagment", "Tagment 600 pg amplified cDNA (Nextera XT, manufacturer's instructions)",
                StepKind.TAGMENTATION,
                ((clean, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste), source=supp)
    graph.mark_terminal(Terminal(waste.id, Disposition.UNRESOLVED, 'Candidate non-library tagmentation products enter PCR and the subsequent AMPure cleanup. Primer/end eligibility does not establish physical removal of these classes.'))
    reviewed = ["puck", "decoder", "decode", "tissue_rna", "hyb_buffer", "capture", "rt_primer", "tso",
                "reverse_transcriptase", "reverse_transcribe", "tissue_digestion", "digested",
                "digest_tissue", "bead_washes", "washed_beads", "wash_waste", "wash_beads",
                "handle_primer", "smart_primer", "kapa", "amplified_cdna", "cdna_pcr",
                "ampure_cdna", "clean_cdna", "ampure_cdna_waste", "cdna_cleanup",
                "transposome", "tn5", "tagment"]
    _finish_pcr(graph, "library", tagged, final, primers, d.endpoint, reviewed,
                label="Slide-seq (Truseq5 + Nextera N7xx, 12 cycles; 0.6× AMPure)")
    graph.annotate(SourceDiscrepancy("tso", "Template-switch oligo versus SMART PCR primer",
        (SourceClaim("Published supplement Table S1 (PMC NIHMS1036108 supplement-1.docx, text line 316); "
                     "bioRxiv v1 Table S1 (media-1 p.22) prints the same value",
                     "TSO AAGCTGGTATCAACGCAGAGTGAATrG+GrG"),
         SourceClaim("same published Table S1 (text line 314)",
                     "SMART PCR primer AAGCAGTGGTATCAACGCAGAGT; the printed TSO lacks the primer's "
                     "AG at positions 5–6 and contains only its 3′ 17 nt"),
         SourceClaim("Stickels et al. 2021 (Slide-seqV2) Supplementary Table 6 and Online Methods "
                     "(same Qiagen order #339414YCO0076714)",
                     "TSO AAGCAGTGGTATCAACGCAGAGTGAATrG+GrG")),
        "Display the TSO exactly as printed; the cDNA PCR is modelled at role level and no SMART-end "
        "bases are inferred. The SMART end is lost from the final library by tagmentation.",
        "The published supplement is a Word document whose text layer gives the value directly; the "
        "preprint rendered PDF page 22 also shows it visibly printed.",
        "Whether the TSO omission is a transcription error is unresolved; the Science-hosted "
        "supplement PDF was not obtained."))
    for name, printed, expected in DECODING_MISMATCHES:
        graph.annotate(SourceDiscrepancy("decoder", f"Decoding primer {name} versus its recessed series",
            (SourceClaim("Published supplement Table S1 (supplement-1.docx text lines 317-321)",
                         f"{name} /5Phos/{printed}"),
             SourceClaim("Computed from the same table: the primer series is the reverse complement "
                         "of the bead TruSeq handle CTACACGACGCTCTTCCGATCT, recessed by n 5′ bases "
                         "(Truseq, -2, -3, -4 all follow it)",
                         f"expected /5Phos/{expected}")),
            "Show the primer exactly as printed; the decoding step is modelled at role level and does "
            "not use the printed bases, so no corrected sequence enters the graph.",
            "The other nine primers match their computed recessed forms exactly.",
            "Whether the order used the printed or the series sequence is not stated; the mismatch "
            "is at the fourth base from the primer 3′ end (C printed, G expected)."))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _slide_graph(_definition, FINAL_LIBRARY, SEQ_PRIMERS)
