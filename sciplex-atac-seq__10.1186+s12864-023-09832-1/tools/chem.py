"""Protocol-owned spatial assay and reaction graph."""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
from protocol_families.spatial import UMI, CELL_A, CELL_B, HASH, SAMPLE
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from pcr import PCRGeometry
from protocol_model import chain, molecule, role, role_library

from types import SimpleNamespace

TITLE = 'sciPlex-ATAC-seq'
NOTES = '01_sciplex_atac_seq.html'
SOURCE = ('Source: <a href="https://doi.org/10.1186/s12864-023-09832-1">Booth et al., <i>BMC Genomics</i> (2023)</a>, '
          'published version of <a href="https://doi.org/10.1101/2023.03.05.531201">bioRxiv 2023.03.05.531201</a>, '
          'with its Supp. Table 4 and the linked author two-level protocols.io procedure.')
SUMMARY = 'Permeabilized nuclei absorb sample-specific ssDNA hashes, then two-level combinatorial ATAC indexing assigns a separate cell identity.'
CAVEAT = ('<b>Published source discrepancy, two-level hash-capture oligo.</b> The published '
          '<i>BMC Genomics</i> Methods ("Co-capture of hash oligo and ATAC profiles with two-level '
          'sci-ATAC") print the capture oligo with + before every third T of the poly(T) tract and define '
          '+ as a locked nucleic acid (LNA). The author two-level protocols.io procedure linked from that '
          'article (pp. 3–4) prints the same oligo with an unmodified T30VN tract. This page displays the '
          'published LNA formulation; the role-level graph does not depend on either form, and which '
          'formulation was used for each published sample remains unresolved. '
          '<b>Hash-read assignment.</b> The published Sequencing section states that hash Read 1 (18 cycles) '
          'reads the 8-nt UMI plus the 10-nt PCR (P5) barcode and that Index 2 reads the capture-oligo '
          'barcode; by the printed sequences, the Nextera Read 1 site of the capture oligo precedes its '
          'N8 UMI and 10-nt well barcode, while Index 2 would read the P5 PCR primer index. The model '
          'therefore does not assign hash-read cycles to a barcode until this is resolved. '
          '<b>Gel windows.</b> Published Methods cut the hash band at 199 bp and ATAC at ~200\u20133000 bp; the '
          'linked protocols.io procedure (step 8) cuts ~200 bp and 300\u20131200 bp. The graph labels the '
          'published windows.')
FINAL_CAPTION = 'Hash reads identify sample origin; the Tn5 (or capture-oligo) well index plus the PCR well indices identify the nucleus; paired genomic reads report accessibility.'
SEQUENCING_INTRO = FINAL_CAPTION
READ_LENGTHS = None
READ_USES = None
SEQUENCING_ENDING = FINAL_CAPTION
SEQUENCING_UNAVAILABLE = ('Published Sequencing section: the two-level ATAC library uses custom primers '
    '(Read 1 GCGATCGAGGACGGCAGATGTGTATAAGAGACAG, Read 2 CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG and '
    'matching index primers) with Read 1 51, Read 2 51, Index 1 8 + 27 dark + 10 and Index 2 10 + 21 dark '
    '+ 8 cycles (T7:P7:P5:T5 index read); the 199-bp hash library is sequenced separately with standard '
    'primers (Read 1 18, Read 2 ≥16, Index 1 10, Index 2 10 cycles). Individual hash, Tn5 and PCR '
    'index choices vary, so one exact primer-landing molecule is not drawn.')
def oligos():
    return (
        "<p>Published Methods, sample hash (5′→3′): GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG<cbc>XXXXXXXXXX</cbc>B" + "A" * 30 + "</p>",
        "<p>Published Methods, two-level capture oligo (5′→3′; + = next base is LNA): TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG<umi>NNNNNNNN</umi><cbc>XXXXXXXXXX</cbc>" + "+".join(["TTT"] * 10) + "VN</p>",
        f"<p>Published Supp. Table 4, P5_i5_1 indexed Tn5 transferred strand (5′→3′): {P5_I5_1_TRANSFERRED[:21]}<cbc>{P5_I5_1_TRANSFERRED[21:29]}</cbc>{P5_I5_1_TRANSFERRED[29:]}</p>",
        f"<p>Published Supp. Table 4, P7_i7_1 indexed Tn5 transferred strand (5′→3′): {P7_I7_1_TRANSFERRED[:27]}<cbc>{P7_I7_1_TRANSFERRED[27:35]}</cbc>{P7_I7_1_TRANSFERRED[35:]}</p>",
        f"<p>Published Supp. Table 4, Tn5ME_Rev loading strand (5′→3′): /5Phos/{TN5ME_REV}</p>",
    )
# Supp. Table 4 (published BMC Genomics MOESM2), representative members as ordered;
# Tn5ME_Rev carries a 5' phosphate.
P5_I5_1_TRANSFERRED = "TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG"
P7_I7_1_TRANSFERRED = "GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG"
TN5ME_REV = "CTGTCTCTTATACACATCT"
TN5_LOADING_STRANDS = (P5_I5_1_TRANSFERRED, P7_I7_1_TRANSFERRED, TN5ME_REV)
FINAL_LIBRARY = role_library(TITLE + " final library roles", (HASH, CELL_A, CELL_B, role("ATAC insert")), index=SAMPLE)
SEQ_PRIMERS = ()
_definition = SimpleNamespace(title=TITLE, doi='10.1186/s12864-023-09832-1', citation='Booth et al., BMC Genomics 24:737 (2023), doi:10.1186/s12864-023-09832-1; bioRxiv 10.1101/2023.03.05.531201', endpoint='Hash reads identify sample origin; the ordered combinatorial indices identify the nucleus; paired genomic reads report accessibility.')


def _sciplex_atac():
    return chain(molecule("permeabilized nuclei separated by sample", (role("accessible chromatin"),)), (
        ("Let nuclei absorb a sample-specific unmodified ssDNA hash oligo; fix, pool and freeze",
         molecule("hashed nucleus", (HASH, role("accessible chromatin"))),
         "The hash records specimen or treatment; it is not a cell barcode."),
        ("Anneal and extend the well-barcoded capture oligo on the hash poly(A); tagment with indexed Tn5",
         molecule("first-round indexed nucleus", (HASH, UMI, CELL_A, role("tagmented accessible DNA"))),
         "The capture-oligo well barcode (hash) and Tn5 N5/N7 indices (chromatin) record the first well."),
        ("Pool, FACS-sort into PCR wells, reverse cross-links and amplify with indexed Nextera primers",
         molecule("two-level cell identity", (HASH, UMI, CELL_A, CELL_B,
                                                role("accessible-DNA fragment"), SAMPLE)),
         "The first-well and PCR-well indices identify a nucleus; the independent hash identifies its sample."),
        ("Gel-separate hash and ATAC libraries",
         molecule("sciPlex-ATAC readouts", (HASH, CELL_A, CELL_B,
                                             role("ATAC insert"), SAMPLE)),
         "Two-level experiments sequence hash and ATAC products separately; three-level variants intermix them."),
    ))

workflow = _sciplex_atac

def _sciplex_graph(d, final):
    source = f"{d.citation}"
    graph = ProtocolGraph(d.title, authority=source)
    nuclei = graph.input("nuclei", "nuclei permeabilized in cold lysis buffer, separated by sample well",
                         MaterialKind.BIOLOGICAL, state=_sciplex_atac().initial, primary=True)
    hash_oligo = graph.input("hash_oligo", "sample-specific ssDNA hash (10-nt hash ID, B, A30); ~0.5 pmol per 1000 cells",
                             MaterialKind.OLIGO)
    hashed = graph.product("hashed", "nuclei carrying absorbed sample hashes",
                           MaterialKind.MOLECULAR,
                           state=molecule("hashed nucleus", (HASH, role("accessible chromatin"))))
    graph.react("hash", "Let nuclei absorb the hash oligo on ice, 5 min",
                StepKind.BINDING,
                ((nuclei, InputRole.TEMPLATE), (hash_oligo, InputRole.STRUCTURAL)), (hashed,))
    fixative = graph.input("fixative", "ice-cold formaldehyde fixation buffer (1% final), 15 min on ice",
                           MaterialKind.CONDITION)
    fixed = graph.product("fixed", "fixed hashed nuclei pooled across samples, washed and frozen",
                          MaterialKind.MOLECULAR, state=hashed.state)
    graph.react("fix_pool", "Fix with formaldehyde, pool all samples, wash and flash-freeze",
                StepKind.CROSS_LINKING,
                ((hashed, InputRole.TEMPLATE), (fixative, InputRole.CONDITION)), (fixed,))
    capture = graph.input("capture_oligo",
                          "well-barcoded capture oligo (Nextera P5 handle, N8 UMI, 10-nt well barcode, LNA-marked T30VN), 25 \u00b5M",
                          MaterialKind.OLIGO)
    hash_polymerase = graph.input("hash_polymerase", "NEBNext High-Fidelity 2X PCR Master Mix",
                                  MaterialKind.CATALYST)
    extended = graph.product("extended", "captured hash extended in each first-round well",
                             MaterialKind.MOLECULAR,
                             state=molecule("captured hash", (HASH, UMI, CELL_A,
                                                              role("accessible chromatin"))))
    graph.react("extend_hash", "Distribute 5000 nuclei/well; anneal capture oligo 55 \u00b0C 5 min, ice; extend 55 \u00b0C 10 min",
                StepKind.EXTENSION,
                ((fixed, InputRole.TEMPLATE), (capture, InputRole.STRUCTURAL),
                 (hash_polymerase, InputRole.CATALYST)), (extended,))
    transposome = graph.input("transposome",
                              "indexed Tn5: one of 96 N5 (8) \u00d7 N7 (12) combinations, loaded with Supp. Table 4 oligos and /5Phos/Tn5ME_Rev "
                              "(P5_i5_1, P7_i7_1 and Tn5ME_Rev strands drawn)",
                              MaterialKind.ADAPTER, state=TN5_LOADING_STRANDS,
                              source="Booth et al. 2023 Supp. Table 4 (oligos for loading Tn5)")
    tn5 = graph.input("tn5", "Tn5 transposase (unloaded, Diagenode C01070010)", MaterialKind.CATALYST)
    tagged = graph.product("tagged", "sample-hashed first-round-indexed nuclei",
                           MaterialKind.MOLECULAR,
                           state=molecule("first-round indexed nucleus", (HASH, UMI, CELL_A,
                               role("tagmented accessible DNA"))))
    waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment in DMF tagmentation buffer + CLB, 55 \u00b0C 15 min; stop with EDTA/spermidine",
                StepKind.TAGMENTATION,
                ((extended, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, waste))
    graph.discard(waste, "does not retain a recoverable hash or ATAC insert")
    sort = graph.input("facs_sort", "DAPI FACS sort of a limited number of nuclei per PCR well",
                       MaterialKind.CONDITION)
    sorted_ = graph.product("sorted", "sorted nuclei in second-round PCR wells",
                            MaterialKind.MOLECULAR, state=tagged.state)
    unsorted = graph.product("unsorted", "nuclei not sorted into PCR wells",
                             MaterialKind.WASTE)
    graph.react("partition", "Pool all wells and FACS-sort into reverse-cross-linking buffer",
                StepKind.PARTITION,
                ((tagged, InputRole.TEMPLATE), (sort, InputRole.CONDITION)),
                (sorted_, unsorted), changes_molecule=False)
    graph.discard(unsorted, "not sorted into a PCR well")
    protk = graph.input("proteinase_k", "proteinase K with 0.04% SDS in EB, 65 \u00b0C 13.5 h",
                        MaterialKind.CATALYST)
    released = graph.product("released", "de-cross-linked hash and tagmented DNA per well",
                             MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("reverse_crosslinks", "Reverse cross-links with proteinase K/SDS at 65 \u00b0C",
                StepKind.DEPROTEINIZATION,
                ((sorted_, InputRole.TEMPLATE), (protk, InputRole.CATALYST)), (released,))
    amplified_state = molecule("two-level cell identity", (HASH, UMI, CELL_A, CELL_B,
                               role("accessible-DNA fragment"), SAMPLE))
    fwd = graph.input("pcr_p5", "indexed Nextera P5 primer (well-specific)", MaterialKind.OLIGO)
    rev = graph.input("pcr_p7", "indexed Nextera P7 primer (well-specific)", MaterialKind.OLIGO)
    pcr_pol = graph.input("pcr_polymerase", "NEBNext High-Fidelity 2X PCR Master Mix",
                          MaterialKind.CATALYST)
    amplified = graph.product("amplified", "co-amplified hash and ATAC products from each well",
                              MaterialKind.MOLECULAR, state=amplified_state)
    graph.react("pcr", "PCR both hash and ATAC products: 72 \u00b0C 5 min, 98 \u00b0C 30 s, 22 cycles",
                StepKind.PCR,
                ((released, InputRole.TEMPLATE), (fwd, InputRole.STRUCTURAL),
                 (rev, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(amplified_state,
                    reason="The indexed Nextera P5/P7 PCR primer sequences are not printed in the published two-level Methods."))
    pooled = graph.product("pooled", "pooled, column-concentrated mixed library",
                           MaterialKind.MOLECULAR, state=amplified_state)
    column_flow = graph.product("column_flow", "unbound PCR reagents and short species",
                                MaterialKind.WASTE)
    graph.react("concentrate", "Pool wells and concentrate (Zymo clean and concentrate, 5X binding buffer)",
                StepKind.CLEANUP, ((amplified, InputRole.TEMPLATE),), (pooled, column_flow),
                changes_molecule=False)
    graph.discard(column_flow, "column flow-through and washes")
    atac_final = role_library("sciPlex-ATAC fragment roles",
                              (HASH, CELL_A, CELL_B, role("ATAC insert")), index=SAMPLE)
    hash_final = role_library("sciPlex-ATAC hash-read roles",
                              (HASH, UMI, CELL_A, CELL_B), index=SAMPLE)
    atac = graph.product("atac_library", "ATAC library (gel ~200\u20133000 bp)", MaterialKind.LIBRARY,
                         state=atac_final)
    hashes = graph.product("hash_library", "hash library (199-bp band)", MaterialKind.LIBRARY,
                           state=hash_final)
    gel_rest = graph.product("gel_rest", "gel regions outside the cut windows", MaterialKind.WASTE)
    graph.react("gel_split", "Separate hash and ATAC libraries on a 1% agarose gel",
                StepKind.SPLIT, ((pooled, InputRole.TEMPLATE),), (atac, hashes, gel_rest),
                changes_molecule=False)
    graph.discard(gel_rest, "left in the gel")
    graph.sequence(atac, SequencingPlan(limitation="ATAC library: " + SEQUENCING_UNAVAILABLE))
    graph.sequence(hashes, SequencingPlan(limitation="Hash library: " + SEQUENCING_UNAVAILABLE))
    reviewed = list(graph.materials) + list(graph.steps)
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def protocol_graph():
    return _sciplex_graph(_definition, FINAL_LIBRARY)
