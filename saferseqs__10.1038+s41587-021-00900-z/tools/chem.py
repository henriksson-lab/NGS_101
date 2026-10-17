"""SaferSeqS copied-duplex-UID and strand-specific target library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from batch_ngs import duplex_fragment, seg
from chemdraw import Construct, MolecularState, Segment, Workflow, feature, revcomp
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from molecular_consensus import CopiedDuplexUID

TITLE = "SaferSeqS — copied duplex UID with strand-specific enrichment"
NOTES = "01_saferseqs.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41587-021-00900-z">'
          'Cohen et al., <i>Nature Biotechnology</i> (2021)</a>.')
SUMMARY = ("A 3′-adapter UID is copied into both strands of each original duplex before "
           "amplification. Separate Watson- and Crick-specific nested PCRs preserve strand "
           "identity, allowing paired family consensuses from the same molecule.")
CAVEAT = ("Supplementary Table 7 publishes the adapters, PCR shells and custom sequencing "
          "primers. Target-specific GSPs and eight-base sample indices vary; the exact "
          "libraries below instantiate the published TP53_342 GSP2 (Supplementary Table 5) "
          "and a representative index. The main endpoint is the product of the Table 7 "
          "<i>Crick</i> primer set; the Table 7 <i>Watson</i> product is drawn separately. "
          "The Methods name the strands and read orientation the other way round (see the "
          "source-discrepancy notice).")

UID = feature("saferseqs_duplex_uid", "umi", "random", group="duplex_uid", part="shared")
I5 = feature("saferseqs_i5", "sample_index", "unknown")
I7 = feature("saferseqs_i7", "sample_index", "unknown")
REPRESENTATIVE_INDEX = "ACGTACGT"
# Supplementary Table 5, amplicon TP53_342 GSP2 (printed CTTGAGTTCCAAGGCCTC*A*T).
TP53_342_GSP2 = "CTTGAGTTCCAAGGCCTCAT"
CRICK_TARGET_HANDLE = "CACACAGGAAACAGCTATGACCATG"
WATSON_TARGET_HANDLE = "CGACGTAAAACGACGGCCAGT"
# 3' N14 Adapter Oligo #1 bases between its 5' phosphate and the N14 UID.
ADAPTER_CONSTANT = "ATAAAACGACGGC"
CRICK_ANCHOR_PRIMER = il.P5 + REPRESENTATIVE_INDEX + il.TRUSEQ_READ1
CRICK_GSP2_PRIMER = il.P7 + REPRESENTATIVE_INDEX + CRICK_TARGET_HANDLE + TP53_342_GSP2
CUSTOM_INDEX_PRIMER = "CATGGTCATAGCTGTTTCCTGTGTG"
CUSTOM_READ1_PRIMER = "GAGATCTACAC" + REPRESENTATIVE_INDEX + WATSON_TARGET_HANDLE
CUSTOM_READ2_PRIMER = "GGCATACGAGAT" + REPRESENTATIVE_INDEX + CRICK_TARGET_HANDLE
FINAL_LIBRARY = Construct([
    Segment("P5 graft", il.P5, "p5"),
    Segment("i5 sample index", REPRESENTATIVE_INDEX, "cbc", feature=I5),
    Segment("5′ adapter (TruSeq Read 1 arm)", il.TRUSEQ_READ1, "r1"),
    seg("14-nt copied duplex UID", "N" * 14, "umi", placeholder=True, feature=UID),
    Segment("copied constant tag", revcomp(ADAPTER_CONSTANT), "umi"),
    seg("strand-specific target amplicon", "X" * 42, placeholder=True),
    Segment("TP53_342 GSP2 site", revcomp(TP53_342_GSP2), "r2"),
    Segment("Crick target-primer handle", revcomp(CRICK_TARGET_HANDLE), "r2"),
    Segment("i7 sample index reverse complement", revcomp(REPRESENTATIVE_INDEX), "cbc",
            feature=I7),
    Segment("P7 graft reverse complement", il.P7_RC, "p7"),
], name="SaferSeqS Table 7 Crick-primer-set representative library")
WATSON_LIBRARY = Construct([
    Segment("P5 graft", il.P5, "p5"),
    Segment("i5 sample index", REPRESENTATIVE_INDEX, "cbc", feature=I5),
    Segment("Watson target-primer handle", WATSON_TARGET_HANDLE, "r1"),
    Segment("TP53_342 GSP2", TP53_342_GSP2, "r1"),
    seg("strand-specific target amplicon", "X" * 42, placeholder=True),
    Segment("3′-adapter constant", ADAPTER_CONSTANT, "umi"),
    seg("14-nt duplex UID", "N" * 14, "umi", placeholder=True, feature=UID),
    Segment("3′ adapter (TruSeq Read 2 arm) reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
    Segment("i7 sample index reverse complement", revcomp(REPRESENTATIVE_INDEX), "cbc",
            feature=I7),
    Segment("P7 graft reverse complement", il.P7_RC, "p7"),
], name="SaferSeqS Table 7 Watson-primer-set representative library")
SEQ_PRIMERS = (
    sp.TRUSEQ["R1"],
    sp.custom("Index 1 (i7)", "SaferSeqS custom index primer",
              CUSTOM_INDEX_PRIMER, "Cohen et al. Supplementary Table 7"),
    sp.TRUSEQ["I2"],
    sp.custom("Read 2", "SaferSeqS custom Read 2 primer",
              CUSTOM_READ2_PRIMER, "Cohen et al. Supplementary Table 7"),
)
WATSON_SEQ_PRIMERS = (
    sp.custom("Read 1", "SaferSeqS custom Read 1 primer",
              CUSTOM_READ1_PRIMER, "Cohen et al. Supplementary Table 7"),
    sp.TRUSEQ["I1"],
    sp.TRUSEQ["R2"],
)
WATSON_I5_LIMITATION = ("Index 2 of the Watson-primer-set library: the standard reverse-complement "
                        "i5 primer has no site (the Watson GSP2 handle, not the TruSeq Read 1 arm, "
                        "follows i5) and Supplementary Table 7 lists no custom i5 primer; the "
                        "Methods state eight-base dual indexing on MiSeq, HiSeq 2500 or HiSeq 4000 "
                        "without describing how this i5 was read.")
# Methods "Library preparation": "Paired-end sequencing was performed with eight-base
# dual indexing", so the Watson library's i5 is read although no primer for it is listed.
WATSON_I5_READ = sp.read_spec("Index 2 (i5)", 8, "i5 sample index (Methods: eight-base dual "
                              "indexing; no i5 primer listed for this library)", "demultiplex")
WATSON_CAPTION = ("Product of the Supplementary Table 7 Watson primer set. The custom Read 1 "
                  "primer starts at the GSP2 end; the standard Read 2 primer reads the 14-nt "
                  "UID and then GCCGTCGTTTTAT. The standard Index 1 primer reads i7. "
                  "Supplementary Table 7 provides this custom primer and the Crick-library "
                  "custom primers, all spiked into the standard Illumina primers.")
READ_USES = {"Read 1": "extract_then_include", "Index 1 (i7)": "demultiplex",
             "Index 2 (i5)": "demultiplex", "Read 2": "include"}
WATSON_READ_USES = {"Read 1": "include", "Index 1 (i7)": "demultiplex",
                    "Read 2": "extract_then_include"}
FINAL_CAPTION = ("Product of the Supplementary Table 7 Crick primer set (Crick GSP1 anchor, "
                 "Crick GSP2 primer, Crick GSP2 anchor/indexing primer). Standard Read 1 "
                 "reads the copied 14-nt UID and then the constant tag GCCGTCGTTTTAT; the "
                 "published custom Read 2 primer enters from the target-specific end.")
FINAL_CAPTIONS = {"crick_library": FINAL_CAPTION, "watson_library": WATSON_CAPTION}
SEQUENCING_INTRO = ("Supplementary Table 7 says the custom read, index and Read 2 primers "
                    "are spiked into the standard Illumina primers. The Crick-primer-set "
                    "molecule uses the standard Read 1 and i5 primers and the published "
                    "custom i7 and Read 2 primers. The Watson-primer-set molecule instead "
                    "uses the custom Read 1 primer and standard index and Read 2 primers. "
                    "Methods: eight-base dual indexing on MiSeq, HiSeq 2500 or HiSeq 4000; "
                    "read lengths are not stated.")
SEQUENCING_UNAVAILABLE = ""


def workflow():
    chemistry = CopiedDuplexUID(14)
    w = Workflow(MolecularState("End-repaired genomic fragment",
                                tuple(duplex_fragment("blunt dephosphorylated DNA").rows())))
    uid_end = AdapterLigationGeometry.unavailable(
        "Supplementary Table 7 publishes both adapter oligos, but their two-ended placement "
        "on the genomic duplex is not yet encoded as an exact ligation geometry.")
    w.ligate_adapter("Ligate the UID-bearing 3′ adapter", uid_end,
            chemistry.three_prime_adapter().rows(), name="3′-adapter-ligated duplex",
            note="The adapter is installed at the two 3′ fragment ends after dephosphorylation and blunting.")
    five_prime_end = AdapterLigationGeometry.unavailable(
        "Supplementary Table 7 publishes the 5′ adapter, but its two-ended placement on the "
        "genomic duplex is not yet encoded as an exact ligation geometry.")
    w.ligate_adapter("Extend the UID and ligate the 5′ adapter", five_prime_end,
            chemistry.copied_barcode().rows(), name="duplex-barcoded library",
            note="Extension copies one random UID into both strands before redundant amplification.")
    w.react("Partition into Watson- and Crick-specific hemi-nested PCRs",
            chemistry.strand_families().rows(), name="oriented target families",
            note="A second nested PCR adds sample indices and further target specificity.")
    return w


def oligos():
    return (
        "3′ N14 Adapter Oligo #1: 5′-/5Phos/ATAAAACGACGGCNNNNNNNNNNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAG*T*C/3SpC3/-3′",
        "3′ Adapter Oligo #2: 5′-GCCGUCGUUUUA/33dT/-3′",
        f"5′ Adapter / library amplification primer #1: 5′-{il.TRUSEQ_READ1}-3′",
        "Library amplification primer #2: 5′-GACTGGAGTTCAGACGTGTGCTCTTCCGATCT-3′",
        "Watson GSP1 anchor primer: 5′-GACTGGAGTTCAGACGT*G*T-3′",
        "Crick GSP1 anchor primer: 5′-ACACTCTTTCCCTACACG*A*C-3′",
        "Watson GSP2 primer shell: 5′-AATGATACGGCGACCACCGAGATCTACACXXXXXXXXCGACGTAAAACGACGGCCAGT[GSP2]-3′",
        "Crick GSP2 primer shell: 5′-CAAGCAGAAGACGGCATACGAGATXXXXXXXXCACACAGGAAACAGCTATGACCATG[GSP2]-3′",
        "Watson GSP2 anchor/indexing primer: 5′-CAAGCAGAAGACGGCATACGAGATXXXXXXXXGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT*C*T-3′",
        "Crick GSP2 anchor/indexing primer: 5′-AATGATACGGCGACCACCGAGATCTACACXXXXXXXXACACTCTTTCCCTACACGACGCTCTTCCGAT*C*T-3′",
        "Custom Read 1 sequencing primer: 5′-GAGATCTACACNNNNNNNNCGACGTAAAACGACGGCCA*G*T-3′",
        "Custom index sequencing primer: 5′-CATGGTCATAGCTGTTTCCTGTG*T*G-3′",
        "Custom Read 2 sequencing primer: 5′-GGCATACGAGATNNNNNNNNCACACAGGAAACAGCTATGACCA*T*G-3′",
        "TP53_342 GSP1 (Supplementary Table 5): 5′-GGGCATCCTTGAGTTCC*A*A-3′",
        "TP53_342 GSP2 portion (Supplementary Table 5): 5′-CTTGAGTTCCAAGGCCTC*A*T-3′",
    )



from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan, SourceClaim,
                            SourceDiscrepancy, StepKind)
from protocol_model import molecule, role
from molecular_consensus import uid_segment

def saferseqs_graph(authority: str, final: Construct, limitation: str, *,
                    primers: tuple = (), read_uses: dict | None = None,
                    watson_final: Construct, watson_primers: tuple = (),
                    watson_read_uses: dict | None = None) -> ProtocolGraph:
    graph = ProtocolGraph("SaferSeqS — copied duplex UID with strand-specific enrichment",
                          authority=authority)
    dna = graph.input("genomic_dna", "end-repaired genomic duplex DNA",
                      MaterialKind.BIOLOGICAL,
                      state=molecule("end-repaired genomic DNA", (role("original Watson strand"),
                                                                   role("original Crick strand"))),
                      primary=True)
    user = graph.input("user_enzyme", "USER enzyme (NEB M5505L, 3 U, 15 min 37 °C)",
                       MaterialKind.CATALYST, source=authority + " Methods, Library preparation")
    excised = graph.product("uracil_excised", "uracil-excised genomic DNA",
                            MaterialKind.MOLECULAR, state=dna.state)
    graph.react("user_pretreat", "Pretreat DNA with USER enzyme to excise uracils",
                StepKind.DIGESTION,
                ((dna, InputRole.TEMPLATE), (user, InputRole.CATALYST)), (excised,))
    phosphatase = graph.input("phosphatase", "Accel-NGS 2S end-repair 1/2 reagents",
                              MaterialKind.CATALYST)
    dephosphorylated = graph.product("dephosphorylated", "blunt dephosphorylated duplex DNA",
                                    MaterialKind.MOLECULAR, state=dna.state)
    graph.react("dephosphorylate", "Accel-NGS 2S end repair 1 and 2: dephosphorylate and repair fragment ends (2.0× and 1.8× SPRI)",
                StepKind.DEPHOSPHORYLATION,
                ((excised, InputRole.TEMPLATE), (phosphatase, InputRole.CATALYST)),
                (dephosphorylated,))
    uid_adapter = graph.input("uid_adapter", "custom 50 µM random-UID-bearing 3-prime adapter (replaces kit Y2)",
                              MaterialKind.ADAPTER,
                              source=authority + " Supplementary Table 7")
    ligase3 = graph.input("three_prime_ligase", "3-prime adapter ligase",
                          MaterialKind.CATALYST)
    uid_ligated = graph.product("uid_ligated", "3-prime UID-adapter-ligated duplex",
                                MaterialKind.MOLECULAR, state=CopiedDuplexUID(14).three_prime_adapter())
    graph.react("ligate_uid", "Ligation 1: ligate the UID-bearing 3-prime adapter (1.2× SPRI)",
                StepKind.LIGATION,
                ((dephosphorylated, InputRole.TEMPLATE), (uid_adapter, InputRole.STRUCTURAL),
                 (ligase3, InputRole.CATALYST)), (uid_ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The published oligo pair is transcribed on the protocol page; its "
                    "two-ended placement on the genomic duplex is not yet encoded as an exact geometry."))
    extension = graph.input("extension_polymerase", "UID-copying DNA polymerase",
                            MaterialKind.CATALYST)
    udg = graph.input("udg", "uracil DNA glycosylase (degrades 3-prime Adapter Oligo #2)",
                      MaterialKind.CATALYST,
                      source=authority + " Supplementary Note, Library construction")
    copied = graph.product("uid_copied", "duplex with the UID copied into both strands",
                           MaterialKind.MOLECULAR, state=CopiedDuplexUID(14).copied_barcode())
    graph.react("copy_uid", "Ligation 2 (concerted): extend the annealed 5-prime adapter through the UID while UDG degrades the blocking oligo",
                StepKind.EXTENSION,
                ((uid_ligated, InputRole.TEMPLATE), (extension, InputRole.CATALYST),
                 (udg, InputRole.CATALYST)),
                (copied,))
    five_adapter = graph.input("five_prime_adapter", "custom 42 µM 5-prime adapter (replaces kit B2)",
                               MaterialKind.ADAPTER,
                               source=authority + " Supplementary Table 7")
    ligase5 = graph.input("five_prime_ligase", "5-prime adapter ligase",
                          MaterialKind.CATALYST)
    sealed = graph.product("sealed", "duplex-barcoded library precursor",
                           MaterialKind.MOLECULAR, state=copied.state)
    graph.react("ligate_five_prime", "Ligation 2 (concerted): seal the extended 5-prime adapter to the 5-prime fragment end (1.05× SPRI)",
                StepKind.LIGATION,
                ((copied, InputRole.TEMPLATE), (five_adapter, InputRole.STRUCTURAL),
                 (ligase5, InputRole.CATALYST)), (sealed,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The published 5-prime adapter is transcribed on the protocol page; its "
                    "two-ended placement on the genomic duplex is not yet encoded as an exact geometry."))
    lib_primer1 = graph.input("library_primer_1", "Library amplification primer #1 (2 µM)",
                              MaterialKind.OLIGO, source=authority + " Supplementary Table 7")
    lib_primer2 = graph.input("library_primer_2", "Library amplification primer #2 (2 µM)",
                              MaterialKind.OLIGO, source=authority + " Supplementary Table 7")
    polymerase = graph.input("nested_polymerase", "NEBNext Ultra II Q5 Master Mix",
                             MaterialKind.CATALYST)
    amplified = graph.product("amplified_library", "redundant copies of both original strands",
                              MaterialKind.MOLECULAR, state=sealed.state)
    graph.react("library_pcr", "Library PCR: 5, 7 or 11 cycles (98 °C 10 s, 65 °C 75 s), 1.8× SPRI",
                StepKind.PCR,
                ((sealed, InputRole.TEMPLATE), (lib_primer1, InputRole.STRUCTURAL),
                 (lib_primer2, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.unavailable(
                    final, reason="The library-amplification primers are transcribed, but the "
                    "two-ended adapter geometry of the genomic library is not encoded."))
    watson = graph.product("watson_template", "library aliquot for the Watson-specific PCR",
                           MaterialKind.MOLECULAR, state=sealed.state)
    crick = graph.product("crick_template", "library aliquot for the Crick-specific PCR",
                          MaterialKind.MOLECULAR, state=sealed.state)
    graph.react("orient", "Partition into Watson- and Crick-specific enrichments",
                StepKind.SPLIT, ((amplified, InputRole.TEMPLATE),), (watson, crick),
                changes_molecule=False)
    reviewed = ["genomic_dna", "user_enzyme", "user_pretreat", "phosphatase",
                "dephosphorylate", "uid_adapter",
                "three_prime_ligase", "ligate_uid", "extension_polymerase", "udg", "copy_uid",
                "five_prime_adapter", "five_prime_ligase", "ligate_five_prime",
                "library_primer_1", "library_primer_2", "library_pcr", "orient", "nested_polymerase"]
    gsp1 = graph.input("gsp1", "GSP1 primer (same for both strands; Supplementary Table 5)",
                       MaterialKind.OLIGO, source=authority + " Supplementary Tables 5 and 7")
    reviewed.append(gsp1.id)
    for prefix, template, library_state, seq_primers, uses in (
            ("watson", watson, watson_final, watson_primers, watson_read_uses),
            ("crick", crick, final, primers, read_uses)):
        anchor = graph.input(prefix + "_gsp1_anchor", prefix.title() + " GSP1 anchor primer (Table 7)",
                             MaterialKind.OLIGO, source=authority + " Supplementary Table 7")
        first = graph.product(prefix + "_gsp1_product", prefix.title() + " first-round hemi-nested product",
                              MaterialKind.MOLECULAR, state=sealed.state)
        graph.react(prefix + "_gsp1_pcr", prefix.title() + " first-round anchored PCR: GSP1 + "
                    + prefix.title() + " GSP1 anchor, 19 cycles",
                    StepKind.PCR,
                    ((template, InputRole.TEMPLATE), (gsp1, InputRole.STRUCTURAL),
                     (anchor, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                    (first,), geometry=PCRGeometry.unavailable(library_state, reason=limitation))
        forward = graph.input(prefix + "_forward", prefix.title() + " GSP2 primer (Table 7 shell + Table 5 GSP2)",
                              MaterialKind.OLIGO, source=authority + " Supplementary Tables 5 and 7")
        reverse = graph.input(prefix + "_reverse", prefix.title() + " GSP2 anchor/indexing primer (Table 7)",
                              MaterialKind.OLIGO, source=authority + " Supplementary Table 7")
        library = graph.product(prefix + "_library", prefix.title() + "-primer-set SaferSeqS library",
                                MaterialKind.LIBRARY, state=library_state)
        graph.react(prefix + "_pcr", prefix.title() + " second-round PCR on 1% of round one: GSP2 + anchor/indexing primer, 17 cycles",
                    StepKind.PCR,
                    ((first, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                     (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                    (library,), geometry=PCRGeometry.unavailable(library_state, reason=limitation))
        graph.sequence(library, SequencingPlan(
            library=library_state if seq_primers else None, primers=seq_primers,
            acquired_roles=(("Read 1", "Index 1 (i7)", "Read 2") if prefix == "watson"
                            else ("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2")),
            read_uses=uses or {},
            unplaced_reads=((WATSON_I5_READ,) if prefix == "watson" and seq_primers else ()),
            limitation=(limitation + " " + WATSON_I5_LIMITATION if prefix == "watson"
                        else limitation)))
        reviewed.extend((anchor.id, prefix + "_gsp1_pcr", first.id,
                         forward.id, reverse.id, prefix + "_pcr", library.id))
    graph.annotate(SourceDiscrepancy(
        "watson_pcr", "Watson/Crick primer and read assignment differs within the source",
        (SourceClaim("Supplementary Table 7 (MOESM3, sheet Supplementary Table 7)",
                     "Watson GSP2 anchor/indexing primer = P7 + index + GTGACTGGAG…; Crick "
                     "GSP2 anchor/indexing primer = P5 + index + ACACTCTTTCC…; Watson GSP2 "
                     "shell carries P5, Crick GSP2 shell carries P7."),
         SourceClaim("Methods, Anchored hemi-nested PCR (first round)",
                     "Watson first round uses the P7 short anchor primer; Crick the P5 short "
                     "anchor primer (consistent with Table 7 GSP1 anchors)."),
         SourceClaim("Methods, Anchored hemi-nested PCR (second round)",
                     "Watson second round substitutes the 'anchor P5 indexing primer'; Crick "
                     "the 'anchor P7 indexing primer'."),
         SourceClaim("Methods, Mutation calling and SaferSeqS analysis pipeline",
                     "Watson-read UID = first 14 bases of Read 1; Crick-read UID = first 14 "
                     "bases of Read 2.")),
        "The model uses the Table 7 sequences and names: the UID-bearing 3-prime adapter "
        "contains the complement of the Watson anchor, so the Watson-primer-set product "
        "places the UID next to P7 (read by Read 2) and the Crick-primer-set product places "
        "the copied UID next to P5 (read by Read 1).",
        "A first-round product made with the P7-side Watson anchor cannot be extended by a "
        "P5-side anchor primer, so the Methods second-round names conflict with the first "
        "round and with Table 7; the Methods read assignment matches those second-round names.",
        "Whether the Methods or Table 7 uses the intended Watson/Crick labels is unresolved. "
        "The molecular products drawn are fixed by the sequences; only the strand label "
        "attached to each is uncertain."))
    graph.declare_reviewed_scope(reviewed, authority)
    return graph.require_publishable()


def protocol_graph():
    return saferseqs_graph(SOURCE, FINAL_LIBRARY, "Target-specific and indexed nested-PCR oligos vary by assay; the page models one published TP53 GSP2 with a representative eight-base index.", primers=SEQ_PRIMERS, read_uses=READ_USES,
                           watson_final=WATSON_LIBRARY, watson_primers=WATSON_SEQ_PRIMERS,
                           watson_read_uses=WATSON_READ_USES)
