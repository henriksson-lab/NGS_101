"""AIRR Community bulk UMI 5-prime-RACE BCR library."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Workflow, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "Bulk UMI 5′-RACE AIRR-seq — B-cell receptor repertoires"
NOTES = "01_airr-umi-5-race.html"
SOURCE = ('Authoritative protocol: <a href="https://doi.org/10.1007/978-1-0716-2115-8_19">Gupta et al., AIRR Community (2022)</a>, '
          'Methods Mol Biol 2453, chapter 19 (open access as NCBI Bookshelf NBK586949), describing the Takara SMARTer Human BCR IgG IgM H/K/L Profiling Kit.')
SUMMARY = "Oligo(dT)-primed SMART reverse transcription places a 12-nt UMI and universal handle at the transcript 5′ end. Semi-nested constant-region PCR then captures complete BCR V(D)J regions and adds dual Illumina indices."
CAVEAT = ("The chapter states that “primers from the vendor are not available”: the dT primer, SMART UMI Oligo, PCR1 and PCR2 primer bases "
          "and the complete adapter arms are not printed. Published molecular roles, lengths, index sequences and Read 1/Read 2 orientation are shown "
          "without inventing bases or primer landing sites. Figure 5 labels the post-UMI linker as 4 bp followed by an unsized "
          "template-switch block, while the text trims 12 + 7 = 19 Read 2 bases; the model keeps the published 7-base linker/template-switch total.")
UMI_LENGTH = 12
SMART_LINKER_LENGTH = 7  # "-len 7": linker + template switch (Subheading 3.10 step 5)
C_PRIMER_SPAN = "27–30 bp"  # Fig. 5: "C-region primer + 7N"
# Section 2.2.2: hBCR PCR2 Universal Forward 1-12 carry these Illumina i7 indexes.
I7_INDEXES = (("D701", "ATTACTCG"), ("D702", "TCCGGAGA"), ("D703", "CGCTCATT"),
              ("D704", "GAGATTCC"), ("D705", "ATTCAGAA"), ("D706", "GAATTCGT"),
              ("D707", "CTGAAGCT"), ("D708", "TAATGCGC"), ("D709", "CGGCTATG"),
              ("D710", "TCCGCGAA"), ("D711", "TCTCGCGC"), ("D712", "AGCGATAG"))
# Section 2.2.3: hBCR PCR2 IgG/IgM/IgK/IgL Reverse 1-4 carry these i5 indexes,
# "as read on a MiSeq instrument" (read strings, not placed as library bases).
I5_INDEXES = (("D501", "TATAGCCT"), ("D502", "ATAGAGGC"), ("D503", "CCTATCCT"),
              ("D504", "GGCTCTGA"))

FINAL_LIBRARY = Construct([
    seg("P5 arm", "[P5 arm]", "p5", placeholder=True),
    seg("eight-base i5 sample index", "J" * 8, "cbc", placeholder=True,
        feature=feature("sample_index_i5", "sample_index", "whitelist",
                        whitelist="chapter section 2.2.3: D501-D504 on hBCR PCR2 Ig reverse 1-4 (MiSeq read strings)")),
    seg("Read 1 entry site", "[Read 1 site]", "r1", placeholder=True),
    seg("C-region primer + 7N (27–30 bp, Fig. 5)", "C" * 30, placeholder=True),
    seg("constant/CDR3/V(D)J/leader region (length variable)", "V" * 300, placeholder=True),
    seg("seven-base linker + template-switch sequence", "L" * SMART_LINKER_LENGTH, placeholder=True),
    seg("twelve-base molecular UMI", "U" * UMI_LENGTH, "umi", placeholder=True,
        feature=feature("bcr_molecule_umi", "umi", "random")),
    seg("Read 2 entry site", "[Read 2 site]", "r2", placeholder=True),
    seg("eight-base i7 sample index", "I" * 8, "cbc", placeholder=True,
        feature=feature("sample_index_i7", "sample_index", "whitelist",
                        whitelist="chapter section 2.2.2: D701-D712 on hBCR PCR2 Universal Forward 1-12")),
    seg("P7 arm", "[P7 arm]", "p7", placeholder=True),
], name="dual-index bulk BCR AIRR-seq library")
SEQ_PRIMERS = ()
FINAL_CAPTION = ("Read 1 starts in the BCR constant-region side (C-region primer + 7N, 27–30 bp in Fig. 5); Read 2 starts at the SMART end, "
                 "reads the 12-base UMI and then the seven-base linker/template-switch sequence before the leader and V region. "
                 "Fig. 5 draws each adapter end (P5 + i5 + Read 1, P7 + i7 + Read 2) as about 70 bp; their bases are not printed.")
SEQUENCING_ENDING = ("The chapter specifies MiSeq Reagent Kit v3, paired-end 2 × 300 with a 10% PhiX spike-in. Paired reads cover constant-region identity and the rearranged variable region. "
                     "The first 19 Read 2 bases are the 12-nt UMI followed by the seven-base linker/template-switch sequence. "
                     "Dual eight-base indexes (D701–D712 i7 on the universal forward primers, D501–D504 i5 on the constant-region reverse primers) identify the sample.")
SEQUENCING_UNAVAILABLE = "The source explicitly says the vendor primers are unavailable. It gives the read orientation and index/UMI lengths, but not exact sequencing-primer sequences from which binding could be computed."

_initial = Scene(); _initial.strand("mRNA", [seg("BCR transcript", "R" * 54, placeholder=True), seg("poly(A)", "A" * 18)], label="polyadenylated BCR mRNA")
INITIAL_ROWS = tuple(_initial.rows()); INITIAL_NAME = "Bulk BCR mRNA"


def workflow():
    cdna = Scene(); cdna.strand("cDNA", [seg("dT primer (sequence not printed)", "X" * 14, placeholder=True),
                                         seg("BCR cDNA", "D" * 54, placeholder=True),
                                         seg("template-switch/linker complement", "L" * SMART_LINKER_LENGTH, placeholder=True),
                                         seg("twelve-base UMI", "U" * UMI_LENGTH, "umi", placeholder=True,
                                             feature=feature("bcr_molecule_umi", "umi", "random")),
                                         seg("SMART universal handle", "L" * 18, placeholder=True)], label="template-switched first-strand cDNA")
    cdna.mark("cDNA", "twelve-base UMI", "one UMI per parental cDNA")
    pcr1 = Scene.duplex([seg("universal SMART end", "L" * 18, placeholder=True),
                         seg("complete V(D)J", "V" * 54, placeholder=True),
                         seg("constant-region segment", "C" * 32, placeholder=True)], label="PCR1 BCR amplicon")
    pcr1.mark("top", "constant-region segment", "IgG / IgM / IgK / IgL reverse primer")
    final = Scene.duplex(list(FINAL_LIBRARY), label="semi-nested PCR2 library",
                         unpaired=tuple(s.name + "'" for s in FINAL_LIBRARY if s.is_role_token()))
    w = Workflow(MolecularState(INITIAL_NAME, INITIAL_ROWS))
    w.react("Prime reverse transcription from poly(A)", cdna.rows(), name="UMI-tagged first-strand cDNA",
            note="SMARTScribe adds non-templated bases; the SMART UMI oligo template-switches and copies a 12-nt UMI plus universal handle.")
    w.react("PCR1 from universal handle to constant region", pcr1.rows(), name="Chain-specific PCR1 amplicon",
            note="Separate IgG, IgM, kappa and lambda reverse primers retain complete V(D)J plus part of the constant region.")
    w.react("Semi-nested PCR2 and dual indexing", final.rows(), name="Sequencing library",
            note="PCR2 uses chain-specific nested reverse primers and indexed universal forward primers to add P5, P7, i5 and i7.")
    return w


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    rna=graph.input("bcr_mrna","polyadenylated BCR mRNA",MaterialKind.MOLECULAR,
                    primary=True)
    oligo_dt=graph.input("oligo_dt","proprietary oligo-dT RT primer",MaterialKind.OLIGO,
                         inferred=True,source=SOURCE+" (vendor sequence unavailable)")
    rt=graph.input("reverse_transcriptase","SMARTScribe reverse transcriptase (42 °C 90 min, 70 °C 10 s)",
                   MaterialKind.CATALYST)
    first=graph.product("first_strand","oligo-dT-primed BCR cDNA",MaterialKind.MOLECULAR,
                        inferred=True,source=SOURCE+" (vendor handle bases unavailable)")
    graph.react("reverse_transcribe","Denature RNA with dT primer and BCR enhancer (72 °C 3 min), then prime RT from poly(A)",
                StepKind.REVERSE_TRANSCRIPTION,
                ((rna,InputRole.TEMPLATE),(oligo_dt,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(first,))
    smart=graph.input("smart_umi","12-nt-UMI SMART template-switch oligo",
                      MaterialKind.OLIGO,inferred=True,
                      source=SOURCE+" (vendor sequence unavailable)")
    switched=graph.product("switched","UMI-tagged template-switched BCR cDNA",
                           MaterialKind.MOLECULAR)
    graph.react("template_switch","Install the UMI and universal 5-prime handle",
                StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(smart,InputRole.STRUCTURAL),
                 (rt,InputRole.CATALYST)),(switched,))
    universal=graph.input("universal_primer","universal SMART PCR1 primer",
                          MaterialKind.OLIGO,inferred=True,
                          source=SOURCE+" (vendor sequence unavailable)")
    constant=graph.input("constant_primer","chain-specific constant-region PCR1 primer pool",
                         MaterialKind.OLIGO,inferred=True,
                         source=SOURCE+" (vendor sequences unavailable)")
    polymerase=graph.input("pcr_polymerase","PrimeSTAR GXL DNA polymerase",
                           MaterialKind.CATALYST)
    pcr1=graph.product("pcr1_product","chain-specific complete-V(D)J PCR1 amplicons",
                       MaterialKind.MOLECULAR)
    graph.react("pcr1","PCR1 (18 cycles, one IgG/IgM/IgK/IgL reverse primer per reaction)",
                StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(universal,InputRole.STRUCTURAL),
                 (constant,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (pcr1,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The authoritative protocol withholds the SMART and gene-specific primer sequences."))
    nested_f=graph.input("nested_forward","indexed universal PCR2 forward primer",
                         MaterialKind.OLIGO,inferred=True,
                         source=SOURCE+" (vendor sequence unavailable)")
    nested_r=graph.input("nested_reverse","semi-nested constant-region PCR2 reverse primer pool",
                         MaterialKind.OLIGO,inferred=True,
                         source=SOURCE+" (vendor sequences unavailable)")
    amplified=graph.product("pcr2_product","indexed PCR2 BCR amplicons",MaterialKind.MOLECULAR,
                          state=FINAL_LIBRARY,inferred=True,
                          source=SOURCE+" (adapter/primer bases not printed)")
    graph.react("pcr2","Semi-nested PCR2 adds P7/i7 (forward) and P5/i5 (reverse); 16–25 cycles per Table 1",StepKind.PCR,
                ((pcr1,InputRole.TEMPLATE),(nested_f,InputRole.STRUCTURAL),
                 (nested_r,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The chapter discloses molecular roles, lengths and index sequences but not the complete PCR2 oligo bases."))
    beads=graph.input("nucleomag","NucleoMag NGS beads: 25 µL then 10 µL to the supernatant",
                      MaterialKind.SELECTION)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY,inferred=True,
                          source=SOURCE+" (adapter/primer bases not printed)")
    first_beads=graph.product("first_bead_fraction","material bound by the first 25 µL bead addition",
                              MaterialKind.WASTE)
    final_sup=graph.product("second_supernatant","supernatant left unbound after the 10 µL bead addition",
                            MaterialKind.WASTE)
    graph.react("size_select","Two-sided bead size selection; elute library in 17 µL",
                StepKind.SELECTION,((amplified,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),
                (library,first_beads,final_sup))
    graph.discard(first_beads,"first bead pellet is discarded; the library stays in its supernatant")
    graph.discard(final_sup,"second-binding supernatant is discarded; the library is bead-bound")
    graph.sequence(library,SequencingPlan(
        limitation=SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope((
        "bcr_mrna","oligo_dt","reverse_transcriptase","reverse_transcribe",
        "smart_umi","template_switch","universal_primer","constant_primer",
        "pcr_polymerase","pcr1","nested_forward","nested_reverse","pcr2",
        "pcr2_product","nucleomag","size_select","library"),SOURCE)
    return graph.require_publishable()
