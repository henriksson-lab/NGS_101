"""10x Chromium Fixed RNA Profiling / Flex probe-pair chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role

seg=Segment

TITLE="10x Fixed RNA Profiling (Flex)"
NOTES="01_flex.html"
SOURCE='Commercial protocol: 10x Genomics <a href="https://www.10xgenomics.com/support/cn/flex-gene-expression/documentation/steps/library-prep/chromium-fixed-rna-profiling-reagent-kits-for-singleplexed-samples-with-feature-barcode-technology-for-protein-using-barcode-oligo-capture">Fixed RNA Profiling user guide CG000674 Rev C</a>.'
SUMMARY="Gene-specific probe pairs hybridize adjacently on fixed RNA. In GEMs, their nick is sealed and the reporter receives a cell barcode and UMI. Antibody-conjugated oligos on the same fixed cells are captured by the gel-bead primer and yield a separate protein-expression library."
CAVEAT="CG000674 Rev C Appendix p. 85 prints both completed library layouts; CG000529 Rev C p. 4 prints the antibody capture oligo pattern, and CG000621 Rev D Tables 7 and 13 print physical custom RNA-probe architecture. The RNA library shell is drawn from p. 85 with the Small RNA Read 2 (Read 2S) arm; its Read 2S and i7 primer sites are derived from that printed constant, not from a printed primer. The LHS/RHS drawing follows the CG000621 custom architecture. 10x publishes commercial WTA target-region sequences through Cell Ranger; full physical LHS/RHS oligos and their barcode-mixture assignment have not been reconstructed, so the generic RNA probe region remains inferred. The protein shell is source-based, but its AB8, antibody barcode and sample index values vary by reagent or run."

def oligos():
    return [
        "<b>Flex gel-bead primer</b> — 5′-CTACACGACGCTCTTCCGATCT-[16-base cell barcode]-[12-base UMI]-TTGCTAGGACCG-3′",
        "<b>Published gene-expression reporter junction</b> — 5′-TTGCTAGGACCG-[BC8]-NN-TACGTGCTAACCGCGT-[ligated probe insert]-TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC-3′",
        "<b>Published custom Flex probe design</b> (<a href=\"https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf\">CG000621 Rev D Table 7</a>; not the commercial whole-transcriptome inventory) — LHS 5′-CCTTGGCACCCGAGAATTCCA-[25-nt target]-3′; RHS /5Phos/-[25-nt target]-ACGCGGTTAGCACGTA-NN-ACTTTAGG-CGGTCCTAGCAA-3′ for singleplex BC001",
        "<b>Official WTA reference example</b> — human v1.0.1 TSPAN6 spliced target-region probe_seq: GGTGACACCACAACAATGCAACGTATTTTGGATCTTGTCTACTGCATGGC (50 nt; reference sequence, not a complete physical oligo or a specimen-specific choice)",
        "<b>Protein antibody capture oligo</b> — <a href=\"https://cdn.10xgenomics.com/image/upload/v1719610805/support-documents/CG000529_DemonstratedProtocol_CellSurfaceProteinLabelingforFixedRNA_RevC.pdf\">CG000529 Rev C p. 4</a> prints /5AmMC12/CGGAGATGTGTATAAGAGACAG-N10-N15-N9-CCCATATAAGAAA for the TotalSeq-C Barcode Oligo Capture branch (N15 is the antibody feature barcode)",
    ]

# Flex probe-pair constants.  CG000621 Rev D Table 7 (custom singleplex probe design)
# prints LHS 5'-CCTTGGCACCCGAGAATTCCA-[target]-3' and RHS
# /5Phos/-[target]-ACGCGGTTAGCACGTA-NN-[BC]-CGGTCCTAGCAA-3'.  CG000674 Rev C Appendix
# p. 85 prints the same constants (as their reverse complements) in the finished
# gene-expression library, so the library below is assembled from these declarations.
LHS_CONSTANT = revcomp(il.TRUSEQ_SMALL_RNA_RA3)          # CCTTGGCACCCGAGAATTCCA
RHS_CONSTANT = "ACGCGGTTAGCACGTA"
CAPTURE_SEQUENCE = "CGGTCCTAGCAA"
GEL_BEAD_TERMINUS = "TTGCTAGGACCG"
if revcomp(CAPTURE_SEQUENCE) != GEL_BEAD_TERMINUS:
    raise ValueError("Flex gel-bead terminus must pair with the RHS capture sequence")
# CG000674 Rev C Appendix p. 85: Read 2S-side constant between the probe insert and i7.
SMALL_RNA_READ2_ARM = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
if not SMALL_RNA_READ2_ARM.startswith(il.TRUSEQ_SMALL_RNA_RA3):
    raise ValueError("Flex Read 2S arm must begin with the LHS-constant complement")
TARGET_ARM_NT = 25   # CG000621 Rev D: 25-nt LHS and RHS target regions (50-nt probe_seq)
PROBE_BARCODE = feature("probe_barcode", "sample_index", "unknown",
                        note="Probe Barcode BC001–BC016; CG000621 Rev D says each commercial "
                             "WTA barcode ID is a mixture of four physical sequences")


def probe_pair_scene() -> Scene:
    """Flex LHS/RHS on target RNA; the Scene enforces probe/RNA pairing.

    Probe architecture follows CG000621 Rev D Table 7 (custom design); the commercial
    WTA physical oligos are not reconstructed, so target arms and BC8 stay placeholders.
    """
    sc = Scene()
    sc.strand("RNA", [seg("RHS-binding target RNA", "X" * TARGET_ARM_NT, placeholder=True,length_bp=TARGET_ARM_NT),
                      seg("LHS-binding target RNA", "X" * TARGET_ARM_NT, placeholder=True,length_bp=TARGET_ARM_NT)],
              label="fixed target RNA")
    sc.anneal("LHS", [seg("LHS constant (Read 2S side)", LHS_CONSTANT, "r2"),
                      seg("LHS target arm", "X" * TARGET_ARM_NT, placeholder=True, inferred=True,length_bp=TARGET_ARM_NT)],
              to="RNA", pair=("LHS target arm", "LHS-binding target RNA"),
              label="LHS probe", unpaired=("LHS constant (Read 2S side)",))
    sc.anneal("RHS", [seg("RHS target arm", "X" * TARGET_ARM_NT, placeholder=True, inferred=True,length_bp=TARGET_ARM_NT),
                      seg("RHS constant", RHS_CONSTANT),
                      seg("RHS NN", "NN", placeholder=True,length_bp=2),
                      seg("probe barcode", "B" * 8, "cbc", placeholder=True, inferred=True,
                          feature=PROBE_BARCODE),
                      seg("capture sequence", CAPTURE_SEQUENCE)],
              to="RNA", pair=("RHS target arm", "RHS-binding target RNA"),
              label="RHS probe", mod5="p",
              unpaired=("RHS constant", "RHS NN", "probe barcode", "capture sequence"))
    sc.arrow("LHS", "LHS 3′ end ligates to the 5′-phosphorylated RHS inside the GEM")
    sc.labels("LHS"); sc.labels("RHS")
    return sc


LIB = Construct([
    seg("P5", il.P5, "p5"),
    seg("i5", "J"*10, "cbc", placeholder=True, inferred=True,
        feature=feature("sample_index_i5", "sample_index", "unknown")),
    seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
    seg("cell barcode", "B"*16, "cbc", placeholder=True, inferred=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="10x Flex whitelist")),
    seg("UMI", "U"*12, "umi", placeholder=True, inferred=True,
        feature=feature("umi", "umi", "random")),
    seg("gel-bead capture terminus", GEL_BEAD_TERMINUS),
    seg("probe barcode (BC8)", "B"*8, "cbc", placeholder=True, inferred=True, feature=PROBE_BARCODE),
    seg("RHS NN complement", "NN", placeholder=True,length_bp=2),
    seg("RHS constant complement", revcomp(RHS_CONSTANT)),
    seg("ligated probe insert", "X"*(2*TARGET_ARM_NT), placeholder=True, inferred=True),
    seg("Small RNA Read 2 arm", SMALL_RNA_READ2_ARM, "r2"),
    seg("i7", "I"*10, "cbc", placeholder=True, inferred=True,
        feature=feature("sample_index_i7", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="Flex gene-expression library")
_R2S_NOTE = ("Sequence derived from the Read 2S-side constant printed in CG000674 Rev C Appendix "
             "p. 85; the guide names Small RNA Read 2 but does not print the primer itself.")
PRIMERS = (sp.TRUSEQ["R1"],
           sp.custom("Index 1 (i7)", "Small RNA Index 1 (derived)", SMALL_RNA_READ2_ARM,
                     "10x CG000674 Rev C Appendix p. 85", _R2S_NOTE),
           sp.TRUSEQ["I2"],
           sp.custom("Read 2", "Small RNA Read 2 (Read 2S, derived)", revcomp(SMALL_RNA_READ2_ARM),
                     "10x CG000674 Rev C Step 6 and Appendix p. 85", _R2S_NOTE))
# CG000674 Rev C Step 6 Sequencing Parameters (p. 73): 28/10/10/90 cycles, both libraries.
READ_LENGTHS = {"Read 1": 28, "Index 1 (i7)": 10, "Index 2 (i5)": 10, "Read 2": 90}
_ge_issues = sp.verify(LIB, PRIMERS, required_roles=tuple(p.role for p in PRIMERS))
if _ge_issues:
    raise ValueError("Flex gene-expression library primer geometry: " + "; ".join(_ge_issues))
ANTIBODY_OLIGO = Construct([
    seg("Nextera partial Read 2", "CGGAGATGTGTATAAGAGACAG", "r2", note="5′ /5AmMC12/ antibody conjugation in CG000529 Rev C"),
    seg("antibody N10", "N"*10, placeholder=True,length_bp=10),
    seg("antibody feature barcode", "N"*15, "cbc", placeholder=True, inferred=True,
        feature=feature("antibody_feature_barcode", "feature_barcode", "unknown")),
    seg("antibody N9", "N"*9, placeholder=True,length_bp=9),
    seg("antibody capture tail", "CCCATATAAGAAA"),
], name="TotalSeq-C antibody oligo for Flex Barcode Oligo Capture")
PROTEIN_LIB = Construct([
    seg("P5", il.P5, "p5"),
    seg("i5", "J"*10, "cbc", placeholder=True, inferred=True,
        feature=feature("protein_i5", "sample_index", "unknown")),
    seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
    seg("GEM barcode", "B"*16, "cbc", placeholder=True, inferred=True,
        feature=feature("protein_cell_barcode", "cell_barcode", "whitelist", whitelist="10x Flex whitelist")),
    seg("UMI", "U"*12, "umi", placeholder=True, inferred=True,
        feature=feature("protein_umi", "umi", "random")),
    seg("gel-bead capture terminus", GEL_BEAD_TERMINUS),
    seg("AB8 antibody multiplex barcode", "A"*8, "cbc", placeholder=True, inferred=True,
        feature=feature("protein_ab8", "sample_index", "unknown",
                        note="Antibody Multiplexing Barcode identifies the sample; CG000674 Rev C p. 20")),
    seg("protein N5", "N"*5, placeholder=True,length_bp=5),
    seg("protein constant 1", "AACCGCGTATTCAGGT"),
    seg("protein constant 2", "TTTCTTATATGGG"),
    seg("protein N9", "N"*9, placeholder=True,length_bp=9),
    seg("antibody feature barcode", "N"*15, "cbc", placeholder=True, inferred=True,
        feature=feature("antibody_feature_barcode", "feature_barcode", "unknown")),
    seg("protein N10", "N"*10, placeholder=True,length_bp=10),
    seg("Nextera Read 2 reverse-complement arm", nx.INDEX1_PRIMER, "r2"),
    seg("i7", "I"*10, "cbc", placeholder=True, inferred=True,
        feature=feature("protein_i7", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="Flex protein-expression library")
PROTEIN_PRIMERS = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"], sp.TRUSEQ["I2"], sp.NEXTERA["R2"])
_protein_issues = sp.verify(PROTEIN_LIB, PROTEIN_PRIMERS,
                            required_roles=tuple(p.role for p in PROTEIN_PRIMERS))
if _protein_issues:
    raise ValueError("Flex protein library primer geometry: " + "; ".join(_protein_issues))
FINAL_LIBRARIES=(("Final probe-derived gene-expression library", LIB, PRIMERS,
                  "INFERRED — The shell follows CG000674 Rev C Appendix p. 85; the 50-nt probe insert, probe barcode, GEM barcode, UMI and selected 10-base indexes are placeholders.",
                  "TruSeq Read 1 and the Small RNA Read 2 (Read 2S) site derived from the printed library constant."),
                 ("Final antibody-derived protein-expression library", PROTEIN_LIB, PROTEIN_PRIMERS,
                  "INFERRED — The guide prints this shell and the antibody oligo pattern, but the AB8 and N15 barcode values and selected 10-base indexes are reagent- or run-dependent.",
                  "TruSeq Read 1 and Nextera Read 2 primer sites are computed on the separate protein library."))

def gel_bead_primer():
    return Construct([seg("partial Read1T",il.TRUSEQ_READ1[-22:],"r1"),
        seg("cell barcode","B"*16,"cbc",placeholder=True,feature=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x Flex whitelist")),
        seg("UMI","U"*12,"umi",placeholder=True,feature=feature("umi","umi","random")),
        seg("gel-bead capture terminus",GEL_BEAD_TERMINUS)],name="CG000674 RevC gel-bead primer")


def barcoded_rna_scene():
    # Single-stranded primer extension product. Complete preamplification/end
    # architectures are not recovered from the finished indexed library.
    sc=Scene();sc.strand("barcoded reporter",[*gel_bead_primer(),
        seg("probe barcode (BC8)","B"*8,"cbc",placeholder=True,inferred=True,feature=PROBE_BARCODE),
        seg("RHS NN complement","NN",placeholder=True,length_bp=2),
        seg("RHS constant complement",revcomp(RHS_CONSTANT)),
        seg("copied probe target region","X"*50,placeholder=True,inferred=True),
        seg("LHS constant complement",revcomp(LHS_CONSTANT),"r2")],
        label="single-stranded barcoded probe reporter (custom-architecture reference)")
    return sc


def reporter_mixture(label, *, amplified):
    sc=Scene()
    sc.strand("probe population",[role("preamplified probe reporters; physical ends unassigned" if amplified else "barcoded probe reporters")],label=label+": probe population")
    sc.strand("protein population",[role("preamplified antibody reporters; physical ends unassigned" if amplified else "barcoded antibody reporters")],label=label+": antibody population")
    return sc


def cell_contents(label, *, antibody_labeled):
    """Separate RNA and surface-associated material in the same cell container."""
    sc=Scene()
    sc.strand("target RNA",[role("target RNA")],label=label+": RNA population")
    sc.strand("surface material",[role("bound antibody oligo" if antibody_labeled else "surface protein")],
              label=label+": surface-associated material")
    return sc


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cell=g.input("cell","unfixed cell carrying target RNA and surface proteins",MaterialKind.BIOLOGICAL,
                 state=cell_contents("cell before Flex labeling",antibody_labeled=False),primary=True)
    antibody=g.input("antibody_oligo","TotalSeq-C antibody-conjugated Barcode Oligo Capture oligo",MaterialKind.SELECTION,
                     state=ANTIBODY_OLIGO, inferred=True,
                     source="10x CG000529 Rev C p. 4; N15 antibody identity depends on the conjugate")
    labeled=g.product("labeled","cell labeled with antibody-oligonucleotide conjugate",MaterialKind.MOLECULAR,
                      state=cell_contents("antibody-labeled cell",antibody_labeled=True))
    unbound_antibody=g.product("unbound_antibody","unbound antibody-oligonucleotide conjugate",MaterialKind.WASTE)
    g.react("label_antibody","Label surface protein with antibody-oligonucleotide conjugate",StepKind.BINDING,
            ((cell,InputRole.TEMPLATE),(antibody,InputRole.SELECTION)),(labeled,unbound_antibody),
            source="10x CG000529 Rev C pp. 4–5; CG000478 Rev C fixation protocol")
    g.discard(unbound_antibody,"unbound antibody conjugates are removed before fixed-cell processing")
    fixative=g.input("fixative","Flex cell fixation condition",MaterialKind.CONDITION,
                     source="10x CG000478 Rev C fixation protocol")
    fixed=g.product("fixed_labeled_cell","fixed antibody-labeled cell with target RNA",MaterialKind.MOLECULAR,
                    state=cell_contents("fixed labeled cell",antibody_labeled=True))
    g.react("fix","Fix antibody-labeled cells",StepKind.CROSS_LINKING,
            ((labeled,InputRole.TEMPLATE),(fixative,InputRole.CONDITION)),(fixed,),
            source="10x CG000478 Rev C; optional protein labeling precedes fixation")
    left=g.input("left_probe","gene-specific left probe half (custom-architecture reference)",MaterialKind.OLIGO,
        state=Construct(probe_pair_scene().strands["LHS"].segs),inferred=True,source=SOURCE+" (panel-specific physical oligo not reconstructed)")
    right=g.input("right_probe","gene-specific right probe half (custom-architecture reference)",MaterialKind.OLIGO,
        state=Construct(probe_pair_scene().strands["RHS"].segs),inferred=True,source=SOURCE+" (panel-specific physical oligo not reconstructed)")
    ab_oligo=g.input("ab_oligo","selected Antibody Multiplexing Barcode oligo AB001–AB016",MaterialKind.OLIGO,
                     inferred=True,source="10x CG000674 Rev C Probe Hybridization step 1.1h; individual physical oligo not reconstructed")
    adjacent=g.product("adjacent","adjacent probe halves hybridized on target RNA",MaterialKind.MOLECULAR,state=probe_pair_scene())
    retained_antibody=g.product("retained_antibody","antibody-conjugated capture oligo on the hybridized cell",MaterialKind.MOLECULAR,
                                state=ANTIBODY_OLIGO)
    retained_ab=g.product("retained_ab","selected Antibody Multiplexing Barcode oligo in probe-hybridized suspension",MaterialKind.MOLECULAR,
                          inferred=True,source="10x CG000674 Rev C step 1.1h; AB001–AB016 sequence choice is sample-specific")
    g.react("hybridize","Hybridize adjacent probe halves to fixed RNA",StepKind.HYBRIDIZATION,
            ((fixed,InputRole.TEMPLATE),(left,InputRole.STRUCTURAL),(right,InputRole.STRUCTURAL),(ab_oligo,InputRole.STRUCTURAL)),
            (adjacent,retained_antibody,retained_ab),
            source="10x CG000674 Rev C Step 1 Probe Hybridization overview p. 15 (overnight; Antibody Multiplexing Barcodes added simultaneously)")
    wash=g.input("post_hyb_wash","post-hybridization wash",MaterialKind.CONDITION,
                 source="10x CG000674 Rev C Step 2 overview p. 15")
    washed=g.product("washed","washed probe-hybridized cell with adjacent probe halves",MaterialKind.MOLECULAR,state=adjacent.state)
    washed_antibody=g.product("washed_antibody","antibody-conjugated capture oligo after post-hybridization wash",MaterialKind.MOLECULAR,
                              state=ANTIBODY_OLIGO)
    washed_ab=g.product("washed_ab","Antibody Multiplexing Barcode oligo retained after post-hybridization wash",MaterialKind.MOLECULAR,
                        inferred=True,source="10x CG000674 Rev C; retention mechanism of the AB oligo through washing not described in the overview")
    unbound_probes=g.product("unbound_probes","unbound LHS/RHS probes in wash supernatant",MaterialKind.WASTE)
    g.react("wash_unbound","Wash off unbound probes after hybridization",StepKind.CLEANUP,
            ((adjacent,InputRole.TEMPLATE),(retained_antibody,InputRole.STRUCTURAL),(retained_ab,InputRole.STRUCTURAL),
             (wash,InputRole.CONDITION)),(washed,washed_antibody,washed_ab,unbound_probes),changes_molecule=False,
            source="10x CG000674 Rev C p. 15: 'After hybridization, the unbound probes are washed off.'")
    g.discard(unbound_probes,"removed in the post-hybridization wash supernatant")
    droplets=g.input("droplets","Flex GEM partitioning condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","partitioned RNA-bound adjacent probe halves",MaterialKind.MOLECULAR,state=adjacent.state)
    occupied_antibody=g.product("occupied_antibody","partitioned antibody-conjugated oligo",MaterialKind.MOLECULAR,
                                state=ANTIBODY_OLIGO)
    occupied_ab=g.product("occupied_ab","partitioned Antibody Multiplexing Barcode oligo",MaterialKind.MOLECULAR,
                          inferred=True,source="10x CG000674 Rev C GEM Generation; physical AB oligo not reconstructed")
    failed=g.product("failed","empty or unresolved GEM material",MaterialKind.MOLECULAR)
    g.react("partition","Partition hybridized RNA probes and antibody conjugates with barcoded gel beads",StepKind.PARTITION,
            ((washed,InputRole.TEMPLATE),(washed_antibody,InputRole.STRUCTURAL),(washed_ab,InputRole.STRUCTURAL),
             (droplets,InputRole.CONDITION)),(occupied,occupied_antibody,occupied_ab,failed),changes_molecule=False)
    g.mark_terminal(Terminal(failed.id,Disposition.UNRESOLVED,
        "GEM products are recovered together. Barcode resolution is computational and does not establish physical removal of empty/multiplet GEM classes."))
    ligase=g.input("ligase","probe-pair ligase in GEM",MaterialKind.CATALYST,source="10x CG000674 Rev C GEM Generation p. 16; enzyme identity not named")
    reporter=g.product("reporter","ligated probe-pair reporter in GEM",MaterialKind.MOLECULAR,state=molecule("ligated Flex reporter",(role("LHS constant"),role("LHS target arm"),role("RHS target arm"),role("probe barcode"),role("capture sequence")),junctions=(("LHS target arm","RHS target arm","probe-pair ligation"),)))
    g.react("seal_probe_pair","Ligate the adjacent probe pair inside the GEM",StepKind.NICK_SEALING,((occupied,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(reporter,),
            source="10x CG000674 Rev C p. 16: ligation seals the LHS/RHS nick while probes remain hybridized")
    bead=g.input("bead_oligo","cell-barcode/UMI gel-bead oligo CTACACGACGCTCTTCCGATCT-N16-N12-TTGCTAGGACCG",MaterialKind.OLIGO,state=gel_bead_primer(),source=SOURCE+" Appendix p. 85")
    pol=g.input("extension_polymerase","barcode-transfer polymerase",MaterialKind.CATALYST)
    barcoded=g.product("barcoded","cell-barcode/UMI-tagged probe reporters",MaterialKind.MOLECULAR,
        state=barcoded_rna_scene(),inferred=True,source="CG000674 RevC gel-bead/final constants and CG000621 RevD custom-probe architecture; selected physical target/barcode values unassigned")
    g.react("barcode","INFERRED — extend cell-barcode primer on a ligated probe reporter reference",StepKind.EXTENSION,((reporter,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(barcoded,),
            source="10x CG000674 Rev C p. 16: Gel Bead primer hybridizes to the capture sequence and is extended")
    protein_joined=g.product("protein_joined","antibody feature oligo attached to selected AB oligo inside GEM",MaterialKind.MOLECULAR,
                             inferred=True,source="10x CG000674 Rev C Inside Individual GEMs; bond and physical AB oligo not reconstructed")
    g.react("attach_ab","Attach Antibody Multiplexing Barcode to protein Feature Barcode inside GEM",StepKind.MAGIC,
            ((occupied_antibody,InputRole.TEMPLATE),(occupied_ab,InputRole.STRUCTURAL)),(protein_joined,),
            reason="CG000674 Rev C says the two oligos attach inside GEMs; their physical end geometry and bond formation have not yet been reconstructed from checked sources.",
            source="10x CG000674 Rev C Inside Individual GEMs, lines 1155–1165")
    protein_barcoded=g.product("protein_barcoded","cell-barcode/UMI-tagged antibody feature reporter",MaterialKind.MOLECULAR,
                               state=molecule("barcoded protein reporter",(role("cell barcode"),role("UMI"),role("antibody feature barcode"))))
    g.react("barcode_protein","Extend gel-bead primer on antibody capture oligo inside the GEM",StepKind.EXTENSION,
            ((protein_joined,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),
            (protein_barcoded,),source="10x CG000674 Rev C GEM Generation and CG000529 Rev C p. 4")
    reviewed=["cell","antibody_oligo","label_antibody","fixative","fix","left_probe","right_probe","hybridize",
              "ab_oligo","ligase","seal_probe_pair","droplets","partition","bead_oligo","extension_polymerase",
              "barcode","attach_ab","barcode_protein","post_hyb_wash","wash_unbound"]
    # Both reporter populations remain together through source Step3.1–3.3.
    recovery=g.input("gem_recovery","125µl Recovery Agent; invert5×,2min; remove pink oil phase (§3.1)",MaterialKind.CONDITION)
    mixed=g.product("recovered_reporters","aqueous RNA/protein reporter mixture",MaterialKind.MOLECULAR,
        state=reporter_mixture("recovered reporter mixture",amplified=False))
    oil=g.product("recovery_oil","Recovery Agent/Partitioning Oil phase",MaterialKind.WASTE)
    g.react("recover_gems","Recover the aqueous reporter populations together",StepKind.EXTRACTION,
        ((barcoded,InputRole.TEMPLATE),(protein_barcoded,InputRole.TEMPLATE),(recovery,InputRole.CONDITION)),
        (mixed,oil),changes_molecule=False)
    g.discard(oil,"pink Recovery Agent/Partitioning Oil aspirated from tube bottom (§3.1e)")
    pre_f=g.input("preamp_forward_role","Pre-Amp Primers C PN-2000953: first binding role of one supplied pool",MaterialKind.OLIGO,inferred=True)
    pre_r=g.input("preamp_reverse_role","Pre-Amp Primers C PN-2000953: opposite binding role of same pool",MaterialKind.OLIGO,inferred=True)
    pre_pol=g.input("preamp_polymerase","Amp Mix PN-2000047",MaterialKind.CATALYST)
    preamp=g.product("preamplified","eight-cycle preamplified reporter mixture",MaterialKind.MOLECULAR,
        state=reporter_mixture("preamplified reporter mixture",amplified=True))
    g.react("preamp_pcr","Pre-amplify recovered reporters8cycles (§3.2)",StepKind.PCR,
        ((mixed,InputRole.TEMPLATE),(pre_f,InputRole.STRUCTURAL),(pre_r,InputRole.STRUCTURAL),
         (pre_pol,InputRole.CATALYST)),(preamp,),geometry=PCRGeometry.unavailable(mixed.state,
            reason="One Pre-Amp PrimersC pool PN-2000953 is named. The graph records two binding roles of that pool, not two separately established oligo species; physical sequences/sites and product ends remain unassigned."))
    spri=g.input("preamp_cleanup","SPRIselect1.8× (126µl to70µl sample),100µl elution (§3.3)",MaterialKind.SELECTION)
    cleaned=g.product("cleaned_reporters","purified preamplified reporter mixture",MaterialKind.MOLECULAR,state=preamp.state)
    spri_waste=g.product("preamp_cleanup_waste","unbound short DNA/reagents in SPRI supernatant and ethanol washes",MaterialKind.WASTE)
    g.react("preamp_cleanup_step","SPRIselect cleanup of shared preamplification",StepKind.CLEANUP,
        ((preamp,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(cleaned,spri_waste),changes_molecule=False)
    g.discard(spri_waste,"supernatant and ethanol washes physically removed (§3.3)")
    rna_input=g.product("flex_index_template","20µl aliquot for TS gene-expression index PCR",MaterialKind.MOLECULAR,state=cleaned.state)
    protein_input=g.product("flex_protein_index_template","20µl aliquot for TN protein-expression index PCR",MaterialKind.MOLECULAR,state=cleaned.state)
    remainder=g.product("other_reporters","remaining cleaned reporter aliquots not followed",MaterialKind.MOLECULAR,state=cleaned.state)
    aliquot_rule=g.input("aliquot_rule","20µl preamp product into each §4.1/§5.1 sample-index PCR",MaterialKind.CONDITION)
    g.react("index_aliquots","Take aliquots for separate gene/protein library construction",StepKind.SAMPLING,
        ((cleaned,InputRole.TEMPLATE),(aliquot_rule,InputRole.CONDITION)),(rna_input,protein_input,remainder),changes_molecule=False)
    g.mark_terminal(Terminal(remainder.id,Disposition.UNRESOLVED,"The model follows20µl for each source library branch; subsequent handling of the remaining reporter mixture is unassigned."))
    # Gene-expression PCR selects the TS-compatible probe-derived reporters.
    ge_f=g.input("flex_forward","Dual Index TS SetA i5 primer (physical sequence unassigned)",MaterialKind.OLIGO)
    ge_r=g.input("flex_reverse","Dual Index TS SetA i7 primer (physical sequence unassigned)",MaterialKind.OLIGO)
    ge_pol=g.input("flex_polymerase","Amp Mix PN-2000103/2000047",MaterialKind.CATALYST)
    ge_product=g.product("flex_pcr_product","TS-indexed gene-expression PCR product",MaterialKind.MOLECULAR,state=LIB)
    g.react("flex_pcr","Sample Index PCR with Dual Index TS SetA (§4.1)",StepKind.PCR,
        ((rna_input,InputRole.TEMPLATE),(ge_f,InputRole.STRUCTURAL),(ge_r,InputRole.STRUCTURAL),
         (ge_pol,InputRole.CATALYST)),(ge_product,),geometry=PCRGeometry.unavailable(rna_input.state,
            reason="Physical TS primers/preamp ends remain unassigned; completed SmallRNA Read2S sites do not supply the preceding template geometry."))
    ge_spri=g.input("flex_spri","SPRIselect1.0×,40µl final eluate (§4.2)",MaterialKind.SELECTION)
    ge_library=g.product("flex_library","Flex probe-derived gene-expression sequencing library",MaterialKind.LIBRARY,state=LIB,inferred=True)
    ge_waste=g.product("flex_spri_waste","gene-expression index-PCR supernatant and washes",MaterialKind.WASTE)
    g.react("flex_cleanup","SPRIselect size-selection cleanup of gene-expression library",StepKind.CLEANUP,
        ((ge_product,InputRole.TEMPLATE),(ge_spri,InputRole.SELECTION)),(ge_library,ge_waste),changes_molecule=False)
    g.discard(ge_waste,"supernatant and ethanol washes removed in §4.2")
    g.sequence(ge_library,SequencingPlan(library=LIB,primers=PRIMERS,
        acquired_roles=tuple(p.role for p in PRIMERS),read_lengths=READ_LENGTHS))
    # Protein PCR has its own TN primer identities and distinct Nextera Read2 arm.
    pr_f=g.input("flex_protein_forward","Dual Index TN SetA i5 primer (physical sequence unassigned)",MaterialKind.OLIGO)
    pr_r=g.input("flex_protein_reverse","Dual Index TN SetA i7 primer (physical sequence unassigned)",MaterialKind.OLIGO)
    pr_pol=g.input("flex_protein_polymerase","Amp Mix PN-2000047/2000103",MaterialKind.CATALYST)
    pr_product=g.product("flex_protein_pcr_product","TN-indexed protein-expression PCR product",MaterialKind.MOLECULAR,state=PROTEIN_LIB)
    g.react("flex_protein_pcr","Sample Index PCR with Dual Index TN SetA (§5.1)",StepKind.PCR,
        ((protein_input,InputRole.TEMPLATE),(pr_f,InputRole.STRUCTURAL),(pr_r,InputRole.STRUCTURAL),
         (pr_pol,InputRole.CATALYST)),(pr_product,),geometry=PCRGeometry.unavailable(protein_input.state,
            reason="Physical TN primers and copied antibody/preamp ends remain unassigned. Final Nextera Read2 sites do not establish the preceding primer input geometry."))
    pr_spri=g.input("flex_protein_spri","SPRIselect1.0×,40µl final eluate (§5.2)",MaterialKind.SELECTION)
    pr_library=g.product("flex_protein_library","Flex antibody-derived protein-expression sequencing library",MaterialKind.LIBRARY,state=PROTEIN_LIB,inferred=True)
    pr_waste=g.product("flex_protein_spri_waste","protein index-PCR supernatant and washes",MaterialKind.WASTE)
    g.react("flex_protein_cleanup","SPRIselect size-selection cleanup of protein library",StepKind.CLEANUP,
        ((pr_product,InputRole.TEMPLATE),(pr_spri,InputRole.SELECTION)),(pr_library,pr_waste),changes_molecule=False)
    g.discard(pr_waste,"supernatant and ethanol washes removed in §5.2")
    g.sequence(pr_library,SequencingPlan(library=PROTEIN_LIB,primers=PROTEIN_PRIMERS,
        acquired_roles=tuple(p.role for p in PROTEIN_PRIMERS),read_lengths=READ_LENGTHS))
    reviewed.extend(("gem_recovery","recover_gems","preamp_forward_role","preamp_reverse_role","preamp_polymerase","preamp_pcr",
        "preamp_cleanup","preamp_cleanup_step","aliquot_rule","index_aliquots","flex_forward","flex_reverse","flex_polymerase","flex_pcr","flex_spri","flex_cleanup","flex_library",
        "flex_protein_forward","flex_protein_reverse","flex_protein_polymerase","flex_protein_pcr","flex_protein_spri","flex_protein_cleanup","flex_protein_library"))
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
