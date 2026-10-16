"""10x Chromium Fixed RNA Profiling / Flex probe-pair chemistry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg, truseq_library
from chemdraw import Construct, feature
from multimodal_spatial import adjacent_probe_scene, finish_pcr_branch
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE="10x Fixed RNA Profiling (Flex)"
NOTES="01_flex.html"
SOURCE='Commercial protocol: 10x Genomics <a href="https://www.10xgenomics.com/support/cn/flex-gene-expression/documentation/steps/library-prep/chromium-fixed-rna-profiling-reagent-kits-for-singleplexed-samples-with-feature-barcode-technology-for-protein-using-barcode-oligo-capture">Fixed RNA Profiling user guide CG000674 Rev C</a>.'
SUMMARY="Gene-specific probe pairs hybridize adjacently on fixed RNA. In GEMs, their nick is sealed and the reporter receives a cell barcode and UMI. Antibody-conjugated oligos on the same fixed cells are captured by the gel-bead primer and yield a separate protein-expression library."
CAVEAT="CG000674 Rev C Appendix p. 85 prints both completed library layouts; CG000529 Rev C p. 4 prints the antibody capture oligo pattern, and CG000621 Rev D Tables 7 and 13 print physical custom RNA-probe architecture. The compact RNA schematic still normalizes platform arms and omits the printed Small RNA Read-2 arm. 10x publishes commercial WTA target-region sequences through Cell Ranger; full physical LHS/RHS oligos and their barcode-mixture assignment have not been reconstructed, so the generic RNA probe region remains inferred. The protein shell is source-based, but its AB8, antibody barcode and sample index values vary by reagent or run."

def oligos():
    return [
        "<b>Flex gel-bead primer</b> — 5′-CTACACGACGCTCTTCCGATCT-[16-base cell barcode]-[12-base UMI]-TTGCTAGGACCG-3′",
        "<b>Published gene-expression reporter junction</b> — 5′-TTGCTAGGACCG-[BC8]-NN-TACGTGCTAACCGCGT-[ligated probe insert]-TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC-3′",
        "<b>Published custom Flex probe design</b> (<a href=\"https://cdn.10xgenomics.com/image/upload/v1729114202/support-documents/CG000621_CustomProbeDesign_RevD.pdf\">CG000621 Rev D Table 7</a>; not the commercial whole-transcriptome inventory) — LHS 5′-CCTTGGCACCCGAGAATTCCA-[25-nt target]-3′; RHS /5Phos/-[25-nt target]-ACGCGGTTAGCACGTA-NN-ACTTTAGG-CGGTCCTAGCAA-3′ for singleplex BC001",
        "<b>Official WTA reference example</b> — human v1.0.1 TSPAN6 spliced target-region probe_seq: GGTGACACCACAACAATGCAACGTATTTTGGATCTTGTCTACTGCATGGC (50 nt; reference sequence, not a complete physical oligo or a specimen-specific choice)",
        "<b>Protein antibody capture oligo</b> — <a href=\"https://cdn.10xgenomics.com/image/upload/v1719610805/support-documents/CG000529_DemonstratedProtocol_CellSurfaceProteinLabelingforFixedRNA_RevC.pdf\">CG000529 Rev C p. 4</a> prints /5AmMC12/CGGAGATGTGTATAAGAGACAG-N10-N15-N9-CCCATATAAGAAA for the TotalSeq-C Barcode Oligo Capture branch (N15 is the antibody feature barcode)",
    ]

LIB, PRIMERS = truseq_library([
    seg("cell barcode", "B"*16, "cbc", placeholder=True, inferred=True,
        feature=feature("cell_barcode", "cell_barcode", "whitelist", whitelist="10x Flex whitelist")),
    seg("UMI", "U"*12, "umi", placeholder=True, inferred=True,
        feature=feature("umi", "umi", "random")),
    seg("ligated probe-pair identifier", "X"*30, placeholder=True, inferred=True)],
    "Flex gene-expression library", inferred_adapters=True)
ANTIBODY_OLIGO = Construct([
    seg("Nextera partial Read 2", "CGGAGATGTGTATAAGAGACAG", "r2", note="5′ /5AmMC12/ antibody conjugation in CG000529 Rev C"),
    seg("antibody N10", "N"*10, placeholder=True),
    seg("antibody feature barcode", "N"*15, "cbc", placeholder=True, inferred=True,
        feature=feature("antibody_feature_barcode", "feature_barcode", "unknown")),
    seg("antibody N9", "N"*9, placeholder=True),
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
    seg("gel-bead capture terminus", "TTGCTAGGACCG"),
    seg("AB8 antibody multiplex barcode", "A"*8, "cbc", placeholder=True, inferred=True,
        feature=feature("protein_ab8", "sample_index", "unknown",
                        note="Antibody Multiplexing Barcode identifies the sample; CG000674 Rev C p. 20")),
    seg("protein N5", "N"*5, placeholder=True),
    seg("protein constant 1", "AACCGCGTATTCAGGT"),
    seg("protein constant 2", "TTTCTTATATGGG"),
    seg("protein N9", "N"*9, placeholder=True),
    seg("antibody feature barcode", "N"*15, "cbc", placeholder=True, inferred=True,
        feature=feature("antibody_feature_barcode", "feature_barcode", "unknown")),
    seg("protein N10", "N"*10, placeholder=True),
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
                  "INFERRED — This compact final duplex preserves read roles but normalizes platform arms; the source-exact Flex reporter constants are transcribed in Key oligos and research notes.",
                  "The declared Illumina primer sites are shown on the completed probe-derived library."),
                 ("Final antibody-derived protein-expression library", PROTEIN_LIB, PROTEIN_PRIMERS,
                  "INFERRED — The guide prints this shell and the antibody oligo pattern, but the AB8 and N15 barcode values and selected 10-base indexes are reagent- or run-dependent.",
                  "TruSeq Read 1 and Nextera Read 2 primer sites are computed on the separate protein library."))

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cell=g.input("cell","unfixed cell carrying target RNA and surface proteins",MaterialKind.BIOLOGICAL,
                 state=molecule("cell before Flex labeling",(role("target RNA"),role("surface protein"))),primary=True)
    antibody=g.input("antibody_oligo","TotalSeq-C antibody-conjugated Barcode Oligo Capture oligo",MaterialKind.SELECTION,
                     state=ANTIBODY_OLIGO, inferred=True,
                     source="10x CG000529 Rev C p. 4; N15 antibody identity depends on the conjugate")
    labeled=g.product("labeled","cell labeled with antibody-oligonucleotide conjugate",MaterialKind.MOLECULAR,
                      state=molecule("antibody-labeled cell",(role("target RNA"),role("bound antibody oligo"))))
    unbound_antibody=g.product("unbound_antibody","unbound antibody-oligonucleotide conjugate",MaterialKind.WASTE)
    g.react("label_antibody","Label surface protein with antibody-oligonucleotide conjugate",StepKind.BINDING,
            ((cell,InputRole.TEMPLATE),(antibody,InputRole.SELECTION)),(labeled,unbound_antibody),
            source="10x CG000529 Rev C pp. 4–5; CG000478 Rev C fixation protocol")
    g.discard(unbound_antibody,"unbound antibody conjugates are removed before fixed-cell processing")
    fixative=g.input("fixative","Flex cell fixation condition",MaterialKind.CONDITION,
                     source="10x CG000478 Rev C fixation protocol")
    fixed=g.product("fixed_labeled_cell","fixed antibody-labeled cell with target RNA",MaterialKind.MOLECULAR,
                    state=molecule("fixed labeled cell",(role("target RNA"),role("bound antibody oligo"))))
    g.react("fix","Fix antibody-labeled cells",StepKind.CROSS_LINKING,
            ((labeled,InputRole.TEMPLATE),(fixative,InputRole.CONDITION)),(fixed,),
            source="10x CG000478 Rev C; optional protein labeling precedes fixation")
    left=g.input("left_probe","gene-specific left probe half",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (panel-specific physical oligo not reconstructed)")
    right=g.input("right_probe","gene-specific right probe half",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (panel-specific physical oligo not reconstructed)")
    ab_oligo=g.input("ab_oligo","selected Antibody Multiplexing Barcode oligo AB001–AB016",MaterialKind.OLIGO,
                     inferred=True,source="10x CG000674 Rev C Probe Hybridization step 1.1h; individual physical oligo not reconstructed")
    adjacent=g.product("adjacent","adjacent probe halves hybridized on target RNA",MaterialKind.MOLECULAR,state=adjacent_probe_scene())
    retained_antibody=g.product("retained_antibody","antibody-conjugated capture oligo on the hybridized cell",MaterialKind.MOLECULAR,
                                state=ANTIBODY_OLIGO)
    retained_ab=g.product("retained_ab","selected Antibody Multiplexing Barcode oligo in probe-hybridized suspension",MaterialKind.MOLECULAR,
                          inferred=True,source="10x CG000674 Rev C step 1.1h; AB001–AB016 sequence choice is sample-specific")
    g.react("hybridize","Hybridize adjacent probe halves to fixed RNA",StepKind.HYBRIDIZATION,
            ((fixed,InputRole.TEMPLATE),(left,InputRole.STRUCTURAL),(right,InputRole.STRUCTURAL),(ab_oligo,InputRole.STRUCTURAL)),
            (adjacent,retained_antibody,retained_ab))
    droplets=g.input("droplets","Flex GEM partitioning condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","partitioned RNA-bound adjacent probe halves",MaterialKind.MOLECULAR,state=adjacent.state)
    occupied_antibody=g.product("occupied_antibody","partitioned antibody-conjugated oligo",MaterialKind.MOLECULAR,
                                state=ANTIBODY_OLIGO)
    occupied_ab=g.product("occupied_ab","partitioned Antibody Multiplexing Barcode oligo",MaterialKind.MOLECULAR,
                          inferred=True,source="10x CG000674 Rev C GEM Generation; physical AB oligo not reconstructed")
    failed=g.product("failed","empty or unresolved GEM material",MaterialKind.WASTE)
    g.react("partition","Partition hybridized RNA probes and antibody conjugates with barcoded gel beads",StepKind.PARTITION,
            ((adjacent,InputRole.TEMPLATE),(retained_antibody,InputRole.STRUCTURAL),(retained_ab,InputRole.STRUCTURAL),
             (droplets,InputRole.CONDITION)),(occupied,occupied_antibody,occupied_ab,failed),changes_molecule=False)
    g.discard(failed,"does not retain a resolvable cell-associated reporter")
    ligase=g.input("ligase","probe-pair ligase in GEM",MaterialKind.CATALYST)
    reporter=g.product("reporter","ligated probe-pair reporter in GEM",MaterialKind.MOLECULAR,state=molecule("ligated Flex reporter",(role("left probe barcode"),role("right probe barcode")),junctions=(("left probe barcode","right probe barcode","probe-pair ligation"),)))
    g.react("seal_probe_pair","Ligate the adjacent probe pair inside the GEM",StepKind.NICK_SEALING,((occupied,InputRole.TEMPLATE),(ligase,InputRole.CATALYST)),(reporter,))
    bead=g.input("bead_oligo","cell-barcode/UMI gel-bead oligo CTACACGACGCTCTTCCGATCT-N16-N12-TTGCTAGGACCG",MaterialKind.OLIGO,source=SOURCE+" Appendix p. 85")
    pol=g.input("extension_polymerase","barcode-transfer polymerase",MaterialKind.CATALYST)
    barcoded=g.product("barcoded","cell-barcode/UMI-tagged probe reporters",MaterialKind.MOLECULAR,state=molecule("barcoded Flex reporter",(role("cell barcode"),role("UMI"),role("probe-pair identifier"))))
    g.react("barcode","Transfer cell barcode and UMI to the ligated reporter",StepKind.EXTENSION,((reporter,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(barcoded,))
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
              "barcode","attach_ab","barcode_protein"]
    finish_pcr_branch(g,"flex",barcoded,LIB,PRIMERS,reviewed,label="Flex probe-derived gene-expression",inferred_primers=True)
    finish_pcr_branch(g,"flex_protein",protein_barcoded,PROTEIN_LIB,PROTEIN_PRIMERS,reviewed,
                      label="Flex antibody-derived protein-expression",inferred_primers=True,
                      limitation="CG000674 Rev C prints the protein-library shell but does not define a single universal protein-index PCR primer pair for every index well.")
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
