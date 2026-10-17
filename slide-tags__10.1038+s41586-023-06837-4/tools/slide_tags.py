"""Slide-tags molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from batch_ngs import seg,truseq_library
from chemdraw import Construct,Row,feature,revcomp
import illumina as il
import nextera as nx
import seqprimers as sp
from multimodal_spatial import finish_pcr_branch, modality_split
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from pcr import PCRGeometry
from protocol_model import molecule, role

TITLE="Slide-tags — spatial barcodes transferred to intact nuclei"
NOTES="01_slide-tags.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41586-023-06837-4">Russell et al., <i>Nature</i> (2024)</a>.'
SUMMARY="Photocleaved oligos from a sequenced bead array diffuse into overlying tissue nuclei; 10x single-nucleus RNA or Multiome processing then links spatial tags to each nucleus barcode."
CAVEAT="The published TAGS/SLAC transfer logic and separate 3′ v3.1 and Multiome branches are shown. Tn5 tagmentation belongs only to Multiome. The cited 10x CG000317 Rev C and CG000338 Rev F guides print fixed gel-bead scaffolds and library structures; variable cell barcodes, UMIs and puck barcodes remain placeholders here. The <a href=\"https://github.com/broadchenf/Slide-tags/blob/main/cell_barcode_matcher.R\">author parser</a> reads a 16-base cell barcode and 12-base feature UMI from its cell-barcode FASTQ input, and reconstructs a 14-base puck barcode from positions 1–8 and 27–32 of its spatial FASTQ input. The script reverses internal R1/R2 variable names relative to file paths; branch-specific run lengths are not asserted here. Derived from the printed TAGS oligo, the published CG000317 3′ feature-library shell and the Supplementary Methods statement that TAGS spatial libraries were prepared as Cell Surface Protein libraries: Read 2 reads the puck J8 block (cycles 1–8), the printed 18-nt linker (9–26) and the J7 block (27–33). The parser uses cycles 1–8 and 27–32, so of the 15 printed <code>J</code> positions it omits the seventh J7 base. The drawn 3′ tag library applies only to TAGS arrays; the SLAC (ligated-puck) mouse-brain branch, amplified with the named STAG_P701_NEX and P5-Truseq Hybrid primers, is not drawn. The Multiome spatial-tag library is a sequencing endpoint in the graph with a visible warning; the puck design used with Multiome is not stated in the inspected Methods, and no strand-level shell is assigned to it. <b>Published discrepancy:</b> the main article Methods (PDF p. 10; PMC10764288 Methods) print <code>TTT_PC_zCCGG...</code> for the SLAC bead oligo, while the published Supplementary Methods section 1.3 and preceding bioRxiv Methods print <code>TTT_PC_GCCGG...</code>. The sequence below uses the supplement's base-only form because <code>z</code> is not a defined nucleotide; its intended meaning in the main article remains uncertain. The <a href=\"https://doi.org/10.1038/s41587-020-0739-1\">cited Slide-seqV2 Methods</a> define <code>J</code> as a split-pool bead barcode base and also print a no-<code>z</code> SLAC precursor."
CELL=feature("slidetags_cell","cell_barcode","whitelist",whitelist="10x Chromium bead barcode")
SPATIAL=feature("slidetags_spatial","spatial_barcode","whitelist",whitelist="sequenced Slide-tags puck")
UMI=feature("slidetags_umi","umi","random")
I5=feature("slidetags_i5","sample_index","unknown")
I7=feature("slidetags_i7","sample_index","unknown")
RNA,RNA_P=truseq_library([seg("cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),seg("RNA UMI","U"*12,"umi",placeholder=True,feature=UMI),seg("nuclear cDNA","X"*40,placeholder=True)],"Slide-tags 3-prime RNA library")
MULTI_RNA,MULTI_RNA_P=truseq_library([seg("cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),seg("RNA UMI","U"*12,"umi",placeholder=True,feature=UMI),seg("nuclear cDNA","X"*40,placeholder=True)],"Slide-tags Multiome RNA library")
# CG000338 Rev F (the Multiome revision the Supplementary Methods 1.5 name),
# Appendix p. 81: ATAC gel-bead oligo and final ATAC library. The 16-base GEM
# barcode and 8-base spacer sit between P5 and the Read 1N (S5 + ME) arm and are
# read by the 24-cycle i5 read (Rev F Sequencing); the only sample index is i7 (N8).
MULTIOME_SPACER="CGCGTCTG"
ATAC_I7=feature("slidetags_atac_i7","sample_index","unknown",
                note="CG000338 Rev F Sample Index N (i7:8); the index set is not printed")
def multiome_atac_library():
    lib=Construct([
        seg("P5",il.P5,"p5"),
        seg("10x cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),
        seg("10x spacer",MULTIOME_SPACER),
        seg("Read 1N S5 arm",nx.S5,"s5"),
        seg("left mosaic end",nx.ME,"me"),
        seg("accessible chromatin","X"*40,placeholder=True),
        seg("right mosaic end reverse complement",nx.ME_RC,"me"),
        seg("Read 2N S7 arm reverse complement",nx.S7_RC,"s7"),
        seg("i7 sample index","I"*8,"cbc",placeholder=True,feature=ATAC_I7),
        seg("P7 reverse complement",il.P7_RC,"p7"),
    ],name="Slide-tags Multiome ATAC library")
    primers=tuple(sp.NEXTERA[k] for k in ("R1","I1","I2","R2"))
    if problems:=sp.verify(lib,primers,required_roles=tuple(p.role for p in primers)):
        raise ValueError("Slide-tags Multiome ATAC library: "+"; ".join(problems))
    return lib,primers
ATAC,ATAC_P=multiome_atac_library()
def bead_primer(name,handle,three_prime):
    """10x gel-bead primer printed as <handle>-N16-N12-<3' end> in the cited guides."""
    return Construct([seg("bead handle",handle,"r1"),
                      seg("cell barcode","N"*16,"cbc",placeholder=True,feature=CELL),
                      seg("UMI","N"*12,"umi",placeholder=True,feature=UMI),
                      *three_prime],name=name)
def phosphorothioate_primer(name,bases,ps_bonds):
    """Custom PCR primer; the last ps_bonds linkages are printed phosphorothioates."""
    return Construct([seg(name,bases,note=f"{ps_bonds} 3'-terminal phosphorothioate linkages (printed *)")],name=name)
STAG_IP7_A1="CAAGCAGAAGACGGCATACGAGATATTTACCGCAGTGACTGGAGTTCAGACGTGT"   # MOESM1 line 89
P5_STAG_IP5_A1="AATGATACGGCGACCACCGAGATCTACACGACAATAAAGACACTCTTTCCCTACACGACGCTC"  # MOESM1 line 89
# Published TAGS bead oligo (Supplementary Methods 1.3, MOESM1 line 73):
# 5'-TTT-PC-<TruSeq Read 2 arm>-J8-<linker>-J7-N7-VV-<capture complement>-3'.
TAGS_R2_ARM="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
TAGS_LINKER="TCTTCAGCGTTCCCGAGA"
TAGS_CAPTURE_COMPLEMENT="GCTTTAAGGCCGGTCCTAGCAA"
TAGS_OLIGO=TAGS_R2_ARM+"J"*8+TAGS_LINKER+"J"*7+"N"*7+"VV"+TAGS_CAPTURE_COMPLEMENT
SLAC_SPLINT_PRINTED_3_TO_5="CTGTGTACCCGCCTCGAAACGATTGC"  # MOESM1 line 70
CAPTURE_SEQUENCE_1="TTGCTAGGACCGGCCTTAAAGC"  # CG000317 Rev C Appendix p. 71
# CG000317 Rev C p. 71 and CG000338 Rev F p. 81 print the same RNA scaffold.
GEM_RNA_PRIMER=bead_primer("10x 3′ v3.1 / Multiome RNA gel-bead primer","CTACACGACGCTCTTCCGATCT",
                           (seg("oligo(dT)30","T"*30),seg("anchor","VN",placeholder=True)))
GEM_FEATURE_PRIMER=bead_primer("10x 3′ v3.1 feature-barcode gel-bead primer","GTCAGATGTGTATAAGAGACAG",
                               (seg("Capture Sequence 1",CAPTURE_SEQUENCE_1,"tso"),))
if revcomp(TAGS_CAPTURE_COMPLEMENT)!=CAPTURE_SEQUENCE_1:
    raise ValueError("TAGS 3' end must be complementary to 10x Capture Sequence 1")
if "A"+il.INDEX1_PRIMER!=revcomp(TAGS_R2_ARM):
    raise ValueError("TAGS 5' arm must yield the TruSeq Read 2/Index 1 arm")
def feature_tag_library():
    # CG000317 Rev C Appendix pp. 71-72 establishes the 3' v3.1 feature-library
    # shell; Supplementary Methods 1.5 says TAGS-array spatial libraries were
    # prepared as Cell Surface Protein libraries. The payload between Capture
    # Sequence 1 and the Read 2 arm is the reverse complement of the published
    # TAGS oligo (derived), so Read 2 reads J8, the 18-nt linker, then J7.
    lib=Construct([
        seg("P5",il.P5,"p5"),
        seg("sample index i5","J"*10,"cbc",placeholder=True,feature=I5),
        seg("Nextera Read 1N",nx.READ1_PRIMER,"r1"),
        seg("cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),
        seg("feature UMI","U"*12,"umi",placeholder=True,feature=UMI),
        seg("Capture Sequence 1",CAPTURE_SEQUENCE_1,"tso"),
        seg("puck VV complement","BB",placeholder=True),
        seg("puck mixed-base N7 complement","N"*7,placeholder=True),
        seg("puck J7 base 7 (Read 2 cycle 33; unused by author parser)","S",
            "cbc",placeholder=True,feature=SPATIAL),
        seg("puck J7 bases 1-6 (Read 2 cycles 27-32)","S"*6,"cbc",placeholder=True,feature=SPATIAL),
        seg("puck linker complement",revcomp(TAGS_LINKER)),
        seg("puck J8 (Read 2 cycles 1-8)","S"*8,"cbc",placeholder=True,feature=SPATIAL),
        seg("Read 2 leading A","A","r2"),
        seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),
        seg("sample index i7 reverse complement","I"*10,"cbc",placeholder=True,feature=I7),
        seg("P7 reverse complement",il.P7_RC,"p7"),
    ],name="Slide-tags 3′ TAGS spatial-tag library")
    primers=(sp.NEXTERA["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
    return lib,primers
TAG,TAG_P=feature_tag_library()
FINAL_LIBRARIES=(
 ("RNA library (3′ branch)",RNA,RNA_P,"10x 3′ v3.1 gene-expression library from tagged nuclei.","The droplet cell barcode and UMI identify each transcript."),
 ("RNA library (Multiome branch)",MULTI_RNA,MULTI_RNA_P,"10x Multiome gene-expression library from tagged nuclei.","The droplet cell barcode and UMI identify each transcript."),
 ("ATAC library (Multiome branch)",ATAC,ATAC_P,"CG000338 Rev F Appendix p. 81 ATAC library from the same spatially tagged nuclei: P5, 16-base GEM barcode and 8-base spacer precede the Read 1N arm; the i7 read gives the 8-base sample index.","Read 1N and Read 2N read the accessible fragment; the 24-cycle i5 read covers the spacer and the GEM barcode (Rev F Sequencing)."),
 ("Spatial-tag library (3′ branch, TAGS arrays)",TAG,TAG_P,"CG000317 Cell Surface Protein library shell with the reverse complement of the published TAGS puck oligo between Capture Sequence 1 and the Read 2 arm (derived).","The 10x guide specifies a 12-base feature UMI. Author code reconstructs a 14-base puck barcode from Read 2 cycles 1–8 (J8) and 27–32 (first six J7 bases)."),
)
READ_LENGTHS={"RNA library (3′ branch)":{},"RNA library (Multiome branch)":{},"ATAC library (Multiome branch)":{},"Spatial-tag library":{}}

def oligos():
 return [
  "<b>SLAC ligated-design bead oligo (Supplementary Methods 1.3)</b> — 5′-TTT_PC_GCCGGTAATACGACTCACTATAGGGCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJNNNNNNNVVGCTCGGACACATGGGCG-3′ (main article instead prints <code>TTT_PC_zCCGG...</code>; see discrepancy above)",
  "<b>Direct-synthesis TAGS puck oligo</b> — 5′-TTT-PC-"+TAGS_OLIGO+"-3′ (PC: photocleavable linker; J: split-pool bead barcode base, defined in cited Slide-seqV2 Methods)",
  "<b>SLAC 10x FB1 extension</b> — 5′-GAGCTTTGCTAACGGTCGAGGCTTTAAGGCCGGTCCTAGCAA-3′",
  "<b>SLAC splint</b> — 5′-"+SLAC_SPLINT_PRINTED_3_TO_5[::-1]+"-3′ (source prints it 3′→5′ as "+SLAC_SPLINT_PRINTED_3_TO_5+")",
  "<b>Poly(A) puck oligo</b> — 5′-TTT-PC-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJNNNNNNNVVA30-3′",
  "<b>10x 3′ v3.1 RNA gel-bead primer (CG000317 Rev C, p. 71)</b> — 5′-CTACACGACGCTCTTCCGATCT-N16-N12-T30VN-3′",
  "<b>10x 3′ v3.1 feature-barcode gel-bead primer (CG000317 Rev C, p. 71)</b> — 5′-GTCAGATGTGTATAAGAGACAG-N16-N12-TTGCTAGGACCGGCCTTAAAGC-3′",
  "<b>10x Multiome ATAC gel-bead primer (CG000338 Rev F, p. 81)</b> — 5′-AATGATACGGCGACCACCGAGATCTACAC-N16-CGCGTCTG-3′",
  "<b>Multiome spatial-tag spike-in</b> — 5′-GTGACTGGAGTTCAGACGT-3′",
  "<b>Multiome spatial-tag P7 primer</b> — 5′-CAAGCAGAAGACGGCATACGAGATATTTACCGCAGTGACTGGAGTTCAGACGT*G*T-3′ (*: phosphorothioate linkage)",
  "<b>Multiome spatial-tag P5 primer</b> — 5′-AATGATACGGCGACCACCGAGATCTACACGACAATAAAGACACTCTTTCCCTACACGACGC*T*C-3′ (*: phosphorothioate linkage)",
 ]

def sections(): return [
 ("Place tissue on a spatially sequenced bead array",[Row(chunks=[("bead — photocleavable linker — spatial barcode — 10x-capture sequence","cbc",False)])],"Each puck barcode already has a measured x/y coordinate."),
 ("Photocleave spatial oligos into the overlying tissue",[Row(chunks=[("UV cleavage → local spatial tags enter nearby nuclei","w1",False)])],"Short diffusion transfers several local bead identities into each intact nucleus."),
 ("Dissociate tagged nuclei and partition with 10x beads",[Row(chunks=[("spatial tag + one GEM cell barcode", "cbc",False)])],"Feature-barcode capture associates puck tags with a nucleus barcode."),
 ("Prepare assay and spatial-tag libraries",modality_split("RNA library","ATAC library in Multiome experiments","spatial-tag library"),"For TAGS arrays the puck payload is the reverse complement of the printed TAGS oligo inside the CG000317 feature-library shell (derived); the Multiome spatial-tag junction remains unresolved."),
]

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 array=g.input("array","spatially sequenced photocleavable bead array",MaterialKind.MOLECULAR,state=molecule("Slide-tags array",(role("mapped spatial barcode"),role("photocleavable linker"),role("10x capture sequence"))),primary=True)
 uv=g.input("uv","photocleavage condition",MaterialKind.CONDITION)
 tags=g.product("released_tags","photocleaved spatial barcode oligos",MaterialKind.MOLECULAR,state=molecule("released spatial tags",(role("spatial barcode"),role("10x capture sequence"))))
 g.react("photocleave","Photocleave spatial oligos from the mapped bead array",StepKind.CHEMICAL_CONVERSION,((array,InputRole.TEMPLATE),(uv,InputRole.CONDITION)),(tags,))
 nuclei=g.input("nuclei","intact tissue nuclei containing RNA and accessible chromatin",MaterialKind.BIOLOGICAL,state=molecule("tissue nuclei",(role("poly(A) RNA"),role("accessible chromatin"))),primary=True)
 tagged=g.product("tagged_nuclei","nuclei retaining local spatial tags",MaterialKind.MOLECULAR,state=molecule("spatially tagged nuclei",(role("spatial barcode tags"),role("poly(A) RNA"),role("accessible chromatin"))))
 unbound=g.product("unbound_tags","unbound or diffusely transferred spatial tags",MaterialKind.WASTE)
 g.react("bind_tags","Transfer local spatial tags into nearby nuclei",StepKind.BINDING,((nuclei,InputRole.TEMPLATE),(tags,InputRole.ANALYTE)),(tagged,unbound)); g.discard(unbound,"is not associated with a recovered nucleus")
 # The paper uses the tagged nuclei in two separate 10x workflows. The 3′
 # assay does not pass through the Multiome ATAC tagmentation step.
 three_prime=g.product("three_prime_nuclei","tagged nuclei assigned to the 3′ assay",MaterialKind.MOLECULAR,state=tagged.state)
 multiome=g.product("multiome_nuclei","tagged nuclei assigned to the Multiome assay",MaterialKind.MOLECULAR,state=tagged.state)
 g.react("assay_split","Route tagged nuclei to the separate 3′ and Multiome assays",StepKind.SPLIT,((tagged,InputRole.TEMPLATE),),(three_prime,multiome),changes_molecule=False)
 transposome=g.input("transposome","Multiome loaded transposome",MaterialKind.ADAPTER,source=SOURCE+"; CG000338 Rev F Appendix p. 81 gives final ATAC library scaffold")
 tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagmented=g.product("tagmented","Multiome nuclei with tagmented accessible chromatin",MaterialKind.MOLECULAR,state=multiome.state)
 twaste=g.product("tn5_waste","non-library tagmentation products",MaterialKind.WASTE)
 g.react("tagment","Tagment accessible chromatin only in the Multiome assay",StepKind.TAGMENTATION,((multiome,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagmented,twaste)); g.discard(twaste,"does not retain an amplifiable accessibility fragment")
 droplets=g.input("droplets","10x GEM partitioning condition",MaterialKind.CONDITION)
 occupied_3p=g.product("occupied_3p","3′ assay cell-barcoded products",MaterialKind.MOLECULAR,state=three_prime.state)
 failed_3p=g.product("failed_3p","empty or unresolved 3′ GEM products",MaterialKind.WASTE)
 g.react("partition_3p","Partition 3′ assay nuclei with barcoded gel beads",StepKind.PARTITION,((three_prime,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied_3p,failed_3p),changes_molecule=False); g.discard(failed_3p,"does not retain a resolvable nucleus barcode")
 occupied_multi=g.product("occupied_multi","Multiome cell-barcoded products",MaterialKind.MOLECULAR,state=tagmented.state)
 failed_multi=g.product("failed_multi","empty or unresolved Multiome GEM products",MaterialKind.WASTE)
 g.react("partition_multi","Partition Multiome nuclei with barcoded gel beads",StepKind.PARTITION,((tagmented,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied_multi,failed_multi),changes_molecule=False); g.discard(failed_multi,"does not retain a resolvable nucleus barcode")
 rna_3p=g.product("rna_3p","3′ RNA branch",MaterialKind.MOLECULAR,state=occupied_3p.state)
 tag_3p=g.product("tag_3p","3′ spatial-tag branch",MaterialKind.MOLECULAR,state=occupied_3p.state)
 g.react("split_3p","Separate 3′ RNA and spatial-tag products",StepKind.SPLIT,((occupied_3p,InputRole.TEMPLATE),),(rna_3p,tag_3p),changes_molecule=False)
 rna_multi=g.product("rna_multi","Multiome RNA branch",MaterialKind.MOLECULAR,state=occupied_multi.state)
 atac=g.product("atac_branch","Multiome ATAC branch",MaterialKind.MOLECULAR,state=occupied_multi.state)
 tag_multi=g.product("tag_multi","Multiome spatial-tag branch",MaterialKind.MOLECULAR,state=occupied_multi.state)
 g.react("split_multi","Separate Multiome RNA, ATAC and spatial-tag products",StepKind.SPLIT,((occupied_multi,InputRole.TEMPLATE),),(rna_multi,atac,tag_multi),changes_molecule=False)
 bead_3p=g.input("bead_primer_3p","3′ v3.1 cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO,state=GEM_RNA_PRIMER,source=SOURCE+"; CG000317 Rev C Appendix p. 71 prints fixed scaffold")
 bead_multi=g.input("bead_primer_multi","Multiome cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO,state=GEM_RNA_PRIMER,source=SOURCE+"; CG000338 Rev F Appendix p. 81 prints fixed scaffold")
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 cdna_3p=g.product("cdna_3p","3′ cell-barcoded nuclear cDNA",MaterialKind.MOLECULAR,state=molecule("Slide-tags 3′ cDNA",(role("cell barcode"),role("RNA UMI"),role("nuclear cDNA"))))
 g.react("reverse_transcribe_3p","Capture and reverse-transcribe 3′ assay nuclear RNA",StepKind.REVERSE_TRANSCRIPTION,((rna_3p,InputRole.TEMPLATE),(bead_3p,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna_3p,))
 cdna_multi=g.product("cdna_multi","Multiome cell-barcoded nuclear cDNA",MaterialKind.MOLECULAR,state=molecule("Slide-tags Multiome cDNA",(role("cell barcode"),role("RNA UMI"),role("nuclear cDNA"))))
 g.react("reverse_transcribe_multi","Capture and reverse-transcribe Multiome nuclear RNA",StepKind.REVERSE_TRANSCRIPTION,((rna_multi,InputRole.TEMPLATE),(bead_multi,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna_multi,))
 capture=g.input("tag_capture_3p","3′ spatial-tag capture primer",MaterialKind.OLIGO,state=GEM_FEATURE_PRIMER,source=SOURCE+"; CG000317 Rev C Appendix p. 71 prints feature-barcode capture scaffold")
 pol=g.input("tag_extension","spatial-tag barcode-transfer polymerase",MaterialKind.CATALYST)
 copied=g.product("copied_tags_3p","3′ cell-barcoded spatial-tag products",MaterialKind.MOLECULAR,state=molecule("Slide-tags 3′ tag precursor",(role("cell barcode"),role("tag UMI"),role("spatial barcode"))))
 g.react("copy_tags_3p","Copy 3′ spatial tags with the nucleus barcode and UMI",StepKind.EXTENSION,((tag_3p,InputRole.TEMPLATE),(capture,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(copied,))
 copied_multi=g.product("copied_tags_multi","Multiome spatial-tag products",MaterialKind.MOLECULAR,state=molecule("Slide-tags Multiome tag precursor",(role("cell barcode"),role("spatial barcode"))))
 g.react("copy_tags_multi","INFERRED — associate Multiome spatial tags with the nucleus barcode",StepKind.MAGIC,((tag_multi,InputRole.TEMPLATE),),(copied_multi,),reason="The published Multiome spike-in and custom PCR are named, but the exact spatial-tag capture junction is not reconstructed from their oligo sequences.")
 reviewed=["array","uv","photocleave","nuclei","bind_tags","assay_split","transposome","tn5","tagment","droplets","partition_3p","partition_multi","split_3p","split_multi","bead_primer_3p","bead_primer_multi","rt","reverse_transcribe_3p","reverse_transcribe_multi","tag_capture_3p","tag_extension","copy_tags_3p","copy_tags_multi"]
 finish_pcr_branch(g,"rna_3p",cdna_3p,RNA,RNA_P,reviewed,label="Slide-tags 3′ RNA",read_lengths=READ_LENGTHS["RNA library (3′ branch)"])
 finish_pcr_branch(g,"rna_multi",cdna_multi,MULTI_RNA,MULTI_RNA_P,reviewed,label="Slide-tags Multiome RNA",read_lengths=READ_LENGTHS["RNA library (Multiome branch)"])
 finish_pcr_branch(g,"atac",atac,ATAC,ATAC_P,reviewed,label="Slide-tags Multiome ATAC",read_lengths=READ_LENGTHS["ATAC library (Multiome branch)"])
 finish_pcr_branch(g,"tag_3p",copied,TAG,TAG_P,reviewed,label="Slide-tags 3′ spatial tag",read_lengths=READ_LENGTHS["Spatial-tag library"])
 # The Multiome branch is sequenced with the paper's custom PCR oligos, but the
 # complete strand-level shell is not established here, so its plan is explicit
 # and limited rather than borrowed from the 3′ Feature Barcode example.
 m_f=g.input("tag_multi_forward","Multiome spatial-tag P5-STAG_ip5_a1 primer",MaterialKind.OLIGO,state=phosphorothioate_primer("P5-STAG_ip5_a1",P5_STAG_IP5_A1,2),source=SOURCE+" Methods: Multiome custom PCR")
 m_r=g.input("tag_multi_reverse","Multiome spatial-tag STAG_iP7_a1 primer",MaterialKind.OLIGO,state=phosphorothioate_primer("STAG_iP7_a1",STAG_IP7_A1,2),source=SOURCE+" Methods: Multiome custom PCR")
 m_pol=g.input("tag_multi_polymerase","Multiome spatial-tag PCR polymerase",MaterialKind.CATALYST)
 m_lib=g.product("tag_multi_library","Multiome spatial-tag sequencing library (junction unresolved)",MaterialKind.LIBRARY,state=molecule("Multiome spatial-tag library",(role("cell barcode"),role("spatial tag"),role("custom P5/P7 PCR arms"))),inferred=True)
 g.react("tag_multi_pcr","Amplify and index the Multiome spatial-tag branch",StepKind.PCR,((copied_multi,InputRole.TEMPLATE),(m_f,InputRole.STRUCTURAL),(m_r,InputRole.STRUCTURAL),(m_pol,InputRole.CATALYST)),(m_lib,),geometry=PCRGeometry.unavailable(Construct([seg("Multiome spatial-tag precursor","X",placeholder=True,inferred=True)],"Multiome spatial-tag PCR substrate (junction unresolved)"),reason="The printed P5/P7 primers do not yet establish the branch-specific precursor junction."))
 g.sequence(m_lib,SequencingPlan(limitation="The published Multiome custom PCR primers are printed, but the complete strand-level library and sequencing-primer sites are not reconstructed."))
 reviewed.extend(("tag_multi_forward","tag_multi_reverse","tag_multi_polymerase","tag_multi_pcr","tag_multi_library"))
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
