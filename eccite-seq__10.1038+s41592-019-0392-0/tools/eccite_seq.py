"""ECCITE-seq 5-prime multimodal tag and direct-guide capture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg,truseq_library
from chemdraw import Construct, feature, revcomp
from pcr import PCRGeometry
from multimodal_spatial import finish_pcr_branch
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE="ECCITE-seq"
NOTES="01_eccite-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41592-019-0392-0">Mimitou et al., <i>Nature Methods</i> (2019)</a>.'
SUMMARY="ECCITE-seq adapts 10x 5′ chemistry so transcript, antibody tag, sample hashtag and directly captured sgRNA products receive a common cell barcode and UMI. Modality-specific handles then support separate libraries, alongside standard immune-receptor enrichment when desired."
CAVEAT=("The schematic follows the Methods: all molecules share one 10x 5′ RT in the GEM; "
        "cDNA amplification carries spiked hashtag, guide-tag and protein-tag additive primers; "
        "0.6X SPRI keeps mRNA cDNA on beads and the short tag/guide products in the supernatant, "
        "which is cleaned twice at 2.0X and split into separate SI-PCR + RPI-x, D7xx_s and "
        "Next_nst_x PCRs (MOESM3 first enriches guides with SI-PCR + gd_add_v4). Tag and guide libraries are assembled from the printed SI-PCR, tag, "
        "RPI1, D701_s and Next_nst_x bases; the 16-nt cell barcode, 10-nt UMI and the last nine "
        "Read 1 bases come from the 10x bead oligo, which ECCITE sources do not print (marked inferred), "
        "and the transcript library's 10x adapters are inferred. The protein-tag library's Index 1 "
        "and Read 2 primers are derived from the small-RNA RPI1 handle; the sources do not name "
        "the sequencing primer mix. "
        "Supplementary Table 3 in the publisher's MOESM1 prints 49 antibody barcodes "
        "of 15 bases and three of 12 bases (lines 219–387); Tables 4–5 print 12-base "
        "mouse/human hashtags (lines 389–460). Representative examples of all three "
        "sizes are drawn here. "
        "The MOESM3 oligo examples (lines 385–396) instead show eight <code>x</code> "
        "characters despite calling the example barcode 12 nt; the identity of those "
        "eight positions is unresolved. They are not substituted for the Table 3–5 barcodes. "
        "Methods name the scaffold RT primer gd_RT_v4, while the MOESM3 oligo list (line 399) "
        "calls the same sequence sg_RT_v4. Immune-receptor (TCR α/β and γ/δ) enrichment follows "
        "the 10x V(D)J branch and is not drawn.")
CELL=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-chromium-5prime")
UMI=feature("umi","umi","random")
ANTIBODY_FEATURE=feature("antibody_barcode","feature_barcode","whitelist",whitelist="ECCITE-seq Supplementary Table 3")
HASHTAG_FEATURE=feature("hashtag_barcode","sample_index","whitelist",whitelist="ECCITE-seq Supplementary Tables 4–5")
GUIDE_FEATURE=feature("guide","guide_barcode","unknown")
ANTIBODY_LONG_BARCODE="TAGAACCGTATCCGT"
ANTIBODY_SHORT_BARCODE="AATAGCGGAGCC"
HASHTAG_A_BARCODE="AGGACCATCCAA"

# Source oligos (MOESM3 "Oligonucleotide sequences", lines 385–429; Methods "ECCITE-seq on
# 10x Genomics instrument").  Phosphorothioate marks (*) and the /5AmMC12/ conjugation
# amine are chemistry, not bases, so only the bases enter the constructs below.
SI_PCR="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC"
TAG_TSO_END="CCCATATAAGAAA"          # 3' end of both tag oligos, anneals to the bead TSO
PROTEIN_TAG_HANDLE="CCTTGGCACCCGAGAATTCCA"
HASHTAG_HANDLE="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
RPI1_P7,RPI1_HANDLE="CAAGCAGAAGACGGCATACGAG","GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
D701S_P7,D701S_HANDLE="CAAGCAGAAGACGGCATACGAGAT","GTGACTGGAGTTCAGACGTGTGC"
NEXT_NST_P7,NEXT_NST_HANDLE,NEXT_NST_SCAFFOLD="CAAGCAGAAGACGGCATACGAGAT","GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG","TATTTCTAGCTCTAAAAC"
I7_TAG=feature("sample_index_i7","sample_index","unknown")
assert il.TRUSEQ_P5_FULL.startswith(SI_PCR)
assert RPI1_HANDLE.endswith(PROTEIN_TAG_HANDLE)          # RPI1 primes on the protein-tag handle
assert HASHTAG_HANDLE.startswith(D701S_HANDLE)           # D701_s primes on the hashtag handle
# Small-RNA RPI1 handle read primers: derived from the source RPI1 primer, not printed as
# sequencing primers by ECCITE sources (the run's primer mix is not stated).
SMALL_RNA_INDEX1=revcomp(RPI1_HANDLE)
SMALL_RNA_READ2=RPI1_HANDLE
_DERIVED="derived from RPI1 primer, Mimitou et al. 2019 MOESM3 line 419; sequencing primer not stated"
ADT_PRIMERS=(sp.TRUSEQ["R1"],sp.custom("Index 1 (i7)","small-RNA RPI handle index primer (derived)",SMALL_RNA_INDEX1,_DERIVED),
             sp.custom("Read 2","small-RNA RPI handle Read 2 primer (derived)",SMALL_RNA_READ2,_DERIVED))
HTO_PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
GDO_PRIMERS=(sp.TRUSEQ["R1"],sp.NEXTERA["I1"],sp.NEXTERA["R2"])
READ_LENGTHS={"Read 1":26,"Index 1 (i7)":8}   # Methods: 26 + 8 cycles; Read 2 39 (HiSeq) or 91 (NovaSeq)

def cell_side():
    """SI-PCR-derived P5/Read 1, bead cell barcode/UMI, and the TSO bases the tags pair with."""
    return [seg("P5 + Read 1 (SI-PCR primer)",SI_PCR,"p5"),
            seg("Read 1 end (bead oligo)",il.TRUSEQ_P5_FULL[len(SI_PCR):],"r1",inferred=True),
            seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL),seg("UMI","U"*10,"umi",placeholder=True,feature=UMI),
            seg("TSO 3′ bases (complement of tag 3′ end)",revcomp(TAG_TSO_END))]

def p7_side(handle,p7,name):
    return [seg(name,revcomp(handle),"r2"),seg("i7 reverse complement","I"*8,"cbc",placeholder=True,feature=I7_TAG),seg("P7 reverse complement",revcomp(p7),"p7")]

def tag_library(name,barcode_name,barcode,feat,handle,primer_handle,p7,primers):
    lib=Construct([*cell_side(),seg(barcode_name,barcode,"cbc",feature=feat),*p7_side(primer_handle,p7,"tag handle / i7 primer arm")],name=name)
    assert revcomp(handle) in lib.top()
    problems=sp.verify(lib,primers,required_roles=tuple(p.role for p in primers))
    if problems: raise ValueError(name+": "+"; ".join(problems))
    return lib,primers

GEX,GEX_P=truseq_library([seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL),seg("UMI","U"*10,"umi",placeholder=True,feature=UMI),
                          seg("TSO 3′ bases (complement of tag 3′ end)",revcomp(TAG_TSO_END)),seg("5-prime cDNA","X"*36,placeholder=True)],
                         "ECCITE-seq transcript library",dual_index=False,inferred_adapters=True)
ADT,ADT_P=tag_library("ECCITE-seq 15-base antibody-tag example","antibody barcode, Table 3 first entry",ANTIBODY_LONG_BARCODE,ANTIBODY_FEATURE,PROTEIN_TAG_HANDLE,RPI1_HANDLE,RPI1_P7,ADT_PRIMERS)
ADT_SHORT,ADT_SHORT_P=tag_library("ECCITE-seq 12-base antibody-tag example","antibody barcode, Table 3 short entry",ANTIBODY_SHORT_BARCODE,ANTIBODY_FEATURE,PROTEIN_TAG_HANDLE,RPI1_HANDLE,RPI1_P7,ADT_PRIMERS)
HTO,HTO_P=tag_library("ECCITE-seq sample-hashtag example","Hashtag_A barcode",HASHTAG_A_BARCODE,HASHTAG_FEATURE,HASHTAG_HANDLE,HASHTAG_HANDLE,D701S_P7,HTO_PRIMERS)
GDO=Construct([*cell_side(),seg("sgRNA protospacer","G"*20,"cbc",placeholder=True,feature=GUIDE_FEATURE),
               seg("scaffold 5′ region (Next_nst_x annealing site)",revcomp(NEXT_NST_SCAFFOLD)),
               *p7_side(NEXT_NST_HANDLE,NEXT_NST_P7,"Nextera Read 2 / i7 primer arm")],name="ECCITE-seq guide-derived oligo library")
GDO_P=GDO_PRIMERS
if sp.verify(GDO,GDO_P,required_roles=tuple(p.role for p in GDO_P)): raise ValueError("guide library primers do not land")
_READS=" Methods: 26 cycles Read 1 (cell barcode + UMI), 8 cycles index, and 39 (HiSeq 2500 rapid) or 91 (NovaSeq S2) cycles Read 2; Read 2 length is not fixed here."
FINAL_LIBRARIES=(
    ("Final transcript library",GEX,GEX_P,"INFERRED — 10x V(D)J-protocol transcript library: 5′ cDNA follows cell barcode, UMI and TSO bases; the 10x fragmentation/ligation adapters are not inspected and are drawn as a standard single-index TruSeq shell.","The transcript library's primer geometry is computed on the inferred duplex."+_READS),
    ("Final antibody-tag library, 15-base example",ADT,ADT_P,"SI-PCR × RPI1 product of a protein tag: TSO bases, the first Supplementary Table 3 barcode TAGAACCGTATCCGT, then the small-RNA RA3/RPI handle and i7. The last nine Read 1 bases (bead oligo) are INFERRED.","Index 1 and Read 2 primers are derived from the RPI1 handle."+_READS),
    ("Final antibody-tag library, 12-base example",ADT_SHORT,ADT_SHORT_P,"The last three Supplementary Table 3 antibody barcodes have 12 bases; AATAGCGGAGCC is one example. Same SI-PCR × RPI1 layout; bead-oligo Read 1 end INFERRED.","This shorter antibody-tag variant uses the same library preparation."+_READS),
    ("Final sample-hashtag library",HTO,HTO_P,"SI-PCR × D701_s product: TSO bases, mouse Hashtag_A AGGACCATCCAA (Supplementary Table 4), then the TruSeq Read 2 handle printed in the hashtag oligo. Bead-oligo Read 1 end INFERRED.","The hashtag library's TruSeq primer geometry is computed on the finished duplex."+_READS),
    ("Final guide-derived oligo library",GDO,GDO_P,"SI-PCR × Next_nst_x (v4) product: TSO bases, the sgRNA protospacer and the scaffold 5′ bases that Next_nst_x anneals to, then the Nextera Read 2 handle and i7. Bead-oligo Read 1 end INFERRED.","The guide library's Nextera index and Read 2 primers are computed on the finished duplex."+_READS),)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("cells","cells carrying poly(A) RNA, antibody/hashtag oligos and sgRNA",MaterialKind.BIOLOGICAL,state=molecule("ECCITE-seq cells",(role("poly(A) RNA"),role("antibody and hashtag oligos"),role("sgRNA"))),primary=True)
    droplets=g.input("droplets","10x 5-prime GEM partitioning condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","occupied GEM products",MaterialKind.MOLECULAR,state=cells.state)
    failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
    g.react("partition","Partition cells with barcoded gel beads",StepKind.PARTITION,((cells,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False)
    g.discard(failed,"does not retain a resolvable cell barcode")
    transcript=g.product("transcript_branch","poly(A) RNA in the GEM",MaterialKind.MOLECULAR,state=occupied.state)
    tag=g.product("tag_branch","antibody and hashtag oligos in the GEM",MaterialKind.MOLECULAR,state=occupied.state)
    guide=g.product("guide_branch","sgRNA in the GEM",MaterialKind.MOLECULAR,state=occupied.state)
    g.react("split","Track the co-encapsulated transcript, tag and sgRNA templates in one RT reaction",StepKind.SPLIT,((occupied,InputRole.TEMPLATE),),(transcript,tag,guide),changes_molecule=False)
    bead=g.input("bead_oligo","10x 5′ gel-bead cell-barcode/UMI template-switch oligo",MaterialKind.OLIGO)
    rt=g.input("rt","10x 5′ RT reagents (template-switching reverse transcriptase, poly(dT) RT primer)",MaterialKind.CATALYST)
    cdna=g.product("cdna","cell-barcoded 5-prime cDNA",MaterialKind.MOLECULAR,state=molecule("ECCITE transcript cDNA",(role("cell barcode"),role("UMI"),role("5-prime cDNA"))))
    g.react("transcript_rt","Reverse-transcribe and template-switch transcript RNA",StepKind.REVERSE_TRANSCRIPTION,((transcript,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
    copied_tag=g.product("copied_tag","cell-barcoded antibody/hashtag products",MaterialKind.MOLECULAR,state=molecule("ECCITE tag precursor",(role("cell barcode"),role("UMI"),role("antibody or hashtag barcode"))))
    g.react("copy_tags","Anneal the tag 3′ end (CCCATATAAGAAA) to the bead TSO and extend during RT",StepKind.EXTENSION,((tag,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied_tag,),
            source=SOURCE+" Results para. 2 and MOESM3 lines 385–396")
    guide_primer=g.input("guide_rt_primer","sgRNA-scaffold RT primer gd_RT_v4 / sg_RT_v4, 12 pmol: AGCAAGTGAGAAGCATCGTGTCAAAGCACCGACTCGGTGCCAC",MaterialKind.OLIGO,source=SOURCE+" Methods (ECCITE-seq on 10x instrument) and MOESM3 lines 104, 399")
    guide_cdna=g.product("guide_cdna","template-switched guide-derived cDNA",MaterialKind.MOLECULAR,state=molecule("ECCITE guide cDNA",(role("cell barcode"),role("UMI"),role("guide protospacer"),role("scaffold"))))
    g.react("guide_rt","Reverse-transcribe sgRNA from its scaffold primer and template-switch",StepKind.REVERSE_TRANSCRIPTION,((guide,InputRole.TEMPLATE),(guide_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(guide_cdna,))
    rt_products=g.product("rt_products","pooled GEM RT products",MaterialKind.MOLECULAR,state=molecule("ECCITE RT products",(role("transcript cDNA"),role("tag products"),role("guide cDNA"))))
    g.react("recover_rt","Recover all barcoded RT products together",StepKind.POOL,((cdna,InputRole.TEMPLATE),(copied_tag,InputRole.TEMPLATE),(guide_cdna,InputRole.TEMPLATE)),(rt_products,),changes_molecule=False)
    additives=g.input("tag_additives","spiked additive primers: hashtag GTGACTGGAGTTCAGACGTGTGCTC (1 pmol), guide-tag AGCAAGTGAGAAGCATCGTGTC (1 pmol), protein-tag CCTTGGCACCCGAGAATTCC (2 pmol)",MaterialKind.OLIGO,source=SOURCE+" Methods (ECCITE-seq on 10x instrument)")
    cdna_primers=g.input("cdna_primers","10x V(D)J cDNA amplification primers (not printed in ECCITE sources)",MaterialKind.OLIGO,inferred=True,source=SOURCE+" (follows the 10x V(D)J protocol)")
    cdna_pol=g.input("cdna_polymerase","10x cDNA amplification master mix",MaterialKind.CATALYST)
    amplified=g.product("amplified","amplified cDNA plus tag and guide products",MaterialKind.MOLECULAR,state=rt_products.state)
    g.react("cdna_amplification","Amplify cDNA with spiked tag and guide additive primers",StepKind.PCR,((rt_products,InputRole.TEMPLATE),(cdna_primers,InputRole.STRUCTURAL),(additives,InputRole.STRUCTURAL),(cdna_pol,InputRole.CATALYST)),(amplified,),
            geometry=PCRGeometry.unavailable(Construct([seg("mixed cDNA, tag and guide products","X"*30,placeholder=True)],name="pooled RT products"),reason="The 10x cDNA primers are not printed in ECCITE sources; the additive primers are printed but their per-species pairing with the 10x forward primer is not drawn here."))
    spri06=g.input("spri_06x","0.6X SPRI beads",MaterialKind.SELECTION)
    large=g.product("large_cdna","large mRNA-derived cDNA (bead-bound)",MaterialKind.MOLECULAR,state=cdna.state)
    small=g.product("small_fraction","protein-tag, hashtag and guide-tag fraction (supernatant)",MaterialKind.MOLECULAR,state=molecule("ECCITE small-product fraction",(role("tag products"),role("guide cDNA"))))
    g.react("spri_size_split","Separate large cDNA (beads) from short tag/guide products (supernatant) with 0.6X SPRI",StepKind.SELECTION,((amplified,InputRole.TEMPLATE),(spri06,InputRole.SELECTION)),(large,small))
    spri20=g.input("spri_20x","SPRI beads to 2.0X, then a second 2.0X round",MaterialKind.SELECTION)
    small_clean=g.product("small_clean","cleaned tag/guide fraction",MaterialKind.MOLECULAR,state=small.state)
    excess=g.product("excess_oligos","excess single-stranded oligos from cDNA amplification",MaterialKind.WASTE)
    g.react("spri_cleanup","Bind short products at 2.0X SPRI twice, ethanol-wash and elute",StepKind.CLEANUP,((small,InputRole.TEMPLATE),(spri20,InputRole.SELECTION)),(small_clean,excess))
    g.discard(excess,"removed by the two 2.0X SPRI rounds")
    copied_adt=g.product("copied_adt","antibody-tag PCR aliquot",MaterialKind.MOLECULAR,
                         state=molecule("ECCITE antibody-tag precursor",(role("cell barcode"),role("UMI"),role("antibody barcode"))))
    copied_hto=g.product("copied_hto","hashtag PCR aliquot",MaterialKind.MOLECULAR,
                         state=molecule("ECCITE hashtag precursor",(role("cell barcode"),role("UMI"),role("12-base hashtag barcode"))))
    guide_aliquot=g.product("guide_aliquot","guide-tag PCR aliquot",MaterialKind.MOLECULAR,state=guide_cdna.state)
    g.react("split_tag_classes","Set up separate protein-tag, hashtag and guide-tag PCRs from the cleaned fraction",StepKind.SPLIT,
            ((small_clean,InputRole.TEMPLATE),),(copied_adt,copied_hto,guide_aliquot),changes_molecule=False,
            source=SOURCE+" Methods: SI-PCR with RPI-x, D7xx_s or Next_nst_x")
    copied_adt_long=g.product("copied_adt_long","antibody-tag products with 15-base barcodes",MaterialKind.MOLECULAR,state=copied_adt.state)
    copied_adt_short=g.product("copied_adt_short","antibody-tag products with 12-base barcodes",MaterialKind.MOLECULAR,state=copied_adt.state)
    g.react("split_antibody_lengths","Resolve published antibody barcode length variants",StepKind.SPLIT,
            ((copied_adt,InputRole.TEMPLATE),),(copied_adt_long,copied_adt_short),changes_molecule=False,
            source=SOURCE+" Supplementary Table 3, final three 12-base barcodes")
    reviewed=["cells","droplets","partition","split","bead_oligo","rt","transcript_rt","copy_tags","guide_rt_primer","guide_rt","recover_rt","tag_additives","cdna_primers","cdna_polymerase","cdna_amplification","spri_06x","spri_size_split","spri_20x","spri_cleanup","split_tag_classes","split_antibody_lengths"]
    finish_pcr_branch(g,"gex",large,GEX,GEX_P,reviewed,label="ECCITE-seq transcript",read_lengths=READ_LENGTHS,limitation="Processed per the 10x Single Cell V(D)J protocol (fragmentation, adapter ligation, sample-index PCR), which was not inspected in this pass.")
    finish_pcr_branch(g,"adt",copied_adt_long,ADT,ADT_P,reviewed,label="ECCITE-seq 15-base antibody-tag",read_lengths=READ_LENGTHS)
    finish_pcr_branch(g,"adt_short",copied_adt_short,ADT_SHORT,ADT_SHORT_P,reviewed,label="ECCITE-seq 12-base antibody-tag",read_lengths=READ_LENGTHS)
    finish_pcr_branch(g,"hto",copied_hto,HTO,HTO_P,reviewed,label="ECCITE-seq sample-hashtag",read_lengths=READ_LENGTHS)
    si_pcr=g.input("guide_enrich_si_pcr","SI-PCR primer "+SI_PCR,MaterialKind.OLIGO,source=SOURCE+" MOESM3 line 405")
    gd_add=g.input("guide_enrich_gd_add","gd_add_v4 (named in the MOESM3 guide PCR; presumed the listed guide-tag additive AGCAAGTGAGAAGCATCGTG*T*C)",MaterialKind.OLIGO,inferred=True,source=SOURCE+" MOESM3 guide-tag library step 1 and oligo list")
    kapa=g.input("guide_enrich_polymerase","2x KAPA HiFi master mix",MaterialKind.CATALYST)
    guide_enriched=g.product("guide_enriched","enriched guide-tag product (~220 bp), 1.6X SPRI purified",MaterialKind.MOLECULAR,state=guide_cdna.state)
    g.react("guide_enrichment","Enrichment PCR of guide tags with SI-PCR and gd_add_v4 (8–10 cycles), then 1.6X SPRI",StepKind.PCR,
            ((guide_aliquot,InputRole.TEMPLATE),(si_pcr,InputRole.STRUCTURAL),(gd_add,InputRole.STRUCTURAL),(kapa,InputRole.CATALYST)),(guide_enriched,),
            geometry=PCRGeometry.unavailable(Construct([seg("guide cDNA before indexing","X"*30,placeholder=True)],name="guide-tag aliquot"),reason="Supplementary protocol (March 2019) step 7.1; Methods call further enrichment optional. Primer landing on the pre-index guide cDNA is not drawn."),
            source=SOURCE+" MOESM3 'Guide-tag library preparation' step 1")
    reviewed.extend(["guide_enrich_si_pcr","guide_enrich_gd_add","guide_enrich_polymerase","guide_enrichment"])
    finish_pcr_branch(g,"guide",guide_enriched,GDO,GDO_P,reviewed,label="ECCITE-seq guide-derived oligo",read_lengths=READ_LENGTHS)
    g.declare_reviewed_scope(reviewed,SOURCE)
    return g.require_publishable()
