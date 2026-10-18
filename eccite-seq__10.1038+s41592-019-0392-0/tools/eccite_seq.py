"""ECCITE-seq 5-prime multimodal tag and direct-guide capture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,feature,revcomp
from pcr import PCRGeometry
from protocol_graph import Disposition,InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind,Terminal
from protocol_model import molecule, role

seg=Segment

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
        "and the transcript library's historical vendor route remains unassigned. The protein-tag library's Index 1 "
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
READ_LENGTHS={"Read 1":26,"Index 1 (i7)":8,"Read 2":39} # selected published HiSeq route; NovaSeq alternative R2=91

def cell_side():
    """SI-PCR-derived P5/Read 1, bead cell barcode/UMI, and the TSO bases the tags pair with."""
    return [seg("P5 + Read 1 (SI-PCR primer)",SI_PCR,"p5"),
            seg("Read 1 end (bead oligo)",il.TRUSEQ_P5_FULL[len(SI_PCR):],"r1",inferred=True),
            seg("cell barcode","B"*16,"cbc",placeholder=True,feature=CELL,inferred=True),seg("UMI","U"*10,"umi",placeholder=True,feature=UMI,inferred=True),
            seg("TSO 3′ bases (complement of tag 3′ end)",revcomp(TAG_TSO_END))]

def p7_side(handle,p7,name):
    return [seg(name,revcomp(handle),"r2"),seg("i7 reverse complement","I"*8,"cbc",placeholder=True,feature=I7_TAG),seg("P7 reverse complement",revcomp(p7),"p7")]

def tag_library(name,barcode_name,barcode,feat,handle,primer_handle,p7,primers):
    lib=Construct([*cell_side(),seg(barcode_name+" reverse complement",revcomp(barcode),"cbc",feature=feat),*p7_side(primer_handle,p7,"tag handle / i7 primer arm")],name=name)
    # DNA tags prime extension on the bead TSO. The Read1-facing strand is
    # the complement of the complete ordered handle–barcode–TSO-end tag.
    # Checking the whole copied interval rejects a mixed-polarity barcode.
    if revcomp(handle+barcode+TAG_TSO_END) not in lib.top():
        raise ValueError(name+": copied interval does not preserve ordered-tag polarity")
    problems=sp.verify(lib,primers,required_roles=tuple(p.role for p in primers))
    if problems: raise ValueError(name+": "+"; ".join(problems))
    return lib,primers

GEX=Construct([role("10x historical V(D)J transcript library; native outer geometry unassigned",inferred=True)],name="ECCITE manufacturer-prepared transcript library (bounded)")
GEX_P=()

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
    ("Final transcript library",GEX,GEX_P,"Historical 10x V(D)J-protocol transcript library; source-independent canonical adapter/index defaults are not assigned.","The historical manufacturer transcript-library route and physical primer geometry remain unassigned."+_READS),
    ("Final antibody-tag library, 15-base example",ADT,ADT_P,"SI-PCR × RPI1 product of a protein tag: TSO bases, the reverse complement of the first ordered Supplementary Table 3 barcode TAGAACCGTATCCGT, then the small-RNA RA3/RPI handle and i7. The last nine Read 1 bases (bead oligo) are INFERRED.","Index 1 and Read 2 primers are derived from the RPI1 handle."+_READS),
    ("Final antibody-tag library, 12-base example",ADT_SHORT,ADT_SHORT_P,"The last three Supplementary Table 3 antibody barcodes have 12 bases; ordered AATAGCGGAGCC is one example, shown reverse-complemented in the copied strand. Same SI-PCR × RPI1 layout; bead-oligo Read 1 end INFERRED.","This shorter antibody-tag variant uses the same library preparation."+_READS),
    ("Final sample-hashtag library",HTO,HTO_P,"SI-PCR × D701_s product: TSO bases, the copy of ordered mouse Hashtag_A AGGACCATCCAA (Supplementary Table4), shown reverse-complemented, then the TruSeq Read 2 handle printed in the hashtag oligo. Bead-oligo Read 1 end INFERRED.","The hashtag library's TruSeq primer geometry is computed on the finished duplex."+_READS),
    ("Final guide-derived oligo library",GDO,GDO_P,"SI-PCR × Next_nst_x (v4) product: TSO bases, the sgRNA protospacer and the scaffold 5′ bases that Next_nst_x anneals to, then the Nextera Read 2 handle and i7. Bead-oligo Read 1 end INFERRED.","The guide library's Nextera index and Read 2 primers are computed on the finished duplex."+_READS),)

READ_SPECS={"Final transcript library":(
    sp.read_spec("Read 1",26,"10x5′ cell barcode/UMI; native primer geometry unassigned","extract"),
    sp.read_spec("Index 1 (i7)",8,"sample index; native species unassigned","demultiplex"),
    sp.read_spec("Read 2",39,"transcript; selected source HiSeq route, NovaSeq alternative91cycles","include"))}
SEQUENCING_UNAVAILABLE={"Final transcript library":"Historical 10xV(D)J manufacturer preparation not reconstructed; source read recipe is retained separately."}


def independent_populations(label,*groups):
    sc=Scene()
    for index,parts in enumerate(groups):
        sc.strand(f"population{index+1}",parts,label=f"{label}: population{index+1}")
    return sc


def identifier_copy(label,content,identifier=None):
    return molecule(label,(role("cell barcode",molecular_feature=CELL),role("UMI",molecular_feature=UMI),role(content,molecular_feature=identifier)))


def short_population(label):
    return independent_populations(label,
        (role("amplified protein-tag population; native ends unassigned",molecular_feature=ANTIBODY_FEATURE),),
        (role("amplified hashtag population; native ends unassigned",molecular_feature=HASHTAG_FEATURE),),
        (role("amplified guide-tag population; native ends unassigned",molecular_feature=GUIDE_FEATURE),))


def primer_state(label,prefix,arm,index=False):
    parts=[Segment(label+"5′part",prefix)] if prefix else []
    if index:parts.append(Segment("i7 sample index","N"*8,placeholder=True,length_bp=8,feature=I7_TAG))
    parts.append(Segment(label+"3′part",arm,note="terminal phosphorothioates preserved in copyable source oligos"))
    return molecule(label,tuple(parts))


SOURCE_OLIGO_EXCEPTIONS={"next_primer":"Published MOESM3 line429 prints the complete Next_nst_x family with x8 and terminal *A*C. Bundle scraper/find does not return that complete long family string; local independent compact-source test checks every fixed base and the complete x8 interval against that exact original row. No altered sequence is substituted."}


def oligos():
    return ("SI-PCR, source phosphorothioates:5′-"+SI_PCR[:-2]+"*T*C-3′",
        "Protein-tag additive:5′-CCTTGGCACCCGAGAATT*C*C-3′",
        "Hashtag additive:5′-GTGACTGGAGTTCAGACGTGTGC*T*C-3′",
        "Guide-tag additive:5′-AGCAAGTGAGAAGCATCGTG*T*C-3′",
        "RPIx source family:5′-"+RPI1_P7+"xxxxxxxx"+RPI1_HANDLE[:-2]+"*C*A-3′",
        "D7xx_s source family:5′-"+D701S_P7+"xxxxxxxx"+D701S_HANDLE[:-2]+"*G*C-3′",
        "Next_nst_x source family:5′-"+NEXT_NST_P7+"xxxxxxxx"+NEXT_NST_HANDLE+NEXT_NST_SCAFFOLD[:-2]+"*A*C-3′")


def protocol_graph():
    source=SOURCE+" published Methods and March2019 MOESM3 protocol"
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input("cells","cells carrying coencapsulated RNA and tag molecules",MaterialKind.BIOLOGICAL,
        state=independent_populations("cell container",(role("poly(A) RNA"),),(role("antibody oligos"),),(role("hashtag oligos"),),(role("sgRNA"),)),primary=True)
    droplets=g.input("droplets","historical10x5′GEM condition",MaterialKind.CONDITION)
    occupied=g.product("occupied","occupied GEM reaction",MaterialKind.MOLECULAR,state=cells.state)
    failed=g.product("failed","empty/ambiguous GEM reaction material",MaterialKind.MOLECULAR)
    g.react("partition","Partition cells and barcoded beads",StepKind.PARTITION,((cells,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False)
    g.mark_terminal(Terminal(failed.id,Disposition.UNRESOLVED,"Unresolved cell assignments are computational; source does not establish physical exclusion of these molecules."))
    bead=g.input("bead_oligo","10x5′ cell-barcode/UMI bead-associated TSO; native species unassigned",MaterialKind.OLIGO)
    rt=g.input("rt","10x5′RT mix with soluble poly(dT) primer",MaterialKind.CATALYST)
    guide_rt=g.input("guide_rt_primer","gd_RT_v4 / listed sg_RT_v4 (12pmol)",MaterialKind.OLIGO,
        state=Construct([Segment("sgRTv4","AGCAAGTGAGAAGCATCGTGTCAAAGCACCGACTCGGTGCCAC")]),source=source)
    cdna=g.product("cdna","single-strand template-switched transcript cDNA",MaterialKind.MOLECULAR,
        state=identifier_copy("transcript cDNA","5′transcript cDNA; distal RT end unassigned"))
    tag=g.product("copied_tag","single-strand protein/hashtag barcoded copies",MaterialKind.MOLECULAR,
        state=independent_populations("tag copies",tuple(identifier_copy("protein copy","protein-tag barcode",ANTIBODY_FEATURE).strands["molecule"].segs),tuple(identifier_copy("hashtag copy","hashtag barcode",HASHTAG_FEATURE).strands["molecule"].segs)))
    guide=g.product("guide_cdna","single-strand template-switched guide cDNA",MaterialKind.MOLECULAR,
        state=identifier_copy("guide cDNA","guide protospacer and copied scaffold",GUIDE_FEATURE))
    g.react("gem_rt","One GEM RT incubation: transcript/guide RT and tag annealing/extension",StepKind.REVERSE_TRANSCRIPTION,
        ((occupied,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(guide_rt,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,tag,guide),source=source)
    # Bookkeeping joins independent molecules into one recovered tube, not a strand.
    recovered=g.product("rt_products","recovered mixed GEM RT products",MaterialKind.MOLECULAR,
        state=independent_populations("recovered RT products",tuple(cdna.state.strands["molecule"].segs),
            *[tuple(st.segs) for st in tag.state.strands.values()],tuple(guide.state.strands["molecule"].segs)))
    g.react("recover_rt","Recover the RT molecular populations together",StepKind.POOL,
        ((cdna,InputRole.TEMPLATE),(tag,InputRole.TEMPLATE),(guide,InputRole.TEMPLATE)),(recovered,),changes_molecule=False)
    cdna_primers=g.input("cdna_primers","historical10xV(D)J cDNA primer mixture; ordered species unassigned",MaterialKind.OLIGO)
    additives=g.input("tag_additives","hashtag1pmol,guide-tag1pmol,protein-tag2pmol additive oligos",MaterialKind.OLIGO,
        state=independent_populations("source additive oligos",(Segment("hashtag","GTGACTGGAGTTCAGACGTGTGCTC",note="terminal phosphorothioates"),),
            (Segment("guide","AGCAAGTGAGAAGCATCGTGTC",note="terminal phosphorothioates"),),
            (Segment("protein","CCTTGGCACCCGAGAATTCC",note="terminal phosphorothioates"),)),source=source)
    cdna_pol=g.input("cdna_polymerase","10x cDNA amplification master mix",MaterialKind.CATALYST)
    amplified=g.product("amplified","amplified transcript and short-product populations",MaterialKind.MOLECULAR,
        state=independent_populations("post-cDNA-PCR mixture",(role("amplified transcript cDNA; native ends unassigned"),),
            *[tuple(st.segs) for st in short_population("short products").strands.values()]))
    g.react("cdna_amplification","Amplify cDNA with source additive primers",StepKind.PCR,
        ((recovered,InputRole.TEMPLATE),(cdna_primers,InputRole.STRUCTURAL),(additives,InputRole.STRUCTURAL),(cdna_pol,InputRole.CATALYST)),(amplified,),
        geometry=PCRGeometry.unavailable(recovered.state,reason="Actual input is independent ss RT populations. Historical10x primer strands and per-species first-cycle geometry are unassigned; no completed library substitutes for them."),source=source)
    spri=g.input("spri_06x","0.6×SPRI",MaterialKind.SELECTION)
    large=g.product("large_cdna","bead-retained large transcript cDNA fraction",MaterialKind.MOLECULAR,state=molecule("large cDNA",(role("amplified transcript fraction; native ends unassigned"),)))
    small=g.product("small_fraction","saved protein/hashtag/guide supernatant",MaterialKind.MOLECULAR,state=short_population("saved short fraction"))
    g.react("spri_size_split","Retain large cDNA on beads and short products in supernatant",StepKind.SELECTION,
        ((amplified,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(large,small),source=source)
    clean_beads=g.input("spri_20x","add1.4×to reach2.0×, then second2.0×SPRI",MaterialKind.SELECTION)
    small_clean=g.product("small_clean","cleaned mixed short-product fraction",MaterialKind.MOLECULAR,state=small.state)
    waste=g.product("excess_oligos","removed2×SPRI supernatants/ethanol and excess oligos",MaterialKind.WASTE)
    g.react("spri_cleanup","Purify the mixed short-product fraction twice at2×SPRI",StepKind.CLEANUP,
        ((small,InputRole.TEMPLATE),(clean_beads,InputRole.SELECTION)),(small_clean,waste),changes_molecule=False,source=source)
    g.discard(waste,"MOESM3 removes supernatants/ethanol washes; not a modality-specific molecular purification")
    adt_in=g.product("adt_aliquot","protein-tag PCR aliquot of mixed short fraction",MaterialKind.MOLECULAR,state=small_clean.state)
    hto_in=g.product("hto_aliquot","hashtag PCR aliquot of mixed short fraction",MaterialKind.MOLECULAR,state=small_clean.state)
    guide_in=g.product("guide_aliquot","guide PCR aliquot of mixed short fraction",MaterialKind.MOLECULAR,state=small_clean.state)
    g.react("split_tag_classes","Aliquot the same mixed fraction for three primer-selected PCRs",StepKind.SPLIT,
        ((small_clean,InputRole.TEMPLATE),),(adt_in,hto_in,guide_in),changes_molecule=False,source=source)
    manufacturer=g.input("gex_manufacturer","historical10xV(D)J transcript library workflow",MaterialKind.CONDITION)
    gex=g.product("gex_library","manufacturer-prepared transcript library",MaterialKind.LIBRARY,state=GEX)
    g.react("gex_library_prep","Prepare transcript library by historical manufacturer protocol",StepKind.MAGIC,
        ((large,InputRole.TEMPLATE),(manufacturer,InputRole.CONDITION)),(gex,),reason="Historical physical fragmentation/adapter/index PCR recipe is not reconstructed; later vendor revisions and canonical endpoints are not substituted.",source=source)
    g.sequence(gex,SequencingPlan(library=GEX,acquired_roles=("Read 1","Index 1 (i7)","Read 2"),read_lengths=READ_LENGTHS,
        limitation=SEQUENCING_UNAVAILABLE["Final transcript library"]))
    si=g.input("si_pcr","source SI-PCR (two terminal phosphorothioates)",MaterialKind.OLIGO,state=primer_state("SI-PCR","",SI_PCR),source=source)
    rpi=g.input("rpi_primer","source RPIx (RPI1 example family; i7=8positions)",MaterialKind.OLIGO,
        state=primer_state("RPIx",RPI1_P7,RPI1_HANDLE,True),source=source)
    d7=g.input("d7_primer","source D7xx_s",MaterialKind.OLIGO,state=primer_state("D7xx_s",D701S_P7,D701S_HANDLE,True),source=source)
    nextp=g.input("next_primer","source Next_nst_x v4",MaterialKind.OLIGO,state=primer_state("Next_nst_x",NEXT_NST_P7,NEXT_NST_HANDLE+NEXT_NST_SCAFFOLD,True),source=source)
    # Local source-defined PCRs have distinct conditions and explicit shared-pool inputs.
    adt_pol=g.input("adt_polymerase","2×KAPA HiFi;100µl,SI/RPIx25pmol each;6–10cycles at60°C",MaterialKind.CATALYST)
    adt_raw=g.product("adt_raw","15-base representative species in shared antibody PCR pool",MaterialKind.MOLECULAR,state=ADT)
    adt_short_raw=g.product("adt_short_raw","12-base representative species in same antibody PCR pool",MaterialKind.MOLECULAR,state=ADT_SHORT)
    g.react("adt_pcr","Protein-tag PCR: SI-PCR/RPIx; representative lengths coexist",StepKind.PCR,
        ((adt_in,InputRole.TEMPLATE),(si,InputRole.STRUCTURAL),(rpi,InputRole.STRUCTURAL),(adt_pol,InputRole.CATALYST)),(adt_raw,adt_short_raw),
        geometry=PCRGeometry.unavailable(adt_in.state,reason="Actual aliquot contains all short-product classes. Native pre-PCR ends are not reconstructed;15/12-base variants are species in one PCR, not physically separated by barcode length."),source=source)
    hto_pol=g.input("hto_polymerase","2×KAPA HiFi;100µl,SI/D7xx_s25pmol each;8–12cycles at64°C",MaterialKind.CATALYST)
    hto_raw=g.product("hto_raw","amplified hashtag library",MaterialKind.MOLECULAR,state=HTO)
    g.react("hto_pcr","Hashtag PCR with SI-PCR/D7xx_s",StepKind.PCR,
        ((hto_in,InputRole.TEMPLATE),(si,InputRole.STRUCTURAL),(d7,InputRole.STRUCTURAL),(hto_pol,InputRole.CATALYST)),(hto_raw,),
        geometry=PCRGeometry.unavailable(hto_in.state,reason="Actual mixed short-product aliquot; native hashtag first-cycle geometry unassigned."),source=source)
    gd_add=g.input("guide_enrich_gd_add","gd_add_v4 named in MOESM3 enrichment; precise listed-member assignment unresolved",MaterialKind.OLIGO)
    enrich_pol=g.input("guide_enrich_polymerase","KAPA HiFi100µl,8–10cycles at64°C",MaterialKind.CATALYST)
    enriched=g.product("guide_enriched","guide-enriched PCR pool; native ends unassigned",MaterialKind.MOLECULAR,
        state=short_population("guide-enrichment product; residual non-guide recovery unassigned"))
    g.react("guide_enrichment","Optional source guide enrichment with SI-PCR/gd_add_v4 (selected)",StepKind.PCR,
        ((guide_in,InputRole.TEMPLATE),(si,InputRole.STRUCTURAL),(gd_add,InputRole.STRUCTURAL),(enrich_pol,InputRole.CATALYST)),(enriched,),
        geometry=PCRGeometry.unavailable(guide_in.state,reason="Actual mixed aliquot; gd_add_v4 species and native first-cycle geometry unassigned."),source=source)
    enrich_beads=g.input("guide_enrich_beads","1.6×SPRI,20µl water",MaterialKind.SELECTION)
    enrich_clean=g.product("guide_enrich_clean","cleaned guide-enriched pool",MaterialKind.MOLECULAR,state=enriched.state)
    enrich_waste=g.product("guide_enrich_waste","1.6×SPRI supernatant and washes",MaterialKind.WASTE)
    g.react("guide_enrich_cleanup","Purify guide enrichment with1.6×SPRI",StepKind.CLEANUP,
        ((enriched,InputRole.TEMPLATE),(enrich_beads,InputRole.SELECTION)),(enrich_clean,enrich_waste),changes_molecule=False,source=source)
    g.discard(enrich_waste,"MOESM3 physical supernatant/ethanol removal")
    guide_pol=g.input("guide_polymerase","KAPA HiFi100µl,SI/Next_nst_x25pmol each;6–10cycles at54°C",MaterialKind.CATALYST)
    guide_raw=g.product("guide_raw","indexed guide library PCR product",MaterialKind.MOLECULAR,state=GDO)
    g.react("guide_pcr","Index guide library using SI-PCR/Next_nst_x v4",StepKind.PCR,
        ((enrich_clean,InputRole.TEMPLATE),(si,InputRole.STRUCTURAL),(nextp,InputRole.STRUCTURAL),(guide_pol,InputRole.CATALYST)),(guide_raw,),
        geometry=PCRGeometry.unavailable(enrich_clean.state,reason="Actual enriched pool, not expected indexed GDO duplex; full native first-cycle geometry unassigned."),source=source)
    # Source gives the same1.6× physical cleanup for the separate tag PCRs.
    for key,inputs,outputs in (
        ("adt",(adt_raw,adt_short_raw),(("adt_library",ADT,ADT_P),("adt_short_library",ADT_SHORT,ADT_SHORT_P))),
        ("hto",(hto_raw,),(("hto_library",HTO,HTO_P),)),
        ("guide",(guide_raw,),(("guide_library",GDO,GDO_P),))):
        cleanup=g.input(key+"_cleanup_beads","source1.6×SPRI and20µl water eluate",MaterialKind.SELECTION)
        products=[g.product(name,"purified "+name+" (representative species)",MaterialKind.LIBRARY,state=lib) for name,lib,_ in outputs]
        discard=g.product(key+"_cleanup_waste","1.6×supernatant and ethanol washes",MaterialKind.WASTE)
        g.react(key+"_cleanup","Purify "+key+" PCR pool with1.6×SPRI",StepKind.CLEANUP,
            tuple((item,InputRole.TEMPLATE) for item in inputs)+((cleanup,InputRole.SELECTION),),(*products,discard),changes_molecule=False,source=source)
        g.discard(discard,"Source physical supernatant/wash removals; no15/12-base species separation")
        for product,(_,lib,primers) in zip(products,outputs):
            g.sequence(product,SequencingPlan(library=lib,primers=primers,acquired_roles=("Read 1","Index 1 (i7)","Read 2"),read_lengths=READ_LENGTHS,
                limitation="Selected published HiSeq26/8/39 recipe; NovaSeq alternative R2=91. Bead-derived endpoint and compatible small-RNA sequencing sites are inferred, not an assigned physical reagent mix."))
    g.declare_reviewed_scope(("cells","droplets","partition","bead_oligo","rt","guide_rt_primer","gem_rt","recover_rt","cdna_primers","tag_additives","cdna_polymerase","cdna_amplification","spri_06x","spri_size_split","spri_20x","spri_cleanup","split_tag_classes","gex_manufacturer","gex_library_prep","gex_library","si_pcr","rpi_primer","d7_primer","next_primer","adt_polymerase","adt_pcr","hto_polymerase","hto_pcr","guide_enrich_gd_add","guide_enrich_polymerase","guide_enrichment","guide_enrich_beads","guide_enrich_cleanup","guide_polymerase","guide_pcr","adt_cleanup_beads","adt_cleanup","adt_library","adt_short_library","hto_cleanup_beads","hto_cleanup","hto_library","guide_cleanup_beads","guide_cleanup","guide_library"),source)
    return g.require_publishable()
