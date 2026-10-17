"""SUM-seq sample pre-indexing and overloaded droplet RNA/ATAC profiling."""
import re
from batch_ngs import seg
from chemdraw import Construct, feature, revcomp
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import (Highlight, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)
TITLE="SUM-seq — sample-indexed, overloaded single-cell RNA and ATAC"
NOTES="01_sum-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41592-025-02700-8">Llorens-Rico et al., <i>Nature Methods</i> (2025)</a>.'
SUMMARY="ATAC fragments and RNA are separately sample-indexed in fixed nuclei, then a second 10x droplet barcode resolves overloaded nuclei; the sample-index–droplet-index pair identifies each cell across both modalities."
CAVEAT="The final libraries are assembled from the Supplementary Table 1 oligos. <inf>Inferred:</inf> the 10x scATAC gel-bead oligo (P5 remainder and 16-nt droplet barcode in front of Nextera Read 1) is not printed in the article or Table 1, and the optional [library_index] length is unprinted (8-nt placeholder). Pre- and final-PCR primer landing is not yet modelled as separate PCR geometry."
ATAC_INDEX_1="CATACGGCAGAAGACGAACACCTCCAACTAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
RNA_INDEX_1="CAGACGTGTGCTCTTCCGATCTNNNNNNNNNNAGTGATTAGCATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN"
TN5_ME_DDC="[Phos]CTGTCTCTTATACACATCT[ddc]"
TN5_BLOCK="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG-[invdT]"
NEXTERA_READ_1="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
PARTIAL_P5_PTO="AATGATACGGCGACCACCGAG*A"
TRUSEQ_R2="CAGACGTGTGCTCTTCCGATCT"
ATAC_SPACER="CATACGGCAGAAGACGAAC"
P7_I7_SPACER="CAAGCAGAAGACGGCATACGAGAT[library_index]CATACGGCAGAAGACGAAC"
PARTIAL_P7_PTO="CAAGCAGAAGACGGCATACGA*G"
P7_I7_TRUSEQ_R2="CAAGCAGAAGACGGCATACGAGAT[library_index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
def oligos():
    return [f"<b>{name}</b> (Supplementary Table 1, 5′→3′): <code>{sequence}</code>"
            for name, sequence in (("ATAC sample-index 1", ATAC_INDEX_1),
                                   ("RNA sample-index 1", RNA_INDEX_1),
                                   ("Tn5 mosaic-end loading oligo", TN5_ME_DDC),
                                   ("Tn5 barcode-hopping blocker", TN5_BLOCK),
                                   ("Nextera Read 1", NEXTERA_READ_1),
                                   ("Partial-P5-PTO", PARTIAL_P5_PTO),
                                   ("TruseqR2", TRUSEQ_R2),
                                   ("ATAC-spacer", ATAC_SPACER),
                                   ("P7-i7-spacer", P7_I7_SPACER),
                                   ("Partial-P7-PTO", PARTIAL_P7_PTO),
                                   ("P7-i7-TruSeqR2", P7_I7_TRUSEQ_R2))]
# --- Final libraries, assembled from the Supplementary Table 1 oligos ---------
# Source scope (PMC12165849 Methods lines 400-427; MOESM3 Table 1b-d): the ATAC
# pre-index oligo is spacer-index-Read2N, the RNA RT primer is TruSeqR2 partial-
# UMI10-index11-T30VN, hybrid tagmentation adds Nextera Read 1, pre-amplification
# uses Partial-P5-PTO with TruseqR2/ATAC-spacer, and the final PCR adds P7 plus an
# optional [library_index].  The 10x scATAC gel-bead oligo (P5 remainder plus the
# 16-nt droplet barcode in front of the Read 1 site) is not printed in the article
# or Table 1, so those regions are inferred.  The [library_index] length is not
# printed either; an 8-nt placeholder is inferred.
SAMPLE_INDEX_NT=11
DROPLET_BARCODE_NT=16
UMI_NT=10
LIBRARY_INDEX_NT=8   # inferred placeholder; Table 1b prints only [library_index]
CELL_PARTS=(("sample index",SAMPLE_INDEX_NT),("droplet barcode",DROPLET_BARCODE_NT))
_P7=P7_I7_SPACER.split("[library_index]")[0]
_READ2N=ATAC_INDEX_1[len(ATAC_SPACER)+SAMPLE_INDEX_NT:]
_TRUSEQ_R2_FULL=P7_I7_TRUSEQ_R2.split("[library_index]")[1]
_POLY_T=RNA_INDEX_1[len(TRUSEQ_R2)+UMI_NT+SAMPLE_INDEX_NT:-2]
if not (ATAC_INDEX_1.startswith(ATAC_SPACER) and _READ2N==TN5_BLOCK.split("-[")[0]
        and _TRUSEQ_R2_FULL.endswith(TRUSEQ_R2) and set(_POLY_T)=={"T"}
        and PARTIAL_P5_PTO.replace("*","")==il.P5[:len(PARTIAL_P5_PTO)-1]):
    raise AssertionError("SUM-seq Table 1 oligo parts do not assemble as described")

def _cell_part(part,length,mark):
    key=part.replace(" ","_")
    return seg(part,mark*length,"cbc",placeholder=True,
               inferred=part=="droplet barcode",
               feature=feature(f"cell_barcode_{key}","cell_barcode","combinatorial",
                               group="cell_barcode",part=part))

def _library_index():
    return seg("library index (reverse complement)","I"*LIBRARY_INDEX_NT,"cbc",
               placeholder=True,inferred=True,
               feature=feature("library_index_i7","sample_index","unknown"))

def _bead_end():
    """Inferred 10x gel-bead contribution upstream of the Nextera Read 1 site."""
    return [seg("P5 (Partial-P5-PTO; remainder from 10x bead, inferred)",il.P5,"p5",inferred=True),
            _cell_part("droplet barcode",DROPLET_BARCODE_NT,"B"),
            seg("Nextera Read 1 (R1N)",NEXTERA_READ_1,"s5")]

ATAC=Construct([*_bead_end(),
    seg("accessible genomic DNA","X"*34,placeholder=True),
    seg("Read2N reverse complement",revcomp(_READ2N),"s7"),
    _cell_part("sample index",SAMPLE_INDEX_NT,"S"),
    seg("ATAC-spacer reverse complement",revcomp(ATAC_SPACER),"r2"),
    _library_index(),
    seg("P7 reverse complement",revcomp(_P7),"p7")],name="SUM-seq ATAC library")
RNA=Construct([*_bead_end(),
    seg("cDNA insert (sense)","X"*34,placeholder=True),
    seg("poly(A) (oligo-dT complement)",revcomp(_POLY_T)),
    _cell_part("sample index",SAMPLE_INDEX_NT,"S"),
    seg("UMI","U"*UMI_NT,"umi",placeholder=True,feature=feature("umi","umi","random")),
    seg("TruSeq Read 2 reverse complement",revcomp(_TRUSEQ_R2_FULL),"r2"),
    _library_index(),
    seg("P7 reverse complement",revcomp(_P7),"p7")],name="SUM-seq RNA library")
AP=(sp.NEXTERA["R1"],sp.NEXTERA["I1"],sp.NEXTERA["I2"],sp.NEXTERA["R2"])
RP=(sp.NEXTERA["R1"],sp.TRUSEQ["I1"],sp.NEXTERA["I2"],sp.TRUSEQ["R2"])
for _lib,_primers in ((ATAC,AP),(RNA,RP)):
    _problems=sp.verify(_lib,_primers,required_roles=tuple(p.role for p in _primers))
    if _problems:
        raise ValueError(f"{_lib.name}: "+"; ".join(_problems))
# Methods 'Sequencing' (line 427): NovaSeq6000/NextSeq2000 run lengths.  The Aviti
# runs (T cell, CRISPRi/a screen) used ATAC 79/11/16/79 and the same RNA lengths.
READ_LENGTHS={"ATAC":{"Read 1":55,"Index 1 (i7)":11,"Index 2 (i5)":16,"Read 2":55},
              "RNA":{"Read 1":95,"Index 1 (i7)":6,"Index 2 (i5)":16,"Read 2":21}}
READ_USES={"ATAC":{"Read 1":"include","Index 1 (i7)":"extract_and_demultiplex",
                   "Index 2 (i5)":"extract","Read 2":"include"},
           "RNA":{"Read 1":"include","Index 1 (i7)":"demultiplex",
                  "Index 2 (i5)":"extract","Read 2":"extract"}}
FINAL_LIBRARIES=(("RNA library",RNA,RP,
    "P5 and droplet barcode from the 10x bead (inferred), hybrid-tagmentation Read 1N, cDNA, RT-primer sample index and UMI, TruSeq Read 2 and the optional P7 library index.",
    "Read 2 (21 cycles) covers the UMI and RNA sample index; Index 2 (16) reads the droplet barcode; Index 1 (6) reads the optional library index."),
    ("ATAC library",ATAC,AP,
    "P5 and droplet barcode from the 10x bead (inferred), Nextera Read 1N, genomic insert, Read2N, ATAC sample index, ATAC-spacer and the optional P7 library index.",
    "Index 1 (11 cycles) from the Nextera Index 1 site reads the ATAC sample index; Index 2 (16) reads the droplet barcode. The optional library index lies beyond the spacer and is not reached by those 11 cycles."))

def _bases(written):
    """Ordering bases of a Table 1 oligo with its terminal/backbone modifications removed."""
    if "[library_index]" in written:
        raise ValueError("an oligo with an unprinted library index has no exact sequence")
    return re.sub(r"\[(Phos|ddc|invdT)\]|-|\*","",written)

def protocol_graph():
    source="Defining source: https://doi.org/10.1038/s41592-025-02700-8"
    m="PMC12165849 Methods"
    graph=ProtocolGraph("SUM-seq",authority=source)
    nuclei=graph.input("nuclei","glyoxal-fixed nuclei, one aliquot per sample (~40,000)",
                       MaterialKind.BIOLOGICAL,primary=True)
    transposome=graph.input("atac_transposome",
        "sample-indexed Tn5: Nextera Read1 + ATAC index oligo annealed 1:1:2 with [Phos]ME-ddC, loaded on in-house Tn5",
        MaterialKind.ADAPTER,state=(NEXTERA_READ_1,ATAC_INDEX_1,_bases(TN5_ME_DDC)))
    tn5=graph.input("tn5","in-house Tn5 transposase (Hennig et al. 2018 protocol)",MaterialKind.CATALYST)
    tagmented=graph.product("tagmented_nuclei","nuclei with sample-indexed ATAC fragments",
                            MaterialKind.MOLECULAR)
    graph.react("atac_tagmentation","Tagment chromatin with the barcoded Tn5 (30 °C, 30 min)",
                StepKind.TAGMENTATION,((nuclei,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                                          (tn5,InputRole.CATALYST)),
                (tagmented,),source=f"{m}, Transposome generation and DNA transposition")
    stop=graph.input("stop_buffer","2× stop buffer (20 mM EDTA, 2% BSA)",MaterialKind.CONDITION)
    stopped=graph.product("stopped_nuclei","stopped, washed transposed nuclei",MaterialKind.MOLECULAR)
    graph.react("stop_atac","Stop tagmentation and wash three times",StepKind.QUENCHING,
                ((tagmented,InputRole.TEMPLATE),(stop,InputRole.CONDITION)),(stopped,),
                changes_molecule=False,source=f"{m}, DNA transposition/In situ RT")
    rt_primer=graph.input("indexed_rt_primer","sample-indexed oligo(dT) RT primer (TruSeqR2–UMI10–index11–T30VN)",
                          MaterialKind.OLIGO,state=RNA_INDEX_1)
    rt=graph.input("reverse_transcriptase","Maxima H Minus RT (12% PEG 8000 in most experiments)",
                   MaterialKind.CATALYST)
    rt_nuclei=graph.product("rt_nuclei","nuclei with sample-indexed RNA–cDNA hybrids and ATAC fragments",
                            MaterialKind.MOLECULAR)
    graph.react("indexed_rt","In situ RT with the sample-indexed oligo(dT), cycled 8–50 °C",
                StepKind.REVERSE_TRANSCRIPTION,
                ((stopped,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),
                (rt_nuclei,),source=f"{m}, In situ RT")
    others=graph.input("other_sample_nuclei","other samples' nuclei, processed identically with different ATAC/RNA index oligos",
                       MaterialKind.MOLECULAR)
    pooled=graph.product("pooled_nuclei","pooled sample-indexed nuclei",MaterialKind.MOLECULAR)
    graph.react("pool_samples","Pool all samples and wash",StepKind.POOL,
                ((rt_nuclei,InputRole.TEMPLATE),(others,InputRole.TEMPLATE)),(pooled,),changes_molecule=False,
                source=f"{m}, In situ RT")
    r1_tn5=graph.input("rna_transposome","Read1N-only Tn5 (Nextera Read1 annealed 1:1 with ME oligo)",
                       MaterialKind.ADAPTER,state=(NEXTERA_READ_1,_bases(TN5_ME_DDC)))
    hybrid=graph.product("hybrid_tagmented","nuclei with Read1N-tagged RNA–cDNA hybrids",
                         MaterialKind.MOLECULAR)
    graph.react("rna_tagmentation","Tagment cDNA–mRNA hybrids in the pooled nuclei (37 °C, 30 min)",
                StepKind.TAGMENTATION,((pooled,InputRole.TEMPLATE),(r1_tn5,InputRole.STRUCTURAL),
                                          (tn5,InputRole.CATALYST)),
                (hybrid,),source=f"{m}, cDNA/mRNA hybrid tagmentation")
    gap_rt=graph.input("gap_fill_rt","Maxima H Minus RT (gap fill)",MaterialKind.CATALYST)
    exo=graph.input("exonuclease_i","Thermolabile Exonuclease I",MaterialKind.CATALYST)
    filled=graph.product("gap_filled","gap-filled nuclei; free single-stranded oligos digested",
                         MaterialKind.MOLECULAR)
    graph.react("gap_fill_exo","Gap fill and ExoI digestion (37 °C, 15 min), wash and strain",
                StepKind.EXTENSION,((hybrid,InputRole.TEMPLATE),(gap_rt,InputRole.CATALYST),
                                    (exo,InputRole.CATALYST)),(filled,),source=f"{m}, Gap fill and ExoI")
    bead=graph.input("droplet_bead_oligo","10x scATAC gel-bead oligo (P5–droplet barcode–Read1N site; not printed)",
                     MaterialKind.ADAPTER)
    enzyme=graph.input("barcoding_enzyme","10x barcoding reagent B / reducing reagent B / barcoding enzyme",
                       MaterialKind.CATALYST)
    blocker=graph.input("hopping_blocker","Blocking invdT oligo (~500 nM; omitted for macrophages)",
                        MaterialKind.OLIGO,state=_bases(TN5_BLOCK))
    droplet=graph.product("droplet_products","sample × droplet barcoded RNA and ATAC fragments",
                          MaterialKind.MOLECULAR)
    graph.react("droplet_barcoding","Overload nuclei into 10x droplets; 72 °C 5 min then 4–12 linear barcoding cycles",
                StepKind.EXTENSION,((filled,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                                    (enzyme,InputRole.CATALYST),(blocker,InputRole.STRUCTURAL)),
                (droplet,),source=f"{m}, GEM generation")
    recovered=graph.product("recovered","Silane- and 1.4× SPRI-purified barcoded fragments",
                            MaterialKind.MOLECULAR)
    gem_waste=graph.product("gem_waste","partitioning oil phase and bead-unbound fraction",
                            MaterialKind.WASTE)
    graph.react("gem_cleanup","Break emulsion; Dynabeads MyOne Silane then 1.4× SPRI",StepKind.CLEANUP,
                ((droplet,InputRole.TEMPLATE),),(recovered,gem_waste),changes_molecule=False,
                source=f"{m}, GEM generation")
    graph.discard(gem_waste,"oil phase removed after recovery reagent; Silane/SPRI supernatants discarded")
    pre_primers=graph.input("preamp_primers","Partial-P5-PTO + TruseqR2 + ATAC-spacer (500 nM each)",
                            MaterialKind.OLIGO,state=(_bases(PARTIAL_P5_PTO),TRUSEQ_R2,ATAC_SPACER))
    pcr=graph.input("pcr_polymerase","NEBNext HF 2× PCR Master Mix",MaterialKind.CATALYST)
    preamp=graph.product("preamplified","pre-amplified RNA and ATAC fragments",MaterialKind.MOLECULAR)
    graph.react("preamp","Pre-amplify both modalities (5 cycles), 1.2× SPRI",StepKind.PCR,
                ((recovered,InputRole.TEMPLATE),(pre_primers,InputRole.STRUCTURAL),(pcr,InputRole.CATALYST)),
                (preamp,),source=f"{m}, Pre- and final library amplification",
                geometry=PCRGeometry.unavailable(RNA,reason="Three-primer pre-amplification of mixed RNA/ATAC fragments is not modelled as separate landing objects.",symmetric_species=True))
    rna_precursor=graph.product("rna_precursor","RNA fraction",MaterialKind.MOLECULAR)
    atac_precursor=graph.product("atac_precursor","ATAC fraction",MaterialKind.MOLECULAR)
    graph.react("library_split","Split into two equal fractions",StepKind.SPLIT,
                ((preamp,InputRole.TEMPLATE),),(rna_precursor,atac_precursor),changes_molecule=False,
                source=f"{m}, Pre- and final library amplification")
    rna_primers=graph.input("rna_pcr_primers","Partial-P5-PTO + P7_i7_TruSeqR2 (only the P5 primer has an exact printed sequence)",
                             MaterialKind.OLIGO,state=(_bases(PARTIAL_P5_PTO),))
    atac_primers=graph.input("atac_pcr_primers","Partial-P5-PTO + P7-i7-spacer (only the P5 primer has an exact printed sequence)",
                              MaterialKind.OLIGO,state=(_bases(PARTIAL_P5_PTO),))
    rna_amp=graph.product("rna_amplified","amplified RNA library",MaterialKind.MOLECULAR)
    atac_amp=graph.product("atac_amplified","amplified ATAC library",MaterialKind.MOLECULAR)
    graph.react("rna_pcr","RNA PCR (5 cycles + qPCR-determined cycles)",StepKind.PCR,
                ((rna_precursor,InputRole.TEMPLATE),(rna_primers,InputRole.STRUCTURAL),(pcr,InputRole.CATALYST)),
                (rna_amp,),source=f"{m}, Pre- and final library amplification",
                geometry=PCRGeometry.unavailable(RNA,reason="Final-PCR primer landing is not yet modelled as separate geometry; the final construct is assembled from the Table 1 oligos.",symmetric_species=True))
    graph.react("atac_pcr","ATAC PCR (5 cycles + qPCR-determined cycles)",StepKind.PCR,
                ((atac_precursor,InputRole.TEMPLATE),(atac_primers,InputRole.STRUCTURAL),(pcr,InputRole.CATALYST)),
                (atac_amp,),source=f"{m}, Pre- and final library amplification",
                geometry=PCRGeometry.unavailable(ATAC,reason="Final-PCR primer landing is not yet modelled as separate geometry; the final construct is assembled from the Table 1 oligos.",symmetric_species=True))
    rna_library=graph.product("rna_library",RNA.name,MaterialKind.LIBRARY,state=RNA)
    atac_library=graph.product("atac_library",ATAC.name,MaterialKind.LIBRARY,state=ATAC)
    for ident,amp,library,ratio in (("rna_cleanup",rna_amp,rna_library,"0.72×"),
                                    ("atac_cleanup",atac_amp,atac_library,"1.2×")):
        waste=graph.product(f"{ident}_waste",f"{ratio} SPRI supernatant",MaterialKind.WASTE)
        graph.react(ident,f"{ratio} SPRI cleanup",StepKind.CLEANUP,((amp,InputRole.TEMPLATE),),
                    (library,waste),changes_molecule=False,
                    source=f"{m}, Pre- and final library amplification")
        graph.discard(waste,f"{ratio} SPRI supernatant discarded")
    for material,library,primers,key in ((rna_library,RNA,RP,"RNA"),(atac_library,ATAC,AP,"ATAC")):
        graph.sequence(material,SequencingPlan(
            library=library,primers=primers,acquired_roles=tuple(p.role for p in primers),
            read_lengths=READ_LENGTHS[key],read_uses=READ_USES[key]))
    graph.annotate(Highlight("droplet_barcoding","Barcode-hopping mitigation",
        "Residual barcoded ATAC pre-index oligos released in multinucleated droplets can prime on ATAC fragments. The 3′-inverted-dT blocking oligo (Read2N) competes for those sites and cannot be extended; fewer linear cycles further reduce hopping. The blocker was omitted in the macrophage experiment."))
    graph.annotate(Highlight("atac_cleanup","Run lengths by platform",
        "NovaSeq6000/NextSeq2000 (benchmark, macrophage): ATAC 55/11/16/55, RNA 95/6/16/21 (Read 1/Index 1/Index 2/Read 2). Aviti (T cell, CRISPRi/a screen): ATAC 79/11/16/79, RNA unchanged. The page uses the Illumina lengths. The ATAC optional library index lies 19 nt beyond the 11 acquired Index 1 cycles; how it is read, if at all, is not stated in the inspected Methods."))
    graph.declare_reviewed_scope((
        "nuclei","atac_transposome","tn5","atac_tagmentation","stop_buffer","stop_atac",
        "indexed_rt_primer","reverse_transcriptase","indexed_rt","other_sample_nuclei","pool_samples",
        "rna_transposome","rna_tagmentation","gap_fill_rt","exonuclease_i","gap_fill_exo",
        "droplet_bead_oligo","barcoding_enzyme","hopping_blocker","droplet_barcoding",
        "gem_cleanup","preamp_primers","pcr_polymerase","preamp","library_split",
        "rna_pcr_primers","atac_pcr_primers","rna_pcr","atac_pcr","rna_cleanup","atac_cleanup",
        "rna_library","atac_library"),source)
    return graph.require_publishable()
