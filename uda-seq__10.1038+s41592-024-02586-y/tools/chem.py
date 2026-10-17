"""UDA-seq universal droplet-first combinatorial indexing.

The graph follows the representative Multiome (sn3'-RNA + ATAC) UDA-seq workflow of the
published Methods (Wang et al., Nature Methods 2025, PMC12165859): "i7-only Tn5
transposome assembly", "Multiome (sn3'-RNA and ATAC-seq) UDA-seq" and its subsections,
and the shared "sc-RNA-seq library construction".  Oligos are transcribed from
Supplementary Table 7 (MOESM8, sheet "primer summary").  10x Genomics kit sequences
(gel-bead oligos, ATAC enzyme B adapters) are not disclosed by the paper and are
marked inferred.
"""
import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (AnnotationSlot, Disposition, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, SourceClaim, SourceDiscrepancy,
                            StepKind, Terminal)

TITLE="UDA-seq — droplet-first combinatorial indexing"
NOTES="01_uda-seq.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41592-024-02586-y">Wang et al., <i>Nature Methods</i> (2025)</a>.'
SUMMARY=("Fixed cells or nuclei are intentionally overloaded into 10x droplets for the first "
         "(droplet) barcode, released intact, redistributed into a 96-well plate and given a "
         "second, well-specific index by PCR. The barcode pair resolves cells for several "
         "droplet-compatible modalities; this page draws the published Multiome RNA + ATAC route.")
CAVEAT=("UDA-seq wraps several 10x assay chemistries (3′/5′ RNA, VDJ, CRISPR, Multiome). This "
        "graph draws the Multiome route. Supplementary Table 7 prints the UDA index primers and "
        "the RNA-route i7-only Tn5 cassette. The 10x gel-bead oligos and ATAC enzyme B adapters "
        "are not disclosed by the paper; their regions are drawn as inferred. Sequencing was on "
        "MGISeq-T7; no library conversion or read primers are described in the inspected text, "
        "so read contents are derived from the drawn construct rather than primer landings.")

AUTH="https://doi.org/10.1038/s41592-024-02586-y"
KIT=AUTH+" (10x Genomics Multiome kit sequence not disclosed; inferred region)"

# Supplementary Table 7, sheet "primer summary".
P5_PRIMER=il.P5[:22]                   # "P5 primer" AATGATACGGCGACCACCGAGA
P7_PRIMER=il.P7                        # "P7 primer"
PARTIAL_TSO=rt.SMART_HANDLE[:-1]       # "Partial TSO" AAGCAGTGGTATCAACGCAGAG
I5_P5_PART=il.P5[:25]                  # 5' part of every Truseq-i5 index primer, as printed


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def well_feature(part_name):
    return feature("cell_well", "cell_barcode", "combinatorial", group="cell_id", part="well")


def droplet_feature():
    return feature("cell_droplet", "cell_barcode", "whitelist",
                   whitelist="10x Multiome gel-bead whitelist", group="cell_id", part="droplet")


def truseq_i5_index_primer():
    """Truseq-i5 index primer (96 versions; well-unique); 3' C-s-T phosphorothioate."""
    return [seg("P5 (as printed, no ACAC)", I5_P5_PART, "p5"),
            seg("well index (round 2)", "N"*8, "cbc", placeholder=True, feature=well_feature("i5")),
            seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1")]


def nextera_i7_index_primer(role="well"):
    feat = (well_feature("i7") if role == "well" else
            feature("sample_i7", "sample_index", "whitelist",
                    whitelist="Supplementary Table 7 Nextera-i7 index set"))
    return [seg("P7", il.P7, "p7"),
            seg("well index (round 2)" if role == "well" else "sample i7 index", "N"*8, "cbc",
                placeholder=True, feature=feat),
            seg("Nextera s7", nx.S7, "s7")]


def tn5_top():
    return [seg("Tn5-top_ME [5Phos]…[3ddC]", nx.ME_RC, "me")]


def tn5_bottom():
    return [seg("Read 2N (s7)", nx.S7, "s7"), seg("mosaic end", nx.ME, "me")]


def oligos():
    def show(name, parts, mods=""):
        text = "".join(p.top for p in parts)
        return f"<b>{name}</b> (Supplementary Table 7, 5′→3′): <code>{mods}{text}</code>"
    return [show("Truseq-i5 index primer (C-s-T 3′ phosphorothioate)", truseq_i5_index_primer()),
            show("Nextare-i7 index primer", nextera_i7_index_primer()),
            show("P5 primer", [seg("", P5_PRIMER)]),
            show("P7 primer", [seg("", P7_PRIMER)]),
            show("Partial TSO (Bio-Partial TSO is 5′-biotinylated)", [seg("", PARTIAL_TSO)]),
            f"<b>Tn5-top_ME</b> (Supplementary Table 7): <code>[5Phos]{nx.ME_RC}[3ddC]</code>",
            show("Tn5-bottom_Read2N", tn5_bottom())]


def rna_library():
    """Representative Multiome GEX UDA-seq library, P5 to P7'."""
    return Construct([
        seg("P5 (Truseq-i5 primer)", I5_P5_PART, "p5"),
        seg("well index (round 2)", "N"*8, "cbc", placeholder=True, feature=well_feature("i5")),
        seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        seg("10x droplet barcode (round 1)", "N"*16, "cbc", placeholder=True, inferred=True,
            feature=droplet_feature()),
        seg("UMI", "N"*12, "umi", placeholder=True, inferred=True,
            feature=feature("umi", "umi", "random")),
        seg("poly(dT)", "T"*30, inferred=True),
        seg("cDNA", "X"*34, placeholder=True),
        seg("mosaic end (rc)", nx.ME_RC, "me"),
        seg("s7 (rc)", nx.S7_RC, "s7"),
        seg("sample i7 index", "N"*8, "cbc", placeholder=True,
            feature=feature("sample_i7", "sample_index", "whitelist",
                            whitelist="Supplementary Table 7 Nextera-i7 index set")),
        seg("P7 (rc)", il.P7_RC, "p7"),
    ], name="UDA-seq Multiome RNA library")


def atac_library():
    """Representative Multiome ATAC UDA-seq library, P5 to P7'."""
    return Construct([
        seg("P5", il.P5, "p5", inferred=True),
        seg("10x droplet barcode (round 1)", "N"*16, "cbc", placeholder=True, inferred=True,
            feature=droplet_feature()),
        seg("s5", nx.S5, "s5", inferred=True),
        seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X"*34, placeholder=True),
        seg("mosaic end (rc)", nx.ME_RC, "me"),
        seg("s7 (rc)", nx.S7_RC, "s7", inferred=True),
        seg("well index (round 2)", "N"*8, "cbc", placeholder=True, feature=well_feature("i7")),
        seg("P7 (rc)", il.P7_RC, "p7"),
    ], name="UDA-seq Multiome ATAC library")


RNA=rna_library()
ATAC=atac_library()
SEQ_LIMIT=("Sequenced on MGISeq-T7 (RNA PE150, ATAC PE100). The inspected Methods describe no "
           "Illumina-to-MGI library conversion and no sequencing primers, and index-read lengths "
           "are not stated; read contents are derived from the drawn construct.")
FINAL_LIBRARIES=(
    ("Representative Multiome RNA library",RNA,(),
     "INFERRED — 10x bead-derived droplet barcode, UMI and poly(dT) regions; well index from the Truseq-i5 primer, sample i7 from the Nextera-i7 primer.",
     "RNA branch"),
    ("Representative Multiome ATAC library",ATAC,(),
     "INFERRED — 10x bead and ATAC-adapter regions; the well index comes from the Nextera-i7 index primer.",
     "ATAC branch"),
)
READ_SPECS={
    "Representative Multiome RNA library":(
        sp.read_spec("Read 1",150,"10x droplet barcode (16) + UMI (12), then poly(dT)/cDNA (derived from construct)","extract"),
        sp.read_spec("Index 2 (i5)",None,"8-bp well index (round-2 'post-index'; FASTQs grouped by it)","extract_and_demultiplex"),
        sp.read_spec("Index 1 (i7)",None,"8-bp library sample index (Nextera-i7 primer)","demultiplex"),
        sp.read_spec("Read 2",150,"cDNA from the Tn5 end","include")),
    "Representative Multiome ATAC library":(
        sp.read_spec("Read 1",100,"accessible genomic DNA","include"),
        sp.read_spec("Index 2 (i5)",None,"16-bp 10x droplet barcode (derived from construct)","extract"),
        sp.read_spec("Index 1 (i7)",None,"8-bp well index (round 2)","extract_and_demultiplex"),
        sp.read_spec("Read 2",100,"accessible genomic DNA","include")),
}


def _pcr(g, ident, label, template, primers, polymerase, product, final, reason):
    g.react(ident, label, StepKind.PCR,
            ((template,InputRole.TEMPLATE),*((p,InputRole.STRUCTURAL) for p in primers),
             (polymerase,InputRole.CATALYST)),(product,),
            geometry=PCRGeometry.unavailable(final, reason=reason))


def _cleanup(g, ident, label, template, kept, waste_label):
    waste=g.product(ident+"_waste",waste_label,MaterialKind.WASTE)
    g.react(ident,label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),),(kept,waste))
    g.discard(waste,"removed by the bead cleanup")


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTH)
    O=MaterialKind.OLIGO
    sample=g.input("sample","isolated nuclei",MaterialKind.BIOLOGICAL,primary=True)
    fixative=g.input("fixative","formaldehyde 1% (cell lines) or 1.5% (tissue), RT 10 min",MaterialKind.CONDITION)
    fixed=g.product("fixed","fixed nuclei",MaterialKind.MOLECULAR)
    g.react("fix","Fix nuclei and wash in PBS-BSA-RNase inhibitor",StepKind.CROSS_LINKING,
            ((sample,InputRole.TEMPLATE),(fixative,InputRole.CONDITION)),(fixed,))

    atac_tn5=g.input("atac_transposome","10x ATAC enzyme B (loaded adapters not disclosed)",
                     MaterialKind.ADAPTER,inferred=True,source=KIT)
    atac_tn5_enzyme=g.input("atac_tn5","Tn5 transposase in ATAC enzyme B",MaterialKind.CATALYST)
    tagged=g.product("tagged","fixed nuclei with tagmented chromatin",MaterialKind.MOLECULAR)
    tag_waste=g.product("tagment_waste","fragments lacking a bead-ligatable end pair",MaterialKind.WASTE)
    g.react("tagment","Tagment chromatin in nuclei before loading (7 µL ATAC buffer B + 3 µL "
            "ATAC enzyme B, 37 °C 60 min)",StepKind.TAGMENTATION,
            ((fixed,InputRole.TEMPLATE),(atac_tn5,InputRole.STRUCTURAL),
             (atac_tn5_enzyme,InputRole.CATALYST)),(tagged,tag_waste))
    g.discard(tag_waste,"not physically separated; carried along without forming the two-ended ATAC library")

    bead=g.input("gel_bead","10x Multiome gel-bead oligos (GEX and ATAC barcodes)",O,
                 inferred=True,source=KIT)
    loading=g.input("overload","overloaded Chromium Controller loading (Multiome kit PN-1000283)",
                    MaterialKind.CONDITION)
    droplets=g.product("droplets","overloaded GEMs with several nuclei per bead",MaterialKind.MOLECULAR)
    empty=g.product("droplet_waste","partitioning oil and unloaded droplets",MaterialKind.WASTE)
    g.react("partition","Overload nuclei into 10x Multiome GEMs",StepKind.PARTITION,
            ((tagged,InputRole.TEMPLATE),(loading,InputRole.CONDITION)),(droplets,empty))
    g.discard(empty,"partitioning oil is removed at droplet release")
    enzymes=g.input("gem_enzymes","10x in-GEM reverse transcriptase and ligase",MaterialKind.CATALYST)
    barcoded=g.product("droplet_barcoded","nuclei with droplet-barcoded cDNA and ATAC fragments",
                       MaterialKind.MOLECULAR,inferred=True,source=KIT)
    g.react("gem_rt_ligation","In-GEM RT and ATAC barcode ligation (37 °C 45 min, 25 °C 30 min)",
            StepKind.MAGIC,((droplets,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),
                            (enzymes,InputRole.CATALYST)),(barcoded,),
            reason="Kit-internal reactions: the 10x bead oligo, splint and TSO structures are not disclosed by the paper.")

    recovery=g.input("recovery_agent","125 µL recovery agent; oil removed; 120 µL PBS added",
                     MaterialKind.CONDITION)
    released=g.product("released","released intact droplet-barcoded nuclei",MaterialKind.MOLECULAR)
    oil=g.product("recovery_waste","recovery agent/partitioning oil phase",MaterialKind.WASTE)
    g.react("release","Break the emulsion and keep the aqueous phase with intact nuclei",
            StepKind.EXTRACTION,((barcoded,InputRole.TEMPLATE),(recovery,InputRole.CONDITION)),
            (released,oil),changes_molecule=False)
    g.discard(oil,"pink oil phase removed from the bottom of the tube and discarded")
    well=g.product("well","representative well (2 µL of nuclei)",MaterialKind.MOLECULAR)
    other=g.product("other_wells","other wells of the 96-well plate",MaterialKind.MOLECULAR)
    g.react("redistribute","Mix and dispense 2 µL per well into a 96-well plate",StepKind.SPLIT,
            ((released,InputRole.TEMPLATE),),(well,other),changes_molecule=False)
    g.mark_terminal(Terminal("other_wells",Disposition.UNRESOLVED,
        "each other well is processed identically with its own well-unique index primers; drawn once"))

    pk=g.input("proteinase_k","1 µL proteinase K + 5 µL EB (55 °C 10 min, 85 °C 5 min)",MaterialKind.CATALYST)
    digested=g.product("digested","deproteinized well contents",MaterialKind.MOLECULAR)
    g.react("deproteinize","Proteinase K digestion and heat inactivation",StepKind.DEPROTEINIZATION,
            ((well,InputRole.TEMPLATE),(pk,InputRole.CATALYST)),(digested,))
    purified=g.product("purified","SILANE-purified nucleic acids (10 µL EB)",MaterialKind.MOLECULAR)
    _cleanup(g,"silane","Dynabeads SILANE purification",digested,purified,"SILANE supernatant and washes")

    i5=g.input("truseq_i5","Truseq-i5 index primer (well-unique)",O,state=tuple(truseq_i5_index_primer()))
    i7w=g.input("nextera_i7_well","Nextare-i7 index primer (well-unique)",O,state=tuple(nextera_i7_index_primer()))
    p5=g.input("p5_primer","P5 primer",O,state=(seg("P5 primer",P5_PRIMER,"p5"),))
    tso=g.input("partial_tso","Partial TSO (Methods: Partial TSO/IS)",O,state=(seg("Partial TSO",PARTIAL_TSO,"tso"),))
    q5=g.input("q5","NEBNext Ultra II Q5 Master Mix",MaterialKind.CATALYST)
    indexed=g.product("well_indexed","well-indexed cDNA and ATAC products",MaterialKind.MOLECULAR,
                      inferred=True,source=KIT)
    g.react("well_pcr","Well-index PCR for both modalities (72 °C 5 min, 98 °C 3 min, 7 × 98 °C 20 s / "
            "63 °C 30 s / 72 °C 1 min, 72 °C 1 min)",StepKind.PCR,
            ((purified,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),(tso,InputRole.STRUCTURAL),
             (p5,InputRole.STRUCTURAL),(i7w,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(indexed,),
            geometry=PCRGeometry.unavailable(RNA,reason="Four primers act on two kit-defined template classes whose 10x end sequences are not disclosed."))
    pooled=g.product("pooled","half of each well's product pooled and purified (1.6× XP, 160 µL EB)",MaterialKind.MOLECULAR)
    unpooled=g.product("unpooled_half","the other half of each well's PCR product",MaterialKind.MOLECULAR)
    xp=g.input("xp_beads","AMPure XP beads",MaterialKind.SELECTION)
    pool_waste=g.product("pool_waste","1.6× XP supernatant",MaterialKind.WASTE)
    g.react("pool_purify","Pool half of the plate's indexed products and purify (1.6× XP)",
            StepKind.SELECTION,((indexed,InputRole.TEMPLATE),(xp,InputRole.SELECTION)),
            (pooled,unpooled,pool_waste))
    g.mark_terminal(Terminal("unpooled_half",Disposition.STORED,"not used in the described library construction"))
    g.discard(pool_waste,"removed by the bead cleanup")
    atac_in=g.product("atac_input","40 µL post-indexing product for ATAC",MaterialKind.MOLECULAR)
    rna_in=g.product("rna_input","40 µL post-indexing product for cDNA",MaterialKind.MOLECULAR)
    g.react("split_modalities","Take 40 µL each for ATAC and cDNA enrichment",StepKind.SPLIT,
            ((pooled,InputRole.TEMPLATE),),(atac_in,rna_in),changes_molecule=False)

    kapa=g.input("kapa","KAPA HiFi HotStart ReadyMix",MaterialKind.CATALYST)
    p7=g.input("p7_primer","P7 primer",O,state=(seg("P7 primer",P7_PRIMER,"p7"),))
    atac_amp=g.product("atac_amplified","amplified ATAC library",MaterialKind.MOLECULAR,inferred=True,source=KIT)
    _pcr(g,"atac_pcr","ATAC enrichment with Partial P5 and P7 primer (98 °C 45 s, 7–9 × 98 °C 20 s / "
         "67 °C 30 s / 72 °C 20 s, 72 °C 1 min)",atac_in,(p5,p7),kapa,atac_amp,ATAC,
         "Methods name a 'Partial P5' primer absent from Supplementary Table 7; the bead-derived P5 end is not disclosed.")
    atac_lib=g.product("atac_library",ATAC.name,MaterialKind.LIBRARY,state=ATAC,inferred=True,source=KIT)
    atac_size_waste=g.product("atac_size_waste","fragments outside 0.5–1.4× XP",MaterialKind.WASTE)
    g.react("atac_size","Double-sided size selection (0.5–1.4× XP, 25.5 µL EB)",StepKind.SELECTION,
            ((atac_amp,InputRole.TEMPLATE),(xp,InputRole.SELECTION)),(atac_lib,atac_size_waste))
    g.discard(atac_size_waste,"removed by double-sided bead size selection")
    g.sequence(atac_lib,SequencingPlan(library=ATAC,limitation=SEQ_LIMIT))

    bio_tso=g.input("bio_partial_tso","Bio-Partial TSO (Methods: Bio-Tso/ISPCR)",O,
                    state=(seg("Bio-Partial TSO (5′ biotin)",PARTIAL_TSO,"tso"),))
    biotin_cdna=g.product("biotin_cdna","5′-biotinylated amplified cDNA",MaterialKind.MOLECULAR,inferred=True,source=KIT)
    _pcr(g,"cdna_bio_pcr","cDNA enrichment with Partial P5 and Bio-Tso/ISPCR (98 °C 3 min, 4–5 × "
         "98 °C 15 s / 63 °C 20 s / 72 °C 1 min)",rna_in,(p5,bio_tso),kapa,biotin_cdna,RNA,
         "The TSO-side 10x template end and the Methods 'Partial P5' identity are not disclosed.")
    c1=g.input("streptavidin_c1","Dynabeads MyOne Streptavidin C1 (B&W buffer, RT 60 min)",MaterialKind.SELECTION)
    captured=g.product("captured_cdna","streptavidin-captured cDNA",MaterialKind.MOLECULAR,inferred=True,source=KIT)
    unbound=g.product("unbound","non-biotinylated products (ATAC fragments, primers)",MaterialKind.WASTE)
    g.react("capture","Capture biotinylated cDNA on streptavidin beads",StepKind.AFFINITY_CAPTURE,
            ((biotin_cdna,InputRole.TEMPLATE),(c1,InputRole.SELECTION)),(captured,unbound))
    g.discard(unbound,"removed with supernatant and B&W/EB washes")
    cdna_amp=g.product("cdna_amplified","on-bead amplified cDNA",MaterialKind.MOLECULAR,inferred=True,source=KIT)
    _pcr(g,"cdna_pcr","On-bead cDNA amplification with Partial P5 and Tso/ISPCR (4–5 cycles)",
         captured,(p5,tso),kapa,cdna_amp,RNA,"As for the biotinylated enrichment PCR.")
    cdna=g.product("cdna","purified cDNA (0.6× XP, 40 µL EB)",MaterialKind.MOLECULAR,inferred=True,source=KIT)
    _cleanup(g,"cdna_cleanup","0.6× XP cleanup",cdna_amp,cdna,"0.6× XP supernatant")

    top=g.input("tn5_top","Tn5-top_ME [5Phos]…[3ddC]",O,state=tuple(tn5_top()))
    bottom=g.input("tn5_bottom","Tn5-bottom_Read2N (Methods: TN5_R2_index)",O,state=tuple(tn5_bottom()))
    cassette=g.product("tn5_cassette","annealed i7-only oligonucleotide cassette",MaterialKind.ADAPTER)
    g.react("anneal_cassette","Anneal 1:1 at 10 µM (75 °C 15 min, 60/50/40 °C 10 min each, 25 °C 30 min)",
            StepKind.HYBRIDIZATION,((bottom,InputRole.TEMPLATE),(top,InputRole.STRUCTURAL)),(cassette,))
    truprep=g.input("trueprep_tn5","TruePrep Tagment Enzyme (Vazyme S601-01) + coupling buffer",MaterialKind.CATALYST)
    i7_tn5=g.product("i7_only_tn5","assembled i7-only Tn5 transposome",MaterialKind.ADAPTER)
    g.react("assemble_tn5","Couple 14 µL cassette with 8 µL Tn5 in 78 µL coupling buffer (30 °C 1 h)",
            StepKind.GENERIC,((cassette,InputRole.TEMPLATE),(truprep,InputRole.CATALYST)),(i7_tn5,),
            changes_molecule=False,reason="Non-covalent loading of the annealed cassette onto Tn5.")
    frag=g.product("tagmented_cdna","i7-only-tagmented cDNA (0.8× XP purified)",MaterialKind.MOLECULAR,inferred=True,source=KIT)
    frag_waste=g.product("cdna_tag_waste","internal fragments with s7 ends only",MaterialKind.WASTE)
    g.react("tagment_cdna","Tagment 50 ng cDNA with i7-only Tn5 (55 °C 10 min; 0.8× XP)",
            StepKind.TAGMENTATION,((cdna,InputRole.TEMPLATE),(i7_tn5,InputRole.STRUCTURAL),
                                   (truprep,InputRole.CATALYST)),(frag,frag_waste))
    g.discard(frag_waste,"not physically separated; lacks the P5-side end needed by the library PCR")
    i7s=g.input("nextera_i7_sample","Nextare-i7 index primer (library sample index)",O,
                state=tuple(nextera_i7_index_primer("sample")))
    rna_amp=g.product("rna_amplified","amplified RNA library",MaterialKind.MOLECULAR,inferred=True,source=KIT)
    _pcr(g,"rna_pcr","Library PCR with Partial P5 and Nextare-i7 index primer (72 °C 5 min, 98 °C 45 s, "
         "8–9 × 98 °C 20 s / 60 °C 30 s / 72 °C 1 min, 72 °C 5 min)",frag,(p5,i7s),q5,rna_amp,RNA,
         "Methods name a 'Partial P5' primer absent from Supplementary Table 7.")
    rna_lib=g.product("rna_library",RNA.name,MaterialKind.LIBRARY,state=RNA,inferred=True,source=KIT)
    rna_size_waste=g.product("rna_size_waste","fragments outside 0.6–0.8× XP",MaterialKind.WASTE)
    g.react("rna_size","Double-sided size selection (0.6–0.8× XP, 25.5 µL EB)",StepKind.SELECTION,
            ((rna_amp,InputRole.TEMPLATE),(xp,InputRole.SELECTION)),(rna_lib,rna_size_waste))
    g.discard(rna_size_waste,"removed by double-sided bead size selection")
    g.sequence(rna_lib,SequencingPlan(library=RNA,limitation=SEQ_LIMIT))

    g.annotate(SourceDiscrepancy("rna_pcr","Primer names differ between Methods and Supplementary Table 7",(
        SourceClaim("Methods, Multiome sc-ATAC/cDNA enrichment and sc-RNA library construction",
                    "'Partial P5', 'Partial TSO/IS', 'Bio-Tso/ISPCR', 'Tso/ISPCR' and 'TN5_R2_index'"),
        SourceClaim("Supplementary Table 7, primer summary",
                    "only 'P5 primer' AATGATACGGCGACCACCGAGA, 'Partial TSO', 'Bio-Partial TSO' and 'Tn5-bottom_Read2N'")),
        model_choice="The Table 7 P5 primer stands for 'Partial P5', Partial TSO for TSO/IS primers, and Tn5-bottom_Read2N for TN5_R2_index.",
        evidence="Table 7 is the only printed oligo list and the Methods cite it for every primer.",
        uncertainty="Whether 'Partial P5' is a distinct, unprinted oligo is not resolved from the inspected sources."))
    g.annotate(Highlight("truseq_i5","Truseq-i5 primer as printed",
        f"Every Truseq-i5 index primer joins {I5_P5_PART} directly to the 8-bp index, without the "
        "ACAC that precedes the i5 index in standard dual-index P5 adapters; the model keeps the printed sequence.",
        slot=AnnotationSlot.SEQUENCING))

    reviewed=[k for k in g.steps]+[k for k in g.materials if k not in g._producer]+["atac_library","rna_library"]
    g.declare_reviewed_scope(tuple(dict.fromkeys(reviewed)),AUTH)
    return g.require_publishable()
