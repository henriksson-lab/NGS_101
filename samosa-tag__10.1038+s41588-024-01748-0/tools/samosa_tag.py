"""SAMOSA-Tag molecular model."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from batch_ngs import joined_scene, seg
from chemdraw import Row, Scene, Segment, feature, revcomp
from longread_footprinting import hairpin_tagmented_smrtbell, marked_fiber, pacbio_entry, smrtbell_rows
from topology import HairpinTopology
import nextera as nx
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "SAMOSA-Tag — hairpin tagmentation of marked chromatin"
NOTES = "01_samosa-tag.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/s41588-024-01748-0">Dennis et al., <i>Nature Genetics</i> (2024)</a>.'
SUMMARY = "EcoGII first records chromatin accessibility as m6dA; hairpin-loaded Tn5 then fragments the marked fibers while installing PacBio-compatible closed ends."
CAVEAT = "Supplementary Table 5 publishes the complete phosphorylated hairpin-Tn5 oligos. A representative barcode member is drawn base by base and folded by its validated 27-bp self-complementary stem. Only the separate PacBio sequencing-primer sequence remains proprietary."
BARCODE_001 = "CGGAAGAA"
HAIRPIN_LOOP = "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT"
HAIRPIN_STEM_NT = len(nx.ME) + len(BARCODE_001)
HAIRPIN_ADAPTER = (revcomp(nx.ME) + revcomp(BARCODE_001) + HAIRPIN_LOOP
                   + BARCODE_001 + nx.ME)
SMRTBELL = hairpin_tagmented_smrtbell(
    "m6dA-marked chromatin fragment", adapter_sequence=HAIRPIN_ADAPTER, stem_length=HAIRPIN_STEM_NT)
FINAL_LIBRARY = None
SEQ_PRIMERS = ()
SEQUENCING_ENDING = pacbio_entry(SMRTBELL, adapter_sequence_known=True)
INITIAL_ROWS=tuple(marked_fiber("EcoGII", target="accessible adenines"))
INITIAL_NAME="Native chromatin fiber"
TOPOLOGY_PANELS = ((
    "Published hairpin-Tn5 oligo",
    HairpinTopology("SMRT-A_bc001", HAIRPIN_ADAPTER, 0,
                    len(HAIRPIN_ADAPTER) - HAIRPIN_STEM_NT,
                    HAIRPIN_STEM_NT).folded_panel(
        "The 5′-phosphorylated 99-nt oligo folds through its reverse-complementary "
        "barcode + mosaic-end arms; Supplementary Table 5 publishes 26 barcode members.")),)


def oligos():
    return [
        f"SMRT-A_bc001 — 5′-/5Phos/{HAIRPIN_ADAPTER}-3′",
        f"sample barcode in SMRT-A_bc001 — 5′-{BARCODE_001}-3′",
        f"Tn5 mosaic end at both stem termini — 5′-{nx.ME}-3′",
    ]


def hairpin_adapter_scene():
    sc = Scene()
    sc.strand("hairpin oligo", [
        Segment("mosaic end reverse complement", revcomp(nx.ME), "me"),
        Segment("barcode reverse complement", revcomp(BARCODE_001), "cbc",
                feature=feature("sample_barcode_001", "sample_index", "whitelist",
                                whitelist="Supplementary Table 5")),
        Segment("PacBio hairpin loop", HAIRPIN_LOOP),
        Segment("sample barcode", BARCODE_001, "cbc",
                feature=feature("sample_barcode_001", "sample_index", "whitelist",
                                whitelist="Supplementary Table 5")),
        Segment("mosaic end", nx.ME, "me"),
    ], label="SMRT-A_bc001", mod5="/5Phos/")
    sc.mark("hairpin oligo", "sample barcode", "8-nt sample barcode")
    return sc

def sections():
    return [
        ("Mark accessible adenines in nuclei", marked_fiber("EcoGII", target="accessible adenines"),
         "EcoGII and SAM deposit m6dA on exposed DNA while protein-bound intervals remain footprints."),
        ("Tagment with hairpin-loaded Tn5", joined_scene([seg("left hairpin-ME","X"*20,"me",placeholder=True),seg("marked chromatin fragment","X"*36,placeholder=True),seg("right ME-hairpin","X"*20,"me",placeholder=True)],(("left hairpin-ME","marked chromatin fragment","Tn5 transfer"),("marked chromatin fragment","right ME-hairpin","Tn5 transfer")),label="hairpin-tagmented chromatin").rows(),
         "Tn5 couples fragmentation to transfer of a hairpin-bearing mosaic-end adapter at both boundaries."),
        ("Gap-repair the tagmented product", smrtbell_rows(SMRTBELL),
         "Polymerase fills the transposition gaps, leaving one covalently closed PacBio template."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    chromatin=graph.input("chromatin","permeabilized nuclei (native chromatin fibers) in Buffer M",
                          MaterialKind.BIOLOGICAL,primary=True)
    ecogii=graph.input("ecogii","EcoGII adenine methyltransferase (250 U, 30 min, 37 °C)",
                       MaterialKind.CATALYST)
    sam=graph.input("sam","S-adenosylmethionine (1 mM; replenished to 1.16 mM at 15 min)",
                    MaterialKind.CONDITION)
    marked=graph.product("marked_fibers","m6dA-marked chromatin fibers",MaterialKind.MOLECULAR,
                         state=INITIAL_ROWS)
    graph.react("mark_accessibility","Mark accessible adenines in nuclei",
                StepKind.CHEMICAL_CONVERSION,
                ((chromatin,InputRole.TEMPLATE),(ecogii,InputRole.CATALYST),
                 (sam,InputRole.CONDITION)),(marked,))
    hairpin_adapter=graph.input("hairpin_transposome","SMRT-A barcoded hairpin oligo, annealed (95 °C 5 min, 25 °C 30 min)",
                                MaterialKind.ADAPTER,state=hairpin_adapter_scene())
    tn5=graph.input("tn5","triple-mutant Tn5 (R27S,E54K,L372P), loaded 55 min at 23 °C",
                    MaterialKind.CATALYST)
    tagmented=graph.product("hairpin_tagmented","hairpin-tagmented marked chromatin",
                            MaterialKind.MOLECULAR)
    graph.react("hairpin_tagment","In situ tagmentation: 9.4–18.8 pmol transposome, Omni-ATAC buffer, 55 °C 45–60 min",
                StepKind.TAGMENTATION,
                ((marked,InputRole.TEMPLATE),(hairpin_adapter,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagmented,))
    deprotein=graph.input("stop_mix","RNase A, then proteinase K / SDS / EDTA (60 °C, 1–2 h)",
                          MaterialKind.CATALYST)
    released=graph.product("released_dna","deproteinized hairpin-tagmented DNA",
                           MaterialKind.MOLECULAR)
    graph.react("terminate","Terminate tagmentation and deproteinize",
                StepKind.DEPROTEINIZATION,
                ((tagmented,InputRole.TEMPLATE),(deprotein,InputRole.CATALYST)),(released,))
    beads=graph.input("spri_tagmented","2× SPRI, then 0.6× SPRI (>500 bp)",MaterialKind.SELECTION)
    purified=graph.product("purified_tagmented","purified hairpin-tagmented DNA (>500 bp enriched)",
                           MaterialKind.MOLECULAR)
    small=graph.product("small_fragments","fragments below the 0.6× SPRI cutoff, protein and buffer",
                        MaterialKind.WASTE)
    graph.react("cleanup_tagmented","Bead cleanup and 0.6× SPRI size enrichment",
                StepKind.CLEANUP,
                ((released,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(purified,small))
    graph.discard(small,"left in the SPRI supernatant / washes")
    repair=graph.input("gap_repair_polymerase","Repair Mix: Phusion-HF (2 U) with 0.8 mM dNTPs",
                       MaterialKind.CATALYST)
    ligase=graph.input("gap_repair_ligase","Taq DNA ligase (80 U) in 1× Taq ligase buffer",
                       MaterialKind.CATALYST)
    repaired=graph.product("gap_repaired","gap-repair products: closed SMRTbell templates plus unrepaired molecules",
                           MaterialKind.MOLECULAR)
    graph.react("gap_repair","Fill both 9-nt transposition gaps and seal the nicks (37 °C, 1 h; 2× SPRI)",
                StepKind.EXTENSION,
                ((purified,InputRole.TEMPLATE),(repair,InputRole.CATALYST),
                 (ligase,InputRole.CATALYST)),(repaired,))
    exo=graph.input("exonuclease_iii","ExoDigest Mix: exonuclease III (100 U per 160 ng), 1× NEBuffer 2",
                    MaterialKind.CATALYST)
    digested=graph.product("exo_digested","closed templates with exonuclease-degraded open molecules",
                           MaterialKind.MOLECULAR)
    graph.react("exo_digest","Exonuclease III digestion enriches covalently closed templates (37 °C, 1 h)",
                StepKind.DIGESTION,
                ((repaired,InputRole.TEMPLATE),(exo,InputRole.CATALYST)),(digested,))
    final_beads=graph.input("spri_final","2× SPRI beads, elution in 12 µl EB",MaterialKind.SELECTION)
    library=graph.product("library",SMRTBELL.name,MaterialKind.LIBRARY,state=SMRTBELL)
    digest_waste=graph.product("digest_waste","exonuclease digestion products of unrepaired molecules",
                               MaterialKind.WASTE)
    graph.react("cleanup_final","Final 2× SPRI cleanup",StepKind.CLEANUP,
                ((digested,InputRole.TEMPLATE),(final_beads,InputRole.SELECTION)),
                (library,digest_waste))
    graph.discard(digest_waste,"left in the SPRI supernatant / washes")
    graph.sequence(library,SequencingPlan(
        limitation="The source publishes the hairpin adapter bases, but the separate PacBio sequencing-primer sequence is proprietary; the closed SMRTbell traversal is known. Sequenced on 8M SMRT Cells (Sequel II, 30 h movies) with polymerase 2.1 or 2.2 (Methods, DNA quality control and PacBio sequencing)."))
    graph.declare_reviewed_scope(("chromatin","ecogii","sam","mark_accessibility",
                                  "hairpin_transposome","tn5","hairpin_tagment",
                                  "stop_mix","terminate","spri_tagmented","cleanup_tagmented",
                                  "gap_repair_polymerase","gap_repair_ligase","gap_repair",
                                  "exonuclease_iii","exo_digest","spri_final","cleanup_final",
                                  "library"),SOURCE)
    return graph.require_publishable()
