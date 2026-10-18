"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from types import SimpleNamespace
from adapter_ligation import AdapterLigationGeometry
from bliss_chem import adapter_a1 as bliss_adapter_a1, A1_BARCODE
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from illumina import TRUSEQ_SMALL_RNA_RA3
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from protocol_graph import SequencingPlan

DEFINITION = SimpleNamespace(
    title="BLISS — breaks labeling in situ and sequencing", doi="10.1038/ncomms15058",
    source_label="Yan et al., Nature Communications (2017)",
    summary="Native DSBs are blunted and ligated in situ to a T7/RA5 adapter carrying a sample barcode and UMI, then linearly amplified by IVT.",
    sequencing="For the representative classic A1 adapter, Read 1 begins with the 8-nt UMI and 8-nt sample barcode, then break-adjacent DNA. Most study libraries were single-end. Exact sequencing primer sites, acquired cycles and the chosen RPI member remain unassigned.",
)



from bliss_chem import oligos as bliss_oligos
from pcr import PCRGeometry
from protocol_graph import SourceClaim, SourceDiscrepancy
from protocol_model import molecule, role

PAPER = "Yan et al. 2017 Nat. Commun. (PMC5437291)"
PROTEX = ("authors' Protocol Exchange method, Research Square revision v2 "
          "(doi:10.21203/rs.2.1448/v2, posted 2019); the article cites "
          "protex.2017.018, whose original revision was not inspected")


def _oligos():
    return (*bliss_oligos(), oligo("RA3 RNA 3′ adapter (custom-synthesized TruSeq Small RNA sequence)",
                                  [Segment("RA3", TRUSEQ_SMALL_RNA_RA3)]))


# Classic A1 is an explicit source-selected variant, not a generic payload size.
FINAL_LIBRARY = Construct([
    role("RP1-derived library end; exact geometry unassigned"),
    role("RA5-side sequencing entry; exact primer unassigned"),
    Segment("8-nt UMI", "N" * 8, "umi", feature=feature("umi", "umi", "random")),
    Segment("A1 sample barcode", A1_BARCODE, "cbc", feature=feature(
        "sample_barcode", "sample_index", "whitelist", whitelist="BLISS Supplementary Data 1 A1")),
    role("break-adjacent DNA"),
    role("RA3/RPI-derived indexed library end; chosen RPI and exact geometry unassigned"),
], name="BLISS classic A1 endpoint (finishing geometry unassigned)")
SPEC = SimpleNamespace(
    TITLE=DEFINITION.title, NOTES="01_bliss.html",
    SOURCE=f'<a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.source_label}</a>; classic adapter A1 from Supplementary Data 1; detailed operations attributed separately to Protocol Exchange 2019 v2.',
    SUMMARY=DEFINITION.summary,
    CAVEAT="Exact A1 break-adapter strands are source-transcribed. The final RP1/RPI-installed geometry is unresolved in this model; no generic Illumina arm or index length is supplied. Detailed conditions below belong to the explicitly identified 2019 Protocol Exchange v2, whose equivalence to the article-cited original revision is unverified.",
    FINAL_LIBRARY=FINAL_LIBRARY, SEQ_PRIMERS=(), FINAL_CAPTION=DEFINITION.sequencing,
    SEQUENCING_INTRO=DEFINITION.sequencing, READ_LENGTHS=None,
    SEQUENCING_ENDING=DEFINITION.sequencing, SEQUENCING_UNAVAILABLE=DEFINITION.sequencing,
    oligos=_oligos,
)


def _state(label, *names, junctions=()):
    parts = []
    for n in names:
        if n == "UMI":
            parts.append(Segment("UMI", "N" * 8, "umi", feature=feature("umi", "umi", "random")))
        elif n in {"sample barcode", "sample barcode complement"}:
            parts.append(Segment(n, A1_BARCODE if n == "sample barcode" else revcomp(A1_BARCODE),
                                 "cbc", feature=feature("sample_barcode", "sample_index", "whitelist",
                                                      whitelist="BLISS Supplementary Data 1 A1")))
        else:
            parts.append(role(n))
    return molecule(label, tuple(parts), junctions=junctions)


def ivt_mixture_scene(dna_label):
    """IVT RNA and retained/degraded DNA are independent populations."""
    scene = Scene()
    scene.strand("rna", list(_state("RNA", "RA5", "UMI", "sample barcode", "break-adjacent RNA").strands["molecule"].segs),
                 label="independent BLISS IVT RNA")
    scene.strand("dna", [role(dna_label)], label=dna_label)
    scene.labels("rna")
    scene.labels("dna")
    return scene


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    protex = f"{source}; {PROTEX}"
    graph = ProtocolGraph(d.title, authority=source)

    dna = graph.input("dna", "formaldehyde-fixed cells, nuclei or tissue section on a coverslip, slide or well",
                      MaterialKind.BIOLOGICAL,
                      state=_state("fixed DNA", "native double-strand break"), primary=True)
    lysis = graph.input("lysis_buffers", "LB1 (1 h, 4 °C) then LB2 with 0.3% SDS (1 h, 37 °C)",
                        MaterialKind.CONDITION, source=protex)
    lysed = graph.product("lysed", "in-situ lysed, still immobilized nuclei",
                          MaterialKind.MOLECULAR,
                          state=_state("accessible fixed DNA", "native double-strand break"))
    graph.react("lyse", "Lyse and permeabilize the immobilized sample in situ",
                StepKind.PERMEABILIZATION,
                ((dna, InputRole.TEMPLATE), (lysis, InputRole.CONDITION)), (lysed,),
                source=protex, changes_molecule=False)
    repair = graph.input("blunting_enzyme", "NEB Quick Blunting Kit enzyme mix with dNTPs (1 h, room temperature)",
                         MaterialKind.CATALYST, source=protex)
    blunted = graph.product("blunted", "blunted native DNA break ends", MaterialKind.MOLECULAR,
                            state=_state("blunted break", "genomic DNA", "blunt DSB end"))
    graph.react("blunt", "Blunt the native break ends in situ", StepKind.END_REPAIR,
                ((lysed, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (blunted,))
    adapter = graph.input("break_adapter", "published UMI/sample-barcode/RA5/T7 break adapter",
                          MaterialKind.ADAPTER, state=bliss_adapter_a1(),
                          source=source + "; Supplementary Data 1 and Supplementary Fig. 1a")
    ligase = graph.input("ligase", "T4 DNA ligase (NEB M0202M), 16–18 h at 16 °C",
                         MaterialKind.CATALYST, source=protex)
    tagged = graph.product("tagged", "adapter-tagged break ends", MaterialKind.MOLECULAR,
                           state=_state("tagged break",
                               "T7 promoter", "RA5", "UMI", "sample barcode", "genomic break",
                               junctions=(("sample barcode", "genomic break", "break-adapter ligation"),)))
    free_adapter = graph.product("ligation_mix", "ligation mix with unligated adapters",
                                 MaterialKind.MOLECULAR,
                                 state=_state("unligated adapter", "free BLISS adapter"))
    graph.react("ligate", "Ligate the barcoded T7 adapter to break ends", StepKind.LIGATION,
                ((blunted, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (tagged, free_adapter),
                geometry=AdapterLigationGeometry.known(
                    "blunt", bliss_adapter_a1(),
                    statement=("Supplementary Fig. 1a marks ligation at the blunt, "
                               "5-prime-phosphorylated sample-barcode end; Supplementary "
                               "Data 1 prints both oligos.")))
    hsw = graph.input("high_salt_wash", "high-salt wash buffer (2 M NaCl), 1 h at 37 °C, repeated 3–5 times",
                      MaterialKind.SELECTION, source=protex)
    washed = graph.product("washed", "immobilized tagged nuclei after high-salt washes",
                           MaterialKind.MOLECULAR, state=tagged.state)
    hsw_waste = graph.product("hsw_waste", "high-salt washes carrying unligated adapters",
                              MaterialKind.WASTE)
    graph.react("wash_adapters", "Remove unligated adapters from the immobilized sample",
                StepKind.CLEANUP,
                ((tagged, InputRole.TEMPLATE), (free_adapter, InputRole.TEMPLATE),
                 (hsw, InputRole.SELECTION)), (washed, hsw_waste), source=protex)
    graph.discard(hsw_waste, "washes aspirated from the coverslip, slide or well")
    protk = graph.input("proteinase_k", "DNA extraction buffer (1% SDS) with proteinase K, 16–18 h at 55 °C",
                        MaterialKind.CATALYST, source=protex)
    lysate = graph.product("lysate", "scraped proteinase K lysate containing tagged genomic DNA",
                           MaterialKind.MOLECULAR, state=tagged.state)
    graph.react("extract", "Scrape and digest the sample to release genomic DNA",
                StepKind.DEPROTEINIZATION,
                ((washed, InputRole.TEMPLATE), (protk, InputRole.CATALYST)), (lysate,),
                source=protex, changes_molecule=False)
    column = graph.input("dna_cleanup", "silica column, or phenol–chloroform and isopropanol precipitation",
                         MaterialKind.SELECTION, source=protex)
    gdna = graph.product("gdna", "purified tagged genomic DNA in TE", MaterialKind.MOLECULAR,
                         state=tagged.state)
    gdna_waste = graph.product("gdna_cleanup_waste", "digested protein and extraction-buffer waste",
                               MaterialKind.WASTE)
    graph.react("purify_dna", "Purify the extracted genomic DNA", StepKind.CLEANUP,
                ((lysate, InputRole.TEMPLATE), (column, InputRole.SELECTION)),
                (gdna, gdna_waste), source=protex)
    graph.discard(gdna_waste, "column flow-through or organic phase/supernatant")
    sonicator = graph.input("sonication", "Bioruptor Plus or Covaris sonication to 300–500 bp mean size",
                            MaterialKind.CONDITION, source=protex)
    sheared = graph.product("sheared", "sheared genomic DNA; tagged fragments end at a random shear point",
                            MaterialKind.MOLECULAR,
                            state=_state("sheared tagged fragment",
                                "T7 promoter", "RA5", "UMI", "sample barcode", "break-adjacent DNA",
                                "random shear end",
                                junctions=(("sample barcode", "break-adjacent DNA", "break-adapter ligation"),)))
    graph.react("shear", "Fragment the extracted DNA by sonication", StepKind.FRAGMENTATION,
                ((gdna, InputRole.TEMPLATE), (sonicator, InputRole.CONDITION)), (sheared,),
                source=protex)
    ampure = graph.input("ampure_concentrate", "AMPure XP 0.8× (option 1; SpeedVac is the stated alternative)",
                         MaterialKind.SELECTION, source=protex)
    concentrated = graph.product("concentrated", "concentrated sheared DNA (8 µl)",
                                 MaterialKind.MOLECULAR, state=sheared.state)
    ampure_waste = graph.product("ampure_supernatant", "AMPure supernatant and ethanol washes",
                                 MaterialKind.WASTE)
    graph.react("concentrate", "Concentrate the sheared DNA on AMPure XP beads", StepKind.CLEANUP,
                ((sheared, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
                (concentrated, ampure_waste), source=protex)
    graph.discard(ampure_waste, "supernatant and washes discarded")
    t7 = graph.input("t7_polymerase", "MEGAscript T7 kit (Thermo AM1334), 14 h at 37 °C",
                     MaterialKind.CATALYST, source=protex)
    ivt = graph.product("ivt", "IVT reaction: break-adjacent RNA plus genomic DNA template",
                        MaterialKind.MOLECULAR,
                        state=ivt_mixture_scene("tagged genomic DNA template"))
    graph.react("transcribe", "Linearly amplify tagged fragments by T7 transcription",
                StepKind.TRANSCRIPTION,
                ((concentrated, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (ivt,))
    dnase = graph.input("dnase", "DNase I, RNase-free (Thermo AM2222), 15 min at 37 °C",
                        MaterialKind.CATALYST, source=protex)
    rna_digest = graph.product("dnase_treated", "IVT RNA with DNase-degraded genomic DNA",
                               MaterialKind.MOLECULAR, state=ivt_mixture_scene("DNase-degraded genomic DNA fragments"))
    graph.react("remove_gdna", "Degrade the genomic DNA template", StepKind.DIGESTION,
                ((ivt, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)), (rna_digest,),
                source=protex)
    rnaclean = graph.input("rnaclean", "RNAClean XP beads, 1.8×", MaterialKind.SELECTION,
                           source=protex)
    rna = graph.product("rna", "purified break-adjacent IVT RNA (6 µl)", MaterialKind.MOLECULAR,
                        state=_state("purified BLISS IVT RNA", "RA5", "UMI", "sample barcode", "break-adjacent RNA"))
    rna_waste = graph.product("rnaclean_supernatant",
                              "supernatant with DNase fragments and nucleotides, plus washes",
                              MaterialKind.WASTE)
    graph.react("rna_cleanup", "Purify the IVT RNA", StepKind.CLEANUP,
                ((rna_digest, InputRole.TEMPLATE), (rnaclean, InputRole.SELECTION)),
                (rna, rna_waste), source=protex)
    graph.discard(rna_waste, "supernatant and ethanol washes discarded")
    ra3 = graph.input("ra3_adapter", "custom TruSeq Small RNA RA3 3′ adapter, 1:5 diluted; denature 2 min at 70 °C",
                      MaterialKind.ADAPTER, source=protex)
    rna_ligase = graph.input("rna_ligase", "T4 RNA ligase 2, truncated (NEB M0242L), 1 h at 28 °C",
                             MaterialKind.CATALYST, source=protex)
    ligated_rna = graph.product("ra3_ligated_rna", "RA3-ligated BLISS IVT RNA",
                               MaterialKind.MOLECULAR,
                               state=_state("RA3-ligated IVT RNA",
                                            "RA5", "UMI", "sample barcode", "break-adjacent RNA", "RA3"))
    graph.react("ligate_ra3", "Ligate RA3 to the IVT RNA 3′ end", StepKind.LIGATION,
                ((rna, InputRole.TEMPLATE), (ra3, InputRole.STRUCTURAL),
                 (rna_ligase, InputRole.CATALYST)), (ligated_rna,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The custom RA3 bases are printed by Illumina, but the custom-synthesis end modifications remain unassigned within the inspected BLISS source scope."))
    rtp = graph.input("rtp", "RTP reverse-transcription primer (custom TruSeq Small RNA sequence)",
                      MaterialKind.OLIGO, source=protex)
    rt = graph.input("rt_enzyme", "SuperScript III reverse transcriptase, 1 h at 50 °C",
                     MaterialKind.CATALYST, source=protex)
    cdna = graph.product("cdna", "first-strand cDNA of the RA3-ligated IVT RNA",
                         MaterialKind.MOLECULAR,
                         state=_state("BLISS cDNA", "RA3 complement", "break-adjacent DNA",
                                      "sample barcode complement", "UMI", "RA5 complement"))
    graph.react("rt", "Reverse-transcribe from RTP", StepKind.REVERSE_TRANSCRIPTION,
                ((ligated_rna, InputRole.TEMPLATE), (rtp, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,), source=protex)
    rp1 = graph.input("rp1", "RP1 forward PCR primer (custom TruSeq Small RNA sequence)",
                      MaterialKind.OLIGO, source=protex)
    rpi = graph.input("rpi", "chosen RPI index primer (member not specified by the source)",
                      MaterialKind.OLIGO, source=protex)
    polymerase = graph.input("pcr_polymerase", "NEBNext High-Fidelity 2X PCR Master Mix; 10–16 cycles",
                             MaterialKind.CATALYST, source=protex)
    pcr = graph.product("pcr_product", "BLISS indexed PCR product", MaterialKind.MOLECULAR,
                        state=final)
    graph.react("pcr", "Index and amplify the BLISS library", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (rp1, InputRole.STRUCTURAL),
                 (rpi, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (pcr,), source=protex,
                geometry=PCRGeometry.unavailable(
                    cdna.state, reason="Actual first-strand cDNA enters PCR. The exact RP1/RPI geometry and chosen RPI member are unassigned in this source-limited model."))
    beads = graph.input("ampure_library", "AMPure XP 0.8×", MaterialKind.SELECTION, source=protex)
    library = graph.product("library", "BLISS sequencing library", MaterialKind.LIBRARY,
                            state=final)
    pcr_waste = graph.product("library_supernatant", "PCR supernatant and ethanol washes",
                              MaterialKind.WASTE)
    graph.react("library_cleanup", "Clean up the amplified library", StepKind.CLEANUP,
                ((pcr, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (library, pcr_waste), source=protex)
    graph.discard(pcr_waste, "supernatant and washes discarded")
    graph.sequence(library, SequencingPlan(library=final, acquisition="unspecified", limitation=d.sequencing))
    graph.annotate(SourceDiscrepancy(
        "break_adapter", "Two classic adapters print the same sample barcode",
        (SourceClaim("Supplementary Data 1, classic adapter A3 (top/bottom rows)",
                     "top ends …NNNNNNNNTGATGATC; bottom begins GATCATCANNNNNNNN…"),
         SourceClaim("Supplementary Data 1, classic adapter A6 (top/bottom rows)",
                     "top ends …NNNNNNNNTGATGATC; bottom begins GATCATCANNNNNNNN…")),
        "The page draws adapter A1 only; no sample is assigned to A3 or A6.",
        "Both rows are in the published supplementary workbook (PMC Cloud copy and Springer copy are identical).",
        "The source does not say whether A3 and A6 were ever used in the same pool or whether one row is a transcription error."))
    reviewed = ["dna", "lysis_buffers", "lyse", "blunting_enzyme", "blunt", "break_adapter",
                "ligase", "ligate", "high_salt_wash", "wash_adapters", "proteinase_k",
                "extract", "dna_cleanup", "purify_dna", "sonication", "shear",
                "ampure_concentrate", "concentrate", "t7_polymerase", "transcribe", "dnase",
                "remove_gdna", "rnaclean", "rna_cleanup", "ra3_adapter", "rna_ligase",
                "ligate_ra3", "rtp", "rt_enzyme", "rt", "rp1", "rpi", "pcr_polymerase",
                "pcr", "pcr_product", "ampure_library", "library_cleanup", "library"]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, SPEC.FINAL_LIBRARY, SPEC.SEQ_PRIMERS)

SPEC.protocol_graph = protocol_graph
