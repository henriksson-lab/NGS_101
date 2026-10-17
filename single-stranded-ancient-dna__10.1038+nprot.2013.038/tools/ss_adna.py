"""Gansauge & Meyer single-stranded ancient-DNA library model."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import illumina as il
import seqprimers as sp
from batch_ngs import seg
from adapter_ligation import AdapterLigationGeometry
from chemdraw import revcomp, Construct, Scene, feature
from pcr import PCRGeometry
from protocol_graph import Disposition, Terminal
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind)

TITLE = "Single-stranded library preparation for ancient or damaged DNA"
NOTES = "01_single-stranded-library.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nprot.2013.038">Gansauge & Meyer, <i>Nature Protocols</i> (2013)</a>.'
PROCEDURE = "Gansauge & Meyer 2013 Procedure"
SUMMARY = "Heat-denatured individual DNA strands receive a biotinylated 3′ adapter, are captured on streptavidin beads, copied into duplexes, blunt-ligated to a second adapter and completed by indexed PCR."
CAVEAT = ("Modified oligos are shown with their published modifications. The damaged insert remains a placeholder because every recovered strand can have a different length and lesion pattern. "
          "Oligos, enzymes and conditions follow the corrected (15 May 2013) Nature Protocols article: Table 1, Reagent setup and Procedure Steps 1–34. "
          "The index-2 sequencing primer (Supplementary Fig. 1 \"Read 4 seq\") is not printed in Table 1 or Step 34, which defers to Illumina's multiplex instructions; its bases are inferred from the Supplementary Fig. 1 adaptor sequence under that arrow. Read lengths are not specified by the protocol.")

# Gansauge & Meyer 2013 Table 1 (corrected-after-print PDF, nprot.2013.038.pdf p. 741),
# also in the authors' 2017 NAR Supplementary Table 2. Bases only; modifications below.
CL78 = "AGATCGGAAG"                                   # Pho-…[C3Spacer]10-TEG-biotin
CL9 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"            # extension primer
CL53 = "CGACGCTCTTCC"                                 # CGACGCTCTTC-ddC (last C dideoxy)
CL73 = "GGAAGAGCGTCGTGTAGGGAAAGAGTGTA"                # Pho-…G*T*G*T*A (phosphorothioates)
# Table 1 "Sequencing primer"; Step 34 replaces the first-read sequencing primer with CL72.
CL72 = "ACACTCTTTCCCTACACGACGCTCTTCC"
# Inferred from 2013 Supplementary Fig. 1: "Read 4 seq" (the second index read) primes back
# over the P5-side adaptor; no oligo is printed in Table 1 or Step 34.
READ4_INDEX2 = "GGAAGAGCGTCGTGTAGGGAAAGAGTGT"
# 2013 Supplementary Table 2 indexing-primer arms 3′ of the 7-nt index window
# (MOESM478 pdf.txt line 9, P5 family; line 329, P7 family).
P5_INDEX_ARM = "ACACTCTTTCCCTACACGACGCTCTT"
P7_INDEX_ARM = "GTGACTGGAGTTCAGACGTGT"
# Table 1 IS4 "P5 primer without index" = P5 joined to the P5-family arm without its
# repeated ACAC (selftest checks the transcription against Table 1).
IS4 = il.P5 + P5_INDEX_ARM[4:]
if CL78[1:] != il.INDEX1_PRIMER[:len(CL78) - 1] or not CL9.endswith(revcomp(CL78)):
    raise ValueError("CL78/CL9 no longer form the P7-side adaptor")
FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"), seg("i5", "J" * 7, "cbc", placeholder=True,
                                  feature=feature("sample_i5", "sample_index", "unknown")),
    seg("CL72 custom Read 1 site", CL72, "r1"),
    seg("ancient DNA insert", "X" * 34, placeholder=True),
    seg("CL78-derived adaptor start", CL78),
    seg("CL9 complement / Read 2 arm", il.INDEX1_PRIMER[len(CL78) - 1:], "r2"),
    seg("i7 reverse complement", "I" * 7, "cbc", placeholder=True,
        feature=feature("sample_i7", "sample_index", "unknown")),
    seg("P7 reverse complement", il.P7_RC, "p7")], name="single-stranded ancient-DNA library")
SEQ_PRIMERS = (sp.custom("Read 1", "CL72", CL72,
                         "Gansauge & Meyer 2013 Table 1 (CL72) and Procedure Step 34",
                         "Step 34: replace the sequencing primer of the first read by the custom primer CL72."),
               sp.TRUSEQ["I1"],
               sp.custom("Index 2 (i5)", "Read 4 seq (index 2)", READ4_INDEX2,
                         "Gansauge & Meyer 2013 Supplementary Fig. 1 (Read 4 seq arrow)",
                         "Inferred reverse complement of the Fig. 1 P5-side adaptor; Table 1 and Step 34 print no index-2 primer."),
               sp.TRUSEQ["R2"])
problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                     required_roles=tuple(p.role for p in SEQ_PRIMERS))
if problems: raise ValueError("ancient-DNA library: " + "; ".join(problems))
FINAL_CAPTION = "Indexed PCR supplies full flow-cell arms around the recovered ancient strand. The P5-side primer site is deliberately truncated and uses the published CL72 custom Read 1 primer."
SEQUENCING_INTRO = "CL72 replaces the standard Read 1 primer because the P5-side library arm is truncated. The read begins directly in the captured ancient DNA strand."

def oligos():
    return [
        "<b>CL78 single-stranded adapter (Table 1, corrected 15 May 2013)</b> — 5′-[Phosphate]" + CL78 + "[C3Spacer]10-[TEG-biotin]-3′",
        "<b>CL9 extension primer (Table 1)</b> — 5′-" + CL9 + "-3′",
        "<b>CL53 double-stranded adapter strand 1 (Table 1)</b> — 5′-" + CL53[:-1] + "-ddC-3′",
        "<b>CL73 double-stranded adapter strand 2 (Table 1)</b> — 5′-[Phosphate]" + CL73[:-4] + "*" + "*".join(CL73[-4:]) + "-3′ (*: phosphorothioate)",
        "<b>IS4 P5 primer without index (Table 1; replaces the P5 indexing primer for single indexing)</b> — 5′-" + IS4 + "-3′",
        "<b>CL72 custom Read 1 sequencing primer (Table 1; Step 34)</b> — 5′-" + CL72 + "-3′",
        "<b>P5 indexing primer family (2013 Suppl. Table 2)</b> — 5′-" + il.P5 + "xxxxxxx" + P5_INDEX_ARM + "-3′ (x: 7-nt index)",
        "<b>P7 indexing primer family (2013 Suppl. Table 2)</b> — 5′-" + il.P7 + "xxxxxxx" + P7_INDEX_ARM + "-3′ (x: 7-nt index)",
    ]


def _cl78_geometry():
    sc=Scene(); sc.strand("CL78",[seg("5′ phosphorylated CL78 ligating end","AGATCGGAAG"),seg("C3 spacers + TEG-biotin","[C3 spacer]10-[TEG-biotin]",placeholder=True)],label="published CL78")
    sc.ligation_end("CL78","5′ phosphorylated CL78 ligating end",kind="single-strand",polarity="5-prime",owner="adapter",text="5′-p single-strand ligating end")
    return AdapterLigationGeometry.known("single-strand",sc,statement="The protocol specifies CircLigase joining the input 3′-OH to the phosphorylated 5′ end of CL78; this is not a duplex overhang.")

def _cl53_cl73_geometry():
    sc=Scene()
    sc.strand("CL73",[seg("CL73 paired region",CL73[:len(CL53)]),
                      seg("CL73 single-stranded tail",CL73[len(CL53):])],
              label="CL73 (last four linkages phosphorothioate)",mod5="p")
    sc.anneal("CL53",[seg("CL53",CL53)],to="CL73",pair=("CL53","CL73 paired region"),
              label="CL53",mod3="ddC")
    sc.ligation_end("CL73","CL73 paired region",kind="blunt",owner="adapter",
                    text="blunt end: CL73 5′-p, CL53 3′-ddC")
    return AdapterLigationGeometry.known("blunt",sc,statement="Reagent setup hybridizes CL53 and CL73 (95 °C 10 s, 0.1 °C/s to 14 °C); CL53 pairs with the 5′ end of 5′-phosphorylated CL73 (computed), giving a blunt end. The article states that only one adapter strand is ligatable while the other is blocked by a 3′ dideoxy modification, so CL73 is joined to the newly synthesized strand.")

def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    dna=graph.input("damaged_dna","purified damaged or ancient DNA",
                    MaterialKind.MOLECULAR,primary=True)
    heat=graph.input("denaturation_heat","95 °C for 2 min, then ice-water bath",
                     MaterialKind.CONDITION,source=PROCEDURE+" Step 4")
    strand_a=graph.product("strand_a","individual denatured ancient-DNA strand A",
                           MaterialKind.MOLECULAR)
    strand_b=graph.product("strand_b","individual denatured ancient-DNA strand B",
                           MaterialKind.MOLECULAR)
    # Procedure Steps 1-4: Endonuclease VIII (+ optional Afu UDG) in CircLigase buffer II
    # with MnCl2, then FastAP, then heat denaturation in the same tube.
    uracil_enzymes=graph.input("uracil_enzymes",
                               "endonuclease VIII (0.12 U/µl) with optional Afu UDG (0.02 U/µl) in CircLigase buffer II/MnCl2",
                               MaterialKind.CATALYST,source=PROCEDURE+" Step 1")
    excised=graph.product("uracil_excised","DNA cleaved at abasic sites (and at deoxyuracils when UDG is included)",
                          MaterialKind.MOLECULAR)
    graph.react("uracil_excision","Excise deoxyuracils (optional UDG) and cleave abasic sites with endonuclease VIII (37 °C, 1 h)",
                StepKind.DIGESTION,
                ((dna,InputRole.TEMPLATE),(uracil_enzymes,InputRole.CATALYST)),(excised,),
                source=PROCEDURE+" Steps 1–2")
    fastap=graph.input("fastap","FastAP heat-labile alkaline phosphatase (1 U)",MaterialKind.CATALYST,
                       source=PROCEDURE+" Step 3")
    treated=graph.product("treated_dna","dephosphorylated DNA",MaterialKind.MOLECULAR)
    graph.react("dephosphorylate","Remove terminal phosphates with FastAP (37 °C, 10 min)",
                StepKind.DEPHOSPHORYLATION,
                ((excised,InputRole.TEMPLATE),(fastap,InputRole.CATALYST)),(treated,),
                source=PROCEDURE+" Steps 3–4")
    graph.react("denature","Heat-denature into individual strands (95 °C, 2 min; ice-water bath)",
                StepKind.STRAND_SEPARATION,
                ((treated,InputRole.TEMPLATE),(heat,InputRole.CONDITION)),
                (strand_a,strand_b),source=PROCEDURE+" Step 4")
    strands=graph.product("single_strands","pool of independently recoverable DNA strands",
                          MaterialKind.MOLECULAR)
    graph.react("collect_strands","Collect both denatured-strand populations",
                StepKind.POOL,
                ((strand_a,InputRole.TEMPLATE),(strand_b,InputRole.TEMPLATE)),
                (strands,),changes_molecule=False)
    cl78=graph.input("cl78","phosphorylated biotinylated CL78 adapter (0.125 µM)",
                     MaterialKind.ADAPTER,source="Table 1; "+PROCEDURE+" Step 5")
    circligase=graph.input("circligase","CircLigase II (5 U/µl) with 20% PEG-4000",MaterialKind.CATALYST,
                           source=PROCEDURE+" Step 5")
    adapted=graph.product("cl78_ligated","CL78-ligated ancient-DNA strands",
                          MaterialKind.MOLECULAR)
    graph.react("ligate_cl78","Ligate CL78 to individual strand 3′ ends (60 °C, 1 h; EDTA stop solution)",
                StepKind.LIGATION,
                ((strands,InputRole.TEMPLATE),(cl78,InputRole.STRUCTURAL),
                 (circligase,InputRole.CATALYST)),(adapted,),geometry=_cl78_geometry(),
                source=PROCEDURE+" Steps 5–7")
    streptavidin=graph.input("streptavidin","MyOne C1 streptavidin beads in bead-binding buffer",
                            MaterialKind.SELECTION,source=PROCEDURE+" Step 8")
    captured=graph.product("captured","bead-captured CL78-bearing strands (with excess CL78)",
                           MaterialKind.MOLECULAR)
    uncaptured=graph.product("uncaptured","unligated, unbiotinylated DNA",MaterialKind.WASTE)
    graph.react("capture","Heat 95 °C 1 min, bind to beads (20 min, RT), wash with buffers A and B",
                StepKind.AFFINITY_CAPTURE,
                ((adapted,InputRole.TEMPLATE),(streptavidin,InputRole.SELECTION)),
                (captured,uncaptured),source=PROCEDURE+" Steps 8–11")
    graph.discard(uncaptured,"removed with the bead supernatant and washes; it lacks the biotinylated CL78 adapter")
    cl9=graph.input("cl9","CL9 extension primer (2 µM)",MaterialKind.OLIGO,state=(CL9,),source="Table 1; "+PROCEDURE+" Step 12")
    copy_polymerase=graph.input("copy_polymerase","Bst 2.0 DNA polymerase (24 U) in isothermal amplification buffer",
                                MaterialKind.CATALYST,source=PROCEDURE+" Step 13")
    duplex=graph.product("copied_duplex","copied ancient-DNA duplex with 3′ overhangs",
                         MaterialKind.MOLECULAR)
    excess_primer=graph.product("excess_cl9","excess CL9, including primer bound to unligated adapter",
                                MaterialKind.WASTE)
    graph.react("copy_strand","Anneal CL9 (65 °C 2 min, ice), extend with Bst 2.0 ramping 15→37 °C at 1 °C/min, 37 °C 5 min; stringency wash at 45 °C",
                StepKind.EXTENSION,
                ((captured,InputRole.TEMPLATE),(cl9,InputRole.STRUCTURAL),
                 (copy_polymerase,InputRole.CATALYST)),(duplex,excess_primer),
                source=PROCEDURE+" Steps 12–15")
    graph.discard(excess_primer,"removed by the 45 °C stringency wash, which would otherwise allow adapter dimers in the second ligation")
    # Introduction: Bst 2.0 lacks 3'-5' exonuclease and leaves 3' overhangs; T4 DNA
    # polymerase removes them (Fig. 1 draws a 3' A).
    blunting=graph.input("overhang_removal","T4 DNA polymerase (5 U) in Buffer Tango with dNTPs",
                         MaterialKind.CATALYST,source=PROCEDURE+" Steps 16–17")
    blunt=graph.product("blunt_duplex","blunt-ended copied ancient-DNA duplex",
                        MaterialKind.MOLECULAR)
    graph.react("remove_overhang","Remove 3′ overhangs with T4 DNA polymerase (25 °C, 15 min); EDTA stop and washes",
                StepKind.END_REPAIR,
                ((duplex,InputRole.TEMPLATE),(blunting,InputRole.CATALYST)),(blunt,),
                source=PROCEDURE+" Steps 16–19")
    second_adapter=graph.input("second_adapter","CL53/CL73 blunt double-stranded adapter (2.5 µM)",
                               MaterialKind.ADAPTER,state=(CL53,CL73),source="Reagent setup; "+PROCEDURE+" Step 20")
    ligase=graph.input("blunt_ligase","T4 DNA ligase (10 U) with 5% PEG-4000",MaterialKind.CATALYST,
                       source=PROCEDURE+" Step 21")
    # No full-duplex state: before indexing PCR the molecule lacks P5/i5/i7/P7, and
    # only CL73 (not CL53, 3′-ddC) is joined; drawing FINAL_LIBRARY here overstated it.
    dual=graph.product("dual_adapter","bead-bound dual-adapter duplex (CL73 joined to the copied strand)",
                       MaterialKind.MOLECULAR)
    excess_adapter=graph.product("excess_adapter","unligated CL53/CL73 adapter",MaterialKind.WASTE)
    graph.react("ligate_second","Blunt-ligate the second adapter (1 h, room temperature); washes A, stringency (45 °C), B",
                StepKind.LIGATION,
                ((blunt,InputRole.TEMPLATE),(second_adapter,InputRole.STRUCTURAL),
                 (ligase,InputRole.CATALYST)),(dual,excess_adapter),geometry=_cl53_cl73_geometry(),
                source=PROCEDURE+" Steps 20–23")
    graph.discard(excess_adapter,"washed away from the beads")
    release_heat=graph.input("release_heat","95 °C for 1 min in 25 µl TET buffer",MaterialKind.CONDITION,
                             source=PROCEDURE+" Steps 24–25")
    released=graph.product("released_library","released library strands (CL9-primed copy carrying CL73)",
                           MaterialKind.MOLECULAR)
    released_templates=graph.product("released_templates","released fraction of biotinylated molecules, including original template strands",
                                     MaterialKind.MOLECULAR)
    retained=graph.product("bead_retained","biotinylated molecules remaining on the beads",MaterialKind.WASTE)
    graph.react("release","Heat-denature and collect the supernatant from the magnet",
                StepKind.STRAND_SEPARATION,
                ((dual,InputRole.TEMPLATE),(release_heat,InputRole.CONDITION)),
                (released,released_templates,retained),source=PROCEDURE+" Steps 24–25")
    graph.discard(retained,"left on the streptavidin beads when the supernatant is transferred")
    graph.mark_terminal(Terminal(released_templates.id, Disposition.UNRESOLVED, 'The source describes original template strands carried into the supernatant and PCR without the library priming sites. Lack of amplification is not physical separation; their later partition is not resolved by this graph.'))
    forward=graph.input("pcr_forward","P5 indexing primer (Supplementary Table 2 family) or IS4 without index (400 nM)",MaterialKind.OLIGO,
                        source="Table 1; Supplementary Table 2; "+PROCEDURE+" Step 30")
    reverse=graph.input("pcr_reverse","P7 indexing primer (Supplementary Table 2 family; 400 nM)",MaterialKind.OLIGO,
                        source="Supplementary Table 2; "+PROCEDURE+" Step 30")
    pcr_polymerase=graph.input("pcr_polymerase","AccuPrime Pfx DNA polymerase",
                               MaterialKind.CATALYST,source=PROCEDURE+" Step 30")
    amplified=graph.product("amplified_library",FINAL_LIBRARY.name+" (unpurified PCR)",MaterialKind.MOLECULAR,
                            state=FINAL_LIBRARY)
    graph.react("index_pcr","Amplify and index: 95 °C 2 min, then qPCR-chosen cycles of 95 °C 15 s / 60 °C 30 s / 68 °C 1 min",
                StepKind.PCR,
                ((released,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
                 (reverse,InputRole.STRUCTURAL),(pcr_polymerase,InputRole.CATALYST)),
                (amplified,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The page models the exact sequencing endpoint, but does not assign one index-primer member as a universal PCR reagent."),
                source=PROCEDURE+" Steps 30–31")
    purifier=graph.input("pcr_cleanup","MinElute PCR purification kit or AMPure XP beads (20 µl TE elution)",
                         MaterialKind.SELECTION,source=PROCEDURE+" Step 32")
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    pcr_waste=graph.product("pcr_waste","primers, nucleotides and polymerase removed by purification",MaterialKind.WASTE)
    graph.react("purify","Purify the amplified library",StepKind.CLEANUP,
                ((amplified,InputRole.TEMPLATE),(purifier,InputRole.SELECTION)),(library,pcr_waste),
                changes_molecule=False,source=PROCEDURE+" Step 32")
    graph.discard(pcr_waste,"removed by column or bead purification")
    graph.sequence(library,SequencingPlan(
        library=FINAL_LIBRARY,primers=tuple(SEQ_PRIMERS)))
    graph.annotate(SourceDiscrepancy(
        "cl78", "Published correction to the CL78 adapter",
        (SourceClaim("Gansauge & Meyer 2013 as initially published, Table 1 CL78 entry (quoted in the 15 May 2013 Change history)",
                     "AGATCGGAAG[C3Spacer]10[TEG-biotin] (no 5′ phosphate)"),
         SourceClaim("Nature Protocols Change history of 15 May 2013 (publisher article page)",
                     "“The correct sequence is [Phosphate]AGATCGGAAG[C3Spacer]10[TEG-biotin]”; corrected in the HTML and PDF versions"),
         SourceClaim("Corrected-after-print article PDF (nprot.2013.038.pdf), Table 1 CL78 row and Corrigendum page",
                     "[Phosphate]AGATCGGAAG[C3Spacer]10[TEG-biotin]"),
         SourceClaim("Gansauge et al. 2017 NAR Supplementary Table 2 (PMC5449542), CL78 row",
                     "Pho-AGATCGGAAG[C3Spacer]10-TEG-biotin, marked as already published in Gansauge & Meyer 2013")),
        "The model uses the 5′-phosphorylated CL78, which CircLigase requires as the donor end.",
        "The authors' later oligo table agrees with the corrected form.",
        "The initially published (uncorrected) Table 1 was not read; its value is known from the Corrigendum quotation."))
    graph.declare_reviewed_scope((
        "damaged_dna","denaturation_heat","denature","collect_strands","cl78",
        "circligase","ligate_cl78","streptavidin","capture","cl9",
        "copy_polymerase","copy_strand","uracil_enzymes","uracil_excision",
        "fastap","dephosphorylate","overhang_removal","remove_overhang","second_adapter","blunt_ligase",
        "ligate_second","release_heat","release","pcr_forward","pcr_reverse","pcr_polymerase",
        "index_pcr","amplified_library","pcr_cleanup","purify","library"),SOURCE)
    return graph.require_publishable()
