"""Oxford Nanopore SQK-PCS114 cDNA-PCR workflow."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Segment, feature, revcomp
from adapter_ligation import AdapterLigationGeometry
from nanopore_library import (CRTA, RAPID_ADAPTER, RTP,
                              SSPII_AS_WRITTEN, SSPII_HANDLE, SSPII_UMI,
                              motor_sequencing_plan)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import Highlight, InputRole, MaterialKind, ProtocolGraph, StepKind
from protocol_model import molecule, role

TITLE = "Oxford Nanopore cDNA-PCR Sequencing V14 — SQK-PCS114"
NOTES = "01_cdna-pcr.html"
SOURCE = 'Commercial protocol: Oxford Nanopore <a href="https://nanoporetech.com/document/pcr-cdna-sequencing-v14-sqk-pcs114">cDNA-PCR Sequencing V14, SQK-PCS114</a>.'
SUMMARY = "Strand-switch reverse transcription adds a UMI and terminal handles, PCR selects and amplifies full-length cDNA while adding rapid-attachment tags, and a motor-loaded Rapid Adapter prepares one cDNA strand for nanopore sequencing."
CAVEAT = ("Modelled from protocol revision PCS_9200_v114_revF_31Mar2026 (the March 2026 "
          "method) and CHTD_500 rev AV. ONT publishes CRTA, RTP, SSPII, both cDNA primers "
          "and the Rapid Adapter top strand. The CRTA overhang strand, the rapid-attachment "
          "chemistry of the cPRM tags and the complete motor-bearing Rapid Adapter assembly "
          "are not disclosed in the inspected documents. Model correction (2026-10-17): the "
          "cPRM reverse primer now includes the published CTTTGTG segment that an earlier "
          "repository transcription omitted.")
# cPRM forward/reverse transcribed from CHTD_500 Rev AV section 15 (local text lines
# 2108 and 2110).  The shared lib reverse constant lacks the 7-nt CTTTGTG tag segment,
# so the published sequences are declared here.
CDNA_PRIMER_FORWARD = ("ATCGCCTACCGTGACAAGAAAGTTGTCGGTGT"
                       "CTTTGTGACTTGCCTGTCGCTCTATCTTC")
CDNA_PRIMER_REVERSE = ("ATCGCCTACCGTGACAAGAAAGTTGTCGGTGT"
                       "CTTTGTGTTTCTGTTGGTGCTGATATTGC")
# Shared 5' tag on both cPRM sequences, computed rather than typed.
CPRM_TAG = CDNA_PRIMER_FORWARD[:next(i for i, (x, y) in enumerate(
    zip(CDNA_PRIMER_FORWARD, CDNA_PRIMER_REVERSE)) if x != y)]
# 3' priming portions: forward matches the RTP 5' end, reverse the SSPII handle.
FORWARD_ANNEAL = next(n for n in range(len(RTP), 0, -1)
                      if CDNA_PRIMER_FORWARD.endswith(RTP[:n]))
REVERSE_ANNEAL = next(n for n in range(len(SSPII_HANDLE), 0, -1)
                      if CDNA_PRIMER_REVERSE.endswith(SSPII_HANDLE[:n]))
_IUPAC = str.maketrans("ACGTVBHDKMRYSWN", "TGCABVDHMKYRSWN")


def _iupac_revcomp(seq: str) -> str:
    return seq.translate(_IUPAC)[::-1]


# RTP pairs with the CRTA 3' end; reverse transcription then copies the remaining
# CRTA 5' region and the poly(A) tail before the transcript (derived, not printed).
CRTA_PAIRED = next(n for n in range(len(RTP), 0, -1) if CRTA.endswith(revcomp(RTP[:n])))
FULL_CDNA = Construct([
    Segment("RT Primer (RTP)", RTP, "r1"),
    Segment("copy of CRTA 5' region", revcomp(CRTA[:len(CRTA) - CRTA_PAIRED]),
            inferred=True),
    Segment("poly(T) copy of poly(A) tail", "T" * 12, placeholder=True),
    Segment("full-length cDNA", "X" * 36, placeholder=True),
    Segment("non-templated C opposite SSPII mGmGmG", "CCC", inferred=True),
    Segment("strand-switch UMI", _iupac_revcomp(SSPII_UMI), "umi", placeholder=True,
            feature=feature("umi", "umi", "random",
                            note="complement of the published SSPII V/T pattern")),
    Segment("SSPII handle complement", revcomp(SSPII_HANDLE), "r2"),
], name="SQK-PCS114 strand-switched cDNA as PCR substrate")
PCR_GEOMETRY = PCRGeometry.exact(
    FULL_CDNA,
    PCRPrimer("cDNA Primer forward", CDNA_PRIMER_FORWARD, anneal_nt=FORWARD_ANNEAL),
    PCRPrimer("cDNA Primer reverse", CDNA_PRIMER_REVERSE, anneal_nt=REVERSE_ANNEAL),
    statement=("Both cPRM sequences are published in ONT CHTD_500 section 15; their "
               "3' ends match the RTP 5' end and the SSPII handle."))
PCR_TEMPLATE = Construct([
    Segment("cPRM forward 5' tag", CDNA_PRIMER_FORWARD[:-FORWARD_ANNEAL], "r1"),
    *FULL_CDNA.segments,
    Segment("cPRM reverse 5' tag complement",
            revcomp(CDNA_PRIMER_REVERSE[:-REVERSE_ANNEAL]), "r2"),
], name="SQK-PCS114 PCR product")
FINAL_LIBRARY = Construct([
    Segment("motor-loaded Rapid Adapter", "[motor-loaded adapter]", "r2",
            placeholder=True, inferred=True),
    *PCR_TEMPLATE.segments,
], name="SQK-PCS114 adapted cDNA")
AMPLIFIED_CDNA = PCR_TEMPLATE
FINAL_CAPTION = (f"The PCR product carries the published cPRM 5' tags; both ends share the "
                 f"same {len(CPRM_TAG)}-nt tag, so both are attachment candidates, and the "
                 "drawing shows one adapted end. The CRTA-copy, non-templated C and "
                 "poly(T) segments are derived. The bracketed motor assembly remains "
                 "unverified pending a component-level source.")
SEQ_PRIMERS = ()
SEQUENCING_MECHANISM = "nanopore_motor"
MOTOR_ADAPTER_NAME = "Rapid Adapter"
SEQUENCING_ENDING = "There is no synthesis sequencing primer. The motor-loaded Rapid Adapter captures one strand of amplified cDNA at the pore and meters that strand through the nanopore; the adapter-bearing end enters first."

def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    rna = g.input("rna", "poly(A)+ RNA (10 ng) or total RNA (500 ng)", MaterialKind.MOLECULAR,
                  state=molecule("poly(A) RNA", (role("full-length transcript"),
                                                  role("poly(A) tail"))), primary=True)
    crta = g.input("crta", f"cDNA RT Adapter (CRTA), double stranded with poly(T) overhang; printed strand {CRTA}",
                   MaterialKind.ADAPTER,
                   source=SOURCE + " (poly(T)-overhang strand bases not printed)")
    anneal_buffer = g.input("annealing_buffer", "Annealing Buffer (AB), 60 °C 5 min then RT 5 min",
                            MaterialKind.CONDITION)
    annealed = g.product("crta_annealed", "CRTA annealed at the end of the poly(A) tail",
                         MaterialKind.MOLECULAR)
    g.react("anneal_crta", "Anneal the CRTA poly(T) overhang to the poly(A) tail end",
            StepKind.HYBRIDIZATION,
            ((rna, InputRole.TEMPLATE), (crta, InputRole.STRUCTURAL),
             (anneal_buffer, InputRole.CONDITION)), (annealed,))
    t3 = g.input("t3_ligase", "T3 DNA Ligase in StickTogether ligase buffer, RNaseOUT",
                 MaterialKind.CATALYST)
    ligated = g.product("crta_ligated", "RNA ligated to CRTA", MaterialKind.MOLECULAR)
    g.react("ligate_crta", "Ligate CRTA to the RNA (10 min RT)", StepKind.LIGATION,
            ((annealed, InputRole.TEMPLATE), (crta, InputRole.STRUCTURAL),
             (t3, InputRole.CATALYST)), (ligated,),
            geometry=AdapterLigationGeometry.unavailable(
                "The protocol says the CRTA poly(T) overhang anneals to the very end of the "
                "poly(A) tail, but does not print the overhang strand, its uracils or which "
                "strand end is joined to the RNA."))
    lexo = g.input("lambda_exo_user", "Lambda Exonuclease and USER enzyme", MaterialKind.CATALYST)
    trimmed = g.product("crta_bottom_removed", "RNA–CRTA with the CRTA bottom strand digested",
                        MaterialKind.MOLECULAR)
    g.react("digest_crta_bottom", "Digest the CRTA bottom strand (37 °C 5 min)",
            StepKind.DIGESTION,
            ((ligated, InputRole.TEMPLATE), (lexo, InputRole.CATALYST)), (trimmed,))
    rnaclean = g.input("rnaclean_xp", "RNAClean XP beads (42 µl), Short Fragment Buffer washes",
                       MaterialKind.SELECTION)
    purified = g.product("crta_purified", "bead-purified CRTA-ligated RNA", MaterialKind.MOLECULAR)
    rnaclean_waste = g.product("rnaclean_supernatant", "RNAClean XP supernatant and SFB washes",
                               MaterialKind.WASTE)
    g.react("purify_crta", "Purify on RNAClean XP beads; elute in 12 µl water", StepKind.CLEANUP,
            ((trimmed, InputRole.TEMPLATE), (rnaclean, InputRole.SELECTION)),
            (purified, rnaclean_waste))
    g.discard(rnaclean_waste, "removed with the bead supernatant and washes")
    rt_primer = g.input("rt_primer", f"RT Primer (RTP): {RTP}",
                        MaterialKind.OLIGO, state=RTP, source=SOURCE)
    rt = g.input("rt", "Maxima H Minus Reverse Transcriptase", MaterialKind.CATALYST)
    first = g.product("first_strand", "first-strand cDNA reaching the RNA 5' end",
                      MaterialKind.MOLECULAR)
    g.react("reverse_transcribe", "Anneal RTP to the CRTA (5 min RT); reverse-transcribe at 42 °C",
            StepKind.REVERSE_TRANSCRIPTION,
            ((purified, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (first,))
    tso = g.input("tso", f"Strand Switching Primer II (SSPII): {SSPII_AS_WRITTEN}",
                  MaterialKind.OLIGO, source=SOURCE)
    full = g.product("full_cdna", "full-length UMI-bearing strand-switched cDNA",
                     MaterialKind.MOLECULAR,
                     source=SOURCE + " (first-strand sequence drawn as the PCR substrate; RNA/SSPII hybrid strands not drawn)")
    g.react("template_switch", "SSPII pairs with the cDNA 5'-end C's; strand switch (42 °C, 30 min total)",
            StepKind.TEMPLATE_SWITCHING,
            ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (full,))
    heat = g.input("rt_heat_kill", "85 °C for 5 min", MaterialKind.CONDITION)
    inactivated = g.product("rt_inactivated", "heat-inactivated RT reaction",
                            MaterialKind.MOLECULAR)
    g.react("heat_inactivate", "Heat-inactivate the reverse transcriptase",
            StepKind.THERMAL_TREATMENT,
            ((full, InputRole.TEMPLATE), (heat, InputRole.CONDITION)), (inactivated,))
    ampure1 = g.input("ampure_rt", "AMPure XP beads (40 µl), 80% ethanol washes",
                      MaterialKind.SELECTION)
    rt_clean = g.product("rt_purified", "bead-purified cDNA (22.5 µl)", MaterialKind.MOLECULAR)
    rt_waste = g.product("ampure_rt_supernatant", "AMPure supernatant and ethanol washes",
                         MaterialKind.WASTE)
    g.react("purify_rt", "Purify the RT reaction on AMPure XP beads", StepKind.CLEANUP,
            ((inactivated, InputRole.TEMPLATE), (ampure1, InputRole.SELECTION)),
            (rt_clean, rt_waste))
    g.discard(rt_waste, "removed with the bead supernatant and washes")
    forward = g.input("forward_primer", f"cDNA Primer (cPRM) forward: {CDNA_PRIMER_FORWARD}",
                      MaterialKind.OLIGO, state=CDNA_PRIMER_FORWARD,
                      source=SOURCE + " (CHTD_500 section 15)")
    reverse = g.input("reverse_primer", f"cDNA Primer (cPRM) reverse: {CDNA_PRIMER_REVERSE}",
                      MaterialKind.OLIGO, state=CDNA_PRIMER_REVERSE,
                      source=SOURCE + " (CHTD_500 section 15; same cPRM tube)")
    polymerase = g.input("polymerase", "TaKaRa Ex Premier DNA Polymerase", MaterialKind.CATALYST)
    amplified = g.product("amplified", "PCR-amplified full-length cDNA",
                          MaterialKind.MOLECULAR, state=AMPLIFIED_CDNA, inferred=True,
                          source=SOURCE + " (published primers; internal derived segments)")
    g.react("pcr", "4 × 50 µl PCRs (5 µl cDNA each): 5 cycles 98/50/68 °C, then 17 cycles 98/68 °C",
            StepKind.PCR,
            ((rt_clean, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
             (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
            (amplified,), geometry=PCR_GEOMETRY)
    exo = g.input("thermolabile_exo1", "Thermolabile Exonuclease I (37 °C 5 min, 80 °C 2 min)",
                  MaterialKind.CATALYST)
    exo_done = g.product("primers_digested", "PCR product after excess-primer digestion",
                         MaterialKind.MOLECULAR, state=AMPLIFIED_CDNA, inferred=True)
    g.react("exo1", "Digest excess single-stranded primers in each PCR", StepKind.DIGESTION,
            ((amplified, InputRole.TEMPLATE), (exo, InputRole.CATALYST)), (exo_done,))
    ampure2 = g.input("ampure_pcr", "AMPure XP beads (140 µl to the 204 µl pool), 80% ethanol washes",
                      MaterialKind.SELECTION)
    pcr_clean = g.product("pcr_purified", "purified amplified cDNA in 15 µl EB",
                          MaterialKind.MOLECULAR, state=AMPLIFIED_CDNA, inferred=True)
    pcr_waste = g.product("ampure_pcr_supernatant", "AMPure supernatant and ethanol washes",
                          MaterialKind.WASTE)
    g.react("purify_pcr", "Pool the four PCRs and purify on AMPure XP beads", StepKind.CLEANUP,
            ((exo_done, InputRole.TEMPLATE), (ampure2, InputRole.SELECTION)),
            (pcr_clean, pcr_waste))
    g.discard(pcr_waste, "removed with the bead supernatant and washes")
    rapid = g.input("rapid_adapter", "Rapid Adapter (RA), diluted in Adapter Buffer (ADB)",
                    MaterialKind.ADAPTER, inferred=True,
                    source=SOURCE + f" (published RA top strand {RAPID_ADAPTER}; motor assembly unverified in inspected sources)")
    library = g.product("library", "motor-adapted amplified-cDNA library",
                        MaterialKind.LIBRARY, state=FINAL_LIBRARY, inferred=True,
                        source=SOURCE + " (complete motor-bearing assembly unverified in inspected sources)")
    g.react("attach", "Add 1 µl diluted RA to 50 fmol amplified cDNA (5 min RT)",
            StepKind.ADAPTER_ATTACHMENT,
            ((pcr_clean, InputRole.TEMPLATE), (rapid, InputRole.STRUCTURAL)), (library,))
    limitation = ("Nanopore sequencing is motor-driven rather than primer-extended; "
                  "the RA top strand is published, but the complete motor-bearing assembly "
                  "and its attachment chemistry to the cPRM tags are not.")
    g.sequence(library, motor_sequencing_plan(FINAL_LIBRARY, limitation=limitation))
    g.annotate(Highlight(
        "template_switch", "Where the non-templated C's sit",
        "SQK-PCS114 Rev F (strand-switching step, local text line 888) says SSPII pairs with deoxycytidines "
        "\"at the 5' end of the first cDNA strand\". Reverse transcriptase adds them at the "
        "growing 3' end, opposite the RNA 5' end; the model places them there and treats "
        "the vendor phrase as referring to the RNA 5' end. Their number is drawn as three "
        "to match the published mGmGmG and is not stated by the vendor."))
    g.declare_reviewed_scope((
        "rna", "crta", "annealing_buffer", "anneal_crta", "t3_ligase", "ligate_crta",
        "lambda_exo_user", "digest_crta_bottom", "rnaclean_xp", "purify_crta", "rt_primer",
        "rt", "reverse_transcribe", "tso", "template_switch", "rt_heat_kill",
        "heat_inactivate", "ampure_rt", "purify_rt", "forward_primer", "reverse_primer", "polymerase",
        "pcr", "thermolabile_exo1", "exo1", "ampure_pcr", "purify_pcr", "rapid_adapter", "attach",
        "library"), SOURCE)
    return g.require_publishable()


def oligos():
    return [
        f"CRTA: 5′-{CRTA}-3′",
        f"RTP: 5′-{RTP}-3′",
        f"SSPII (as published): 5′-{SSPII_AS_WRITTEN}-3′",
        f"cDNA Primer forward: 5′-{CDNA_PRIMER_FORWARD}-3′",
        f"cDNA Primer reverse: 5′-{CDNA_PRIMER_REVERSE}-3′",
        f"Rapid Adapter (RA) top strand / RAP T: 5′-{RAPID_ADAPTER}-3′",
    ]
