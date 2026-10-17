"""NEBNext Ultra II Directional RNA library preparation, E7760."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Row, Segment
from adapter_ligation import reported_adapter_end
from endprep import dA_tailed_scene, repair_and_dA_tail
from nebnext import STANDARD_PRIMERS, indexed_library
from pcr import PCRGeometry
from protocol_graph import AnnotationSlot, Highlight, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
TITLE="NEBNext Ultra II Directional RNA"
NOTES="01_nebnext-ultra-ii-directional-rna.html"
SOURCE='Commercial protocol: <a href="https://www.neb.com/en-us/products/next-generation-sequencing-library-preparation/library-preparation-for-illumina/rna-library-prep-for-illumina">NEB E7760/E7765 manual v7.0_4/25</a> (Section 4, purified mRNA or rRNA-depleted RNA; other sections differ in RNA enrichment and fragmentation only).'
SUMMARY="Random-primed fragmented RNA is copied into cDNA in the presence of Strand Specificity Reagent, and the second strand is made with a dUTP mix. USER Enzyme is added to the adaptor ligation; removal of the dU-marked second strand before PCR is the inferred basis of strand specificity (the manual names the reagents but does not state the mechanism)."
RNA=Construct([Segment("fragmented RNA", "R"*42, placeholder=True)], name="RNA fragment")
FIRST=Construct([Segment("first-strand cDNA", "X"*42, placeholder=True)], name="RNA:cDNA product")
SECOND=[Row(chunks=[("first-strand cDNA  3′-xxxxxxxxxxxxxxxx-5′",None,False)]), Row(chunks=[("second strand      5′-xxxUxxxxUxxxxxxx-3′", "w1",False)])]
USER=[Row(chunks=[("USER cuts the dUTP-labelled second strand → it cannot seed library PCR",None,False)]), Row(chunks=[("the unmarked first strand defines the sequenced orientation",None,False)])]
END_PREP=repair_and_dA_tail(Construct([Segment("dUTP-marked cDNA duplex", "X"*42, placeholder=True)]))
ADAPTOR_END=reported_adapter_end(
    "NEBNext adaptor", kind="cohesive", polarity="3-prime", terminal_sequence="T",
    substrate_scenes=(dA_tailed_scene(END_PREP),),
    statement="The Ultra II directional workflow dA-tails the cDNA before ligating the complementary 3′-dT NEBNext adaptor; the selected oligo-set body varies by kit configuration.")
STEPS=(("Fragment RNA and synthesize the first cDNA strand", FIRST, "Random priming copies enriched or rRNA-depleted RNA fragments."),
       ("Make a dUTP-marked second strand", SECOND, "dUTP substitutes for dT during second-strand synthesis; the mark belongs only to strand two."),
       ("End-prep, ligate adaptor and remove strand two", USER, "USER Enzyme is added to the ligation mixture; excision of the dU-containing strand is inferred.", ADAPTOR_END),
       ("PCR-complete the surviving strand", Construct([Segment("orientation-preserving cDNA insert", "X"*42, placeholder=True)], name="directional insert"), "Primers from a separately purchased NEBNext oligo kit amplify the surviving strand."))
SEQ_PRIMERS=STANDARD_PRIMERS; REQUIRED_ROLES=tuple(p.role for p in SEQ_PRIMERS)
READOUT="Read 1 and Read 2 traverse opposite ends of the cDNA insert; their orientation is interpretable because the dUTP-labelled second strand was removed before amplification."
CAVEAT="INFERRED: the kit supplies no adaptor or PCR primers, so the P5/i5/Read 1 and Read 2/i7/P7 arms of the final library are drawn as one dual-index NEBNext oligo-kit configuration; the 8-nt index windows and arm sequences are not printed in the E7760 manual."
FINAL_CAPTION="INFERRED — adaptor arms and indexes depend on the separately purchased oligo kit. Dual-index directional RNA library. The molecule has ordinary TruSeq-compatible sequencing sites; directionality was imposed earlier by selective second-strand destruction."
# The prep kit contains no adaptor or PCR primers (oligo kit not selected), so every
# platform segment of the drawn layout is inferred.
def final_library(): return indexed_library(Construct([Segment("directional cDNA insert", "X"*42, placeholder=True)]), "Ultra II Directional RNA library", inferred_adaptors=True)

MANUAL="NEB E7760/E7765 manual v7.0_4/25"
FRAG_MINUTES_INTACT=(8,15)   # Table 4.1.3A, 94 C, ~200 nt inserts, RIN >= 7
RT_PROGRAM=((25,10),(42,15),(70,15))   # Step 4.2.3, minutes
SECOND_STRAND=(16,60)
END_PREP_PROGRAM=((20,30),(65,30))
LIGATION=(20,15)
USER_DIGEST=(37,15)
PCR_CYCLES=(6,13)   # Table 4.8.3A/B, input-dependent


def _cleanup(g, reviewed, ident, label, material, ratio, waste_text, kept_label, *, state=None, kind=MaterialKind.MOLECULAR):
    """A manual-specified SPRIselect / NEBNext Sample Purification Bead cleanup."""
    beads=g.input(ident+"_beads", f"{ratio} SPRIselect or NEBNext Sample Purification Beads", MaterialKind.SELECTION)
    kept=g.product(ident+"_kept", kept_label, kind, state=material.state if state is None else state, inferred=material.inferred, source=material.source)
    waste=g.product(ident+"_waste", waste_text, MaterialKind.WASTE)
    g.react(ident, label, StepKind.CLEANUP, ((material,InputRole.TEMPLATE),(beads,InputRole.SELECTION)), (kept,waste), changes_molecule=False)
    g.discard(waste, "bead supernatant removed on the magnet")
    reviewed.extend((beads.id, ident, kept.id, waste.id))
    return kept


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    reviewed=[]
    rna=g.input("rna","DNA-free poly(A)-enriched, rRNA-depleted or total RNA",MaterialKind.BIOLOGICAL,state=Construct([Segment("input RNA", "R"*42, placeholder=True)], name="input RNA"),primary=True)
    fs_buffer=g.input("fs_buffer","NEBNext First Strand Synthesis Reaction Buffer (fragmentation at 94 C)",MaterialKind.CONDITION)
    fragmented=g.product("fragmented","fragmented RNA in first-strand buffer with random primers",MaterialKind.MOLECULAR,state=RNA)
    g.react("fragment","Fragment RNA at 94 C in First Strand Synthesis Reaction Buffer (Section 4.1)",StepKind.FRAGMENTATION,((rna,InputRole.TEMPLATE),(fs_buffer,InputRole.CONDITION)),(fragmented,))
    rtprimer=g.input("rt_primer","kit Random Primers (added in the fragmentation and priming mix)",MaterialKind.OLIGO,inferred=True,source=MANUAL+" lists Random Primers as a kit component; sequence and length not printed")
    ssr=g.input("strand_specificity_reagent","NEBNext Strand Specificity Reagent (composition not printed)",MaterialKind.CONDITION)
    rt=g.input("rt","NEBNext First Strand Synthesis Enzyme Mix",MaterialKind.CATALYST)
    first=g.product("first_strand","RNA:first-strand cDNA",MaterialKind.MOLECULAR,state=FIRST)
    g.react("reverse_transcribe","Synthesize first-strand cDNA: 25 C 10 min, 42 C 15 min, 70 C 15 min",StepKind.REVERSE_TRANSCRIPTION,((fragmented,InputRole.TEMPLATE),(rtprimer,InputRole.STRUCTURAL),(ssr,InputRole.CONDITION),(rt,InputRole.CATALYST)),(first,))
    second_primer=g.input("second_primer","second-strand priming (mechanism not described by the manual)",MaterialKind.OLIGO,inferred=True,source=MANUAL+" lists only the Second Strand Synthesis Enzyme Mix and Reaction Buffer with dUTP Mix; no primer is named")
    dutp=g.input("dutp_buffer","NEBNext Second Strand Synthesis Reaction Buffer with dUTP Mix",MaterialKind.CONDITION)
    second_pol=g.input("second_polymerase","NEBNext Second Strand Synthesis Enzyme Mix",MaterialKind.CATALYST)
    duplex=g.product("duplex","dUTP-marked cDNA duplex",MaterialKind.MOLECULAR,state=END_PREP.insert)
    g.react("second_strand","Synthesize a dUTP-marked second strand: 16 C 1 h",StepKind.EXTENSION,((first,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),(dutp,InputRole.CONDITION),(second_pol,InputRole.CATALYST)),(duplex,))
    duplex=_cleanup(g,reviewed,"ds_cleanup","Purify double-stranded cDNA (1.8X beads, elute 0.1X TE)",duplex,"1.8X","second-strand reaction supernatant (enzymes, nucleotides, primers)","purified dUTP-marked cDNA duplex")
    endbuf=g.input("endprep_buffer","NEBNext Ultra II End Prep Reaction Buffer",MaterialKind.CONDITION)
    endmix=g.input("endprep_mix","NEBNext Ultra II End Prep Enzyme Mix",MaterialKind.CATALYST)
    repaired=g.product("repaired","end-repaired cDNA duplex",MaterialKind.MOLECULAR,state=END_PREP.insert)
    g.react("end_repair","Repair cDNA duplex ends (End Prep, 20 C 30 min)",StepKind.END_REPAIR,((duplex,InputRole.TEMPLATE),(endbuf,InputRole.CONDITION),(endmix,InputRole.CATALYST)),(repaired,))
    tailed=g.product("tailed","end-repaired, dA-tailed cDNA duplex",MaterialKind.MOLECULAR,state=dA_tailed_scene(END_PREP))
    g.react("da_tail","Add one 3-prime dA to each cDNA end (End Prep, then 65 C 30 min)",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(tailed,))
    adapter=g.input("adapter","NEBNext Adaptor (non-indexed; from a separately purchased NEBNext oligo kit)",MaterialKind.ADAPTER,inferred=True,source=MANUAL+" (adaptor supplied in a separate oligo kit; its structure is not printed in this manual)")
    enhancer=g.input("ligation_enhancer","NEBNext Ligation Enhancer",MaterialKind.CONDITION)
    ligase=g.input("ligase","NEBNext Ultra II Ligation Master Mix",MaterialKind.CATALYST)
    adapted=g.product("adapted","adapter-ligated dUTP-marked cDNA",MaterialKind.MOLECULAR,state=final_library(),inferred=True,source=MANUAL+" (adaptor body varies by oligo kit; drawn as the PCR-completed layout)")
    g.react("ligate","Ligate NEBNext Adaptor (20 C 15 min)",StepKind.LIGATION,((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(enhancer,InputRole.CONDITION),(ligase,InputRole.CATALYST)),(adapted,),geometry=ADAPTOR_END)
    user=g.input("user","NEBNext USER Enzyme (3 ul added to the ligation mixture)",MaterialKind.CATALYST)
    selected=g.product("selected","adapter-ligated first strand (second strand cleaved at dU)",MaterialKind.MOLECULAR,state=adapted.state,inferred=True,source=MANUAL+" adds USER to the ligation mixture but does not state the strand-selection mechanism; dU-strand excision is inferred from the dUTP Mix and USER components")
    destroyed=g.product("destroyed_second","USER-cleaved dU-containing second-strand fragments",MaterialKind.WASTE)
    g.react("destroy_second","USER digestion in the ligation mixture (37 C 15 min)",StepKind.DIGESTION,((adapted,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(selected,destroyed))
    g.discard(destroyed,"cleaved at dU sites, so it cannot template full-length library PCR (degradation, not a represented physical separation)")
    selected=_cleanup(g,reviewed,"ligation_cleanup","Purify the ligation reaction (0.9X beads, elute 0.1X TE)",selected,"0.9X","ligation supernatant (unligated adaptor, small fragments)","purified USER-treated ligation product")
    q5=g.input("pcr_polymerase","NEBNext Ultra II Q5 Master Mix",MaterialKind.CATALYST)
    fwd=g.input("pcr_forward","NEBNext Universal PCR Primer / i5 primer (oligo kit)",MaterialKind.OLIGO,inferred=True,source=MANUAL+" (oligo kit not selected; sequence not printed)")
    rev=g.input("pcr_reverse","NEBNext Index (X) Primer / i7 primer (oligo kit)",MaterialKind.OLIGO,inferred=True,source=MANUAL+" (oligo kit not selected; sequence not printed)")
    amplified=g.product("amplified","PCR-amplified directional library",MaterialKind.MOLECULAR,state=final_library())
    g.react("pcr","PCR enrichment: 98 C 30 s; 6-13 cycles 98 C 10 s / 65 C 75 s; 65 C 5 min",StepKind.PCR,((selected,InputRole.TEMPLATE),(fwd,InputRole.STRUCTURAL),(rev,InputRole.STRUCTURAL),(q5,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(final_library(),reason="The manual does not print the selected oligo kit's primer sequences."))
    library=_cleanup(g,reviewed,"pcr_cleanup","Purify the PCR reaction (0.9X beads, elute 0.1X TE)",amplified,"0.9X","PCR supernatant (primers, adaptor dimer below the cut-off)","directional dual-index RNA library",kind=MaterialKind.LIBRARY)
    g.sequence(library,SequencingPlan(library=final_library(),primers=tuple(SEQ_PRIMERS),acquired_roles=REQUIRED_ROLES))
    g.annotate(Highlight("destroy_second","USER acts in the ligation tube",
        "Every protocol section of the E7760 manual is written for a non-indexed NEBNext Adaptor and adds 3 ul USER Enzyme directly to the ligation mixture (e.g. Step 4.6.5), 37 C 15 min. "
        "The manual does not describe the adaptor structure or say what USER cleaves; this page draws only the inferred removal of the dUTP-marked second strand and does not draw any USER action on the separately supplied adaptor."))
    g.annotate(Highlight("pcr","Oligo kit not selected",
        "The library prep kit contains no adaptor or PCR primers. Single-index kits use a Universal PCR Primer (no i5 index read); dual-index kits use i5/i7 primers. "
        "The four-read dual-index layout drawn below is one oligo-kit configuration, not a requirement of E7760. Read 1/Read 2 adapter-trimming sequences are printed on manual p. 3.",
        slot=AnnotationSlot.SEQUENCING))
    g.declare_reviewed_scope(["rna","fs_buffer","fragment","rt_primer","strand_specificity_reagent","rt","reverse_transcribe","second_primer","dutp_buffer","second_polymerase","second_strand","endprep_buffer","endprep_mix","end_repair","da_tail","adapter","ligation_enhancer","ligase","ligate","user","destroy_second","pcr_forward","pcr_reverse","pcr_polymerase","pcr",*reviewed],SOURCE)
    return g.require_publishable()
