"""CIRCLE-seq source-owned chemistry; historical primer/lot limits stay explicit."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Scene, Segment, circle_rows, complement_segments
from circular import CircularizationGeometry
from hairpin import Hairpin
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceClaim, SourceDiscrepancy, StepKind, Disposition, Terminal)
from protocol_model import Drawing, molecule, role
from topology import HairpinTopology

TITLE = "CRISPR CIRCLE-seq — in-vitro off-target cleavage"
NOTES = "01_crispr-circle-seq.html"
SOURCE = 'Defining source: <a href="https://doi.org/10.1038/nmeth.4278">Tsai et al., <i>Nature Methods</i> (2017)</a>; published Supplementary Protocol.'
SUMMARY = "Random genomic fragments are converted to covalently closed duplex circles. Cas9–guide cleavage re-creates ends that accept sequencing adapters after residual linear DNA has been digested."
CAVEAT = ('This protocol is distinct from naturally occurring eccDNA Circle-Seq. Source-printed adapter sequences are reference reagents; the historical physical E7601A kit lot is unverified. '
          '<b>Published adapter discrepancy:</b> Supplementary Protocol step 21 prints a 64-nt E7601A sequence with ACACTCTTTCCTACACGACG; the later NEB E7600 Set 1 manual prints 65 nt with ACACTCTTTCCCTACACGACG and records a November 2017 sequence correction. '
          'This model uses the defining protocol’s printed 64-nt reference, with physical ligation geometry unresolved. The differing C is in the loop, outside the paired stem. Exact historical i5/i7 primers, chosen indexes, index lengths and final primer sites remain unassigned.')
# Supplementary Protocol pp.1,6; terminal * is the final internucleotide linkage.
STEM_LOOP = "CGGTGGACCGATGATCUATCGGTCCACCGT"
PUBLISHED_LIBRARY_ADAPTER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCUACACTCTTTCCTACACGACGCTCTTCCGATCT"
STEM_MODEL = Hairpin("oSQT1288", STEM_LOOP, stem_len=12, overhang_3_len=1)
LIBRARY_MODEL = Hairpin("source-printed E7601A reference", PUBLISHED_LIBRARY_ADAPTER, stem_len=12, overhang_3_len=1)
STEM_TOPOLOGY = HairpinTopology("oSQT1288", STEM_LOOP, 0, len(STEM_LOOP)-13, 12)
LIBRARY_TOPOLOGY = HairpinTopology("source-printed E7601A reference", PUBLISHED_LIBRARY_ADAPTER, 0, len(PUBLISHED_LIBRARY_ADAPTER)-13, 12)


def hairpin_reagent(model):
    scene = Scene()
    scene.strand("hairpin", [Segment(model.name, model.sequence, "r1",
                 note="U denotes /ideoxyU/; final C*T has a phosphorothioate linkage")],
                 label=model.name + " (one folded oligo)", mod5="P")
    scene.note("hairpin", "5′-/5Phos/; internal /ideoxyU/; final C*T phosphorothioate linkage")
    scene.labels("hairpin")
    return scene


def oligos():
    def written(sequence):
        modified = sequence.replace("U", "/ideoxyU/")
        return "/5Phos/" + modified[:-1] + "*" + modified[-1]
    return ("oSQT1288 — 5′-" + written(STEM_LOOP) + "-3′",
            "E7601A source-printed reference (physical lot unverified) — 5′-" + written(PUBLISHED_LIBRARY_ADAPTER) + "-3′")


def ds_state(label, *regions):
    scene = molecule(label + " (representative strand of duplex)", tuple(role(name) for name in regions))
    scene.note("molecule", "Source-described duplex DNA; complementary strand and base-level geometry unassigned.")
    return scene


def tailed_state(label, insert):
    # Exact dA placement is intentionally not rebuilt as a falsely paired duplex.
    scene = ds_state(label, insert)
    scene.note("molecule", "Source reports 3′ dA-tailing on duplex fragment ends; terminal base geometry is not reconstructed.")
    return scene


def circle_state():
    # Two independently closed complementary strands, never ssDNA circularization.
    top = Construct([role("genomic insert"), role("USER/PNK intramolecular closure region; exact sequence unassigned")], name="circle strand 1")
    opposite = Construct(list(reversed(complement_segments(list(top)))), name="circle strand 2")
    return Drawing(tuple(circle_rows(top, "closed genomic strand 1") + circle_rows(opposite, "closed genomic strand 2")))


FINAL_LIBRARY = Construct([
    role("NEBNext i5-primer-derived library end; exact primer/index unassigned"),
    role("Cas9-cut-to-circularization-junction flank A"),
    role("original circularization junction; exact sequence unassigned"),
    role("circularization-junction-to-Cas9-cut flank B"),
    role("NEBNext i7-primer-derived library end; exact primer/index unassigned"),
], name="CIRCLE-seq source-limited dual-index endpoint")
SEQ_PRIMERS = ()
FINAL_CAPTION = "The two adapter-facing ends flank the same Cas9 cleavage site in the original circle. Their insert path passes through the original intramolecular circularization junction. Exact physical library arms are unassigned."
SEQUENCING_INTRO = "The defining Methods state 150-bp paired-end MiSeq sequencing. Index acquisition and exact sequencing-primer sites remain unassigned."
SEQUENCING_ENDING = SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE = "Source-declared paired 150-bp reads are retained without inventing exact historical primer sites or acquired index cycles."
READ_LENGTHS = {"Read 1": 150, "Read 2": 150}


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    dna = g.input("genomic_dna", "purified U2OS or HEK293 genomic DNA", MaterialKind.BIOLOGICAL,
                  state=ds_state("genomic DNA", "genomic DNA"), primary=True)
    shear = g.input("shearing", "Covaris S200, mean fragment length about 300 bp", MaterialKind.CONDITION)
    fragments = g.product("fragments", "random genomic DNA fragments", MaterialKind.MOLECULAR,
                          state=ds_state("sheared genomic fragments", "genomic insert"))
    g.react("shear", "Shear genomic DNA", StepKind.FRAGMENTATION, ((dna,InputRole.TEMPLATE),(shear,InputRole.CONDITION)), (fragments,))
    repair = g.input("repair_mix", "KAPA PCR-free library kit end-repair enzyme and buffer", MaterialKind.CATALYST)
    repaired = g.product("repaired", "end-repaired genomic fragments", MaterialKind.MOLECULAR, state=fragments.state)
    g.react("end_repair", "Repair genomic fragment ends", StepKind.END_REPAIR, ((fragments,InputRole.TEMPLATE),(repair,InputRole.CATALYST)), (repaired,))
    tail = g.input("tailing_polymerase", "KAPA A-tailing enzyme and buffer, 30 °C 30 min", MaterialKind.CATALYST)
    datp = g.input("datp", "dATP supplied for A-tailing", MaterialKind.MOLECULAR)
    tailed = g.product("da_tailed", "3′-dA-tailed genomic fragments", MaterialKind.MOLECULAR, state=tailed_state("dA-tailed genomic fragments", "genomic insert"))
    g.react("da_tail", "Add 3′ dA ends", StepKind.TAILING, ((repaired,InputRole.TEMPLATE),(datp,InputRole.STRUCTURAL),(tail,InputRole.CATALYST)), (tailed,))
    stem = g.input("stem_loop", "oSQT1288: /5Phos/, /ideoxyU/, final C*T", MaterialKind.ADAPTER,
                   state=hairpin_reagent(STEM_MODEL), chemical_state=STEM_TOPOLOGY)
    ligase = g.input("stem_ligase", "KAPA T4 DNA ligase and ligation buffer, 20 °C 1 h", MaterialKind.CATALYST)
    stemmed = g.product("stem_ligated", "genomic fragments capped with un-opened oSQT1288 hairpins", MaterialKind.MOLECULAR,
                       state=molecule("one hairpin-capped DNA molecule (fold geometry unassigned)", (role("oSQT1288 cap"),role("genomic DNA path"),role("opposite oSQT1288 cap"))))
    g.react("ligate_stem", "Ligate the USER stem-loop", StepKind.LIGATION, ((tailed,InputRole.TEMPLATE),(stem,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)), (stemmed,),
            geometry=reported_adapter_end("oSQT1288 source-printed adapter", kind="cohesive", polarity="3-prime", terminal_sequence="T", substrate_scenes=(tailed.state,),
                statement="Published oSQT1288 has a 3′T and the source explicitly dA-tails the duplex substrate. Full capped-duplex geometry is not reconstructed."))
    exo = g.input("lambda_exo", "Lambda Exonuclease plus E. coli Exonuclease I, 37 °C 1 h", MaterialKind.CATALYST)
    selected = g.product("exo_selected", "hairpin-closed molecules surviving exonuclease", MaterialKind.MOLECULAR, state=stemmed.state)
    waste = g.product("open_end_waste", "hydrolysis products of molecules with unprotected ends", MaterialKind.MOLECULAR)
    g.react("exo_select", "Digest molecules retaining an unprotected end", StepKind.DIGESTION, ((stemmed,InputRole.TEMPLATE),(exo,InputRole.CATALYST)), (selected,waste))
    g.mark_terminal(Terminal(waste.id,Disposition.UNRESOLVED,"Lambda Exonuclease/E. coli Exonuclease I degrade unprotected DNA. Hydrolysis does not establish physical removal; downstream handling of fragments/nucleotides is not followed in this graph, whose intermediate cleanups remain deferred."))
    user = g.input("user_enzyme", "USER and T4 PNK in T4 DNA ligase buffer, 37 °C 1 h", MaterialKind.CATALYST)
    opened = g.product("opened_stem", "USER-opened phosphorylated duplex circularization substrate", MaterialKind.MOLECULAR,
                       state=ds_state("opened stem-loop duplex substrate", "genomic insert with USER/PNK-processed ends; exact scar unassigned"))
    g.react("open_stem", "Excise dU and phosphorylate to expose palindromic cohesive ends", StepKind.DIGESTION, ((selected,InputRole.TEMPLATE),(user,InputRole.CATALYST)), (opened,))
    ligase2 = g.input("circle_ligase", "T4 DNA ligase, 5 ng/µl DNA, 16 °C 16 h", MaterialKind.CATALYST)
    circles = g.product("circles", "covalently closed duplex genomic circles", MaterialKind.MOLECULAR, state=circle_state())
    g.react("circularize", "Circularize the duplex fragments intramolecularly", StepKind.CIRCULARIZATION, ((opened,InputRole.TEMPLATE),(ligase2,InputRole.CATALYST)), (circles,),
            geometry=CircularizationGeometry.known("double-strand intramolecular ligation", statement="Figure 1 and the published protocol specify 4-nt palindromic overhangs, USER/PNK treatment and dilute T4 ligation; exact scar sequence is unassigned."))
    ps = g.input("exonuclease", "Plasmid-Safe ATP-dependent DNase with ATP", MaterialKind.CATALYST)
    survivors = g.product("selected_circles", "exonuclease-resistant circles", MaterialKind.MOLECULAR, state=circles.state)
    linear = g.product("linear_waste", "hydrolysis products of residual linear DNA", MaterialKind.MOLECULAR)
    g.react("digest_linear", "Digest residual linear DNA", StepKind.DIGESTION, ((circles,InputRole.TEMPLATE),(ps,InputRole.CATALYST)), (survivors,linear))
    g.mark_terminal(Terminal(linear.id,Disposition.UNRESOLVED,"Plasmid-Safe degrades residual linear DNA. Degraded products remain reaction material until a physical separation; their handling through source intermediate cleanup is outside this graph, not an asserted discard."))
    guide = g.input("guide_rna", "90 nM in-vitro-transcribed study-specific gRNA", MaterialKind.OLIGO)
    cas9 = g.input("cas9", "90 nM SpCas9 in Cas9 buffer on 250 ng circular DNA, 37 °C 1 h", MaterialKind.CATALYST)
    cut = g.product("cas9_linearized", "Cas9-linearized circles with two fresh ends", MaterialKind.MOLECULAR,
                    state=ds_state("Cas9-cleaved duplex circle", "cut-to-closure flank A", "original circularization junction", "closure-to-cut flank B"))
    uncut = g.product("uncut_circles", "uncleaved circles potentially remaining in reaction", MaterialKind.MOLECULAR, state=survivors.state)
    g.react("cas9_cut", "Cut target-bearing circles with Cas9–guide", StepKind.DIGESTION, ((survivors,InputRole.TEMPLATE),(guide,InputRole.STRUCTURAL),(cas9,InputRole.CATALYST)), (cut,uncut))
    g.mark_terminal(Terminal(uncut.id, Disposition.UNRESOLVED, "Uncleaved circles are not assigned a physical removal; they lack the fresh ends selected by adapter ligation."))
    tail2 = g.input("post_cut_tailing", "KAPA A-tailing enzyme/buffer, 30 °C 30 min", MaterialKind.CATALYST)
    datp2 = g.input("post_cut_datp", "dATP for post-cut A-tailing", MaterialKind.MOLECULAR)
    cut_tailed = g.product("cut_da_tailed", "dA-tailed Cas9-created ends", MaterialKind.MOLECULAR, state=tailed_state("dA-tailed Cas9-cleaved circle", "Cas9-cut duplex containing original circle-closure junction"))
    g.react("tail_cut", "Add 3′ dA to Cas9-created ends", StepKind.TAILING, ((cut,InputRole.TEMPLATE),(datp2,InputRole.STRUCTURAL),(tail2,InputRole.CATALYST)), (cut_tailed,))
    adapter = g.input("library_adapters", "source-printed 64-nt E7601A reference; physical historical lot unverified", MaterialKind.ADAPTER,
                      state=hairpin_reagent(LIBRARY_MODEL), chemical_state=LIBRARY_TOPOLOGY)
    ligase3 = g.input("library_ligase", "KAPA T4 DNA ligase/buffer, 20 °C 1 h", MaterialKind.CATALYST)
    ligated = g.product("adapter_ligated", "hairpin-adapted Cas9-cut DNA before USER", MaterialKind.MOLECULAR,
                        state=molecule("hairpin-capped cleavage product; historical geometry unassigned", (role("un-opened E7601A hairpin"),role("Cas9-cut circle-derived DNA"),role("opposite un-opened E7601A hairpin"))))
    from adapter_ligation import AdapterLigationGeometry
    g.react("ligate_library", "Ligate the named NEBNext hairpin adapter", StepKind.LIGATION, ((cut_tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase3,InputRole.CATALYST)), (ligated,),
            geometry=AdapterLigationGeometry.unavailable("dA substrate and 3′T adapter are source-described, but the 64/65-nt published discrepancy leaves historical physical adapter identity/complete geometry unresolved."))
    opener = g.input("post_ligation_user", "USER enzyme, supplementary protocol step 23", MaterialKind.CATALYST)
    addressable = g.product("user_opened_library", "USER-opened adapter-bearing DNA entering PCR", MaterialKind.MOLECULAR,
                           state=ds_state("actual USER-opened pre-PCR DNA", "opened NEBNext adapter-derived end", "Cas9-cut circle-derived insert", "opposite opened adapter-derived end"))
    g.react("open_library_adapter", "Open adapter dU with USER before PCR", StepKind.DIGESTION, ((ligated,InputRole.TEMPLATE),(opener,InputRole.CATALYST)), (addressable,))
    f = g.input("pcr_forward", "NEBNext Dual Index Set 1 i5 primer, member and physical sequence unassigned", MaterialKind.OLIGO)
    r = g.input("pcr_reverse", "NEBNext Dual Index Set 1 i7 primer, member and physical sequence unassigned", MaterialKind.OLIGO)
    pol = g.input("pcr_polymerase", "KAPA HiFi HotStart ReadyMix, 22 cycles", MaterialKind.CATALYST)
    amplified = g.product("pcr_product", "dual-index CIRCLE-seq amplicons", MaterialKind.MOLECULAR, state=FINAL_LIBRARY)
    g.react("library_pcr", "Amplify USER-opened DNA with i5/i7 primers", StepKind.PCR, ((addressable,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)), (amplified,),
            geometry=PCRGeometry.unavailable(addressable.state, reason="The actual USER-opened DNA is the PCR template; historical physical i5/i7 primer identities and final junction geometry are unassigned in this source scope."))
    beads = g.input("post_pcr_beads", "0.7× AMPure XP cleanup, supplementary protocol step 26", MaterialKind.SELECTION)
    library = g.product("sequencing_library", FINAL_LIBRARY.name, MaterialKind.LIBRARY, state=FINAL_LIBRARY)
    wash = g.product("post_pcr_waste", "post-PCR cleanup supernatant and washes", MaterialKind.WASTE)
    g.react("post_pcr_cleanup", "Purify PCR library on AMPure beads", StepKind.CLEANUP, ((amplified,InputRole.TEMPLATE),(beads,InputRole.SELECTION)), (library,wash)); g.discard(wash,"removed during stated bead cleanup")
    g.sequence(library, SequencingPlan(library=FINAL_LIBRARY, acquired_roles=("Read 1","Read 2"), read_lengths=READ_LENGTHS, limitation=SEQUENCING_UNAVAILABLE))
    g.annotate(SourceDiscrepancy("library_adapters", "Historical E7601A adapter sequence is unresolved",
        (SourceClaim("Published Supplementary Protocol step 21", "64 nt; ACACTCTTTCCTACACGACG"),
         SourceClaim("Later NEB E7600 Set 1 manual, adaptor overview/revision history", "65 nt; ACACTCTTTCCCTACACGACG; correction dated November 2017")),
        "Draw the defining protocol’s 64-nt printed reference; do not assign exact historical physical lot or PCR-end geometry.",
        "Discrepancy retained from prior source review; the vendor manual was not re-read in this bounded pass.",
        "The physical kit lot used in the study is unverified."))
    reviewed = ("genomic_dna","shearing","shear","repair_mix","end_repair","tailing_polymerase","datp","da_tail","stem_loop","stem_ligase","ligate_stem","lambda_exo","exo_select","user_enzyme","open_stem","circle_ligase","circularize","exonuclease","digest_linear","guide_rna","cas9","cas9_cut","post_cut_tailing","post_cut_datp","tail_cut","library_adapters","library_ligase","ligate_library","post_ligation_user","open_library_adapter","pcr_forward","pcr_reverse","pcr_polymerase","library_pcr","pcr_product","post_pcr_beads","post_pcr_cleanup","sequencing_library")
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
