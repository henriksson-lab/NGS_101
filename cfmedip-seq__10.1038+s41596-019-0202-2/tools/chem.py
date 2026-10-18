"""Source-owned protocol definition and reaction graph."""
from __future__ import annotations

from types import SimpleNamespace
from chemdraw import Construct, Scene, feature
from endprep import repair_and_dA_tail, dA_tailed_scene
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole
from protocol_graph import MaterialKind
from protocol_graph import ProtocolGraph
from protocol_graph import StepKind
from pcr import PCRGeometry
from protocol_graph import Disposition, SequencingPlan, Terminal
from protocol_model import molecule, role

DEFINITION = SimpleNamespace(
    title='cfMeDIP-seq — cell-free DNA methylome enrichment',
    doi='10.1038/s41596-019-0202-2',
    source_label='Shen et al., Nature Protocols (2019)',
    summary='Adapter-ligated low-input cfDNA is supplemented with controls and filler, denatured, and enriched with a 5mC antibody without bisulfite conversion.',
    sequencing="Illumina reads report the sequence of immunoprecipitated cfDNA fragments. Procedure Step 158 sequences IP libraries paired-end 75 bp with at least 30 million reads (input-control libraries optionally); the 'Sequencing' section reports HiSeq and NextSeq runs with single-end or paired-end 50–125-bp reads (≥30–35 million per sample), paired-end when fragment size is wanted.",
    caveat='Native NEBNext kit architecture remains unassigned; independent filler and controls are not covalently joined to the cfDNA library.')

TITLE = DEFINITION.title
NOTES = "01_cfmedip_seq.html"
SOURCE = f'Defining source: <a href="https://doi.org/{DEFINITION.doi}">{DEFINITION.source_label}</a>.'
SUMMARY = DEFINITION.summary
CAVEAT = ("The protocol names NEBNext Multiplex Oligos E7335L or equivalent non-methylated adapters. "
          "Study-specific kit revision, physical adapter and primer sequences/sites remain unassigned here. "
          "A generic dual-index TruSeq shell cannot establish those molecules. Filler and spike-in controls "
          "remain independent populations; source PCR enriches adapter-bearing templates.")
FINAL_LIBRARY = Construct([role("NEB universal-primer-derived end; physical sequence unassigned", inferred=True),
    role("cfDNA insert from IP library"),
    role("NEB unique index-primer-derived end; physical sequence unassigned", inferred=True,
         molecular_feature=feature("cfmedip_sample_index","sample_index","unknown"))],
    name="cfMeDIP NEBNext library endpoint roles")
SEQ_PRIMERS = ()
FINAL_CAPTION = "Source-named NEBNext universal/unique index-primer endpoint; physical strand geometry unassigned."
SEQUENCING_INTRO = DEFINITION.sequencing
READ_LENGTHS = {"Read 1":75,"Read 2":75}
SEQUENCING_UNAVAILABLE = "Source Step158 specifies paired75 reads; native sequencing-primer sites and index acquisition remain unassigned."
SEQUENCING_ENDING = "The broader Sequencing section reports single- or paired-end50–125bp runs; this graph selects Step158 paired75."

def oligos():
    return ()

def library_mixture(label, *, controls=False, single=False, enriched=False):
    sc=Scene()
    populations=[("cfDNA library", "5mC-enriched adapter-bearing cfDNA" if enriched else "adapter-bearing cfDNA"),
                 ("filler", "methylated/unmethylated lambda filler; not adapter-ligated")]
    if controls: populations.extend((("methylated control","methylated A. thaliana spike-in"),
                                     ("unmethylated control","unmethylated A. thaliana spike-in")))
    for key,description in populations:
        if single:
            for polarity in ("one polarity","complementary polarity"):
                sc.strand(key+polarity,[role(description+"; "+polarity)],label=label+": "+key+" strands")
        else:
            sc.strand(key,[role(description+"; native duplex architecture unresolved")],label=label+": "+key)
    return sc


def _local_graph(d, final, primers):
    source = f"{d.source_label}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    dna = graph.input("cfdna", "naturally fragmented plasma cell-free DNA",
                      MaterialKind.BIOLOGICAL,
                      state=Scene.duplex([role("naturally fragmented cfDNA")],unpaired=("naturally fragmented cfDNA'",)), primary=True)
    repair = graph.input("repair_mix", "cfDNA end-repair enzyme mix", MaterialKind.CATALYST)
    repaired = graph.product("repaired", "end-repaired cfDNA", MaterialKind.MOLECULAR,
                             state=dna.state)
    graph.react("end_repair", "End-repair cfDNA fragments", StepKind.END_REPAIR,
                ((dna, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("tailing_enzyme", "dA-tailing enzyme", MaterialKind.CATALYST)
    tailed = graph.product("tailed", "3-prime-dA-tailed cfDNA", MaterialKind.MOLECULAR,
                           state=dA_tailed_scene(repair_and_dA_tail(Construct([role("cfDNA")])) ))
    graph.react("tail", "Add one 3-prime dA to each fragment end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("adapters", "NEBNext Multiplex Oligos for Illumina adaptor (sequence and end geometry not printed)",
                           MaterialKind.ADAPTER, inferred=True,
                           source=source + " (Procedure Step 28 names the kit; sequences not printed)")
    ligase = graph.input("ligase", "adapter DNA ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "adapter-ligated cfDNA", MaterialKind.MOLECULAR,
                            state=molecule("adapter-ligated cfDNA; native geometry unresolved", (
                                role("NEBNext adapter-bearing cfDNA; ligating strands and branches unassigned", inferred=True),)))
    graph.react("ligate", "INFERRED adaptor sequence — ligate the NEBNext adaptor to dA-tailed cfDNA", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The article names the NEBNext Multiplex Oligos adaptor (Step 28) but prints neither its sequence nor its end geometry, and the kit revision is not stated. dA/dT ligation is the summary-level description."))
    lig_beads=graph.input("ligation_beads","88µL AMPure XP into110µL ligation; two80%ethanol washes;27µL elution",MaterialKind.SELECTION)
    lig_clean=graph.product("ligation_clean","purified adapter-ligated cfDNA before USER",MaterialKind.MOLECULAR,state=adapted.state)
    lig_waste=graph.product("ligation_waste","post-ligation bead supernatant and ethanol washes",MaterialKind.WASTE)
    graph.react("ligation_cleanup","Post-ligation0.8× AMPure cleanup (Steps31–41)",StepKind.CLEANUP,
        ((adapted,InputRole.TEMPLATE),(lig_beads,InputRole.SELECTION)),(lig_clean,lig_waste),changes_molecule=False)
    graph.discard(lig_waste,"Steps35–38 supernatant/ethanol physically removed")
    # Procedure Step 42: the NEBNext adaptor "must be cleaved with USER enzyme, in order
    # to ensure proper separation of the DNA strand under denaturation conditions";
    # 37 C 15 min, then MinElute purification. The adaptor bases stay unprinted.
    user = graph.input("user_enzyme", "USER enzyme (NEBNext Multiplex Oligos kit) with T4 DNA ligase buffer",
                       MaterialKind.CATALYST, source=source + " (Procedure Step 42)")
    opened = graph.product("user_opened", "USER-cleaved adapter-ligated cfDNA", MaterialKind.MOLECULAR,
                           state=adapted.state)
    graph.react("user_cleave", "Cleave the NEBNext adaptor with USER enzyme (37 °C, 15 min)",
                StepKind.DIGESTION, ((lig_clean, InputRole.TEMPLATE), (user, InputRole.CATALYST)),
                (opened,), source=source + " (Procedure Step 42)")
    minelute=graph.input("minelute","MinElute PCR Purification Kit;40µL EB elution (Step44)",MaterialKind.SELECTION)
    opened_clean=graph.product("opened_clean","purified USER-opened library",MaterialKind.MOLECULAR,state=opened.state)
    minelute_waste=graph.product("minelute_waste","MinElute flowthrough/washes",MaterialKind.WASTE)
    graph.react("user_cleanup","Purify USER-treated sample with MinElute",StepKind.CLEANUP,
        ((opened,InputRole.TEMPLATE),(minelute,InputRole.SELECTION)),(opened_clean,minelute_waste),changes_molecule=False)
    graph.discard(minelute_waste,"excluded by source column purification; exact component recovery unassigned")
    # Procedure Step 45 and 'Filler DNA generation': six lambda PCR amplicons (Table 1),
    # 1/5/10/15/20L-CpG M.SssI-methylated and 20S-CpG unmethylated, 50:50 w/w, topping the
    # DNA up to 100 ng. Filler is not adapter-ligated and yields no reads.
    filler = graph.input("filler", "λ filler DNA: M.SssI-methylated and unmethylated λ PCR amplicons, 50:50 w/w, to 100 ng total (not adapter-ligated)",
                         MaterialKind.MOLECULAR,
                         source=source + " (Procedure Step 45; Table 1; 'Filler DNA generation')")
    with_filler = graph.product("with_filler", "adapter-ligated cfDNA topped up with λ filler",
                                MaterialKind.MOLECULAR,
                                state=library_mixture("library and filler"))
    graph.react("add_filler", "Top up to 100 ng with λ filler DNA", StepKind.POOL,
                ((opened_clean, InputRole.TEMPLATE), (filler, InputRole.TEMPLATE)), (with_filler,))
    # Procedure Step 48: 1.5 ul each of methylated and unmethylated A. thaliana control
    # DNA (0.2 ng/ul; MagMeDIP kit) with MagBuffer A and B, before heat denaturation.
    controls = graph.input("controls", "methylated and unmethylated A. thaliana spike-in DNA, 0.3 ng each (Diagenode MagMeDIP kit)",
                           MaterialKind.MOLECULAR,
                           source=source + " (Procedure Step 48; 'QC of cfMeDIP-seq')")
    mixture = graph.product("mixture", "low-input cfDNA library mixture",
                            MaterialKind.MOLECULAR,
                            state=library_mixture("library mixture",controls=True))
    graph.react("pool_controls", "Add methylated and unmethylated A. thaliana spike-ins",
                StepKind.POOL,
                ((with_filler, InputRole.TEMPLATE), (controls, InputRole.TEMPLATE)), (mixture,))
    denaturing = graph.input("denaturing_condition", "95 °C for 10 min, then snap-cooling in ice-water slush (Steps 51–54)",
                             MaterialKind.CONDITION)
    single = graph.product("single_strand", "denatured cfDNA library strands of both polarities",
                           MaterialKind.MOLECULAR, state=library_mixture("denatured mixture",controls=True,single=True))
    plus = graph.product("plus_strands", "denatured library strands, one polarity", MaterialKind.MOLECULAR)
    minus = graph.product("minus_strands", "denatured library strands, complementary polarity", MaterialKind.MOLECULAR)
    graph.react("denature", "Denature the low-input library mixture",
                StepKind.STRAND_SEPARATION,
                ((mixture, InputRole.TEMPLATE), (denaturing, InputRole.CONDITION)),
                (plus, minus))
    graph.react("retain_both_polarities", "Account for both strand polarities in the same tube (no transfer)",
                StepKind.POOL, ((plus, InputRole.TEMPLATE), (minus, InputRole.TEMPLATE)), (single,))
    # Steps 62-64: 7.9 ul (10%) of the denatured sample becomes the input control (IC),
    # 79 ul the IP sample. The IC is purified side by side (Steps 89, 96-116), is the
    # reference for QC1/QC2 qPCR (Steps 117-128) and may optionally be sequenced (Step 158).
    ip_sample = graph.product("ip_sample", "denatured IP sample (79 µL)", MaterialKind.MOLECULAR,
                              state=single.state)
    input_control = graph.product("input_control", "10% input control (IC, 7.9 µL)", MaterialKind.QC,state=single.state)
    graph.react("split_ic", "Set aside 10% of the denatured sample as the input control (IC; Steps 62–64)",
                StepKind.SPLIT, ((single, InputRole.TEMPLATE),), (ip_sample, input_control),
                source=source + " (Procedure Steps 62–64)")
    graph.mark_terminal(Terminal(input_control.id, Disposition.QC,
        "IC is purified with the IPure kit beside the IP sample and used as the qPCR recovery and cycle-number reference (Steps 89, 96–128); Step 158 says IC samples can also be sequenced when copy-number or enrichment metrics are needed, which is not modelled as a separate library here"))
    antibody = graph.input("antibody", "Diagenode MagMeDIP-kit 5-mC monoclonal antibody 33D3 (C15200081), 0.16 µg per reaction, with MagBuffer C and washed MagMeDIP Magbeads",
                           MaterialKind.SELECTION,
                           source=source + " (published Reporting Summary, 'Antibodies used'; Procedure Steps 56–71)")
    enriched = graph.product("enriched", "bead-bound 5mC-enriched adapter-ligated cfDNA strands",
                             MaterialKind.MOLECULAR, state=library_mixture("IP-enriched mixture; recovery of each population unresolved",controls=True,single=True,enriched=True))
    depleted = graph.product("depleted", "unbound supernatant and MagWash buffer 1/2 washes",
                             MaterialKind.WASTE)
    graph.react("immunoprecipitate", "Immunoprecipitate 5mC-containing strands (4 °C, 17 h, 10 r.p.m.), then wash 3× MagWash Buffer-1 and 1× MagWash Buffer-2 on an ice-cold magnet",
                StepKind.AFFINITY_CAPTURE,
                ((ip_sample, InputRole.TEMPLATE), (antibody, InputRole.SELECTION)),
                (enriched, depleted), source=source + " (Procedure Steps 64–86)")
    graph.discard(depleted, "supernatant and washes removed from the bead-bound methylated fraction (Steps 79–86)")
    ipure = graph.input("ipure", "Diagenode IPure Kit v2: IP-EB elution buffer (buffer A + B), carrier, isopropanol, magnetic beads, wash buffers 1/2, buffer C",
                        MaterialKind.CONDITION, source=source + " (Procedure Steps 73–75, 87–116)")
    eluted = graph.product("eluted", "purified single-stranded IP DNA in 50 µL buffer C",
                           MaterialKind.MOLECULAR, state=enriched.state)
    ipure_waste = graph.product("ipure_waste", "MagMeDIP beads and IPure bead washes",
                                MaterialKind.WASTE)
    graph.react("elute_purify", "Elute twice with IP-EB (2 × 50 µL, room temperature), then bind to IPure beads with carrier and isopropanol, wash, and elute in buffer C (25 + 23 µL, topped to 50 µL)",
                StepKind.CLEANUP, ((enriched, InputRole.TEMPLATE), (ipure, InputRole.CONDITION)),
                (eluted, ipure_waste), source=source + " (Procedure Steps 87–116)")
    graph.discard(ipure_waste, "beads and wash supernatants removed during IP elution and IPure purification (Steps 92–113)")
    qc_rule=graph.input("qc_aliquot_rule","QC1 four1µL aliquots; QC2 two1µL aliquots; retain44µL (Steps117–129)",MaterialKind.CONDITION)
    production=graph.product("production_template","44µL purified IP DNA retained for library PCR",MaterialKind.MOLECULAR,state=eluted.state)
    qc_aliquots=graph.product("qc_aliquots","six1µL aliquots consumed by QC1/QC2 qPCR",MaterialKind.QC,state=eluted.state)
    graph.react("qc_sampling","Sample duplicate recovery/specificity qPCR and cycle-number qPCR",StepKind.SAMPLING,
        ((eluted,InputRole.TEMPLATE),(qc_rule,InputRole.CONDITION)),(production,qc_aliquots),changes_molecule=False,
        source=source+" Steps117–129; QC1 two target pairs in duplicate, QC2 duplicate")
    graph.mark_terminal(Terminal(qc_aliquots.id,Disposition.QC,"QC1 compares methylated/unmethylated spike-in recovery; QC2 uses pooled indexing primers only for analytical cycle determination. These aliquots do not return to production PCR."))
    reviewed = ["ligation_beads","ligation_cleanup","minelute","user_cleanup","qc_aliquot_rule","qc_sampling","qc_aliquots","production_template","cfdna", "repair_mix", "end_repair", "tailing_enzyme", "tail", "adapters",
                "ligase", "ligate", "user_enzyme", "user_cleave", "filler", "add_filler",
                "controls", "pool_controls", "denaturing_condition",
                "denature", "retain_both_polarities", "ip_sample", "input_control", "split_ic", "antibody", "enriched",
                "depleted", "immunoprecipitate", "ipure", "eluted", "ipure_waste", "elute_purify"]
    # Steps 124-131: PCR (inlined so the source post-PCR cleanups can follow it).
    primer_source = source + " (Procedure Steps 124–131)"
    forward = graph.input("library_forward", "NEB universal primer (NEBNext Multiplex Oligos for Illumina; sequence not printed)",
                          MaterialKind.OLIGO, source=primer_source, inferred=True)
    reverse = graph.input("library_reverse", "NEB index primer, one unique index per sample (NEBNext Multiplex Oligos; sequence not printed)",
                          MaterialKind.OLIGO, source=primer_source, inferred=True)
    polymerase = graph.input("library_polymerase", "KAPA HiFi HotStart ReadyMix", MaterialKind.CATALYST)
    amplified = graph.product("amplified", "amplified double-stranded cfMeDIP-seq library (2 × 55 µL)",
                              MaterialKind.MOLECULAR, state=final)
    residual=graph.product("coip_nonlibrary","co-IP filler/control material outside amplified cfDNA endpoint",MaterialKind.MOLECULAR)
    graph.react("library_pcr", "INFERRED primer sequences — amplify and index the 44 µL eluate with the NEB universal and a unique NEB index primer (KAPA HiFi HotStart; 95 °C 3 min, then the QC2 qPCR-determined number of 98 °C 20 s / 65 °C 15 s / 72 °C 30 s cycles, ≤15 for IP samples; 72 °C 1 min)",
                StepKind.PCR,
                ((production, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,residual), source=primer_source,
                geometry=PCRGeometry.unavailable(production.state, reason=
                    "The actual input is the purified single-stranded IP mixture. Physical NEBNext primer sites and first-copy geometry remain undeclared; completed library endpoints do not supply those coordinates."))
    graph.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,"Lambda filler is not adapter-ligated and is not a NEBNext library-PCR target. Non-amplification does not prove physical removal; control/filler recovery through later cleanup is not quantified here."))
    beads = graph.input("ampure", "AMPure XP beads, 80% ethanol, low TE", MaterialKind.CONDITION,
                        source=source + " (Procedure Steps 132–156)")
    cleaned = graph.product("cleaned", "1× AMPure-purified library (100 µL low TE)",
                            MaterialKind.MOLECULAR, state=final)
    cleanup_waste = graph.product("cleanup_waste", "1× AMPure supernatant (primers, unincorporated reagents)",
                                  MaterialKind.WASTE)
    graph.react("ampure_cleanup", "Pool the two PCR tubes and purify with 1× AMPure XP (Steps 132–143)",
                StepKind.CLEANUP, ((amplified, InputRole.TEMPLATE), (beads, InputRole.CONDITION)),
                (cleaned, cleanup_waste), source=source + " (Procedure Steps 132–143)")
    graph.discard(cleanup_waste, "1× AMPure supernatant removed")
    size_beads = graph.input("size_beads", "AMPure XP beads at 0.6× then 1.2× of the 100 µL sample volume",
                             MaterialKind.SELECTION, source=source + " (Procedure Steps 144–148)")
    library = graph.product("library_library", "cfMeDIP-seq sequencing library (20 µL low TE)",
                            MaterialKind.LIBRARY, state=final)
    large = graph.product("large_fragments", "large fragments bound at 0.6× AMPure", MaterialKind.WASTE)
    small = graph.product("small_fragments", "small fragments left unbound at 1.2× AMPure", MaterialKind.WASTE)
    graph.react("dual_size_selection", "Dual size selection: 0.6× AMPure XP (discard beads), then 1.2× AMPure XP on the supernatant; elute in 20 µL low TE (Steps 144–156)",
                StepKind.SELECTION, ((cleaned, InputRole.TEMPLATE), (size_beads, InputRole.SELECTION)),
                (library, large, small), source=source + " (Procedure Steps 144–156)")
    graph.discard(large, "0.6× AMPure beads with the larger unwanted fragments are discarded (Step 147)")
    graph.discard(small, "supernatant of the 1.2× AMPure binding is removed (Step 152)")
    plan = SequencingPlan(library=final,acquired_roles=("Read 1","Read 2"),read_lengths=READ_LENGTHS,
            limitation=SEQUENCING_UNAVAILABLE+" "+SEQUENCING_ENDING)
    graph.sequence(library, plan)
    reviewed.extend(("library_forward", "library_reverse", "library_polymerase", "amplified",
                     "library_pcr", "ampure", "cleaned", "cleanup_waste", "ampure_cleanup",
                     "size_beads", "library_library", "large_fragments", "small_fragments",
                     "dual_size_selection"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    return _local_graph(DEFINITION, FINAL_LIBRARY, SEQ_PRIMERS)
