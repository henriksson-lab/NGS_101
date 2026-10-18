"""NG Capture-C native source-owned graph; accepted-manuscript Methods and published probe workbook."""
from types import SimpleNamespace
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, oligo
from adapter_ligation import AdapterLigationGeometry
from endprep import EndPreparedDuplex, dA_tailed_scene
from pcr import PCRGeometry
from protocol_model import molecule, role
from protocol_graph import (SourceClaim, SourceDiscrepancy, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal, Disposition)

TITLE="NG Capture-C — next-generation Capture-C"
SOURCE='Davies et al., <a href="https://doi.org/10.1038/nmeth.3664">Nature Methods (2016)</a>; PMC4724891 accepted manuscript and Online Methods; published Supplementary Data probe workbook.'
SUMMARY="DpnII 3C contacts are sonicated and converted with NEBNext reagents and TruSeq indexed adapters, then selected by two rounds of biotin-probe hybridization, streptavidin capture and PCR from beads."
CAVEAT="Native TruSeq adapter strands, index width/placement and PCR primer species remain unassigned. The expected sequencing library is distinct from the actual pre-index-PCR ligation product and denatured captured templates. The published Hba-1 probe is shown as an example; capture pools vary by viewpoint."
SEQUENCING="Mostly MiSeq150nt paired-end, with one HiSeq100nt paired-end run. Index acquisition and cycles remain unassigned. DpnII digestion and identical-coordinate PCR-duplicate removal in analysis are computational, not physical molecular disposal."
READ_SPECS=(sp.read_spec("Read 1","unassigned","150nt MiSeq or100nt HiSeq source alternatives; viewpoint/reporter fragment","include"),sp.read_spec("Read 2","unassigned","Opposite fragment end,150nt MiSeq or100nt HiSeq alternatives","include"),sp.read_spec("Index 1 (i7)","unassigned","TruSeq sample index; native acquisition recipe unassigned","undetermined"))
HBA_PROBE="GATCCATGGTAGCACAGGGCAGCTAAGATGCAAGTCTGAAGGAGGAGTCTGGCGAGCTGCTCCTGCAGTTCCCTGGACCCAGAAGGATGAGCTAGCAGATTCACTTGAGCCAAAGGATTC"

def contact_parts():
 return (role("captured viewpoint fragment; native boundaries unassigned"),Segment("DpnII contact junction","GATC"),role("reporter fragment; native boundaries unassigned"))

def contact_scene():
 parts=contact_parts()
 sc=Scene.duplex(parts,unpaired=tuple(x.name+"'" for x in parts if x.is_role_token()))
 sc.junction("top",parts[0].name,parts[1].name,"proximity ligation; DpnII GATC junction")
 return sc

def independent_regions(label):
 sc=Scene()
 for name in ("viewpoint","proximal reporter"):
  sc.strand(name,[role(name+" chromatin DNA; native ends unassigned")],label=label)
 return sc

FINAL_LIBRARY=Construct([role("P5-containing TruSeq left end; native geometry unassigned"),*contact_parts(),role("indexed TruSeq right end; native index placement unassigned",molecular_feature=feature("ngcapture_sample_index","sample_index","unknown",note="Sample index carried within unresolved TruSeq end; native position, width and bases unassigned"))],name="NG Capture-C source-limited indexed target library")
SEQ_PRIMERS=()

def probe_pool():
 sc=Scene()
 sc.strand("Hba-1 example",[Segment("source Hba-1 capture probe",HBA_PROBE)],label="published pool member example; viewpoint-specific",mod5="Biosg")
 sc.strand("other viewpoints",[role("other independently ordered pool members",length_bp=120)],label="separate probe population",mod5="Biosg")
 return sc

def captured_state():
 sc=Scene()
 sc.strand("captured library strand",list(FINAL_LIBRARY),label="denatured target strand; captured orientation/extent unassigned")
 sc.strand("biotin capture probe",[role("pool member complementary to chosen viewpoint",length_bp=120)],label="independent 5′-biotinylated probe; no covalent attachment",mod5="Biosg")
 sc.note("biotin capture probe","Hybridization retains target molecules; the generic pool does not assign a base-resolved binding coordinate")
 return sc

SOURCE_METHODS = ("Online Methods of the accepted author manuscript, "
                  "PMC4724891 supplementary file NIHMS65861-supplement-01")
RESTRICTION_ENZYME = "DpnII"


def protocol_graph() -> ProtocolGraph:
    """NG Capture-C graph from the Online Methods (3C, library prep, two captures)."""
    source=SOURCE
    final=FINAL_LIBRARY
    graph=ProtocolGraph(TITLE,authority=source)
    fixed = graph.input("fixed_chromatin", "2% formaldehyde-fixed cells, lysed to nuclei",
                        MaterialKind.BIOLOGICAL,
                        state=independent_regions("formaldehyde-fixed chromatin; independent DNA regions"), primary=True)
    restriction = graph.input("primary_restriction_enzyme",
                              f"{RESTRICTION_ENZYME} (3 x 500 U, 37 °C, 16–24 h, after SDS/Triton X-100)",
                              MaterialKind.CATALYST)
    digested = graph.product(
        "restriction_fragments", "DpnII-digested fixed chromatin",
        MaterialKind.MOLECULAR,
        state=independent_regions("independent DpnII fragments in fixed chromatin"))
    graph.react("primary_digest", f"Digest fixed chromatin with {RESTRICTION_ENZYME}",
                StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE),
                 (restriction, InputRole.CATALYST)), (digested,))
    ligase = graph.input("proximity_ligase",
                         "high-concentration T4 DNA ligase (Thermo Scientific), 16 °C overnight",
                         MaterialKind.CATALYST)
    contact_state=contact_scene()
    contact = graph.product("contact_product", "proximity-ligated 3C contact",
                            MaterialKind.MOLECULAR, state=contact_state)
    graph.react("proximity_ligation",
                "Heat-inactivate DpnII (65 °C) and proximity-ligate restriction fragments",
                StepKind.PROXIMITY_LIGATION,
                ((digested, InputRole.TEMPLATE),
                 (ligase, InputRole.CATALYST)), (contact,))
    proteinase = graph.input("crosslink_reversal",
                             "Proteinase K at 65 °C overnight, RNase, phenol-chloroform and ethanol precipitation",
                             MaterialKind.CONDITION)
    deproteinized = graph.product("deproteinized_3c_dna", "deproteinized 3C library DNA",
                                  MaterialKind.MOLECULAR, state=contact_state)
    graph.react("reverse_crosslinks", "Reverse cross-links and recover 3C DNA",
                StepKind.DEPROTEINIZATION,
                ((contact, InputRole.TEMPLATE),
                 (proteinase, InputRole.CONDITION)), (deproteinized,))
    sonication = graph.input("sonication",
                             "Covaris S220 sonication of 5 µg 3C DNA to about 200 bp",
                             MaterialKind.CONDITION)
    sonicated = graph.product("sonicated_3c_dna", "sonicated 3C DNA",
                              MaterialKind.MOLECULAR, state=deproteinized.state)
    graph.react("sonicate", "Sonicate recovered 3C DNA", StepKind.FRAGMENTATION,
                ((deproteinized, InputRole.TEMPLATE),
                 (sonication, InputRole.CONDITION)), (sonicated,))
    repair = graph.input("end_repair_enzymes", "NEBNext end-repair reagents",
                         MaterialKind.CATALYST)
    repaired = graph.product("end_repaired_dna", "end-repaired sonicated 3C DNA",
                             MaterialKind.MOLECULAR, state=sonicated.state)
    graph.react("end_repair", "End-repair sonicated 3C DNA", StepKind.END_REPAIR,
                ((sonicated, InputRole.TEMPLATE),
                 (repair, InputRole.CATALYST)), (repaired,))
    tailing = graph.input("a_tailing_enzyme", "NEBNext dA-tailing reagents",
                          MaterialKind.CATALYST)
    tailed_state=dA_tailed_scene(EndPreparedDuplex(Construct(contact_parts(),name="end-repaired 3C fragment"),five_prime_phosphate=True,three_prime_overhang="A"),"dA-tailed source 3C fragment; terminal phosphate derived")
    tailed = graph.product("da_tailed_dna", "dA-tailed 3C DNA",
                           MaterialKind.MOLECULAR, state=tailed_state)
    graph.react("da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE),
                 (tailing, InputRole.CATALYST)), (tailed,))
    adapters = graph.input("illumina_adapters", "Illumina TruSeq sequencing adapters",
                           MaterialKind.ADAPTER)
    adapter_ligase = graph.input("adapter_ligase", "NEBNext adapter-ligation reagents",
                                 MaterialKind.CATALYST)
    adapted = graph.product("adapter_ligated_dna", "adapter-ligated 3C DNA",
                            MaterialKind.MOLECULAR, state=molecule("actual ligated pre-index-PCR fragment",(role("native ligated TruSeq left adapter end"),*contact_parts(),role("native ligated TruSeq right adapter end; index incorporation stage unassigned"))))
    graph.react("adapter_ligation", "Ligate TruSeq adapters (AMPure XP 1:1.8 cleanups)",
                StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE),
                 (adapters, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (adapted,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The Online Methods name Illumina TruSeq indexed adapters ligated with NEBnext "
                    "reagents (E6000/E6040/E7335/E7500; E7370 for low cell numbers) but print "
                    "no adapter strands or end geometry."))
    polymerase = graph.input("pcr_polymerase", "Herculase II PCR kit (Agilent), 6–8 cycles",
                             MaterialKind.CATALYST)
    forward = graph.input("index_forward_primer", "library PCR primer 1 (sequence not printed)",
                          MaterialKind.OLIGO)
    reverse = graph.input("index_reverse_primer", "library PCR primer2; index incorporation stage unassigned",
                          MaterialKind.OLIGO)
    indexed = graph.product("indexed_library", "indexed 3C library",
                            MaterialKind.MOLECULAR, state=final)
    graph.react("index_pcr", "Index and amplify the 3C library", StepKind.PCR,
                ((adapted, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (indexed,),
                geometry=PCRGeometry.unavailable(
                    adapted.state, reason=(
                        "The Online Methods say libraries were indexed with TruSeq indexed adapters "
                        "using NEBNext reagents for adaptor ligation and PCR indexing, and that 6–8 "
                        "Herculase II cycles added the TruSeq indices; whether the index enters on "
                        "the adapter or the PCR primer, and the primer sequences, are not printed.")))

    beads = graph.input("streptavidin_beads",
                        "M270 streptavidin beads with SeqCap EZ wash buffers (47 °C stringent washes)",
                        MaterialKind.SELECTION)
    lm_polymerase = graph.input(
        "post_capture_pcr_mix",
        "SeqCap EZ Post-Capture LM-PCR master mix (18 cycles) or KAPA mix from SeqCap EZ accessory kit v2 (14 cycles)",
        MaterialKind.CATALYST)
    lm_forward = graph.input("post_capture_primer_1", "SeqCap EZ Post LM-PCR oligo, first of pair (sequence not printed)",
                             MaterialKind.OLIGO)
    lm_reverse = graph.input("post_capture_primer_2", "SeqCap EZ Post LM-PCR oligo, second of pair (sequence not printed)",
                             MaterialKind.OLIGO)
    template = indexed
    for n, hours in ((1, "64–72 h"), (2, "24 h")):
        probes = graph.input(f"capture_oligo_round{n}",
                             f"5′-biotinylated 120-nt viewpoint capture oligo pool, round {n} (13 pmol)",
                             MaterialKind.OLIGO,state=probe_pool())
        blockers = graph.input(
            f"hybridization_blockers_round{n}",
            "COT DNA, NimbleGen HE universal and index-specific blocking oligos, NimbleGen hybridization buffer",
            MaterialKind.CONDITION)
        hybrid = graph.product(f"capture_hybrid_round{n}",
                               f"round-{n} library hybridized to biotinylated probes",
                               MaterialKind.MOLECULAR, state=captured_state())
        substrate = ("1.5–2 µg indexed library, or pooled indexed libraries"
                     if n == 1 else "75% of round-1 amplified material, up to 2 µg")
        graph.react(f"capture_hybridization_round{n}",
                    f"Denature and hybridize {substrate} to capture oligos at 47 °C for {hours}",
                    StepKind.HYBRIDIZATION,
                    ((template, InputRole.TEMPLATE),
                     (probes, InputRole.STRUCTURAL),
                     (blockers, InputRole.CONDITION)), (hybrid,))
        bound = graph.product(f"capture_round{n}", f"round-{n} streptavidin-bead-bound library",
                              MaterialKind.MOLECULAR, state=hybrid.state)
        reject = graph.product(f"capture_round{n}_reject", f"uncaptured round-{n} library",
                               MaterialKind.WASTE)
        graph.react(f"capture_round{n}_step",
                    f"Round-{n} streptavidin pull-down and SeqCap EZ washes",
                    StepKind.AFFINITY_CAPTURE,
                    ((hybrid, InputRole.TEMPLATE),
                     (beads, InputRole.SELECTION)), (bound, reject))
        graph.discard(reject, "not bound to a biotinylated probe; removed with the "
                              "bead-binding and wash supernatants")
        last = n == 2
        amplified = graph.product(
            "sequencing_library" if last else f"capture_round{n}_amplified",
            TITLE + " captured library" if last else "round-1 amplified captured library",
            MaterialKind.LIBRARY if last else MaterialKind.MOLECULAR, state=final)
        residual=graph.product(f"capture_round{n}_residual",f"round-{n} original bead/probe/template material after PCR; recovery unassigned",MaterialKind.MOLECULAR,state=bound.state)
        graph.react(f"post_capture_pcr_round{n}",
                    f"Amplify round-{n} captured material from the beads (two reactions, AMPure XP cleanup)",
                    StepKind.PCR,
                    ((bound, InputRole.TEMPLATE),
                     (lm_forward, InputRole.STRUCTURAL),
                     (lm_reverse, InputRole.STRUCTURAL),
                     (lm_polymerase, InputRole.CATALYST)), (amplified,residual),
                    geometry=PCRGeometry.unavailable(
                        bound.state, reason=(
                            "Fig. 1a legend: PCR from beads uses the adapter P5 and P7 sequences; "
                            "the SeqCap EZ post-LM-PCR oligo sequences are vendor reagents not printed in the source.")))
        graph.mark_terminal(Terminal(residual.id,Disposition.UNRESOLVED,"PCR copies captured templates without consuming their probe-bound original strands; their exact partition through AMPure recovery is unassigned."))
        template = amplified
    graph.sequence(template,SequencingPlan(final,acquired_roles=("Read 1","Read 2"),limitation=SEQUENCING))
    graph.annotate(SourceDiscrepancy(
        target="capture_round2_step",
        title="Reported double-capture enrichment",
        claims=(SourceClaim("Results paragraph on single versus double capture (PMC4724891 JATS)",
                            "two sequential capture steps resulted in up to 1,000,000-fold enrichment"),
                SourceClaim("Figure 1a and Figure 2a legends (PMC4724891 JATS)",
                            "double capture produced enrichments up to 3,000,000-fold")),
        model_choice="The model encodes no fold-enrichment value; the two capture rounds are drawn without a yield claim.",
        evidence="Both values are printed in the open-access accepted manuscript; the paywalled publisher version of the main text was not inspected.",
        uncertainty="Which maximum belongs to which experiment, and whether the published version harmonizes the numbers, is unresolved."))
    graph.declare_reviewed_scope(("fixed_chromatin","primary_restriction_enzyme","primary_digest","proximity_ligase","proximity_ligation","crosslink_reversal","reverse_crosslinks","sonication","sonicate","end_repair_enzymes","end_repair","a_tailing_enzyme","da_tail","illumina_adapters","adapter_ligase","adapter_ligation","pcr_polymerase","index_forward_primer","index_reverse_primer","index_pcr","streptavidin_beads","post_capture_pcr_mix","post_capture_primer_1","post_capture_primer_2","capture_oligo_round1","hybridization_blockers_round1","capture_hybridization_round1","capture_round1_step","post_capture_pcr_round1","capture_oligo_round2","hybridization_blockers_round2","capture_hybridization_round2","capture_round2_step","post_capture_pcr_round2","sequencing_library"),source)
    return graph.require_publishable()


def oligos():
 return (oligo("Published Hba-1 viewpoint probe example (Supplementary Data, oligo sheet row94)",[Segment("Hba-1 probe",HBA_PROBE)],five="5′-/5Biosg/-"),)

SPEC=SimpleNamespace(TITLE=TITLE,NOTES="01_ng_capture_c.html",SOURCE=SOURCE,SUMMARY=SUMMARY,CAVEAT=CAVEAT,FINAL_LIBRARY=FINAL_LIBRARY,SEQ_PRIMERS=SEQ_PRIMERS,READ_LENGTHS=None,ACQUIRED_ROLES=("Read 1","Read 2"),SHOW_READ_CONTENT=False,READ_SPECS=READ_SPECS,SEQUENCING_INTRO=SEQUENCING,SEQUENCING_ENDING=SEQUENCING,SEQUENCING_UNAVAILABLE="Native TruSeq sequencing-primer and index geometry remain unassigned.",FINAL_CAPTION="Role-level native TruSeq ends around the DpnII contact fragment; no canonical single-index shell asserted.",oligos=oligos,protocol_graph=protocol_graph)
