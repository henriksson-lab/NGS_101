"""Oudelaar2018 Tri-C: native contact and capture declarations, no shared endpoint recipe."""
from types import SimpleNamespace
from chemdraw import Construct, Scene, Segment, feature, complement_segments
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind,InputRole,StepKind,SequencingPlan,Terminal,Disposition,SourceClaim,SourceDiscrepancy
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from endprep import EndPreparedDuplex,dA_tailed_scene
import seqprimers as sp
TITLE='Tri-C — targeted multiway chromatin contacts'
NOTES='01_tri_c.html'
SOURCE='Oudelaar et al., Nature Genetics2018, <a href="https://doi.org/10.1038/s41588-018-0253-2">doi:10.1038/s41588-018-0253-2</a>; defining Methods and cited2018 detailed protocol.'
SUMMARY='Native NlaIII concatemers are gently sheared, adapted and indexed before two viewpoint-capture rounds. Captured DNA is amplified from streptavidin beads; capture probes remain separate molecules.'
CAVEAT="""Published source boundaries: article Methods uses2.9nM total capture-probe pool; cited2018 ProtocolExchange v1 Capture enrichment p.7 says2.9nM per oligonucleotide. The2018 materials list p.2 names E6040, while its library-prep heading p.5 says E7370. The2019 v2 revision names E6040 and presents E7645 UltraII as a later low-input option; it does not identify the historical kit lot. Both detailed versions state10nmol in3.46µl gives2.9µM; arithmetic gives about2.9mM. The graph uses the source-listed end-repair/dA/ligation operations with actual kit identity and probe concentration unassigned. The article0.7:1 post-sonication cleanup differs from the2018 protocol1.8× instruction; the graph follows the article for that step. Superseded preprint values are not current authority. Native adapters, selected index members and PCR/USER end boundaries remain unassigned."""
SHEARING_ID=feature('sonication_endpoint','umi','random',note='Physical genomic shearing coordinate used computationally for deduplication; no synthesized nucleotide barcode.')
SAMPLE_ID=feature('tri_c_sample','sample_index','unknown',note='Differential library index; native member, count and width unassigned.')
ENDPOINT=role('physical sonication endpoint; coordinate identifier, no added barcode bases',molecular_feature=SHEARING_ID)
def contact_parts():
 return [role('viewpoint fragment'),Segment('first NlaIII ligation junction','CATG'),role('reporter fragment1'),Segment('second NlaIII ligation junction','CATG'),role('reporter fragment2')]
FINAL_LIBRARY=Construct([role('native Illumina end A'),ENDPOINT,*contact_parts(),role('native indexed Illumina end B',molecular_feature=SAMPLE_ID)],name='Tri-C source-limited captured endpoint')
SEQ_PRIMERS=()
FINAL_CAPTION='Illustrative multiway endpoint, with physical shearing-coordinate identity and native adapter/index roles. It is not an assumed canonical TruSeq library.'
SEQUENCING_INTRO='First erythroid experiment MiSeq v2 PE250; second multiplexed erythroid/ES experiment NextSeq v2 PE150. Index acquisition remains unassigned.'
SEQUENCING_ENDING=SEQUENCING_INTRO
READ_SPECS=(sp.read_spec('Read 1','250 first experiment;150 second','viewpoint/reporter contacts; computational duplicate filtering by physical coordinates','include'),sp.read_spec('Read 2','250 first experiment;150 second','viewpoint/reporter contacts','include'),sp.read_spec('Index 1 (i7)','unassigned','Native sample index acquisition unassigned','undetermined'),sp.read_spec('Index 2 (i5)','unassigned','Native index acquisition unassigned','undetermined'))
DEFINITION=SimpleNamespace(title=TITLE,doi='10.1038/s41588-018-0253-2',citation='Oudelaar et al., Nature Genetics2018')

def duplex(parts,label):
 return Scene.duplex(parts,label=label,unpaired=tuple(p.name+"'" for p in parts if p.is_role_token()))
def fixed_state():
 sc=Scene();sc.strand('chromatin',[role('chromosomal DNA; unresolved locus topology before digestion')],label='fixed chromatin population');return sc
def fragment_state():
 sc=Scene()
 for i,label in enumerate(('viewpoint fragment','reporter fragment1','reporter fragment2')):
  body=role(label+' body; locus sequence/length unassigned')
  opposite=complement_segments([body])[0]
  tail=Segment(label+' top3′CATG','CATG');bottom_tail=Segment(label+' bottom3′CATG','CATG')
  sc.strand('top'+str(i),[body,tail],label='NlaIII fragment; cleavage-derived5′P',mod5='p')
  sc.anneal('bottom'+str(i),[opposite,bottom_tail],to='top'+str(i),pair=(opposite.name,body.name),unpaired=(opposite.name,),label='NlaIII fragment; cleavage-derived5′P',mod5='p')
  sc.overhang('top'+str(i),tail.name,polarity='3-prime',owner='substrate',text='NlaIII3′CATG')
  sc.overhang('bottom'+str(i),bottom_tail.name,polarity='3-prime',owner='substrate',text='NlaIII3′CATG')
 return sc

def molecule(label,parts,**kw):
 sc=Scene();sc.strand('DNA',[p if isinstance(p,Segment) else role(p) for p in parts],label=label);return sc

def multiway_contact(enzyme,viewpoint,reporters):
 if enzyme!='NlaIII' or len(reporters)<2:raise ValueError('Local Tri-C example requires source NlaIII and at least two reporters')
 return SimpleNamespace(scene=lambda label:duplex(contact_parts(),label))
def da_tailed(label,insert_name):
 return dA_tailed_scene(EndPreparedDuplex(Construct([ENDPOINT,*contact_parts()],name=insert_name),five_prime_phosphate=False,three_prime_overhang='A'),label=label+';5′phosphate assignment unassigned')
def adapted_state(*,opened=False):
 return duplex([role('native adaptor end A; '+('USER end boundary unassigned' if opened else 'ligated topology unassigned')),ENDPOINT,*contact_parts(),role('native adaptor end B; '+('USER end boundary unassigned' if opened else 'ligated topology unassigned'))],'native adapted 3C precursor')
def post_pcr_state(*,capture_round=0):
 sc=duplex(list(FINAL_LIBRARY),'expected amplified multiway library population')
 sc.strand('native templates',[role('possible co-carried native adapted/captured templates; abundance, integrity and partition unassigned')],label='independent residual template lineage; no physical separation asserted')
 if capture_round:
  sc.strand('capture probes',[role('possible co-carried biotinylated capture oligos through round'+str(capture_round)+'; integrity/partition unassigned')],label='separate capture-probe molecules')
 return sc

def capture_state(n):
 sc=post_pcr_state(capture_round=n)
 sc.strand('hybridizing probes',[role('round'+str(n)+' biotinylated120nt viewpoint probes; exact sequence/site mapping unassigned')],label='separate hybridized capture oligos; no covalent library join')
 return sc

ARTICLE = ("Oudelaar et al., Nature Genetics 50:1744 (2018), doi:10.1038/s41588-018-0253-2, "
           "Online Methods 'Tri-C - Experimental procedure'")
PROTEX = ("Oudelaar et al. 2018 Tri-C Protocol Exchange v1, doi:10.1038/protex.2018.113, "
          "'Sonication and addition of Illumina sequencing adaptors' and 'Capture enrichment' "
          "pp. 5-8, Tables 3 and 5 pp. 10-11 (cited as defining-article reference 49)")


def _tri_c_graph(definition, final: Construct,
                 primers: tuple[sp.SeqPrimer, ...]) -> ProtocolGraph:
    """Tri-C graph from the defining Methods and the directly cited detailed method.

    The defining article says 3C libraries were "prepared as previously described"
    (its ref. 41, NG Capture-C) and capture used NimbleGen SeqCap EZ reagents "as
    previously described". Steps whose detail comes only from the cited 2018
    Protocol Exchange method carry that source explicitly; nothing is taken from
    the predecessor's own text.
    """
    source = f"{definition.citation}; doi:{definition.doi}"
    graph = ProtocolGraph(definition.title, authority=source)
    fixed = graph.input("fixed_chromatin", "cells fixed in 2% formaldehyde for 10 min, "
                        "lysed, Dounce-homogenized, SDS/Triton X-100 treated",
                        MaterialKind.BIOLOGICAL,
                        state=fixed_state(),
                        primary=True, source=ARTICLE)
    restriction = graph.input("primary_restriction_enzyme",
                              "NlaIII, 3 x 600 U over 16-24 h at 37 °C", MaterialKind.CATALYST,
                              source=ARTICLE)
    digested = graph.product(
        "restriction_fragments", "NlaIII-digested fixed chromatin", MaterialKind.MOLECULAR,
        state=fragment_state(),
        source=ARTICLE)
    graph.react("primary_digest", "Digest fixed chromatin with NlaIII", StepKind.DIGESTION,
                ((fixed, InputRole.TEMPLATE), (restriction, InputRole.CATALYST)),
                (digested,), source=ARTICLE)
    heat = graph.input("digest_heat_inactivation", "65 °C for 20 min, then ice",
                       MaterialKind.CONDITION, source=ARTICLE)
    inactivated = graph.product("inactivated_digest", "heat-inactivated NlaIII digest",
                                MaterialKind.MOLECULAR, state=digested.state, source=ARTICLE)
    graph.react("inactivate_digest", "Heat-inactivate the digestion", StepKind.THERMAL_TREATMENT,
                ((digested, InputRole.TEMPLATE), (heat, InputRole.CONDITION)),
                (inactivated,), source=ARTICLE, changes_molecule=False)
    ligase = graph.input("proximity_ligase", "T4 DNA ligase, 720 U, overnight at 16 °C",
                         MaterialKind.CATALYST, source=ARTICLE)
    contact_model = multiway_contact(
        "NlaIII", "viewpoint fragment", ("reporter fragment 1", "reporter fragment 2"))
    contact_state = contact_model.scene("multiway 3C concatemer")
    contact = graph.product("contact_product", "multiway proximity-ligation concatemer",
                            MaterialKind.MOLECULAR, state=contact_state, source=ARTICLE)
    graph.react("proximity_ligation", "Proximity-ligate NlaIII fragments in nuclei",
                StepKind.PROXIMITY_LIGATION,
                ((inactivated, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)),
                (contact,), source=ARTICLE)
    proteinase = graph.input("crosslink_reversal",
                             "de-crosslinking with Proteinase K (3 U) overnight at 65 °C",
                             MaterialKind.CONDITION, source=ARTICLE)
    deproteinized = graph.product("deproteinized_3c_dna", "de-crosslinked 3C DNA",
                                  MaterialKind.MOLECULAR, state=contact_state, source=ARTICLE)
    graph.react("reverse_crosslinks", "Reverse cross-links and digest protein",
                StepKind.DEPROTEINIZATION,
                ((contact, InputRole.TEMPLATE), (proteinase, InputRole.CONDITION)),
                (deproteinized,), source=ARTICLE)
    extraction = graph.input("dna_extraction",
                             "RNase (7.5 mU), phenol-chloroform, sodium acetate/ethanol "
                             "precipitation at -80 °C, 70% ethanol wash",
                             MaterialKind.CONDITION, source=ARTICLE)
    library_3c = graph.product("purified_3c_library", "purified 3C library DNA",
                               MaterialKind.MOLECULAR, state=contact_state, source=ARTICLE)
    graph.react("extract_3c_dna", "RNase-treat and extract 3C DNA", StepKind.EXTRACTION,
                ((deproteinized, InputRole.TEMPLATE), (extraction, InputRole.CONDITION)),
                (library_3c,), source=ARTICLE, changes_molecule=False)

    sonication = graph.input("sonication",
                             "Covaris S220, one 55-s cycle, 10% duty, intensity 4, "
                             "200 cycles/burst, target ~450 bp", MaterialKind.CONDITION,
                             source=ARTICLE)
    sonicated_state=duplex([ENDPOINT,*contact_parts()],'sonicated multiway DNA fragment; endpoint coordinate identifier')
    sonicated = graph.product("sonicated_3c_dna", "~450-bp sonicated 3C DNA",
                              MaterialKind.MOLECULAR, state=sonicated_state, source=ARTICLE)
    graph.react("sonicate", "Sonicate gently to retain multiple junctions",
                StepKind.FRAGMENTATION,
                ((library_3c, InputRole.TEMPLATE), (sonication, InputRole.CONDITION)),
                (sonicated,), source=ARTICLE)
    beads_sonication = graph.input("sonication_ampure", "AMPure XP, 0.7:1 bead-to-sample",
                                   MaterialKind.SELECTION, source=ARTICLE)
    sized = graph.product("size_selected_3c_dna", "bead-recovered sonicated 3C DNA",
                          MaterialKind.MOLECULAR, state=sonicated_state, source=ARTICLE)
    sized_waste = graph.product("sonication_cleanup_supernatant",
                                "bead supernatant (fragments and solutes not bound at this ratio)",
                                MaterialKind.WASTE, source=ARTICLE)
    graph.react("sonication_cleanup", "Clean up and size-select after sonication",
                StepKind.CLEANUP,
                ((sonicated, InputRole.TEMPLATE), (beads_sonication, InputRole.SELECTION)),
                (sized, sized_waste), source=ARTICLE, changes_molecule=False)
    graph.discard(sized_waste, "supernatant removed from the magnetically held AMPure beads")
    graph.annotate(SourceDiscrepancy(
        "sonication_cleanup", "Post-sonication bead ratio",
        (SourceClaim("Defining article Online Methods, Tri-C experimental procedure",
                     "Clean-up and size selection after sonication used AMPure XP at 0.7:1."),
         SourceClaim("2018 Protocol Exchange method, Sonication step I.3 and note E (p. 5)",
                     "220 µl beads are added to the 120-µl sonicate (about 1.8x); "
                     "all clean-ups use 1.8 volumes and size selection is called unnecessary.")),
        "The graph uses the article's 0.7:1 size-selecting clean-up for the published experiment.",
        "The defining article describes the libraries actually sequenced; the detailed method "
        "is a general protocol.",
        "Which ratio each sequenced library received is not stated beyond the article's sentence."))

    repair = graph.input("end_repair_enzymes",
                         "NEBNext End Repair buffer and enzyme mix, 20 °C 30 min",
                         MaterialKind.CATALYST, source=PROTEX)
    repaired = graph.product("end_repaired_dna", "end-repaired sonicated 3C DNA",
                             MaterialKind.MOLECULAR, state=sonicated_state, source=PROTEX)
    graph.react("end_repair", "End-repair sonicated 3C DNA", StepKind.END_REPAIR,
                ((sized, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (repaired,), source=PROTEX)
    tailing = graph.input("a_tailing_enzyme", "NEBNext dA-tailing buffer and Klenow, 37 °C 30 min",
                          MaterialKind.CATALYST, source=PROTEX)
    tailed = graph.product("da_tailed_dna", "dA-tailed 3C DNA", MaterialKind.MOLECULAR,
                           state=da_tailed("dA-tailed 3C fragment", "viewpoint-contact fragment"),
                           source=PROTEX)
    graph.react("da_tail", "Add one 3-prime dA to each end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (tailing, InputRole.CATALYST)),
                (tailed,), source=PROTEX)
    adapters = graph.input("illumina_adapters",
                           "Illumina TruSeq-compatible adaptor supplied with the NEBNext kit",
                           MaterialKind.ADAPTER, inferred=True, source=PROTEX)
    adapter_ligase = graph.input("adapter_ligase",
                                 "Quick Ligation buffer and T4 DNA ligase, 20 °C 15 min",
                                 MaterialKind.CATALYST, source=PROTEX)
    ligated = graph.product("adapter_ligated_dna", "adapter-ligated 3C DNA",
                            MaterialKind.MOLECULAR, state=adapted_state(), inferred=False, source=PROTEX)
    graph.react("adapter_ligation", "Ligate Illumina adaptors", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (adapter_ligase, InputRole.CATALYST)), (ligated,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The article names Illumina TruSeq adaptors added with NEBNext reagents and "
                    "the detailed method names the kit 'Adaptor'; neither prints the adaptor "
                    "strands, their uracil position or the ligation-ready end."),
                source=PROTEX)
    user = graph.input("user_enzyme", "USER Enzyme Mix (3 µl), 37 °C 15 min",
                       MaterialKind.CATALYST, source=PROTEX)
    user_treated = graph.product("user_treated_dna", "USER-treated adapter-ligated 3C DNA",
                                 MaterialKind.MOLECULAR, state=adapted_state(opened=True),
                                 source=PROTEX)
    user_products=graph.product('user_products','USER cleavage products; physical partition unassigned',MaterialKind.MOLECULAR)
    graph.react("user_treatment", "Treat the ligation with USER enzyme", StepKind.DIGESTION,
                ((ligated, InputRole.TEMPLATE), (user, InputRole.CATALYST)),
                (user_treated,user_products), source=PROTEX,
                reason="The detailed method adds USER directly to the ligation. Its molecular "
                       "target (a uracil in the kit adaptor) is not printed in the inspected sources.")
    graph.mark_terminal(Terminal(user_products.id,Disposition.UNRESOLVED,reason='Chemical USER cleavage does not determine physical recovery of its fragments.'))
    beads_ligation = graph.input("ligation_ampure", "AMPure XP, 1.8:1 bead-to-sample",
                                 MaterialKind.SELECTION, source=ARTICLE)
    cleaned = graph.product("cleaned_adapter_ligated_dna", "bead-recovered adapter-ligated DNA",
                            MaterialKind.MOLECULAR, state=adapted_state(opened=True), source=ARTICLE)
    cleaned_waste = graph.product("ligation_cleanup_supernatant",
                                  "bead supernatant after ligation/USER",
                                  MaterialKind.WASTE, source=ARTICLE)
    graph.react("ligation_cleanup", "Bead clean-up after ligation and USER treatment",
                StepKind.CLEANUP,
                ((user_treated, InputRole.TEMPLATE), (beads_ligation, InputRole.SELECTION)),
                (cleaned, cleaned_waste), source=ARTICLE, changes_molecule=False)
    graph.discard(cleaned_waste, "supernatant removed from the magnetically held AMPure beads")

    polymerase = graph.input("pcr_polymerase", "Agilent Herculase II Fusion DNA polymerase",
                             MaterialKind.CATALYST, source=ARTICLE)
    forward = graph.input("index_forward_primer", "NEB Universal primer", MaterialKind.OLIGO,
                          source=PROTEX)
    reverse = graph.input("index_reverse_primer", "NEB Index primer (E7335/E7500 sets)",
                          MaterialKind.OLIGO, source=PROTEX)
    indexed = graph.product("indexed_library", "indexed adapter-ligated 3C library",
                            MaterialKind.MOLECULAR, state=post_pcr_state(), source=ARTICLE)
    graph.react("index_pcr", "Seven-cycle indexing PCR", StepKind.PCR,
                ((cleaned, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (indexed,),
                geometry=PCRGeometry.unavailable(
                    cleaned.state, reason="Article: 7-cycle Herculase II PCR; Protocol Exchange Table 3 "
                                  "names NEB Universal and Index primers without bases or the "
                                  "index chosen for each library."),
                source=ARTICLE)
    others = graph.input("other_indexed_libraries",
                         "other differentially indexed Tri-C libraries",
                         MaterialKind.MOLECULAR, source=ARTICLE)
    pooled = graph.product("multiplexed_library", "multiplexed indexed libraries, "
                           "1.5 µg per library", MaterialKind.MOLECULAR, state=post_pcr_state(),
                           source=ARTICLE)
    graph.react("pool_libraries", "Pool libraries for one multiplexed capture", StepKind.POOL,
                ((indexed, InputRole.TEMPLATE), (others, InputRole.TEMPLATE)),
                (pooled,), source=ARTICLE, changes_molecule=False)

    blockers = graph.input("capture_blockers",
                           "COT DNA plus NimbleGen TS-HE Universal and TS-HE Index blocking oligos",
                           MaterialKind.CONDITION, source=PROTEX)
    beads = graph.input("streptavidin_beads", "M-270 streptavidin beads and SeqCap washes",
                        MaterialKind.SELECTION, source=PROTEX)
    capture_pcr = graph.input("post_capture_pcr_mix",
                              "KAPA HiFi HotStart ReadyMix", MaterialKind.CATALYST, source=PROTEX)
    post_lm_1 = graph.input("post_lm_pcr_oligo_1", "NimbleGen Post-LM-PCR Oligo 1",
                            MaterialKind.OLIGO, source=PROTEX)
    post_lm_2 = graph.input("post_lm_pcr_oligo_2", "NimbleGen Post-LM-PCR Oligo 2",
                            MaterialKind.OLIGO, source=PROTEX)
    current = pooled
    for n, hours in ((1, "64-72 h"), (2, "24 h")):
        probes = graph.input(f"capture_oligo_round{n}",
                             "biotinylated 120-nt viewpoint capture oligos "
                             "(Supplementary Table 5 coordinates)",
                             MaterialKind.OLIGO, source=ARTICLE)
        hybrid = graph.product(f"capture_hybrids_round{n}",
                               f"round-{n} library/probe hybridization mix",
                               MaterialKind.MOLECULAR, state=capture_state(n), source=PROTEX)
        graph.react(f"hybridize_round{n}",
                    f"Denature at 95 °C and hybridize with probes at 47 °C for {hours}",
                    StepKind.HYBRIDIZATION,
                    ((current, InputRole.TEMPLATE), (probes, InputRole.STRUCTURAL),
                     (blockers, InputRole.CONDITION)), (hybrid,), source=PROTEX,
                    changes_molecule=False)
        bound = graph.product(f"captured_round{n}", f"round-{n} bead-bound captured library",
                              MaterialKind.MOLECULAR, state=capture_state(n), source=PROTEX)
        reject = graph.product(f"capture_round{n}_reject", f"round-{n} unbound/washed-off library",
                               MaterialKind.WASTE, source=PROTEX)
        graph.react(f"capture_round{n}_step", f"Round-{n} streptavidin capture and stringent washes",
                    StepKind.AFFINITY_CAPTURE,
                    ((hybrid, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                    (bound, reject), source=PROTEX, changes_molecule=False)
        graph.discard(reject, "removed in bead supernatants and 47 °C stringent washes")
        last = n == 2
        amplified = graph.product(
            "sequencing_library" if last else f"amplified_capture_round{n}",
            definition.title + " captured library" if last else
            f"round-{n} amplified captured library",
            MaterialKind.LIBRARY if last else MaterialKind.MOLECULAR, state=post_pcr_state(capture_round=n),
            source=PROTEX)
        graph.react(f"post_capture_pcr_round{n}",
                    f"Amplify round-{n} capture off the beads (KAPA HiFi, Table 5)",
                    StepKind.PCR,
                    ((bound, InputRole.TEMPLATE), (post_lm_1, InputRole.STRUCTURAL),
                     (post_lm_2, InputRole.STRUCTURAL), (capture_pcr, InputRole.CATALYST)),
                    (amplified,),
                    geometry=PCRGeometry.unavailable(
                        bound.state, reason="Protocol Exchange Table 5 names Post-LM-PCR Oligos 1 and 2 "
                                      "(no bases) and a 98/60/72 °C programme with 'go to step 2, "
                                      "14 cycles'; DNA is amplified on the beads, not eluted."),
                    source=PROTEX)
        current = amplified
    graph.sequence(current, SequencingPlan(library=final,primers=(),limitation='MiSeq PE250 first experiment or NextSeq PE150 second experiment; native primer/index acquisition unassigned.'))
    graph.declare_reviewed_scope([
        "fixed_chromatin", "primary_restriction_enzyme", "primary_digest",
        "digest_heat_inactivation", "inactivate_digest", "proximity_ligase",
        "proximity_ligation", "crosslink_reversal", "reverse_crosslinks", "dna_extraction",
        "extract_3c_dna", "sonication", "sonicate", "sonication_ampure", "sonication_cleanup",
        "end_repair_enzymes", "end_repair", "a_tailing_enzyme", "da_tail",
        "illumina_adapters", "adapter_ligase", "adapter_ligation", "user_enzyme",
        "user_treatment", "ligation_ampure", "ligation_cleanup", "pcr_polymerase",
        "index_forward_primer", "index_reverse_primer", "index_pcr",
        "other_indexed_libraries", "pool_libraries", "capture_blockers",
        "streptavidin_beads", "post_capture_pcr_mix", "post_lm_pcr_oligo_1",
        "post_lm_pcr_oligo_2", "capture_oligo_round1", "hybridize_round1", "capture_round1_step", "post_capture_pcr_round1",
        "capture_oligo_round2", "hybridize_round2", "capture_round2_step", "post_capture_pcr_round2", "sequencing_library",
    ], f"{ARTICLE}; {PROTEX}")
    return graph.require_publishable()


def protocol_graph():
 return _tri_c_graph(DEFINITION,FINAL_LIBRARY,())
SPEC=SimpleNamespace(TITLE=TITLE,NOTES=NOTES,SOURCE=SOURCE,SUMMARY=SUMMARY,CAVEAT=CAVEAT,FINAL_LIBRARY=FINAL_LIBRARY,SEQ_PRIMERS=(),FINAL_CAPTION=FINAL_CAPTION,SEQUENCING_INTRO=SEQUENCING_INTRO,SEQUENCING_ENDING=SEQUENCING_ENDING,READ_SPECS=READ_SPECS,SHOW_READ_CONTENT=False,protocol_graph=protocol_graph)
