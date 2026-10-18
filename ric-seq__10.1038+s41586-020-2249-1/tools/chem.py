"""Native2020 RIC-seq molecular populations, labelled junction and USER precursor."""
from types import SimpleNamespace
from chemdraw import Construct,Scene,Segment,feature,oligo,revcomp,complement_segments
from protocol_model import molecule,role
from adapter_ligation import AdapterLigationGeometry
from endprep import EndPreparedDuplex,dA_tailed_scene
from pcr import PCRGeometry
from protocol_graph import ProtocolGraph,InputRole,MaterialKind,StepKind,SequencingPlan,Terminal,Disposition
import seqprimers as sp
TITLE='RIC-seq — RNA in situ conformation sequencing'
SOURCE='Cai et al., Nature(2020), <a href="https://doi.org/10.1038/s41586-020-2249-1">doi:10.1038/s41586-020-2249-1</a>; original Methods2067–2195, published SupplementaryTable9, GEO GSM3629915.'
NOTES='01_ric_seq.html'
SUMMARY='Formaldehyde-fixed RNP neighbourhoods undergo MNase/FastAP/pCp–biotin/FastAP/PNK preparation, proximity ligation, labelled-RNA enrichment and dUTP strand-specific library preparation.'
CAVEAT='Selected rRNA-depleted2020 route; totalRNA without depletion is a separate source option. Biotin selection retains both chimeras and labelled singletons, not only contacts.2021 PEI adapter sequences are not substituted. Complete2020 Y-adapter/PCR geometry and read/index acquisition remain unassigned; Table9 preserves its7 eight-base and5 six-base index members.'
PE1='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
INDEXES=('AACGTGAT','ACAGGAGC','ACCGAGCT','ACCTTGAA','ACACATCG','ACGCCTAA','ACTGGTCA','CGAAAC','GATCTG','TTGACT','GGAACT','TGACAT')
INDEX_FEATURE=feature('ric_sample','sample_index','whitelist',whitelist='2020Table9 indexprimers1–12;7x8nt/5x6nt;selectedmemberunassigned')
def reverse_primer(index):return Construct([Segment('P7',il.P7),Segment('source index',index,feature=INDEX_FEATURE),Segment('source reverse annealing arm','GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')],name='2020 index primer'+str(INDEXES.index(index)+1))
import illumina as il
FINAL_LIBRARY=Construct([Segment('source PE1 PCR end',PE1),role('RNA-derived DNA; chimera or singleton'),role('opposite2020 PCR end; index member/orientation unassigned',molecular_feature=INDEX_FEATURE)],name='RIC-seq2020 source-limited strand-specific library')
SEQ_PRIMERS=()
READ_LENGTHS=None
SEQUENCING_ENDING='2020 source HiSeqXTen paired-end; biological read lengths and index acquisition remain unassigned.'
SEQUENCING_INTRO=SEQUENCING_ENDING
SEQUENCING_UNAVAILABLE='Native2020 Y-adapter/strand occupancy and exact sequencing-primer sites are not reconstructed.'
READ_SPECS=(sp.read_spec('Read 1',None,'RNA-derived sequence; exact primer geometry unassigned','include'),sp.read_spec('Read 2',None,'opposite end; exact primer geometry unassigned','include'),sp.read_spec('Index 1 (i7)',None,'selected6/8nt physical index family; acquired cycles unassigned','undetermined'))
FINAL_CAPTION='RNA-derived duplex PCR products; source-known PCR oligos do not supply the unknown pre-PCR Y-adapter geometry.'
def oligos():return (oligo('2020Table9 PE1.0',[Segment('PE1.0',PE1)]),*(oligo(x.name,list(x)) for x in map(reverse_primer,INDEXES)))
def rnp_scene(stage='native'):
 sc=Scene()
 for name in ('A','B'):
  parts=[role('independent RNA'+name+' fragment; length unassigned')]
  if stage in ('labelled','cp_oh','kinased'):parts.append(Segment('RNA'+name+' added biotinC','C',note='pCp–biotin cytidine; attachment atom unassigned'))
  sc.strand('RNA'+name,parts,label='fixed RNP RNA'+name,mod5='P' if stage=='kinased' else 'OH/heterogeneous',mod3='P' if stage in ('mnase','labelled') else ('OH' if stage in ('dephos','cp_oh','kinased') else 'native'))
  sc.note('RNA'+name,'Shared formaldehyde-fixed RNP neighbourhood; proximity does not establish RNA–RNA base pairing.')
 return sc

def labelled_products():
 sc=Scene();parts=[role('RNA A fragment'),Segment('internal biotinC at ligation junction','C',note='retained pCp–biotin junction base'),role('RNA B fragment'),Segment('terminal biotinC','C',note='labelled fragmentB terminalC before further fragmentation')]
 sc.strand('chimera',parts,label='proximity-ligated RNA chimera',mod5='P',mod3='OH');sc.junction('chimera',parts[1].name,parts[2].name,'RNA proximity ligation');sc.mark('chimera',parts[1].name,'covalent biotin label; affinity handle')
 sc.strand('singleton',[role('unligated labelled RNA'),Segment('singleton biotinC','C')],label='independent labelled singleton',mod3='OH')
 return sc

def fragmented_products():
 sc=Scene();sc.strand('chimera',[role('fragmented RNA A interval'),Segment('retained internal biotinC','C'),role('fragmented RNA B interval')],label='labelled junction-bearing target; ends variable')
 sc.strand('singleton',[role('fragmented labelled singleton RNA'),Segment('retained singleton biotinC','C')],label='independent terminal-labelled singleton; source may lose terminalC later',mod3='P/OH/cyclicP; heterogeneous')
 return sc

def dna_parts():return [Segment('randomN6 primer','N'*6,placeholder=True,length_bp=6),role('RNA B-derived cDNA interval'),Segment('copied internalC complement',revcomp('C')),role('RNA A-derived cDNA interval')]
def dna_scene(*,duplex=False,adapted=False,surviving=False):
 sc=Scene()
 for name,parts in [('chimera',dna_parts()),('singleton',[Segment('singleton randomN6','N'*6,placeholder=True,length_bp=6),role('singleton RNA-derived cDNA')])]:
  if adapted:parts=[role(name+' nativeY adapter left; geometry unassigned'),*parts,role(name+' nativeY adapter right; geometry unassigned')]
  sc.strand(name,parts,label=name+' first-strand target',mod3='OH')
  if duplex:
   bottom=complement_segments(parts);sc.anneal(name+'_second',bottom,to=name,pair=(bottom[-1].name,parts[0].name),unpaired=tuple(s.name for s in bottom if s.is_role_token()),label='dUTP-marked second-strand population')
   sc.note(name+'_second','80% replacement of dTTP with dUTP; positions/distribution unassigned. NativeY-adapter opposite strands are unassigned.' if adapted else '80%dTTP replacement by dUTP; exact uracil positions unassigned.')
  if surviving:sc.note(name,'Surviving original first-strand target after USER; original dU strand is hydrolysed, not physically discarded here.')
 return sc

def a_tailed_targets():
 # Draw each independent target duplex separately, keeping explicit source dA overhangs.
 from protocol_model import Drawing
 rows=[]
 for label,parts in [('chimera',dna_parts()),('singleton',[role('singleton RNA-derived DNA')])]:
  rows.extend(dA_tailed_scene(EndPreparedDuplex(Construct(parts,name=label),five_prime_phosphate=True,three_prime_overhang='A'),label+' dA-tailed duplex').rows())
 return Drawing(tuple(rows))

def pcp_oligo():
 sc=Scene();sc.strand('pcp',[Segment('biotin-cytidine','C',note='Biotin attachment atom unassigned')],label='pCp–biotin reagent',mod5='P',mod3='P');return sc

def pcp_geometry():
 sc=Scene();sc.strand('RNA',[role('MNase RNA substrate')],label='RNA3′OH',mod3='OH');sc.strand('pcp',[Segment('pCp–biotin','C')],label='5′P cytidine-biotin3′P',mod5='P',mod3='P');sc.ligation_end('pcp','pCp–biotin',kind='single-strand',polarity='5-prime',owner='adapter',text='source5′P joins RNA3′OH; biotinC3′P retained')
 return AdapterLigationGeometry.known('single-strand',sc,statement='2020Methods specify pCp–biotin ligation to dephosphorylated RNA3′OH with T4RNAligase; unknown RNA terminal base does not erase known end chemistry.')

def index_inventory():
 sc=Scene()
 for i,index in enumerate(INDEXES,1):sc.strand('index'+str(i),list(reverse_primer(index)),label='source Table9 index primer'+str(i))
 return sc

def protocol_graph():
    source = SOURCE
    final = FINAL_LIBRARY
    graph = ProtocolGraph(TITLE, authority=source)
    cells = graph.input("cells", "cells containing spatially proximal RNA–protein assemblies",
                        MaterialKind.BIOLOGICAL,
                        state=rnp_scene(), primary=True)
    formaldehyde = graph.input("crosslinker", "1% formaldehyde in PBS, 10 min room temperature; 0.125 M glycine quench "
                               "(Cai 2020 Methods; GEO GSE127188 GSM3629915 treatment protocol)",
                               MaterialKind.CONDITION)
    fixed = graph.product("fixed", "crosslinked RNA–protein neighbourhoods",
                          MaterialKind.MOLECULAR, state=rnp_scene())
    graph.react("crosslink", "Crosslink RNA–protein assemblies in situ",
                StepKind.CROSS_LINKING,
                ((cells, InputRole.TEMPLATE), (formaldehyde, InputRole.CONDITION)), (fixed,))
    permeabilization = graph.input("permeabilization", "permeabilization buffer (NP-40, Triton X-100, "
                                   "Tween 20, SUPERase In), 15 min on ice; PNK-buffer washes "
                                   "(GEO GSE127188 GSM3629915 treatment protocol)", MaterialKind.CONDITION)
    permeable = graph.product("permeabilized", "permeabilized crosslinked cells",
                              MaterialKind.MOLECULAR, state=rnp_scene())
    graph.react("permeabilize", "Permeabilize the crosslinked cells", StepKind.PERMEABILIZATION,
                ((fixed, InputRole.TEMPLATE), (permeabilization, InputRole.CONDITION)), (permeable,))
    fragmentation = graph.input("mnase", "micrococcal nuclease 0.03 U/µl, 10 min 37 °C; stopped by "
                                "EGTA washes (GEO GSE127188 GSM3629915 treatment protocol)",
                                MaterialKind.CATALYST)
    fragments = graph.product("fragments", "fragmented fixed RNA neighbourhoods",
                              MaterialKind.MOLECULAR, state=rnp_scene("mnase"))
    graph.react("fragment", "Fragment RNA in situ with micrococcal nuclease", StepKind.FRAGMENTATION,
                ((permeable, InputRole.TEMPLATE), (fragmentation, InputRole.CATALYST)), (fragments,))
    meth = "Cai 2020 Methods, 'pCp–biotin labelling and proximity ligation'"
    fastap1 = graph.input("fastap_1", f"FastAP alkaline phosphatase 10 U, 37 °C 15 min; PNK+EGTA, high-salt "
                          f"and PNK washes ({meth})", MaterialKind.CATALYST)
    dephos = graph.product("dephosphorylated", "fixed RNA fragments with 3′-OH ends (3′ phosphates removed)",
                           MaterialKind.MOLECULAR, state=rnp_scene("dephos"))
    graph.react("dephosphorylate_3p", "Remove MNase-generated 3′ phosphates", StepKind.DEPHOSPHORYLATION,
                ((fragments, InputRole.TEMPLATE), (fastap1, InputRole.CATALYST)), (dephos,))
    pcp = graph.input("pcp_biotin", "pCp–biotin 4 µl of 1 mM (Thermo Fisher 20160)", MaterialKind.OLIGO, state=pcp_oligo())
    pcp_ligase = graph.input("pcp_ligase", f"T4 RNA ligase 100 U, RNA ligase buffer, RNase inhibitor, "
                             f"15% PEG 20000, overnight 16 °C; PNK washes ({meth})", MaterialKind.CATALYST)
    labelled = graph.product("labelled", "fixed RNA fragments with 3′-terminal pCp–biotin (3′ phosphate)",
                             MaterialKind.MOLECULAR, state=rnp_scene("labelled"))
    graph.react("label_ends", "Ligate pCp–biotin to RNA 3′ ends in situ (omitted in −pCp controls, ED Fig. 1b,c)",
                StepKind.LIGATION,
                ((dephos, InputRole.TEMPLATE), (pcp, InputRole.STRUCTURAL), (pcp_ligase, InputRole.CATALYST)),
                (labelled,), geometry=pcp_geometry())
    fastap2 = graph.input("fastap_2", f"FastAP alkaline phosphatase 10 U, 15 min ({meth})", MaterialKind.CATALYST)
    cp_oh = graph.product("cp_biotin_oh", "fixed RNA fragments ending in Cp–biotin with a 3′-OH",
                          MaterialKind.MOLECULAR, state=rnp_scene("cp_oh"))
    graph.react("dephosphorylate_pcp", "Remove the 3′ phosphate of the added pCp–biotin",
                StepKind.DEPHOSPHORYLATION,
                ((labelled, InputRole.TEMPLATE), (fastap2, InputRole.CATALYST)), (cp_oh,))
    pnk = graph.input("pnk", f"T4 polynucleotide kinase (EK0032) with 1.5 mM ATP in PNK buffer A, 37 °C 45 min; "
                      f"PNK+EGTA and PNK washes ({meth})", MaterialKind.CATALYST)
    kinased = graph.product("kinased", "fixed RNA fragments with 5′-phosphate and 3′ Cp–biotin-OH ends",
                            MaterialKind.MOLECULAR, state=rnp_scene("kinased"))
    graph.react("phosphorylate_5p", "Phosphorylate RNA 5′ ends", StepKind.PHOSPHORYLATION,
                ((cp_oh, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (kinased,))
    ligase = graph.input("proximity_ligase", f"T4 RNA ligase 100 U (Thermo Fisher EL0021), RNA ligase buffer, "
                         f"RNase inhibitor, 0.1 mg/ml BSA in 200 µl, overnight 16 °C ({meth})", MaterialKind.CATALYST)
    chimera = graph.product("chimera", "in-situ proximity-ligated RNA chimera with internal C–biotin at the junction",
                            MaterialKind.MOLECULAR, state=labelled_products())
    graph.react("ligate", "Proximity-ligate labelled RNA ends in situ",
                StepKind.PROXIMITY_LIGATION,
                ((kinased, InputRole.TEMPLATE), (ligase, InputRole.CATALYST)), (chimera,))
    reversal = graph.input("reversal", "proteinase K in 0.5% SDS buffer, 37 °C 60 min then 56 °C 15 min; "
                           "TRIzol LS/chloroform extraction (GEO GSE127188 GSM3629915 extract protocol)",
                           MaterialKind.CONDITION)
    purified = graph.product("purified", "de-crosslinked RNA chimeras",
                             MaterialKind.MOLECULAR, state=labelled_products())
    graph.react("reverse_crosslinks", "Reverse crosslinks and recover RNA chimeras",
                StepKind.DEPROTEINIZATION,
                ((chimera, InputRole.TEMPLATE), (reversal, InputRole.CONDITION)), (purified,))
    dnase = graph.input("dnase", "RQ1 DNase I, 37 °C 20 min; acid phenol:chloroform, ethanol "
                        "precipitation (GEO GSE127188 GSM3629915 extract protocol)", MaterialKind.CATALYST)
    dnased = graph.product("dnase_treated", "DNase-treated total RNA including chimeras",
                           MaterialKind.MOLECULAR, state=labelled_products())
    dna_waste = graph.product("dna_digested", "digested DNA, including RNA–DNA ligation by-products",
                              MaterialKind.MOLECULAR)
    graph.react("dnase_digest", "Digest DNA and RNA–DNA ligation by-products", StepKind.DIGESTION,
                ((purified, InputRole.TEMPLATE), (dnase, InputRole.CATALYST)), (dnased, dna_waste))
    graph.mark_terminal(Terminal(dna_waste.id,Disposition.UNRESOLVED,"DNase produces molecular fragments; their partition through source acid-phenol/ethanol recovery is unresolved, not physical discard by digestion."))
    depletion = graph.input("rrna_depletion", "Ribo-off rRNA Depletion Kit (Vazyme N406-01); optional "
                            "per GEO — GSE127188 holds rRNA_depleted and Total libraries",
                            MaterialKind.SELECTION)
    depleted = graph.product("depleted", "rRNA-depleted total RNA including chimeras",
                             MaterialKind.MOLECULAR, state=labelled_products())
    rrna = graph.product("rrna_removed", "depleted rRNA", MaterialKind.WASTE)
    graph.react("deplete_rrna", "Optionally deplete rRNA (ED Fig. 1a, top; GEO extract protocol)",
                StepKind.SELECTION,
                ((dnased, InputRole.TEMPLATE), (depletion, InputRole.SELECTION)), (depleted, rrna))
    graph.discard(rrna, "removed by rRNA depletion")
    refrag = graph.input("rna_fragmentation", "94 °C for 5 min in 1× first-strand buffer (Mg2+), then ice "
                         "(GEO GSE127188 GSM3629915 extract protocol)", MaterialKind.CONDITION)
    sized = graph.product("sized", "fragmented rRNA-depleted RNA", MaterialKind.MOLECULAR,
                          state=fragmented_products())
    graph.react("fragment_rna", "Fragment the recovered RNA (ED Fig. 1a, middle)",
                StepKind.FRAGMENTATION,
                ((depleted, InputRole.TEMPLATE), (refrag, InputRole.CONDITION)), (sized,))
    beads = graph.input("c1_beads", "MyOne Streptavidin C1 beads in TWB, 30 min room temperature; "
                        "eluted 3× at 95 °C in PK buffer (GEO GSE127188 GSM3629915 extract protocol)",
                        MaterialKind.SELECTION)
    enriched = graph.product("enriched", "pCp–biotin-labelled RNA fragments after C1 pull-down",
                             MaterialKind.MOLECULAR, state=fragmented_products())
    unbound = graph.product("unbound", "unlabelled RNA fragments in C1 flow-through",
                            MaterialKind.WASTE)
    graph.react("c1_pulldown", "Enrich biotin-labelled RNA on C1 beads (ED Fig. 1a, bottom)",
                StepKind.AFFINITY_CAPTURE,
                ((sized, InputRole.TEMPLATE), (beads, InputRole.SELECTION)), (enriched, unbound))
    graph.discard(unbound, "not retained by streptavidin C1 beads")
    rt_primer = graph.input("library_rt_primer", "N6 random primer (Cai 2020 Methods; GEO GSE127188 GSM3629915)",
                            MaterialKind.OLIGO,state=molecule("randomN6 primer",(Segment("N6","N"*6,placeholder=True,length_bp=6),)))
    rt = graph.input("library_rt_enzyme", "SuperScript II reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product("library_cdna", "RIC-seq chimera cDNA", MaterialKind.MOLECULAR,
                         state=dna_scene())
    rt_rna=graph.product("rt_template_rna","RNA templates after RT",MaterialKind.MOLECULAR,state=fragmented_products())
    graph.react("library_rt", "Reverse-transcribe the enriched RNA", StepKind.REVERSE_TRANSCRIPTION,
                ((enriched, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,rt_rna))
    second = graph.input("library_second_strand", "RNase H, E. coli DNA polymerase I and dNTPs with "
                         "80% of dTTP replaced by dUTP, 16 °C 2 h (GEO extract protocol)",
                         MaterialKind.CATALYST)
    ds = graph.product("library_ds_cdna", "double-stranded cDNA with dUTP-marked second strand",
                       MaterialKind.MOLECULAR,
                       state=dna_scene(duplex=True))
    rna_debris=graph.product("rt_rna_debris","RNaseH RNA hydrolysis products",MaterialKind.MOLECULAR)
    graph.react("library_second_strand_synthesis", "Synthesize the dUTP-marked second strand",
                StepKind.EXTENSION, ((cdna, InputRole.TEMPLATE),(rt_rna,InputRole.TEMPLATE), (second, InputRole.CATALYST)), (ds,rna_debris))
    graph.mark_terminal(Terminal(rna_debris.id,Disposition.UNRESOLVED,"RNaseH produces fragments; precise recovery through subsequent bead cleanups is unassigned."))
    ds_clean_rule=graph.input("ds_clean_rule","1.8×AMPure recovery after second-strand synthesis",MaterialKind.SELECTION)
    ds_clean=graph.product("ds_clean","1.8×AMPure recovery after second-strand synthesis targetDNA",MaterialKind.MOLECULAR,state=ds.state)
    ds_clean_waste=graph.product("ds_clean_waste","source bead supernatant/wash fractions",MaterialKind.WASTE)
    graph.react("ds_clean_step","1.8×AMPure recovery after second-strand synthesis",StepKind.CLEANUP,((ds,InputRole.TEMPLATE),(ds_clean_rule,InputRole.SELECTION)),(ds_clean,ds_clean_waste));graph.discard(ds_clean_waste,"Physically excluded source bead supernatant and wash fractions; detailed hydrolysate partition remains unassigned.")
    repair = graph.input("library_end_repair", "T4 DNA polymerase, Klenow DNA polymerase and T4 PNK "
                         "after 1.8× AMPure (GEO extract protocol)", MaterialKind.CATALYST)
    repaired = graph.product("library_repaired", "end-repaired ds cDNA", MaterialKind.MOLECULAR,
                             state=ds.state)
    graph.react("library_end_repair_step", "End-repair the ds cDNA", StepKind.END_REPAIR,
                ((ds_clean, InputRole.TEMPLATE), (repair, InputRole.CATALYST)), (repaired,))
    repair_clean_rule=graph.input("repair_clean_rule","1.8×AMPure recovery after end repair",MaterialKind.SELECTION)
    repair_clean=graph.product("repair_clean","1.8×AMPure recovery after end repair targetDNA",MaterialKind.MOLECULAR,state=repaired.state)
    repair_clean_waste=graph.product("repair_clean_waste","source bead supernatant/wash fractions",MaterialKind.WASTE)
    graph.react("repair_clean_step","1.8×AMPure recovery after end repair",StepKind.CLEANUP,((repaired,InputRole.TEMPLATE),(repair_clean_rule,InputRole.SELECTION)),(repair_clean,repair_clean_waste));graph.discard(repair_clean_waste,"Physically excluded source bead supernatant and wash fractions; detailed hydrolysate partition remains unassigned.")
    klenow = graph.input("library_a_tailing", "Klenow exo− (GEO extract protocol)", MaterialKind.CATALYST)
    tailed = graph.product("library_a_tailed", "3′-A-tailed ds cDNA", MaterialKind.MOLECULAR,
                           state=a_tailed_targets())
    graph.react("library_a_tail", "Add 3′ adenines", StepKind.TAILING,
                ((repair_clean, InputRole.TEMPLATE), (klenow, InputRole.CATALYST)), (tailed,))
    adaptors = graph.input("library_adapters", "Illumina Y-shaped adaptors, 2 pmol (Cai 2020 Methods; "
                           "bases and ends not printed there)", MaterialKind.ADAPTER)
    ligase_dna = graph.input("library_adapter_ligase", "Quick T4 DNA ligase, 20 °C 15 min; AMPure ×2 "
                             "(GEO extract protocol)", MaterialKind.CATALYST)
    adapted = graph.product("library_adapted", "adapter-flanked chimera cDNA", MaterialKind.MOLECULAR,
                            state=dna_scene(duplex=True,adapted=True))
    graph.react("library_adapter_install", "Ligate paired-end adapters to the A-tailed cDNA",
                StepKind.LIGATION, ((tailed, InputRole.TEMPLATE), (adaptors, InputRole.STRUCTURAL),
                                    (ligase_dna, InputRole.CATALYST)),
                (adapted,), geometry=AdapterLigationGeometry.unavailable(
                    "The Cai 2020 Methods name 'Illumina Y-shaped adaptors (2 pmol)' after "
                    "Klenow exo- A-tailing but print no adaptor bases or ends. The PEI adapter A/B "
                    "printed in Cao 2021 Nat Protoc Table 1 belong to that later protocol and are "
                    "not assumed for the 2020 libraries."))
    adapted_clean_rule=graph.input("adapted_clean_rule","two consecutiveAMPure recoveries after adapter ligation",MaterialKind.SELECTION)
    adapted_clean=graph.product("adapted_clean","two consecutiveAMPure recoveries after adapter ligation targetDNA",MaterialKind.MOLECULAR,state=adapted.state)
    adapted_clean_waste=graph.product("adapted_clean_waste","source bead supernatant/wash fractions",MaterialKind.WASTE)
    graph.react("adapted_clean_step","two consecutiveAMPure recoveries after adapter ligation",StepKind.CLEANUP,((adapted,InputRole.TEMPLATE),(adapted_clean_rule,InputRole.SELECTION)),(adapted_clean,adapted_clean_waste));graph.discard(adapted_clean_waste,"Physically excluded source bead supernatant and wash fractions; detailed hydrolysate partition remains unassigned.")
    user = graph.input("library_user", "USER enzyme, 37 °C 15 min in the PCR mix to digest the dUTP "
                       "strand (GEO extract protocol)", MaterialKind.CATALYST)
    forward = graph.input("library_forward", "Illumina PE1.0 primer (Supplementary Table 9)",
                          MaterialKind.OLIGO,state=PE1)
    reverse = graph.input("library_reverse", "Index primer 1–12 (Supplementary Table 9)",
                          MaterialKind.OLIGO,state=index_inventory())
    polymerase = graph.input("library_polymerase", "Platinum Pfx DNA polymerase, 11–13 cycles (GEO extract protocol)", MaterialKind.CATALYST)
    surviving=graph.product("user_surviving","surviving first-strand targets after USER",MaterialKind.MOLECULAR,state=dna_scene(adapted=True,surviving=True))
    user_debris=graph.product("user_debris","degraded uracil-strand products",MaterialKind.MOLECULAR)
    graph.react("user_digest","USER phase within PCR incubation: digest dU strand37°C15min",StepKind.DIGESTION,((adapted_clean,InputRole.TEMPLATE),(user,InputRole.CATALYST)),(surviving,user_debris))
    graph.mark_terminal(Terminal(user_debris.id,Disposition.UNRESOLVED,"Hydrolysed dU-strand products remain molecular; subsequent gel recovery is not mapped molecule by molecule."))
    amplified=graph.product("library_amplified","PCR-amplified RIC-seq library",MaterialKind.MOLECULAR,state=final)
    graph.react("library_pcr","Amplify surviving first-strand targets11–13cycles",StepKind.PCR,((surviving,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(surviving.state,reason="Actual postUSER first-strand target population, not finished library.2020Table9 primers are known, but native Y-end/first-copy geometry remains unassigned; later2021 adapters are not substituted."))
    gel = graph.input("gel", "agarose gel, 200–450 bp products, MinElute Gel Extraction (Qiagen 28606) "
                      "(Cai 2020 Methods; ED Fig. 1b,c)",
                      MaterialKind.SELECTION)
    library = graph.product("library_library", "RIC-seq sequencing library", MaterialKind.LIBRARY,
                            state=final)
    off_size = graph.product("gel_excluded", "library fragments outside the excised gel region",
                             MaterialKind.WASTE)
    graph.react("gel_purify", "Gel-purify the library (ED Fig. 1b,c)", StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (gel, InputRole.SELECTION)), (library, off_size))
    graph.discard(off_size, "outside excised gel region")
    graph.sequence(library,SequencingPlan(library=final,acquired_roles=('Read 1','Read 2'),limitation=SEQUENCING_ENDING+' '+SEQUENCING_UNAVAILABLE))
    graph.declare_reviewed_scope(('cells','crosslinker','crosslink','permeabilization','permeabilize','mnase','fragment','fastap_1','dephosphorylate_3p','pcp_biotin','pcp_ligase','label_ends','fastap_2','dephosphorylate_pcp','pnk','phosphorylate_5p','proximity_ligase','ligate','reversal','reverse_crosslinks','dnase','dnase_digest','rrna_depletion','deplete_rrna','rna_fragmentation','fragment_rna','c1_beads','c1_pulldown','library_rt_primer','library_rt_enzyme','library_rt','library_second_strand','library_second_strand_synthesis','ds_clean_rule','ds_clean_step','library_end_repair','library_end_repair_step','repair_clean_rule','repair_clean_step','library_a_tailing','library_a_tail','library_adapters','library_adapter_ligase','library_adapter_install','adapted_clean_rule','adapted_clean_step','library_user','user_digest','library_forward','library_reverse','library_polymerase','library_pcr','gel','gel_purify','library_library'),SOURCE)
    return graph.require_publishable()

SPEC=SimpleNamespace(**{k:v for k,v in globals().copy().items() if k.isupper()},oligos=oligos,protocol_graph=protocol_graph)
