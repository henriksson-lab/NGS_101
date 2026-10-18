"""GAM: independent nuclear-profile DNA with source-owned WGA/HT boundaries."""
from types import SimpleNamespace
from chemdraw import Construct, Scene, feature
from protocol_model import molecule, role
from pcr import PCRGeometry
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind,Disposition,Terminal
import seqprimers as sp

TITLE="GAM — Genome Architecture Mapping by nuclear-profile co-segregation"
NOTES="01_gam.html"
SOURCE='<a href="https://doi.org/10.1038/nature21411">Beagrie et al., Nature (2017)</a>; Methods WGA/library preparation and Supplementary Table1.'
SUMMARY="One220nm nuclear profile is microdissected and amplified as its own library. Association is statistical co-segregation across profiles; independent loci are never a covalent contact chimera."
CAVEAT=("Selected route: original TruSeq DNA HT preparation, including Pippin300–500nt selection. The source’s Nano route uses a different350nt bead-selection branch and is not silently merged. WGA4 end adapters occupy30bp; their sequences, native attachment geometry and PCR primer species remain unassigned. The printed single-end100bp run includes30 dark cycles whose relation to the100 is unresolved; no imaged genomic span is inferred. Published positive-control discrepancy: Isolation names E14mESC DNA, while WGA names human genomic DNA in an empty lid. Neither value changes the NP libraries.")


def nuclear_profile_scene(label, *, amplified=False):
 sc=Scene()
 for n in range(1,4):
  parts=[]
  if amplified:parts.append(role("5-prime WGA4 adapter "+str(n),length_bp=30))
  parts.append(role("independent genomic region "+str(n)+"; molecular extent unassigned"))
  if amplified:parts.append(role("3-prime WGA4 adapter "+str(n),length_bp=30))
  sc.strand("region"+str(n),parts,label=label+": region"+str(n))
 return sc


FINAL_LIBRARY=Construct([
 role("TruSeq HT outer end; native geometry unassigned"),
 role("5-prime WGA4 adapter",length_bp=30),
 role("one profile-derived genomic interval; length unassigned"),
 role("3-prime WGA4 adapter",length_bp=30),
 role("opposite TruSeq HT end; native geometry unassigned")],name="GAM source-limited TruSeq HT endpoint")
SEQ_PRIMERS=()
READ_LENGTHS=None
FINAL_CAPTION="One independent profile-derived molecular species; WGA4 adapters are30bp, but native TruSeq ends/index placement remain unassigned. No contact junction exists."
SEQUENCING_INTRO="Methods reports single-end100bp HiSeq rapid-run with the first30bp not imaged. Whether those dark cycles are inside or additional to the100 is unresolved. Table1 gives8-base i7/i5 index sequences for selected NP HT libraries; acquired index cycles remain unassigned."
SEQUENCING_ENDING=SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE="Exact physical sequencing-primer sites and imaged genomic coverage are unassigned."
READ_SPECS=(sp.read_spec("Read 1","100bp as printed;30 dark cycles, relation unresolved","profile genomic sequence after WGA4 end; imaged genomic extent unassigned","include"),)


def oligos():return ()


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 nucleus=g.input("nucleus","fixed sucrose-embedded nucleus",MaterialKind.BIOLOGICAL,state=nuclear_profile_scene("nucleus"),primary=True)
 section=g.input("cryosection","approximately220nm random cryosection",MaterialKind.CONDITION)
 profile=g.product("profile","one thin nuclear profile",MaterialKind.MOLECULAR,state=nuclear_profile_scene("profile"))
 other=g.product("excluded_nucleus","nuclear material outside selected plane",MaterialKind.MOLECULAR)
 g.react("section","Sample a thin nuclear plane",StepKind.SAMPLING,((nucleus,InputRole.TEMPLATE),(section,InputRole.CONDITION)),(profile,other))
 g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"Outside this sampled profile; remaining section/nuclear material may persist and is not assigned a physical disposal."))
 cap=g.input("capture_cap","adhesive laser-microdissection cap",MaterialKind.SELECTION)
 isolated=g.product("isolated","individually captured nuclear profile",MaterialKind.MOLECULAR,state=profile.state)
 uncollected=g.product("uncaptured","section material not collected in this cap",MaterialKind.MOLECULAR)
 g.react("microdissect","Laser-microdissect one profile into its own cap",StepKind.SELECTION,((profile,InputRole.TEMPLATE),(cap,InputRole.SELECTION)),(isolated,uncollected))
 g.mark_terminal(Terminal(uncollected.id,Disposition.UNRESOLVED,"Not collected in this cap; physical disposal of other section material is unassigned."))
 lysis=g.input("wga4_lysis","WGA4 ProteinaseK/lysis-fragmentation buffer;50°C4h,99°C4min",MaterialKind.CATALYST)
 lysed=g.product("lysed","deproteinized profile DNA",MaterialKind.MOLECULAR,state=profile.state)
 g.react("lyse","Lyse the profile in its cap",StepKind.DEPROTEINIZATION,((isolated,InputRole.TEMPLATE),(lysis,InputRole.CATALYST)),(lysed,))
 wga=g.input("wga4_reagent","WGA4 library buffer,stabilisation solution,library-preparation enzyme",MaterialKind.CONDITION)
 prepared=g.product("prepared","WGA4 library-prepared profile DNA",MaterialKind.MOLECULAR,state=nuclear_profile_scene("prepared",amplified=True))
 g.react("prepare","Prepare WGA4 templates at95°C4min then16/24/37°C20min each,75°C5min",StepKind.MAGIC,((lysed,InputRole.TEMPLATE),(wga,InputRole.CONDITION)),(prepared,),reason="The source names WGA4 reagents/times and30bp end adapters; native adapter attachment and oligo geometry are unassigned.")
 f=g.input("wga_forward","WGA4 forward primer role; physical pool/species unassigned",MaterialKind.OLIGO)
 r=g.input("wga_reverse","WGA4 opposing primer role; physical pool/species unassigned",MaterialKind.OLIGO)
 pol=g.input("wga_polymerase","WGA DNA polymerase plus amplification master mix; supplier PCR program",MaterialKind.CATALYST)
 amp=g.product("amplified","WGA-amplified independent profile loci",MaterialKind.MOLECULAR,state=nuclear_profile_scene("amplified",amplified=True))
 g.react("wga","Amplify actual WGA4-prepared profile DNA",StepKind.PCR,((prepared,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(prepared.state,reason="Actual independent WGA4-prepared profile molecules enter PCR; physical primer species/ends are unassigned, not a completed sequencing library."))
 minelute=g.input("wga_cleanup","Qiagen MinElute PCR Purification;elute50µl",MaterialKind.SELECTION)
 clean=g.product("clean_wga","purified WGA profile DNA",MaterialKind.MOLECULAR,state=amp.state)
 wash=g.product("wga_cleanup_waste","MinElute flowthrough/washes",MaterialKind.WASTE)
 g.react("purify_wga","Purify WGA DNA before sequencing-library preparation",StepKind.CLEANUP,((amp,InputRole.TEMPLATE),(minelute,InputRole.SELECTION)),(clean,wash));g.discard(wash,"Excluded by the source column transfer; component partition remains unassigned.")
 kit=g.input("ht_kit","TruSeq DNA HT Sample Prep Kit;whole WGA yield up to1.1µg in55µl",MaterialKind.CONDITION)
 index_inventory=Scene()
 for name in ("i7","i5"):
  index_inventory.strand(name,[role(name+" sample-index oligo inventory",length_bp=8,molecular_feature=feature("gam_"+name,"sample_index","unknown"))],label="source HT NP index inventory; native placement unassigned")
 indexes=g.input("ht_indexes","HT NP indexes;8base i7/i5 sequences in SupplementaryTable1",MaterialKind.OLIGO,state=index_inventory)
 product=g.product("ht_product","TruSeq HT profile library before size selection",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
 g.react("ht_library","Build the selected TruSeq DNA HT profile library",StepKind.MAGIC,((clean,InputRole.TEMPLATE),(kit,InputRole.CONDITION),(indexes,InputRole.STRUCTURAL)),(product,),reason="The source delegates native end processing, ligation and internal PCR to the HT manufacturer protocol; no generic outer shell or terminal PCR on WGA input is substituted.")
 pippin=g.input("pippin","Pippin Prep300–500nt,EtBr-free1.5%agarose cassette",MaterialKind.SELECTION)
 library=g.product("library","size-selected indexed HT profile library",MaterialKind.LIBRARY,state=FINAL_LIBRARY)
 offsize=g.product("off_size","library products outside selected Pippin window",MaterialKind.WASTE)
 g.react("select_library","Select the source HT300–500nt library window",StepKind.SELECTION,((product,InputRole.TEMPLATE),(pippin,InputRole.SELECTION)),(library,offsize));g.discard(offsize,"Not included in the selected Pippin output fraction.")
 g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,acquired_roles=("Read 1",),limitation=SEQUENCING_INTRO))
 g.declare_reviewed_scope(("nucleus","cryosection","section","capture_cap","microdissect","wga4_lysis","lyse","wga4_reagent","prepare","wga_forward","wga_reverse","wga_polymerase","wga","wga_cleanup","purify_wga","ht_kit","ht_indexes","ht_library","pippin","select_library","library"),SOURCE)
 return g.require_publishable()

SPEC=SimpleNamespace(**{k:v for k,v in globals().copy().items() if k.isupper()},oligos=oligos,protocol_graph=protocol_graph)
