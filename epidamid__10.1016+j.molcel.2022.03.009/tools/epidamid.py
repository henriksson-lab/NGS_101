"""EpiDamID / scDam&T-seq molecular model.

Source scope.  Rang et al. 2022 (PMC9153956) prints AdRt/AdRb/AdR_PCR (attributed to
Vogel et al. 2007, "for m6A-PCR"), RandomhexRT and the RNA PCR primers.  Its STAR
Methods state that scDam&T-seq "was performed as previously described in detail
(Markodimitraki et al., 2020)" with halved volumes, and its Key Resources attribute the
"Adapters for DamID2" and the SORT-seq/CEL-Seq2 primers to Markodimitraki et al. 2020
without printing them.  The DamID2 adapter below is therefore the example printed in
that cited protocol (PMC7779467 Reagents, DamID adapter); the well-specific barcode set
used for each EpiDamID plate is not assigned.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct,Scene,Segment,feature,oligo,revcomp
from adapter_ligation import AdapterLigationGeometry
from damid import DamMarkedGATC
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)

# Rang et al. 2022 Key Resources (PMC9153956 JATS text lines 574-593; mmc4 2206-2246)
ADRT="CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA"
ADRB="TCCTCGGCCGCG"
ADR_PCR="GGTCGCGGCCGAGGATC"
RANDOMHEX_RT="GCCTTGGCACCCGAGAATTCCANNNNNN"
RNA_PCR_1="AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"
RNA_PCR_INDEX="CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
RNA_PCR_INDEX_EXAMPLE="CGTGAT"   # boundary inferred: Key Resources says "Follow Illumina design"; not delimited in source

# Markodimitraki et al. 2020 Reagents, "Example of one DamID adapter" (PMC7779467 JATS 284-288)
DAMID2_FORK_TOP="GGTGAT"
DAMID2_CORE="CCGGTAATACGACTCACTATAGGGGTTCAGAGTTCTACAGTCCGACGATC"
DAMID2_BC1,DAMID2_BC2="TGCA","TATG"
DAMID2_TOP=DAMID2_FORK_TOP+DAMID2_CORE+"NNN"+DAMID2_BC1+"NNN"+DAMID2_BC2+"GA"
DAMID2_BOTTOM="TCCATANNNTGCANNNGATCGTCGGACTGTAGAACTCTGAACCCCTATAGTGAGTCGTATTACCGGGAGCTT"
DAMID2_FORK_BOTTOM=DAMID2_BOTTOM[-6:]

MARK=DamMarkedGATC("a reader-domain-bound histone modification")
SAMPLE_INDEX=feature("rna_pcr_index","sample_index","fixed",
                     note="6-nt boundary inferred from the Illumina RNA PCR index-primer layout; the source prints the primer without delimiting the index")

def _part(name):
 return "2" if "part 2" in name else "1"
def _umi(name):
 return Segment(name,"NNN","umi",placeholder=True,
                feature=feature("damid_umi_"+_part(name),"umi","random",group="damid_umi",
                                part="UMI part "+_part(name)))
def _bc(name,bases):
 return Segment(name,bases,"cbc",feature=feature("damid_cell_barcode_"+_part(name),"cell_barcode",
                "fixed",group="damid_cell_barcode",part="barcode part "+_part(name)))

def mark_scene():
 sc=Scene.duplex(list(MARK.locus()),label="living-cell chromatin record");sc.note("top","targeting domain recruits Dam; m6A accumulates at nearby GATC sites");sc.labels("top");return sc

def damid2_adapter_scene():
 top=[Segment("top six-base fork",DAMID2_FORK_TOP),Segment("T7 / RA5 constant core",DAMID2_CORE,"t7"),
      _umi("UMI part 1"),_bc("barcode part 1",DAMID2_BC1),_umi("UMI part 2"),
      _bc("barcode part 2",DAMID2_BC2),Segment("terminal GA","GA")]
 bottom=[Segment("terminal TC","TC"),_bc("barcode part 2 bottom",revcomp(DAMID2_BC2)),
         _umi("UMI part 2 bottom"),_bc("barcode part 1 bottom",revcomp(DAMID2_BC1)),
         _umi("UMI part 1 bottom"),Segment("constant core bottom",revcomp(DAMID2_CORE),"t7"),
         Segment("bottom six-base fork",DAMID2_FORK_BOTTOM)]
 if "".join(s.top for s in bottom)!=DAMID2_BOTTOM:
  raise ValueError("DamID2 bottom oligo assembly does not match the printed example")
 sc=Scene();sc.strand("DamID2 top",top,label="DamID2 top oligo")
 sc.anneal("DamID2 bottom",bottom,to="DamID2 top",pair=("terminal TC","terminal GA"),
           label="DamID2 bottom oligo",mod5="/5Phos/",unpaired=("bottom six-base fork",))
 sc.mark("DamID2 top","top six-base fork","non-complementary fork prevents concatemers")
 sc.ligation_end("DamID2 top","terminal GA",kind="blunt",owner="adapter",text="blunt GA/TC end")
 sc.labels("DamID2 top")
 return sc

def dpni_fragment_scene():
 sc=Scene.duplex([Segment("left DpnI half-site","TC"),
                  Segment("methylated-GATC genomic fragment","X"*30,placeholder=True),
                  Segment("right DpnI half-site","GA")],label="DpnI fragment, blunt half-GATC ends")
 sc.ligation_end("top","left DpnI half-site",kind="blunt",owner="substrate",text="blunt DpnI end")
 sc.labels("top")
 return sc

def adapter_ligation_geometry():
 return AdapterLigationGeometry.known("blunt",damid2_adapter_scene(),dpni_fragment_scene(),
  statement="DpnI cuts methylated GATC bluntly; the cited Markodimitraki 2020 protocol prints both DamID2 adapter strands, with a 5'-phosphorylated bottom oligo and a terminal GA that restores GATC with the fragment TC.")

TAGGED=Construct([Segment("top six-base fork",DAMID2_FORK_TOP),
                  Segment("T7 / RA5 constant core",DAMID2_CORE,"t7"),
                  _umi("UMI part 1"),_bc("barcode part 1",DAMID2_BC1),_umi("UMI part 2"),
                  _bc("barcode part 2",DAMID2_BC2),Segment("adapter half-site GA","GA"),
                  Segment("genomic half-site TC","TC"),
                  Segment("Dam-marked genomic insert","X"*30,placeholder=True),
                  Segment("opposite DpnI end adapter","X"*18,placeholder=True,inferred=True)],
                 name="DamID2-adapter-ligated fragment (top strand)")

def _final_library():
 overlap=RNA_PCR_1[-21:]
 i=DAMID2_TOP.find(overlap)
 if i<0:raise ValueError("RNA PCR primer 1 does not overlap the DamID2 adapter")
 after=DAMID2_TOP[i+len(overlap):DAMID2_TOP.index("NNN")]
 tail=RANDOMHEX_RT[1:-6]   # 3'-terminal 21 nt of the RandomhexRT tail
 if not RNA_PCR_INDEX.endswith(tail):
  raise ValueError("RNA PCR index primer does not end in the RandomhexRT tail")
 anneal=RNA_PCR_INDEX[len("CAAGCAGAAGACGGCATACGAGAT")+len(RNA_PCR_INDEX_EXAMPLE):]
 return Construct([
  Segment("P5 / RNA PCR primer 1",RNA_PCR_1,"p5"),
  Segment("adapter RA5 remainder",after,"r1"),
  _umi("UMI part 1"),_bc("barcode part 1",DAMID2_BC1),_umi("UMI part 2"),
  _bc("barcode part 2",DAMID2_BC2),Segment("restored GATC","GATC"),
  Segment("Dam-marked genomic insert","X"*30,placeholder=True),
  Segment("random-hexamer copy","N"*6,placeholder=True),
  Segment("RandomhexRT / index-primer arm (reverse complement)",revcomp(anneal),"r2"),
  Segment("example index (reverse complement; boundary inferred)",revcomp(RNA_PCR_INDEX_EXAMPLE),"cbc",feature=SAMPLE_INDEX,inferred=True),
  Segment("P7 (reverse complement)",revcomp("CAAGCAGAAGACGGCATACGAGAT"),"p7"),
 ],name="EpiDamID DamID-branch library (derived from printed primers and cited adapter)")
FINAL_LIBRARY=_final_library()
SEQ_PRIMERS=()

TITLE="EpiDamID — histone-targeted adenine recording with scDam&T-seq"
NOTES="01_epidamid.html"
SOURCE='<a href="https://doi.org/10.1016/j.molcel.2022.03.009">Rang and de Luca et al., <i>Molecular Cell</i> (2022)</a>; DamID2 adapter from the cited <a href="https://doi.org/10.1038/s41596-020-0314-8">Markodimitraki et al. (2020)</a> protocol.'
SUMMARY="A histone-binding domain or mintbody recruits Dam in living cells, leaving m6A at nearby GATC sites. In each sorted cell, barcoded CEL-Seq2 primers copy mRNA; DpnI then cuts methylated GATC and barcoded T7 DamID2 adapters are ligated, so both branches are amplified together by T7 transcription and finished as one CEL-Seq2 library."
CAVEAT=("<b>Model correction.</b> AdRt/AdRb/AdR_PCR printed in the EpiDamID Key Resources are the Vogel 2007 adapters “for m6A-PCR”, used to check Dam methylation by gel during clone screening; an earlier version of this page wrongly drew them as the scDam&amp;T-seq DamID2 adapter. "
        "The DamID2 adapter is not printed in the EpiDamID paper; the drawn adapter is the single example in the cited Markodimitraki 2020 protocol (barcode TGCA-TATG), and the plate barcode sets are not assigned. "
        "INFERRED — the final DamID library is assembled by sequence overlap of the printed RNA PCR primers, RandomhexRT and that adapter; the source prints no complete final strand or sequencing primers. "
        "The transcript (CEL-Seq2) branch shares the same outer arms but is not drawn base by base here.")
FINAL_CAPTION="DamID-branch library derived from the printed RNA PCR primer 1 (3′ 21 nt matching the adapter RA5 region), the cited DamID2 adapter example and the RandomhexRT tail contained in the RNA PCR index primer."
SEQUENCING_ENDING="Libraries were sequenced on NextSeq 500 (75-base reads) or NextSeq 2000 (100-base reads); reads are demultiplexed by the barcode at the start of R1."
SEQUENCING_UNAVAILABLE="The EpiDamID paper does not name the sequencing primers or paired-end configuration; primer landing is therefore not computed. Read 1 is stated to begin with the UMI/barcode."

def oligos():
 yield oligo("DamID2 adapter top (Markodimitraki 2020 example)",[Segment("DamID2 top",DAMID2_TOP,"t7")])
 yield oligo("DamID2 adapter bottom (Markodimitraki 2020 example)",[Segment("DamID2 bottom",DAMID2_BOTTOM,"t7")],mods="/5Phos/")
 yield oligo("RandomhexRT",[Segment("random-hexamer RT primer",RANDOMHEX_RT,"r2")])
 yield oligo("RNA PCR primer 1",[Segment("RNA PCR 1",RNA_PCR_1,"p5")],mods="terminal G*A phosphorothioate")
 yield oligo("RNA PCR index primer (example)",[Segment("indexed RNA PCR",RNA_PCR_INDEX,"p7")],mods="terminal C*A phosphorothioate")
 yield oligo("AdRt (Vogel 2007; m6A-PCR QC only)",[Segment("AdRt",ADRT)])
 yield oligo("AdRb (Vogel 2007; m6A-PCR QC only)",[Segment("AdRb",ADRB)])
 yield oligo("AdR_PCR (Vogel 2007; m6A-PCR QC only)",[Segment("AdR_PCR",ADR_PCR)])

def protocol_graph():
 source="Defining source: https://doi.org/10.1016/j.molcel.2022.03.009 (STAR Methods, Single-cell Dam&T-seq)"
 cited=source+"; DamID2 adapter example from cited Markodimitraki et al. 2020 (doi:10.1038/s41596-020-0314-8)"
 derived=cited+" (final strand assembled by sequence overlap; not printed)"
 graph=ProtocolGraph("EpiDamID",authority=source)
 chromatin=graph.input("chromatin","sorted single cell: chromatin and mRNA",MaterialKind.BIOLOGICAL,
                       primary=True,state=mark_scene())
 dam=graph.input("targeted_dam","histone-targeted Dam methylase",MaterialKind.CATALYST)
 marked=graph.product("marked_chromatin","m6A-marked chromatin record",MaterialKind.MOLECULAR,
                      state=mark_scene())
 graph.react("m6a_recording","Record nearby GATC sites through targeted Dam methylation",
             StepKind.CHEMICAL_CONVERSION,
             ((chromatin,InputRole.TEMPLATE),(dam,InputRole.CATALYST)),(marked,))
 cel_primer=graph.input("cel_primer","barcoded T7 CEL-Seq2 poly(dT) primer (pre-dispensed)",MaterialKind.OLIGO,source=cited)
 lysis=graph.input("lysis","lysis mix 65 °C 5 min",MaterialKind.CONDITION)
 rt1=graph.input("superscript_ii","SuperScript II, 42 °C 2 h",MaterialKind.CATALYST)
 first=graph.product("first_strand","well with mRNA first-strand cDNA and m6A-marked genomic DNA",MaterialKind.MOLECULAR)
 graph.react("mrna_rt","Lyse and reverse-transcribe mRNA with the barcoded CEL-Seq2 primer",
             StepKind.REVERSE_TRANSCRIPTION,
             ((marked,InputRole.TEMPLATE),(cel_primer,InputRole.STRUCTURAL),
              (lysis,InputRole.CONDITION),(rt1,InputRole.CATALYST)),(first,))
 ss=graph.input("second_strand","E. coli DNA polymerase I, DNA ligase and RNase H, 16 °C 2 h",MaterialKind.CATALYST)
 dscdna=graph.product("ds_cdna","well with double-stranded T7 cDNA and genomic DNA",MaterialKind.MOLECULAR)
 graph.react("second_strand_synthesis","Synthesize second-strand cDNA",StepKind.EXTENSION,
             ((first,InputRole.TEMPLATE),(ss,InputRole.CATALYST)),(dscdna,))
 pk=graph.input("proteinase_k","Proteinase K, 50 °C 10 h then 80 °C",MaterialKind.CATALYST)
 deprot=graph.product("deproteinized","deproteinized well contents",MaterialKind.MOLECULAR)
 graph.react("deproteinize","Digest proteins",StepKind.DEPROTEINIZATION,
             ((dscdna,InputRole.TEMPLATE),(pk,InputRole.CATALYST)),(deprot,))
 dpni=graph.input("dpni","DpnI, 37 °C 6 h then 80 °C",MaterialKind.CATALYST)
 fragments=graph.product("dpni_fragments","DpnI-cut methylated GATC fragments (with cDNA)",
                         MaterialKind.MOLECULAR,state=dpni_fragment_scene())
 graph.react("dpni_digest","Digest methylated GATC with DpnI",StepKind.DIGESTION,
             ((deprot,InputRole.TEMPLATE),(dpni,InputRole.CATALYST)),(fragments,))
 adapter=graph.input("damid2_adapter","barcoded T7 DamID2 adapter (cited 2020 example)",
                     MaterialKind.ADAPTER,state=damid2_adapter_scene(),source=cited)
 ligase=graph.input("ligase","T4 DNA ligase (Roche), 16 °C 16 h",MaterialKind.CATALYST)
 tagged=graph.product("adapter_ligated",TAGGED.name,MaterialKind.MOLECULAR,state=TAGGED,source=cited)
 graph.react("adapter_ligation","Ligate the barcoded DamID2 adapter",StepKind.LIGATION,
             ((fragments,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),
              (ligase,InputRole.CATALYST)),(tagged,),geometry=adapter_ligation_geometry(),source=cited)
 beads=graph.input("pool_beads","pool wells; 0.8 vol PEG/NaCl magnetic beads",MaterialKind.SELECTION)
 pooled=graph.product("pooled","pooled, bead-purified templates",MaterialKind.MOLECULAR)
 bead_waste=graph.product("pool_waste","bead supernatant (small fragments and free adapter fraction)",MaterialKind.WASTE)
 graph.react("pool_cleanup","Pool wells and bead-purify",StepKind.CLEANUP,
             ((tagged,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(pooled,bead_waste))
 graph.discard(bead_waste,"supernatant removed after 0.8-volume bead binding")
 t7=graph.input("t7_polymerase","MEGAscript T7, 37 °C 14 h",MaterialKind.CATALYST)
 rna=graph.product("ivt_rna","amplified RNA from DamID and CEL-Seq2 templates",MaterialKind.MOLECULAR)
 graph.react("ivt","Amplify T7 templates by in-vitro transcription",
             StepKind.TRANSCRIPTION,
             ((pooled,InputRole.TEMPLATE),(t7,InputRole.CATALYST)),(rna,))
 fragmentation=graph.input("fragmentation_condition","0.25 vol fragmentation buffer, 94 °C 90 s; EDTA stop; bead cleanups",
                           MaterialKind.CONDITION)
 fragmented=graph.product("fragmented_rna","fragmented aRNA",MaterialKind.MOLECULAR)
 graph.react("fragmentation","Fragment the bead-purified aRNA",StepKind.FRAGMENTATION,
             ((rna,InputRole.TEMPLATE),(fragmentation,InputRole.CONDITION)),(fragmented,))
 rt_primer=graph.input("rt_primer","RandomhexRT primer",MaterialKind.OLIGO,source=source,state=RANDOMHEX_RT)
 rt=graph.input("reverse_transcriptase","reverse transcriptase (CEL-Seq2 protocol)",MaterialKind.CATALYST)
 cdna=graph.product("cdna","random-primed aRNA cDNA",MaterialKind.MOLECULAR)
 graph.react("reverse_transcription","Random-prime reverse transcription of fragmented aRNA",
             StepKind.REVERSE_TRANSCRIPTION,
             ((fragmented,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),
              (rt,InputRole.CATALYST)),(cdna,))
 p5=graph.input("p5_primer","RNA PCR primer 1",MaterialKind.OLIGO,source=source,state=RNA_PCR_1)
 p7=graph.input("p7_primer","RNA PCR index primer",MaterialKind.OLIGO,source=source,state=RNA_PCR_INDEX)
 polymerase=graph.input("pcr_polymerase","library PCR (8–10 cycles; polymerase per CEL-Seq2)",MaterialKind.CATALYST)
 library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY,
                       source=derived,inferred=True)
 graph.react("library_pcr","PCR-complete the CEL-Seq2-style library",StepKind.PCR,
             ((cdna,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
              (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
             geometry=PCRGeometry.unavailable(
                 FINAL_LIBRARY,reason="The pre-PCR cDNA duplex is not drawn; the final DamID strand is assembled by sequence overlap of printed primers and the cited adapter."),
             source=derived)
 graph.sequence(library,SequencingPlan(
     library=FINAL_LIBRARY,limitation=SEQUENCING_UNAVAILABLE,
     acquired_roles=("Read 1",),read_lengths={"Read 1":75},
     read_uses={"Read 1":"include"}))
 graph.declare_reviewed_scope((
     "chromatin","targeted_dam","m6a_recording","cel_primer","lysis","superscript_ii",
     "mrna_rt","second_strand","second_strand_synthesis","proteinase_k","deproteinize",
     "dpni","dpni_digest","damid2_adapter","ligase","adapter_ligation","pool_beads",
     "pool_cleanup","t7_polymerase","ivt","fragmentation_condition","fragmentation",
     "rt_primer","reverse_transcriptase","reverse_transcription","p5_primer","p7_primer",
     "pcr_polymerase","library_pcr","library"),source)
 return graph.require_publishable()
