"""Protocol-owned native spatial chemistry; molecular populations and source limits explicit."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from dataclasses import replace
from chemdraw import Construct,Scene,Segment,feature,oligo,revcomp,complement_segments
from protocol_model import role
from protocol_graph import ProtocolGraph,MaterialKind,InputRole,StepKind,SequencingPlan,Terminal,Disposition,SourceClaim,SourceDiscrepancy
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
import seqprimers as sp

def add_duplex(sc,key,parts,label):
 sc.strand(key,parts,label=label);bottom=complement_segments(parts)
 sc.anneal(key+" complement",bottom,to=key,pair=(parts[0].name+"'",parts[0].name),unpaired=tuple(p.name for p in bottom if p.is_role_token()),label=label)
def rc_parts(parts):
 return [role(p.name+" copied complement",molecular_feature=p.feature,inferred=p.inferred) if p.is_role_token() else replace(p,name=p.name+" copied complement",top=revcomp(p.top)) for p in reversed(parts)]
def tools(g):
 def inp(i,label,kind=MaterialKind.CATALYST,state=None):return g.input(i,label,kind,state=state)
 def out(i,label,state=None,kind=MaterialKind.MOLECULAR):return g.product(i,label,kind,state=state)
 def op(i,label,kind,inputs,state,geometry=None,reason=""):
  product=out(i+"_out",label,state);g.react(i,label,kind,inputs,(product,),geometry=geometry,**({"reason":reason} if reason else {}));return product
 def clean(i,label,x,state=None):
  selection=inp(i+"_selection",label,MaterialKind.SELECTION);kept=out(i+"_out",label,x.state if state is None else state,x.kind);waste=out(i+"_waste","physically removed wash/flowthrough fraction",None,MaterialKind.WASTE)
  g.react(i,label,StepKind.CLEANUP,((x,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(kept,waste));g.discard(waste,"Source physical recovery fraction; exact molecular partition/purity unassigned.");return kept
 return inp,out,op,clean
TRUSEQ_PCR_HANDLE = "CTACACGACGCTCTTCCGATCT"
BEAD_LINKER = "TCTTCAGCGTTCCCGAGA"
DECODING_PRIMERS = (
    ("Truseq", "AGATCGGAAGAGCGTCGTGTAG"), ("Truseq -1", "GATCGGAAGAGCGTCGTCTAG"),
    ("Truseq -2", "ATCGGAAGAGCGTCGTGTAG"), ("TruSeq-3", "TCGGAAGAGCGTCGTGTAG"),
    ("TruSeq-4", "CGGAAGAGCGTCGTGTAG"), ("UP", "TCTCGGGAACGCTGAAGA"),
    ("UP-1", "CTCGGGAACGCTGAAGA"), ("UP-2", "TCGGGAACGCTGAAGA"),
    ("UP-3", "CGGGAACGCTGAAGA"), ("UP-4", "GGGAACGCTGAAGA"))


def _recessed_series_mismatches():
    """Printed decoding primers that are not the 5'-recessed form of their parent."""
    parents = {"Truseq": revcomp(TRUSEQ_PCR_HANDLE), "UP": revcomp(BEAD_LINKER)}
    out = []
    for name, seq in DECODING_PRIMERS:
        family = "UP" if name.startswith("UP") else "Truseq"
        n = int(name.rsplit("-", 1)[1]) if "-" in name else 0
        expected = parents[family][n:]
        if seq != expected:
            out.append((name, seq, expected))
    return tuple(out)


DECODING_MISMATCHES = _recessed_series_mismatches()

TITLE="Slide-seq v1 — native Sequence2 bead chemistry"
NOTES="01_slide_seq.html"
SOURCE='<a href="https://doi.org/10.1126/science.aaw1219">Rodriques et al., Science2019</a>; PMC6927209 original author-manuscript supplement.'
SUMMARY="Sequence2beads retain their source5′linker/T8 prefix through RT and bead recovery. Source-known SOLiD decoding, template switching, bead-bound preamplification and NexteraXT tagmentation precede native mixed-handle library PCR."
CAVEAT="Selected native branch is source BeadSequence2 (J8+J7,UMI8,T30). Sequence1 is an independently manufactured PC/T7/J6+J7 bead variant, not an aliquot physically split from Sequence2 and not a reason to invent photorelease here. Optical decoding reports14bases whereas Sequence2physicalbarcode has15; exact optical coverage is unassigned. The published author-manuscript supplement is readable; identity with the gated Science-hostedPDF remains unchecked. Printed TSO/SMART and decoderTruseq-1 discrepancies stay public. Native NexteraXT adapter/primer-member geometry and molecular N700index extent remain unassigned;8acquiredindexcycles do not supply those strands."
BARCODE=feature("slide.bead","spatial_barcode","whitelist",group="slide.bead",whitelist="decoded Slide-seq beadpool")
UMI=feature("slide.umi","umi","random")
INDEX=feature("slide.index","sample_index","unknown",note="N700series member/physical index extent unassigned in capped pass")
SMART="AAGCAGTGGTATCAACGCAGAGT"
TSO_DNA="AAGCTGGTATCAACGCAGAGTGAAT"
TSO_WRITTEN=TSO_DNA+"rG+GrG"
from illumina import TRUSEQ_P5_FULL
TRUSEQ5=TRUSEQ_P5_FULL
def tag_parts():return [Segment("native TruSeq handle",TRUSEQ_PCR_HANDLE,"r1"),Segment("barcode block1","J"*8,"cbc",placeholder=True,feature=BARCODE),Segment("source bead linker",BEAD_LINKER),Segment("barcode block2","J"*7,"cbc",placeholder=True,feature=BARCODE),Segment("UMI","N"*8,"umi",placeholder=True,feature=UMI),Segment("poly(dT)30","T"*30)]
def bead_parts():return [Segment("source T8prefix","T"*8),*tag_parts()]
def bead_scene(extended=False):
 sc=Scene();parts=bead_parts()+([role("first-strand transcript cDNA"),role("native template-switch C/G junction; boundary unassigned"),Segment("copied source TSO DNAprefix",revcomp(TSO_DNA),inferred=True)] if extended else [])
 sc.strand("native bead DNA",parts,mod5="native bead linker; chemistry unassigned",label="actual bead-bound ssDNA; full5′T8 retained until PCR")
 if extended:sc.strand("original RNA",[role("independent captured RNA template; physical U")],label="RNA not consumed by RT")
 return sc
def tso_scene():
 sc=Scene();sc.strand("TSO",[Segment("source TSO DNAprefix",TSO_DNA),Segment("terminal RNA G1","G"),Segment("terminal LNA G2","G"),Segment("terminal RNA G3","G")],label="source-written rG+GrG:2RNA G and middle LNA G; no prefix repair")
 return sc
def captured_scene():
 sc=bead_scene();sc.strand("RNA",[role("independent tissue poly(A)RNA; physical U")],label="hybridized native RNA; not covalently joined to DNA")
 return sc
def decoding_scene(stage):
 sc=bead_scene();sc.strand("decoder strand",[role("source5′phosphorylated decoding primer with "+stage)],label="independent primer/probe strand, not joined to bead DNA")
 return sc
def preamp_parts():return [*tag_parts(),role("transcript cDNA"),role("native SMART-primer-derived opposite end; printedTSOfootprint unresolved",inferred=True)]
def preamp_scene():
 sc=Scene();add_duplex(sc,"expected preamp",preamp_parts(),"expected two-primer bounded preamplified DNA")
 sc.strand("prior populations",[role("possible bead-linked originals, RNA, TSO/decoder/primer products and other copies; recovery unresolved")],label="independent original/copy populations")
 return sc
def tagged_scene():
 sc=Scene();sc.strand("barcoded tagged fragment",[*tag_parts(),role("native transcript fragment"),role("native NexteraXTn7-tagged opposite end; gap/strand geometry unassigned")],label="actual tagged precursor; no completed P5/P7shell")
 sc.strand("other populations",[role("possible prior amplicons, bead-linked originals/RNA/oligos, nonamplifying tagged/offtarget fragments; physical partition unresolved")],label="independent possible co-carried species; nonamplification is not discard")
 return sc
FINAL_LIBRARY=Construct([role("Truseq5-derived native Read1 end"),*tag_parts()[1:],role("native transcript fragment"),role("N700indexed opposite end; native extent unassigned",molecular_feature=INDEX)],name="Slide-seq source-limited native endpoint")
SEQ_PRIMERS=();FINAL_CAPTION="Sequence2physicalJ8/linker18/J7/UMI8/T30; native outer/index ends unassigned."
SEQUENCING_ENDING="Published supplement NovaSeqS2:R1=42,I1=8,R2=50. Barcodephysical15vs opticaldecoded14 stay separate; index member/strand extent/primer entry unassigned."
READ_SPECS=(sp.read_spec("Read 1",42,"Sequence2beadbarcode blocks/linker and UMI; physical41bases before firstT","extract"),sp.read_spec("Index 1 (i7)",8,"N700sample index; native member/extent unassigned","demultiplex"),sp.read_spec("Read 2",50,"transcript","include"),sp.read_spec("Index 2 (i5)",None,"not assigned in inspected run recipe","undetermined"))
def oligos():
 yield oligo("Sequence2native bead DNA",bead_parts(),mods="native5′linker")
 yield "Printed TSO:5′-"+TSO_WRITTEN+"-3′ (2RNA G,middleLNA G; no repair)"
 yield oligo("SMARTPCR",[Segment("SMART",SMART)])
 yield oligo("TruSeqhandlePCR",[Segment("source handle",TRUSEQ_PCR_HANDLE)])
 yield oligo("Truseq5",[Segment("Truseq5",TRUSEQ5)])
 for n,s in DECODING_PRIMERS:yield oligo("source SOLiDdecoder "+n,[Segment("source bases",s)],mods="/5Phos/")
def decoder_scene():
 sc=Scene()
 for i,(name,bases) in enumerate(DECODING_PRIMERS):sc.strand("decoder"+str(i),[Segment("source "+name,bases)],mod5="p",label="independent source5′phosphorylated "+name)
 return sc

def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE);graph=g;inp,out,op,clean=tools(g)
 bead=g.input("bead","Sequence2native bead puck; independently synthesized variant",MaterialKind.BIOLOGICAL,state=bead_scene(),primary=True)
 decoder=inp("decoder","10source5′phosphorylated SOLiDdecoding primers",MaterialKind.OLIGO,decoder_scene())
 dec=op("decoder_hybridize","Hybridize decodingprimer4×SSC;20min",StepKind.HYBRIDIZATION,((bead,InputRole.TEMPLATE),(decoder,InputRole.STRUCTURAL)),decoding_scene("annealed constant handle"))
 probe=inp("solid_probes","SOLiDSR75 sequencing probes; exact probe species unassigned",MaterialKind.OLIGO);lig=inp("decoding_ligase","T4DNAligase;2ligations/primer for10primers")
 geo=AdapterLigationGeometry("unavailable","unavailable",(decoding_scene("annealed constant handle"),),"Native primer/probe ligation is separate from bead DNA; source probe geometry not assigned.")
 dec=op("decode_ligate","Ligate SOLiDprobe onto annealed decoder strand; image",StepKind.LIGATION,((dec,InputRole.TEMPLATE),(probe,InputRole.STRUCTURAL),(lig,InputRole.CATALYST)),decoding_scene("ligated SOLiDprobe"),geometry=geo)
 buffer=inp("cleave_buffers","SOLiDbuffersC/B cleave fluorophore; bufferF washes",MaterialKind.CONDITION)
 dec=op("decode_cleave","Cleave fluorophore; repeat probe ligation for optical decoding",StepKind.CHEMICAL_CONVERSION,((dec,InputRole.TEMPLATE),(buffer,InputRole.CONDITION)),decoding_scene("probe cleavage products; exact residue state unassigned"))
 dec=clean("decode_strip","80%formamide strip10min and instrument-buffer flow; repeat primers",dec,bead_scene())
 rna=inp("tissue_rna","independent tissue RNA;6×SSC15min capture",MaterialKind.MOLECULAR)
 captured=op("capture_rna","Hybridize RNA to decoded native bead poly(dT)",StepKind.HYBRIDIZATION,((dec,InputRole.TEMPLATE),(rna,InputRole.STRUCTURAL)),captured_scene())
 primer=inp("bead_primer","actual native bead capture primer; same molecular population",MaterialKind.OLIGO,Construct(bead_parts()));rt=inp("rt","MaximaH−42°C1h");tso=inp("tso","printed TSO DNA/rG+GrG reagent",MaterialKind.OLIGO,tso_scene())
 cdna=op("rt_switch","INFERRED junction—RT and source TSO copying; no repaired TSO prefix",StepKind.REVERSE_TRANSCRIPTION,((captured,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),bead_scene(True))
 proteinase=inp("proteinase","ProteinaseK/SDS tissue digestion37°C40min")
 cdna=op("digest_tissue","Digest tissue before physical bead recovery",StepKind.DEPROTEINIZATION,((cdna,InputRole.TEMPLATE),(proteinase,InputRole.CATALYST)),cdna.state)
 cdna=clean("recover_beads","Vigorouspipetting;removeglass;pellet/washbeads repeatedly thenwater",cdna)
 f=inp("pre_forward","TruseqhandlePCRprimer",MaterialKind.OLIGO,TRUSEQ_PCR_HANDLE);r=inp("pre_reverse","sourceSMARTPCRprimer",MaterialKind.OLIGO,SMART);p=inp("pre_pol","KAPAHiFiHotStart;4cycles65°C+9cycles67°C")
 amp=out("preamp","native preamplified DNA plus possible originals",preamp_scene())
 g.react("pre_pcr","INFERRED SMART junction—amplify actual bead-bound ssDNA",StepKind.PCR,((cdna,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(p,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(cdna.state,reason="Actual full-prefix beadssDNA provided; source printedTSO/SMART discrepancy leaves opposite-end/firstcycle geometry unresolved."))
 amp=clean("pre_cleanup","0.6×AMPure;10µlrecovery;600pgDNA selected forNexteraXT",amp)
 transposome=inp("transposome","native NexteraXT loadedTn5 adapter; full strands unassigned",MaterialKind.ADAPTER);tn5=inp("tn5","NexteraXT tagmentation mix; native member/geometry unassigned")
 tagged=op("tagment","Tagment actual preamplified DNA with NexteraXT",StepKind.TAGMENTATION,((amp,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),tagged_scene())
 boundary=inp("kit_conditions","native NexteraXT post-tagmentation handoff; detailed recovery not read",MaterialKind.CONDITION)
 tagged=op("native_handoff","Bounded native post-tagmentation kit handoff before indexingPCR",StepKind.MAGIC,((tagged,InputRole.TEMPLATE),(boundary,InputRole.CONDITION)),tagged.state,reason="Supplement delegates NexteraXT tagmentation to manufacturer; intervening native neutralization/recovery detail remains unassigned within cap.")
 f=inp("library_forward","sourceTruseq5",MaterialKind.OLIGO,TRUSEQ5);r=inp("library_reverse","nativeN700indexprimer; exactmember/fullstrand unassigned",MaterialKind.OLIGO);p=inp("library_pol","native indexingPCR;72°C3min then12cycles55°Cannealing")
 tagged=op("initial_extension","Native72°C3min pre-cycle extension; completion geometry unassigned",StepKind.EXTENSION,((tagged,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(p,InputRole.CATALYST)),tagged_scene())
 final=tagged_scene();add_duplex(final,"expected library",list(FINAL_LIBRARY),"expected native daughter; prior populations retained")
 library=out("library","Slide-seq native indexed library",final,MaterialKind.LIBRARY)
 g.react("library_pcr","Amplify actual native tagged PCR precursor",StepKind.PCR,((tagged,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(p,InputRole.CATALYST)),(library,),geometry=PCRGeometry.unavailable(tagged.state,reason="Actual native tagged population supplied; N700member and early-cycle adapter geometry remain unassigned, no finished library used as template."))
 library=clean("final_cleanup","0.6×AMPure;10µlrecovery;4nMnormalization",library)
 g.sequence(library,SequencingPlan(acquired_roles=("Read 1","Index 1 (i7)","Read 2"),read_lengths={"Read 1":42,"Index 1 (i7)":8,"Read 2":50},limitation=SEQUENCING_ENDING))
 g.declare_reviewed_scope(("bead","decoder","decoder_hybridize","solid_probes","decoding_ligase","decode_ligate","cleave_buffers","decode_cleave","decode_strip_selection","decode_strip","tissue_rna","capture_rna","bead_primer","rt","tso","rt_switch","proteinase","digest_tissue","recover_beads_selection","recover_beads","pre_forward","pre_reverse","pre_pol","pre_pcr","pre_cleanup_selection","pre_cleanup","transposome","tn5","tagment","kit_conditions","native_handoff","library_forward","library_reverse","library_pol","initial_extension","library_pcr","final_cleanup_selection","final_cleanup","final_cleanup_out"),SOURCE)
 graph.annotate(SourceDiscrepancy("tso", "Template-switch oligo versus SMART PCR primer",
     (SourceClaim("Published supplement Table S1 (PMC NIHMS1036108 supplement-1.docx, text line 316)",
                  "TSO AAGCTGGTATCAACGCAGAGTGAATrG+GrG"),
      SourceClaim("same published Table S1 (text line 314)",
                  "SMART PCR primer AAGCAGTGGTATCAACGCAGAGT; the printed TSO lacks the primer's "
                  "AG at positions 5–6 and contains only its 3′ 17 nt"),
      SourceClaim("Stickels et al. 2021 (Slide-seqV2) Supplementary Table 6 and Online Methods "
                  "(same Qiagen order #339414YCO0076714)",
                  "TSO AAGCAGTGGTATCAACGCAGAGTGAATrG+GrG")),
     "Display the TSO exactly as printed; the cDNA PCR is modelled at role level and no SMART-end "
     "bases are inferred. The SMART end is lost from the final library by tagmentation.",
     "The published author-manuscript supplement is a Word document whose text layer gives the value directly.",
     "Whether the TSO omission is a transcription error is unresolved; the Science-hosted "
     "supplement PDF was not obtained."))
 for name, printed, expected in DECODING_MISMATCHES:
     graph.annotate(SourceDiscrepancy("decoder", f"Decoding primer {name} versus its recessed series",
         (SourceClaim("Published supplement Table S1 (supplement-1.docx text lines 317-321)",
                      f"{name} /5Phos/{printed}"),
          SourceClaim("Computed from the same table: the primer series is the reverse complement "
                      "of the bead TruSeq handle CTACACGACGCTCTTCCGATCT, recessed by n 5′ bases "
                      "(Truseq, -2, -3, -4 all follow it)",
                      f"expected /5Phos/{expected}")),
         "Show the primer exactly as printed; the decoding step is modelled at role level and does "
         "not use the printed bases, so no corrected sequence enters the graph.",
         "The other nine primers match their computed recessed forms exactly.",
         "Whether the order used the printed or the series sequence is not stated; the mismatch "
         "is at the fourth base from the primer 3′ end (C printed, G expected)."))
 return g.require_publishable()
