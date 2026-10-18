"""Slide-tags molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
from chemdraw import Construct, Scene, Segment, complement_segments, feature,revcomp
seg=Segment
import illumina as il
import nextera as nx
import seqprimers as sp
from protocol_graph import Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from pcr import PCRGeometry
from protocol_model import role
from adapter_ligation import AdapterLigationGeometry
from endprep import EndPreparedDuplex, dA_tailed_scene

TITLE="Slide-tags — spatial barcodes transferred to intact nuclei"
NOTES="01_slide-tags.html"
SOURCE='Defining source: <a href="https://doi.org/10.1038/s41586-023-06837-4">Russell et al., <i>Nature</i> (2024)</a>.'
SUMMARY="Photocleaved oligos from a sequenced bead array diffuse into overlying tissue nuclei; 10x single-nucleus RNA or Multiome processing then links spatial tags to each nucleus barcode."
CAVEAT="The published TAGS/SLAC transfer logic and separate 3′ v3.1 and Multiome branches are shown. Tn5 tagmentation belongs only to Multiome. The cited 10x CG000317 Rev C and CG000338 Rev F guides print fixed gel-bead scaffolds and library structures; variable cell barcodes, UMIs and puck barcodes remain placeholders here. The <a href=\"https://github.com/broadchenf/Slide-tags/blob/main/cell_barcode_matcher.R\">author parser</a> reads a 16-base cell barcode and 12-base feature UMI from its cell-barcode FASTQ input, and reconstructs a 14-base puck barcode from positions 1–8 and 27–32 of its spatial FASTQ input. The script reverses internal R1/R2 variable names relative to file paths; branch-specific run lengths are not asserted here. Derived from the printed TAGS oligo, the published CG000317 3′ feature-library shell and the Supplementary Methods statement that TAGS spatial libraries were prepared as Cell Surface Protein libraries: Read 2 reads the puck J8 block (cycles 1–8), the printed 18-nt linker (9–26) and the J7 block (27–33). The parser uses cycles 1–8 and 27–32, so of the 15 printed <code>J</code> positions it omits the seventh J7 base. The drawn 3′ tag library applies only to TAGS arrays; the SLAC (ligated-puck) mouse-brain branch, amplified with the named STAG_P701_NEX and P5-Truseq Hybrid primers, is represented as a separate native-endpoint boundary. The Multiome spatial-tag library is a sequencing endpoint in the graph with a visible warning; the puck design used with Multiome remains unassigned in this scoped reconstruction, and no strand-level shell is assigned to it. <b>Published discrepancy:</b> the main article Methods (PDF p. 10; PMC10764288 Methods) print <code>TTT_PC_zCCGG...</code> for the SLAC bead oligo, while the published Supplementary Methods section 1.3 print <code>TTT_PC_GCCGG...</code>. The sequence below uses the published supplement's base-only form because <code>z</code> is not a defined nucleotide; its intended meaning in the main article remains uncertain. The <a href=\"https://doi.org/10.1038/s41587-020-0739-1\">cited Slide-seqV2 Methods</a> define <code>J</code> as a split-pool bead barcode base and also print a no-<code>z</code> SLAC precursor."
CELL=feature("slidetags_cell","cell_barcode","whitelist",whitelist="10x Chromium bead barcode")
SPATIAL=feature("slidetags_spatial","spatial_barcode","whitelist",whitelist="sequenced Slide-tags puck")
UMI=feature("slidetags_umi","umi","random")
I5=feature("slidetags_i5","sample_index","unknown")
I7=feature("slidetags_i7","sample_index","unknown")
RNA=Construct([role("nativeCG317GEX end A; strand assembly unassigned"),role("nuclear cDNA fragment"),role("cell identity; native placement unassigned",molecular_feature=CELL),role("RNAUMI; native placement unassigned",molecular_feature=UMI),role("nativeGEX opposite/index region; widths unassigned")],name="Slide-tags3′RNA native endpoint")
SLAC_RNA=Construct(RNA.segments,name="Slide-tagsSLACbrainRNA native endpoint")
MULTI_RNA=Construct([role("nativeCG338GEX end A; strand assembly unassigned"),role("nuclear cDNA fragment"),role("cell identity; native placement unassigned",molecular_feature=CELL),role("RNAUMI; native placement unassigned",molecular_feature=UMI),role("nativeGEX opposite/index region; widths unassigned")],name="Slide-tagsMultiomeRNA native endpoint")
RNA_P=MULTI_RNA_P=()
# CG000338 Rev F (the Multiome revision the Supplementary Methods 1.5 name),
# Appendix p. 81: ATAC gel-bead oligo and final ATAC library. The 16-base GEM
# barcode and 8-base spacer sit between P5 and the Read 1N (S5 + ME) arm and are
# read by the 24-cycle i5 read (Rev F Sequencing); the only sample index is i7 (N8).
MULTIOME_SPACER="CGCGTCTG"
ATAC_I7=feature("slidetags_atac_i7","sample_index","unknown",
                note="CG000338 Rev F Sample Index N (i7:8); the index set is not printed")
def multiome_atac_library():
    lib=Construct([
        seg("P5",il.P5,"p5"),
        seg("10x cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),
        seg("10x spacer",MULTIOME_SPACER),
        seg("Read 1N S5 arm",nx.S5,"s5"),
        seg("left mosaic end",nx.ME,"me"),
        seg("accessible chromatin","X"*40,placeholder=True),
        seg("right mosaic end reverse complement",nx.ME_RC,"me"),
        seg("Read 2N S7 arm reverse complement",nx.S7_RC,"s7"),
        seg("i7 sample index","I"*8,"cbc",placeholder=True,feature=ATAC_I7),
        seg("P7 reverse complement",il.P7_RC,"p7"),
    ],name="Slide-tags Multiome ATAC library")
    primers=tuple(sp.NEXTERA[k] for k in ("R1","I1","I2","R2"))
    if problems:=sp.verify(lib,primers,required_roles=tuple(p.role for p in primers)):
        raise ValueError("Slide-tags Multiome ATAC library: "+"; ".join(problems))
    return lib,primers
ATAC,ATAC_P=multiome_atac_library()
def bead_primer(name,handle,three_prime):
    """10x gel-bead primer printed as <handle>-N16-N12-<3' end> in the cited guides."""
    return Construct([seg("bead handle",handle,"r1"),
                      seg("cell barcode","N"*16,"cbc",placeholder=True,feature=CELL),
                      seg("UMI","N"*12,"umi",placeholder=True,feature=UMI),
                      *three_prime],name=name)
def phosphorothioate_primer(name,bases,ps_bonds):
    """Custom PCR primer; the last ps_bonds linkages are printed phosphorothioates."""
    return Construct([seg(name,bases,note=f"{ps_bonds} 3'-terminal phosphorothioate linkages (printed *)")],name=name)
STAG_IP7_A1="CAAGCAGAAGACGGCATACGAGATATTTACCGCAGTGACTGGAGTTCAGACGTGT"   # MOESM1 line 89
P5_STAG_IP5_A1="AATGATACGGCGACCACCGAGATCTACACGACAATAAAGACACTCTTTCCCTACACGACGCTC"  # MOESM1 line 89
# Published TAGS bead oligo (Supplementary Methods 1.3, MOESM1 line 73):
# 5'-TTT-PC-<TruSeq Read 2 arm>-J8-<linker>-J7-N7-VV-<capture complement>-3'.
TAGS_R2_ARM="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
TAGS_LINKER="TCTTCAGCGTTCCCGAGA"
TAGS_CAPTURE_COMPLEMENT="GCTTTAAGGCCGGTCCTAGCAA"
TAGS_OLIGO=TAGS_R2_ARM+"J"*8+TAGS_LINKER+"J"*7+"N"*7+"VV"+TAGS_CAPTURE_COMPLEMENT
SLAC_SPLINT_PRINTED_3_TO_5="CTGTGTACCCGCCTCGAAACGATTGC"  # MOESM1 line 70
CAPTURE_SEQUENCE_1="TTGCTAGGACCGGCCTTAAAGC"  # CG000317 Rev C Appendix p. 71
# CG000317 Rev C p. 71 and CG000338 Rev F p. 81 print the same RNA scaffold.
GEM_RNA_PRIMER=bead_primer("10x 3′ v3.1 / Multiome RNA gel-bead primer","CTACACGACGCTCTTCCGATCT",
                           (seg("oligo(dT)30","T"*30),seg("anchor","VN",placeholder=True)))
GEM_FEATURE_PRIMER=bead_primer("10x 3′ v3.1 feature-barcode gel-bead primer","GTCAGATGTGTATAAGAGACAG",
                               (seg("Capture Sequence 1",CAPTURE_SEQUENCE_1,"tso"),))
if revcomp(TAGS_CAPTURE_COMPLEMENT)!=CAPTURE_SEQUENCE_1:
    raise ValueError("TAGS 3' end must be complementary to 10x Capture Sequence 1")
if "A"+il.INDEX1_PRIMER!=revcomp(TAGS_R2_ARM):
    raise ValueError("TAGS 5' arm must yield the TruSeq Read 2/Index 1 arm")
def feature_tag_library():
    # CG000317 Rev C Appendix pp. 71-72 establishes the 3' v3.1 feature-library
    # shell; Supplementary Methods 1.5 says TAGS-array spatial libraries were
    # prepared as Cell Surface Protein libraries. The payload between Capture
    # Sequence 1 and the Read 2 arm is the reverse complement of the published
    # TAGS oligo (derived), so Read 2 reads J8, the 18-nt linker, then J7.
    lib=Construct([
        seg("P5",il.P5,"p5"),
        seg("sample index i5","J"*10,"cbc",placeholder=True,feature=I5),
        seg("Nextera Read 1N",nx.READ1_PRIMER,"r1"),
        seg("cell barcode","C"*16,"cbc",placeholder=True,feature=CELL),
        seg("feature UMI","U"*12,"umi",placeholder=True,feature=UMI),
        seg("Capture Sequence 1",CAPTURE_SEQUENCE_1,"tso"),
        seg("puck VV complement","BB",placeholder=True),
        seg("puck mixed-base N7 complement","N"*7,placeholder=True),
        seg("puck J7 base 7 (Read 2 cycle 33; unused by author parser)","S",
            "cbc",placeholder=True,feature=SPATIAL),
        seg("puck J7 bases 1-6 (Read 2 cycles 27-32)","S"*6,"cbc",placeholder=True,feature=SPATIAL),
        seg("puck linker complement",revcomp(TAGS_LINKER)),
        seg("puck J8 (Read 2 cycles 1-8)","S"*8,"cbc",placeholder=True,feature=SPATIAL),
        seg("Read 2 leading A","A","r2"),
        seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),
        seg("sample index i7 reverse complement","I"*10,"cbc",placeholder=True,feature=I7),
        seg("P7 reverse complement",il.P7_RC,"p7"),
    ],name="Slide-tags 3′ TAGS spatial-tag library")
    primers=(sp.NEXTERA["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
    return lib,primers
TAG,TAG_P=feature_tag_library()
FINAL_LIBRARIES=(
 ("RNA library (3′ branch)",RNA,RNA_P,"10x 3′ v3.1 gene-expression library from tagged nuclei.","The droplet cell barcode and UMI identify each transcript."),
 ("RNA library (Multiome branch)",MULTI_RNA,MULTI_RNA_P,"10x Multiome gene-expression library from tagged nuclei.","The droplet cell barcode and UMI identify each transcript."),
 ("ATAC library (Multiome branch)",ATAC,ATAC_P,"CG000338 Rev F Appendix p. 81 ATAC library from the same spatially tagged nuclei: P5, 16-base GEM barcode and 8-base spacer precede the Read 1N arm; the guide reference identifies a physical8-base sample index, without assigning study index acquisition.","Read 1N and Read 2N read the accessible fragment; the guide reference assigns24cycles to spacer+GEMbarcode, but this does not establish study index acquisition."),
 ("Spatial-tag library (3′ branch, TAGS arrays)",TAG,TAG_P,"CG000317 Cell Surface Protein library shell with the reverse complement of the published TAGS puck oligo between Capture Sequence 1 and the Read 2 arm (derived).","The 10x guide specifies a 12-base feature UMI. Author code reconstructs a 14-base puck barcode from Read 2 cycles 1–8 (J8) and 27–32 (first six J7 bases)."),
)
READ_LENGTHS={"RNA library (3′ branch)":{},"RNA library (Multiome branch)":{},"ATAC library (Multiome branch)":{},"Spatial-tag library":{}}

def oligos():
 return [
  "<b>SLAC ligated-design bead oligo (Supplementary Methods 1.3)</b> — 5′-TTT_PC_GCCGGTAATACGACTCACTATAGGGCTACACGACGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJNNNNNNNVVGCTCGGACACATGGGCG-3′ (main article instead prints <code>TTT_PC_zCCGG...</code>; see discrepancy above)",
  "<b>Direct-synthesis TAGS puck oligo</b> — 5′-TTT-PC-"+TAGS_OLIGO+"-3′ (PC: photocleavable linker; J: split-pool bead barcode base, defined in cited Slide-seqV2 Methods)",
  "<b>SLAC 10x FB1 extension</b> — 5′-GAGCTTTGCTAACGGTCGAGGCTTTAAGGCCGGTCCTAGCAA-3′",
  "<b>SLAC splint</b> — 5′-"+SLAC_SPLINT_PRINTED_3_TO_5[::-1]+"-3′ (source prints it 3′→5′ as "+SLAC_SPLINT_PRINTED_3_TO_5+")",
  "<b>Poly(A) puck oligo</b> — 5′-TTT-PC-GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTJJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJNNNNNNNVVA30-3′",
  "<b>10x 3′ v3.1 RNA gel-bead primer (CG000317 Rev C, p. 71)</b> — 5′-CTACACGACGCTCTTCCGATCT-N16-N12-T30VN-3′",
  "<b>10x 3′ v3.1 feature-barcode gel-bead primer (CG000317 Rev C, p. 71)</b> — 5′-GTCAGATGTGTATAAGAGACAG-N16-N12-TTGCTAGGACCGGCCTTAAAGC-3′",
  "<b>10x Multiome ATAC gel-bead primer (CG000338 Rev F, p. 81)</b> — 5′-AATGATACGGCGACCACCGAGATCTACAC-N16-CGCGTCTG-3′",
  "<b>Multiome spatial-tag spike-in</b> — 5′-GTGACTGGAGTTCAGACGT-3′",
  "<b>Multiome spatial-tag P7 primer</b> — 5′-CAAGCAGAAGACGGCATACGAGATATTTACCGCAGTGACTGGAGTTCAGACGT*G*T-3′ (*: phosphorothioate linkage)",
  "<b>Multiome spatial-tag P5 primer</b> — 5′-AATGATACGGCGACCACCGAGATCTACACGACAATAAAGACACTCTTTCCCTACACGACGC*T*C-3′ (*: phosphorothioate linkage)",
 ]

MULTI_TAG=Construct([role("nativecustomSTAG indexing end A; placement unassigned",molecular_feature=I5),role("cell identity; native width/placement unassigned",molecular_feature=CELL),role("spatialbarcode; Multiomepuckdesign unassigned",molecular_feature=SPATIAL),role("nativecustomSTAG indexing end B; junction unassigned",molecular_feature=I7)],name="Slide-tagsMultiome spatial native boundary")
SLAC_TAG=Construct([role("nativeSTAG_P701_NEX/P5TruseqHybrid end A"),role("cell identity; native assembly unassigned",molecular_feature=CELL),role("SLACspatialbarcode; native capture/junction unassigned",molecular_feature=SPATIAL),role("nativeSLAC end B")],name="Slide-tagsSLACbrain spatial native boundary")
FINAL_LIBRARIES=(*FINAL_LIBRARIES,("RNA library (SLAC brain)",SLAC_RNA,(),"Native3′GEX endpoint for separateSLACbrain experiment; internalCG317orderunassigned.","Study runallocationunassigned."),("Spatial-tag library (SLAC brain)",SLAC_TAG,(),"Named10cyclecustomPCR from10µlcleaned3′supernatant; nativeprimerstrands unassigned.","NoTAGS orMultiome shell substituted."),("Spatial-tag library (Multiome branch)",MULTI_TAG,(),"Source15cyclemodifiedcustomPCR from5µlcleaned4.3supernatant.","Puckdesign/capturejunction unassigned."))
SHOW_READ_CONTENT=False
READ_SPECS={name:tuple(sp.read_spec(r,"unassigned","NextSeq1000P2 100cyclekit or someNovaSeqSP reported; branch-specific acquired cycle allocation unassigned","undetermined") for r in ("Read 1","Read 2","Index 1 (i7)","Index 2 (i5)")) for name,*_ in FINAL_LIBRARIES}
READ_SPECS_NOTE="HistoricalCG317TAGS/CSP andCG338ATAC primer geometry is a source-qualified reference, not evidence of study index/read acquisition. Published article names instruments/kit but no branch run mask is assigned here."


def _tag_parts(design):
 if design=="tags":return [Segment("TAGS5′arm",TAGS_R2_ARM),Segment("puckJ8","N"*8,"cbc",placeholder=True,feature=SPATIAL),Segment("pucklinker",TAGS_LINKER),Segment("puckJ7","N"*7,"cbc",placeholder=True,feature=SPATIAL),Segment("mixedN7VV","N"*7+"VV",placeholder=True),Segment("capturecomplement",TAGS_CAPTURE_COMPLEMENT)]
 if design=="slac":return [role("publishedSLAC fixedarm/capture assembly; native ligated extension unassigned"),Segment("SLACJ8","N"*8,"cbc",placeholder=True,feature=SPATIAL),Segment("SLACpucklinker",TAGS_LINKER),Segment("SLACJ7","N"*7,"cbc",placeholder=True,feature=SPATIAL),role("SLACFB1 capture-extension junction unassigned")]
 return [role("Multiome spatialpuck identity; design/width unassigned",molecular_feature=SPATIAL),role("native capture context unassigned")]


def mixed(design,stage):
 sc=Scene()
 sc.strand("RNA",[role("native nuclear RNA; fate through recovery unassigned")],label="independent RNA population")
 dna="transposed chromatin; donor/gap geometry unassigned" if design=="multi" else "genomic DNA; not a separateTapestri assay"
 sc.strand("DNA",[role(dna)],label="independent DNA population")
 if stage!="tissue":sc.strand("tag",_tag_parts(design),label="independent spatialtag molecules; ambient/residual fraction unassigned")
 if stage in ("rt","preamp","cdna_amp","native"):
  sc.strand("cDNA",[*GEM_RNA_PRIMER.segments,role("first-strand nuclear cDNA; native extent unassigned")],label="ssRT cDNA" if stage=="rt" else "native amplified cDNA population; full ends unassigned")
  if stage!="rt":sc.strand("cDNAopposite",[role("native complementary cDNA copies; terminal extent unassigned")],label="separate PCR-copy population")
  if design=="tags":
   upstream=_tag_parts(design)[:-1]
   sc.strand("featurecopy",[*GEM_FEATURE_PRIMER.segments,*complement_segments(upstream)],label="ssTAGS copy; nativefirst-cycle extent unassigned")
  else:sc.strand("tagcopy",[role("nucleus identity; native tagcopy arrangement unassigned",molecular_feature=CELL),role("spatialidentity; native product/junction unassigned",molecular_feature=SPATIAL)],label="native spatialcopy class; not a completed endpoint")
 return sc


def _clean(g,key,template,label):
 reagent=g.input(key+"_reagent",label,MaterialKind.SELECTION)
 keep=g.product(key+"_retained","source-recovered "+template.label,MaterialKind.MOLECULAR,state=template.state)
 waste=g.product(key+"_waste","physical non-retained fraction; composition unassigned",MaterialKind.WASTE)
 g.react(key,label,StepKind.CLEANUP,((template,InputRole.TEMPLATE),(reagent,InputRole.SELECTION)),(keep,waste),changes_molecule=False)
 g.discard(waste,"Physical supernatant/washes or removedfraction; no species-specific complete depletion claimed.")
 return keep


def _pcr(g,key,template,state,label,forward=None,reverse=None,enzyme_label="native kitPCR polymerase; identity unassigned"):
 f=forward or g.input(key+"_forward","nativePCR primerA; selectedstrand unassigned",MaterialKind.OLIGO)
 r=reverse or g.input(key+"_reverse","nativePCR primerB; relationship unassigned",MaterialKind.OLIGO)
 enzyme=g.input(key+"_polymerase",enzyme_label,MaterialKind.CATALYST)
 product=g.product(key+"_product",label+" products",MaterialKind.MOLECULAR,state=state)
 other=g.product(key+"_residual","co-carried RNA/other polymers; class-specific recovery unassigned",MaterialKind.MOLECULAR,state=template.state)
 g.react(key,label,StepKind.PCR,((template,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(enzyme,InputRole.CATALYST)),(product,other),geometry=PCRGeometry.unavailable(template.state,reason="Actual preceding mixed/native state is the PCRtemplate. Complete first-cycle/nativeend and primerlanding geometry unassigned; endpoint is separate."))
 g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,"PCR eligibility is not physical modality purification; co-carriedclass fate remains unassigned."))
 return product


def _library(g,key,template,endpoint,primers=()):
 lib=g.product(key,"source-qualified sequencing endpoint",MaterialKind.LIBRARY,state=endpoint)
 g.react(key+"_retain","Retain source library for sequencing",StepKind.GENERIC,((template,InputRole.TEMPLATE),),(lib,),changes_molecule=False)
 g.sequence(lib,SequencingPlan(library=endpoint,primers=tuple(primers),acquisition="unspecified",limitation="Study branch-specific run allocation unassigned. Historicalknownprimergeometry is reference scope only."))
 return lib


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 routes={}
 sections=g.input("experiment_sections","collection of separate tissue sections allocated toTAGS3′,SLACbrain3′ andMultiomeexperiments; notsame nuclei processedthreeways",MaterialKind.BIOLOGICAL,primary=True)
 tissues={key:g.product(key+"_tissue","separate tissue section for "+key+" experiment",MaterialKind.BIOLOGICAL,state=mixed(key,"tissue")) for key in ("tags","slac","multi")}
 g.react("allocate_experiments","Represent separately allocated tissue sections forsource experimentvariants",StepKind.SPLIT,((sections,InputRole.TEMPLATE),),tuple(tissues.values()),changes_molecule=False)
 for design,exposure,diffusion in (("tags","365nm30s","7.5min"),("slac","365nm3min","5min"),("multi","design-specific365nmexposure unassigned","design-specificdiffusion unassigned")):
  array_sc=Scene();array_sc.strand("array",[Segment("source5′TTT","TTT"),role("photocleavable bead tether; nativeend chemistry unassigned"),*_tag_parts(design)],label=design+" array; Multiomepuckdesign unassigned" if design=="multi" else design+" published design")
  array=g.input(design+"_array","separate experiment: mapped "+design+" spatial array",MaterialKind.MOLECULAR,state=array_sc,primary=True)
  uv=g.input(design+"_uv",exposure+"; spatialdiffusion "+diffusion,MaterialKind.CONDITION)
  freed_sc=Scene();freed_sc.strand("tag",_tag_parts(design),label="released tags; cleavage-terminalchemistry unassigned")
  freed=g.product(design+"_freed","photoreleased spatialoligos",MaterialKind.MOLECULAR,state=freed_sc)
  remnant=g.product(design+"_array_remnant","post-cleavage bead/tether remnant",MaterialKind.MOLECULAR)
  g.react(design+"_photocleave","Photorelease "+design+" tags and diffuse into tissue",StepKind.CHEMICAL_CONVERSION,((array,InputRole.TEMPLATE),(uv,InputRole.CONDITION)),(freed,remnant))
  g.mark_terminal(Terminal(remnant.id,Disposition.UNRESOLVED,"Arrayremnant not recovered as nuclear library substrate; disposal notassigned."))
  tissue=tissues[design]
  tagged=g.product(design+"_tagged","tagged tissue/nuclear mixture including possible free spatialtags",MaterialKind.MOLECULAR,state=mixed(design,"initial"))
  g.react(design+"_transfer","Transfer local spatialtags into tissue nuclei",StepKind.BINDING,((tissue,InputRole.TEMPLATE),(freed,InputRole.ANALYTE)),(tagged,))
  recovered=_clean(g,design+"_pellet600",tagged,"Nucleusisolation/600gpellet; physicalsupernatant removed")
  recovered=_clean(g,design+"_filter40",recovered,"40µmfiltration; retainfiltrate andremove filterretentate")
  dapi=g.input(design+"_dapi","DAPI staining",MaterialKind.CONDITION)
  stained=g.product(design+"_stained","DAPI-stained nuclear mixture",MaterialKind.MOLECULAR,state=recovered.state)
  g.react(design+"_stain","Stain isolated nuclei withDAPI",StepKind.GENERIC,((recovered,InputRole.TEMPLATE),(dapi,InputRole.CONDITION)),(stained,),changes_molecule=False)
  routes[design]=_clean(g,design+"_pellet200",stained,"200gpellet; physicalsupernatant removed")
 # Separate3′TAGS andSLAC experiments. NoRNA/tag separation precedesGEM/RT orcDNA amplification.
 for design in ("tags","slac"):
  partition=g.input(design+"_gem_condition","PN1000268 / CG000317RevC3′v3.1 GEMpartition",MaterialKind.CONDITION)
  gem=g.product(design+"_gem","mixed3′GEM contents",MaterialKind.MOLECULAR,state=mixed(design,"initial"))
  others=g.product(design+"_other_gems","other partition contents/recovery unassigned",MaterialKind.MOLECULAR)
  g.react(design+"_partition","Partition nuclei with3′gelbeads",StepKind.PARTITION,((routes[design],InputRole.TEMPLATE),(partition,InputRole.CONDITION)),(gem,others),changes_molecule=False)
  g.mark_terminal(Terminal(others.id,Disposition.UNRESOLVED,"Empty/noncalledGEMs are not a physically discarded molecularclass."))
  bead=g.input(design+"_rna_bead","historicalCG317RevC RNAbead scaffold reference",MaterialKind.OLIGO,state=GEM_RNA_PRIMER)
  feature=g.input(design+"_feature_bead","historicalCG317RevC featurebead scaffold reference",MaterialKind.OLIGO,state=GEM_FEATURE_PRIMER)
  rt=g.input(design+"_rt","native3′GEM RT/extension mix; identities unassigned",MaterialKind.CATALYST)
  copied=g.product(design+"_rt_mixture","mixedss nuclearcDNA/spatialcopy reaction",MaterialKind.MOLECULAR,state=mixed(design,"rt"))
  g.react(design+"_rt_capture","Capture RNA andspatialtags in mixed3′RT reaction",StepKind.REVERSE_TRANSCRIPTION,((gem,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(feature,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
  recovery=g.input(design+"_gem_recovery_kit","CG317 nativepostGEM recovery; internalreagents/order unassigned",MaterialKind.CONDITION)
  recovered=g.product(design+"_gem_recovered","nativepostGEM recovered mixedss/copy material",MaterialKind.MOLECULAR,state=copied.state)
  g.react(design+"_gem_recovery","Bounded nativeCG317 postGEM recovery before cDNAamplification",StepKind.MAGIC,((copied,InputRole.TEMPLATE),(recovery,InputRole.CONDITION)),(recovered,),reason="HistoricalCG317 component sourcepass didnotinspect postGEM recovery operations. MixedRT material is carried through without assigning oil/reagent or molecularclass removal from CG338.")
  amplified=_pcr(g,design+"_cdna_pcr",recovered,mixed(design,"preamp"),"Source3′cDNA amplification; nativeGEMrecovery/primerdetails unassigned")
  spri=g.input(design+"_spri","source3′SPRI separation at§2.3; native ratio unassigned",MaterialKind.SELECTION)
  cdna=g.product(design+"_cdna_fraction","retained cDNA-size fraction; class-specific purity unassigned",MaterialKind.MOLECULAR,state=amplified.state)
  short=g.product(design+"_tag_supernatant","retained short spatialcopy-containing SPRIsupernatant§2.3",MaterialKind.MOLECULAR,state=amplified.state)
  g.react(design+"_spri_split","Physically retain cDNA fraction andspatialtag supernatant",StepKind.SELECTION,((amplified,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(cdna,short),changes_molecule=False)
  native_kit=g.input(design+"_gex_kit","CG317 historical3′GEX construction reagents; internalorderunassigned",MaterialKind.CONDITION)
  rna_endpoint=RNA if design=="tags" else SLAC_RNA
  rna=g.product(design+"_rna_native","native3′GEX endpoint",MaterialKind.MOLECULAR,state=rna_endpoint)
  g.react(design+"_gex_construct","Bounded nativeCG3173′GEX construction after cDNA recovery",StepKind.MAGIC,((cdna,InputRole.TEMPLATE),(native_kit,InputRole.CONDITION)),(rna,),reason="PriorCG317source scope inspectedAppendix71–72/scaffolds/CSPendpoint only, not3′GEX internalreactionorder. CG338Multiome order is notprojected intoCG317.")
  _library(g,design+"_rna_library",rna,rna_endpoint)
  short=_clean(g,design+"_supernatant_cleanup",short,"Source cleaned spatialtag supernatant; nativecleanupimplementation unassigned")
  if design=="slac":
   amount=g.input("slac_input_amount","10µl cleaned§2.3supernatant",MaterialKind.CONDITION)
   chosen=g.product("slac_tag_input","10µl mixed spatialtag supernatant",MaterialKind.MOLECULAR,state=short.state)
   rest=g.product("slac_tag_rest","other supernatant outsidePCR aliquot",MaterialKind.MOLECULAR,state=short.state)
   g.react("slac_take","Take10µl cleanedsupernatant",StepKind.SAMPLING,((short,InputRole.TEMPLATE),(amount,InputRole.CONDITION)),(chosen,rest),changes_molecule=False)
   g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,"Unfollowed supernatant; no physicaldiscard claimed."))
   f=g.input("slac_forward","namedP5TruseqHybrid; strand unassigned",MaterialKind.OLIGO)
   r=g.input("slac_reverse","namedSTAG_P701_NEX; strand unassigned",MaterialKind.OLIGO)
   product=_pcr(g,"slac_tag_pcr",chosen,SLAC_TAG,"SLACbrain customspatialPCR10cycles",f,r,"NEBNextHighFidelity")
   _library(g,"slac_tag_library",product,SLAC_TAG)
  else:
   f=g.input("tags_tag_forward","nativeCG317CSPindexPCRprimerA; memberunassigned",MaterialKind.OLIGO)
   r=g.input("tags_tag_reverse","nativeCG317CSPindexPCRprimerB; memberunassigned",MaterialKind.OLIGO)
   product=_pcr(g,"tags_tag_pcr",short,TAG,"TAGS spatialCSPindexPCR:7cycles5.5mmpuck /9cycles3mmpuck",f,r)
   _library(g,"tags_tag_library",product,TAG,TAG_P)
 # Multiome nativeCG338RevF order is positively read, includingmixedpreamp fractionation.
 trans=g.input("multi_transposome","nativeMultiome loaded adapters; donorstrandsunassigned",MaterialKind.ADAPTER)
 tn5=g.input("multi_tn5","Multiome transposition reagents; exactconditionsunassigned",MaterialKind.CATALYST)
 tagged=g.product("multi_transposed","transposed nuclei withRNA andspatialtags retained",MaterialKind.MOLECULAR,state=mixed("multi","initial"))
 off=g.product("multi_other_ends","otherATACend classes",MaterialKind.MOLECULAR)
 g.react("multi_transpose","Tagment accessibleDNA inMultiome nuclei",StepKind.TAGMENTATION,((routes["multi"],InputRole.TEMPLATE),(trans,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,off))
 g.mark_terminal(Terminal(off.id,Disposition.UNRESOLVED,"Endeligibility doesnot establish physicalremoval orzeroamplification."))
 gem_condition=g.input("multi_gem_condition","CG338RevF§2.2 GEMmastermix/transposednuclei partition",MaterialKind.CONDITION)
 gem=g.product("multi_gem","mixedMultiomeGEM contents",MaterialKind.MOLECULAR,state=tagged.state)
 other_gems=g.product("multi_other_gems","otherGEM contents; recoveryunassigned",MaterialKind.MOLECULAR)
 g.react("multi_partition","Partition transposed nuclei inMultiomeGEMs",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(gem_condition,InputRole.CONDITION)),(gem,other_gems),changes_molecule=False)
 g.mark_terminal(Terminal(other_gems.id,Disposition.UNRESOLVED,"Empty/noncalledGEMs are computationalclasses, not asserted physicalwaste."))
 bead=g.input("multi_rna_bead","CG338 RNAbead scaffold historicalcomponentreference",MaterialKind.OLIGO,state=GEM_RNA_PRIMER)
 rt=g.input("multi_rt","nativeMultiomeGEM RT/extension enzymes; componentidentityunassigned",MaterialKind.CATALYST)
 copied=g.product("multi_rt_mixture","mixednativeGEM cDNA/ATAC/spatialcopy material",MaterialKind.MOLECULAR,state=mixed("multi","rt"))
 g.react("multi_rt_capture","NativeGEM RNA/barcode extension withspatialtags present",StepKind.REVERSE_TRANSCRIPTION,((gem,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(copied,))
 recovery=g.input("multi_recovery_agent","CG338§3.1 RecoveryAgent; removeoil/recoveryphase",MaterialKind.SELECTION)
 aqueous=g.product("multi_aqueous","recovered aqueousmixedGEM products",MaterialKind.MOLECULAR,state=copied.state)
 oil=g.product("multi_oil","physicallyremoved oil/RecoveryAgent phase",MaterialKind.WASTE)
 g.react("multi_break","BreakGEMs andretain aqueousphase",StepKind.EXTRACTION,((copied,InputRole.TEMPLATE),(recovery,InputRole.SELECTION)),(aqueous,oil),changes_molecule=False)
 g.discard(oil,"Source physicaloil/RecoveryAgent fraction; no modality-specific moleculardiscard inferred.")
 recovered=_clean(g,"multi_dynabeads",aqueous,"CG338§3.1 Dynabeadsrecovery ofmixedGEM products")
 spike=g.input("multi_spike","spatialtagspikeGTGACTGGAGTTCAGACGT:1µl0.329µM atstep4.1",MaterialKind.OLIGO,state=Construct([Segment("STAGspike","GTGACTGGAGTTCAGACGT")],name="MultiomeSTAGspike"))
 pre=g.input("multi_preamp_primers","CG338PreAmpPrimers plusSTAGspike",MaterialKind.OLIGO)
 premix=g.input("multi_preamp_mix","CG338PreAmpMix;7cycles",MaterialKind.CATALYST)
 amplified=g.product("multi_preamplified","mixedcDNA+ATAC+spatialcopy preampmaterial",MaterialKind.MOLECULAR,state=mixed("multi","preamp"))
 g.react("multi_preamp","Step4.1 mixedpreamplification7cycles",StepKind.PCR,((recovered,InputRole.TEMPLATE),(pre,InputRole.STRUCTURAL),(spike,InputRole.STRUCTURAL),(premix,InputRole.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(recovered.state,reason="Actualmixedss/nativeGEMrecovery material; nativePreAmpPrimer species andfirst-cyclegeometry unassigned."))
 spri=g.input("multi_spri160","CG338§4.3 160µlSPRI; recover160µleluate",MaterialKind.SELECTION)
 eluate=g.product("multi_eluate160","160µlrecovered mixedcDNA/ATAC eluate",MaterialKind.MOLECULAR,state=amplified.state)
 short=g.product("multi_tag_supernatant","retained spatialcopy-containing SPRIsupernatant§4.3",MaterialKind.MOLECULAR,state=amplified.state)
 g.react("multi_preamp_spri","Step4.3 recovereluate andretain spatialtagsupernatant",StepKind.SELECTION,((amplified,InputRole.TEMPLATE),(spri,InputRole.SELECTION)),(eluate,short),changes_molecule=False)
 atac=g.product("multi_atac40","40µl mixedeluate assignedATACindexPCR",MaterialKind.MOLECULAR,state=eluate.state)
 cdna=g.product("multi_cdna35","35µl mixedeluate assignedcDNAPCR",MaterialKind.MOLECULAR,state=eluate.state)
 stored=g.product("multi_stored85","remaining85µl sourceeluate",MaterialKind.MOLECULAR,state=eluate.state)
 g.react("multi_split_eluate","Aliquot160µleluate:40ATAC+35cDNA+85stored",StepKind.SPLIT,((eluate,InputRole.TEMPLATE),),(atac,cdna,stored),changes_molecule=False)
 g.mark_terminal(Terminal(stored.id,Disposition.STORED,"Source remainingeluate stored; mixed aliquots do notpurify modalities."))
 atac_amp=_pcr(g,"multi_atac_index_pcr",atac,ATAC,"CG338§5.1 nativeATACindexPCR")
 atac_clean=_clean(g,"multi_atac_double_spri",atac_amp,"CG338§5.2 doublesidedSPRI; exactretentionratiosunassigned")
 _library(g,"multi_atac_library",atac_clean,ATAC,ATAC_P)
 cdna_amp=_pcr(g,"multi_cdna_pcr",cdna,mixed("multi","cdna_amp"),"CG338§6.1 cDNAamplification")
 cdna_clean=_clean(g,"multi_cdna_spri",cdna_amp,"CG338§6.2 SPRIrecovery ofcDNAamplicons")
 frag_condition=g.input("multi_fragment_reagent","CG338§7.1 fragmentationreagents; identitiesunassigned",MaterialKind.CATALYST)
 fragment_state=Scene();fragment_state.strand("target",[role("cDNAfragment withcellidentity; nativeboundaryunassigned",molecular_feature=CELL),role("UMI-bearing fragmentclass; nativeboundaryunassigned",molecular_feature=UMI),role("RNA-derivedfragmentbody")],label="eligiblecDNAfragmentreference; otherfragmentclasses unassigned")
 fragments=g.product("multi_gex_fragments","nativefragmentedGEX DNA",MaterialKind.MOLECULAR,state=fragment_state)
 g.react("multi_fragment","CG338§7.1 fragmentcDNA",StepKind.FRAGMENTATION,((cdna_clean,InputRole.TEMPLATE),(frag_condition,InputRole.CATALYST)),(fragments,))
 er=g.input("multi_endrepair_reagent","CG338§7.1 endrepairreagents; identitiesunassigned",MaterialKind.CATALYST)
 repaired=g.product("multi_gex_repaired","nativeend-repairedGEX DNA",MaterialKind.MOLECULAR,state=fragment_state)
 g.react("multi_endrepair","CG338§7.1 endrepair",StepKind.END_REPAIR,((fragments,InputRole.TEMPLATE),(er,InputRole.CATALYST)),(repaired,))
 at=g.input("multi_at_reagent","CG338§7.1 Ataillingreagents; identitiesunassigned",MaterialKind.CATALYST)
 prepared=EndPreparedDuplex(Construct(fragment_state.strands["target"].segs,name="nativeGEX fragmentreference"),five_prime_phosphate=False,three_prime_overhang="A",overhang_inferred=False)
 tailed=g.product("multi_gex_tailed","3′A-tailedGEX fragments;5′chemistryunassigned",MaterialKind.MOLECULAR,state=dA_tailed_scene(prepared,label="source3′A;5′statusunassigned"))
 g.react("multi_atail","CG338§7.1 Atailling",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(at,InputRole.CATALYST)),(tailed,))
 selected=_clean(g,"multi_gex_double_spri",tailed,"CG338§7.2 doublesidedSPRI; nativecutoffsunassigned")
 adapter=g.input("multi_gex_adapter","CG338§7.3 nativeadapter; strandgeometryunassigned",MaterialKind.ADAPTER)
 ligase=g.input("multi_gex_ligase","CG338§7.3 ligationenzymes; identitiesunassigned",MaterialKind.CATALYST)
 adapted_state=Scene()
 adapted_state.strand("adapted_fragment",[role("native adapter-derived endA; junction unassigned"),*fragment_state.strands["target"].segs,role("native adapter-derived endB; junction unassigned")],label="actual adapted cDNAfragment class; paired extent unassigned")
 adapted_state.strand("opposite",[role("opposite adapted DNAstrand; native geometry unassigned")],label="separate complementary-copy role")
 adapted_state.strand("other",[role("co-carried RNA/gDNA/shorttags and other fragmentclasses; recovery unassigned")],label="other independent molecular populations")
 adapted=g.product("multi_gex_adapted","nativeadapter-ligatedGEX fragments",MaterialKind.MOLECULAR,state=adapted_state)
 g.react("multi_gex_ligate","CG338§7.3 adapterligation",StepKind.LIGATION,((selected,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),geometry=AdapterLigationGeometry.unavailable("Knownnative3′A substrate; sourcecapdidnotreadexactadapterstrands/junctiongeometry."))
 adapted=_clean(g,"multi_gex_ligation_spri",adapted,"CG338§7.4 SPRIafterligation")
 indexed=_pcr(g,"multi_gex_index_pcr",adapted,MULTI_RNA,"CG338§7.5 nativeGEXindexPCR")
 indexed=_clean(g,"multi_gex_final_spri",indexed,"CG338§7.6 doublesidedSPRIafterindexPCR")
 _library(g,"multi_rna_library",indexed,MULTI_RNA)
 short=_clean(g,"multi_tag_supernatant_cleanup",short,"Sourcecleaned4.3spatialtag supernatant; nativecleanupimplementationunassigned")
 amount=g.input("multi_tag_amount","5µl cleaned4.3supernatant",MaterialKind.CONDITION)
 chosen=g.product("multi_tag_input","5µl nativeSPRI spatialtaginput",MaterialKind.MOLECULAR,state=short.state)
 rest=g.product("multi_tag_rest","otherSPRI supernatant",MaterialKind.MOLECULAR,state=short.state)
 g.react("multi_tag_take","Take5µlcleaned4.3supernatant",StepKind.SAMPLING,((short,InputRole.TEMPLATE),(amount,InputRole.CONDITION)),(chosen,rest),changes_molecule=False)
 g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,"SupernatantoutsidecustomPCRaliquot; fateunassigned."))
 mf=g.input("multi_tag_forward","P5-STAG_ip5_a1;2terminalPSbonds",MaterialKind.OLIGO,state=phosphorothioate_primer("P5-STAG_ip5_a1",P5_STAG_IP5_A1,2))
 mr=g.input("multi_tag_reverse","STAG_iP7_a1;2terminalPSbonds",MaterialKind.OLIGO,state=phosphorothioate_primer("STAG_iP7_a1",STAG_IP7_A1,2))
 product=_pcr(g,"multi_tag_pcr",chosen,MULTI_TAG,"MultiomecustomspatialPCR15cycles",mf,mr)
 _library(g,"multi_tag_library",product,MULTI_TAG)
 g.declare_reviewed_scope(REVIEWED,SOURCE)
 return g.require_publishable()

REVIEWED=("tags_gem_recovery_kit","tags_gem_recovery","slac_gem_recovery_kit","slac_gem_recovery","experiment_sections","allocate_experiments","tags_array","tags_uv","tags_photocleave","tags_tissue","tags_transfer","tags_pellet600_reagent","tags_pellet600","tags_filter40_reagent","tags_filter40","tags_dapi","tags_stain","tags_pellet200_reagent","tags_pellet200","slac_array","slac_uv","slac_photocleave","slac_tissue","slac_transfer","slac_pellet600_reagent","slac_pellet600","slac_filter40_reagent","slac_filter40","slac_dapi","slac_stain","slac_pellet200_reagent","slac_pellet200","multi_array","multi_uv","multi_photocleave","multi_tissue","multi_transfer","multi_pellet600_reagent","multi_pellet600","multi_filter40_reagent","multi_filter40","multi_dapi","multi_stain","multi_pellet200_reagent","multi_pellet200","tags_gem_condition","tags_partition","tags_rna_bead","tags_feature_bead","tags_rt","tags_rt_capture","tags_cdna_pcr_forward","tags_cdna_pcr_reverse","tags_cdna_pcr_polymerase","tags_cdna_pcr","tags_spri","tags_spri_split","tags_gex_kit","tags_gex_construct","tags_rna_library","tags_tag_forward","tags_tag_reverse","tags_tag_pcr_polymerase","tags_tag_pcr","tags_tag_library","slac_gem_condition","slac_partition","slac_rna_bead","slac_feature_bead","slac_rt","slac_rt_capture","slac_cdna_pcr_forward","slac_cdna_pcr_reverse","slac_cdna_pcr_polymerase","slac_cdna_pcr","slac_spri","slac_spri_split","slac_gex_kit","slac_gex_construct","slac_rna_library","slac_input_amount","slac_take","slac_forward","slac_reverse","slac_tag_pcr_polymerase","slac_tag_pcr","slac_tag_library","multi_transposome","multi_tn5","multi_transpose","multi_gem_condition","multi_partition","multi_rna_bead","multi_rt","multi_rt_capture","multi_recovery_agent","multi_break","multi_dynabeads_reagent","multi_dynabeads","multi_spike","multi_preamp_primers","multi_preamp_mix","multi_preamp","multi_spri160","multi_preamp_spri","multi_split_eluate","multi_atac_index_pcr_forward","multi_atac_index_pcr_reverse","multi_atac_index_pcr_polymerase","multi_atac_index_pcr","multi_atac_double_spri_reagent","multi_atac_double_spri","multi_atac_library","multi_cdna_pcr_forward","multi_cdna_pcr_reverse","multi_cdna_pcr_polymerase","multi_cdna_pcr","multi_cdna_spri_reagent","multi_cdna_spri","multi_fragment_reagent","multi_fragment","multi_endrepair_reagent","multi_endrepair","multi_at_reagent","multi_atail","multi_gex_double_spri_reagent","multi_gex_double_spri","multi_gex_adapter","multi_gex_ligase","multi_gex_ligate","multi_gex_ligation_spri_reagent","multi_gex_ligation_spri","multi_gex_index_pcr_forward","multi_gex_index_pcr_reverse","multi_gex_index_pcr_polymerase","multi_gex_index_pcr","multi_gex_final_spri_reagent","multi_gex_final_spri","multi_rna_library","multi_tag_amount","multi_tag_take","multi_tag_forward","multi_tag_reverse","multi_tag_pcr_polymerase","multi_tag_pcr","multi_tag_library","multi_atac_library_retain","multi_rna_library_retain","multi_tag_library_retain","multi_tag_supernatant_cleanup","multi_tag_supernatant_cleanup_reagent","slac_rna_library_retain","slac_supernatant_cleanup","slac_supernatant_cleanup_reagent","slac_tag_library_retain","tags_rna_library_retain","tags_supernatant_cleanup","tags_supernatant_cleanup_reagent","tags_tag_library_retain")
