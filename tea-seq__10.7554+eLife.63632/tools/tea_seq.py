"""TEA-seq molecular model."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"lib"))
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Row, Segment, feature, oligo, revcomp
from pcr import PCRGeometry, PCRPrimer
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role
CELL=feature("tea_cell","cell_barcode","whitelist",whitelist="10x Multiome bead barcode")
RNA_I5=feature("tea_rna_i5","sample_index","unknown")
RNA_I7=feature("tea_rna_i7","sample_index","unknown")
ATAC_I7=feature("tea_atac_i7","sample_index","unknown")
ADT_I7=feature("tea_adt_i7","sample_index","unknown")
RNA_UMI=feature("tea_rna_umi","umi","random")
ADT_UMI=feature("tea_adt_umi","umi","random")
ANTIBODY=feature("tea_antibody","feature_barcode","whitelist",whitelist="TotalSeq-A antibody-tag panel")
PARTIAL_R1=il.TRUSEQ_READ1[11:]
ATAC_SPACER="CGCGTCTG"
# Table 5 supplies this complete 3-prime arm, including terminal phosphorothioate A.
ADT_I7_ARM="GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
ADT_CAPTURE=ADT_I7_ARM[12:]
ADT_P5_PRIMER=il.NEBNEXT_UNIVERSAL_PRIMER[:49]


def rna_bead_primer():
    return Construct([Segment("partial TruSeq Read 1",PARTIAL_R1,"r1"),
        Segment("cell barcode","N"*16,"cbc",placeholder=True,feature=CELL),
        Segment("RNA UMI","N"*12,"umi",placeholder=True,feature=RNA_UMI),
        Segment("poly(dT)","T"*30),Segment("VN anchor","VN",placeholder=True,length_bp=2)])


def rna_library():
    # CG000338 Rev A p80 prints the complete dual-index GEX product.
    return Construct([Segment("P5 flank",il.P5[:-4],"p5"),
        Segment("i5 sample index","N"*10,"cbc",placeholder=True,feature=RNA_I5),
        Segment("TruSeq Read 1",il.TRUSEQ_READ1,"r1"),
        *rna_bead_primer().segments[1:],
        Segment("cDNA insert","X"*40,placeholder=True),
        Segment("dA junction","A"),Segment("TruSeq Read 2 arm",il.INDEX1_PRIMER,"r2"),
        Segment("i7 reverse complement","N"*10,"cbc",placeholder=True,feature=RNA_I7),
        Segment("P7 reverse complement",il.P7_RC,"p7")],name="TEA-seq dual-index gene-expression library")


def atac_library():
    # The cell barcode is at the P5 end and is acquired in forward-workflow I2.
    return Construct([Segment("P5",il.P5,"p5"),
        Segment("shared cell barcode","N"*16,"cbc",placeholder=True,feature=CELL),
        Segment("ATAC spacer",ATAC_SPACER),Segment("Nextera Read 1",nx.READ1_PRIMER,"r1"),
        Segment("accessible genomic DNA","X"*42,placeholder=True),
        Segment("Nextera Read 2 arm",nx.INDEX1_PRIMER,"r2"),
        Segment("i7 reverse complement","N"*8,"cbc",placeholder=True,feature=ATAC_I7),
        Segment("P7 reverse complement",il.P7_RC,"p7")],name="TEA-seq Multiome ATAC library")


def adt_library():
    # Derived from the Multiome bead and Table 5 SI-PCR-Oligo/ADT-i7 pair.
    return Construct([Segment("P5 flank",il.P5[:-4],"p5"),
        Segment("TruSeq Read 1",il.TRUSEQ_READ1,"r1"),
        Segment("shared cell barcode","N"*16,"cbc",placeholder=True,feature=CELL),
        Segment("ADT UMI","N"*12,"umi",placeholder=True,feature=ADT_UMI),
        Segment("poly(dT)","T"*30),Segment("VN anchor","VN",placeholder=True,length_bp=2),
        Segment("antibody barcode reverse complement","N"*15,"cbc",placeholder=True,feature=ANTIBODY),
        Segment("small-RNA Read 2 arm",revcomp(ADT_I7_ARM),"r2"),
        Segment("i7 reverse complement","N"*8,"cbc",placeholder=True,feature=ADT_I7),
        Segment("P7 reverse complement",il.P7_RC,"p7")],name="TEA-seq antibody-tag PCR endpoint")


RNA=rna_library()
ATAC=atac_library()
ADT=adt_library()
FORWARD_INDEX_SOURCE="10x CG000373 Rev A pp2–3; author TEA bench protocol step80 (forward indexing)"
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],
    sp.custom("Index 2 (i5)","P5-facing GEX forward-index priming span",il.P5[:-4],FORWARD_INDEX_SOURCE,
        "The printed P5-side binding span locates the forward-workflow index entry; no added custom sequencing oligo is asserted."),
    sp.TRUSEQ["R2"])
ATAC_P=(sp.NEXTERA["R1"],sp.NEXTERA["I1"],
    sp.custom("Index 2 (i5)","P5-facing ATAC forward-index priming span",il.P5,FORWARD_INDEX_SOURCE,
        "Forward I2 enters the cell barcode from P5; it does not traverse the eight-base spacer."),
    sp.NEXTERA["R2"])
ADT_P=()
for _library,_primers in ((RNA,SEQ_PRIMERS),(ATAC,ATAC_P)):
    _problems=sp.verify(_library,_primers,required_roles=tuple(p.role for p in _primers))
    if _problems:raise ValueError("; ".join(_problems))
if sp.verify(ADT,(sp.TRUSEQ["R1"],),required_roles=("Read 1",)):
    raise ValueError("TEA-seq ADT must retain the Table 5 TruSeq Read 1 site")
TITLE="TEA-seq — transcripts, epitopes and accessibility"
NOTES="01_tea-seq.html"
SOURCE='<a href="https://doi.org/10.7554/eLife.63632">Swanson et al., <i>eLife</i> (2021)</a>.'
MULTIOME_GUIDE='<a href="https://assets.ctfassets.net/an68im79xiti/5tnYerYC4FdDwYK0hlcfWj/cc46120ef02327249b5c42f24f44d699/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevA.pdf">10x CG000338 Rev A, Sequences appendix, p. 80</a>'
SUMMARY="Oligo-tagged antibodies are combined with a 10x Single Cell Multiome workflow so RNA, surface epitopes and accessible chromatin share one cell identity but become three sequencing libraries."
CAVEAT="The cited "+MULTIOME_GUIDE+" prints both gel-bead scaffolds and complete GEX/ATAC library layouts. The authors' forward-indexing joint run acquires 50/10/16/90 cycles; the per-library tables below use their demultiplexing masks. This workflow is confirmed by their bench protocol step 80 and 10x CG000373 Rev A. The ADT endpoint is derived from the original Table 5 PCR primers; exact loaded ADT Read 2 and Index 1 oligos are not established here. Individual barcode members and selected TotalSeq-A tag sequences remain unresolved."
CAVEAT += (" <strong>Apparent vendor-guide inconsistency:</strong> CG000338 Rev A p.80 prints "
           "the GEX P5-side prefix AATGATACGGCGACCACCGAGATCT (25 bases) before i5. "
           "<a href='https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000317_ChromiumNextGEMSingleCell3_v3.1_CellSurfaceProtein_RevE.pdf'>CG000317 Rev E p.72</a> "
           "prints AATGATACGGCGACCACCGAGATCTACAC (29 bases) in both the TT forward primer and GEX product. "
           "Both guides name TT kit PN-1000215 and plate PN-3000431. This model follows the cited Multiome Rev A layout; "
           "the actual historical TEA-seq prefix remains unverified. The difference does not establish a reagent variant.")
CAVEAT += (" <strong>Article vs bench protocol:</strong> the article Methods (TEA-seq library preparation) use "
           "5 µL Pre-Amp Primers per 100 µL preamplification; the authors' protocols.io TEA-seq V.2 step 46 lists 4 µL. "
           "The model does not depend on the volume and does not choose between them. The article calls the final ADT cleanup "
           "a \"dual-sided 1.6×\" SPRIselect cleanup but gives one ratio; protocols.io step 64 gives a single 1.6× cleanup, "
           "which the model labels.")

FINAL_LIBRARIES=(("RNA library",RNA,SEQ_PRIMERS,"Poly(A) RNA captured by the Multiome bead.","Read 1 carries cell barcode and UMI; Read 2 reports cDNA."),("ATAC library",ATAC,ATAC_P,"Tagmented accessible chromatin carrying the same cell identity.","Paired genomic reads use Nextera sites."),("Antibody-tag library",ADT,ADT_P,"Polyadenylated antibody tags selectively amplified from the shared product.","Barcode and UMI identify the epitope observation."))
READ_LENGTHS={"RNA library":{"Read 1":28,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":90},"ATAC library":{"Read 1":50,"Index 1 (i7)":8,"Index 2 (i5)":16,"Read 2":50},"Antibody-tag library":{"Read 1":28,"Index 1 (i7)":8,"Read 2":90}}
ADT_LIMITATION="Table 5 establishes the PCR-derived small-RNA read arm, but the exact loaded ADT Read 2/Index 1 sequencing oligos have not been verified. Their binding sites and barcode-cycle coordinates are not asserted."
READ_SPECS={"Antibody-tag library":(
    sp.read_spec("Read 1",28,"16-base cell barcode and 12-base UMI","extract"),
    sp.read_spec("Index 1 (i7)",8,"ADT sample index","demultiplex"),
    sp.read_spec("Index 2 (i5)",None,"not retained by the ADT demultiplexing mask","not_applicable"),
    sp.read_spec("Read 2",90,"antibody-tag read; exact entry site unresolved","extract"))}
SEQUENCING_UNAVAILABLE={"Antibody-tag library":ADT_LIMITATION}
def oligos():
    """Author Table 5 primers and the cited 10x Rev A gel-bead scaffolds."""
    yield oligo("SI-PCR-Oligo for ADT", [Segment("P5 ADT primer", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC")])
    yield oligo("ADT i7 primer", [Segment("P7 and index flank", "CAAGCAGAAGACGGCATACGAGAT"), Segment("8-base i7 index", "NNNNNNNN", placeholder=True), Segment("ADT-binding arm", "GTGACTGGAGTTCCTTGGCACCCGAGAATTCC"), Segment("terminal phosphorothioate A", "*A", placeholder=True)])
    yield oligo("10x Multiome RNA gel-bead primer (CG000338 Rev A)", [Segment("Read 1 arm", "CTACACGACGCTCTTCCGATCT"), Segment("cell barcode", "N"*16, placeholder=True), Segment("UMI", "N"*12, placeholder=True), Segment("oligo(dT) capture", "T"*30+"VN", placeholder=True)])
    yield oligo("10x Multiome ATAC gel-bead primer (same guide)", [Segment("P5 and index flank", "AATGATACGGCGACCACCGAGATCTACAC"), Segment("cell barcode", "N"*16, placeholder=True), Segment("ATAC spacer", "CGCGTCTG")])
def sections():return [("Stain permeabilized cells",[Row(chunks=[("epitope — antibody — feature barcode — poly(A)","cbc",False)])],"TotalSeq-A oligos encode surface epitopes."),("Tagment and partition",[Row(chunks=[("accessible DNA + poly(A) RNA + antibody tags → one barcoded GEM","me",False)])],"The Multiome workflow associates all three analytes with the same cell."),("Split by analyte",[Row(chunks=[("shared product ├─ RNA\n               ├─ ATAC\n               └─ antibody tags",None,False)])],"Modality-specific PCRs create three libraries.")]


def adt_precursor():
    return Construct([Segment("partial Read 1",PARTIAL_R1,"r1"),
        *ADT.segments[2:7],Segment("copied ADT capture handle",revcomp(ADT_CAPTURE),"r2")],
        name="TEA-seq pre-index ADT cDNA")


def _index_library(g,prefix,template,final,primers,lengths,reviewed,*,geometry,
                   forward_state=None,reverse_state=None,limitation="",conditions=""):
    """Declare this protocol's indexing step with explicit endpoint and geometry."""
    forward=g.input(prefix+"_forward",prefix+" indexing forward primer",MaterialKind.OLIGO,
                    state=forward_state,inferred=forward_state is None)
    reverse=g.input(prefix+"_reverse",prefix+" indexing reverse primer",MaterialKind.OLIGO,
                    state=reverse_state,inferred=reverse_state is None)
    polymerase=g.input(prefix+"_polymerase",prefix+" indexing polymerase",MaterialKind.CATALYST)
    final_material=g.product(prefix+"_library",final.name,MaterialKind.LIBRARY,state=final)
    g.react(prefix+"_pcr","Index and amplify the "+prefix+" library"+(" ("+conditions+")" if conditions else ""),StepKind.PCR,
        ((template,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),
         (reverse,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(final_material,),geometry=geometry)
    g.sequence(final_material,SequencingPlan(library=final,primers=tuple(primers),
        acquired_roles=tuple(lengths),read_lengths=lengths,limitation=limitation))
    reviewed.extend((forward.id,reverse.id,polymerase.id,prefix+"_pcr",final_material.id))


def protocol_graph():
 g=ProtocolGraph(TITLE,authority=SOURCE)
 cells=g.input("cells","permeabilized cells containing RNA, chromatin and surface epitopes",MaterialKind.BIOLOGICAL,state=molecule("TEA-seq cells",(role("poly(A) RNA"),role("accessible chromatin"),role("surface epitopes"))),primary=True)
 antibodies=g.input("antibodies","poly(A)-tailed TotalSeq-A antibodies",MaterialKind.SELECTION,inferred=True,source=SOURCE+" (selected antibody-tag member sequences unverified)")
 stained=g.product("stained","cells retaining epitope-specific antibody tags",MaterialKind.MOLECULAR,state=molecule("stained cells",(role("poly(A) RNA"),role("accessible chromatin"),role("antibody tags"))))
 unbound=g.product("unbound","unbound antibody tags",MaterialKind.WASTE)
 g.react("stain","Stain permeabilized cells with oligo-tagged antibodies",StepKind.BINDING,((cells,InputRole.TEMPLATE),(antibodies,InputRole.SELECTION)),(stained,unbound)); g.discard(unbound,"does not report a cell-associated epitope")
 transposome=g.input("transposome","ATAC Enzyme B PN-2000265 loaded transposome",MaterialKind.ADAPTER,inferred=True,source=SOURCE+" (source names the Multiome reagent; its loaded-end molecular state is undeclared in this model)")
 tn5=g.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
 tagged=g.product("tagged","stained cells retaining tagmented chromatin and RNA",MaterialKind.MOLECULAR,state=stained.state)
 twaste=g.product("tn5_waste","non-library tagmentation products",MaterialKind.WASTE)
 g.react("tagment","Tagment accessible chromatin",StepKind.TAGMENTATION,((stained,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,twaste)); g.discard(twaste,"does not retain an amplifiable accessibility fragment")
 droplets=g.input("droplets","10x Multiome GEM partitioning condition",MaterialKind.CONDITION)
 occupied=g.product("occupied","shared GEM contents before RT and barcoding",MaterialKind.MOLECULAR,state=tagged.state)
 failed=g.product("failed","empty or unresolved GEM products",MaterialKind.WASTE)
 g.react("partition","Partition cells with barcoded gel beads",StepKind.PARTITION,((tagged,InputRole.TEMPLATE),(droplets,InputRole.CONDITION)),(occupied,failed),changes_molecule=False); g.discard(failed,"does not retain a resolvable cell barcode")
 bead=g.input("bead_primer","cell-barcode/UMI oligo(dT) bead primer",MaterialKind.OLIGO,state=rna_bead_primer(),source=MULTIOME_GUIDE)
 rt=g.input("rt","reverse transcriptase",MaterialKind.CATALYST)
 tso=g.input("tso","Template Switch Oligo PN-3000228",MaterialKind.OLIGO,inferred=True,
     source=SOURCE+" (named reagent; exact TSO state undeclared in this graph)")
 atac_bead=g.input("atac_bead_primer","Multiome P5/cell-barcode ATAC bead primer",MaterialKind.OLIGO,
     state=Construct([Segment("P5",il.P5),Segment("cell barcode","N"*16,placeholder=True,feature=CELL),Segment("ATAC spacer",ATAC_SPACER)]),source=MULTIOME_GUIDE)
 barcoding=g.input("barcoding_enzyme","Multiome Barcoding Enzyme Mix PN-2000266",MaterialKind.CATALYST)
 pool_state=molecule("barcoded TEA-seq mixture",(role("cell-barcoded transcript cDNA"),role("cell-barcoded ATAC DNA"),role("cell-barcoded ADT cDNA")))
 barcoded=g.product("barcoded_pool","barcoded GEX, ATAC and ADT molecules in one pool",MaterialKind.MOLECULAR,state=pool_state)
 g.react("reverse_transcribe","GEM reverse transcription, template switching and ATAC barcoding (37°C 45 min, 25°C 30 min; then Quenching Agent)",StepKind.REVERSE_TRANSCRIPTION,
     ((occupied,InputRole.TEMPLATE),(bead,InputRole.STRUCTURAL),(atac_bead,InputRole.STRUCTURAL),
      (tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST),(barcoding,InputRole.CATALYST)),(barcoded,))
 cleanup=g.input("gem_cleanup","SILANE and 2.0x SPRI cleanup conditions",MaterialKind.SELECTION)
 recovered=g.product("recovered_pool","recovered barcoded TEA-seq molecules",MaterialKind.MOLECULAR,state=pool_state)
 spent=g.product("gem_cleanup_waste","removed barcoding reagents and unretained molecules",MaterialKind.WASTE)
 g.react("recover","Recover GEM products and remove barcoding reagents",StepKind.CLEANUP,
     ((barcoded,InputRole.TEMPLATE),(cleanup,InputRole.SELECTION)),(recovered,spent),changes_molecule=False)
 g.discard(spent,"does not enter the retained preamplification mixture")
 preamp=g.input("preamp_primers","Pre-Amp Primers PN-2000271",MaterialKind.OLIGO,inferred=True,
     source=SOURCE+" (named reagent; exact primer states undeclared in this model)")
 adt_amp=g.input("adt_rev_amp","ADT-Rev-AMP",MaterialKind.OLIGO,
     state=Construct([Segment("ADT-Rev-AMP",ADT_CAPTURE[:-1])]))
 preamp_pol=g.input("preamp_polymerase","Pre-Amp Mix PN-2000270",MaterialKind.CATALYST)
 amplified=g.product("preamp_products","preamplified ATAC, GEX and ADT mixture before cleanup",MaterialKind.MOLECULAR,state=pool_state)
 g.react("preamp","Preamplify the shared molecules with added ADT-Rev-AMP (7 cycles)",StepKind.PCR,
     ((recovered,InputRole.TEMPLATE),(preamp,InputRole.STRUCTURAL),(adt_amp,InputRole.STRUCTURAL),
      (preamp_pol,InputRole.CATALYST)),(amplified,),
     geometry=PCRGeometry.unavailable(pool_state,reason="The source describes shared ATAC/cDNA preamplification; the commercial primer mixture and multiple precursor duplexes are not assembled in this model."))
 preamp_spri=g.input("preamp_spri","2.0× SPRIselect cleanup",MaterialKind.SELECTION)
 preamplified=g.product("preamplified_pool","purified preamplified ATAC, GEX and ADT mixture",MaterialKind.MOLECULAR,state=pool_state)
 preamp_waste=g.product("preamp_spri_waste","unbound primers and reagents",MaterialKind.WASTE)
 g.react("preamp_cleanup","Purify the preamplified mixture (2.0× SPRIselect)",StepKind.CLEANUP,
     ((amplified,InputRole.TEMPLATE),(preamp_spri,InputRole.SELECTION)),(preamplified,preamp_waste),changes_molecule=False)
 g.discard(preamp_waste,"removed in the 2.0× SPRIselect supernatant")
 atac=g.product("atac_branch","aliquot for ATAC sample-index PCR",MaterialKind.MOLECULAR,state=pool_state)
 cdna_aliquot=g.product("cdna_aliquot","aliquot for cDNA and ADT amplification",MaterialKind.MOLECULAR,state=pool_state)
 g.react("allocate","Allocate preamplified material to ATAC and cDNA workflows",StepKind.SPLIT,
     ((preamplified,InputRole.TEMPLATE),),(atac,cdna_aliquot),changes_molecule=False)
 cdna_primers=g.input("cdna_primers","cDNA Primers PN-2000089",MaterialKind.OLIGO,inferred=True,
     source=SOURCE+" (exact primer pair undeclared in this model)")
 cdna_pol=g.input("cdna_polymerase","cDNA amplification mix",MaterialKind.CATALYST)
 expanded=g.product("expanded_cdna_pool","amplified transcript and ADT cDNA",MaterialKind.MOLECULAR,
     state=molecule("amplified cDNA mixture",(role("transcript cDNA"),role("ADT cDNA"))))
 g.react("amplify_cdna","Amplify cDNA and ADTs with added ADT-Rev-AMP",StepKind.PCR,
     ((cdna_aliquot,InputRole.TEMPLATE),(cdna_primers,InputRole.STRUCTURAL),
      (adt_amp,InputRole.STRUCTURAL),(cdna_pol,InputRole.CATALYST)),(expanded,),
     geometry=PCRGeometry.unavailable(pool_state,reason="The source names the cDNA primer mixture and ADT-Rev-AMP; commercial primer-pair binding geometry is not assembled in this model."))
 selection=g.input("cdna_size_selection","0.6x followed by 2.0x SPRI size separation",MaterialKind.SELECTION)
 cdna=g.product("cdna","large transcript-cDNA fraction",MaterialKind.MOLECULAR,
     state=molecule("GEX cDNA",(role("cell barcode"),role("UMI"),role("transcript cDNA"))))
 copied=g.product("copied_tags","small ADT-cDNA fraction",MaterialKind.MOLECULAR,state=adt_precursor())
 g.react("separate_cdna","Separate large cDNA and small ADT products by size",StepKind.SELECTION,
     ((expanded,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(cdna,copied),changes_molecule=False)
 reviewed=["cells","antibodies","stain","transposome","tn5","tagment","droplets","partition",
     "bead_primer","rt","tso","atac_bead_primer","barcoding_enzyme","reverse_transcribe","gem_cleanup","recover",
     "preamp_primers","adt_rev_amp","preamp_polymerase","preamp","preamp_spri","preamp_cleanup","allocate","cdna_primers",
     "cdna_polymerase","amplify_cdna","cdna_size_selection","separate_cdna"]
 fragmentation=g.input("rna_fragmentation_enzyme","GEX fragmentation enzyme",MaterialKind.CATALYST)
 fragments=g.product("rna_fragments","retained GEX fragments",MaterialKind.MOLECULAR,state=cdna.state)
 g.react("rna_fragment","Fragment amplified GEX cDNA (joint fragmentation/end-repair/A-tailing reaction: 32°C 5 min, 65°C 30 min)",StepKind.FRAGMENTATION,
         ((cdna,InputRole.TEMPLATE),(fragmentation,InputRole.CATALYST)),(fragments,))
 endrepair=g.input("rna_end_repair","GEX end-repair enzyme mixture",MaterialKind.CATALYST)
 repaired=g.product("rna_repaired","end-repaired GEX fragments",MaterialKind.MOLECULAR,state=fragments.state)
 g.react("rna_repair","End-repair GEX fragments",StepKind.END_REPAIR,
         ((fragments,InputRole.TEMPLATE),(endrepair,InputRole.CATALYST)),(repaired,))
 tailase=g.input("rna_a_tailing","GEX A-tailing enzyme mixture",MaterialKind.CATALYST)
 tailed=g.product("rna_tailed","A-tailed GEX fragments",MaterialKind.MOLECULAR,state=repaired.state)
 g.react("rna_tail","A-tail GEX fragments",StepKind.TAILING,
         ((repaired,InputRole.TEMPLATE),(tailase,InputRole.CATALYST)),(tailed,))
 adapter=g.input("rna_adapter","GEX Adapter Oligos PN-2000094",MaterialKind.ADAPTER)
 ligase=g.input("rna_ligase","GEX DNA ligase",MaterialKind.CATALYST)
 adapted=g.product("rna_adapted","adapter-ligated GEX fragments",MaterialKind.MOLECULAR,state=tailed.state)
 g.react("rna_ligate","Ligate the GEX sequencing adapter (after 0.6×/0.8× SPRI; 20°C 15 min)",StepKind.LIGATION,
         ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),
         geometry=AdapterLigationGeometry.unavailable("The article names PN-2000094 and the guide prints the finished library; its ligation-ready end geometry is not assembled in this model."))
 reviewed.extend(("rna_fragmentation_enzyme","rna_fragment","rna_end_repair","rna_repair",
                  "rna_a_tailing","rna_tail","rna_adapter","rna_ligase","rna_ligate"))
 _index_library(g,"rna",adapted,RNA,SEQ_PRIMERS,READ_LENGTHS["RNA library"],reviewed,
     conditions="after 0.8× SPRI; Dual Index TT Set A, 14 cycles; 0.6×/0.8× SPRI",
     geometry=PCRGeometry.unavailable(adapted.state,reason="The completed GEX product is printed in CG000338 Rev A; a selected TT primer pair and its precursor sites have not been declared in this graph."))
 _index_library(g,"atac",atac,ATAC,ATAC_P,READ_LENGTHS["ATAC library"],reviewed,
     conditions="SI-PCR Primer B + Sample Index Plate N Set A, 9 cycles; 0.6×/1.6× SPRI",
     geometry=PCRGeometry.unavailable(atac.state,reason="The completed ATAC product is printed in CG000338 Rev A; the barcoded pre-index duplex and selected sample-index primer pair are not yet assembled in this graph."))
 adt_forward=PCRPrimer("SI-PCR-Oligo",ADT_P5_PRIMER,anneal_nt=13)
 adt_reverse=PCRPrimer("ADT-i7 primer",il.P7+"N"*8+ADT_I7_ARM,anneal_nt=len(ADT_CAPTURE))
 _index_library(g,"adt",copied,ADT,ADT_P,READ_LENGTHS["Antibody-tag library"],reviewed,
     conditions="KAPA HiFi, 15 cycles; 1.6× SPRI",
     forward_state=Construct([Segment(adt_forward.name,adt_forward.sequence)]),
     reverse_state=Construct([Segment("P7",il.P7),Segment("i7 index","N"*8,placeholder=True,feature=ADT_I7),Segment("ADT binding arm",ADT_I7_ARM)]),
     geometry=PCRGeometry.exact(copied.state,adt_forward,adt_reverse,
         statement="TEA-seq Supplementary file 5 prints SI-PCR-Oligo and ADT-i7 primer; the precursor uses the printed Multiome bead scaffold."),
     limitation=ADT_LIMITATION)
 g.declare_reviewed_scope(reviewed,SOURCE)
 return g.require_publishable()
