"""MATQ-seq source-owned mixed priming, tailing, conversion and library graph."""
from chemdraw import Construct,Scene,Segment,oligo,revcomp,feature
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind,Disposition,Terminal,SourceDiscrepancy,SourceClaim
from protocol_model import molecule,role
import seqprimers as sp
from endprep import repair_and_dA_tail,dA_tailed_scene
seg=Segment
GAT27 = "GTGAGTGATGGTTGAGGATGTGTGGAG"
GAT27_DT = GAT27 + "N" * 5 + "T" * 20
GAT27_5N3G = GAT27 + "N" * 5 + "GGG"
GAT27_5N3T = GAT27 + "N" * 5 + "TTT"
GAT21_6N3G = "GATGGTTGAGGATGTGTGGAG" + "N" * 6 + "GGG"
THREE_N_GAT24 = "N" * 3 + "AGTGATGGTTGAGGATGTGTGGAG"

TITLE="MATQ-seq — multiple annealing and dC-tailing single-cell RNA-seq"
NOTES="01_matq-seq.html"
SOURCE='Sheng et al., Nature Methods(2017), <a href="https://doi.org/10.1038/nmeth.4145">doi:10.1038/nmeth.4145</a>; published supplementary protocol.'
SUMMARY="Mixed oligo-dT/random RT primers repeatedly sample total RNA; variable dC tails support second-strand synthesis before GAT27 amplification and3NGAT24 conversion."
CAVEAT="INFERRED library ends — the dC tract length and complete multi-priming product geometries are unassigned, not fixed12nt tails or annealed-primer views used as synthesized products. Native TruSeq adapter/index species remain unassigned. Prior article2×75 versus SupplementaryFig21 85×2 read statements remain distinct below."
FINAL_LIBRARY=Construct([role("TruSeq-family left end; physical species unassigned",inferred=True),role("sheared MATQ amplified cDNA"),role("TruSeq-family opposite end/index architecture unassigned",inferred=True,molecular_feature=feature("matq_sample_index","sample_index","unknown"))],name="MATQ source-limited TruSeq endpoint")
SEQ_PRIMERS=()
FINAL_CAPTION="INFERRED library ends — historical TruSeq endpoint roles; native primer/adapter geometry and index widths remain unassigned."
SEQUENCING_INTRO="Prior defining OnlineMethods review records NextSeq5002×75; published SupplementaryFig21 says85×2. The graph selects2×75 without merging the alternatives."
SEQUENCING_ENDING=SEQUENCING_INTRO
SEQUENCING_UNAVAILABLE="Physical sequencing primers and acquired index reads remain unassigned."
READ_LENGTHS={"Read 1":75,"Read 2":75}
READ_SPECS=(
    sp.read_spec("Read 1",75,"MATQ cDNA; physical primer site and exact entry unassigned","include"),
    sp.read_spec("Read 2",75,"Opposite cDNA end; physical primer site unassigned","include"),
    sp.read_spec("Index 1 (i7)","unassigned","Index acquisition and cycles unassigned","undetermined"),
    sp.read_spec("Index 2 (i5)","unassigned","Index acquisition and cycles unassigned","undetermined"),
)
READ_SPECS_NOTE="Selected prior defining2×75 run; SupplementaryFig21 85×2 remains a separate claim. Unknown index acquisition is not a negative acquisition finding."

def oligos():
    yield oligo("GAT27dT", [seg("GAT27", GAT27, "tso"),
                             seg("N5", "N" * 5, placeholder=True), seg("dT20", "T" * 20)])
    yield oligo("GAT27-5N3G", [seg("GAT27", GAT27, "tso"),
                                seg("N5", "N" * 5, placeholder=True), seg("3G", "GGG")])
    yield oligo("GAT27-5N3T", [seg("GAT27", GAT27, "tso"),
                                seg("N5", "N" * 5, placeholder=True), seg("3T", "TTT")])
    yield oligo("GAT21-6N3G", [seg("GAT21", "GATGGTTGAGGATGTGTGGAG", "tso"),
                                seg("N6", "N" * 6, placeholder=True), seg("3G", "GGG")])
    yield oligo("GAT27 PCR", [seg("GAT27", GAT27, "tso")])
    yield oligo("3NGAT24", [seg("N3", "N" * 3, placeholder=True),
                             seg("GAT24", "AGTGATGGTTGAGGATGTGTGGAG", "tso")])


def rt_primers():
    sc=Scene()
    for key,seq in (("dT",GAT27_DT),("3G",GAT27_5N3G),("3T",GAT27_5N3T)):
        sc.strand(key,[Segment("GAT27 "+key,seq)],label="source RT primer "+key)
    return sc

def target_copies(*,tailed=False):
    sc=Scene()
    for key,primer in rt_primers().strands.items():
        parts=[*primer.segs,role("RNA-derived first-strand cDNA; variable extent")]
        if tailed:parts.append(role("TdT dC tract; length unassigned"))
        sc.strand(key,parts,label="independent "+key+"-primed cDNA target class")
    return sc

def duplex_products(label):
    sc=Scene()
    sc.strand("representative strand",[role(label+" cDNA product; native handle/loop extents unassigned")],label=label+" representative strand")
    sc.strand("opposite strand",[role(label+" complementary product; native geometry unassigned")],label=label+" opposite strand")
    return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    rna=g.input("rna","single-cell total RNA",MaterialKind.BIOLOGICAL,state=molecule("total RNA",(role("poly(A) and other RNA templates"),)),primary=True)
    primers=g.input("mixed_rt_primers","GAT27dT1.5µM/GAT27-5N3G5µM/GAT27-5N3T5µM mix",MaterialKind.OLIGO,state=rt_primers())
    rt=g.input("reverse_transcriptase","SuperScriptIII;10stepped8–50°Ccycles,then50°C15min",MaterialKind.CATALYST)
    first=g.product("first_strand","three priming classes of first-strand cDNA",MaterialKind.MOLECULAR,state=target_copies())
    template=g.product("rt_template_rna","RNA templates remaining after RT",MaterialKind.MOLECULAR,state=rna.state)
    g.react("reverse_transcribe","Repeated mixed-primer RT",StepKind.REVERSE_TRANSCRIPTION,((rna,InputRole.TEMPLATE),(primers,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,template))
    cleanup=g.input("cleanup_enzymes","T4DNApolymerase primer digestion;RNaseH/RNaseI37°C15min,72°C15min",MaterialKind.CATALYST)
    cleaned=g.product("cleaned","cDNA targets after source primer/RNA digestion",MaterialKind.MOLECULAR,state=first.state)
    debris=g.product("digestion_products","RNA/unused-primer hydrolysis products; not physically removed by digestion",MaterialKind.MOLECULAR)
    g.react("cleanup","Digest unused primers and RNA templates in the same tube",StepKind.DIGESTION,((first,InputRole.TEMPLATE),(template,InputRole.TEMPLATE),(cleanup,InputRole.CATALYST)),(cleaned,debris))
    g.mark_terminal(Terminal(debris.id,Disposition.UNRESOLVED,"Source enzyme digestion/heat treatment changes molecules; physical handling of small products through later preparation is not reconstructed."))
    tdt=g.input("tdt","TdT withdCTP,37°C15min,72°C15min",MaterialKind.CATALYST)
    tailed=g.product("tailed","variable dC-tailed cDNA target populations",MaterialKind.MOLECULAR,state=target_copies(tailed=True))
    g.react("tail","Add terminal dC tract; do not assign tail length",StepKind.TAILING,((cleaned,InputRole.TEMPLATE),(tdt,InputRole.CATALYST)),(tailed,))
    secondp=g.input("second_strand_primer","GAT21-6N3G",MaterialKind.OLIGO,state=GAT21_6N3G)
    vent=g.input("deep_vent","DeepVent exo−;10×48°C20s/72°C1min",MaterialKind.CATALYST)
    second=g.product("duplex","second-strand extension/amplification products",MaterialKind.MOLECULAR,state=duplex_products("post-GAT21 synthesis"))
    g.react("second_strand","GAT21-6N3G priming and second-strand synthesis",StepKind.EXTENSION,((tailed,InputRole.TEMPLATE),(secondp,InputRole.STRUCTURAL),(vent,InputRole.CATALYST)),(second,))
    p=g.input("gat27_primer","single source GAT27 PCR primer",MaterialKind.OLIGO,state=GAT27)
    wta=g.product("wta_product","GAT27-amplified cDNA",MaterialKind.MOLECULAR,state=duplex_products("GAT27 amplification"))
    g.react("wta","DeepVent24-cycle GAT27 amplification",StepKind.PCR,((second,InputRole.TEMPLATE),(p,InputRole.STRUCTURAL),(vent,InputRole.CATALYST)),(wta,),geometry=PCRGeometry.unavailable(second.state,symmetric_species=True,reason="Actual post-second-strand product classes enter this single-primer amplification. Native loop/handle first-cycle geometry is unassigned."))
    column=g.input("wta_column","QiagenPCR purification",MaterialKind.SELECTION)
    cleanw=g.product("wta_clean","column-purified WTA product",MaterialKind.MOLECULAR,state=wta.state)
    cw=g.product("wta_column_waste","column flowthrough/washes",MaterialKind.WASTE)
    g.react("wta_cleanup","Purify WTA product",StepKind.CLEANUP,((wta,InputRole.TEMPLATE),(column,InputRole.SELECTION)),(cleanw,cw),changes_molecule=False);g.discard(cw,"source column exclusion")
    conv=g.input("gat24_primer","3NGAT24 conversion primer",MaterialKind.OLIGO,state=THREE_N_GAT24)
    converted=g.product("converted","3NGAT24 converted duplex target population",MaterialKind.MOLECULAR,state=duplex_products("post-3NGAT24 conversion"))
    g.react("conversion","200ng input;source conversion lists two20-cycle blocks with intervening95°C",StepKind.PCR,((cleanw,InputRole.TEMPLATE),(conv,InputRole.STRUCTURAL),(vent,InputRole.CATALYST)),(converted,),geometry=PCRGeometry.unavailable(cleanw.state,symmetric_species=True,reason="The published protocol lists a single3NGAT24 species and two20-cycle blocks; conversion-end geometry remains undeclared rather than borrowing the final sequencing library."))
    cb=g.input("conversion_beads","1.2×AMPureXP",MaterialKind.SELECTION)
    amplified=g.product("amplified","purified converted MATQ cDNA",MaterialKind.MOLECULAR,state=converted.state)
    cwater=g.product("conversion_waste","AMPure supernatant/washes",MaterialKind.WASTE)
    g.react("conversion_cleanup","Purify converted cDNA",StepKind.CLEANUP,((converted,InputRole.TEMPLATE),(cb,InputRole.SELECTION)),(amplified,cwater),changes_molecule=False);g.discard(cwater,"source1.2×bead cleanup")
    shear=g.input("shear_condition","CovarisS220~200bp;prior defining OnlineMethods scope retained",MaterialKind.CONDITION)
    fragments=g.product("fragments","sheared cDNA fragments",MaterialKind.MOLECULAR,state=molecule("fragment target population",(role("sheared MATQ cDNA; inherited handle classes not retained whole"),)))
    g.react("shear","Shear converted amplicons",StepKind.FRAGMENTATION,((amplified,InputRole.TEMPLATE),(shear,InputRole.CONDITION)),(fragments,))
    er=g.input("end_repair_mix","NEBNext EndRepairModule",MaterialKind.CATALYST)
    repaired=g.product("repaired","end-repaired fragments",MaterialKind.MOLECULAR,state=fragments.state)
    g.react("end_repair","Repair fragment ends",StepKind.END_REPAIR,((fragments,InputRole.TEMPLATE),(er,InputRole.CATALYST)),(repaired,))
    k=g.input("klenow_exo","Klenow exo− dA-tailing",MaterialKind.CATALYST)
    tailed2=g.product("da_tailed","dA-tailed fragments",MaterialKind.MOLECULAR,state=dA_tailed_scene(repair_and_dA_tail(Construct([role("sheared MATQ cDNA")]))))
    g.react("da_tail","dA-tail repaired fragments",StepKind.TAILING,((repaired,InputRole.TEMPLATE),(k,InputRole.CATALYST)),(tailed2,))
    adapter=g.input("illumina_adapter","historical TruSeq adapter; selected species unassigned",MaterialKind.ADAPTER)
    ligase=g.input("adapter_ligase","T4DNA QuickLigase",MaterialKind.CATALYST)
    ligated=g.product("ligated","adapter-ligated fragments",MaterialKind.MOLECULAR,state=molecule("actual pre-index fragment",(role("TruSeq-adapter-bearing fragment; native strands/ends unassigned"),)))
    g.react("ligate_adapter","Ligate source-named TruSeq adapters",StepKind.LIGATION,((tailed2,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.unavailable("Historical TruSeq variant/physical adapter ends remain unassigned in prior article source scope; native geometry is not supplied by a canonical shell."))
    f=g.input("library_forward","TruSeqPCR primer mix first binding role; species unassigned",MaterialKind.OLIGO)
    r=g.input("library_reverse","same supplied primer mix opposite binding role; species unassigned",MaterialKind.OLIGO)
    pol=g.input("library_polymerase","KAPAHiFi HotStart;5cycles from prior defining OnlineMethods review",MaterialKind.CATALYST)
    lib=g.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,state=FINAL_LIBRARY)
    g.react("library_pcr","INFERRED ends — amplify sequencing library for5cycles",StepKind.PCR,((ligated,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(lib,),geometry=PCRGeometry.unavailable(ligated.state,reason="Actual adapter-ligated fragment is the input, with historical primer/end geometry unassigned."))
    g.sequence(lib,SequencingPlan(library=FINAL_LIBRARY,acquired_roles=("Read 1","Read 2"),read_lengths=READ_LENGTHS,limitation=SEQUENCING_INTRO+" "+SEQUENCING_UNAVAILABLE))
    g.annotate(SourceDiscrepancy("library","Source run alternatives remain unresolved",(SourceClaim("Prior defining OnlineMethods cost paragraph","NextSeq5002×75 midoutput kit"),SourceClaim("Published SupplementaryFig21","paired85×2 bases")),"Select2×75; keep85×2 as separate source claim.","Existing article/facet source review retained; Supplementary figure source inspected in this round.","Exact relation between those statements remains unassigned."))
    g.declare_reviewed_scope(("rna","mixed_rt_primers","reverse_transcriptase","reverse_transcribe","cleanup_enzymes","cleanup","tdt","tail","second_strand_primer","deep_vent","second_strand","gat27_primer","wta","wta_column","wta_cleanup","gat24_primer","conversion","conversion_beads","conversion_cleanup","shear_condition","shear","end_repair_mix","end_repair","klenow_exo","da_tail","illumina_adapter","adapter_ligase","ligate_adapter","library_forward","library_reverse","library_polymerase","library_pcr","library"),SOURCE)
    return g.require_publishable()
