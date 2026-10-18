"""CG000239 Rev F owns the Visium fresh-frozen molecular route and readout."""
from dataclasses import replace
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp
from adapter_ligation import AdapterLigationGeometry
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import DetailView, Disposition, InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind, Terminal
from protocol_model import molecule, role
TITLE = "10x Visium Spatial Gene Expression — array capture"
NOTES = "01_visium.html"
SOURCE = ('Commercial protocol: 10x Genomics <a href="https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000239_Visium_Spatial_Gene_Expression_User_Guide_Rev_F.pdf">CG000239 Rev F</a>. '
          'The <a href="https://cdn.10xgenomics.com/raw/upload/v1655151898/support/in-line%20documents/Dual_Index_Kit_TT_Set_A.json">Dual Index TT Set A table</a> publishes the 96 well-specific i7/i5 pairs.')

SUMMARY = "Array-bound oligo-dT captures local RNA, installing spatial barcode and UMI before second-strand synthesis, release and dual-index library construction."
CAVEAT = ("INFERRED — KOH-released strand identity: the guide prints both cDNA strands but does not "
          "establish the slide attachment chemistry. The model follows transfer of the complementary "
          "second strand. Printed oligo/library patterns are known; individual spot, UMI and index "
          "members vary. Fragment classes and native PCR binding geometry remain bounded below.")
SLIDE_READ1 = il.TRUSEQ_READ1[-22:]
TSO = "AAGCAGTGGTATCAACGCAGAGTACATrGrGrG"
SECOND_STRAND_PRIMER = "AAGCAGTGGTATCAACGCAGAG"
CDNA_PCR_FORWARD = SLIDE_READ1
CDNA_PCR_REVERSE = SECOND_STRAND_PRIMER
ADAPTER_TOP = il.INDEX1_PRIMER
ADAPTER_BOTTOM = il.STEM_COMPLEMENT + "T"
INDEX_I5_PREFIX = il.P5
INDEX_I5_SUFFIX = "ACACTCTTTCCCTACACGACGCTC"
INDEX_I7_PREFIX = il.P7
INDEX_I7_SUFFIX = "GTGACTGGAGTTCAGACGTGT"
SPATIAL = feature("spatial_barcode", "spatial_barcode", "whitelist", whitelist="10x-visium-spot-barcode-set")
UMI_FEATURE = feature("umi", "umi", "random")
I5_FEATURE = feature("sample_index_i5", "sample_index", "unknown")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
READ_LENGTHS = {"Read 1":28, "Index 1 (i7)":10, "Index 2 (i5)":10, "Read 2":90}
SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
SEQUENCING_INTRO = "Rev F: R1 28 cycles (spatial barcode16, UMI12), i7/i5 each10, R2 90 cycles (transcript)."
SEQUENCING_ENDING = "Exact spot and TT well members are run-specific; the guide prints the complete pattern."


def capture_parts():
    return [Segment("partial Read 1", SLIDE_READ1, "r1"),
            Segment("spatial barcode", "N"*16, "cbc", placeholder=True, length_bp=16, feature=SPATIAL),
            Segment("UMI", "N"*12, "umi", placeholder=True, length_bp=12, feature=UMI_FEATURE),
            Segment("poly(dT)30", "T"*30), Segment("anchoring V", "V", placeholder=True, length_bp=1, bottom="b"), Segment("anchoring N", "N", placeholder=True, length_bp=1, bottom="n")]


def full_cdna_parts():
    return [*capture_parts(), Segment("captured cDNA", "X"*24, placeholder=True, note="Variable transcript; glyph count is drawing width, not length"), Segment("untemplated CCC", "CCC"),
            Segment("copied TSO", revcomp(TSO[:-6]), "tso")]


def reversed_parts(parts, *, inferred=False):
    return [replace(s, name=s.name+" complement", top=(s.bottom_text() if s.name in ("anchoring V", "anchoring N") else s.top) if s.placeholder else revcomp(s.top),
                    bottom=None, inferred=inferred or s.inferred) for s in reversed(parts)]


def tso_state():
    return Construct([Segment("DNA part", TSO[:-6], "tso"),
                      Segment("terminal rGrGrG", "GGG", "tso", note="Three RNA guanosines")], name="PN-3000228 TSO")


def printed(name, *parts):
    return Construct([Segment(part, bases) for part, bases in parts], name=name)


def capture_scene():
    sc = Scene(); sc.strand("RNA", [role("transcript"), Segment("poly(A) docking span", "A"*30)], label="tissue mRNA")
    sc.anneal("primer", capture_parts(), to="RNA", pair=("poly(dT)30", "poly(A) docking span"),
              label="source slide-primer pattern", unpaired=("anchoring V", "anchoring N"))
    sc.arrow("primer", "on-slide reverse transcription")
    return sc


def adapter_scene():
    sc = Scene(); sc.strand("top", [Segment("Read 2 adapter", ADAPTER_TOP, "r2")], label="PN-2000094 top")
    sc.anneal("bottom", [Segment("paired stem", ADAPTER_BOTTOM[:-1], "r2"), Segment("3′ T", "T", "r2")],
              to="top", pair=("paired stem", "Read 2 adapter"), label="PN-2000094 bottom")
    sc.ligation_end("bottom", "3′ T", kind="cohesive", polarity="3-prime", owner="adapter", text="source 3′-T overhang")
    return sc


def pre_index_parts():
    return [*capture_parts(), Segment("captured cDNA", "X"*24, placeholder=True, note="Variable transcript; glyph count is drawing width, not length"), Segment("dA junction", "A"),
            Segment("adapter paired stem", ADAPTER_TOP[:12], "r2"),
            Segment("adapter fork arm", ADAPTER_TOP[12:], "r2")]


def pre_index_scene():
    # Source p71 ligation product has a long top arm and only the short adapter
    # strand below it. It has no P5, P7 or sample index before index PCR.
    sc = Scene(); sc.strand("top", pre_index_parts(), label="adapter-ligated cDNA, before indexing")
    bottom = [Segment("bottom adapter stem", ADAPTER_BOTTOM[:-1], "r2"), Segment("bottom dT junction", "T"),
              *reversed_parts([*capture_parts(), Segment("captured cDNA", "X"*24, placeholder=True, note="Variable transcript; glyph count is drawing width, not length")])]
    sc.anneal("bottom", bottom, to="top", pair=("bottom adapter stem", "adapter paired stem"), label="source short adapter strand")
    sc.junction("top", "dA junction", "adapter paired stem", "adapter ligation")
    return sc


def index_primer(which):
    if which=="i5":
        return Construct([Segment("P5 arm", INDEX_I5_PREFIX, "p5"),
            Segment("i5", "N"*10, "cbc", placeholder=True, length_bp=10, feature=I5_FEATURE),
            Segment("i5 3′ footprint", INDEX_I5_SUFFIX, "r1")], name="TT i5 primer")
    return Construct([Segment("P7 arm", INDEX_I7_PREFIX, "p7"),
        Segment("i7", "N"*10, "cbc", placeholder=True, length_bp=10, feature=I7_FEATURE),
        Segment("i7 3′ footprint", INDEX_I7_SUFFIX, "r2")], name="TT i7 primer")


def final_library():
    con = Construct([Segment("P5", il.P5, "p5"),
        Segment("i5", "N"*10, "cbc", placeholder=True, length_bp=10, feature=I5_FEATURE),
        Segment("Read 1 arm", il.TRUSEQ_READ1, "r1"), *capture_parts()[1:], Segment("captured cDNA", "X"*24, placeholder=True, note="Variable transcript; glyph count is drawing width, not length"),
        Segment("dA junction", "A"), Segment("Read 2 arm", ADAPTER_TOP, "r2"),
        Segment("i7 reverse complement", "N"*10, "cbc", placeholder=True, length_bp=10, feature=I7_FEATURE),
        Segment("P7 reverse complement", il.P7_RC, "p7")], name="CG000239 Rev F indexed Visium library")
    problems = sp.verify(con, SEQ_PRIMERS, required_roles=tuple(p.role for p in SEQ_PRIMERS))
    if problems: raise ValueError("invalid source Visium library: "+"; ".join(problems))
    return con

FINAL_LIBRARY = final_library()
FINAL_CAPTION = "Rev F Appendix p71 indexed library pattern; the pre-index molecule lacks its P5/P7/index additions."


def oligos():
    return [oligo("slide capture primer pattern", capture_parts()),
        "template-switch oligo: 5′-"+TSO+"-3′",
        oligo("second-strand/cDNA PCR reverse", [Segment("primer", SECOND_STRAND_PRIMER)]),
        oligo("cDNA PCR forward", [Segment("primer", CDNA_PCR_FORWARD)]),
        oligo("ligation adapter top", [Segment("adapter", ADAPTER_TOP)]),
        oligo("ligation adapter bottom, ordering orientation", [Segment("adapter", ADAPTER_BOTTOM)]),
        oligo("TT i5 primer", list(index_primer("i5"))), oligo("TT i7 primer", list(index_primer("i7")))]


def _cleanup(g, key, template, label, *, final=False, double=False):
    beads=g.input(key+"_beads", label, MaterialKind.SELECTION, source=SOURCE)
    kept=g.product(key+"_product", label+" retained cDNA", MaterialKind.LIBRARY if final else MaterialKind.MOLECULAR, state=template.state)
    wastes=[g.product(key+"_small", "removed supernatant/short-product fraction", MaterialKind.WASTE)]
    if double: wastes.append(g.product(key+"_large", "first-cut bead fraction not transferred with supernatant", MaterialKind.WASTE))
    g.react(key, label, StepKind.CLEANUP, ((template,InputRole.TEMPLATE),(beads,InputRole.SELECTION)), (kept,*wastes), changes_molecule=False)
    for waste in wastes:g.discard(waste,"physically excluded by the source bead/supernatant transfers")
    return kept


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    tissue=g.input("tissue","imaged tissue section",MaterialKind.BIOLOGICAL,state=molecule("tissue",(role("local poly(A) RNA"),)),primary=True)
    permeabilization=g.input("permeabilization","Rev F tissue permeabilization",MaterialKind.CONDITION)
    released=g.product("released_rna","locally released RNA",MaterialKind.MOLECULAR,state=molecule("RNA",(role("transcript"),role("poly(A)"))))
    g.react("permeabilize","Permeabilize tissue over the array",StepKind.EXTRACTION,((tissue,InputRole.TEMPLATE),(permeabilization,InputRole.CONDITION)),(released,))
    surface=g.input("surface_oligo","source slide-primer pattern; attachment chemistry unassigned",MaterialKind.OLIGO,
                    state=molecule("capture primer",tuple(capture_parts())),source=SOURCE+" Appendix p70")
    tso=g.input("tso","PN-3000228 template-switch oligo",MaterialKind.OLIGO,state=tso_state(),source=SOURCE+" Appendix p70")
    rt=g.input("rt","on-slide reverse transcription mix",MaterialKind.CATALYST)
    first_state=molecule("first-strand cDNA on slide",tuple(full_cdna_parts()))
    first=g.product("surface_cdna","single-stranded spatially barcoded cDNA",MaterialKind.MOLECULAR,state=first_state)
    g.react("reverse_transcribe","Reverse-transcribe and template-switch on the slide",StepKind.REVERSE_TRANSCRIPTION,
            ((released,InputRole.TEMPLATE),(surface,InputRole.STRUCTURAL),(tso,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(first,))
    second_primer=g.input("second_primer","PN-2000217 second-strand primer",MaterialKind.OLIGO,state=printed("second-strand primer",("primer",SECOND_STRAND_PRIMER)),source=SOURCE+" Appendix p70")
    second_pol=g.input("second_polymerase","Second Strand Enzyme; 65 °C15 min",MaterialKind.CATALYST)
    duplex_state=Scene.duplex(full_cdna_parts())
    duplex=g.product("surface_duplex","double-stranded cDNA on slide",MaterialKind.MOLECULAR,state=duplex_state)
    g.react("second_strand","Synthesize the second cDNA strand on the slide",StepKind.EXTENSION,
            ((first,InputRole.TEMPLATE),(second_primer,InputRole.STRUCTURAL),(second_pol,InputRole.CATALYST)),(duplex,))
    koh=g.input("koh","0.08 M KOH10 minRT; transfer into1 M Tris pH7.0",MaterialKind.CONDITION)
    second_state=molecule("INFERRED released second strand",tuple(reversed_parts(full_cdna_parts(),inferred=True)))
    cdna=g.product("released_cdna","INFERRED second-strand cDNA transferred from slide",MaterialKind.MOLECULAR,
                   state=second_state,inferred=True,source=SOURCE+" Step2.2 transfer; strand identity conditional on unassigned slide attachment")
    residue=g.product("slide_residue","discarded slide and retained material",MaterialKind.WASTE)
    g.react("denature","INFERRED released strand — denature and transfer cDNA",StepKind.STRAND_SEPARATION,
            ((duplex,InputRole.TEMPLATE),(koh,InputRole.CONDITION)),(cdna,residue))
    g.discard(residue,"Step2.2h permits discarding slide/cassette after transfer; exact retained strand chemistry unresolved")
    aliquoting=g.input("qpcr_aliquot_condition","Step3.1:1 µL qPCR aliquot; ~35 µL retained for amplification",MaterialKind.CONDITION)
    main=g.product("main_cdna","cDNA retained for production amplification",MaterialKind.MOLECULAR,state=cdna.state)
    qc=g.product("qpcr_aliquot","qPCR aliquot for Cq determination",MaterialKind.QC,state=cdna.state)
    g.react("qpcr_sampling","Take the source qPCR aliquot",StepKind.SAMPLING,((cdna,InputRole.TEMPLATE),(aliquoting,InputRole.CONDITION)),(main,qc),changes_molecule=False)
    g.mark_terminal(Terminal(qc.id,Disposition.QC,"Step3.1 KAPA SYBR FAST qPCR/cDNA Primers; rounded Cq sets production PCR cycles. Molecular QC reaction not reconstructed."))
    wf=g.input("wta_forward","PN-2000089 cDNA forward",MaterialKind.OLIGO,state=printed("forward",("primer",CDNA_PCR_FORWARD)),source=SOURCE+" Appendix p70")
    wr=g.input("wta_reverse","PN-2000089 cDNA reverse",MaterialKind.OLIGO,state=printed("reverse",("primer",CDNA_PCR_REVERSE)),source=SOURCE+" Appendix p70")
    wp=g.input("wta_polymerase","Amp Mix PN-2000047/2000103; cycles rounded from Cq",MaterialKind.CATALYST)
    amp=g.product("amplified","amplified double-stranded spatial cDNA",MaterialKind.MOLECULAR,state=duplex_state)
    g.react("wta","Amplify the transferred cDNA with the printed primers",StepKind.PCR,
            ((main,InputRole.TEMPLATE),(wf,InputRole.STRUCTURAL),(wr,InputRole.STRUCTURAL),(wp,InputRole.CATALYST)),(amp,),
            geometry=PCRGeometry.unavailable(main.state,reason="The actual template is the conditional transferred single strand. "
                "Printed primer identities/scaffolds are known; native first-cycle geometry is not declared by substituting the amplified duplex."))
    clean=_cleanup(g,"cdna_cleanup",amp,"Step3.3 SPRIselect cleanup")
    fragmentase=g.input("fragmentation","Step4.1 Fragmentation Mix/Enzyme",MaterialKind.CATALYST)
    fragments=g.product("fragments","fragmented cDNA population (end classes unresolved)",MaterialKind.MOLECULAR,
                         state=molecule("fragment population",(role("variable cDNA fragments; barcode-retention classes not reconstructed"),)))
    g.react("fragment","Fragment amplified cDNA",StepKind.FRAGMENTATION,((clean,InputRole.TEMPLATE),(fragmentase,InputRole.CATALYST)),(fragments,))
    endmix=g.input("endprep","Step4.1 end-repair/A-tailing components in same tube",MaterialKind.CATALYST)
    tailed=g.product("tailed","end-repaired dA-tailed fragments",MaterialKind.MOLECULAR,
                      state=molecule("dA-tailed population",(role("cDNA duplex fragments with3′dA ends"),)))
    g.react("end_prep","End-repair and A-tail in the Step4.1 reaction",StepKind.END_REPAIR,((fragments,InputRole.TEMPLATE),(endmix,InputRole.CATALYST)),(tailed,))
    sized=_cleanup(g,"frag_size_select",tailed,"Step4.2 double-sided SPRIselect selection",double=True)
    adapter=g.input("adapter","PN-2000094 top/bottom adapter",MaterialKind.ADAPTER,state=adapter_scene(),source=SOURCE+" Appendix p71")
    ligase=g.input("ligase","Step4.3 DNA ligase",MaterialKind.CATALYST)
    adapted=g.product("adapted","source barcode-bearing adapter-ligated precursor (no sample indexes)",MaterialKind.MOLECULAR,state=pre_index_scene())
    g.react("ligate","Ligate source adapters; show barcode-bearing precursor",StepKind.LIGATION,
            ((sized,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(adapted,),
            geometry=AdapterLigationGeometry.known("cohesive",adapter_scene(),statement="RevF Appendix p71 prints the3′T adapter and ligated product. Fragment-class fates remain unresolved."))
    ligclean=_cleanup(g,"ligation_cleanup",adapted,"Step4.4 post-ligation SPRIselect cleanup")
    forward=g.input("visium_forward","Dual Index TT i5 primer pattern",MaterialKind.OLIGO,state=index_primer("i5"),source=SOURCE+" Appendix p71")
    reverse=g.input("visium_reverse","Dual Index TT i7 primer pattern",MaterialKind.OLIGO,state=index_primer("i7"),source=SOURCE+" Appendix p71")
    pol=g.input("visium_polymerase","Step4.5 Amp Mix; cycles input-dependent per guide table",MaterialKind.CATALYST)
    indexed=g.product("visium_amplified","sample-indexed Visium PCR product",MaterialKind.MOLECULAR,state=FINAL_LIBRARY)
    g.react("visium_pcr","Add P5/P7 and both10-base sample indexes by PCR",StepKind.PCR,
            ((ligclean,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(indexed,),
            geometry=PCRGeometry.unavailable(ligclean.state,reason="Actual input is the source pre-index forked molecule, without P5/P7/indexes. "
                "Both printed primer patterns are supplied; native fork/first-cycle binding geometry is not implemented by the duplex PCR validator."))
    library=_cleanup(g,"index_size_select",indexed,"Step4.6 double-sided SPRIselect0.6X/0.8X;35 µL elution",final=True,double=True)
    g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,primers=SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in SEQ_PRIMERS),read_lengths=READ_LENGTHS))
    g.declare_reviewed_scope(("tissue","permeabilization","permeabilize","surface_oligo","tso","rt","reverse_transcribe",
        "second_primer","second_polymerase","second_strand","koh","denature","qpcr_aliquot_condition","qpcr_sampling","wta_forward","wta_reverse","wta_polymerase","wta",
        "cdna_cleanup_beads","cdna_cleanup","fragmentation","fragment","endprep","end_prep","frag_size_select_beads","frag_size_select",
        "adapter","ligase","ligate","ligation_cleanup_beads","ligation_cleanup","visium_forward","visium_reverse","visium_polymerase","visium_pcr",
        "index_size_select_beads","index_size_select","index_size_select_product"),SOURCE)
    g.annotate(DetailView("reverse_transcribe","Slide-primer docking","Sequence pattern is printed; slide attachment chemistry remains unassigned.",state=capture_scene()))
    return g.require_publishable()
