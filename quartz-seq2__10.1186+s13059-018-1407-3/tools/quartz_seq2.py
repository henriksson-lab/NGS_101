"""Molecular construct model for Quartz-Seq2 (Sasagawa et al., 2018)."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp, complement_segments

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="quartz-seq2-v3.2-rt-primer-set")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE, "i7": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


M = "TATAGAATTCGCGGCCGCTCGCGAT"
TAGGING_PRIMER = M + "T" * 23
GM_PRIMER = "G" + M
V32_BARCODE = "ATCAATCCTATCTGC"       # eMDRT0001
RT_PRIMER_1536_FIRST = M + "AC" + V32_BARCODE + "N" * 8 + "T" * 24
RYSHAPE_P5 = "GATCGGAAGAGCGTCGTGTA"
LT06_OLIGO_INDEX = "ATTGGC"
RYSHAPE_P7_LT06 = il.P7 + LT06_OLIGO_INDEX + il.TRUSEQ_READ2
P5_GMAC = il.P5[:28] + "TT" + GM_PRIMER + "AC"
READ1_DROPQUARTZ = "ACATT" + GM_PRIMER + "AC"
TPC2 = il.P7[:22]


def rt_primer(barcode_nt: int = 15) -> list[Segment]:
    return [seg("M handle", M, "tso"), seg("AC", "AC"),
            seg("cell barcode", V32_BARCODE if barcode_nt == 15 else "B" * barcode_nt, "cbc", placeholder=barcode_nt != 15),
            seg("UMI", "N" * 8, "umi", placeholder=True), seg("dT24", "T" * 24)]


def tagging_primer() -> list[Segment]:
    return [seg("M handle", M, "tso"), seg("dT23", "T" * 23)]


def mrna_rt_scene() -> Scene:
    mrna = [seg("RNA body", "X" * 28, placeholder=True), seg("poly(A)", "A" * 24)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT primer", rt_primer(), to="mRNA", pair=("dT24", "poly(A)"),
              label="barcoded RT primer")
    sc.arrow("RT primer", "reverse transcriptase")
    return sc


def first_strand() -> Construct:
    return Construct([
        *rt_primer(), seg("cDNA", "X" * 28, placeholder=True),
    ], name="Quartz-Seq2 barcoded first-strand cDNA")


def tailed_first_strand() -> Construct:
    return Construct([
        *first_strand(),
        seg("TdT poly(A)", "[variable TdT-added poly(dA)]", placeholder=True),
    ], name="poly(A)-tailed Quartz-Seq2 first strand")


def tagging_scene() -> Scene:
    target = tailed_first_strand()
    sc = Scene()
    sc.strand("first strand", list(target), label="first strand")
    sc.anneal("tagging primer", tagging_primer(), to="first strand",
              pair=("dT23", "TdT poly(A)"), label="tagging primer", unpaired=("dT23", "M handle"))
    sc.arrow("tagging primer", "second-strand synthesis")
    return sc


def wta_amplicon() -> Construct:
    """Long cDNA product with the single gM primer facing inward at both ends."""
    return Construct([
        seg("tagging-end gM", GM_PRIMER, "tso"),
        seg("tagging poly(T)", "T" * 23),
        seg("cDNA", "X" * 32, placeholder=True),
        seg("RT-primer poly(A)", "A" * 24),
        seg("UMI reverse complement", "N" * 8, "umi", placeholder=True, feature=UMI_FEATURE),
        seg("cell barcode reverse complement", revcomp(V32_BARCODE), "cbc", feature=CELL_BARCODE),
        seg("AC reverse complement", "GT"),
        seg("RT-end gM reverse complement", revcomp(GM_PRIMER), "tso"),
    ], name="Quartz-Seq2 suppression-amplified cDNA")


def truncated_adapter_scene() -> Scene:
    """Source-listed truncated Y adapter; Scene enforces its 12-base stem."""
    long = [seg("P7", il.P7, "p7"), seg("i7", LT06_OLIGO_INDEX, "cbc"),
            seg("Read 2 arm", il.TRUSEQ_READ2[:-13], "r2"),
            seg("stem complement", il.TRUSEQ_READ2[-13:-1], "r2"), seg("T", "T")]
    short = [seg("stem", RYSHAPE_P5[:12], "r1"),
             seg("unpaired Read 1 arm", RYSHAPE_P5[12:], "r1")]
    sc = Scene()
    sc.strand("long", long, label="P7 strand")
    sc.anneal("short", short, to="long", pair=("stem", "stem complement"),
              label="short strand", unpaired=("unpaired Read 1 arm",))
    sc.overhang("long", "T", polarity="3-prime", owner="adapter",
                text="adapter 3′-T ligation overhang")
    return sc


READ1_PRIMER = sp.custom("Read 1", "Read1DropQuartz", READ1_DROPQUARTZ,
                         "Sasagawa et al. 2018, Supplementary Table S4")
SEQ_PRIMERS = (READ1_PRIMER, sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def selective_pcr_template(insert_nt: int = 36) -> Construct:
    """Actual LT-strand ligation product; adapter is not copied onto the RT strand yet."""
    rt_fragment=[seg("gM primer",GM_PRIMER,"tso"),seg("AC","AC"),
                 seg("cell barcode",V32_BARCODE,"cbc"),
                 seg("UMI","N"*8,"umi",placeholder=True),seg("poly(T)","T"*24),
                 seg("cDNA fragment","X"*insert_nt,placeholder=True)]
    return Construct([seg("native rYshapeP7LT06",RYSHAPE_P7_LT06,"p7"),
                      *complement_segments(rt_fragment),seg("original opposite cDNA3′dA","A",note="source end-prep-added dA on original opposite strand; not copied adapter sequence")],
                     name="native LT strand covalently ligated to opposite cDNA5′ end")


def final_library(insert_nt: int = 36) -> Construct:
    lib = Construct([
        seg("truncated P5", il.P5[:28], "p5"),
        seg("TT", "TT"),
        seg("gM primer", GM_PRIMER, "tso"),
        seg("AC", "AC"),
        seg("cell barcode", V32_BARCODE, "cbc"),
        seg("UMI", "N" * 8, "umi", placeholder=True),
        seg("poly(T)", "T" * 24),
        seg("3-prime cDNA fragment", "X" * insert_nt, placeholder=True),
        seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7", revcomp(LT06_OLIGO_INDEX), "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="Quartz-Seq2 sequencing library")
    for primer in SEQ_PRIMERS:
        if sp.locate(lib, primer) is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
    return lib


def primer_landings(lib: Construct | None = None):
    lib = lib or final_library()
    out = []
    for primer in SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
        out.append((primer, hit))
    return out


def read1_layout() -> list[tuple[str, str]]:
    return [(x.cycles, x.feature.label) for x in
            sp.feature_spans(final_library(), SEQ_PRIMERS, {"Read 1": 23})]


# Native source-owned process declarations. Independent polymers are Scene rows;
# role labels bound geometry rather than supplying a canonical kit shell.
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (ProtocolGraph,MaterialKind,InputRole,StepKind,SequencingPlan,Terminal,Disposition)


def role(name,text):
    return seg(name,"["+text+"]",placeholder=True)


def native_state(stage):
    sc=Scene()
    if stage=="sample":
        sc.strand("rna",[role("RNA","cellular RNA; polyadenylated transcripts selected by primer")])
        sc.strand("ercc",[role("ERCC","independent ERCC spike-in RNA; not covalently joined to cell RNA")])
        sc.strand("other",[role("other lysate","other cellular polymers; recovery not class-resolved")])
        return sc
    if stage in ("rt","exo","tail"):
        parts=list(first_strand())
        if stage=="tail":parts += [role("TdT tail","variable 3′ poly(dA); no fixed tail length assigned")]
        sc.strand("first",parts,label="actual first-strand DNA")
    elif stage=="second":
        # First cycle has tagging M, not the extra5′G introduced by suppression primer.
        parts=list(wta_amplicon())
        parts[0]=seg("tagging-end M",M,"tso")
        parts[-1]=seg("RT-end M reverse complement",revcomp(M),"tso")
        sc.strand("second",parts,label="actual tagging-primer-extended ss DNA")
        sc.strand("original_first",list(tailed_first_strand()),label="original tailed first strand; pairing extent not fixed")
    elif stage=="wta":
        parts=list(wta_amplicon());sc.strand("wta",parts,label="suppression-PCR daughter",mod5="")
        sc.anneal("wta_complement",complement_segments(parts),to="wta",pair=(parts[0].name+"'",parts[0].name),mod5="")
        sc.strand("originals",[role("originals","possible residual native first/second strands; chemical modifications remain lineage-specific")])
        sc.strand("short_byproducts",[role("short byproducts","possible tailed-primer-derived short DNA; suppression is nonamplification, not physical discard")])
    elif stage in ("fragment","repaired","tailed_fragments"):
        sc.strand("body",[role("fragment","sheared cDNA fragment; random endpoints and possible original WTA termini")],mod5="5′P" if stage in ("repaired","tailed_fragments") else "")
        sc.strand("opposite",[role("opposite","opposite cDNA strand; random endpoints")],mod5="5′P" if stage in ("repaired","tailed_fragments") else "")
        if stage=="tailed_fragments":
            # Both physical strands have3′A, not one linear A appended to a duplex shell.
            sc.strands["body"].segs += [seg("3′A","A")]
            sc.strands["opposite"].segs += [seg("opposite 3′A","A")]
    elif stage=="ligated":
        sc.strand("native_ligated_opposite",list(selective_pcr_template()),label="actual LT3′T→opposite cDNA5′ ligation product; no copied adapter shell")
        sc.strand("native_rt_end",[seg("upstream native LT strand",RYSHAPE_P7_LT06,"p7",note="native LT3′T ligated to RT-end fragment5′; lies outside internal gM selective-PCR boundary"),seg("gM",GM_PRIMER),seg("AC","AC"),seg("cell barcode",V32_BARCODE,"cbc"),seg("UMI","N"*8,"umi",placeholder=True),seg("dT24","T"*24),role("native body","RT-end cDNA fragment"),seg("3′A","A")],label="RT-end fragment with source LT strand upstream of gM; no copied adapter or P5-gMac extension")
        sc.strand("short_adapter",[seg("short rYshapeP5",RYSHAPE_P5)],label="source short adapter strand; covalent5′-ligation status not assigned from absent phosphate notation")
        sc.strand("nonselected",[role("nonselected","internal / tagging-end fragments, incompletely ligated RT-end classes and other adapter products; nonamplification is not removal")])
    elif stage=="library":
        parts=list(final_library());sc.strand("daughter",parts,label="library PCR daughter; primer-derived phosphorothioate bonds retained, no extra nucleotide")
        sc.anneal("daughter_complement",complement_segments(parts),to="daughter",pair=(parts[0].name+"'",parts[0].name))
        sc.strand("native_residuals",[role("native residuals","possible co-carried original ligated molecules, other fragment classes and WTA originals; class-specific clearance unassigned")])
    sc.strand("rna_residual",[role("RNA residual","original RNA / partial RNA hydrolysates; RNaseH activity" if stage not in ("rt","exo") else "original RNA templates; pairing extent not assigned")])
    if stage != "rt":sc.strand("other_well_lineages",[role("other well lineages","independent other-well DNA copies retain their own cell barcodes and N8UMIs through the same native operation; no covalent concatenation")])
    sc.strand("ercc_residual",[role("ERCC residual","independent ERCC RNA/copy lineage; corresponding DNA copies undergo same primer operations")])
    sc.strand("other_residual",[role("other residual","possible co-carried other lysate polymers; no complete class clearance assigned")])
    if stage not in ("rt",):sc.strand("oligo_residual",[role("oligo residual","possible surviving oligo / primer hydrolysates; class-specific fate unassigned")])
    return sc


def protocol_graph():
    source="https://doi.org/10.1186/s13059-018-1407-3; published Methods232–258 and Additionalfile5 TableS4; selected v3.2/ligation-library branch"
    g=ProtocolGraph("Quartz-Seq2 native v3.2 ligation branch",authority=source)
    M,R,K=MaterialKind,InputRole,StepKind
    def inp(i,label,kind=M.CATALYST,state=None):return g.input(i,label,kind,state=state)
    def prod(i,label,stage,kind=M.MOLECULAR):return g.product(i,label,kind,state=native_state(stage))
    def cleanup(i,label,x):
        sel=inp(i+"_medium",label,M.SELECTION)
        y=g.product(i+"_retained",label+" retained mixture",x.kind,state=x.state)
        w=g.product(i+"_waste","physically unretained fraction; composition unresolved",M.WASTE)
        g.react(i,label,K.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(y,w));g.discard(w,"physically unretained fraction; no complete class-specific RNA/oligo/template removal asserted")
        return y
    sample=g.input("sample","single-cell lysate plus independent ERCC RNA",M.BIOLOGICAL,state=native_state("sample"),primary=True)
    primer=inp("rt_primer","source RT primer v3.2 selected eMDRT0001; N8DNA UMI",M.OLIGO,Construct(rt_primer()))
    rt=inp("rt","SuperScriptIII");first=prod("first","native ss cDNA plus independent RNA templates","rt")
    g.react("reverse_transcription","Reverse transcription 35°C5min/50°C50min;70°C15min",K.REVERSE_TRANSCRIPTION,((sample,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST)),(first,))
    peer_state=Scene();peer_state.strand("peer_copies",[role("other-well copies","independent first-strand DNA molecules with their own source well-barcode members and N8UMIs")]);peer_state.strand("peer_rna",[role("peer RNA","independent original cell/ERCC RNA templates and lysate residuals")])
    peers=inp("peer_wells","other independently barcoded wells; not a covalent concatenation",M.MOLECULAR,peer_state)
    pooled=prod("pooled","physically pooled independent well mixtures","rt")
    pooled.state.strand("peer_copies",[role("peer copies","other wells retain their own barcode/UMI identities; no barcode substitution or covalent joining")])
    g.react("pool_wells","Collect pooled wells by inverted-plate centrifugation",K.POOL,((first,R.TEMPLATE),(peers,R.TEMPLATE)),(pooled,),changes_molecule=False)
    current=cleanup("pooled_column","DNA Clean & Concentrator5;20µl elution per column",pooled)
    tdt=inp("tdt","terminal deoxynucleotidyl transferase");rnase=inp("rnase_h","RNaseH; source tailing-mixture activity, not a separate prior incubation");datp=inp("datp","dATP",M.CONDITION)
    tailed=prod("tailed","native poly(dA)-tailed first DNA plus RNA hydrolysates/residuals","tail")
    g.react("tail_and_rnase_h","TdT plus RNaseH in one source mixture:37°C75s;65°C10min",K.TAILING,((current,R.TEMPLATE),(tdt,R.CATALYST),(rnase,R.CATALYST),(datp,R.CONDITION)),(tailed,))
    # Source physically distributes45µl into four approximately11µl reactions.
    aliquots=[]
    for n in range(4):aliquots.append(prod("aliquot_"+str(n+1),"one of four physical TdT-mixture aliquots","tail"))
    g.react("split_four","Distribute45µl into four approximately11µl aliquots",K.SPLIT,((tailed,R.TEMPLATE),),tuple(aliquots),changes_molecule=False)
    tag=inp("tagging_primer","source M-dT23 tagging primer",M.OLIGO,Construct(tagging_primer()));pol=inp("wta_polymerase","MightyAmp DNA polymerase")
    suppression_state=Scene();suppression_state.strand("primer",[seg("suppression primer",GM_PRIMER)],mod5="")
    suppression=inp("suppression_primer","gM primer; no source amine modification assigned",M.OLIGO,suppression_state)
    amps=[]
    for n,aliquot in enumerate(aliquots,1):
        suffix="_"+str(n)
        second=prod("second"+suffix,"actual tagging-primer-extended second DNA plus first originals","second")
        g.react("second_strand"+suffix,"Tagging:98°C130s/40°C1min;0.2°C/s increment to68°C;68°C5min",K.EXTENSION,((aliquot,R.TEMPLATE),(tag,R.STRUCTURAL),(pol,R.CATALYST)),(second,))
        amp=prod("wta"+suffix,"suppression-PCR daughters plus unresolved originals/byproducts","wta")
        g.react("suppression_pcr"+suffix,"Suppression PCR11cycles;98°C10s/65°C15s/68°C5min",K.PCR,((second,R.TEMPLATE),(suppression,R.STRUCTURAL),(pol,R.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(second.state,reason="Actual native second strands and first-strand originals enter single-species suppression PCR; variable tail and pan-like byproduct folding/geometry remain bounded.",symmetric_species=True))
        amps.append(amp)
    pooled_amp=prod("pooled_wta","pooled four PCR mixtures","wta")
    g.react("pool_wta","Pool all four PCR solutions",K.POOL,tuple((a,R.TEMPLATE) for a in amps),(pooled_amp,),changes_molecule=False)
    current=cleanup("wta_column","MinElute column;40µl water elution",pooled_amp)
    current=cleanup("wta_ampure","32µl AMPureXP;32µl final cDNA",current)
    shear=inp("covaris","CovarisLE220;5–10ng amplified cDNA",M.CONDITION)
    fragmented=prod("fragmented","sheared native WTA DNA mixture","fragment")
    g.react("fragment","Physical cDNA shearing",K.FRAGMENTATION,((current,R.TEMPLATE),(shear,R.CONDITION)),(fragmented,))
    current=cleanup("fragment_cleanup","DNA Clean & Concentrator5;10µl water",fragmented)
    mix=inp("endprep_mix","KAPA End repair & A-tailing Enzyme");tailed_frag=prod("tailed_fragments","native end-repaired3′A-tailed fragments","tailed_fragments")
    g.react("endprep","Source combined end-repair/A-tail mix;37°C60min then65°C30min",K.END_REPAIR,((current,R.TEMPLATE),(mix,R.CATALYST)),(tailed_frag,))
    adapter=inp("adapter","source truncated rYshapeP5/rYshapeP7LT06",M.ADAPTER,truncated_adapter_scene());ligase=inp("ligase","KAPA DNA ligase")
    ligated=prod("ligated","native ligated precursor strands, not completed PCR product","ligated")
    g.react("ligate","Adapter ligation20°C15min",K.LIGATION,((tailed_frag,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.known("cohesive",truncated_adapter_scene(),statement="Source strand sequences give opposed adapter3′T / substrate3′A. Native variable-end ligation/fill extent is bounded independently from the eventual amplified endpoint."))
    current=cleanup("ligation_cleanup","18µl AMPureXP per22µl ligation;20µl water",ligated)
    fstate=Construct([seg("P5-gMac",P5_GMAC,note="source two3′-proximal phosphorothioate bonds; no extra nucleotide")]);f=inp("p5_gmac","P5-gMac_hybrid",M.OLIGO,fstate)
    r=inp("tpc2","TPC2 with terminal phosphorothioate bond",M.OLIGO,Construct([seg("TPC2",TPC2,note="source terminal phosphorothioate bond")]))
    libpol=inp("library_polymerase","KAPA HiFi ready mix");library=prod("amplified_library","native primer-selected library daughters plus unresolved originals","library",M.LIBRARY)
    g.react("library_pcr","P5-gMac/TPC2 selective PCR8cycles;98°C45s initial denaturation",K.PCR,((current,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(libpol,R.CATALYST)),(library,),geometry=PCRGeometry.unavailable(current.state,reason="Actual asymmetric native adapter-ligated strands enter PCR, not a completed amplified duplex. Exact first-cycle nick/gap/fill geometry across heterogeneous sheared endpoints remains unassigned."))
    final=cleanup("library_cleanup","40µl AMPureXP per50µl PCR;20–30µl water",library)
    g.sequence(final,SequencingPlan(library=final_library(),primers=SEQ_PRIMERS,acquired_roles=("Read 1","Index 1 (i7)","Read 2"),read_lengths={"Read 1":23,"Index 1 (i7)":6,"Read 2":63},read_uses={"Read 1":"extract","Read 2":"include","Index 1 (i7)":"demultiplex"},limitation="Selected source v3.2 mask23/6/63; v3.1 uses22/6/64–118 and is not merged. Cell barcode is well identity, N8UMI labels original RNA-derived RT molecules; ERCC and cell RNA remain independent."))
    reviewed=('endprep_mix','fragment_cleanup_medium','library_cleanup_medium','library_cleanup_retained','ligation_cleanup_medium','pooled_column_medium','rt','wta_ampure_medium','wta_column_medium','sample','rt_primer','reverse_transcription','peer_wells','pool_wells','pooled_column','tdt','rnase_h','datp','tail_and_rnase_h','split_four','tagging_primer','wta_polymerase','suppression_primer','second_strand_1','suppression_pcr_1','second_strand_2','suppression_pcr_2','second_strand_3','suppression_pcr_3','second_strand_4','suppression_pcr_4','pool_wta','wta_column','wta_ampure','covaris','fragment','fragment_cleanup','endprep','adapter','ligase','ligate','ligation_cleanup','p5_gmac','tpc2','library_polymerase','library_pcr','library_cleanup')
    g.declare_reviewed_scope(reviewed,source)
    return g.require_publishable()
