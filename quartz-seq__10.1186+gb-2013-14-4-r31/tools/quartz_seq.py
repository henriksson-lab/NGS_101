"""Molecular construct model for Quartz-Seq (Sasagawa et al., 2013)."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp, complement_segments


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


M = "TATAGAATTCGCGGCCGCTCGCGAT"
RT_MIDDLE = "AATACGACTCACTATAGGGCG"
RT_PRIMER = M + RT_MIDDLE + "T" * 24
TAGGING_PRIMER = M + "T" * 23
SUPPRESSION_PRIMER = "G" + M
TRSU = il.TRUSEQ_P5_FULL
TRSI2_INDEX = "CGATGT"
TRSI2 = il.INDEX1_PRIMER + TRSI2_INDEX + il.P7_RC
TPC1 = il.P5[:21]
TPC2 = il.P7[:22]


def rt_primer() -> list[Segment]:
    return [seg("M handle", M, "tso"), seg("T7/GGCG", RT_MIDDLE),
            seg("dT24", "T" * 24)]


def tagging_primer() -> list[Segment]:
    return [seg("M handle", M, "tso"), seg("dT23", "T" * 23)]


def mrna_rt_scene() -> Scene:
    mrna = [seg("RNA body", "X" * 28, placeholder=True), seg("poly(A)", "A" * 24)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT primer", rt_primer(), to="mRNA", pair=("dT24", "poly(A)"),
              label="RT primer")
    sc.arrow("RT primer", "reverse transcriptase")
    return sc


def first_strand() -> Construct:
    return Construct([
        *rt_primer(), seg("cDNA", "X" * 28, placeholder=True),
    ], name="Quartz-Seq first-strand cDNA")


def tailed_first_strand() -> Construct:
    return Construct([
        *first_strand(),
        seg("TdT poly(A)", "[variable TdT-added poly(dA)]", placeholder=True),
    ], name="poly(A)-tailed Quartz-Seq first strand")


def tagging_scene() -> Scene:
    target = tailed_first_strand()
    sc = Scene()
    sc.strand("first strand", list(target), label="first strand")
    sc.anneal("tagging primer", tagging_primer(), to="first strand",
              pair=("dT23", "TdT poly(A)"), label="tagging primer", unpaired=("dT23", "M handle"))
    sc.arrow("tagging primer", "second-strand synthesis")
    return sc


def adapter_scene() -> Scene:
    """Source-printed TruSeq ligating stem and its single 3′-T overhang."""
    bottom = [seg("stem complement", TRSU[-13:-1], "r1"),
              seg("3′ T", TRSU[-1], "r1")]
    top = [seg("stem", TRSI2[:12], "r2")]
    sc = Scene(); sc.strand("TRSU", bottom, label="TRSU ligating end")
    sc.anneal("TRSI-2", top, to="TRSU", pair=("stem", "stem complement"),
              label="5′-phosphorylated TRSI-2", mod5="p")
    sc.overhang("TRSU", "3′ T", polarity="3-prime", owner="adapter",
                text="adapter 3′-T overhang")
    return sc


def wta_amplicon() -> Construct:
    return Construct([
        seg("tagging-end gM", SUPPRESSION_PRIMER, "tso"),
        seg("tagging poly(T)", "T" * 23),
        seg("cDNA", "X" * 32, placeholder=True),
        seg("RT-primer poly(A)", "A" * 24),
        seg("T7/GGCG reverse complement", revcomp(RT_MIDDLE)),
        seg("RT-end gM reverse complement", revcomp(SUPPRESSION_PRIMER), "tso"),
    ], name="Quartz-Seq WTA product")


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def final_library(insert_nt: int = 36) -> Construct:
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("TruSeq Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        seg("cDNA fragment", "X" * insert_nt, placeholder=True),
        seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7", TRSI2_INDEX, "cbc",
            feature=feature("sample_i7", "sample_index", "whitelist",
                            whitelist="published TRSI-2 index set")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="Quartz-Seq sequencing library")
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
        parts=list(wta_amplicon());sc.strand("wta",parts,label="suppression-PCR daughter",mod5="NH2")
        sc.anneal("wta_complement",complement_segments(parts),to="wta",pair=(parts[0].name+"'",parts[0].name),mod5="NH2")
        sc.strand("originals",[role("originals","possible residual native first/second strands; chemical modifications remain lineage-specific")])
        sc.strand("short_byproducts",[role("short byproducts","possible tailed-primer-derived short DNA; suppression is nonamplification, not physical discard")])
    elif stage in ("fragment","repaired","tailed_fragments"):
        sc.strand("body",[role("fragment","internal sheared cDNA fragment; random endpoints")],mod5="5′P" if stage in ("repaired","tailed_fragments") else "")
        sc.strand("opposite",[role("opposite","opposite internal cDNA strand; random endpoints")],mod5="5′P" if stage in ("repaired","tailed_fragments") else "")
        sc.strand("blocked_end_fragment",[role("source WTA end","separate WTA-terminal fragment from NH2-primer-derived strand; 5′NH2 remains after shearing/end repair")],mod5="NH2")
        sc.strand("blocked_end_opposite",[role("opposite blocked-end fragment","opposite terminal-fragment strand; cut/repair endpoint and modification lineage unassigned")])
        if stage=="tailed_fragments":
            # Both physical strands have3′A, not one linear A appended to a duplex shell.
            sc.strands["body"].segs += [seg("3′A","A")]
            sc.strands["opposite"].segs += [seg("opposite 3′A","A")]
    elif stage=="ligated":
        sc.strand("native",[seg("TRSU",TRSU,note="source terminal phosphorothioate bond"),role("body","sheared cDNA"),seg("substrate3′A","A"),seg("TRSI-2",TRSI2,note="source5′phosphate consumed at ligation; terminal phosphorothioate bond retained")],label="actual covalent ligated strand; unfilled Y/nick geometry not completed")
        sc.strand("blocked_terminal_residuals",[role("blocked terminal residuals","possible NH2-bearing original WTA-end fragments and partially ligated end classes; not all end-repaired fragments acquire5′P")])
        sc.strand("opposite",[seg("opposite TRSU",TRSU),role("opposite body","opposite sheared cDNA"),seg("opposite substrate3′A","A"),seg("opposite TRSI-2",TRSI2)],label="opposite native strand; source adapter arms remain asymmetric")
    elif stage=="library":
        parts=list(final_library());sc.strand("daughter",parts,label="library PCR daughter; primer-derived phosphorothioate bonds retained, no extra nucleotide")
        sc.anneal("daughter_complement",complement_segments(parts),to="daughter",pair=(parts[0].name+"'",parts[0].name))
        sc.strand("native_residuals",[role("native residuals","possible co-carried original ligated molecules, other fragment classes and WTA originals; class-specific clearance unassigned")])
    sc.strand("rna_residual",[role("RNA residual","original RNA / partial RNA hydrolysates; RNaseH activity" if stage not in ("rt","exo") else "original RNA templates; pairing extent not assigned")])
    sc.strand("other_residual",[role("other residual","possible co-carried other lysate polymers; no complete class clearance assigned")])
    if stage not in ("rt",):sc.strand("oligo_residual",[role("oligo residual","possible surviving oligo / primer hydrolysates; class-specific fate unassigned")])
    return sc


def protocol_graph():
    source="https://doi.org/10.1186/gb-2013-14-4-r31; published WTA Methods/Table1/LIMprep; RNA-seq branch, not Quartz-Chip IVT"
    g=ProtocolGraph("Quartz-Seq native LIMprep",authority=source)
    M,R,K=MaterialKind,InputRole,StepKind
    def inp(i,label,kind=M.CATALYST,state=None):return g.input(i,label,kind,state=state)
    def prod(i,label,stage,kind=M.MOLECULAR):return g.product(i,label,kind,state=native_state(stage))
    def cleanup(i,label,x):
        sel=inp(i+"_medium",label,M.SELECTION)
        y=g.product(i+"_retained",label+" retained mixture",x.kind,state=x.state)
        w=g.product(i+"_waste","physically unretained fraction; composition unresolved",M.WASTE)
        g.react(i,label,K.CLEANUP,((x,R.TEMPLATE),(sel,R.SELECTION)),(y,w));g.discard(w,"physically unretained fraction; no complete class-specific RNA/oligo/template removal asserted")
        return y
    sample=g.input("sample","single-cell lysate",M.BIOLOGICAL,state=native_state("sample"),primary=True)
    primer=inp("rt_primer","source RT primer",M.OLIGO,Construct(rt_primer()))
    rt=inp("rt","SuperScriptIII");first=prod("first","native ss cDNA plus independent RNA templates","rt")
    g.react("reverse_transcription","Reverse transcription 35°C5min/45°C20min;70°C10min",K.REVERSE_TRANSCRIPTION,((sample,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST)),(first,))
    exo=inp("exoi","ExonucleaseI");current=prod("primer_digested","first-strand mixture plus digested and surviving RT primer","exo")
    g.react("digest_rt_primer","ExonucleaseI37°C30min;80°C10min",K.DIGESTION,((first,R.TEMPLATE),(exo,R.CATALYST)),(current,))
    tdt=inp("tdt","terminal deoxynucleotidyl transferase");rnase=inp("rnase_h","RNaseH; source tailing-mixture activity, not a separate prior incubation");datp=inp("datp","dATP",M.CONDITION)
    tailed=prod("tailed","native poly(dA)-tailed first DNA plus RNA hydrolysates/residuals","tail")
    g.react("tail_and_rnase_h","TdT plus RNaseH in one source mixture:37°C50s;65°C10min",K.TAILING,((current,R.TEMPLATE),(tdt,R.CATALYST),(rnase,R.CATALYST),(datp,R.CONDITION)),(tailed,))
    aliquots=[tailed]
    tag=inp("tagging_primer","source M-dT23 tagging primer",M.OLIGO,Construct(tagging_primer()));pol=inp("wta_polymerase","MightyAmp DNA polymerase")
    suppression_state=Scene();suppression_state.strand("primer",[seg("suppression primer",SUPPRESSION_PRIMER)],mod5="NH2")
    suppression=inp("suppression_primer","5′NH2 suppression primer",M.OLIGO,suppression_state)
    amps=[]
    for n,aliquot in enumerate(aliquots,1):
        suffix="_"+str(n)
        second=prod("second"+suffix,"actual tagging-primer-extended second DNA plus first originals","second")
        g.react("second_strand"+suffix,"Tagging:98°C130s/40°C1min;68°C5min",K.EXTENSION,((aliquot,R.TEMPLATE),(tag,R.STRUCTURAL),(pol,R.CATALYST)),(second,))
        amp=prod("wta"+suffix,"suppression-PCR daughters plus unresolved originals/byproducts","wta")
        g.react("suppression_pcr"+suffix,"Suppression PCR21cycles;98°C10s/65°C15s/68°C5min",K.PCR,((second,R.TEMPLATE),(suppression,R.STRUCTURAL),(pol,R.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(second.state,reason="Actual native second strands and first-strand originals enter single-species suppression PCR; variable tail and pan-like byproduct folding/geometry remain bounded.",symmetric_species=True))
        amps.append(amp)
    current=cleanup("wta_cleanup","MinElute purification (source alternative:AMPureXP)",amps[0])
    shear=inp("covaris","CovarisS220;20ng amplified cDNA",M.CONDITION)
    fragmented=prod("fragmented","sheared native WTA DNA mixture","fragment")
    g.react("fragment","Physical cDNA shearing",K.FRAGMENTATION,((current,R.TEMPLATE),(shear,R.CONDITION)),(fragmented,))
    current=cleanup("fragment_cleanup","DNA Clean & Concentrator5;10µl water",fragmented)
    repair=inp("end_repair_mix","KAPA end-repair enzymes");repaired=prod("repaired","native repaired cDNA fragments","repaired")
    g.react("end_repair","End repair20°C30min",K.END_REPAIR,((current,R.TEMPLATE),(repair,R.CATALYST)),(repaired,))
    current=cleanup("repair_cleanup","DNA Clean & Concentrator5;12.5µl EB1/10",repaired)
    adenylase=inp("a_tail_mix","KAPA A-tailing enzyme");tailed_frag=prod("tailed_fragments","native3′A-tailed cDNA fragments","tailed_fragments")
    g.react("a_tail","A-tail30°C30min",K.TAILING,((current,R.TEMPLATE),(adenylase,R.CATALYST)),(tailed_frag,))
    tailed_frag=cleanup("a_tail_cleanup","DNA Clean & Concentrator5;12.5µl EB1/10",tailed_frag)
    adapter=inp("adapter","source TRSU/TRSI-2; TRSI5′phosphate and terminal S bonds",M.ADAPTER,adapter_scene());ligase=inp("ligase","KAPA DNA ligase")
    ligated=prod("ligated","native ligated precursor strands, not completed PCR product","ligated")
    g.react("ligate","Adapter ligation20°C15min",K.LIGATION,((tailed_frag,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)),(ligated,),geometry=AdapterLigationGeometry.known("cohesive",adapter_scene(),statement="Source strand sequences give opposed adapter3′T / substrate3′A. Native variable-end ligation/fill extent is bounded independently from the eventual amplified endpoint."))
    current=cleanup("ligation_cleanup_1","First high-salt AMPureXP adapter-dimer depletion",ligated)
    current=cleanup("ligation_cleanup_2","Repeat high-salt AMPureXP adapter-dimer depletion;20µl EB1/10",current)
    f=inp("tpc1","TPC1 with terminal phosphorothioate bond",M.OLIGO,Construct([seg("TPC1",TPC1,note="source terminal phosphorothioate bond")]))
    r=inp("tpc2","TPC2 with terminal phosphorothioate bond",M.OLIGO,Construct([seg("TPC2",TPC2,note="source terminal phosphorothioate bond")]))
    libpol=inp("library_polymerase","KAPA HiFi ready mix");library=prod("amplified_library","native primer-selected library daughters plus unresolved originals","library",M.LIBRARY)
    g.react("library_pcr","TPC1/TPC2 PCR10–12cycles;98°C45s initial denaturation",K.PCR,((current,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(libpol,R.CATALYST)),(library,),geometry=PCRGeometry.unavailable(current.state,reason="Actual asymmetric native adapter-ligated strands enter PCR, not a completed amplified duplex. Exact first-cycle nick/gap/fill geometry across heterogeneous sheared endpoints remains unassigned."))
    final=cleanup("library_cleanup","AMPureXP purification after library PCR",library)
    g.sequence(final,SequencingPlan(library=final_library(),primers=SEQ_PRIMERS,acquired_roles=("Read 1","Index 1 (i7)","Read 2"),read_lengths={"Read 1":50,"Index 1 (i7)":None,"Read 2":50},read_uses={"Read 1":"include","Read 2":"include","Index 1 (i7)":"demultiplex"},limitation="Published PE50 comparison cohort selected. Native TRSI-2 index is6nt but this capped pass does not establish acquired index-cycle count; computational WTA-primer trimming is not physical fragment removal. Quartz-Chip IVT branch is not this RNA-seq preparation."))
    reviewed=('rt','wta_cleanup_medium','fragment_cleanup_medium','end_repair_mix','repair_cleanup_medium','a_tail_mix','a_tail_cleanup_medium','ligation_cleanup_1_medium','ligation_cleanup_2_medium','library_cleanup_medium','library_cleanup_retained','sample','rt_primer','reverse_transcription','exoi','digest_rt_primer','tdt','rnase_h','datp','tail_and_rnase_h','tagging_primer','wta_polymerase','suppression_primer','second_strand_1','suppression_pcr_1','wta_cleanup','covaris','fragment','fragment_cleanup','end_repair','repair_cleanup','a_tail','a_tail_cleanup','adapter','ligase','ligate','ligation_cleanup_1','ligation_cleanup_2','tpc1','tpc2','library_polymerase','library_pcr','library_cleanup')
    g.declare_reviewed_scope(reviewed,source)
    return g.require_publishable()
