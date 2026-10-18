"""Published standard plate RamDA-seq; source-owned reactions and kit boundaries."""
from chemdraw import Construct, Scene, Segment, complement_segments, feature, oligo, revcomp
from protocol_model import role
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, Highlight, InputRole, MaterialKind,
                            ProtocolGraph, SequencingPlan, StepKind)
import seqprimers as sp

TITLE="RamDA-seq — single-cell full-length total RNA"
NOTES="01_ramda-seq.html"
SOURCE='Hayashi et al., <a href="https://doi.org/10.1038/s41467-018-02866-0">Nature Communications (2018)</a>, published Methods standard plate preparation, Supplementary Note1 and Supplementary Data3.'
SUMMARY="Not-so-random and oligo(dT) primers initiate a coupled RT/nicking/displacement reaction. Second-strand synthesis, bead recovery and scaled Nextera XT preparation produce the standard plate library."
CAVEAT="This graph selects the original standard plate method:13 library PCR cycles. Advanced plate and C1 are separate variants. Exact historical Nextera XT loaded strands, native first-cycle assembly, selected index members and physical sequencing-primer sites remain unassigned here. INFERRED native end roles do not substitute a canonical completed library for the PCR input."
FIRST_NSR_1="AAAGCA"
SECOND_NSR_1=revcomp(FIRST_NSR_1)
I5=feature("ramda_xt_i5","sample_index","unknown",note="Native index position/width/selected member unassigned")
I7=feature("ramda_xt_i7","sample_index","unknown",note="Native index position/width/selected member unassigned")
FINAL_LIBRARY=Construct([role("native Nextera XT end region; arrangement unassigned",inferred=True,molecular_feature=I5),role("RamDA cDNA fragment"),role("opposite native Nextera XT end region; arrangement unassigned",inferred=True,molecular_feature=I7)],name="standard plate RamDA library")
SEQ_PRIMERS=()
FINAL_CAPTION="INFERRED native XT end roles — standard plate RamDA library; barcode/index arrangements and physical primer sites unassigned."
SEQUENCING_ENDING="Published standard plate libraries were pooled200fmol/well and sequenced single-read76cycles on NextSeq500. Exact physical read/index sites and index acquisition/cycles remain unassigned in this source scope."
READ_SPECS=(sp.read_spec("Read 1",76,"cDNA fragment; physical entry site unassigned","include"),sp.read_spec("Index 1 (i7)","unassigned","Index acquisition/cycles unassigned","undetermined"),sp.read_spec("Index 2 (i5)","unassigned","Index acquisition/cycles unassigned","undetermined"))


def oligos():
    return (oligo("1st_MmNSR_1 (ordering example, not complete pool)",[Segment("first NSR",FIRST_NSR_1)]),oligo("2nd_MmNSR_1 (ordering example)",[Segment("second NSR",SECOND_NSR_1)]),oligo("oligo(dT)18",[Segment("dT18","T"*18)]))


def population(*names):
    sc=Scene()
    for i,name in enumerate(names):sc.strand('population_'+str(i),[role(name)],label=name)
    return sc


def primer_pool(name,example):
    sc=Scene();sc.strand('example',[Segment(name,example)],label=name+' source ordering example')
    sc.strand('remaining_pool',[role('other distinct members of source408-member pool')],label='independent primer species; complete pool not reproduced')
    return sc


def duplex():
    body=[role('second-strand-copied RamDA cDNA')]
    return Scene.duplex(body,unpaired=tuple(s.name+"'" for s in body),label='double-stranded cDNA')


def second_mixture():
    sc=duplex()
    for name in ('sample RNA templates','ERCC RNA templates','carried genomic-DNA/nuclease fragments'):
        sc.strand(name,[role(name)],label='independent carried '+name+'; later partition unassigned')
    return sc


def cleanup(g,key,template,reagent,state=None,library=False):
    selection=g.input(key+'_selection',reagent,MaterialKind.SELECTION)
    out=g.product(key+'_retained','physically recovered '+template.label,MaterialKind.LIBRARY if library else MaterialKind.MOLECULAR,state=state if state is not None else template.state)
    waste=g.product(key+'_waste','removed supernatant/washes',MaterialKind.WASTE)
    g.react(key,'Recover source-selected DNA fraction',StepKind.CLEANUP,((template,InputRole.TEMPLATE),(selection,InputRole.SELECTION)),(out,waste),changes_molecule=False)
    g.discard(waste,'Source bead supernatant/wash transfer; identities of individual RNA/hydrolysis species in those fractions remain unassigned.')
    return out


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    sample=g.input('sample','1µl standard plate single-cell lysate or10pgRNA',MaterialKind.BIOLOGICAL,state=population('total RNA','possible contaminating genomic DNA'),primary=True)
    heat=g.input('initial_heat','70°C90s; hold4°C',MaterialKind.CONDITION)
    denatured=g.product('denatured','heat-treated standard plate sample',MaterialKind.MOLECULAR,state=sample.state)
    g.react('denature','Denature sample before genomic-DNA digestion',StepKind.THERMAL_TREATMENT,((sample,InputRole.TEMPLATE),(heat,InputRole.CONDITION)),(denatured,))
    dnase=g.input('dnase','DNase I Amplification Grade; initial30°C5min, retained in subsequent RT mixture',MaterialKind.CATALYST)
    ercc=g.input('ercc','ERCC RNA spike-in,1:5000000 in gDNA digestion mix',MaterialKind.BIOLOGICAL,state=population('independent ERCC RNA molecules'))
    cleaned=g.product('cleaned','RNA reaction mixture after contaminating-DNA digestion',MaterialKind.MOLECULAR,state=population('sample RNA','independent ERCC RNA','digested genomic-DNA fragments retained in reaction mixture'))
    g.react('digest_gdna','Digest contaminating genomic DNA and add ERCC spike-in',StepKind.DIGESTION,((denatured,InputRole.TEMPLATE),(dnase,InputRole.CATALYST),(ercc,InputRole.SPIKE_IN)),(cleaned,))
    nsr=g.input('first_nsr','source408-member first-strand NSR pool',MaterialKind.OLIGO,state=primer_pool('1st_MmNSR_1',FIRST_NSR_1))
    dt=g.input('oligo_dt','oligo(dT)18',MaterialKind.OLIGO,state='T'*18)
    rt=g.input('rt','RNase-H-minus PrimeScript enzyme mix',MaterialKind.CATALYST)
    gp32=g.input('gp32','T4 gene32 single-strand-binding protein',MaterialKind.CONDITION)
    program=g.input('rt_program','25°C10min,30°C10min,37°C60min,50°C5min,94°C5min',MaterialKind.CONDITION)
    amplified=g.product('ramda_product','coupled RT-RamDA product mixture',MaterialKind.MOLECULAR,state=population('retained sample RNA templates','independent ERCC RNA templates','nicked cDNA strands; RNA-hybrid pairing unassigned','independent displaced single-strand cDNA copies','carried genomic-DNA/nuclease product populations'))
    g.react('rt_ramda','Coupled reverse transcription, cDNA nicking and displacement',StepKind.REVERSE_TRANSCRIPTION,((cleaned,InputRole.TEMPLATE),(nsr,InputRole.STRUCTURAL),(dt,InputRole.STRUCTURAL),(rt,InputRole.CATALYST),(dnase,InputRole.CATALYST),(gp32,InputRole.CONDITION),(program,InputRole.CONDITION)),(amplified,))
    second=g.input('second_nsr','source complementary second-strand NSR pool',MaterialKind.OLIGO,state=primer_pool('2nd_MmNSR_1',SECOND_NSR_1))
    kl=g.input('klenow','Klenow Fragment3′→5′exo−; ramp4→37°C1°C/min,37°C30min',MaterialKind.CATALYST)
    ds=g.product('ds_cdna','second-strand-copied RamDA cDNA plus carried reaction components',MaterialKind.MOLECULAR,state=second_mixture())
    g.react('second_strand','Synthesize complementary cDNA in standard plate conditions',StepKind.EXTENSION,((amplified,InputRole.TEMPLATE),(second,InputRole.STRUCTURAL),(kl,InputRole.CATALYST)),(ds,))
    recovered=cleanup(g,'pre_tagment',ds,'15µl AMPureXP; washed bead-bound cDNA directly eluted in6µl1×TagmentDNA Buffer',duplex())
    loaded=g.input('xt_transposome','historical Nextera XT loaded transposomes; native strands unassigned',MaterialKind.ADAPTER)
    tn5=g.input('xt_tn5','Nextera XT Tn5;2/5volume standard plate preparation',MaterialKind.CATALYST)
    tagged=g.product('tagged','actual pre-PCR tagmented cDNA strand classes',MaterialKind.MOLECULAR,state=population('native XT-tagged cDNA; transferred strands, gaps and complete ends unassigned'))
    other=g.product('other_ends','other tagmentation end classes',MaterialKind.MOLECULAR)
    g.react('tagment','Tagment recovered double-stranded RamDA cDNA',StepKind.TAGMENTATION,((recovered,InputRole.TEMPLATE),(loaded,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,other))
    g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Amplification compatibility is not physical removal; individual recovery through subsequent beads unassigned.'))
    f=g.input('xt_forward','native XT indexing primer first binding role; member/strand unassigned',MaterialKind.OLIGO)
    r=g.input('xt_reverse','native XT indexing primer opposite binding role; member/strand unassigned',MaterialKind.OLIGO)
    pol=g.input('xt_pcr','kit PCR polymerase;13cycles standard plate',MaterialKind.CATALYST)
    amp=g.product('amplified_library','standard plate indexed library',MaterialKind.MOLECULAR,state=FINAL_LIBRARY,inferred=True)
    g.react('index_pcr','Complete native XT library with13cycles',StepKind.PCR,((tagged,InputRole.TEMPLATE),(f,InputRole.STRUCTURAL),(r,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(tagged.state,reason='Actual pre-PCR tagged strands, not final duplex, enter PCR. Historical XT primers/first-cycle geometry are source-delegated and unassigned here.'))
    library=cleanup(g,'post_pcr',amp,'1.2×AMPureXP;24µlTE elution',library=True)
    g.sequence(library,SequencingPlan(library=FINAL_LIBRARY,acquired_roles=('Read 1',),read_lengths={'Read 1':76},limitation=SEQUENCING_ENDING))
    g.annotate(Highlight('rt_ramda','One coupled amplification reaction','Supplementary Note1 describes RT, DNaseI nicking and gp32-assisted strand displacement within the same reaction; graph ordering does not invent separate nicking/RT reagent additions. NativeRNA and nuclease products are not physically discarded by catalysis.'))
    g.annotate(Highlight('index_pcr','Other source variants remain separate','Advanced plate uses14PCRcycles, altered gp32/incubation and second-strand conditions; C1 is a separate IFC method. Neither is merged into this standard13cycle recipe. The printed C1 gDNA-digestion130°C wording is outside this selected variant and remains a deferred original-PDF check.'))
    g.declare_reviewed_scope(('sample','initial_heat','denature','dnase','ercc','digest_gdna','first_nsr','oligo_dt','rt','gp32','rt_program','rt_ramda','second_nsr','klenow','second_strand','pre_tagment_selection','pre_tagment','xt_transposome','xt_tn5','tagment','xt_forward','xt_reverse','xt_pcr','index_pcr','post_pcr_selection','post_pcr','post_pcr_retained'),SOURCE)
    return g.require_publishable()
