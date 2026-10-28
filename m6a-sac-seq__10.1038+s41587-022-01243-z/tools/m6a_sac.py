"""Source-owned m6A-SAC-seq chemical branches; modified RNA is not a basecall."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene, Segment, feature, revcomp
import illumina as il
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition
SOURCE='Hu et al.2022 doi:10.1038/s41587-022-01243-z, accepted manuscript Methods211–227; publisher SupplementaryTable1; VOR procedure applicability unverified'
RNA_ADAPTER_WRITTEN='rAppNNNNNATCACGAGATCGGAAGAGCACACGTC/iBiodT//3SpC3/'
CDNA_ADAPTER_WRITTEN='Phos-'+'N'*5+revcomp(il.TRUSEQ_SMALL_RNA_RA5)+'/3SpC3/'
RT_PRIMER=il.TRUSEQ_READ2[-21:]

@dataclass(frozen=True)
class AdenosinePath:
    original: str
    stage: str='native'
    fto_erased: bool=False
    def __post_init__(self):
        if self.original not in {'A','m6A'} or self.stage not in {'native','allyl','cyclized'}:
            raise ValueError('unsupported source adenosine pathway')
        if self.fto_erased and self.original!='m6A':raise ValueError('FTO-erased path must originate at m6A')
    @property
    def substrate(self):return 'A' if self.fto_erased else self.original
    def erase_selected(self):
        if self.original!='m6A' or self.stage!='native' or self.fto_erased:raise ValueError('selected FTO path requires native m6A')
        return replace(self,fto_erased=True)
    def allylate(self):
        if self.stage!='native':raise ValueError('only native A/m6A substrate enters this labeling event')
        return replace(self,stage='allyl')
    def cyclize(self):
        if self.stage!='allyl':raise ValueError('iodine cyclization requires the allylated substrate')
        return replace(self,stage='cyclized')
    @property
    def species(self):
        if self.stage=='native':return self.substrate
        label='a6m6A' if self.substrate=='m6A' else 'a6A'
        return ('cyclized '+label) if self.stage=='cyclized' else label
    def as_canonical_rna_base(self):
        if self.stage!='native' or self.substrate!='A':
            raise ValueError('a covalent RNA adduct/modification is not a canonical sequencing base')
        return 'A'

@dataclass(frozen=True)
class RNAPopulation:
    control: str
    selected_reacted_paths: tuple[AdenosinePath,...]
    residual_native_material: str='unassigned fractions retained; not complete conversion'
    soluble_material: str='present until physical washes; exact composition unassigned'


def published_adapter_inventory_scene():
    sc=Scene()
    sc.strand('stock adapter',[
        Segment('random adapter identifier','NNNNN',placeholder=True,inferred=True,
          feature=feature('sac_random_identifier','umi','random',note='UMI role inferred; actual read-processing assignment remains unreviewed')),
        Segment('printed index','ATCACG',feature=feature('sac_printed_index','sample_index','fixed')),
        Segment('adapter DNA','AGATCGGAAGAGCACACGTC'),
        Segment('internal biotin dT','T',note='Source /iBiodT/ retained'),
        Segment('terminal C3','[C3 linker]',placeholder=True)],label='publisher-SI stock adapter; not assigned to manuscript branch',mod5='rApp',mod3='C3 block')
    sc.note('stock adapter','Publisher-SI inventory only: 5′ activation, internal biotin-dT and terminal C3 remain explicit. Manuscript-version applicability unverified.')
    return sc


def graph_for_control(control):
    if control not in {'input','FTO-','FTO+'}:raise ValueError('source controls are input,FTO-,FTO+')
    g=ProtocolGraph('m6A-SAC-seq — '+control,authority=SOURCE,representation='published_instance')
    population=RNAPopulation(control,(AdenosinePath('A'),AdenosinePath('m6A')))
    prior=g.input('bound_rna','prepared RNA3′adapter–biotin/C1-bead population',MaterialKind.MOLECULAR,primary=True,
                  chemical_state=population)
    scope=['bound_rna']
    if control=='FTO+':
        fto=g.input('fto','FTO and source demethylation buffer',MaterialKind.CATALYST)
        population=replace(population,selected_reacted_paths=(*population.selected_reacted_paths,AdenosinePath('m6A').erase_selected()))
        treated=g.product('fto_treated','FTO-treated RNA; residual m6A and native A fractions unassigned',MaterialKind.MOLECULAR,chemical_state=population)
        g.react('demethylate','Source FTO background control; no complete erasure assumed',StepKind.CHEMICAL_CONVERSION,
                ((prior,InputRole.TEMPLATE),(fto,InputRole.CATALYST)),(treated,));prior=treated;scope+=['fto','demethylate']
        washed_fto=g.product('fto_washed','RNA after source FTO bead washes; residual partition unquantified',MaterialKind.MOLECULAR,chemical_state=population)
        fto_waste=g.product('fto_waste','FTO wash fractions',MaterialKind.WASTE)
        g.react('wash_fto','Physical source washing after FTO',StepKind.CLEANUP,((prior,InputRole.TEMPLATE),),(washed_fto,fto_waste));g.discard(fto_waste,'Actual FTO wash separation');prior=washed_fto;scope+=['wash_fto']
    if control!='input':
        enzyme=g.input('mjdim1','MjDim1',MaterialKind.CATALYST)
        sam=g.input('allyl_sam','allylic-SAM',MaterialKind.MOLECULAR)
        population=replace(population,selected_reacted_paths=tuple(p.allylate() for p in population.selected_reacted_paths))
        labeled=g.product('allylated','selected a6m6A and background a6A pathways; residual native RNA unassigned',MaterialKind.MOLECULAR,chemical_state=population)
        g.react('allylate','On-bead RNA labeling:50°C1h then six fresh20min rounds (Methods219)',StepKind.CHEMICAL_CONVERSION,
                ((prior,InputRole.TEMPLATE),(enzyme,InputRole.CATALYST),(sam,InputRole.STRUCTURAL)),(labeled,));prior=labeled;scope+=['mjdim1','allyl_sam','allylate']
        labeling_wash=g.product('labeling_wash','source wash fractions after labeling',MaterialKind.WASTE)
        washed_label=g.product('washed_label','bead-bound allylated RNA after source wash before iodine',MaterialKind.MOLECULAR,chemical_state=population)
        g.react('wash_label','Wash labeling reagents before iodine treatment',StepKind.CLEANUP,((prior,InputRole.TEMPLATE),),(washed_label,labeling_wash));g.discard(labeling_wash,'Actual labeling-stage wash');prior=washed_label;scope+=['wash_label']
        iodine=g.input('iodine','I2 treatment, room temperature/dark1h',MaterialKind.CONDITION)
        population=replace(population,selected_reacted_paths=tuple(p.cyclize() for p in population.selected_reacted_paths))
        cyclized=g.product('cyclized','covalent cyclized RNA adduct pathways plus reaction material; no RNA-base rewriting',MaterialKind.MOLECULAR,chemical_state=population)
        g.react('cyclize','Iodine-induced cyclization, then source quench; reaction products retained until wash',StepKind.CHEMICAL_CONVERSION,
                ((prior,InputRole.TEMPLATE),(iodine,InputRole.CONDITION)),(cyclized,));prior=cyclized;scope+=['iodine','cyclize']
        clean=g.product('washed_rna','washed bead-bound modified RNA; soluble partition unquantified',MaterialKind.MOLECULAR,
                        chemical_state=replace(population,soluble_material='partition after physical washing unquantified'))
        wash=g.product('wash_reject','source bead wash fractions',MaterialKind.WASTE)
        g.react('wash','Actual magnetic-bead washing after chemical treatment',StepKind.CLEANUP,((prior,InputRole.TEMPLATE),),(clean,wash));g.discard(wash,'Physical wash separation, not automatic chemical waste');prior=clean;scope+=['wash']
    g.mark_terminal(Terminal(prior.id,Disposition.UNRESOLVED,'HIV-RT copying, RNaseH/elution, second adapter ligation, exact PCR/index/read geometry remain unmodelled. Source reports HiSeqXTen PE150; final VOR procedure not reconciled.'))
    g.declare_reviewed_scope(scope,SOURCE)
    return g

def protocol_graph():return graph_for_control('FTO-')

def selftest():
    from checks import Check
    c=Check()
    for origin,species in [('A','cyclized a6A'),('m6A','cyclized a6m6A')]:
        result=AdenosinePath(origin).allylate().cyclize()
        c(origin+' keeps its chemical origin through labeling/cyclization',(result.original,result.species),(origin,species))
    try:AdenosinePath('m6A').cyclize()
    except ValueError:blocked=True
    else:blocked=False
    c('no cyclization before allyl installation',blocked,True)
    try:AdenosinePath('m6A').allylate().cyclize().as_canonical_rna_base()
    except ValueError:blocked=True
    else:blocked=False
    c('mutational readout cannot replace covalent RNA adduct with a base',blocked,True)
    erased=AdenosinePath('m6A').erase_selected().allylate().cyclize()
    c('FTO-erased selected m6A retains origin but follows A-background chemistry',(erased.original,erased.species),('m6A','cyclized a6A'))
    controlled=graph_for_control('FTO+').materials['fto_treated'].chemical_state.selected_reacted_paths
    c('FTO control keeps erased and residual m6A-origin paths without assigning fractions',sorted(x.fto_erased for x in controlled if x.original=='m6A'),[False,True])
    c('untreated input does not acquire experimental chemistry','allylate' in graph_for_control('input').steps,False)
    from checks import Source, have
    from storage import sources_root
    original=Source(sources_root()/'m6a-sac-seq__10.1038+s41587-022-01243-z'/'s41587-022-01243-z_41587_2022_1243_MOESM1_ESM.pdf.txt',
                    'publisher supplementary Table1, doi:10.1038/s41587-022-01243-z')
    if have(c,original,label='publisher-SI oligo transcription only'):
        text=''.join(original.path.read_text().split())
        c('publisher inventory retains rApp/biotin-dT/C3 tokens',RNA_ADAPTER_WRITTEN in text,True)
        c('publisher inventory RT matches source-assigned canonical suffix',RT_PRIMER in text,True)
    c.report()
if __name__=='__main__':selftest()
