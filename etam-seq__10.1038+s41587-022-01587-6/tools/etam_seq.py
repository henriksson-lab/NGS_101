"""eTAM-seq selected chemical/copying paths with independent RNA and reported-base identity."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct, Scene, Segment, feature
import illumina as il
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition
SOURCE='Xiao et al.2023 doi:10.1038/s41587-022-01587-6, accepted manuscript Methods203–211 and publisher SI Note1; final VOR procedure unverified'
RNA_ADAPTER='/5rApp/'+il.INDEX2_PRIMER_RC[:19]+'/3Bio/'
RT_PRIMER=il.TRUSEQ_READ1[-20:]
CDNA_ADAPTER='/5Phos/'+'N'*6+il.TRIM_SEEN_IN_READ1[:22]+'/3SpC3/'

@dataclass(frozen=True)
class SitePath:
    original: str
    rna_residue: str
    first_strand_dna: str | None=None
    fto_erased: bool=False
    def __post_init__(self):
        allowed={('m6A','A'),('m6A','I')} if self.fto_erased else {('A','A'),('A','I'),('m6A','m6A')}
        if (self.original,self.rna_residue) not in allowed:
            raise ValueError('only declared A-to-I or protected m6A paths supported')
        if self.first_strand_dna is not None and self.first_strand_dna!={'A':'T','I':'C','m6A':'T'}[self.rna_residue]:
            raise ValueError('this selected complementary-copy path cannot assert a different DNA base')
    def erase_selected(self):
        if self.rna_residue!='m6A' or self.first_strand_dna is not None:raise ValueError('FTO acts on the selected native RNA m6A site')
        return replace(self,rna_residue='A',fto_erased=True)
    def deaminate_selected_site(self):
        if self.first_strand_dna is not None:raise ValueError('deamination operates on RNA before copying')
        return replace(self,rna_residue='I' if self.rna_residue in {'A','I'} else 'm6A')
    def copy(self):
        return replace(self,first_strand_dna={'A':'T','I':'C','m6A':'T'}[self.rna_residue])
    def reported_reference_base(self):
        if self.first_strand_dna is None:raise ValueError('reported DNA base requires copying first')
        return {'A':'T','C':'G','G':'C','T':'A'}[self.first_strand_dna]

@dataclass(frozen=True)
class SampleState:
    control: str
    selected_site_paths: tuple[SitePath,...]
    residual_conversion_and_fto: str='unassigned fractions; residual A is not by itself proof of m6A'
    rna_retained: bool=True
    stage: str='RNA adapter ligated / bead bound'


def cdna_adapter_inventory():
    sc=Scene();sc.strand('cDNA adapter',[
        Segment('random identifier','NNNNNN',placeholder=True,inferred=True,
          feature=feature('etam_random_identifier','umi','random',note='identifier role inferred pending complete read-processing review')),
        Segment('adapter DNA',il.TRIM_SEEN_IN_READ1[:22]),
        Segment('C3 linker','[C3 linker]',placeholder=True)],mod5='phosphate',mod3='C3 block')
    return sc


def graph_for_control(control):
    if control not in {'FTO-','FTO+','IVT'}:raise ValueError('source controls are FTO-,FTO+,IVT')
    g=ProtocolGraph('eTAM-seq — '+control,authority=SOURCE,representation='published_instance')
    state=SampleState(control,(SitePath('A','A'),) if control=='IVT' else (SitePath('A','A'),SitePath('m6A','m6A')))
    prior=g.input('bound_rna','RNA3′adapter-ligated/biotin-C1 population; '+control+' target (calibration probes separate)',MaterialKind.MOLECULAR,primary=True,chemical_state=state)
    scope=['bound_rna']
    if control=='FTO+':
        fto=g.input('fto','FTO in specified demethylation buffer37°C1h',MaterialKind.CATALYST)
        state=replace(state,selected_site_paths=(*state.selected_site_paths,SitePath('m6A','m6A').erase_selected()))
        treated=g.product('fto_treated','FTO-treated target RNA; residual modification fractions unassigned',MaterialKind.MOLECULAR,chemical_state=state)
        g.react('demethylate','Source FTO control; no complete demethylation asserted',StepKind.CHEMICAL_CONVERSION,
                ((prior,InputRole.TEMPLATE),(fto,InputRole.CATALYST)),(treated,));prior=treated;scope+=['fto','demethylate']
        washed_fto=g.product('fto_washed','RNA after source FTO bead washing',MaterialKind.MOLECULAR,chemical_state=state)
        fto_waste=g.product('fto_waste','FTO wash fractions',MaterialKind.WASTE)
        g.react('wash_fto','Source FTO wash before TadA treatment',StepKind.CLEANUP,((prior,InputRole.TEMPLATE),),(washed_fto,fto_waste));g.discard(fto_waste,'Physical FTO wash reject');prior=washed_fto;scope+=['wash_fto']
    tada=g.input('tada','TadA8.20',MaterialKind.CATALYST)
    state=replace(state,selected_site_paths=tuple(x.deaminate_selected_site() for x in state.selected_site_paths),stage='selected reacted RNA paths; native/unconverted material retained as unassigned fractions')
    changed=g.product('deaminated','selected RNA A→I path and resistant m6A path; no DNA-base substitution in RNA',MaterialKind.MOLECULAR,chemical_state=state)
    g.react('deaminate','On beads:53°C1h, then two fresh44°C1h rounds',StepKind.CHEMICAL_CONVERSION,
            ((prior,InputRole.TEMPLATE),(tada,InputRole.CATALYST)),(changed,));scope+=['tada','deaminate']
    washed=g.product('washed','washed bead-bound RNA after deamination',MaterialKind.MOLECULAR,chemical_state=state)
    reject=g.product('wash_waste','deamination supernatants and washes',MaterialKind.WASTE)
    g.react('wash','Actual source bead washing; soluble material not automatically waste before separation',StepKind.CLEANUP,
            ((changed,InputRole.TEMPLATE),),(washed,reject));g.discard(reject,'Physical wash fractions');scope+=['wash']
    primer=g.input('rt_primer','source RT primer',MaterialKind.OLIGO,state=Construct([Segment('RT primer',RT_PRIMER)]))
    rt=g.input('rt','Maxima H− RT50°C1h',MaterialKind.CATALYST)
    copied=replace(state,selected_site_paths=tuple(x.copy() for x in state.selected_site_paths),stage='RNA–cDNA hybrid; source variable RNA sequence and mutation/background fractions unassigned')
    cdna=g.product('hybrid','conditional complementary first-strand path: RNA I pairs DNA C; reference-oriented report G',MaterialKind.MOLECULAR,chemical_state=copied,inferred=True)
    g.react('reverse_transcribe','Source Maxima RT; local selected-site complement relation only',StepKind.REVERSE_TRANSCRIPTION,
            ((washed,InputRole.TEMPLATE),(primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,));scope+=['rt_primer','rt','reverse_transcribe']
    g.mark_terminal(Terminal(cdna.id,Disposition.UNRESOLVED,'Source cDNA release/cleanup, C3-blocked adapter ligation, E6440S unique-dual-index PCR and acquired read plan remain unmodelled. Full RNA sequence/background conversion and VOR applicability unresolved.'))
    g.declare_reviewed_scope(scope,SOURCE)
    return g

def protocol_graph():return graph_for_control('FTO-')

def selftest():
    from checks import Check
    c=Check();p=SitePath('A','A').deaminate_selected_site().copy()
    c('RNA chemistry, complementary first strand and reference-oriented report remain distinct',
      (p.rna_residue,p.first_strand_dna,p.reported_reference_base()),('I','C','G'))
    q=SitePath('m6A','m6A').deaminate_selected_site().copy()
    c('protected m6A path remains chemically methylated after copying',(q.rna_residue,q.reported_reference_base()),('m6A','A'))
    erased=SitePath('m6A','m6A').erase_selected().deaminate_selected_site().copy()
    c('FTO-erased selected m6A is deaminatable while retaining source identity',(erased.original,erased.rna_residue,erased.reported_reference_base()),('m6A','I','G'))
    controlled=graph_for_control('FTO+').materials['deaminated'].chemical_state.selected_site_paths
    c('FTO control distinguishes resistant and erased/deaminated m6A-origin paths',sorted((x.rna_residue,x.fto_erased) for x in controlled if x.original=='m6A'),[('I',True),('m6A',False)])
    c('IVT target control does not invent endogenous methylation',
      [x.original for x in graph_for_control('IVT').materials['bound_rna'].chemical_state.selected_site_paths],['A'])
    c('RNA template is retained in nascent cDNA hybrid',protocol_graph().materials['hybrid'].chemical_state.rna_retained,True)
    from checks import Source, have
    from storage import sources_root
    original=Source(sources_root()/'etam-seq__10.1038+s41587-022-01587-6'/'s41587-022-01587-6_PMC10625715.xml.txt',
                    'accepted-manuscript Methods205/211, doi:10.1038/s41587-022-01587-6')
    if have(c,original,label='accepted-manuscript oligo transcription'):
        text=''.join(original.path.read_text().split())
        for label,sequence in [('RNA adapter',RNA_ADAPTER),('RT primer',RT_PRIMER),('cDNA adapter',CDNA_ADAPTER)]:
            c('source '+label+' sequence/modifications match',sequence in text,True)
    c.report()
if __name__=='__main__':selftest()
