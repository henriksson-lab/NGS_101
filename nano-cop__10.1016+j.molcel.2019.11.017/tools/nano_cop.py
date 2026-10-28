"""nano-COP native-RNA tail and supportive-cDNA states, qualified to inspected manuscript."""
from dataclasses import dataclass, replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph, MaterialKind, InputRole, StepKind, Terminal, Disposition
from adapter_ligation import AdapterLigationGeometry
SOURCE='Drexler et al.2020 doi:10.1016/j.molcel.2019.11.017, accepted manuscript JATS54–56,160,176; VOR applicability unverified'

@dataclass(frozen=True)
class NativeRNA:
    artificial_tail: str='none'
    supporting_cdna: bool=False
    splint: str='none'
    def __post_init__(self):
        if self.artificial_tail not in {'none','A','I'}: raise ValueError('native tail identity must remain A or inosine, never generic DNA')
        if self.splint not in {'none','T10','C10'}: raise ValueError('unknown source splint')
        if self.splint=='C10' and self.artificial_tail!='I': raise ValueError('C10 branch requires the explicitly selected poly(I) RNA')
    @property
    def sequenced_polymer(self): return 'RNA'
    def tailed(self,kind):
        if self.artificial_tail!='none' or kind not in {'A','I'}: raise ValueError('select one explicit source tail route')
        return replace(self,artificial_tail=kind)
    def ligated(self):
        if self.artificial_tail=='none': raise ValueError('unassigned endogenous tail in no-tailing control is not a universally ligatable substrate')
        return replace(self,splint='C10' if self.artificial_tail=='I' else 'T10')
    def reverse_transcribed(self):
        if self.splint=='none': raise ValueError('source RT follows adapter ligation')
        return replace(self,supporting_cdna=True)


def graph_for_tail(tail):
    if tail not in {'A','I','none'}: raise ValueError('source branches are A,I,none')
    g=ProtocolGraph('nano-COP — '+tail+' tail branch',authority=SOURCE,representation='published_instance')
    rna=g.input('labelled_chr_rna','chromatin-associated RNA after8min4sU pulse and fractionation',MaterialKind.MOLECULAR,primary=True)
    reagent=g.input('biotin_hpdp','EZ-Link Biotin-HPDP reaction24°C1.5h',MaterialKind.CONDITION)
    tagged=g.product('biotin_mixture','biotinylated RNA plus soluble reaction products',MaterialKind.MOLECULAR)
    g.react('biotinylate','4sU-dependent biotinylation; retain reaction mixture before separation',StepKind.CHEMICAL_CONVERSION,
            ((rna,InputRole.TEMPLATE),(reagent,InputRole.CONDITION)),(tagged,))
    clean=g.product('clean_tagged_rna','RNA after chloroform/isopropanol recovery',MaterialKind.MOLECULAR)
    waste=g.product('tagging_waste','extraction phases and precipitation supernatant',MaterialKind.WASTE)
    g.react('recover_tagged','Actual extraction and isopropanol precipitation',StepKind.CLEANUP,((tagged,InputRole.TEMPLATE),),(clean,waste));g.discard(waste,'Physical separation, not reaction-generated automatic waste')
    beads=g.input('streptavidin','μMACS streptavidin beads',MaterialKind.SELECTION)
    bound=g.product('bound_rna','biotinylated4sU RNA on streptavidin',MaterialKind.MOLECULAR)
    unbound=g.product('unbound','affinity flow-through and washes',MaterialKind.WASTE)
    g.react('capture','Streptavidin affinity capture and source washes',StepKind.AFFINITY_CAPTURE,
            ((clean,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(bound,unbound));g.discard(unbound,'Physical affinity reject')
    dtt=g.input('dtt','0.1M DTT elution',MaterialKind.CONDITION)
    eluate=g.product('released_rna','DTT-eluted RNA-containing mixture; no RNA-base rewrite',MaterialKind.MOLECULAR)
    g.react('release','DTT release of selected RNA',StepKind.CHEMICAL_CONVERSION,
            ((bound,InputRole.TEMPLATE),(dtt,InputRole.CONDITION)),(eluate,))
    purified=g.product('purified_rna','miRNeasy-purified RNA with on-column DNase',MaterialKind.MOLECULAR,chemical_state=NativeRNA())
    reject=g.product('purification_waste','column and cleanup reject',MaterialKind.WASTE)
    g.react('purify','miRNeasy cleanup after elution',StepKind.CLEANUP,((eluate,InputRole.TEMPLATE),),(purified,reject));g.discard(reject,'Physical column reject')
    kit=g.input('rrna_kit','source Ribo-Zero Gold or RiboMinus v2; per-sample assignment in gated TableS1',MaterialKind.SELECTION)
    depleted=g.product('depleted','rRNA-depleted RNA population; kit molecular mechanics unreviewed',MaterialKind.MOLECULAR,chemical_state=NativeRNA())
    removed=g.product('rrna_material','depletion-rejected molecular material; physical final fate unreviewed',MaterialKind.MOLECULAR)
    g.react('rrna_deplete','Source kit rRNA depletion; molecular implementation unresolved',StepKind.MAGIC,
            ((purified,InputRole.TEMPLATE),(kit,InputRole.SELECTION)),(depleted,removed),reason='The original names two rRNA depletion kits; molecular steps and rejected-material fate remain unreviewed.')
    g.mark_terminal(Terminal(removed.id,Disposition.UNRESOLVED,'Source rRNA kit internals and rejected-material fate unreviewed.'))
    endpoint=depleted
    if tail!='none':
        enzyme=g.input('pap','yeast poly(A) polymerase37°C30min' if tail=='I' else 'E.coli PAP/NEB37°C1h (unperturbed branch)',MaterialKind.CATALYST)
        ntp=g.input('ntp','ITP' if tail=='I' else 'ATP required for poly(A) synthesis; exact buffer not curated',MaterialKind.MOLECULAR,inferred=tail=='A')
        tailed_state=NativeRNA().tailed(tail)
        tailed=g.product('tailed','native RNA with artificial poly('+tail+') tail; tail length unassigned',MaterialKind.MOLECULAR,chemical_state=tailed_state)
        g.react('tail','Source enzymatic RNA tailing',StepKind.TAILING,
                ((depleted,InputRole.TEMPLATE),(enzyme,InputRole.CATALYST),(ntp,InputRole.STRUCTURAL)),(tailed,))
        adapter=g.input('adapter','custom C10 splinted adapter (TableS4)' if tail=='I' else 'SQK-RNA001/002 T10 splinted adapter',MaterialKind.ADAPTER)
        ligase=g.input('ligase','T4 DNA ligase',MaterialKind.CATALYST)
        joined_state=tailed_state.ligated()
        joined=g.product('adapter_rna','RNA with source splinted3′adapter; full adapter sequence unreviewed',MaterialKind.MOLECULAR,chemical_state=joined_state)
        g.react('ligate','Source3′ adapter ligation15min at room temperature',StepKind.LIGATION,
                ((tailed,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(ligase,InputRole.CATALYST)),(joined,),
                geometry=AdapterLigationGeometry.unavailable('Reported C10/poly(I) or T10/poly(A) recognition is known; full source oligos/end geometry not reviewed (custom TableS4 gated).'))
        rt=g.input('rt','SuperScript III',MaterialKind.CATALYST)
        rt_primer=g.input('rt_primer','source RT primer role; exact kit oligo identity unreviewed',MaterialKind.OLIGO,inferred=True)
        endpoint=g.product('supported_rna','native RNA with supportive cDNA; RNA is the sequenced polymer',MaterialKind.MOLECULAR,chemical_state=joined_state.reverse_transcribed())
        g.react('supportive_rt','Source supportive RT; resulting cDNA is not sequenced',StepKind.REVERSE_TRANSCRIPTION,
                ((joined,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(endpoint,))
    g.mark_terminal(Terminal(endpoint.id,Disposition.UNRESOLVED,'Source SQK-RNA001/002 motor-adapter attachment, complete chemical state and acquisition geometry remain unreviewed; no current ONT kit substituted.'))
    g.declare_reviewed_scope(('labelled_chr_rna','biotin_hpdp','biotinylate','recover_tagged','streptavidin','capture','dtt','release','purify','rrna_kit','rrna_deplete'),SOURCE)
    return g

def protocol_graph(): return graph_for_tail('I')

def selftest():
    from checks import Check
    c=Check();s=NativeRNA().tailed('I').ligated().reverse_transcribed()
    c('inosine tail retains C10 recognition and native RNA readout after supportive RT',
      (s.artificial_tail,s.splint,s.supporting_cdna,s.sequenced_polymer),('I','C10',True,'RNA'))
    try: NativeRNA().ligated()
    except ValueError: rejected=True
    else: rejected=False
    c('no-tail control does not invent universal adapter-compatible native ends',rejected,True)
    c('poly(A) branch retains separate T10 state',graph_for_tail('A').materials['supported_rna'].chemical_state.splint,'T10')
    c.report()
if __name__=='__main__': selftest()
