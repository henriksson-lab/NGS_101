"""ac4C-seq selected RNA chemistry and complementary DNA report remain distinct."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Scene, Segment

def unreviewed_geometry(reason):
    sc=Scene()
    sc.strand("unassembled end",[Segment("source geometry not assembled","[unreviewed molecular geometry]",placeholder=True)],label="Local reading/assembly boundary")
    sc.note("unassembled end",reason+" This is not a claim that the original source fails to report an end.")
    return AdapterLigationGeometry("unavailable","unavailable",(sc,),reason)
SOURCE='Sas-Chen et al.2020 DOI10.1038/s41586-020-2418-2; accepted manuscript PMC8130014 Results100, library Methods228–230; VOR/SI applicability unverified'

@dataclass(frozen=True)
class Site:
    origin: str='ac4C'
    residue: str='ac4C'
    first_dna: str|None=None
    acetyl_products: bool=False
    def __post_init__(self):
        if self.origin not in {'ac4C','C'} or self.residue not in {'ac4C','C','N4-acetyltetrahydrocytidine'}:raise ValueError('selected cytidine-path classes only')
        if self.origin=='C' and self.residue!='C':raise ValueError('unmodified origin cannot acquire acetylation in this path')
        if self.first_dna is not None:
            allowed={'A','G'} if self.residue=='N4-acetyltetrahydrocytidine' else {'G'}
            if self.first_dna not in allowed:raise ValueError('outside selected C>T/nonmutation copying paths; background outcomes unmodelled')
    def deacetylate(self):
        if self.first_dna is not None:raise ValueError('RNA pretreatment precedes copying')
        return replace(self,residue='C',acetyl_products=self.origin=='ac4C')
    def reduce_selected(self):
        if self.first_dna is not None:raise ValueError('reduction precedes copying')
        return replace(self,residue='N4-acetyltetrahydrocytidine' if self.residue=='ac4C' else self.residue)
    def copy_selected(self,base='A'):return replace(self,first_dna=base)
    @property
    def reported(self):
        if self.first_dna is None:raise ValueError('reference report requires actual copying path')
        return {'A':'T','G':'C'}[self.first_dna]

def graph_for_control(control='reduced'):
    if control not in {'reduced','alkali_then_reduced','mock'}:raise ValueError('source chemical-control branches only')
    g=ProtocolGraph('ac4C-seq — '+control+' RNA/copying boundary',authority=SOURCE,representation='published_instance')
    site=Site();prior=g.input('rna','Source RNA containing caller-selected ac4C site; sequence/population fractions unassigned',K.MOLECULAR,primary=True,chemical_state=site)
    scope=['rna']
    if control=='alkali_then_reduced':
        alk=g.input('alkali','Source mild-alkali pretreatment; sequencing conditions in unreviewed SI',K.CONDITION,inferred=True)
        site=site.deacetylate();out=g.product('deacetylated','Selected RNA ac4C→C path plus acetyl-derived reaction products; no physical removal yet',K.MOLECULAR,chemical_state=site)
        g.react('deacetylate','Source deacetylation control; complete conversion not predicted',S.CHEMICAL_CONVERSION,((prior,R.TEMPLATE),(alk,R.CONDITION)),(out,));prior=out;scope+=['deacetylate']
    condition=g.input('reduction_condition','AcidicNaCNBH3' if control!='mock' else 'Source mock chemical treatment; sequencing branch buffer/workup unreviewed',K.CONDITION)
    if control!='mock':site=site.reduce_selected()
    changed=g.product('treated','Selected RNA residue '+site.residue+'; free reagents/other products retained, rates unassigned',K.MOLECULAR,chemical_state=site)
    g.react('treat','Source reduction/mock chemistry; validation timing not substituted for NGS branch',S.CHEMICAL_CONVERSION,((prior,R.TEMPLATE),(condition,R.CONDITION)),(changed,));scope+=['reduction_condition','treat']
    ready=g.product('ready','Source-treated RNA entering library preparation; actual chemical workup unassembled',K.MOLECULAR,chemical_state=site)
    g.react('unreviewed_workup','Unreviewed source chemical cleanup/handoff',S.MAGIC,((changed,R.TEMPLATE),),(ready,),reason='Do not infer an uninspected quench/purification or disappear reaction products.')
    ap=g.input('fastap','FastAP thermosensitive alkaline phosphatase',K.CATALYST)
    repaired=g.product('repaired','FastAP-treated RNA; source end geometry unassembled',K.MOLECULAR,chemical_state=site)
    g.react('fastap_treat','Source library preparation starts with FastAP',S.DEPHOSPHORYLATION,((ready,R.TEMPLATE),(ap,R.CATALYST)),(repaired,));scope+=['fastap','fastap_treat']
    adapter=g.input('rna_adapter','Source3′RNA adapter; exact original polymer/modifications unreviewed',K.ADAPTER)
    ligase=g.input('t4_ligase','Source T4 ligase; exact kit/type/end conditions unreviewed',K.CATALYST)
    tagged=g.product('tagged','Source adapter-ligated RNA class; other RNA populations unassigned',K.MOLECULAR,chemical_state=site)
    g.react('rna_ligation','Source3′RNA adapter ligation',S.LIGATION,((repaired,R.TEMPLATE),(adapter,R.STRUCTURAL),(ligase,R.CATALYST)),(tagged,),geometry=unreviewed_geometry('Source adapter end/polymer and complete ligation-ready RNA geometry not reviewed.'));scope+=['rna_adapter','t4_ligase','rna_ligation']
    primer=g.input('rt_primer','RT-priming reagent of cited strand-specific protocol; original sequence/version unreviewed',K.OLIGO)
    rt=g.input('tgirt','TGIRT-III',K.CATALYST)
    copied=site.copy_selected('A' if site.residue=='N4-acetyltetrahydrocytidine' else 'G')
    hybrid=g.product('hybrid','Selected complementary DNA path with original RNA retained; C>T is a reference-oriented report, not RNA thymidine',K.MOLECULAR,chemical_state={'site':copied,'RNA_retained':True},inferred=True)
    g.react('rt','Source TGIRT copying; selected mismatch/nonmutation class only',S.REVERSE_TRANSCRIPTION,((tagged,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST)),(hybrid,));scope+=['tgirt','rt']
    g.mark_terminal(Terminal(hybrid.id,Disposition.UNRESOLVED,'RNA/cDNA separation, cDNA3′adapter ligation/PCR/read geometry and exact chemical conditions remain unassembled; final VOR/SI reconciliation incomplete.'))
    g.declare_reviewed_scope(scope,SOURCE)
    g.add_source_item(SourceItem('unreviewed_chemistry_library','Publisher chemistry SI, exact treatment/adapter/RT-primer versions, predecessors39/40 and read plan','Published DOI10.1038/s41586-020-2418-2 and cited strand-specific library protocols'))
    return g

def protocol_graph():return graph_for_control()
def selftest():
    from checks import Check
    c=Check();p=Site().reduce_selected().copy_selected('A')
    c('reduced RNA, complementary DNA and reported base remain distinct',(p.residue,p.first_dna,p.reported),('N4-acetyltetrahydrocytidine','A','T'))
    q=Site().deacetylate().reduce_selected().copy_selected('G')
    c('deacetylated origin survives control and preserves released-material bookkeeping',(q.origin,q.residue,q.reported,q.acetyl_products),('ac4C','C','C',True))
    c('reduced template permits selected nonmutation copy',Site().reduce_selected().copy_selected('G').reported,'C')
    c('RNA retained after selected RT',protocol_graph().materials['hybrid'].chemical_state['RNA_retained'],True)
    c.report()
if __name__=='__main__':selftest()
