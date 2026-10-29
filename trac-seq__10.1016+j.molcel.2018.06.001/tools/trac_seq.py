"""TRAC-seq selected m7G cleavage with retained site/upstream material lineage."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
SOURCE='Lin et al.2018 DOI10.1016/j.molcel.2018.06.001; accepted manuscript PMC6086580 Results101–103, Discussion137, Methods351–353; VOR/SI applicability unverified'

@dataclass(frozen=True)
class SelectedSite:
    rna: str
    position: int
    selected_m7g: bool
    reduced: bool=False
    cleaved: bool=False
    def __post_init__(self):
        if not self.rna or set(self.rna)-set('ACGU') or type(self.position) is not int or not 0<=self.position<len(self.rna)-1 or self.rna[self.position]!='G':raise ValueError('explicit selected internal G and downstream RNA required')
        if not self.selected_m7g and (self.reduced or self.cleaved):raise ValueError('this selected m7G-reaction relation cannot infer modification from sequence')
        if self.cleaved and not self.reduced:raise ValueError('source reduction precedes selected cleavage')
    def reduce(self):
        if not self.selected_m7g:raise ValueError('caller must explicitly designate selected m7G substrate')
        return replace(self,reduced=True)
    def cleave(self):
        if not self.reduced:raise ValueError('source reduction before aniline')
        return replace(self,cleaved=True)
    @property
    def downstream(self):
        if not self.cleaved:raise ValueError('no released downstream fragment before cleavage')
        return self.rna[self.position+1:]
    @property
    def non_downstream_lineage(self):
        if not self.cleaved:raise ValueError('no assigned product partition before cleavage')
        return self.rna[:self.position+1]
    @property
    def downstream5(self):return 'phosphate' if self.cleaved else 'not assigned'

def protocol_graph():
    g=ProtocolGraph('TRAC-seq — reduction/aniline RNA prefix',authority=SOURCE,representation='published_instance')
    prior=g.input('alkb_rna','Source WT/D135S-AlkB-treated smallRNA; upstream actual reactions/other-modification fractions unassembled',K.MOLECULAR,primary=True)
    hydride=g.input('hydride','0.1MNaBH4,30min on ice/in dark',K.CONDITION)
    carrier=g.input('free_m7gtp','1mM free m7GTP methylation carrier; separate from covalent RNA residues',K.MOLECULAR)
    reduced=g.product('reduced','Reduced RNA-site paths plus free carrier/reduction products and unreacted material; no residue deletion',K.MOLECULAR,chemical_state={'RNA_site':'reduction-treated selected m7G; detailed reduced/intermediate structure unassigned','free_carrier':'separate molecular population; not RNA incorporation'})
    g.react('reduce','Source RNA reduction in free m7GTP-containing solution',S.CHEMICAL_CONVERSION,((prior,R.TEMPLATE),(hydride,R.CONDITION),(carrier,R.STRUCTURAL)),(reduced,))
    def precipitate(key,x):
        reagent=g.input(key+'_reagents','Source sodium acetate/coldethanol precipitation',K.SELECTION)
        keep=g.product(key+'_rna','Physically recovered RNA-containing population; soluble-product carryover unassigned',K.MOLECULAR)
        waste=g.product(key+'_waste','Physical precipitation supernatant/wash fractions',K.WASTE)
        g.react(key,'Source ethanol precipitation',S.CLEANUP,((x,R.TEMPLATE),(reagent,R.SELECTION)),(keep,waste));g.discard(waste,'Physical source precipitation fractions');return keep
    ready=precipitate('first_precipitation',reduced)
    aniline=g.input('aniline','Water/glacialaceticacid/aniline7:3:1,roomtemperature/dark2h',K.CONDITION)
    cleaved=g.product('cleaved','Selected downstreamRNA5′p after m7G, plus upstream/site-derived and uncut material retained in reaction',K.MOLECULAR,chemical_state={'downstreamRNA':'starts after selected m7G; source5′phosphate','other_material':'upstream and chemically changed site residue lineage retained; detailed structures/ends not assigned'})
    g.react('cleave','Source aniline-induced cleavage; not computational deletion of methylated G',S.CHEMICAL_CONVERSION,((ready,R.TEMPLATE),(aniline,R.CONDITION)),(cleaved,))
    final=precipitate('second_precipitation',cleaved)
    g.mark_terminal(Terminal(final.id,Disposition.UNRESOLVED,'Source cDNA library handoff; actual AlkB front-end, historical smallRNA adapter/end/PCR/read architecture and final VOR/SI applicability unassembled.'))
    g.declare_reviewed_scope(['alkb_rna','hydride','free_m7gtp','reduce','first_precipitation_reagents','first_precipitation','aniline','cleave','second_precipitation_reagents','second_precipitation'],SOURCE)
    g.add_source_item(SourceItem('unreviewed_full_chemistry','Defining supplementary originals, older cleavage chemistry, AlkB/historical kit and final read layout','Defining DOI10.1016/j.molcel.2018.06.001 and cited source methods'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=SelectedSite('ACGUA',2,True).reduce().cleave()
    c('released3′fragment starts after selectedG with source5′p',(p.downstream,p.downstream5),('UA','phosphate'))
    c('selected site is retained in non-downstream residue lineage',p.non_downstream_lineage+p.downstream,p.rna)
    c('non-RAGGU context is not falsely excluded when caller supplies selectedm7G',p.selected_m7g,True)
    try:SelectedSite('ACGUA',2,False).reduce();bad=False
    except ValueError:bad=True
    c('sequenceG alone does not imply reactive methylation',bad,True)
    original=Source(sources_root()/'trac-seq__10.1016+j.molcel.2018.06.001'/'j.molcel.2018.06.001_PMC6086580.xml.txt','accepted manuscript Methods353')
    if have(c,original,label='source freecarrier/cleavage check'):
        text=' '.join(original.path.read_text().split());c('source explicitly identifies free m7GTP carrier','1 mM free m7GTP' in text,True)
    c.report()
if __name__=='__main__':selftest()
