"""CeU-seq source-limited selected probe chemistry; not an assembled library route."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition
SOURCE='Li et al.2015 DOI10.1038/nchembio.1836; publisher SI FigS2(b), S3 caption, Supplementary NoteS21–S22; main procedure gated'

@dataclass(frozen=True)
class ProbeSite:
    residue: str
    n3_cmc: bool=False
    biotin_clicked: bool=False
    rna_backbone_retained: bool=True
    def __post_init__(self):
        if self.residue not in {'Psi','U','I'}:raise ValueError('only source probe classes represented')
        if self.n3_cmc and self.residue!='Psi':raise ValueError('selected FigS2 probe path has no U/I cross-reactive adduct')
        if self.biotin_clicked and not self.n3_cmc:raise ValueError('click requires azido-CMC-bearing probe in this selected path')
        if not self.rna_backbone_retained:raise ValueError('label/click does not delete RNA residue')
    def label_selected(self):return replace(self,n3_cmc=self.residue=='Psi')
    def click_selected(self):
        if not self.n3_cmc:raise ValueError('selected RNA-associated click requires N3-CMC adduct')
        return replace(self,biotin_clicked=True)

def protocol_graph():
    g=ProtocolGraph('CeU-seq — published probe-labeling event only',authority=SOURCE,representation='published_instance')
    prior=g.input('rna_probe','Published model RNA probe classes Ψ/U/I; exact sequence/procedure not reviewed',K.MOLECULAR,primary=True,chemical_state=tuple(ProbeSite(x) for x in ('Psi','U','I')))
    condition=g.input('unreviewed_conditions','Actual FigS2 experimental conditions are unreviewed; no buffer, time or temperature assigned',K.CONDITION,inferred=True)
    reagent=g.input('n3_cmc','N3-CMC azido carbodiimide; synthesized in source Supplementary Note',K.MOLECULAR)
    labelled=g.product('labelled','Selected FigS2 Ψ–N3-CMC adduct; U/I probe classes unmodified in this illustrated relation',K.MOLECULAR,chemical_state=tuple(ProbeSite(x).label_selected() for x in ('Psi','U','I')))
    g.react('label','Source FigS2 probe-labeling relation; no unreviewed buffer/rate/yield assignment',S.CHEMICAL_CONVERSION,((prior,R.TEMPLATE),(reagent,R.STRUCTURAL),(condition,R.CONDITION)),(labelled,))
    biotin=g.input('dbco_biotin','DBCO–biotin (printed FigS2b)',K.MOLECULAR)
    clicked=g.product('clicked','Selected Ψ–CMC–biotin adduct; RNA backbone and residue identity retained; unreacted fractions unassigned',K.MOLECULAR,chemical_state=(ProbeSite('Psi').label_selected().click_selected(),ProbeSite('U'),ProbeSite('I')))
    g.react('click','Source FigS2 DBCO–biotin conjugation; no copper catalyst inferred',S.CHEMICAL_CONVERSION,((labelled,R.TEMPLATE),(biotin,R.STRUCTURAL),(condition,R.CONDITION)),(clicked,))
    g.mark_terminal(Terminal(clicked.id,Disposition.UNRESOLVED,'Main article chemistry/workup, reaction order/conditions, RNA enrichment, actual RT termination/cDNA adaptation and sequencing/read layout require defining full text. FigS3 supports pulldown/stop signal but is not sufficient to reconstruct that library.'))
    g.declare_reviewed_scope(['rna_probe','n3_cmc','label','dbco_biotin','click'],SOURCE);return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=ProbeSite('Psi').label_selected().click_selected()
    c('click appends biotin while preserving labelled RNA residue',(p.residue,p.n3_cmc,p.biotin_clicked,p.rna_backbone_retained),('Psi',True,True,True))
    c('selected source U/I controls preserve original identity',tuple(ProbeSite(x).label_selected().n3_cmc for x in ('U','I')),(False,False))
    try:ProbeSite('Psi').click_selected();bad=False
    except ValueError:bad=True
    c('biotin is not attached before selected azido handle',bad,True)
    original=Source(sources_root()/'ceu-seq__10.1038+nchembio.1836'/'nchembio.1836_41589_2015_BFnchembio1836_MOESM157_ESM.pdf.txt','publisher FigS2 caption; image separately inspected')
    if have(c,original,label='source figure caption'):
        text=' '.join(original.path.read_text().split())
        c('source distinguishes click conjugation from RNA deletion','Biotin can be conjugated by click chemistry.' in text,True)
    c.report()
if __name__=='__main__':selftest()
