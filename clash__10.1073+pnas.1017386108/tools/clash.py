"""CLASH original: source-limited selected RNA-arm joining, no inferred library recipe."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
from adapter_ligation import AdapterLigationGeometry
SOURCE='Kudla et al.2011 DOI10.1073/pnas.1017386108; published PMC3116431 Results1129–1151 and Methods1218–1220; SI and CRAC predecessor unreviewed'

@dataclass(frozen=True)
class SelectedJoin:
    first: str
    second: str
    donor5: str
    acceptor3: str
    joined: bool=False
    def __post_init__(self):
        if not self.first or not self.second or set(self.first+self.second)-set('ACGU'):raise ValueError('two explicit RNA arms required')
        if self.donor5 not in {'phosphate','unassigned','consumed'} or self.acceptor3 not in {'OH','unassigned','consumed'}:raise ValueError('invalid end state')
        if self.joined and (self.donor5,self.acceptor3)!=('consumed','consumed'):raise ValueError('joined junction consumes selected end functionalities')
    def ligate(self):
        if (self.donor5,self.acceptor3)!=('phosphate','OH') or self.joined:raise ValueError('selected ligation requires explicitly compatible ends; source workup unreviewed')
        return replace(self,donor5='consumed',acceptor3='consumed',joined=True)
    @property
    def product(self):
        if not self.joined:raise ValueError('unjoined arms do not constitute a chimeric molecule')
        return self.first+self.second

def protocol_graph():
    g=ProtocolGraph('CLASH — nickel-immobilized hybrid ligation boundary',authority=SOURCE,representation='published_instance')
    rna=g.input('nickel_complex','Source UV-crosslinked His-tagged protein/RNA complex on nickel resin after IgG/TEV/6Mguanidine purification; actual RNA end workup unreviewed',K.MOLECULAR,primary=True,chemical_state={'RNA_arms':'possibly base-paired; sequences and5′/3′ends unassigned','covalent_RNP':'retained','residual_background':'unassigned'})
    condition=g.input('ligation_conditions','Reported linker-ligation reactions16°C or25°C; SI ligase/reagent recipe not inspected',K.CONDITION,inferred=True)
    ligase=g.input('unreviewed_ligase','Source ligation enzyme identity/conditions unreviewed in gated SI; no specific ligase inferred',K.CATALYST,inferred=True)
    hybrid=g.product('hybrid','Source-supported hybrid-RNA class: direct arm junction and two remaining ends available for linkers; unligated/off-target populations unresolved',K.MOLECULAR,chemical_state={'junction':'conditional direct RNAarm1–RNAarm2','adapter_sequences':'unreviewed','physical_other_products':'remain unassigned in reaction'},inferred=True)
    g.react('hybrid_ligation','Selected hybrid-forming ligation; do not substitute an RT template-switching event',S.PROXIMITY_LIGATION,((rna,R.TEMPLATE),(condition,R.CONDITION),(ligase,R.CATALYST)),(hybrid,),geometry=AdapterLigationGeometry.unavailable('Original SI/predecessor end preparation and ligase recipe not inspected; no phosphates or precise stem-end placement are assigned.'))
    g.mark_terminal(Terminal(hybrid.id,Disposition.UNRESOLVED,'Detailed CRAC purification/end repair/ligase and adapter chemistry, protein release/RT/PCR and acquired50nt Illumina–Solexa read geometry remain unassembled. Source supports the hybrid class, not this complete molecular reconstruction.'))
    g.declare_reviewed_scope(['nickel_complex','hybrid_ligation'],SOURCE)
    g.add_source_item(SourceItem('unreviewed_chemistry','SI Materials/Methods and cited CRAC predecessor, including ligase, RNA ends and linkers','PMC3116431 SI link; DOI10.1073/pnas.0901997106'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=SelectedJoin('ACGU','UGCA','phosphate','OH').ligate()
    c('selected direct ligation retains both complete RNA arms',p.product,'ACGUUGCA')
    c('selected junction consumes compatible end groups',(p.donor5,p.acceptor3),('consumed','consumed'))
    try:SelectedJoin('ACGU','UGCA','unassigned','unassigned').ligate();bad=False
    except ValueError:bad=True
    c('unknown source ends cannot silently ligate in constructor',bad,True)
    original=Source(sources_root()/'clash__10.1073+pnas.1017386108'/'pnas.1017386108_PMC3116431.html.txt','published Methods1218–1220')
    if have(c,original,label='published library method summary'):
        text=' '.join(original.path.read_text().split())
        c('source retains hybrid ligation on immobilized protein complexes','Protein–RNA complexes were immobilized on nickel beads' in text,True)
    c.report()
if __name__=='__main__':selftest()
