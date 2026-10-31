"""Private source-qualified partial preparation; no complete molecular/absence claim."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
class B:
    def __init__(self,title,authority,label,kind=K.MOLECULAR):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',label,kind,primary=True)
    def stock(self,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cat(self,key,label):return self.stock(key,label,K.CATALYST,R.CATALYST)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind=S.MAGIC,inputs=(),product='Actual source preparation; whole occupied substrates/products/side products and conversion unassembled',**kw):
        y=self.g.product(key+'_product',product,K.MOLECULAR)
        if kind==S.MAGIC:kw.setdefault('reason','Positive source operation; whole occupied molecules, end distributions, detailed molecular products and conversion unassembled, not a source-absence claim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Actual retained physical source fraction, carryover/recovery unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual excluded workup fraction',K.WASTE)
        self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Actual source workup fraction excluded; no quantitative purification claim.');self.x=y;return y
    def aliquot(self,key,label,kept='Selected source aliquot',rest='Untransferred stock, amount/fate unassigned',qc=False):
        y=self.g.product(key+'_selected',kept,K.QC if qc else K.MOLECULAR);r=self.g.product(key+'_remainder',rest,K.MOLECULAR)
        self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),(y,r))
        if qc:self.g.mark_terminal(Terminal(y.id,Disposition.QC,kept));self.x=r;return r
        self.g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,rest));self.x=y;return y


AUTHORITY='Published main Nature HTML abstract/ExtendedData4 caption; actual author/publisher corrections read. Full main Methods and all SI originals uninspected. Narrow in-vitro TetrahymenaRNA5min branch only, no full kit or successor equivalence.'
def protocol_graph():
    b=B('PORE-cupine selected chemical-probed native-RNA preparation',AUTHORITY,'Source in-vitro transcribed full-length Tetrahymena RNA; synthesis/whole sequence/ends unassembled');g=b.g
    g.add_source_item(SourceItem('author_correction','Published AuthorCorrection10.1038/s41587-020-00755-w actual body707: missing Lambert affiliation and acknowledgment.','Retain stated correction scope only; does not establish complete article/recipe equivalence.'))
    g.add_source_item(SourceItem('publisher_correction','Published PublisherCorrection10.1038/s41587-021-00889-5 actual body713: initially attached ReportingSummary belonged to another paper.','Use current main caption only for selected preparation; corrected summary/SI chemistry not read or inherited.'))
    g.add_source_item(SourceItem('selected_scope','Published main ExtendedData4f–i caption: NAI-N3 100mM5min and25min comparison with Bioanalyzer/mapping/readlength results.','Model selects source5min condition; no inference of complete adduct occupancy, RNA ends, library recipe or modification-to-signal specificity. Unmodified and25min comparisons are not pooled into this branch.'))
    b.op('probe','Source NAI-N3 100mM5min chemical RNA probing',inputs=(b.stock('nai_n3','Actual source NAI-N3 chemical contributor; molecular adduct populations unassembled'),b.cond('probe_conditions','Source100mM5min; other detailed reaction conditions unassembled')),product='Source chemically modified Tetrahymena RNA preparation; whole adduct locations/stereochemistry and unreacted/off-product populations unassembled')
    b.aliquot('bioanalyzer','Actual source Bioanalyzer analytical allocation; amount unassigned',qc=True)
    b.op('native_finish','Source direct native-RNA nanopore library preparation; detailed finishing order unassembled',inputs=(b.stock('library_reagents','Historical source direct-RNA library reagent inventory; complete oligos, protein, nucleotide supplies and packaging unassembled'),),product='Source native modified-RNA sequencing preparation; whole RNA/adapter/motor architecture unassembled, no cDNA sequencing product substituted')
    b.aliquot('loading','Source selected native-RNA sequencing allocation; amount unassigned')
    lib=g.product('library','Source selected modified native-RNA nanopore library; whole molecular geometry unassembled',K.LIBRARY)
    g.react('load','Actual selected native-RNA sequencing stock',S.GENERIC,((b.x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Published main abstract/ExtendedData4 native long-read RNA sequencing. Source bioinformatic length>50percent annotated transcript and100modified/200unmodified read thresholds are computational filters, not physical RNA size selection. Flowcell/motor/run/read landing unassembled; no successor or kit defaults.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Probe chemical contributes molecular inventory',g.steps['probe'].inputs[1].role,R.STRUCTURAL)
    c('QC sample is separate physical allocation',g.terminals['bioanalyzer_selected'].disposition,Disposition.QC)
    c('Source Bioanalyzer analytical fraction has QC material kind',g.materials['bioanalyzer_selected'].kind,K.QC)
    c('Analytical fraction does not enter native RNA finishing',g.steps['native_finish'].inputs[0].material,'bioanalyzer_remainder')
    c('Actual selected source and both correction gates remain open',all(any(repr(k) in p for p in g.publication_problems()) for k in ('author_correction','publisher_correction','selected_scope')),True)
    c('Loading stock remainder retained','loading_remainder' in g.terminals,True)
    c('No copied cDNA transition imported',all(x.kind!=S.REVERSE_TRANSCRIPTION for x in g.steps.values()),True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'s41587-020-0712-z.html.txt'
    if have(c,Source(p,'https://www.nature.com/articles/s41587-020-0712-z')):
        t=''.join(p.read_text().split());c('Selected literal probe example in published caption','NAI-N3(100mMfor5minsor25mins)modifiedTetrahymenaRNA' in t,True);c('Filter remains computational','longerthan50%ofannotatedlengths' in t,True)
    c('Private incomplete chemistry',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
