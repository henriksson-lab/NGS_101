"""Holo-seq: source-local post-library carrier-DNA digestion branches."""
from dataclasses import dataclass
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SourceDiscrepancy,SourceClaim
SOURCE='Xiao et al.2018 DOI10.1186/s13059-018-1553-7; published PMC6193298 Methods181–197; Additionalfile4 step-by-step library recipes; FigS1 original carrier-template image'

@dataclass(frozen=True)
class SelectedReferenceCut:
    """Conditional reference partition only; NotI stagger/end chemistry not reconstructed."""
    dna: str
    position: int
    origin: str
    def __post_init__(self):
        if not self.dna or set(self.dna)-set('ACGT') or type(self.position) is not int or not 0<self.position<len(self.dna):raise ValueError('explicit DNA with internal selected reference boundary required')
        if self.origin not in {'carrier_derived','endogenous','late_dna_carrier'}:raise ValueError('declared molecular origin required')
    @property
    def reference_pieces(self):return self.dna[:self.position],self.dna[self.position:]
    @property
    def inherited_outer_handles(self):return (frozenset({'left'}),frozenset({'right'}))

def graph_for_branch(branch='mrna'):
    if branch not in {'mrna','total','small'}:raise ValueError('source mRNA/total/small branches only')
    g=ProtocolGraph('Holo-seq — '+branch+' carrier-DNA removal boundary',authority=SOURCE,representation='published_instance')
    prior=g.input('library_dna',('Post-PCR small-RNA library' if branch=='small' else 'Source adapter-ligated '+branch+' cDNA')+' containing native and carrier-derived DNA; full upstream kit architecture unassembled',K.MOLECULAR,primary=True)
    scope=['library_dna']
    def user(key,x):
        enzyme=g.input(key+'_enzyme','Source USER enzyme37°C15min; exact historical substrate/cleavage geometry unassembled',K.CATALYST)
        out=g.product(key+'_product','Source USER-treated DNA mixture; cleavage/degradation products and unreacted material not physically separated yet',K.MOLECULAR)
        g.react(key,'Source USER treatment at branch-specific position',S.DIGESTION,((x,R.TEMPLATE),(enzyme,R.CATALYST)),(out,));scope.extend([key+'_enzyme',key]);return out
    def cleanup(key,x,late=False,ratio='1×'):
        inputs=[(x,R.TEMPLATE)]
        if late:
            carrier=g.input(key+'_dna_carrier','Source100ng~2kb DNA carrier; sequence/ends/recovery unreviewed',K.MOLECULAR)
            inputs.append((carrier,R.STRUCTURAL));scope.append(key+'_dna_carrier')
        beads=g.input(key+'_beads',ratio+'AMPureXP; source physical purification',K.SELECTION)
        keep=g.product(key+'_kept','Physically recovered DNA mixture; native/remaining carrier/late DNA fractions unquantified',K.MOLECULAR)
        waste=g.product(key+'_waste','Physical bead-binding/wash rejects; exact fragment recovery unassigned',K.WASTE)
        g.react(key,'Source bead cleanup; no species disappears before separation',S.CLEANUP,tuple(inputs)+((beads,R.SELECTION),),(keep,waste));g.discard(waste,'Physical purification rejects only');scope.extend([key+'_beads',key]);return keep
    if branch=='mrna':prior=user('user_before_noti',prior)
    if branch!='small':prior=cleanup('ligation_cleanup',prior,late=True)
    enzyme=g.input('noti','NotI2µl in source CutSmart50µl,37°C2h',K.CATALYST)
    cut=g.product('noti_digest','NotI-treated DNA containing cut reference pieces, uncut molecules and unassigned carrier/native recovery; no direct RNA excision',K.MOLECULAR,chemical_state={'substrate':'library DNA','cut_products':'both physical fragments retained until cleanup','specificity_scope':'native NotI-site DNA can also be affected'})
    g.react('digest_noti','Source carrier-library cleavage; do not assert perfect selective depletion',S.DIGESTION,((prior,R.TEMPLATE),(enzyme,R.CATALYST)),(cut,));scope+=['noti','digest_noti']
    prior=cleanup('digest_cleanup',cut,late=branch=='small')
    if branch=='total':prior=user('user_after_noti',prior)
    g.mark_terminal(Terminal(prior.id,Disposition.UNRESOLVED,'Source post-removal PCR'+('/USER mixture' if branch=='total' else '')+' or small-library sequencing handoff; exact historical full molecular/USER/adapter/read chemistry unassembled. Carrier depletion efficiency and residual fractions not predicted.'))
    g.declare_reviewed_scope(scope,SOURCE)
    g.add_source_item(SourceItem('unreviewed_full_library','Historical E7490/E7530/E7420/E7300 polymers, upstream reactions, USER substrate topology and final read plan','Published main Methods177–197 and Additionalfile4 cited kit procedures'))
    if branch=='small':
        g.annotate(SourceDiscrepancy('library_dna','Published small-RNA PCR counts differ',
           (SourceClaim('Main Methods193','15cycles before NotI'),SourceClaim('Additionalfile4 small-RNA PCR block','22cycles before NotI')),
           'Model begins after PCR and assigns no cycle count.',
           'Both original published locations inspected.',
           'Applicability to individual datasets/procedure revisions unresolved; no numerical repair inferred.'))
    return g

def protocol_graph():return graph_for_branch()
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=SelectedReferenceCut('ACGTGCGGCCGCTT',6,'carrier_derived')
    c('selected reference cut preserves all DNA bases',''.join(p.reference_pieces),p.dna)
    c('cut pieces retain opposite original outer handles',p.inherited_outer_handles,(frozenset({'left'}),frozenset({'right'})))
    q=SelectedReferenceCut(p.dna,p.position,'endogenous');c('origin label does not protect a susceptible endogenous fragment',p.reference_pieces,q.reference_pieces)
    c('mRNA USER precedes carrier digest',tuple(graph_for_branch('mrna').steps)[:3],('user_before_noti','ligation_cleanup','digest_noti'))
    c('total-RNA USER follows actual carrier-removal cleanup',tuple(graph_for_branch('total').steps)[-2:],('digest_cleanup','user_after_noti'))
    c('small branch lateDNA carrier is a declared material','digest_cleanup_dna_carrier' in graph_for_branch('small').materials,True)
    original=Source(sources_root()/'holo-seq__10.1186+s13059-018-1553-7'/'s13059-018-1553-7_supp_13059_2018_1553_MOESM4_ESM.docx.txt','Additionalfile4')
    if have(c,original,label='source latecarrier and smallbranch checks'):
        text=''.join(original.path.read_text().split())
        c('source lateDNA carrier retained in recipe','Add100ng~2KbDNAcarrier' in text,True)
        c('source detailed smallRNA count kept distinct','22cyclesof' in text,True)
    c.report()
if __name__=='__main__':selftest()
