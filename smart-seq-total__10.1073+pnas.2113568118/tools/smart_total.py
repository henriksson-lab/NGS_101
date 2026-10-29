"""Smart-seq-total v1/v2: artificial RNA tails and typed degradable TSO chemistry."""
from dataclasses import dataclass,replace
from pathlib import Path
import re,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Construct,Segment,feature
import nextera as nx
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem,SourceDiscrepancy,SourceClaim,Highlight
UDG_CATALYSIS='NEB M0280 official current product description,2026-10-29: monofunctional glycosylase hydrolyzes deoxyuridine N-glycosidic bond; current formulation/units not substituted for2021 reagent'
SOURCE='Isakova et al.2021 DOI10.1073/pnas.2113568118; published PMC8713755 Methods111–133, SI FigS2/TableS2 and Supplementary Protocol'
TSO_WRITTEN={1:'dUCGdUCGGCAGCGdUCAGdUdUGdUAdUCAACdUCAGACAdUrGrG+G',2:'AdUGGCdUCGGAGAdUGdUGdUAdUAAGAGACAGdUCdUrGrG+G'}

@dataclass(frozen=True)
class Residue:
    base: str
    sugar: str
    abasic: bool=False
    def __post_init__(self):
        if self.base not in 'ACGTU' or self.sugar not in {'DNA','RNA','LNA'}:raise ValueError('invalid nucleotide')
        if self.abasic and (self.base,self.sugar)!=('U','DNA'):raise ValueError('only selected dU glycosylase product represented')
    def udg(self):return replace(self,abasic=True) if (self.base,self.sugar)==('U','DNA') else self

@dataclass(frozen=True)
class TSO:
    residues: tuple[Residue,...]
    biotin5: bool=True
    backbone_retained: bool=True
    def __post_init__(self):
        if not self.residues or not self.backbone_retained:raise ValueError('TSO backbone cannot disappear through base excision')
    @classmethod
    def source(cls,version):
        raw=TSO_WRITTEN[version];tokens=re.findall(r'dU|rG|\+G|[ACGT]',raw)
        if ''.join(tokens)!=raw:raise ValueError('unparsed source modification')
        return cls(tuple(Residue('U','DNA') if t=='dU' else Residue('G','RNA') if t=='rG' else Residue('G','LNA') if t=='+G' else Residue(t,'DNA') for t in tokens))
    def udg(self):return replace(self,residues=tuple(r.udg() for r in self.residues))
    @property
    def abasic_positions(self):return tuple(i for i,r in enumerate(self.residues) if r.abasic)
    def complementary_copy(self):
        if self.abasic_positions:raise ValueError('abasic-template traversal requires a separate specified path')
        comp={'A':'T','C':'G','G':'C','T':'A','U':'A'}
        return ''.join(comp[r.base] for r in reversed(self.residues))

@dataclass(frozen=True)
class TailPath:
    rna: str
    endogenous_tail: int=0
    added_tail: int=0
    def __post_init__(self):
        if not self.rna or set(self.rna)-set('ACGU') or self.endogenous_tail<0 or self.added_tail<0:raise ValueError('explicit RNA and nonnegative tail lengths required')
    def extend_selected(self,number):
        if not isinstance(number,int) or number<=0:raise ValueError('selected positive addition length required; source yield distribution unknown')
        return replace(self,added_tail=self.added_tail+number)
    @property
    def tailed_rna(self):return self.rna+'A'*(self.endogenous_tail+self.added_tail)

def oligo_dt(version):
    if version==1:
        return Construct([Segment('v1 handle','CATA'+nx.ADAPTOR_S7),Segment('dT30','T'*30),Segment('anchor','VN',placeholder=True)])
    if version!=2:raise ValueError('source versions1/2 only')
    return Construct([Segment('v2 handle',nx.S5+'AGTTGTATCAACTCAGAC'),Segment('constrained UMI','NNNHNNNHNNHNNHNN',placeholder=True,feature=feature('smart_total_v2_umi','umi','degenerate')),
       Segment('local linker','GCTGTAAC'),Segment('dT22','T'*22),Segment('anchor','VN',placeholder=True)])

def graph_for_version(version=2):
    if version not in {1,2}:raise ValueError('versions1/2 only')
    g=ProtocolGraph('Smart-seq-total v'+str(version)+' — RNA/TSO prefix',authority=SOURCE,representation='published_instance')
    prior=g.input('cell_rna','Single-cell total RNA after source72°C3min lysis; original RNA ends/circular species retained as unassigned populations',K.MOLECULAR,primary=True)
    polya=g.input('polya','E.coli poly(A) polymerase, source version-specific activity,37°C15min',K.CATALYST)
    atp=g.input('atp','ATP source version-specific concentration; incorporation length unassigned',K.MOLECULAR)
    tailed=g.product('tailed','Selected eligible linear RNA obtains an artificial3′poly(A) extension; endogenous tail origin no longer identifiable from total tail alone',K.MOLECULAR,chemical_state={'added_tail':'unassigned length','native_RNA_retained':True,'v2_capping':'reported simultaneous treatment; substrate-specific cap chemistry unassembled' if version==2 else 'not part of v1 treatment'})
    cap_inputs=()
    if version==2:
        cap=g.input('capping_enzyme','Vaccinia capping enzyme0.1U in source simultaneous tail/cap mix',K.CATALYST)
        donors=g.input('cap_donors','GTP0.05mM and SAM2.5nM; source cap chemistry/eligible RNA ends unassembled',K.MOLECULAR)
        cap_inputs=((cap,R.CATALYST),(donors,R.STRUCTURAL))
    g.react('tail','Source polyadenylation'+(' with simultaneous capping treatment; cap topology/yields unassigned' if version==2 else '')+'; no guarantee all RNA/circular species acquire an extendable tail',S.TAILING,((prior,R.TEMPLATE),(polya,R.CATALYST),(atp,R.STRUCTURAL))+cap_inputs,(tailed,))
    primer=g.input('oligo_dt','Source biotinylated anchored oligo-dT'+(' with constrained16baseUMI' if version==2 else ''),K.OLIGO,state=oligo_dt(version))
    rt=g.input('ssii','SuperScriptII,42°C90min; source buffer/betaine/MgCl2/dNTPs',K.CATALYST)
    cdna=g.product('firststrand','Selected oligo-dT-primed cDNA with RNA template retained; actual priming sites and residual RNA classes unassigned',K.MOLECULAR,chemical_state={'RNA_retained':True,'primer':oligo_dt(version)})
    g.react('rt','Source first-strand synthesis; no undocumented RNaseH removal',S.REVERSE_TRANSCRIPTION,((tailed,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST)),(cdna,))
    tso=g.input('tso','Source5′biotin TSO; TableS2 dU/RNA/LNA identity retained',K.OLIGO,chemical_state=TSO.source(version))
    switched=g.product('switched','Conditional TSO-templated cDNA extension; original RNA and free/extended TSO populations retained',K.MOLECULAR,chemical_state={'RNA_retained':True,'tso':TSO.source(version),'copied_TSO_DNA':TSO.source(version).complementary_copy(),'copy_scope':'selected complete TSO-templated path, efficiency/truncation unassigned'})
    g.react('switch','Source template switching; typed TSO is template, not a DNA molecule silently substituted by its sequence',S.TEMPLATE_SWITCHING,((cdna,R.TEMPLATE),(tso,R.STRUCTURAL),(rt,R.CATALYST)),(switched,))
    udg=g.input('udg','UDG '+('1U37°C30min' if version==1 else '2U37°C60min'),K.CATALYST)
    excised=g.product('abasic_tso','Selected dU base-excision products; TSO backbone retained with AP sites and liberated uracil; cDNA/RNA remain in reaction',K.MOLECULAR,chemical_state={'tso':TSO.source(version).udg(),'liberated_bases':'uracil from selected dU sites; not yet physically removed','cDNA':'copied canonical DNA retained','RNA_retained':True})
    g.react('udg_excision','Source calls this TSO digestion; local model resolves base excision only, backbone breakdown not inferred',S.DIGESTION,((switched,R.TEMPLATE),(udg,R.CATALYST)),(excised,))
    g.mark_terminal(Terminal(excised.id,Disposition.UNRESOLVED,'Remaining TSO/AP processing and source PCR/index/tagmentation branches, physical cleanup, DASH and acquired read plan unassembled. v2 simultaneous capping not yet reconstructed; no universal5′end eligibility assumed.'))
    g.declare_reviewed_scope(['cell_rna','polya','atp','tail','oligo_dt','ssii','rt','tso','switch','udg','udg_excision'],SOURCE)
    if version==2:g.declare_reviewed_scope(['capping_enzyme','cap_donors'],SOURCE)
    g.add_source_item(SourceItem('unreviewed_remaining_preparation','Version-specific capping, cDNA/PCR/index branches, DASH and acquired read geometry','Published Methods123–133; SI Supplementary Protocol and TableS2'))
    if version==1:
        g.annotate(Highlight('tso','Source TSO notation: explicit dU table used',
            'Published main Methods117/PDFp6 prints UCGUCGGCAGCGUCAGUUGUAUCAACUCAGACAUrGrG+G, whereas SI TableS2p37 explicitly marks internal deoxyU as dU. Bare U in the main oligo string does not by itself specify RNA sugar and can abbreviate dU; no RNA-versus-DNA conflict or publication error is established. Local ordered residues follow the explicit TableS2 notation, without claiming an independently verified order sheet.'))
    if version==2:
        g.annotate(SourceDiscrepancy('abasic_tso','Published v2 custom Read1 values differ',
            (SourceClaim('Published main Methods133/PDFp7','TCGTCGGCAGCGTCAGTTGTATCAACTCAGAC'),SourceClaim('Published SI TableS2 originalp37','TCGGCAGCGTCAGTTGTATCAACTCAGACATGGG')),
            'No final read plan or primer geometry is declared; neither value is silently substituted.',
            'Both original PDF locations inspected; v2 row in SI repeats the v1 custom Read1 value.',
            'Actual acquired primer/reagent identity and branch applicability remain unresolved.'))
    return g

def protocol_graph():return graph_for_version(2)
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();t=TSO.source(2);u=t.udg()
    c('UDG preserves backbone residue positions and terminal RNA/LNA', (len(u.residues),u.residues[-3:]),(len(t.residues),t.residues[-3:]))
    c('only DNA uracils become abasic',all((r.base,r.sugar)==('U','DNA') for r in u.residues if r.abasic),True)
    c('RNA uracil is not silently made UDG-susceptible',Residue('U','RNA').udg().abasic,False)
    try:u.complementary_copy();bad=False
    except ValueError:bad=True
    c('AP traversal cannot masquerade as ordinary complete copying',bad,True)
    a=TailPath('CGU',2).extend_selected(3);b=TailPath('CGU',0).extend_selected(5)
    c('same acquired tail cannot determine endogenous polyA status',(a.tailed_rna==b.tailed_rna,a.endogenous_tail!=b.endogenous_tail),(True,True))
    c('v2 UMI keeps constrainedH synthesis',next(s.top for s in oligo_dt(2) if s.feature),'NNNHNNNHNNHNNHNN')
    original=Source(sources_root()/'smart-seq-total__10.1073+pnas.2113568118'/'pnas.2113568118_supp_pnas.2113568118.sapp.pdf.txt','SI TableS2, rendered originalp37 inspected')
    if have(c,original,label='original oligos'):
        text=''.join(original.path.read_text().split())
        for v in (1,2):c('sourceTSOv'+str(v)+' full modification spelling',TSO_WRITTEN[v] in text,True)
        c('sourcev2 anchored UMI primer sequence',''.join(s.top for s in oligo_dt(2)) in text,True)
    c.report()
if __name__=='__main__':selftest()
