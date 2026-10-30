"""Provisional RBS-seq main/SI model: physical RNA versus complementary copying."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
SOURCE='Khoddami et al.2019 DOI10.1073/pnas.1817334116; published main/SI Methodspp3–6/9–11; linked correction10.1073/pnas.1905628116 CONTENTS UNREAD, applicability provisional'

@dataclass(frozen=True)
class RNA:
    residues: tuple[str,...]
    origins: tuple[str,...]=()
    heat_mg_conditioned: bool=False
    def __post_init__(self):
        if not self.origins:object.__setattr__(self,'origins',self.residues)
        allowed={'A':{'A'},'C':{'C','U'},'G':{'G'},'U':{'U'},'m5C':{'m5C'},'m1A':{'m1A','m6A'},'m6A':{'m6A'},'Psi':{'Psi','PsiBS'}}
        if not self.residues or len(self.residues)!=len(self.origins) or any(o not in allowed or n not in allowed[o] for o,n in zip(self.origins,self.residues)):raise ValueError('explicit retained RNA residue lineage and selected source chemistry required')
    def bs_processed(self):
        changes={'C':'U','m1A':'m6A','Psi':'PsiBS'}
        return replace(self,residues=tuple(changes.get(r,r) for r in self.residues))
    def condition_heat_mg(self):return replace(self,heat_mg_conditioned=True)
    def copy_selected(self,*,skip=frozenset(),native_m1a_choices=None):
        if not isinstance(skip,frozenset):raise ValueError('explicit immutable skipped template positions required')
        if any(type(i) is not int or not 0<=i<len(self.residues) or self.residues[i]!='PsiBS' for i in skip):raise ValueError('this selected bypass class concerns explicit adduct positions only')
        if skip and not self.heat_mg_conditioned:raise ValueError('selected source bypass explanation requires post-BS heat/Mg qualification; other assay regimes outside model')
        choices={} if native_m1a_choices is None else native_m1a_choices
        comp={'A':'T','C':'G','G':'C','U':'A','m5C':'G','m6A':'T','Psi':'A','PsiBS':'A'};bases=[]
        for i,r in enumerate(self.residues):
            if i in skip:continue
            if r=='m1A':
                if choices.get(i) not in {'A','C','G','T'}:raise ValueError('native m1A first-strand choice must be explicit, not predicted')
                bases.append(choices[i])
            else:bases.append(comp[r])
        return Copy(self,''.join(reversed(bases)),skip)

@dataclass(frozen=True)
class Copy:
    template: RNA
    dna5to3: str
    skipped: frozenset[int]
    def __post_init__(self):
        if set(self.dna5to3)-set('ACGT') or len(self.dna5to3)+len(self.skipped)!=len(self.template.residues):raise ValueError('selected complementary copy must retain template-coordinate accounting')
    @property
    def reference_oriented(self):return revcomp(self.dna5to3)

def unreviewed_adapter_end(which):
    sc=Scene();sc.strand('adapter',[Segment('unassembled source end','[unreviewed '+which+' end]',placeholder=True)])
    return AdapterLigationGeometry('unavailable','unavailable',(sc,),'Historical source adapter sequence/modifications/end geometry unassembled; no assertion that source did not report them.')

def graph_for_sample(fraction='polyA',bs=True):
    if fraction not in {'polyA','ribo_depleted','small'}:raise ValueError('three source RNA fractions only')
    g=ProtocolGraph('RBS-seq — '+fraction+(' BS' if bs else ' NBS'),authority=SOURCE,representation='published_instance')
    prior=g.input('prepared_rna','Source separately prepared '+fraction+' RNA fraction; source fragmenting/selection front-end not fully assembled',K.MOLECULAR,primary=True)
    scope=['prepared_rna']
    def clean(key,x,label,procedure,kind=K.MOLECULAR):
        reagent=g.input(key+'_reagent',procedure,K.SELECTION);out=g.product(key+'_kept',label,kind);waste=g.product(key+'_waste','Physical '+key+' rejects; exact composition/recovery unassigned',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(reagent,R.SELECTION)),(out,waste));g.discard(waste,'Actual physical purification fractions');scope.extend([key+'_reagent',key]);return out
    if bs:
        formamide=g.input('formamide_heat','45µlRNA+240µlformamide,95°C5min then ice',K.CONDITION)
        den=g.product('denatured','Source denatured RNA; refolding/survival fractions unassigned',K.MOLECULAR)
        g.react('denature','Source RNA denaturation before BS',S.THERMAL_TREATMENT,((prior,R.TEMPLATE),(formamide,R.CONDITION)),(den,));scope+=['formamide_heat','denature']
        bisulfite=g.input('bisulfite','312µl5Mbisulfite stock+3µl100mMhydroquinone,600µltotal;50°C16h',K.CONDITION)
        changed=g.product('bs_reaction','RNA and reaction products; sulfonated/deaminated/adduct intermediates and survival fractions unresolved',K.MOLECULAR)
        g.react('bs_treat','Source BS reaction; no physical RNA deletion assigned',S.CHEMICAL_CONVERSION,((den,R.TEMPLATE),(bisulfite,R.CONDITION)),(changed,));scope+=['bisulfite','bs_treat']
        prior=clean('nap10',changed,'NAP-10-recovered RNA population','Source equilibrated NAP-10 column and water elution')
        tris=g.input('tris','Mix0.5mlsample+0.5ml2MTris-stock pH9;37°C2h',K.CONDITION)
        des=g.product('desulfonated','Selected surviving RNA C→U,m5Cretained,m1A→m6A,PsiBS classes; RNA backbone remains, rates/unreacted products unassigned',K.MOLECULAR)
        g.react('desulfonate','Source complete BS/desulfonation workflow; exact site reaction timing not overassigned',S.CHEMICAL_CONVERSION,((prior,R.TEMPLATE),(tris,R.CONDITION)),(des,));scope+=['tris','desulfonate']
        prior=clean('ethanol_recovery',des,'Ethanol-recovered source RNA','Source sodiumacetate/GlycoBlue/ethanol;≥30min−80°C, pellet and70%ethanol wash')
    ap=g.input('antarctic_phosphatase','Source Antarctic phosphatase37°C30min,65°C5min',K.CATALYST)
    repaired=g.product('dephosphorylated','Source RNA after phosphatase/heat; exact end/site class fractions unassigned',K.MOLECULAR)
    g.react('repair3','Source end repair starts with phosphatase',S.DEPHOSPHORYLATION,((prior,R.TEMPLATE),(ap,R.CATALYST)),(repaired,));scope+=['antarctic_phosphatase','repair3']
    pnk=g.input('pnk','Source PNK37°C60min',K.CATALYST);atp=g.input('atp','Source ATP phosphate donor',K.MOLECULAR)
    phosphorylated=g.product('phosphorylated','Source end-repaired RNA with reagent/carryover material retained',K.MOLECULAR)
    g.react('repair5','Source PNK+ATP treatment',S.PHOSPHORYLATION,((repaired,R.TEMPLATE),(pnk,R.CATALYST),(atp,R.STRUCTURAL)),(phosphorylated,));scope+=['pnk','atp','repair5']
    prior=clean('rna_cleanup',phosphorylated,'Source RNeasyMinElute recovered RNA','Historical RNeasyMinElute protocol')
    a3=g.input('adapter3','Illumina v1.5sRNA3′adapter; exact historical polymer/modifications unreviewed',K.ADAPTER);lig2=g.input('ligase2','T4RNA ligase2 truncated; sourcebuffer+MgCl2,22°C1h',K.CATALYST)
    linked3=g.product('linked3','Source3′adapter-ligated RNA and residual substrates; no invented RNA deletion',K.MOLECULAR)
    g.react('ligate3','Source70°C2min preanneal then3′ligation',S.LIGATION,((prior,R.TEMPLATE),(a3,R.STRUCTURAL),(lig2,R.CATALYST)),(linked3,),geometry=unreviewed_adapter_end('v1.5sRNA3′'));scope+=['adapter3','ligase2','ligate3']
    a5=g.input('adapter5','Source Illumina5′adapter; historical sequence/ends unreviewed',K.ADAPTER);lig1=g.input('ligase1','Source T4RNA ligase withATP,20°C1h',K.CATALYST)
    linked5=g.product('linked5','Source doubly-adapter-ligated RNA/carryover population',K.MOLECULAR)
    g.react('ligate5','Source separately70°C2min-heated5′adapter then ligation',S.LIGATION,((linked3,R.TEMPLATE),(a5,R.STRUCTURAL),(lig1,R.CATALYST)),(linked5,),geometry=unreviewed_adapter_end('5′adapter'));scope+=['adapter5','ligase1','ligate5']
    primer=g.input('rt_primer','Source barcodedRTprimer; exact member, barcode/end/priming geometry unassembled',K.OLIGO);rt=g.input('ssii','SuperScriptII per source supplier procedure',K.CATALYST)
    hybrid=g.product('hybrid','Source cDNA with original RNA retained; conditional Psi bypass changes cDNA only, native m1A copy distribution unassigned',K.MOLECULAR,chemical_state={'RNA_retained':True,'copy_scope':'selected classes; no RNA residue deletion','Psi_mechanism':'heat/Mg-qualified bypass explanation; exact library exposure/conformation unassembled'})
    g.react('rt','Source first-strand synthesis; copying outcomes not physical RNA base substitution',S.REVERSE_TRANSCRIPTION,((linked5,R.TEMPLATE),(primer,R.STRUCTURAL),(rt,R.CATALYST)),(hybrid,));scope+=['rt_primer','ssii','rt']
    f=g.input('pcr_forward','Source Illumina forward primer; historical sequence unreviewed',K.OLIGO);r=g.input('pcr_reverse','Source Illumina reverse primer; historical sequence unreviewed',K.OLIGO);poly=g.input('phusion','Source PhusionHF,15cycles',K.CATALYST)
    amplified=g.product('amplified','Source amplified cDNA plus original/carryover molecular material',K.MOLECULAR)
    g.react('pcr','Source15cycle PCR; exact template/primer geometry unassembled',S.PCR,((hybrid,R.TEMPLATE),(f,R.STRUCTURAL),(r,R.STRUCTURAL),(poly,R.CATALYST)),(amplified,),geometry=PCRGeometry.unavailable(Scene(),reason='Historical original full adapters, barcodedRT/PCR/run-primer sites unassembled.'));scope+=['pcr_forward','pcr_reverse','phusion','pcr']
    purified=clean('ampure',amplified,'Source AMPure-recovered library','Source AMPureXP cleanup',K.LIBRARY)
    alkali=g.input('naoh','Source NaOH flowcell-preparation denaturation; detailed mixture unassembled',K.CONDITION)
    final=g.product('library','Source denatured DNA-library population, both strand/carryover fates unassigned',K.LIBRARY)
    g.react('flowcell_denature','Source library denaturation before separate-lane sequencing',S.GENERIC,((purified,R.TEMPLATE),(alkali,R.CONDITION)),(final,),changes_molecule=False,reason='Source NaOH-mediated noncovalent strand denaturation; no heat or covalent transformation assigned.');scope+=['naoh','flowcell_denature']
    g.sequence(final,SequencingPlan(acquired_roles=('Read 1','Read 2'),read_lengths={'Read 1':101,'Read 2':101},limitation='Published main/SI report HiSeq2000PE101; full original library/run-primer/index geometry and correction applicability unassembled.'))
    g.declare_reviewed_scope(scope,SOURCE)
    g.add_source_item(SourceItem('unread_correction','Actual linked correction contents remain unread','https://doi.org/10.1073/pnas.1905628116'))
    g.add_source_item(SourceItem('unreviewed_library_geometry','Original historical full oligos, DatasetS10 and exact heat/Mg/copy/read geometry','Published main/SI and directly cited historical protocols'))
    return g

def protocol_graph():return graph_for_sample()
def selftest():
    from checks import Check
    c=Check();rna=RNA(('C','m5C','m1A','Psi','G'));treated=rna.bs_processed().condition_heat_mg();copied=treated.copy_selected(skip=frozenset({3}))
    c('one RNA can retain distinct treated chemical classes',treated.residues,('U','m5C','m6A','PsiBS','G'))
    c('physical RNA stays intact while selected cDNA is shorter',(len(copied.template.residues),len(copied.dna5to3)),(5,4))
    c('complementary firstDNA and reference-oriented report differ',(copied.dna5to3,copied.reference_oriented),('CTGA','TCAG'))
    c('RNA modification origin preserved through copying',copied.template.origins,rna.origins)
    try:rna.bs_processed().copy_selected(skip=frozenset({3}));bad=False
    except ValueError:bad=True
    c('selected bypass path requires explicit heat/Mg qualification',bad,True)
    try:rna.copy_selected();bad=False
    except ValueError:bad=True
    c('native m1A copying outcome is not silently predicted',bad,True)
    c('explicit native m1A choice permits a selected mismatch',RNA(('m1A',)).copy_selected(native_m1a_choices={0:'C'}).reference_oriented,'G')
    c.report()
if __name__=='__main__':selftest()
