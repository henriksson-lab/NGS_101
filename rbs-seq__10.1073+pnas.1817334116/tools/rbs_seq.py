"""Provisional RBS-seq main/SI model: physical RNA versus complementary copying."""
from dataclasses import dataclass,replace
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
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
        formamide=g.input('formamide','Actual240µldeionizedformamide molecular contribution',K.MOLECULAR)
        den=g.product('denatured','Source denatured RNA; refolding/survival fractions unassigned',K.MOLECULAR)
        g.react('denature','Source RNA denaturation before BS',S.THERMAL_TREATMENT,((prior,R.TEMPLATE),(formamide,R.STRUCTURAL),(g.materials['formamide_heat'],R.CONDITION)),(den,));scope+=['formamide_heat','formamide','denature']
        bisulfite=g.input('bisulfite','312µl5Mbisulfite stock+3µl100mMhydroquinone,600µltotal;50°C16h',K.CONDITION)
        bs_reagent=g.input('bs_reagent','Actual312µlfresh5Msodiumbisulfite stock,pH5',K.MOLECULAR)
        hydroquinone=g.input('hydroquinone','Actual3µl100mMhydroquinone',K.MOLECULAR)
        changed=g.product('bs_reaction','RNA and reaction products; sulfonated/deaminated/adduct intermediates and survival fractions unresolved',K.MOLECULAR)
        g.react('bs_treat','Source BS reaction; no physical RNA deletion assigned',S.CHEMICAL_CONVERSION,((den,R.TEMPLATE),(bisulfite,R.CONDITION),(bs_reagent,R.STRUCTURAL),(hydroquinone,R.STRUCTURAL)),(changed,));scope+=['bisulfite','bs_reagent','hydroquinone','bs_treat']
        prior=clean('nap10',changed,'NAP-10-recovered RNA population','Source equilibrated NAP-10 column and water elution')
        portions=tuple(g.product('des_portion_'+str(i),'Actual0.5mlportion of1mlNAP10eluate for separate desulfonation',K.MOLECULAR)for i in range(2))
        g.react('des_allocate','Actual1ml eluate splits into two0.5ml tubes BEFORE desulfonation',S.SPLIT,((prior,R.TEMPLATE),),portions)
        recovered=[]
        for i,portion in enumerate(portions):
            k='tube_'+str(i)
            tris=g.input(k+'_tris','Actual0.5ml2MTris pH9 stock added to this0.5mlRNA portion',K.MOLECULAR)
            conditions=g.input(k+'_condition','Actual37°C2h desulfonation',K.CONDITION)
            ds=g.product(k+'_desulfonated','Source surviving RNA C→U,m5Cretained,m1A→m6A,ΨBS selected classes; actual rates/residuals unassigned',K.MOLECULAR)
            g.react(k+'_desulfonate','Actual independent desulfonation tube',S.CHEMICAL_CONVERSION,((portion,R.TEMPLATE),(tris,R.STRUCTURAL),(conditions,R.CONDITION)),(ds,))
            acetate=g.input(k+'_acetate','Actual111µl3Msodiumacetate pH5.5',K.MOLECULAR)
            glycoblue=g.input(k+'_glycoblue','Actual20µl15mg/mlGlycoBlue carrier',K.MOLECULAR)
            mix=g.product(k+'_mix','Actual1mlRNA/Tris+111µlacetate+20µlGlycoBlue preparation',K.MOLECULAR)
            g.react(k+'_mix_prepare','Actual precipitation reagents added before three-way split',S.GENERIC,((ds,R.TEMPLATE),(acetate,R.STRUCTURAL),(glycoblue,R.STRUCTURAL)),(mix,),changes_molecule=False)
            precip=tuple(g.product(k+'_part_'+str(j),'Actual377µlprecipitation portion, not a cloned reaction product',K.MOLECULAR)for j in range(3))
            g.react(k+'_split','Actual independent377µl portions; three tubes from each desulfonated preparation',S.SPLIT,((mix,R.TEMPLATE),),precip)
            for j,part in enumerate(precip):
                q=k+'_precip_'+str(j)
                ethanol=g.input(q+'_ethanol','Actual950µlabsoluteethanol then1ml70%ethanolwash for this tube',K.MOLECULAR)
                environment=g.input(q+'_condition','Actual≥30min−80°C;21,000g20min4°C;wash/airdry;resuspend10µlwater',K.CONDITION)
                keep=g.product(q+'_kept','Actual10µl recoveredRNA from this independent precipitation tube',K.MOLECULAR)
                waste=g.product(q+'_waste','Actual precipitation/washsupernatants from this tube',K.WASTE)
                g.react(q,'Actual independent source precipitation and recovery',S.CLEANUP,((part,R.TEMPLATE),(ethanol,R.STRUCTURAL),(environment,R.CONDITION)),(keep,waste));g.discard(waste,'Actual removed precipitation/wash fractions');recovered.append(keep)
        prior=g.product('ethanol_pooled','Actual six10µl recoveredRNA contributions pooled to60µl stock',K.MOLECULAR)
        g.react('ethanol_pool','Pool actual six physically produced RNA fractions after recovery',S.POOL,tuple((x,R.TEMPLATE)for x in recovered),(prior,))
        g.declare_reviewed_scope(TUBE_TARGETS,SOURCE)
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

# Independent curator inventory: actual two desulfonation and six precipitation source tubes.
TUBE_TARGETS=tuple("""des_allocate des_portion_0 des_portion_1
 tube_0_tris tube_0_condition tube_0_desulfonate tube_0_desulfonated tube_0_acetate tube_0_glycoblue tube_0_mix_prepare tube_0_mix tube_0_split tube_0_part_0 tube_0_part_1 tube_0_part_2
 tube_1_tris tube_1_condition tube_1_desulfonate tube_1_desulfonated tube_1_acetate tube_1_glycoblue tube_1_mix_prepare tube_1_mix tube_1_split tube_1_part_0 tube_1_part_1 tube_1_part_2
 tube_0_precip_0_ethanol tube_0_precip_0_condition tube_0_precip_0 tube_0_precip_0_kept tube_0_precip_0_waste
 tube_0_precip_1_ethanol tube_0_precip_1_condition tube_0_precip_1 tube_0_precip_1_kept tube_0_precip_1_waste
 tube_0_precip_2_ethanol tube_0_precip_2_condition tube_0_precip_2 tube_0_precip_2_kept tube_0_precip_2_waste
 tube_1_precip_0_ethanol tube_1_precip_0_condition tube_1_precip_0 tube_1_precip_0_kept tube_1_precip_0_waste
 tube_1_precip_1_ethanol tube_1_precip_1_condition tube_1_precip_1 tube_1_precip_1_kept tube_1_precip_1_waste
 tube_1_precip_2_ethanol tube_1_precip_2_condition tube_1_precip_2 tube_1_precip_2_kept tube_1_precip_2_waste ethanol_pool ethanol_pooled""".split())

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
    g=protocol_graph()
    c('Actual two RNA portions before desulfonation',tuple(g.steps['des_allocate'].outputs),('des_portion_0','des_portion_1'))
    c('Six physically produced10µl fractions pooled',len(g.steps['ethanol_pool'].inputs),6)
    c('Each pooled component has its own source precipitation producer',len({g.producer_of(x.material).id for x in g.steps['ethanol_pool'].inputs}),6)
    c('Bisulfite is a molecular reactant',g.materials['bs_reagent'].kind,K.MOLECULAR)
    c('Actual correction source item remains uncovered','unread_correction' in g.source_coverage,False)
    c('Actual full geometry source item remains uncovered','unreviewed_library_geometry' in g.source_coverage,False)
    c.report(exit_on_fail=True)
if __name__=='__main__':selftest()
