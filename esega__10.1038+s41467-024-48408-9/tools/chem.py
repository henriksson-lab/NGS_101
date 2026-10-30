"""Source-owned private partial ESEGA ss-ZP1 deamination/five-triphosphate route.

Synthetic Z/P remain chemical states: no ACGT-normalization or supported-topology claim.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature,revcomp
from illumina import TRUSEQ_READ1,NEBNEXT_ARM_READ2
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from pcr import PCRGeometry

TITLE='ESEGA selected ss-ZP1 deaminated five-triphosphate route'
AUTHORITY='Published PMC11094070 main129–145/233–251/277–295 and published SI Table1 p5, Table2 p6; two originals. Ref/control, ds-separation and upstream six-triphosphate synthesis are separate branches.'
ZP1='TAAGATGAGAGTTGAGGAGAGTTACGTGZACGCPTPGTCAZCACAGTATAGTAGTGTAAGTAGATAGTGGA'
A1='ATGATGAAGCCTGAAGTTCCAGATTCCCATCAGTGATTCTGCGCTGATCATGCTGTATGCCATGTGCATCCTCACCTAAGATGAGAGTTGAGGAGAGTTA'
B1='GGTGTTAGTCTTGCTGTTGCTGTGACCTCCGAGAAGTGGGTGGCTTATCCAATCGGTAGTAGCGGTGACCTTCTGCTCCACTATCTACTTACACTACTAT'
A_CODES=('GAGTCCAT','TAGTCGAG','CCGATACG','TTGCACTG','GAAGATTC')
B_CODES=('CTTCACGA','TGAAGCTC','GAGATGCT','GCAGTATC','AACTCGTC')
A_TAIL='ATGATGAAGCCTGAAGTTCCA'
B_TAIL='GGTGTTAGTCTTGCTGTTGCT'
B_HANDLE=revcomp(NEBNEXT_ARM_READ2)+'T'
def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[source molecular geometry unassembled]',placeholder=True)]);return sc
def primer_scene(name):
    sc=Scene()
    if name in {'A1','B1'}:parts=[Segment(name,{'A1':A1,'B1':B1}[name])]
    else:
        handle,tail=(TRUSEQ_READ1,A_TAIL) if name=='A1-Bod' else (B_HANDLE,B_TAIL)
        parts=[Segment('source handle',handle),Segment('sample index','N'*8,placeholder=True,feature=feature('esega_'+name.lower(),'sample_index','unknown'),note='Table1 five actual barcodes; selected experiment member unassigned, not a random UMI'),Segment('source complete binding tail',tail)]
    sc.strand('oligo',parts);return sc
def barcode_members(which):
    h,t,c=(TRUSEQ_READ1,A_TAIL,A_CODES)if which=='A'else(B_HANDLE,B_TAIL,B_CODES)
    return tuple(h+b+t for b in c)

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('zp1','Source synthetic single-stranded ZP1 supplied stock; synthesis outside selected route',K.MOLECULAR,primary=True,chemical_state={'polymer':'DNA with synthetic Z/P','alphabet':'ACGTZP','source_ordering':ZP1})
    def op(key,label,kind,stocks=(),product='Source reaction mixture; incompletely processed substrate and free reagents retained until physical workup'):
        nonlocal x
        ins=[(x,R.TEMPLATE)]
        for sid,text,mk,role,cs,state in stocks:ins.append((g.input(sid.lower().replace("-","_"),text,mk,chemical_state=cs,state=state),role))
        kw={}
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Whole source-specific Z/P conversion and primer landing/product geometry unassembled; not a source absence claim.')
        if kind==S.GENERIC:kw['changes_molecule']=False
        if kind==S.MAGIC:kw['reason']='Source service finish/whole molecular geometry unassembled.'
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,ins,(y,),**kw);x=y
    def cond(key,text):return(key,text,K.CONDITION,R.CONDITION,None,None)
    def enz(key,text):return(key,text,K.CATALYST,R.CATALYST,None,None)
    def split(key,label,rest,disp=Disposition.UNRESOLVED):
        nonlocal x
        a=g.product(key+'_used',label,K.MOLECULAR);b=g.product(key+'_rest',rest,K.QC if disp==Disposition.QC else K.MOLECULAR)
        g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,b));g.mark_terminal(Terminal(b.id,disp,rest));x=a
    op('denature','Source 20ul 1uM DNA+4ul formamide,85C10min then ice',S.THERMAL_TREATMENT,(cond('denature_condition','Source formamide/temperature programme'),))
    op('deaminate','Source add100ul premix containing APOBEC,37C3h',S.CHEMICAL_CONVERSION,(enz('apobec','Source APOBEC cytidine deaminase'),cond('deam_condition','Source premix formulation not assembled here')),product='Source C-to-U treated synthetic DNA plus residual C/Z/P molecules; no quantitative completion assumed')
    op('heat_stop','Selected source direct-use alternative:90C60sec heat inactivation, without bead/column cleanup',S.QUENCHING,(cond('heat_condition','Source heat alternative; purification routes separate'),))
    op('dilute','Source dilute treated DNA before amplification; factor unassigned',S.GENERIC,(cond('dilution','Source aqueous dilution'),))
    split('template_aliquot','Source1ul diluted template used in50ul length-extension PCR','Untransferred diluted template; exact volume/fate unassigned')
    nts=('five_nts','Source Results five-triphosphate pool:0.2mM each A/C/G/T plus0.05mMdZTP, no dPTP; added dZTP volume in main279 generic premix unassigned',K.MOLECULAR,R.STRUCTURAL,{'nucleotides':('dATP','dCTP','dGTP','dTTP','dZTP'),'concentration_mM':(.2,.2,.2,.2,.05),'scope':'Selected Results route; not silently equated to generic four-standard premix'},None)
    ps=tuple((n,'Source5ul '+n+'-100 500nM',K.OLIGO,R.STRUCTURAL,{'polymer':'DNA'},primer_scene(n))for n in('A1','B1'))
    op('length_pcr','Source50ul:5ul each A1/B1,1ul template,25ul2xTaKaRaTaq pH8.9,14ul water;95C1min,10x95C15s/60C10s/72C30s. Selected five-triphosphate addition bookkeeping unassigned',S.PCR,ps+(enz('length_taq','Source fresh TaKaRaTaq pH8.9'),nts),product='Source length-extended/transliterated copy populations plus original Z/P templates and off-products; intended U->T/Z->C/P->G readout is not guaranteed complete conversion')
    split('barcode_aliquot','Source5ul from50ul length-extension product transferred to new barcoding tube','Remaining45ul length-extension mixture; fate unassigned')
    bp=tuple((n,'Source '+n+' indexed family,0.5uM; selected family member unassigned',K.OLIGO,R.STRUCTURAL,{'polymer':'DNA','source_members':barcode_members(w)},primer_scene(n))for n,w in(('A1-Bod','A'),('B1-Bod','B')))
    op('barcode_pcr','Source new TaKaRaTaq master mix;95C1min,12–20x95C15s/58C10s/72C30s, EvaGreen monitoring; second-PCR nucleotide formulation unassigned',S.PCR,bp+(enz('barcode_taq','Source fresh second-tube TaKaRaTaq master mix'),),product='Source barcoded products, residual original templates/primers and off-products before physical workup; no exact whole library asserted')
    split('page_qc','Production DNA retained','Source separate analytical15percentPAGE-urea/TBE aliquot; allocation unassigned',Disposition.QC)
    cc=g.input('gel_condition','Source2percentagarose/QIAgel extraction',K.CONDITION);kept=g.product('gel_kept','Source gel-recovered production DNA',K.MOLECULAR);w=g.product('gel_waste','Physically separated gel/workup fractions',K.WASTE)
    g.react('gel_extract','Source physical gel extraction after PCR',S.CLEANUP,((x,R.TEMPLATE),(cc,R.CONDITION)),(kept,w));g.discard(w,'Source physical extraction, no perfect molecular purity assumed.');x=kept
    op('service','Source GENEWIZ Amplicon-EZ NGS; historical service finishing, P5/P7/index acquisition geometry and instrument/read lengths unassembled',S.MAGIC)
    lib=g.product('library','Source sequenced amplicon preparation; complete library geometry unassembled',K.LIBRARY);g.react('loading','Source service loading preparation',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='Published Fig3 compares untreated/deaminated five-triphosphate transliteration. Z:P expected to read C:G; ordinary C:G deamination expected T:A. Conversion probabilities, original-template carryover, service read/index geometry remain unassembled.'))
    g.add_source_item(SourceItem('pending','Historical Amplicon-EZ finishing and exact Z/P-containing whole product/primer geometry unassembled; main279 premix vs Results five-nucleotide supplement bookkeeping unassigned.','Selected source-positive private scope; no missing-information claim or full linked-source audit.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph()
    c('Raw six-letter source substrate preserved',g.materials['zp1'].chemical_state['source_ordering'],ZP1)
    c('Synthetic bases not normalized',set(ZP1)-set('ACGT'),{'Z','P'})
    c('Heat alternative directly follows reaction',g.steps['heat_stop'].inputs[0].material,'deaminate_product')
    c('No purification invented before direct heat branch',sum(s.kind==S.CLEANUP for s in g.steps.values()),1)
    c('Five-triphosphate pool excludes dPTP',g.materials['five_nts'].chemical_state['nucleotides'],('dATP','dCTP','dGTP','dTTP','dZTP'))
    c('Second PCR consumes actual5ul transfer',g.steps['barcode_pcr'].inputs[0].material,'barcode_aliquot_used')
    c('Untransferred45ul explicit',g.terminals['barcode_aliquot_rest'].disposition,Disposition.UNRESOLVED)
    c('Two fresh source PCR catalyst additions',[(s.id,[p.material for p in s.inputs if p.role==R.CATALYST])for s in g.steps.values()if s.kind==S.PCR],[('length_pcr',['length_taq']),('barcode_pcr',['barcode_taq'])])
    c('Complete A1 family binding tail preserved',barcode_members('A')[0][-len(A_TAIL):],A_TAIL)
    c('Complete B1 family binding tail preserved',barcode_members('B')[0][-len(B_TAIL):],B_TAIL)
    c('Barcode declared as sample identifier',primer_scene('A1-Bod').strands['oligo'].segs[1].feature.role,'sample_index')
    path=sources_root()/Path(__file__).resolve().parents[1].name/'s41467-024-48408-9_supp_41467_2024_48408_MOESM1_ESM.pdf.txt'
    if have(c,Source(path,'https://doi.org/10.1038/s41467-024-48408-9')):
        compact=''.join(path.read_text().split())
        for name,seq in [('ZP1',ZP1),('A1',A1),('B1',B1)]+[(f'{w}Bod{i+1}',seq)for w in('A','B')for i,seq in enumerate(barcode_members(w))]:c('Whole source ordering '+name,seq in compact,True)
    c('Private source scope remains incomplete',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
