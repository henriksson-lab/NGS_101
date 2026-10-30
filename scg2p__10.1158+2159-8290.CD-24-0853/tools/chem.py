from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition

def unknown_scene(label='Unassembled whole molecular product'):
    sc=Scene();sc.strand('unassembled',[Segment('unassembled','[unassembled]',placeholder=True)]);return sc

def add(g,key,label,kind=K.MOLECULAR,role=R.STRUCTURAL,**kw):
    return (g.input(key,label,kind,**kw),role)
def condition(g,key,label):return add(g,key,label,K.CONDITION,R.CONDITION)
def opaque(g,key,label,x,ins=(),kind=S.MAGIC):
    y=g.product(key+'_product',label+'; whole molecular identities/junctions unassembled',K.MOLECULAR)
    kw={'reason':'Source-positive operation; whole kit composition, substrate context or product geometry unassembled, not publication absence.'} if kind==S.MAGIC else {}
    g.react(key,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);return y
def cleanup(g,key,label,x):
    y=g.product(key+'_kept','Source retained preparation; recovery/carryover unassigned',K.MOLECULAR);w=g.product(key+'_removed','Source workup fraction',K.WASTE)
    g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),condition(g,key+'_condition',label)),(y,w));g.discard(w,'Physical workup, not perfect molecular selectivity.');return y
def aliquot(g,key,label,x,remainder_label='Untransferred source stock; later fate unassigned'):
    y=g.product(key+'_chosen',label,K.MOLECULAR);r=g.product(key+'_remainder',remainder_label,K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,remainder_label));return y
def qc(g,key,label,x):
    a=g.product(key+'_analytical',label,K.MOLECULAR);k=g.product(key+'_kept','Retained production stock',K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,k));g.mark_terminal(Terminal(a.id,Disposition.QC,'Analytical aliquot not returned; amount unassigned.'));return k
TITLE='scG2P: associated-manuscript targeted RNA and DNA co-capture'
AUTHORITY='Associated NIHMS2135808/PMC12874418.1 Methods189–217 and Results97/105. Not verified equivalent to final published body; one original. No modern Tapestri constructor inherited.'
FORWARD_HANDLE='GTACTCGCAGTAGTC'
REVERSE_HANDLE='GTGATACACGACTATGAGCGCTA'
CAPTURE='/5BioTinTEG/'+REVERSE_HANDLE+'/3C6/'
def capture_oligo():
    sc=Scene();sc.strand('capture',[Segment('rna_amplicon_handle_probe',REVERSE_HANDLE)],mod5='BioTinTEG',mod3='C6');return sc
def capture_pair():
    sc=capture_oligo();sc.anneal('rna_derived_amplicon',[Segment('handle_complement',revcomp(REVERSE_HANDLE),inferred=True)],to='capture',pair=('handle_complement','rna_amplicon_handle_probe'));return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    stock=g.input('cell_stock','Source cells/nuclei in wash buffer already containing1percentBSA,1xDPBS,ProtectorRNaseinhibitor and1mMDTT;1000–3000perul,min35ul; sample unassigned',K.BIOLOGICAL,primary=True,chemical_state={'stock_BSA':'Source1percentBSA co-present in original cell/nucleus wash stock; no new addition','source_analytes':'Cells/nuclei carry RNA and genomicDNA; whole construct unassembled'})
    x=opaque(g,'emulsion_rt','Source first Tapestri v2 encapsulation; lysis and targeted RT50C60min/80C10min',stock,(add(g,'rna_reverse_pool','Source14.6uM RNA reverse-primer pool, each5prime '+REVERSE_HANDLE+' plus target-specific3prime arm; whole pool/orderings unassembled',K.OLIGO),add(g,'ssiv','Source10U/ul SuperScriptIV',K.CATALYST,R.CATALYST),add(g,'rt_nt','Source0.5mM EACH dNTP',K.MOLECULAR,R.STRUCTURAL),add(g,'proteinase','Source0.04mg/ml ProteinaseK,added immediately before loading; concurrent proteolysis/RT compatibility not proven here',K.CATALYST,R.CATALYST),condition(g,'rt_buffer','Source DTT/buffer/RNaseinhibitor/1percentNP40 and encapsulation oil; excess oil removal part of compound boundary')))
    x=opaque(g,'second_emulsion','Source second encapsulation,293.8ul BarcodingMasterMixv2 and UVrelease of bead-bound barcoded forward primers; repeated description205/207 is one operation',x,(add(g,'barcoding_stock','Source barcoding master mix, DNA forward0.21uM/reverse2.1uM,RNA forward0.42uM and barcode beads; whole components/barcodes unassembled'),condition(g,'barcoding_oil','Source barcoding oil/excess oil removed')))
    x=opaque(g,'target_pcr','Source targeted PCR:10 cycles72/61C regime then12 cycles72/51C regime; same-emulsion primer/enzyme/nucleotide components retained',x)
    x=cleanup(g,'break_emulsion','Source15ul extraction agent then cleanup buffer/enzyme37C60min350rpm',x);x=cleanup(g,'ampure1','Source0.72xAMPureXP first cleanup',x);x=cleanup(g,'ampure2','Source0.72xAMPureXP second cleanup',x)
    probe=add(g,'biotin_capture_oligo','Source0.7uM RNA Biotin Oligo; complete modified ordering '+CAPTURE,K.OLIGO,R.STRUCTURAL,state=capture_oligo())
    x=opaque(g,'capture_hybridize','Source DNA/cDNA amplicon mixture plus biotin oligo99C5min→ice5min',x,(probe,),S.HYBRIDIZATION)
    bound=g.product('rna_amplicon_bound','Source bead-bound RNA-derived DNA amplicon enriched fraction; recovery/background unassigned',K.MOLECULAR);sup=g.product('dna_amplicon_supernatant','Source DNA amplicon enriched supernatant; background unassigned',K.MOLECULAR)
    g.react('magnetic_capture','Source M270streptavidin RT30min then magnetic separation; RNA and DNA refer to amplicon origins, not RNA versus DNA polymers',S.AFFINITY_CAPTURE,((x,R.TEMPLATE),add(g,'m270','Source M270streptavidin beads',K.SELECTION,R.SELECTION),condition(g,'capture_buffer','Source printed10MTris-HCl/1mMEDTA/2MNaCl; literal10M retained, concentration not repaired')),(bound,sup))
    sup=cleanup(g,'dna_cleanup','Source0.72xAMPureXP DNA fraction cleanup',sup)
    for label,part in [('dna',sup),('rna',bound)]:
        part=opaque(g,label+'_index_pcr','Source ENTIRE '+label+'fraction, LibraryMix and '+('MissionBio kit indexes'if label=='dna'else'custom RNA indices')+';15cycles98/62/72C',part,(add(g,label+'_index_stock','Source LibraryMix plus '+label+' index-primer preparation; complete components unassembled'),))
        lib=g.product(label+'_library','Source '+label+'-origin amplicon library',K.LIBRARY);g.react(label+'_load','Source library sequencing; no unreported pooling assumed',S.GENERIC,((part,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='Associated manuscript NovaSeq6000 paired150bp. Complete kit/custom index primers, barcode identities/read landings and final published-version equivalence unassembled.'))
    g.add_source_item(SourceItem('unassembled','Actual linked supplementary workbook gated by PMC/TDM; full final article ordinary unverified. Complete source-linked vendor v2 components unassembled.','No absence/full-source review.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('Two source-origin library endpoints',len(g.sequencing),2)
    c('Second encapsulation recorded once',sum(s.id=='second_emulsion'for s in g.steps.values()),1)
    c('Explicit RT nucleotide contribution',any(p.material=='rt_nt'and p.role==R.STRUCTURAL for p in g.steps['emulsion_rt'].inputs),True)
    c('Post-PCR capture operates on mixed DNA products',list(g.steps).index('target_pcr')<list(g.steps).index('magnetic_capture'),True)
    c('Capture pair candidate remains inferred',capture_pair().strands['rna_derived_amplicon'].segs[0].inferred,True)
    p=sources_root()/Path(__file__).resolve().parents[1].name/'2159-8290.CD-24-0853_PMC12874418.xml.txt'
    if have(c,Source(p,'https://doi.org/10.1158/2159-8290.CD-24-0853')):
        t=p.read_text()
        for name,seq in [('forward handle',FORWARD_HANDLE),('reverse handle',REVERSE_HANDLE),('whole modified biotin probe',CAPTURE)]:c('Original '+name,seq in t,True)
    c('Private graph remains incomplete',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
