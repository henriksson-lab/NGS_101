"""Private source-qualified m7G-seq pulldown branch; associated manuscript."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from illumina import TRUSEQ_READ1,TRUSEQ_READ2
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
TITLE='m7G-seq: Zhang2019 selected enriched mRNA branch'
AUTHORITY='Associated AM PMC6588483.1/NIHMS1526211 main211–215, Fig3B97; proteinase reference207. Final published body and cited MethodS1 not checked for complete equivalence.'
RNA_BODY=revcomp(TRUSEQ_READ1)[:19]
DNA_BODY=revcomp(TRUSEQ_READ2)[:22]
RT_PRIMER=TRUSEQ_READ1[13:]
SOURCE_ORDERINGS={'rna_adapter':'5rApp-'+RNA_BODY+'-3SpC3','dna_adapter':'5Phos-'+'N'*10+DNA_BODY+'-3SpC3','rt_primer':RT_PRIMER}
def adapter(which):
    s=Scene();parts=[]
    if which=='dna':parts.append(Segment('source_random10','N'*10,placeholder=True,note='Source N10 degeneracy; purpose and read landing unassigned'))
    parts.append(Segment('source_body',RNA_BODY if which=='rna'else DNA_BODY));s.strand('adapter',parts,mod5='rApp'if which=='rna'else'P',mod3='SpC3');return s
def ends(which):
    s=adapter(which);s.strand('substrate',[Segment('substrate','[RNA fragment]'if which=='rna'else'[cDNA]',placeholder=True)],mod3='OH')
    s.ligation_end('adapter','source_random10'if which=='dna'else'source_body',through='source_body',kind='single-strand',polarity='5-prime',owner='adapter',text='Source5rApp'if which=='rna'else'Source5P; blocked3C3')
    s.ligation_end('substrate','substrate',kind='single-strand',polarity='3-prime',owner='substrate',text='RNA after PNK/rSAP or nascent cDNA; full junction unassembled');return s

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance');x=g.input('mrna','Source4–6ug human mRNA after two polyA purification rounds; cell line unassigned',K.MOLECULAR,primary=True)
    def inp(k,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(g.input(k,label,kind,**kw),role)
    def op(k,label,kind,ins=(),geometry=None,product='Source molecular mixture, residual substrates and reagents retained until workup'):
        nonlocal x
        y=g.product(k+'_product',product,K.MOLECULAR);kw={}
        if geometry is not None:kw['geometry']=geometry
        if kind==S.MAGIC:kw['reason']='Source operation positive; complete reagents, molecular geometry or reaction allocation unassembled, not confirmed absent.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        g.react(k,label,kind,((x,R.TEMPLATE),)+tuple(ins),(y,),**kw);x=y
    def clean(k,label):
        nonlocal x
        c=g.input(k+'_condition',label,K.CONDITION);y=g.product(k+'_kept','Source retained preparation; recovery/carryover unassigned',K.MOLECULAR);w=g.product(k+'_waste','Source physical workup fractions',K.WASTE);g.react(k,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Physical workup, not perfect substrate specificity.');x=y
    def nt(k,label):return inp(k,label,K.MOLECULAR,R.STRUCTURAL)
    def cond(k,label):return inp(k,label,K.CONDITION,R.CONDITION)
    op('fragment','SourceAM8740 fragmentation into50–100nt; exact reaction unassembled',S.MAGIC)
    op('decap','SourceTobaccoDecappingPlus2 standard procedure',S.MAGIC,(inp('decapping_stock','SourceTobaccoDecappingPlus2 preparation'),))
    clean('decap_clean','SourceOligoCleanConcentrator')
    op('pnk','Source5ulPNK+5ul10xbuffer50ul37C1h,3prime dephosphorylation',S.DEPHOSPHORYLATION,(inp('pnk_enzyme','SourceT4PNK EK0032'),cond('pnk_buffer','Source B0201S buffer; no ATP addition specified here')))
    op('rsap','Source5ulrSAP50ulCutSmart37C1.5h;5prime/3prime dephosphorylation',S.DEPHOSPHORYLATION,(inp('rsap_enzyme','Source rSAP'),))
    a=g.product('rna20','Source20ul repaired/dephosphorylated RNA',K.MOLECULAR);r=g.product('rna_remainder','Untransferred repaired RNA stock; amount/fate unassigned',K.MOLECULAR);g.react('rna_aliquot','Source20ul selected from preceding50ul preparation',S.SPLIT,((x,R.TEMPLATE),),(a,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Untransferred stock.'));x=a
    ad=g.input('rna_adapter','Source1.6ul100uM3linker '+SOURCE_ORDERINGS['rna_adapter'],K.ADAPTER,state=adapter('rna'),chemical_state={'source_ordering':SOURCE_ORDERINGS['rna_adapter'],'five_prime':'rApp','three_prime':'C3spacer'})
    op('rna_anneal','Source70C2min/immediateice withRNAadapter',S.HYBRIDIZATION,((ad,R.STRUCTURAL),),product='Same-tube RNA and adapter component before ligase addition')
    t=g.product('rna_template_component','Same-tube RNA component',K.MOLECULAR);p=g.product('retained_rna_adapter','Same-tube retained adapter, not a new stock addition',K.ADAPTER,state=adapter('rna'));g.react('rna_components','Same-tube components, no physical fractionation',S.GENERIC,((x,R.TEMPLATE),),(t,p),changes_molecule=False);x=t
    op('rna_ligate','Source4ulRnl2truncatedKQ/5ul10xbuffer/15ul50percentPEG,50ul25C2h16C8hormore',S.LIGATION,((p,R.STRUCTURAL),inp('rnl2','Source T4RNA ligase2 truncatedKQ')),AdapterLigationGeometry.known('single-strand',ends('rna'),statement='Source-known adapter5rApp/RNA3OH interface; sealed whole product unassembled.'))
    op('deadenylate','Source dilute97ul+3ul5Deadenylase30C1h',S.CHEMICAL_CONVERSION,(inp('deadenylase','Source5Deadenylase'),))
    op('recjf','Source1.5ulRecJf37C1h',S.DIGESTION,(inp('recjf_enzyme','SourceRecJf'),),product='Source ligated RNA plus residual adapter/substrate classes; no perfect depletion claim')
    clean('ligated_clean','SourceRNAcleanconcentrator')
    op('reduce','Source1MNaBH4 stock as5x,25C45–60min; conversion incomplete',S.CHEMICAL_CONVERSION,(nt('nabh4','Source freshly aqueous NaBH4 reducing reagent; final convention from5xstock'),cond('reduction_condition','Source25C45–60minaqueous reaction'),))
    clean('reduction_workup','Source7volRNA binding buffer to quench;RNAcolumn recovery')
    op('biotinylate','Source35ulRNA+5ul1MMESpH4.5+10ul50mMhydrazidebiotininDMSO55C1h',S.CHEMICAL_CONVERSION,(nt('hydrazide_biotin','SourceEZLinkHydrazideBiotin, incorporated tag contribution'),cond('mes','SourceMES/DMSO')),product='Source biotinylated abasic RNA and unconverted/reduced/unlabelled classes; no stoichiometric conversion asserted')
    clean('biotin_clean','SourcequickRNApurification')
    op('capture','Source20ulMyOneC1resinwashed2xB&W/resuspend2xB&W equal sample volume4C20min',S.BINDING,(inp('streptavidin','SourceMyOneStreptavidinC1',K.SELECTION,R.SELECTION),))
    clean('capture_wash','Source5washes50mMTrispH7.4/300mMKCl/.05percentNP40;retain beads')
    op('pk','SourceproteinaseK55C30min, referenced buffer207; one30min selected here, no extra second digestion inferred',S.DEPROTEINIZATION,(inp('pk_enzyme','SourceProteinaseK'),))
    clean('elute','Source save flowthrough/RNAcleanconcentrator; residual bead fraction physical')
    primer=g.input('rt_primer','Source1ul2uM '+RT_PRIMER,K.OLIGO,chemical_state={'source_ordering':RT_PRIMER})
    op('rt_anneal','Source65C2min/immediateice',S.HYBRIDIZATION,((primer,R.STRUCTURAL),))
    t=g.product('rt_template','Same-tube RNA component',K.MOLECULAR);p=g.product('retained_rt_primer','Same-tube retained RT primer',K.OLIGO);g.react('rt_components','Same-tube components, not fractions/new stock',S.GENERIC,((x,R.TEMPLATE),),(t,p),changes_molecule=False);x=t
    op('rt','Source1–2ulHIVRT+2ul10mMdNTP/AMVbuffer/RNaseOUT37C1h75C5min',S.REVERSE_TRANSCRIPTION,((p,R.STRUCTURAL),inp('hiv_rt','Source Recombinant HIV RT Worthington'),nt('rt_nt','Source2ul10mMdNTP mix; perbase/total convention unassigned')))
    clean('cdna_clean','SourceOligoCleanConcentrator')
    ad=g.input('dna_adapter','Source1ul50uM '+SOURCE_ORDERINGS['dna_adapter'],K.ADAPTER,state=adapter('dna'),chemical_state={'source_ordering':SOURCE_ORDERINGS['dna_adapter']})
    op('dna_anneal','Source75C2min withcDNAlinker',S.HYBRIDIZATION,((ad,R.STRUCTURAL),))
    t=g.product('dna_template_component','Same-tube cDNA component',K.MOLECULAR);p=g.product('retained_dna_adapter','Same-tube retained cDNA adapter',K.ADAPTER,state=adapter('dna'));g.react('dna_components','Same-tube components not new stock',S.GENERIC,((x,R.TEMPLATE),),(t,p),changes_molecule=False);x=t
    op('dna_ligate','Source2ulbuffer/10ul50percentPEG/2ul10mMATP/1ulRnl125C12h',S.LIGATION,((p,R.STRUCTURAL),inp('rnl1','SourceT4RNA ligase1 high concentration'),cond('atp','SourceATP/PEG/buffer')),AdapterLigationGeometry.known('single-strand',ends('dna'),statement='Known cDNA3OH/adapter5P local interface; sealed whole junction unassembled.'))
    op('pcr','Sourceuniversal/indexedNEBNextMultiplexOligos PCR; exact whole primers/cycles unassembled',S.MAGIC,(nt('pcr_preparation','Source PCR preparation including named universal/indexed primer stocks; complete components/catalyst identity unassembled'),))
    lib=g.product('library','Source enriched m7G-seq library',K.LIBRARY);g.react('load','Source sequencerloading',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False);g.sequence(lib,SequencingPlan(limitation='SourceNextSeq500 singleend80bp. Whole primer landing/read masks unassembled. N10 literal degeneracy has unassigned purpose and read landing; no molecule-UMI meaning assigned.'))
    g.add_source_item(SourceItem('pending','Cited MethodS1 actual supplement1 PDF PMC gated; associated AM recipe not claimed entirely equivalent to final published body. Input/beforepulldown companion branches and complete sealed constructs unassembled.','No source absence claim.'));return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();g=protocol_graph();c('AdapterC3doesnotaddcytidine',''.join(x.top for x in adapter('rna').strands['adapter'].segs),RNA_BODY)
    c('RNAadapterretainedfromsameaddition',g.producer_of('retained_rna_adapter').id,'rna_components')
    c('RTretainedprimerrole',any(p.material=='retained_rt_primer'and p.role==R.STRUCTURAL for p in g.steps['rt'].inputs),True)
    c('Biotinisstructuralcontribution',any(p.material=='hydrazide_biotin'and p.role==R.STRUCTURAL for p in g.steps['biotinylate'].inputs),True)
    c('RNAremainderunresolved',g.terminals['rna_remainder'].disposition,Disposition.UNRESOLVED)
    path=sources_root()/Path(__file__).resolve().parents[1].name/'j.molcel.2019.03.036_PMC6588483.xml.txt'
    if have(c,Source(path,'https://doi.org/10.1016/j.molcel.2019.03.036')):
        text=path.read_text();c('OriginalRNAwholeordering',"5’rApp-"+RNA_BODY+'-3SpC3' in text,True);c('OriginalRTprimer',RT_PRIMER in text,True);c('OriginalcDNAwholeordering',"5’ Phos-"+'N'*10+DNA_BODY[:10]+' '+DNA_BODY[10:]+'-3SpC3' in text,True)
    c('Privateincompletegraph',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
