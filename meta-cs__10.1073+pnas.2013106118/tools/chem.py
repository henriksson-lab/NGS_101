"""Published selected preparation; private partial molecular model, not a complete audit."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from chemdraw import Scene,Segment
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from circular import CircularizationGeometry

TITLE='META-CS selected single-cell preparation'
AUTHORITY='Published2021 PMC7923680 main71/173–195; DatasetS1wholeoligotable; SIAppendixFigS1–S4captions. Three originals closed; no predecessor orDatasetS2 reading.'
INITIAL='One source mouth-pipetted/FACS-sorted single cell; source cell identity/ploidy unassigned'
READOUT='Source long fractionHiSeq2500PE250 and middleHiSeq4000PE150; shortPE150 onlynamed340-series samples. Source20percentPhiX added to each loading pool. Computationalstrandagreement is not physical removal or guaranteed error correction.'
DEFERRED='DatasetS2 run allocation; source-specific full Tn5/genomic ends, exact historical NEB multiplex/index primers and read geometry; citedMETA/DipC/LIANTI components not blindly inherited.'

def inventory_state(label):return None

from nextera import ME
from illumina import TRUSEQ_READ1,TRUSEQ_READ2
from chemdraw import revcomp,feature
TAG_PREFIXES=('GGCACCGAAAA', 'CTCGGCGATAAA', 'GGTGGAGCATAA', 'CGAGCGCATTAA', 'AGCCCGGTTATA', 'TCGGCACCAATA', 'GCCTGTGGATTA', 'GCGACCCTTTTA', 'GCATGCGGTAAT', 'GCGTTGCCATAT', 'GGCCGCATTTAT', 'ACCGCCTCTATT', 'CCGTGCCAAAAT', 'TCTCCGGGAATT', 'CCGCGCTTATTT', 'CTGAGCTCGTTTT')
def oligos(kind):
    tags=tuple(t+ME for t in TAG_PREFIXES)
    if kind=='tag':return tags
    arm=TRUSEQ_READ1 if kind=='adp1' else TRUSEQ_READ2[2:]
    return tuple(arm+t[:-6] for t in tags)
def inventory_state(label):
    kind='tag' if 'equal16' in label else 'adp1' if label.startswith('Adp1') else 'adp2' if label.startswith('Adp2') else None
    if kind:
        sc=Scene()
        for i,seq in enumerate(oligos(kind)):
            prefix=TAG_PREFIXES[i];start=seq.index(prefix);parts=[]
            if start:parts.append(Segment('source primer arm',seq[:start]))
            parts.extend([Segment('META finite tag',prefix,feature=feature('meta_tag_'+str(i+1),'molecular_complex_barcode','fixed')),Segment('source mosaic-end interval',seq[start+len(prefix):])]);sc.strand(kind+str(i+1),parts)
        return sc
    if label.startswith('Carrier'):
        sc=Scene();sc.strand('carrier',[Segment('carrier DNA','TCAGGTTTTCCTGAA')]);return sc
    return None

def unknown_scene():
    sc=Scene();sc.strand('unassembled',[Segment('geometry','[source geometry unassembled]',placeholder=True)]);return sc

def protocol_graph():
    g=ProtocolGraph(TITLE,authority=AUTHORITY,representation='published_instance')
    x=g.input('input',INITIAL,K.BIOLOGICAL,primary=True)
    def op(key,label,kind,reagent,product,oligos=()):
        nonlocal x
        enzyme=kind in {S.DIGESTION,S.END_REPAIR,S.EXTENSION,S.TAILING,S.CIRCULARIZATION,S.LIGATION,S.PCR,S.DEPROTEINIZATION,S.TAGMENTATION,S.CHEMICAL_CONVERSION}
        c=g.input(key+'_stock',reagent,K.CATALYST if enzyme else K.CONDITION)
        inputs=[(x,R.TEMPLATE),(c,R.CATALYST if enzyme else R.CONDITION)]
        for i,a in enumerate(oligos):
            o=g.input(key+'_oligo_'+str(i),a,K.OLIGO,state=inventory_state(a));inputs.append((o,R.STRUCTURAL))
        kw={}
        if kind==S.MAGIC:kw['reason']='Actual published component; full source-specific molecular transition unassembled.'
        if kind==S.CIRCULARIZATION:kw['geometry']=CircularizationGeometry.unavailable('Source intramolecular closure known; full substrate and closure geometry unassembled, not a claim of source absence.')
        if kind==S.LIGATION:kw['geometry']=AdapterLigationGeometry('unavailable','unavailable',(unknown_scene(),),'Source ligation operation known; full ends and sealed junction unassembled.')
        if kind==S.PCR:kw['geometry']=PCRGeometry.unavailable(unknown_scene(),reason='Source primer inventories/whole amplicon geometry not assembled in this selected model.')
        y=g.product(key+'_product',product,K.MOLECULAR);g.react(key,label,kind,inputs,(y,),**kw);x=y
    def clean(key,label):
        nonlocal x
        c=g.input(key+'_condition',label,K.CONDITION);y=g.product(key+'_kept','Physical source recovered preparation; carryover and recovery unquantified',K.MOLECULAR);w=g.product(key+'_waste','Actual source physically separated fractions',K.WASTE)
        g.react(key,label,S.CLEANUP,((x,R.TEMPLATE),(c,R.CONDITION)),(y,w));g.discard(w,'Physical workup, not a perfect molecular specificity assertion.');x=y
    def qc(key,label):
        nonlocal x
        y=g.product(key+'_retained','Production preparation retained after independent analytical allocation',K.MOLECULAR);q=g.product(key+'_assay',label,K.QC);g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(y,q));g.mark_terminal(Terminal(q.id,Disposition.QC,'Analytical aliquot; not returned to production.'));x=y
    op('lysis','2ul50C1h/70C15min',S.DEPROTEINIZATION,'METAlysisbuffer20mMTris/20mMNaCl/0.15percentTriton/25mMDTT/1mMEDTA/Qiagenprotease1.5mgml','Source cell lysate with genomic DNA/proteolysis products/reagents; optional-80C storage not imposed',('Carrier ssDNA500nM TCAGGTTTTCCTGAA',))
    op('tag','Add8ul;55C10min',S.TAGMENTATION,'META16-tagTn5transposome0.38nM/10mMTAPS/5mMMgCl2/8percentPEG8000','Tn5-tagged genomic fragments plus unprocessed DNA; randomtagpair combinations, off-products and exact junctions unassembled',('Source equal16 META tags annealed5uM each thenpooled;Tn5equimolar30minRT/125nMstock-80C',))
    op('strip','Add2ulstopbuffer50C30min/70C15min',S.DEPROTEINIZATION,'Source printed300nMNaCl/45mMEDTA/0.01percentTriton/Qiagenprotease1mgml','Proteolyzed source transposition preparation; cleavage products/remnants retained')
    op('adp1','Add13ulQ5mix;72C3min/98C30s/62C5min/72C1min',S.EXTENSION,'Q5polymerase0.25ul/buffer/GCenhancer/Mg/dNTP/BSA','Source first tagged-strand extension products plus templates/unused primers',('Adp1 source16primer pool0.8ul50uMtotal',))
    op('exo1','Add1ulExoI37C30min/80C20min',S.DIGESTION,'ExonucleaseI','Source first-extension preparation with cleaved primer fragments; no physical wash or perfect ssDNA specificity asserted')
    op('adp2','Add4ulQ5mix;98C30s/62C5min/72C1min',S.EXTENSION,'Q5polymerase0.05ul/buffer/GCenhancer/dNTP','Source second tagged-strand extension products plus earlier species/reagents',('Adp2 source16primer pool1ul50uMtotal',))
    op('exo2','Add1ulExoI37C30min/80C20min',S.DIGESTION,'ExonucleaseI','Source doublytagged strand preparations plus primer fragments/reagents')
    op('pcr','Add19ulQ5mix;98C20s;10cycles diploid/11haploid98C10s72C2min;72C2minfinal',S.PCR,'Q5polymerase0.2ul/buffer/GCenhancer/dNTP','Source amplified libraries; no fragment recovery/duplex fidelity guarantee',('NEBMultiplexUniversalPrimer5ul exacthistoricalidentityunassembled','NEBIndexPrimer5ulE7335S/E7500S/E7710S/E7730S exactidentityunassembled'))
    clean('column','DCC5columncleanup')
    qc('bioanalyzer','SourceBioanalyzer amplificationefficiency QC')
    other=g.input('other_cells','Other independently prepared single-cell libraries; identities/masses unassigned',K.LIBRARY,primary=True);pool=g.product('cell_pool','Physicalsource single-cell library pool',K.MOLECULAR);g.react('pool_cells','Sourcepool BEFORE sizefractionation',S.POOL,((x,R.TEMPLATE),(other,R.TEMPLATE)),(pool,));x=pool
    c=g.input('beads','SourceAMPure0.6x→further0.15x(total0.75x)→further0.15x(total0.9x)',K.SELECTION)
    long=g.product('long','Source first0.6x bead-retained long fraction',K.LIBRARY);mid=g.product('middle','Source further0.15x middle fraction from firstsupernatant',K.LIBRARY);short=g.product('short','Source finalfurther0.15x short fraction',K.LIBRARY);rest=g.product('final_supernatant','Physical unassigned residual after lastselection',K.MOLECULAR)
    g.react('fractionate','Sequentialbeadselection; actualfractions, notthreeequalaliquots',S.PARTITION,((x,R.TEMPLATE),(c,R.SELECTION)),(long,mid,short,rest));g.mark_terminal(Terminal(rest.id,Disposition.UNRESOLVED,'Lastsupernatant sourcefate unassigned; no idealized completeDNArecovery.'));g.mark_terminal(Terminal(short.id,Disposition.UNRESOLVED,'Shortfraction sequencedonlysource340-2/3/4/5/13/14/15; selectedcell/run unassigned.'))
    for key,part,plan in [('long',long,'HiSeq2500PE250'),('middle',mid,'HiSeq4000PE150')]:
        phix=g.input(key+'_phix','Source20percentPhiX loadingaddition; absolutevolumeunassigned',K.LIBRARY,primary=True);out=g.product(key+'_load','Physicalsource library/PhiXpool',K.LIBRARY);g.react(key+'_pool','SourcePhiXaddition',S.POOL,((part,R.TEMPLATE),(phix,R.TEMPLATE)),(out,));g.sequence(out,SequencingPlan(limitation=plan+';exact historicalprimer/index/rungeometryunassembled'))
    g.add_source_item(SourceItem('pending',DEFERRED,'Selected original sections read; linked-source/whole molecular geometry review incomplete. No accepted absence claim.'))
    return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('ExoI followsfirsttagging',g.steps['exo1'].inputs[0].material,'adp1_product')
    c('Secondtagging retainsexonuclease reactionmixture',g.steps['adp2'].inputs[0].material,'exo1_product')
    c('PCR followssecondExoI',g.steps['pcr'].inputs[0].material,'exo2_product')
    c('Pooling precedessizefractionation',g.steps['fractionate'].inputs[0].material,'cell_pool')
    c('Unassignedshortfractionretained',g.terminals['short'].disposition,Disposition.UNRESOLVED)
    c('All48sourceoligosassembled',sum(len(oligos(k)) for k in ('tag','adp1','adp2')),48)
    c('Private partial model cannot be publicly promoted',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
