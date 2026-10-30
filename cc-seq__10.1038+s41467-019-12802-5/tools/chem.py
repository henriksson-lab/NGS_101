"""Private source-owned CC-seq selected human post-extraction component."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,revcomp
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,SourceItem,SequencingPlan,Terminal,Disposition
from adapter_ligation import AdapterLigationGeometry
class B:
    def __init__(self,title,authority,start):
        self.g=ProtocolGraph(title,authority=authority,representation='published_instance');self.x=self.g.input('source',start,K.MOLECULAR,primary=True)
    def stock(self,key,label,kind=K.CATALYST,role=R.CATALYST,**kw):return(self.g.input(key,label,kind,**kw),role)
    def cond(self,key,label):return self.stock(key,label,K.CONDITION,R.CONDITION)
    def op(self,key,label,kind,inputs=(),geometry=None,product='Source copresent reaction mixture; incomplete/unreacted classes and carryover unassigned'):
        y=self.g.product(key+'_product',product,K.MOLECULAR);kw={}
        if kind==S.MAGIC:kw['reason']='Positive source operation; historical component formulation/whole molecular products unassembled, not source absence.'
        if kind==S.GENERIC:kw['changes_molecule']=False
        if geometry is not None:kw['geometry']=geometry
        self.g.react(key,label,kind,((self.x,R.TEMPLATE),)+tuple(inputs),(y,),**kw);self.x=y;return y
    def split(self,key,label,branches):
        yy=tuple(self.g.product(key+'_'+name,txt,K.MOLECULAR)for name,txt,disp in branches);self.g.react(key,label,S.SPLIT,((self.x,R.TEMPLATE),),yy)
        for y,(_,txt,disp)in zip(yy,branches):
            if disp is not None:self.g.mark_terminal(Terminal(y.id,disp,txt))
        self.x=yy[0];return yy
    def clean(self,key,label):
        y=self.g.product(key+'_kept','Source retained physical fraction; recovery/carryover unassigned',K.MOLECULAR);w=self.g.product(key+'_waste','Actual separated wash/workup fraction',K.WASTE);self.g.react(key,label,S.CLEANUP,((self.x,R.TEMPLATE),self.cond(key+'_condition',label)),(y,w));self.g.discard(w,'Source separated workup; no perfect specificity/yield.');self.x=y
    def capture(self,key,label,inputs):
        y=self.g.product(key+'_bound','Source bead/membrane bound mixture; purity unassigned',K.MOLECULAR);u=self.g.product(key+'_unbound','Unretained physical capture fraction',K.MOLECULAR);self.g.react(key,label,S.AFFINITY_CAPTURE,((self.x,R.TEMPLATE),)+tuple(inputs),(y,u));self.g.mark_terminal(Terminal(u.id,Disposition.UNRESOLVED,'No later recovery assigned.'));self.x=y
    def finish(self,limitation):
        self.split('loading_split','Source library aliquot for loading/pooling; amounts unassembled',[('aliquot','Selected library loading aliquot',None),('remainder','Untransferred library stock',Disposition.UNRESOLVED)])
        other=self.g.input('other_library_aliquots','Other source multiplexed library aliquots; identities/amounts unassembled',K.MOLECULAR);lib=self.g.product('library','Actual source library pool; whole products unassembled',K.LIBRARY);self.g.react('pool','Source multiplexed sequencing pool',S.POOL,((self.x,R.TEMPLATE),(other,R.TEMPLATE)),(lib,));self.g.sequence(lib,SequencingPlan(limitation=limitation));self.g.add_source_item(SourceItem('pending','Whole molecular products/history-specific kit internals/read landing unassembled.','Selected bounded pass only; no confirmed absence/full chemistry proof.'));return self.g

from illumina import TRUSEQ_READ1,NEBNEXT_ARM_READ2
P5TOP=TRUSEQ_READ1
P5BOTTOM=revcomp(P5TOP[:-1])+'T'
P7TOP=NEBNEXT_ARM_READ2+'ACT'
P7BOTTOM=revcomp(P7TOP[:-1])+'T'
SOURCE_ORDERINGS={'p5top':P5TOP,'p5bottom':'5P-'+P5BOTTOM[:-1]+'-3invdT','p7top':'5P-'+P7TOP[:-1]+'-3invdT','p7bottom':P7BOTTOM}
SOURCE_OLIGO_EXCEPTIONS={'p5':'PublishedSI actualTable5p17 originalrender read: P5upper ACACTCTTTCCCTACACGACGCTCTTCCGATCT; lower drawn3to5 TiTGTGAGAAAGGGATGTGCTGCGAGAAGGCTAGp, actual5to3 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT-invdT. Parser cannot match additional inverted terminalT in reverse displayed orientation. Local wholeoriginalordering+modification source fixtures retain every fixedbase andTi, with32bpcomputedpair and3Tcohesiveupperend; no hidden mismatch repaired.'}
def adapter(n):
    s=Scene()
    if n==5:
        s.strand('top',[Segment('core',P5TOP[:-1],note='Source firstAterminalnucleoside3OH; linkagechemistryunassembled'),Segment('cohesive_T','T')],mod5='A nucleoside3OH',mod3='OH')
        s.anneal('bottom',[Segment('core_rc',P5BOTTOM[:-1]),Segment('inverted_dT','T',note='Source inverted dT with3prime–3prime bond, not conventional terminalT')],to='top',pair=('core_rc','core'),unpaired=('inverted_dT',),mod5='P',mod3='inverted dT3prime–3prime')
    else:
        s.strand('top',[Segment('core',P7TOP[:-1]),Segment('inverted_dT','T',note='Source inverted dT with3prime–3prime bond')],mod5='P',mod3='inverted dT3prime–3prime')
        s.anneal('bottom',[Segment('core_rc',P7BOTTOM[:-1]),Segment('cohesive_T','T')],to='top',pair=('core_rc','core'),unpaired=('cohesive_T',),mod5='G nucleoside3OH',mod3='OH')
    return s
def end(n):
    s=adapter(n)
    if n==7:
        s.anneal('substrate',[Segment('insert','[sourcefragment]',placeholder=True),Segment('dA','A')],to='bottom',pair=('dA','cohesive_T'),unpaired=('insert',),mod3='OH')
        s.ligation_end('bottom','cohesive_T',kind='cohesive',polarity='3-prime',owner='adapter',text='Source P7 3T/adjacent5P ligated end')
        s.ligation_end('substrate','dA',kind='cohesive',polarity='3-prime',owner='substrate',text='Source dA-tailed fragment; other ends/junctions unassembled')
    else:
        s.anneal('substrate',[Segment('insert','[sourcefragment]',placeholder=True),Segment('dA','A')],to='top',pair=('dA','cohesive_T'),unpaired=('insert',),mod3='OH')
        s.ligation_end('top','cohesive_T',kind='cohesive',polarity='3-prime',owner='adapter',text='Table5P5source3Tcohesiveend/adjacentbottom5P')
        s.ligation_end('substrate','dA',kind='cohesive',polarity='3-prime',owner='substrate',text='SourceinitialdA-tailretainedthroughTDP2; completeproduct/unconvertedclassesunassembled')
    return s
def protocol_graph():
    b=B('CC-seq selected human component','Publishedmain205–225/261 plusSI Table5originalp17;2originals. Initial nonproteolysed extraction described in notes, outside this selected component.','Source pooled1ml nonproteolysed/RNase-treated humanDNA; individual aliquot allocations/leftovers unassembled');g=b.g
    b.op('sonicate','Source Covaris300–400bp 15min',S.FRAGMENTATION,(b.cond('covaris','10percentduty75W200cycles/burst'),))
    b.op('binding_buffer','Source1mlTE dilution plus.3MNaCl/.2percentTriton/.1percentsarcosyl',S.GENERIC,(b.cond('binding_conditions','Source protein-binding buffer'),))
    b.split('column_aliquot','Source divided acrosssilica columns approximately20ugDNAeach',[('selected','Selected approximately20ug column input',None),('others','Other source column input aliquots for actualparallelcapture/elution',None)])
    b.capture('silica','Source QIAGEN silica protein-binding capture, flowthrough reapplied', (b.stock('silica_membrane','Source silica-fibre membrane',K.SELECTION,R.SELECTION),))
    b.clean('wash','Source6x600ulTEN1min, retainmembrane')
    y=g.product('tes_eluate','Source100ulTES eluate',K.MOLECULAR);r=g.product('silica_residual','Membrane/residual fraction',K.MOLECULAR);g.react('elute','Source100ulTES .5percentSDS releases CCs',S.CLEANUP,((b.x,R.TEMPLATE),b.cond('tes','10mMTris/1mMEDTA/.5percentSDS')),(y,r));g.mark_terminal(Terminal(r.id,Disposition.UNRESOLVED,'Residualmembrane'));b.x=y
    others=g.product('other_column_eluates','Other actualparallelcolumns retained eluates for500ulpool; allocationsunassembled',K.MOLECULAR);other_residual=g.product('other_column_residual','Othercolumns separatedunretained/wash/membranefractions; recoveryunassigned',K.MOLECULAR);g.react('parallel_columns','Otheractualcolumns highsaltcapture/flowthroughreapply/6TENwashes/100ulTESelution; detailedpercolumnallocationunassembled',S.MAGIC,((g.materials['column_aliquot_others'],R.TEMPLATE),b.stock('parallel_column_membranes','Sourceotheractualsilicacolumnmembranes',K.SELECTION,R.SELECTION),b.cond('parallel_column_conditions','Same sourcebinding/wash/TESrecipe, quantitiesperpoolunassembled')),(others,other_residual),reason='Sourceparallel-columnworkuppositive; exactallocation/internalintermediategeometryunassembled.');g.mark_terminal(Terminal(other_residual.id,Disposition.UNRESOLVED,'Physicalothercolumnworkupfractions.'));pool=g.product('pooled_eluates','Source500ul pooledTES eluates',K.MOLECULAR);g.react('elution_pool','Source pool columneluates to500ul',S.POOL,((b.x,R.TEMPLATE),(others,R.TEMPLATE)),(pool,));b.x=pool
    b.op('proteinase','Source1mg/mlPK30min60C',S.DEPROTEINIZATION,(b.stock('pk','SourceProteinaseK'),))
    b.op('carrier','Source.2mg/mlglycogen/200mMNaOAc/1.41mlethanol−80Covernight',S.GENERIC,(b.stock('glycogen','Sourcephysicalglycogencarrier',K.MOLECULAR,R.STRUCTURAL),b.cond('precip_conditions','Sourceethanol/acetate/temperature')))
    b.clean('precipitate','Source20,000g1h4C/1.5ml70percentethanol/airdry10min/resuspend52ulTris')
    b.split('qubit1','Source52ul=50ulproduction+2ulQubit',[('production','Source50ulproduction',None),('qc','Source2ulQubit',Disposition.QC)])
    # A single opaque source kitpreparation contributes molecular and catalytic components; no extra fresh doses inferred.
    st=g.input('kit_stock','HistoricalUltraII endprep/ligation preparation; actualcomposition/aliquotsunassembled',K.MOLECULAR);cat=g.product('kit_endprep_enzyme','Endprep catalyst component of source preparedreaction',K.CATALYST);nt=g.product('kit_nucleotides','Source endprep nucleotide component includingdATP; actualamountsunassembled',K.MOLECULAR);lg=g.product('kit_ligase','Ligase component allocated fromsource preparation',K.CATALYST);g.react('kit_components','Represent components ofsource preparation, not separateinventedstock additions',S.GENERIC,((st,R.TEMPLATE),),(cat,nt,lg),changes_molecule=False)
    b.op('endprep','SourceUltraII endrepair/dA-tail; source-specific formulationunassembled',S.END_REPAIR,((cat,R.CATALYST),(nt,R.STRUCTURAL)))
    p7=g.input('p7','SourceTable5preparedcustomP7 selectedreactiondose; stockpreparationoutsideboundary',K.ADAPTER,state=adapter(7),chemical_state={'source_ordering':SOURCE_ORDERINGS['p7top']+' ; '+SOURCE_ORDERINGS['p7bottom'],'nonstandard_ends':'invdT3prime–3prime andGterminalnucleoside3OH'})
    b.op('ligate_p7','SourcecustomP7 ligation BEFORETDP2',S.LIGATION,((p7,R.STRUCTURAL),(lg,R.CATALYST)),AdapterLigationGeometry.known('cohesive',end(7),statement='Table5wholemodifiedadapter/source dA-tail knownlocalA:T pair; whole sealedproductsunassembled.'))
    b.clean('p7_clean','SourceAMPure78:90/elute50ulTris')
    b.op('tdp2','Source3ul10uMhumanTDP2+50ulreactionbuffer/50ulDNA1h37C',S.CHEMICAL_CONVERSION,(b.stock('tdp2_enzyme','Source recombinant humanTDP2'),b.stock('bsa','SourcephysicalBSAinTDP2buffer',K.MOLECULAR,R.STRUCTURAL),b.cond('tdp2_buffer','Source100mMTrisOAc/100mMNaOAc/2mMMgOAc/2mMDTT/200ug/mlBSA addedbuffer; final dilution sourcequalified')))
    b.clean('tdp2_clean','SourceAMPure103:103/elute52ulTris')
    p5=g.input('p5','SourceTable5customP5 selecteddose; wholekitdoseunassembled',K.ADAPTER,state=adapter(5),chemical_state={'source_ordering':SOURCE_ORDERINGS['p5top']+' ; '+SOURCE_ORDERINGS['p5bottom'],'nonstandard_ends':'invdT3prime–3prime andAterminalnucleoside3OH'})
    b.op('ligate_p5','SourceP5ligation WITHOUTsecondendrepair',S.LIGATION,((p5,R.STRUCTURAL),b.stock('p5_ligase','Source secondligation preparation catalyst; exacthistoricalstockunassembled')),AdapterLigationGeometry.known('cohesive',end(5),statement='Table5P5source3Tcohesiveend; initialdA-tail, completeTDP2substratepopulation/wholejunctionunassembled.'))
    b.clean('p5_clean','SourceAMPure78:90/elute17ul')
    b.split('qubit2','Source17ul=15ulPCR+2ulQubit',[('production','Source15ulPCRtemplate',None),('qc','Source2ulQubit',Disposition.QC)])
    b.op('pcr','SourceUltraII PCR/universalP5/indexedP7; exactoligos/cyclesunassembled',S.MAGIC,(b.stock('pcr_stock','ActualhistoricalPCRpreparation molecularcontributor',K.MOLECULAR,R.STRUCTURAL),))
    b.clean('pcr_clean','Source50ulTrisdilution/AMPure84:100/elute30ul1mMTris8.1')
    b.clean('size_select','SourceBluePippin200–600bpactualphysicalselection')
    b.split('bioanalyzer','SourceBioanalyzermolarityQC aliquot/remainingstock; volumesunassigned',[('production','Remainingproductionlibrary',None),('qc','ActualBioanalyzerQCaliquot',Disposition.QC)])
    return b.finish('PublishedsourceMiSeqv3PE75 orNextSeq500v2PE42; wholePCRprimers/index/readlandingunassembled. Narrowpartialcomponent, notcompleteCCvalidation.')
def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    import re
    c=Check();g=protocol_graph();c('P7 before unmask',g.steps['tdp2'].inputs[0].material,'p7_clean_kept');c('No secondendprep',g.steps['ligate_p5'].inputs[0].material,'tdp2_clean_kept');c('Qubitactualaliquot',g.terminals['qubit1_qc'].disposition,Disposition.QC);c('Qubitnotwholeproduction',g.steps['endprep'].inputs[0].material,'qubit1_production');c('SecondQubitnotPCR',g.steps['pcr'].inputs[0].material,'qubit2_production');c('Silicaunretainedopen',g.terminals['silica_unbound'].disposition,Disposition.UNRESOLVED);c('Actualcolumnpool',g.steps['elution_pool'].kind,S.POOL);c('OthercolumneluatesderivefromsameinitialDNA',g.steps['parallel_columns'].inputs[0].material,'column_aliquot_others')
    c('P5wholeupper',P5TOP,'ACACTCTTTCCCTACACGACGCTCTTCCGATCT');c('P5wholelowerinclinvT',P5BOTTOM,'GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTT');c('P7wholeupperinclinvT',P7TOP,'GATCGGAAGAGCACACGTCTGAACTCCAGTCACT');c('P7wholelower',P7BOTTOM,'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')
    c('Invertedbaseinventoryretained',sum(len(z.top)for z in adapter(7).strands['top'].segs),len(P7TOP));c('Privatepartial',bool(g.publication_problems()),True)
    path=sources_root()/Path(__file__).resolve().parents[1].name/'s41467-019-12802-5_supp_41467_2019_12802_MOESM1_ESM.pdf.txt'
    if have(c,Source(path,'https://doi.org/10.1038/s41467-019-12802-5')):
        text=re.sub(r'\s+','',path.read_text())
        c('SourceP5wholeupper',P5TOP in text,True)
        c('SourceP5loweroriginalorientationmods',('Ti'+P5BOTTOM[:-1][::-1]+'p') in text,True)
        c('SourceP7upperwholemods',('p'+P7TOP[:-1]+'Ti') in text,True)
        c('SourceP7loweroriginalorientation',P7BOTTOM[::-1] in text,True)
    c.report()
if __name__=='__main__':selftest()
