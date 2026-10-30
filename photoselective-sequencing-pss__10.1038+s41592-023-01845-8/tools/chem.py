"""Published PSS oligo inventory and local interfaces; no ordered protocol graph."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'lib'))
from chemdraw import Scene,Segment,feature
from nextera import ADAPTOR_S5,ADAPTOR_S7,ME_RC
from illumina import P5,P7
TITLE='Photoselective sequencing (PSS) source components'
SOURCE='PublishedSI MOESM1 originalp2 oligoTable, ExtendedDataFig1A; mainHTMLsubscriptionpreview. ActualcoreMethodsnotreceived; no wholepreparationgraph.'
SOURCE_ORDERINGS={
'PSS1':'/5Alex546N//iSpPC/TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG',
'PSS2':'/5Alex546N//iSpPC/GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG',
'PSS3':'/5phos/C*T*G*T*C*T*C*T*T*A*T*A*C*A*/3ddC/',
'PSS4':'AATGATACGGCGACCACCGAGATCTACACCTCTCTAT',
'PSS5':'CAAGCAGAAGACGGCATACGAGATBBBBBBBBCGACTGACTCT',
'PSS6':'TGCCGACGAATAGAGAGGTGTAG/3ddC/',
'PSS7':'CCCACGAGACAGAGTCAGTCG/3ddC/',
'PSS8':'AATGATACGGCGACCACCGAGAT',
'PSS9':'CAAGCAGAAGACGGCATACGAGAT'}

def seq(name):
    return {'PSS1':ADAPTOR_S5,'PSS2':ADAPTOR_S7,'PSS3':ME_RC[:15],'PSS4':P5+'CTCTCTAT','PSS5':P7+'B'*8+'CGACTGACTCT','PSS6':'TGCCGACGAATAGAGAGGTGTAGC','PSS7':'CCCACGAGACAGAGTCAGTCGC','PSS8':P5[:23],'PSS9':P7}[name]

def inventory(name):
    sc=Scene();parts=[] if name=='PSS5' else [Segment(name,seq(name))]
    if name=='PSS5':parts=[Segment('P7',P7),Segment('source_known_index_placeholder','B'*8,placeholder=True,feature=feature('pss_sample','sample_index','unknown'),note='Bmeansknownsamplebarcode,notIUPACnonA;exactindexnotselected'),Segment('source_tail','CGACTGACTCT')]
    sc.strand(name,parts,mod5='Alexa546/photocleavable linker'if name in {'PSS1','PSS2'}else'P'if name=='PSS3'else'',mod3='ddC'if name in {'PSS3','PSS6','PSS7'}else'')
    return sc

def source_blocked_adapter(name):
    """Inventory pairing only; no claim of actual loaded Tn5 stoichiometry."""
    if name not in {'PSS1','PSS2'}:raise ValueError('Only two source tagmentation oligos')
    sc=Scene();body=seq(name);sc.strand('long',[Segment('handle',body[:-19]),Segment('ME',body[-19:-15]),Segment('ME_paired',body[-15:])],mod5='Alexa546/photocleavable linker')
    sc.anneal('short',[Segment('source_short',seq('PSS3'),note='First14sourcepositionsPS;lastCis/3ddC/modifiednucleotide')],to='long',pair=('source_short','ME_paired'),mod5='P',mod3='ddC')
    return sc

def ligation_interface(side):
    """Source oligo complementarity only; no ordered preparation/product claim.

    Terminal ddC is an explicit extra base outside the two tested splint arms;
    its pairing/occupancy is not asserted against an upstream handle/barcode.
    """
    if side not in {'P5','P7'}:raise ValueError('P5 or P7')
    substrate,adapter,splint,n=('PSS1','PSS4','PSS6',9)if side=='P5'else('PSS2','PSS5','PSS7',10)
    a=seq(adapter);spl=seq(splint);tail=len(spl)-n-1
    sc=Scene();sc.strand('splint',[Segment('substrate_arm',spl[:n]),Segment('adapter_arm',spl[n:-1]),Segment('terminal_ddC','C',note='Source/3ddC/;outsideassertedpairedarms')],mod3='ddC')
    sc.anneal('substrate',[Segment('cleaved_handle',seq(substrate)[:n]),Segment('handle_remainder',seq(substrate)[n:]),Segment('DNA_remainder','[DNA context unassembled]',placeholder=True)],to='splint',pair=('cleaved_handle','substrate_arm'),mod5='P')
    parts=([Segment('P7',P7),Segment('sample_code','B'*8,placeholder=True,feature=feature('pss_sample','sample_index','unknown'))]if side=='P7'else[Segment('adapter_prefix',a[:-tail])])+[Segment('ligation_tail',a[-tail:])]
    sc.anneal('adapter',parts,to='splint',pair=('ligation_tail','adapter_arm'),unpaired=('P7','sample_code')if side=='P7'else('adapter_prefix',))
    if sc.strands['substrate'].span('cleaved_handle')[1]!=sc.strands['adapter'].span('ligation_tail')[0]:raise ValueError('Source complementary arms must meet at adjacent nick')
    return sc

def selftest():
    from checks import Check
    c=Check()
    c('WholePSS1','/5Alex546N//iSpPC/'+seq('PSS1'),SOURCE_ORDERINGS['PSS1'])
    c('WholePSS2','/5Alex546N//iSpPC/'+seq('PSS2'),SOURCE_ORDERINGS['PSS2'])
    c('WholemodifiedshortME','/5phos/'+'*'.join(seq('PSS3')[:-1])+'*/3ddC/',SOURCE_ORDERINGS['PSS3'])
    c('WholePSS4',seq('PSS4'),SOURCE_ORDERINGS['PSS4'])
    c('WholePSS5pattern',seq('PSS5'),SOURCE_ORDERINGS['PSS5'])
    c('WholePSS6',seq('PSS6')[:-1]+'/3ddC/',SOURCE_ORDERINGS['PSS6'])
    c('WholePSS7',seq('PSS7')[:-1]+'/3ddC/',SOURCE_ORDERINGS['PSS7'])
    c('WholePSS8',seq('PSS8'),SOURCE_ORDERINGS['PSS8'])
    c('WholePSS9',seq('PSS9'),SOURCE_ORDERINGS['PSS9'])
    for n in ('PSS1','PSS2'):c(n+'sourceblockedshortpairs',source_blocked_adapter(n).strands['short'].mod3,'ddC')
    for side in ('P5','P7'):
        sc=ligation_interface(side)
        c(side+'sourcecomplementaryarmsadjacent',sc.strands['substrate'].span('cleaved_handle')[1],sc.strands['adapter'].span('ligation_tail')[0])
        c(side+'ddCretainedoutsidepairedarms',(sc.strands['splint'].segs[-1].top,sc.strands['splint'].mod3),('C','ddC'))
    c('Barcodeplaceholderstructured',inventory('PSS5').strands['PSS5'].segs[1].feature.role,'sample_index');c.report()
if __name__=='__main__':selftest()
