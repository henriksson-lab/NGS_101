"""Published PRIM-seq FigS9 component only; full NGS preparation remains uncollected."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'lib'))
from chemdraw import Scene,Segment,revcomp
from adapter_ligation import AdapterLigationGeometry
from protocol_graph import ProtocolGraph,MaterialKind as K,InputRole as R,StepKind as S,Terminal,Disposition,SourceItem
BOTTOM_AS_PRINTED_3_TO_5='GGTTACCGCGGCCCGGAAAGAAATACAAAAACCGCAGAACCAGT'
BOTTOM_5_TO_3=BOTTOM_AS_PRINTED_3_TO_5[::-1]
TOP='TGACCAATGGCGCCGGGCCTTTCTTTATGTTTTTGGCGTCTTGG'
AUTHORITY='Published2025/issue2026PRIMmainpreviewidentity/abstract/ED1–2 only; publishedSI1Fig9actualPDFp13/printed12. Componentonly, notwholeSMART-display/REILISrecipe.'
def bottom_scene(activated=True):
    sc=Scene();sc.strand('bottom',[Segment('published44ntDNA',BOTTOM_5_TO_3)],mod5='App'if activated else'');sc.note('bottom','Sourceinternalbiotinmarkpresent; exactcovalentattachmentpositionnotatom-mapped');
    if activated:sc.ligation_end('bottom','published44ntDNA',kind='single-strand',polarity='5-prime',text='Knownsource5App donor in FREEactivatedlinker; afterligationactivationAMPisreleased,notbackbone')
    return sc
def top_scene():
    sc=Scene();sc.strand('top',[Segment('published44ntDNA',TOP)],mod5='P');return sc
def duplex_scene(joined=False,inferred=True):
    if TOP[3:]!=revcomp(BOTTOM_5_TO_3[3:]):raise ValueError('Source41bplinkercoredoesnotpair')
    sc=Scene();sc.strand('top',[Segment('source5overhang','TGA'),Segment('pairedtopDNA',TOP[3:])],mod5='P')
    parts=([Segment('INFERREDupstreamRNA','[sourceRNA]',placeholder=True,inferred=inferred,note='FullnativeRNA/endprocessingunassembled; intended3OHacceptorviewinferred') ]if joined else[])+[Segment('sourcebottom5overhang',BOTTOM_5_TO_3[:3]),Segment('pairedbottomDNA',BOTTOM_5_TO_3[3:])]
    sc.anneal('bottom',parts,to='top',pair=('pairedbottomDNA','pairedtopDNA'),mod5=''if joined else'App');sc.note('bottom','Sourceinternalbiotinpresent. PostligationpartialviewhasRNA-DNAjunctionandNOretained5App; fullRNA/productgeometryunassembled'if joined else'Sourcefreelinker5App/internalbiotin; RNArecipientnotincluded');return sc

def protocol_graph():
    g=ProtocolGraph('PRIM-seq publishedRNA-linkerCOMPONENTONLY; wholeNGSpreparationunassembled',authority=AUTHORITY,representation='published_instance')
    dna=g.input('bottom_dna','ActualsourceFig9bottomDNA/internalbiotin; preactivation5endchemistryunassigned',K.OLIGO,state=bottom_scene(False),primary=True,chemical_state={'source_3_to_5_ordering':BOTTOM_AS_PRINTED_3_TO_5,'five_to_three_transcription':'reverseprintedcharacters, notcomplement','biotin':'Sourceinternalmarkpresent; exactatompositionunassembled'})
    supplies=g.input('adenylation_supplies','Actualsource5DNAadenylationreagentinventory; identities/doses/packagingunassembled',K.MOLECULAR)
    activated=g.product('activated_donor','Source5App-activatedbottomDNA; completeadenylationmixture/byproductfatesunassembled',K.ADAPTER,state=bottom_scene(),chemical_state={'five_prime_activation':'Appactivatingadenylate, notDNAbackboneresidue','DNA_length':44,'biotin':'Internalmarkpresent'})
    g.react('adenylate','SourceFig9a→b explicitly5primeDNAadenylation',S.MAGIC,((dna,R.TEMPLATE),(supplies,R.STRUCTURAL)),(activated,),reason='Knownsource5App productcomponent; exactreactioninventory/chemistry/coproducts/occupancyunassembled.')
    rna=g.input('rna','SourceRNArecipientinFig9; exactsequence/quantities/endprocessingunassembled, intended3OHclassderived',K.BIOLOGICAL,primary=True)
    lig_sup=g.input('ligation_supplies','ActualsourceRNA-linkerligationreagents; identity/dose/packagingunassembled',K.MOLECULAR)
    joined=g.product('joined_mixture','SourceRNA-linkedbottomDNA plusactivationAMP/coproducts/unreactedclassesCOPRESENT; noimaginedAMPphysicalseparation',K.MOLECULAR,chemical_state={'intended_junction':'RNA3O toDNA5P derivedfromsource5Appdonor/intendedligation; wholeoccupiedmoleculesunassembled','activation_AMP':'Releasedfromactivateddonorbyintendedligation; remainscopresentmixture, notDNABACKBONE; downstreamfateunassigned'})
    g.react('rna_ligate','SourceFig9b+c→d RNA-linkerbottomligation',S.MAGIC,((rna,R.TEMPLATE),(activated,R.STRUCTURAL),(lig_sup,R.STRUCTURAL)),(joined,),geometry=AdapterLigationGeometry.known('single-strand',bottom_scene(),statement='Knownsource5AppDNA donor; actualRNArecipient3OH/intendedjunctionviewinferred, wholeproductunassembled.'),reason='SourcepositiveRNA-linkerligation, noactualenzymes/doses/quantitativeconversionorwholeRNAinvented.')
    top=g.input('top_dna','Actualsource5PhostopDNA',K.OLIGO,state=top_scene())
    partial=g.product('annealed_partial','SourceFig9e partialRNA-linkedDNAduplexmixture includingcopresentAMP/unreactedclasses; wholePRIMreporterunassembled',K.MOLECULAR,state=duplex_scene(True))
    g.react('anneal','SourceFig9d+e→RNA-linkedduplexannealing',S.HYBRIDIZATION,((joined,R.TEMPLATE),(top,R.STRUCTURAL)),(partial,))
    g.mark_terminal(Terminal(partial.id,Disposition.UNRESOLVED,'Component-onlypartialmixture. FullSMART-display/REILISandactualNGSlibrarypreparationnotread; no fabricatedsequencingendpoint.'))
    g.add_source_item(SourceItem('source_orientation','PublishedSI1Fig9printedbottom3→5; AppatRIGHT5end, captionexplicit5DNAadenylation; top5Phos.', 'Scraper3Apporientationlabelisnotauthority. Reverseprintedbottombodyto5→3; activationAMPnotDNAbackbone.41bpduplexderived; wholeRNA/productinferred/unassembled.'))
    return g

def selftest():
    from checks import Check,Source,have
    from storage import sources_root
    c=Check();p=sources_root()/HERE.parent.name/'s41587-025-02780-z_41587_2025_2780_MOESM1_ESM.pdf.txt'
    if have(c,Source(p,'https://doi.org/10.1038/s41587-025-02780-z')):
        t=p.read_text();c('Wholeprinted3to5bottomsource',BOTTOM_AS_PRINTED_3_TO_5 in t,True);c('Whole5Phostopsource','/Phos/'+TOP in t,True);c('Sourceexplicit5DNAadenylation','5´ DNA adenylation'in t or'5′ DNA adenylation'in t,True)
    g=protocol_graph();sc=duplex_scene(True)
    c('Reverseprintedbottomnotrevcomp',BOTTOM_5_TO_3[::-1],BOTTOM_AS_PRINTED_3_TO_5)
    c('RNAviewinferreddefault',sc.strands['bottom'].segs[0].inferred,True)
    c('PostligationviewnoactivationApp',sc.strands['bottom'].mod5,'')
    c('FreedonorhasactivationApp',bottom_scene().strands['bottom'].mod5,'App')
    c('Annealingusespostligationmixture',g.steps['anneal'].inputs[0].material,'joined_mixture')
    c('NoimaginedAMPseparation',all(m.id!='released_amp'for m in g.materials.values()),True)
    c('NoNGSlibraryinvented',g.sequencing,{})
    c('Componentisnotcompleteprotocol',bool(g.publication_problems()),True);c.report()
if __name__=='__main__':selftest()
