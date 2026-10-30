"""Protocol-owned source-qualified partial chemistry. No public schematic endpoint."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib'))
from chemdraw import Scene, Segment
from protocol_graph import ProtocolGraph, MaterialKind as K, InputRole as R, StepKind as S, SourceItem, SequencingPlan, Terminal, Disposition
from pcr import PCRGeometry
from adapter_ligation import AdapterLigationGeometry

def scene(name, text):
    sc = Scene(); sc.strand(name, [Segment(name, text, placeholder=True)]); return sc

def op(g, key, x, label, kind=S.MAGIC, *, reagent='Source mixture; internal chemistry unassembled', product=None, structural=()):
    catalyst = kind in {S.DIGESTION, S.END_REPAIR, S.EXTENSION, S.TAILING, S.REVERSE_TRANSCRIPTION, S.PCR, S.LIGATION, S.DEPHOSPHORYLATION, S.PHOSPHORYLATION}
    stock = g.input(key + '_reagent', reagent, K.CATALYST if catalyst else K.CONDITION)
    inputs = [(x, R.TEMPLATE), (stock, R.CATALYST if catalyst else R.CONDITION)]
    inputs.extend((a, R.STRUCTURAL) for a in structural)
    kwargs = {}
    if kind == S.MAGIC:
        kwargs['reason'] = 'Source operation established; exact molecular transition/kit internals unassembled. No accepted absence claim.'
    if kind == S.PCR:
        kwargs['geometry'] = PCRGeometry.unavailable(Scene(), reason='Exact template/primer endpoint geometry unassembled; selected oligos do not establish whole amplicon.')
        if len(structural) < 2:
            for arm in ('f', 'r'):
                inputs.append((g.input(key + '_' + arm, 'Source PCR primer; exact identity/landing unassembled', K.OLIGO, state=scene(key + '_' + arm, '[source primer]')), R.STRUCTURAL))
    if kind == S.REVERSE_TRANSCRIPTION and not structural:
        inputs.append((g.input(key + '_primer', 'Source RT primer; geometry unassembled', K.OLIGO, state=scene(key + '_primer', '[source RT primer]')), R.STRUCTURAL))
    if kind == S.LIGATION:
        kwargs['geometry'] = AdapterLigationGeometry('unavailable', 'unavailable', (scene(key + '_ends', '[source ligation; complete substrate/end geometry unassembled]'),), 'This incomplete model does not claim publication absence.')
        if not structural:
            inputs.append((g.input(key + '_adapter', 'Source ligating adapter; full geometry unassembled', K.ADAPTER, state=scene(key + '_adapter', '[source adapter]')), R.STRUCTURAL))
    y = g.product(key + '_product', product or label + '; physical companions retained until separation', K.MOLECULAR)
    g.react(key, label, kind, inputs, (y,), **kwargs); return y

def clean(g, key, x, label):
    condition = g.input(key + '_condition', label, K.CONDITION)
    kept = g.product(key + '_kept', 'Source physically recovered fraction; carryover/yield unquantified', K.MOLECULAR)
    waste = g.product(key + '_waste', 'Actual separated cleanup fractions; no perfect classification assumed', K.WASTE)
    g.react(key, label, S.CLEANUP, ((x, R.TEMPLATE), (condition, R.CONDITION)), (kept, waste))
    g.discard(waste, 'Actual physical separation rather than computational removal.'); return kept

def qc(g, key, x, label):
    kept = g.product(key + '_kept', 'Source retained preparation', K.MOLECULAR)
    assay = g.product(key + '_aliquot', 'Separate analytical aliquot; no return to library assumed', K.QC)
    g.react(key, label, S.SPLIT, ((x, R.TEMPLATE),), (kept, assay))
    g.mark_terminal(Terminal(assay.id, Disposition.QC, 'Source analytical aliquot.')); return kept

SOURCE='2016detailedprotocol associatedacceptedmanuscriptPMC4852555 Introduction/Proceduresteps11–66;2015definingmain/SIunavailable. Noequivalenceorversionmergingclaimed.'
def split(g,key,x,label,used,remaining):
    a=g.product(key+'_used',used,K.MOLECULAR);b=g.product(key+'_remainder',remaining,K.MOLECULAR)
    g.react(key,label,S.SPLIT,((x,R.TEMPLATE),),(a,b));g.mark_terminal(Terminal(b.id,Disposition.UNRESOLVED,'Actualremainingfraction;fatenotassigned'));return a

def protocol_graph():
    g=ProtocolGraph('cP-RNA-seq —2016associateddetailedprotocol, not2015reconstruction',authority=SOURCE,representation='published_instance')
    x=g.input('total_rna','SourceextractedtotalRNA;targetsizeclassesstillmixedbyterminalchemistry',K.MOLECULAR,primary=True)
    x=clean(g,'rna_page',x,'Source15percentdenaturingPAGEtargetlengthwindow/elution/isopropanolprecipitation;50ulexample')
    x=split(g,'cip_input',x,'20ulof50ulexamplepurifiedRNA','20ulCIPinput','30ulremainingRNA')
    x=op(g,'cip',x,'CIP37C40min;removesaccessible5prime/3primephosphates,not3primecP',S.DEPHOSPHORYLATION,reagent='CIP/CutSmart/RNasin',product='RNAterminalclasseswithCIP-removedphosphates;3primecPclassesretained; no physicalclassdisappearance')
    x=clean(g,'cip_clean',x,'Acidphenol-chloroform/ethanolprecipitation;45ulelution')
    x=op(g,'periodate',x,'FreshNaIO4final10mM40minoniceindark',S.CHEMICAL_CONVERSION,reagent='NaIO4',product='Oxidized3primecis-diolRNAclasses plusprotectedcP/otherresistantends;notperfectspecificityorphysicaldepletion')
    x=clean(g,'periodate_clean',x,'Ethanolprecipitation;10ulelution')
    x=op(g,'pnk',x,'T4PNK/ATP37C40min restores5primeP/opens3primecP',S.END_REPAIR,reagent='T4PNK/ATP/RNasin',product='LigationcompetentcP-derivedRNAclasses plusoxidized/resistant/nontargetRNAclasses/reagents;efficiencyunpredicted')
    x=clean(g,'pnk_clean',x,'Phenol-chloroform/ethanolprecipitation;10ulelution')
    x=split(g,'ligation_input',x,'5ulof10ulPNK-treatedRNA','5ulRA3ligationinput','5ulremainingPNK-treatedRNA')
    x=op(g,'ra3',x,'RA3anneal70C2minice1min;truncatedT4RNL2 28C1h/4Covernight',S.LIGATION,reagent='HistoricalTruSeqSmallRNA RA3/HML/RNaseinhibitor/truncatedT4RNL2')
    x=op(g,'stop',x,'STP28C15min;nointerveningcleanup',S.QUENCHING,reagent='SourceSTPstopbuffer',product='StoppedRA3reactionincludingunligatedRNAs/adapters/reagents')
    x=op(g,'ra5',x,'SeparatelyconditionedRA5/ATP/T4RNL mix;3uladded;28C1h/4Covernight',S.LIGATION,reagent='RA5/ATP/T4RNAligase',product='Dual-adaptedRNAclassesplusunligatedRNA/adapters/reagents')
    x=split(g,'rt_input',x,'Source6ulofadapterligationreaction','6ulRTinput','Remainingligationreaction;allocationunassigned')
    x=op(g,'rt',x,'RTPprimer70C2minthenSSIII50C60min',S.REVERSE_TRANSCRIPTION,reagent='SSIII/DTT/dNTP/RNaseinhibitor')
    x=op(g,'pcr',x,'Add37.5ulPCRmixtotheRTreaction;11–15cycles',S.PCR,reagent='TruSeqPML/RP1/RPIX;RPIXissampleindex,notUMI')
    x=clean(g,'native_page',x,'8percentnativePAGE;imageandselecttargetcDNAbands;finalrecoverydetailunassembled')
    lib=g.product('library','SelectedcP-RNA-seqcDNAclasses;exactkitadapters/junctionsnotassembled',K.LIBRARY)
    g.react('library_identity','Partialphysicalsequencinginterface',S.GENERIC,((x,R.TEMPLATE),),(lib,),changes_molecule=False)
    g.sequence(lib,SequencingPlan(limitation='2016detailedprotocolcallsIlluminaNGSbutacquisition/bioinformaticsoutsidescope;2015actualrununreviewed. Noinstrument/readcyclesinferred.'))
    g.add_source_item(SourceItem('pending','2015original/SI;final2016article;TruSeqcomponentrevision;nativePAGEfinalrecovery/acquisition','AM-qualifiedpartialmodel,notpublishedabsence.'));return g

def selftest():
    from checks import Check
    c=Check();g=protocol_graph()
    c('PeriodateactsafterCIPcleanup',g.steps['periodate'].inputs[0].material,'cip_clean_kept')
    c('PNKactsafterperiodatecleanup',g.steps['pnk'].inputs[0].material,'periodate_clean_kept')
    c('RA3usesphysicalhalfinput',g.steps['ra3'].inputs[0].material,'ligation_input_used')
    c('RA5retainsstoppedRA3reaction',g.steps['ra5'].inputs[0].material,'stop_product')
    c('RTuses6ulphysicalaliquot',g.steps['rt'].inputs[0].material,'rt_input_used')
    c('PCRdirectlyreceivesRTreaction',g.steps['pcr'].inputs[0].material,'rt_product')
    c('RemainingPNKRNAexplicit',g.terminals['ligation_input_remainder'].disposition,Disposition.UNRESOLVED)
    c('Incompletegraphnotpublishable',bool(g.publication_problems()),True)
    c.report()
if __name__=='__main__':selftest()
