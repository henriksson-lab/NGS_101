"""Published2022 MALBAC-DT: source-owned original and automated plate routes."""
from chemdraw import Construct,Scene,feature
from protocol_model import role
from protocol_graph import (ProtocolGraph,StepKind,MaterialKind,InputRole,SequencingPlan,
                            Terminal,Disposition,Highlight)
from pcr import PCRGeometry
import seqprimers as sp

TITLE='MALBAC-DT — high-detectability single-cell transcriptomics'
NOTES='01_malbac_dt.html'
SOURCE='<a href="https://doi.org/10.1073/pnas.2206938119">Chapman, Lee et al., PNAS2022, PMC9907105</a>, published Methods original/automated MALBAC-DT and Sequence Processing. The directory names the superseded preprint; it is not active recipe authority.'
SUMMARY='Barcoded poly(T) reverse transcription, excess-primer digestion and MALBAC-style amplification precede Tru2-primer completion and source-specific library preparation. Original Nextera and automated TruePrep routes remain distinct.'
CAVEAT='INFERRED region order — inspected published passages name primer species and barcode/UMI roles but their native strands, internal order, junctions and selected kit/index geometry remain unassigned. No full canonical Nextera shell,8nt index, P7 end or exact customRead2 primer is inferred from names. Published Sequence Processing gives a20nt consensus UMI pattern; this does not establish the physical primer region order. Relevant linked protocol identity remains unresolved; no missing-information claim is accepted.'

WELL='malbac_dt_well'
INNER=feature('malbac_inner','cell_barcode','unknown',group=WELL,part='inner RT3-An row/well barcode')
OUTER=feature('malbac_outer','cell_barcode','unknown',group=WELL,part='outer Tru2-Gn-RT column barcode; automated uses one member')
UMI=feature('malbac_umi','umi','degenerate',note='Published consensus(HBDV)5, excluding RT-Bn(VDBH)5 similarity; physical primer pattern/orientation unassigned')
I5=feature('malbac_kit_index','sample_index','unknown',note='Original S5XX named; automated TruePrep native index architecture unassigned')
IX=feature('malbac_ix','sample_index','unknown',note='Original Ix-Tru2 named; selected member/physical site unassigned')
UMI_LENGTH=len('HBDV'*5)


def identity_regions(outer=False):
    regions=[role('RT-An inner barcode; native position/width unassigned',inferred=True,molecular_feature=INNER),role('UMI; physical primer location unassigned',length_bp=UMI_LENGTH,inferred=True,molecular_feature=UMI)]
    if outer:regions.append(role('Tru2-associated outer barcode; native position unassigned',inferred=True,molecular_feature=OUTER))
    return regions


def state(label,*,outer=False,rna=False,other=False,tagged=False):
    sc=Scene()
    parts=[role(label+' native cDNA/handle assembly unassigned',inferred=True),*identity_regions(outer)]
    if tagged:parts.append(role('source-kit transferred end/gaps; native arrangement unassigned',inferred=True))
    sc.strand('cDNA',parts,label='INFERRED display order; identifier-bearing cDNA lineage')
    if rna:sc.strand('RNA',[role('independent RNA templates; hybrid geometry unassigned')])
    if other:sc.strand('other',[role('independent excess primers, hydrolysis products and non-target copies; partition unassigned')])
    return sc


def final(name,auto=False):
    return Construct([role(('native TruePrep end region' if auto else 'native S5XX-associated Nextera end region')+'; full strand/index geometry unassigned',inferred=True,molecular_feature=I5),role('transcript-derived cDNA'),*identity_regions(True),role(('native automated Tru2-associated opposite end' if auto else 'native Ix-Tru2-associated opposite end')+'; selected primer/end geometry unassigned',inferred=True,molecular_feature=IX)],name=name)

ORIGINAL=final('original U2OS/HEK293T plate library')
AUTOMATED=final('automated K562 plate library',True)
FINAL_LIBRARY=ORIGINAL
SEQ_PRIMERS=()
FINAL_LIBRARIES=(('Original plate library',ORIGINAL,(),'INFERRED native regions — original Nextera/S5XX plus Ix-Tru2 route. Internal order and full primer junctions unassigned.','Original acquired run layout unassigned.'),('Automated plate library',AUTOMATED,(),'INFERRED native regions — distinct VazymeTruePrep route. Kit-end geometry is not inherited from original Nextera.','Published NovaSeq paired2×150 with customRead2 primer; physical sites unassigned.'))
SEQUENCING_UNAVAILABLE='Native kit/read/index primer geometry remains unassigned; custom Read2 is source-positive for automated plate but no sequence is supplied to this model.'
READ_SPECS={name:tuple(sp.read_spec(read,150 if name=='Automated plate library' and read in ('Read 1','Read 2') else 'unassigned','Transcript signal' if read=='Read 1' else 'Tru2-side cell barcode/UMI signal; exact entry unassigned' if read=='Read 2' else 'Index acquisition/cycles unassigned','include' if read=='Read 1' and name=='Automated plate library' else 'extract' if read=='Read 2' and name=='Automated plate library' else 'undetermined') for read in ('Read 1','Read 2','Index 1 (i7)','Index 2 (i5)')) for name,*_ in FINAL_LIBRARIES}
READ_SPECS_NOTE='Only automated plate2×150 is assigned. Kit/index presence does not assign acquired index cycles or original-plate read configuration.'

def oligos():return ()


def protocol_graph():
    g=ProtocolGraph(TITLE,authority=SOURCE)
    cells=g.input('cells','sorted cells configured for original or automated plate experiment',MaterialKind.BIOLOGICAL,primary=True)
    original=g.product('original_cells','original U2OS/HEK293T experiment cells',MaterialKind.MOLECULAR)
    automated=g.product('automated_cells','automated K562 experiment cells',MaterialKind.MOLECULAR)
    g.react('source_routes','Represent separately configured published experiments',StepKind.SPLIT,((cells,InputRole.TEMPLATE),),(original,automated),changes_molecule=False)
    for key,sample,endpoint,heading in (('original',original,ORIGINAL,'Original plate library'),('automated',automated,AUTOMATED,'Automated plate library')):
        rt_primer=g.input(key+'_rt_primer','RT-An/RT3-An;8row members original,96well members automated; native species unassigned',MaterialKind.OLIGO)
        rt=g.input(key+'_rt','SuperScriptIV;72°C3min beforeenzymeaddition,55°C10minRT' if key=='original' else 'source automatedRTmix; enzymeidentityunassigned,55°C10min',MaterialKind.CATALYST)
        cdna=g.product(key+'_cdna','actual first-strand cDNA/RNA reaction mixture',MaterialKind.MOLECULAR,state=state('first-strand',rna=True,other=True))
        g.react(key+'_reverse_transcribe','Reverse-transcribe source sorted-cell RNA',StepKind.REVERSE_TRANSCRIPTION,((sample,InputRole.TEMPLATE),(rt_primer,InputRole.STRUCTURAL),(rt,InputRole.CATALYST)),(cdna,))
        exo=g.input(key+'_exo','ExonucleaseI;37°C30min,80°C20min' if key=='original' else 'source automatedExomix; enzymeidentityunassigned,37°C30min/80°C20min',MaterialKind.CATALYST)
        bn=g.input(key+'_rt_bn','RT-Bn added in exonuclease mix; precise molecular role unassigned' if key=='original' else 'automatedExomix other components unassigned; original RT-Bn inheritance unproven',MaterialKind.OLIGO if key=='original' else MaterialKind.CONDITION)
        digested=g.product(key+'_digested','cDNA plus carried RNA/primer hydrolysis products',MaterialKind.MOLECULAR,state=cdna.state)
        g.react(key+'_exo_digest','Digest excess RT primers without a physical separation',StepKind.DIGESTION,((cdna,InputRole.TEMPLATE),(exo,InputRole.CATALYST),(bn,InputRole.STRUCTURAL if key=='original' else InputRole.CONDITION)),(digested,))
        random=g.input(key+'_gat7n','GAT-7N random/common primer role; native species unassigned' if key=='original' else 'automatedPCRmix first primer role; speciesunassigned',MaterialKind.OLIGO)
        common=g.input(key+'_gatcom','GAT-COM common primer role; native species unassigned' if key=='original' else 'automatedPCRmix opposite primer role; speciesunassigned',MaterialKind.OLIGO)
        pol=g.input(key+'_pol','DeepVentexo− in original; automated master mix source-delegated',MaterialKind.CATALYST)
        amp=g.product(key+'_malbac_product','MALBAC-style amplified cDNA population',MaterialKind.MOLECULAR,state=state('MALBAC-derived',rna=True,other=True))
        if key=='original':
            looped=g.product('original_looped','source10loopingcycle product population; native loop/copy geometry unassigned',MaterialKind.MOLECULAR,state=state('looping-cycle',rna=True,other=True))
            g.react('original_looping','Source10loopingcycles before PCR',StepKind.EXTENSION,((digested,InputRole.TEMPLATE),(random,InputRole.STRUCTURAL),(common,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(looped,))
        else:looped=digested
        g.react(key+'_malbac','Original subsequent PCR: printed Goto13 seventeen times' if key=='original' else 'Automated15PCRcycles',StepKind.PCR,((looped,InputRole.TEMPLATE),(random,InputRole.STRUCTURAL),(common,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(amp,),geometry=PCRGeometry.unavailable(looped.state,reason='Actual incoming looping-cycle population (original) or Exomix-treated cDNA (automated). Native primer/copy/loop geometry unassigned; no completed duplex substituted.'))
        tru2=g.input(key+'_tru2','Tru2-Gn-RT12column members original; one Tru2-G-RT automated; native species unassigned',MaterialKind.OLIGO)
        finished=g.product(key+'_finished','Tru2-completed cDNA population',MaterialKind.MOLECULAR,state=state('Tru2-primer-derived',outer=True,rna=True,other=True))
        g.react(key+'_tru2_pcr','Add source Tru2 primer and amplify5morecycles',StepKind.PCR,((amp,InputRole.TEMPLATE),(tru2,InputRole.STRUCTURAL),(common,InputRole.STRUCTURAL),(pol,InputRole.CATALYST)),(finished,),geometry=PCRGeometry.unavailable(amp.state,reason='Actual incoming MALBAC cDNA population retained. Native Tru2/GAT binding and first-copy assembly are unassigned; completed product is not the template.'))
        beads=g.input(key+'_beads','Original1µl/well pooling and0.8×AMPure; automated equivalent cleanup not established by inspected summary',MaterialKind.SELECTION if key=='original' else MaterialKind.CONDITION)
        recovered=g.product(key+'_recovered','pooled cDNA for source-kit library preparation',MaterialKind.MOLECULAR,state=finished.state)
        if key=='original':
            waste=g.product(key+'_bead_waste','physically removed bead supernatant/washes',MaterialKind.WASTE)
            g.react(key+'_pool','Pool and recover original plate cDNA',StepKind.CLEANUP,((finished,InputRole.TEMPLATE),(beads,InputRole.SELECTION)),(recovered,waste),changes_molecule=False);g.discard(waste,'Source bead fraction transfer; individual RNA/hydrolysate recovery unassigned.')
        else:g.react(key+'_pool','Automated input recovery remains source-delegated',StepKind.MAGIC,((finished,InputRole.TEMPLATE),(beads,InputRole.CONDITION)),(recovered,),reason='Automated summary names downstreamTruePrep but does not establish original0.8×pooling cleanup equivalence.')
        adapter=g.input(key+'_transposome','VazymeTruePrep native loaded adapter' if key=='automated' else 'original Nextera native loaded adapter',MaterialKind.ADAPTER)
        tn5=g.input(key+'_tn5','source-kit Tn5',MaterialKind.CATALYST)
        tagged=g.product(key+'_tagged','native pre-PCR tagged barcode-retaining cDNA classes',MaterialKind.MOLECULAR,state=state('source-kit tagged',outer=True,rna=True,other=True,tagged=True))
        other=g.product(key+'_other_ends','internal/nonselected tagmentation classes',MaterialKind.MOLECULAR)
        g.react(key+'_tagment','Prepare source-specific tagged cDNA',StepKind.TAGMENTATION,((recovered,InputRole.TEMPLATE),(adapter,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),(tagged,other))
        g.mark_terminal(Terminal(other.id,Disposition.UNRESOLVED,'Lacking Tru2 amplification compatibility does not physically discard molecules; later recovery is unassigned.'))
        forward=g.input(key+'_library_f','nativeTruePrep forward primer unassigned' if key=='automated' else 'NexteraS5XX; selected member/index width unassigned',MaterialKind.OLIGO)
        reverse=g.input(key+'_library_r','nativeautomated Tru2-side kit/index completion primer unassigned' if key=='automated' else 'Ix-Tru2 replacesS7XX; native strand/selected member unassigned',MaterialKind.OLIGO)
        pcr=g.input(key+'_library_pol','source-kit PCR polymerase/cycles unassigned',MaterialKind.CATALYST)
        library=g.product(key+'_library',heading,MaterialKind.LIBRARY,state=endpoint,inferred=True)
        carried=g.product(key+'_carried_species','carried RNA/primer/nontarget populations after libraryPCR; native recoveryunassigned',MaterialKind.MOLECULAR,state=state('carried nonselected species',outer=True,rna=True,other=True))
        g.react(key+'_library_pcr','Complete source-specific library from actual native input',StepKind.PCR,((tagged,InputRole.TEMPLATE),(forward,InputRole.STRUCTURAL),(reverse,InputRole.STRUCTURAL),(pcr,InputRole.CATALYST)),(library,carried),geometry=PCRGeometry.unavailable(tagged.state,reason='Actual tagged cDNA input retained. Native loaded adapter/PCR/index strands and full end assembly remain unassigned; no canonical finished duplex or default index width substituted.'))
        g.mark_terminal(Terminal(carried.id,Disposition.UNRESOLVED,'PCR does not itself physically remove carriedRNA, primers or nonselected copies. Subsequent kit recovery is not reconstructed from inspected source scope.'))
        g.sequence(library,SequencingPlan(library=endpoint,acquired_roles=('Read 1','Read 2') if key=='automated' else (),read_lengths={'Read 1':150,'Read 2':150} if key=='automated' else {},limitation=SEQUENCING_UNAVAILABLE))
    g.annotate(Highlight('original_malbac','Printed loop/PCR schedule remains source-qualified','Original program says Goto2 ten times and Goto13 seventeen times; the exact iteration interpretation is unresolved rather than silently choosing17or18PCRcycles. Automated15+5cycles is a distinct positive source statement.'))
    g.annotate(Highlight('automated_pool','Historical source boundary','Published2022 source supersedes preprint. Original Nextera/S7XXreplacement is not transferred to the automatedTruePrep/customRead2 route. Native RNA and hydrolysates are carried without declaring physical disappearance; computational UMI/barcode filters are not separations.'))
    g.declare_reviewed_scope(REVIEWED,SOURCE)
    return g.require_publishable()

REVIEWED=('cells','source_routes','original_rt_primer','original_rt','original_reverse_transcribe','original_exo','original_rt_bn','original_exo_digest','original_gat7n','original_gatcom','original_pol','original_looping','original_malbac','original_tru2','original_tru2_pcr','original_beads','original_pool','original_transposome','original_tn5','original_tagment','original_library_f','original_library_r','original_library_pol','original_library_pcr','original_library','automated_rt_primer','automated_rt','automated_reverse_transcribe','automated_exo','automated_rt_bn','automated_exo_digest','automated_gat7n','automated_gatcom','automated_pol','automated_malbac','automated_tru2','automated_tru2_pcr','automated_beads','automated_pool','automated_transposome','automated_tn5','automated_tagment','automated_library_f','automated_library_r','automated_library_pol','automated_library_pcr','automated_library')
