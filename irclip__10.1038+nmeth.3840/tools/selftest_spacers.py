"""Synthetic topology tests plus optional transcription check against the cached source."""
from irclip import (rt_primer_backbone, linear_cdna_backbone, conditional_first_extension,
                    RT_WRITTEN, protocol_graph)
from spacer_backbone import DNA, C18, Backbone, Landing, StopHypothesis, extend_to_linker
from checks import Check, Source, have
from chemdraw import Construct
from storage import sources_root


def rejects(fn):
    try: fn()
    except ValueError: return True
    return False


def run(check):
    policy=StopHypothesis(frozenset({'/iSp18/'}),'Synthetic test assumption: stop at the explicit linker.')
    template=Backbone((DNA('head','ACGT'),C18('barrier'),DNA('after','CCGA'),DNA('tail','ATGC'))).close()
    # Primer AC pairs template head[2:4]=GT. Copy head AC, wrap, copy ATGC then CCGA.
    result=extend_to_linker(template,Landing(0,2,4,'AC'),hypothesis=policy)
    check('synthetic extension copies reverse complements through circle wraparound',
          result.product.exact_sequence(),'ACGTGCATTCGG')
    check('synthetic traversal retains exact visited intervals and linker stop',
          ([(x.node_index,x.template_start,x.template_end) for x in result.copied_path],
           result.stopped_at,result.closure_crossings), ([(0,0,2),(3,0,4),(2,0,4)],1,1))
    ambiguous=Backbone((DNA('head','ACGT'),C18('barrier'),DNA('tail','ARYK'))).close()
    check('copying complements IUPAC classes rather than only reversing their letters',
          extend_to_linker(ambiguous,Landing(0,2,4,'AC'),hypothesis=policy).product.exact_sequence(),
          'ACGTMRYT')
    immediate=extend_to_linker(template,Landing(2,0,4,'TCGG'),hypothesis=policy)
    check('polymerase does not add bases when linker is immediately before primer',
          (immediate.product.exact_sequence(),immediate.copied_path,immediate.stopped_at),('TCGG',(),1))
    check('no unspecified linker traversal or bypass',rejects(lambda:extend_to_linker(
        template,Landing(0,2,4,'AC'),hypothesis=StopHypothesis(frozenset(),'No linker behavior assigned.'))),True)
    no_linker=Backbone((DNA('head','ACGT'),DNA('tail','ATGC'))).close()
    check('no linker-free rolling-circle product is invented',rejects(lambda:extend_to_linker(
        no_linker,Landing(0,2,4,'AC'),hypothesis=policy)),True)
    check('linear molecule cannot silently become circular traversal',rejects(lambda:extend_to_linker(
        Backbone((DNA('head','ACGT'),DNA('tail','ATGC'))),Landing(0,2,4,'AC'),hypothesis=policy)),True)
    check('invalid primer landing rejected before copying',rejects(lambda:extend_to_linker(
        template,Landing(0,2,4,'AA'),hypothesis=policy)),True)
    check('5-prime OH cannot silently close as phosphorylated DNA',
          rejects(lambda:Backbone((DNA('unphosphorylated','ACGT'),),five_prime='OH').close()),True)

    primer=rt_primer_backbone();linear=linear_cdna_backbone();circle=linear.close()
    check('source primer written form preserves both non-nucleotide linkers',primer.written(),RT_WRITTEN)
    check('ordered cDNA retains every primer node and both linker objects',
          (linear.nodes[:len(primer.nodes)]==primer.nodes,
           sum(isinstance(x,C18) for x in circle.nodes),circle.nodes==linear.nodes),(True,2,True))
    check('non-nucleotide backbone cannot be passed off as a DNA sequence',
          rejects(circle.exact_sequence),True)
    extension=conditional_first_extension(circle)
    check('source-construct conditional copy stops at second C18 after one closure crossing',
          (circle.nodes[extension.stopped_at].name,extension.closure_crossings),('second C18',1))
    check('RNA-derived block changes copy orientation without inventing sequence/length',
          [(x.origin,x.reverse_complemented) for x in extension.product.nodes
           if hasattr(x,'origin')],[('captured RNA',False)])
    original_features={x.feature.id:x.feature for x in primer.nodes if isinstance(x,DNA) and x.feature}
    copied_features={x.feature.id:x.feature for x in extension.product.nodes if isinstance(x,DNA) and x.feature}
    umi_parts=sorted((f.id,f.group,f.part,f.role) for f in original_features.values() if f.role=='umi')
    check('multipart UMI has distinct stable part IDs and one shared group',umi_parts,[
          ('irclip_degenerate_identifier_n3','irclip_degenerate_identifier','N3','umi'),
          ('irclip_degenerate_identifier_n5','irclip_degenerate_identifier','N5','umi'),
          ('irclip_degenerate_identifier_w3','irclip_degenerate_identifier','W3','umi')])
    check('conditional copying preserves complete identifier definitions and grouping',copied_features,original_features)
    represented=Construct(primer.display_segments(),name='source primer annotation with C18 retained')
    check('source primer forms a valid Construct without deleting either C18 placeholder',
          [(s.name,s.top,s.placeholder) for s in represented if s.name.endswith('C18')],
          [('first C18','[C18 linker]',True),('second C18','[C18 linker]',True)])
    check('conditional product is not an exact sequence when RNA identity is unassigned',
          rejects(extension.product.exact_sequence),True)
    g=protocol_graph()
    check('graph retains actual spacer circle and labels first extension inferred',
          (g.materials['cdna_circle'].chemical_state.backbone.nodes==circle.nodes,
           g.materials['conditional_first_strand'].inferred),(True,True))

    original=Source(sources_root()/'irclip__10.1038+nmeth.3840'/'nmeth.3840_PMC5477425.xml.txt',
                    'python3 tools/get_sources.py --slugs irclip__10.1038+nmeth.3840; accepted manuscript Methods71')
    if have(check,original,label='irCLIP modified-primer source transcription'):
        normalized=''.join(original.path.read_text().split()).replace('/iSpl8/','/iSp18/')
        check('source prints phosphorylated primer with both carbon-spacer tokens',
              RT_WRITTEN in normalized,True)
        check('source states linker function without proving this precise stop mechanism',
              'requiredforgenerationoflinearPCRproductusingcircularizedcDNAastemplate' in normalized,True)

if __name__=='__main__':
    c=Check();run(c);c.report()
