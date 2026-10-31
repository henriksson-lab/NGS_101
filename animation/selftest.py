"""Optional authoring checks: python3 -m animation.selftest [--native]."""
from dataclasses import replace
import json
import subprocess
import tempfile
from pathlib import Path
import sys
import unittest
from .model import Movie, Event
from .examples import reference, smartseq


class MolecularTests(unittest.TestCase):
    def test_copy_rejects_wrong_base(self):
        m=reference();m.shots[0].strands[1].sequence='A'*len(m.shots[0].strands[1].sequence)
        with self.assertRaisesRegex(ValueError,'complement'):m.validate()

    def test_copy_rejects_blocked_end(self):
        m=reference();m.shots[0].strands[1].three_prime='ddC'
        with self.assertRaisesRegex(ValueError,'extendable'):m.validate()

    def test_mapping_rejects_wrong_coordinates(self):
        m=reference();m.shots[0].events[0].template_start=999
        with self.assertRaisesRegex(ValueError,'leaves template'):m.validate()

    def test_unknown_is_not_exact(self):
        m=reference();m.shots[0].strands[0].schematic=True
        with self.assertRaisesRegex(ValueError,'unknown classes'):m.validate()

    def test_speech_extends_hold_not_action(self):
        shot=reference().shots[0];short=shot.timing();long=replace(shot,narration_seconds=50).timing()
        self.assertEqual(short[:2],long[:2]);self.assertEqual(long[2],51)

    def test_duration_does_not_truncate(self):
        shot=replace(reference().shots[0],narration_seconds=50,duration=10)
        with self.assertRaisesRegex(ValueError,'truncate'):shot.timing()

    def test_roundtrip_and_seek(self):
        m=smartseq();r=Movie.from_dict(json.loads(json.dumps(m.to_dict())))
        for t in [0,5,25,10,0,10000]:
            self.assertEqual(m.at(t)[0].id,r.at(t)[0].id)
        self.assertIn('not SMART-seq2',r.provenance['original_protocol'])
        self.assertTrue(all(s.schematic for shot in r.shots for s in shot.strands))

    def test_invalid_speed_and_timing(self):
        for speed in [0,-1,float('nan')]:
            with self.assertRaises(ValueError):reference().duration(speed=speed)
        with self.assertRaises(ValueError):reference().duration(wpm=0)

    def test_speech_estimate_uses_chosen_rate_at_export(self):
        m=reference();shot=m.shots[0];shot.narration=' '.join(['word']*50);shot.duration=12
        m.validate()
        self.assertEqual(m.duration(wpm=300,speed=1),12)
        with self.assertRaisesRegex(ValueError,'truncate'):m.duration(wpm=100)

    def test_version_is_checked(self):
        d=reference().to_dict();d['schema']='future/99'
        with self.assertRaisesRegex(ValueError,'schema'):Movie.from_dict(d)

    def test_ligation_not_exact_without_ends(self):
        m=reference();m.shots[0].strands[1].nicks=[8];m.shots[0].events=[Event('seal','daughter',bond=8,exact=True)]
        with self.assertRaisesRegex(ValueError,'two-end'):m.validate()


class ReactionTests(unittest.TestCase):
    def test_closure_checks_actual_end_modifications(self):
        from .protocols import crispr_mip
        m=crispr_mip();m.shots[-1].strands[0].five_prime='OH'
        with self.assertRaisesRegex(ValueError,'terminal modifications'):m.validate()

    def test_blocked_junction_cannot_ligate(self):
        from .protocols import smallseq
        m=smallseq();next(e for e in m.shots[1].events if e.kind=='seal').junction.acceptor_3prime='ddC'
        with self.assertRaisesRegex(ValueError,'3primeOH'):m.validate()

    def test_long_probe_backbone_is_retained_not_forced_to_pair(self):
        from .protocols import crispr_mip
        m=crispr_mip();shot=m.shots[1];event=shot.events[0]
        self.assertTrue(event.exact)
        self.assertEqual(event.primer_binding_length,19)
        self.assertGreater(event.primer_length,event.primer_binding_length)
        sequence=shot.strands[0].sequence
        shot.strands[0].sequence=sequence[:event.primer_length]+'A'+sequence[event.primer_length+1:]
        with self.assertRaisesRegex(ValueError,'complement'):m.validate()

    def test_fill_must_precede_closure(self):
        from .protocols import crispr_mip,closure_event
        m=crispr_mip();m.shots[1].events.append(closure_event('probe','Regression end fixture',start=0))
        with self.assertRaisesRegex(ValueError,'before extension'):m.validate()

    def test_bending_does_not_create_a_covalent_circle(self):
        from .evaluate import evaluate
        from .protocols import crispr_mip
        s=crispr_mip().shots[0]
        states,_=evaluate(s,s.establish+s.action)
        self.assertFalse(states['probe'].closed)

    def test_sealing_preserves_polymer_and_closes_exactly_one_nick(self):
        from .evaluate import evaluate
        from .protocols import crispr_mip
        s=crispr_mip().shots[-1];before,_=evaluate(s,0);after,_=evaluate(s,s.establish+s.action)
        self.assertEqual(before['probe'].strand.sequence,after['probe'].strand.sequence)
        self.assertFalse(before['probe'].closed);self.assertTrue(after['probe'].closed)
        self.assertEqual(after['probe'].nicks[0],0)

    def test_circle_requires_ordered_closed_route(self):
        from .protocols import boloramis
        m=boloramis();m.shots[-1].strands[0].path[-1]=(len(m.shots[-1].strands[0].sequence),1,1)
        with self.assertRaisesRegex(ValueError,'geometrically closed'):m.validate()

    def test_rca_cannot_copy_linear_template_as_a_circle(self):
        from .protocols import boloramis
        m=boloramis();m.shots[-1].strands[0].circular=False
        with self.assertRaisesRegex(ValueError,'covalently closed'):m.validate()

    def test_rca_copied_repeat_and_origin_are_checked(self):
        from .protocols import boloramis
        m=boloramis();s=m.shots[-1].strands[1]
        s.sequence=s.sequence[:-1]+('A'if s.sequence[-1]!='A'else'C')
        with self.assertRaisesRegex(ValueError,'complement'):m.validate()

    def test_rca_threeprime_front_stays_near_template_polymerase(self):
        from .protocols import boloramis
        from .evaluate import evaluate
        import math
        shot=boloramis().shots[-1]
        for f in [0,.2,.6,1]:
            states,actors=evaluate(shot,shot.establish+shot.action*f)
            child=states['concatemer'];front=child.point(child.visible[0][1])
            self.assertLess(math.dist(front,actors[0].position),20)

    def test_retained_identity_cannot_silently_change_bases(self):
        from .protocols import smallseq
        m=smallseq();shot=m.shots[1]
        shot.inventory['ra3_input']='A'*len(shot.inventory['ra3_input'])
        with self.assertRaisesRegex(ValueError,'original molecule'):m.validate()

    def test_pairing_rejects_parallel_or_noncomplementary_arms(self):
        from .protocols import crispr_mip
        m=crispr_mip();m.shots[0].pairings[0].b_step=1
        with self.assertRaisesRegex(ValueError,'antiparallel'):m.validate()
        m=crispr_mip();m.shots[0].pairings[0].b_start+=1
        with self.assertRaisesRegex(ValueError,'complementary'):m.validate()

    def test_conversion_retains_protected_c_and_physical_dna_u(self):
        from .protocols import scrrbs
        from .evaluate import evaluate
        m=scrrbs();shot=m.shots[1];states,actors=evaluate(shot,1000)
        self.assertEqual(set(states['original'].converted.values()),{'U'})
        self.assertNotIn(3,states['original'].converted)
        self.assertEqual(states['original'].strand.polymer,'DNA')
        self.assertEqual(actors,[])  # Bisulfite is not a protein enzyme.

    def test_source_missing_slots_never_become_exact_claims(self):
        from .examples import EXAMPLES
        m=EXAMPLES['smartseq']()
        self.assertEqual(m.provenance['original_protocol'],'Ramskold2012; not SMART-seq2')
        self.assertIn('index',m.provenance['steps_used'])
        self.assertFalse(any(e.exact for s in m.shots for e in s.events))

    def test_ambiguous_base_family_does_not_establish_exact_copy(self):
        m=reference();m.shots[0].strands[0].sequence='N'*len(m.shots[0].strands[0].sequence)
        with self.assertRaisesRegex(ValueError,'unknown classes'):m.validate()

    def test_ordinary_model_use_does_not_require_animation_libraries(self):
        code="import sys;from animation.examples import EXAMPLES;EXAMPLES['smartseq']().validate();assert 'cairo'not in sys.modules"
        subprocess.run([sys.executable,'-S','-c',code],check=True)

    def test_all_protocols_have_explicit_source_bindings_and_roundtrip(self):
        from .examples import EXAMPLES
        for name,fn in EXAMPLES.items():
            with self.subTest(example=name):
                m=fn();r=Movie.from_dict(json.loads(json.dumps(m.to_dict())))
                self.assertEqual(json.dumps(m.to_dict(),sort_keys=True),json.dumps(r.to_dict(),sort_keys=True))
                if name!='reference':
                    self.assertEqual(set(s.id for s in m.shots),set(m.provenance['steps_used']))
                    self.assertTrue(m.provenance['module_sha256'])


class PresentationTests(unittest.TestCase):
    def test_default_speed_and_seek_share_multiplier(self):
        m=reference()
        self.assertAlmostEqual(m.duration()*1.5,m.duration(speed=1))
        for t in [0,1,3,5]:self.assertEqual(m.at(t),m.at(t*1.5,speed=1))

    def test_poly_a_and_primer_handle_are_physically_separate(self):
        m=smartseq();shot=m.shots[0]
        rna,primer=shot.strands;pair=shot.pairings[0]
        # Check the displayed representative member, independently of schematic body.
        from .model import complement
        for j in range(pair.length):
            self.assertEqual(primer.sequence[pair.a_start+j],complement(rna.sequence[pair.b_start-j]))
        self.assertTrue(primer.regions[0].end<=pair.a_start)
        self.assertEqual(primer.sequence[pair.a_start:pair.a_start+30],'T'*30)

    def test_rna_and_cdna_persist_through_template_switch(self):
        shots={s.id:s for s in smartseq().shots}
        for id in ['rt','terminal','switch']:
            rna=next(s for s in shots[id].strands if s.id=='rna')
            self.assertEqual(rna.sequence,shots['anchor'].strands[0].sequence)
            self.assertEqual(rna.start,shots['anchor'].strands[0].start)
        first=next(s for s in shots['firstcopy'].strands if s.id=='reverse')
        self.assertLess(first.initially_visible,len(first.sequence))
        self.assertTrue(any(e.kind=='extend'and e.strand=='reverse'for e in shots['firstcopy'].events))

    def test_transposition_has_dimers_and_retains_cleaved_pieces(self):
        from .protocols import tagmentation_shot
        from .evaluate import evaluate
        shot=tagmentation_shot();shot.source='Schematic test';shot.validate()
        state,actors=evaluate(shot,shot.establish+shot.action)
        self.assertEqual(sum(a.complex=='homodimer'for a in actors),2)
        top=state['top']
        self.assertEqual(top.visible,[(0.,140.)])
        self.assertLess(top.segment_point(30,29)[0],top.segment_point(30,30)[0])
        self.assertGreater(top.segment_point(105,106)[0],top.segment_point(105,104)[0])

    def test_index_overhang_is_not_forced_to_pair(self):
        shot=next(s for s in smartseq().shots if s.id=='index')
        event=next(e for e in shot.events if e.kind=='extend')
        self.assertLess(event.primer_binding_length,event.primer_length)
        self.assertTrue(any(r.label=='Sample index'for r in shot.strands[1].regions))
        self.assertEqual(len([e for e in shot.events if e.kind=='extend']),2)

    def test_functional_regions_and_initial_pieces_roundtrip(self):
        from .protocols import smallseq
        from .evaluate import evaluate
        m=smallseq();r=Movie.from_dict(json.loads(json.dumps(m.to_dict())))
        shot=next(s for s in r.shots if s.id=='lig3')
        start,_=evaluate(shot,0);end,_=evaluate(shot,shot.timing()[2])
        junction=shot.strands[0].nicks[0]
        self.assertGreater(start['rna_adapter'].segment_point(junction,junction)[0],start['rna_adapter'].segment_point(junction,junction-1)[0])
        self.assertEqual(end['rna_adapter'].segment_point(junction,junction),end['rna_adapter'].segment_point(junction,junction-1))
        self.assertEqual(end['rna_adapter'].nicks[junction],0)


class NativeTests(unittest.TestCase):
    def test_deterministic_pixels_and_stills(self):
        from .render import Renderer,Style
        r=Renderer(Style(640,360));m=reference()
        a=bytes(r.image(m,4).get_data());r.image(m,9);b=bytes(r.image(m,4).get_data())
        self.assertEqual(a,b)
        with tempfile.TemporaryDirectory()as tmp:
            png=Path(tmp)/'frame.png';svg=Path(tmp)/'frame.svg'
            r.png(png,m,4);r.svg(svg,m,4)
            self.assertTrue(png.read_bytes().startswith(b'\x89PNG'))
            self.assertIn('<svg',svg.read_text())

    def test_native_video_and_failure(self):
        from .render import Renderer,Style
        from .export import export_mp4
        r=Renderer(Style(320,180));m=reference()
        # Shorten presentation for smoke test without altering the molecular event.
        s=m.shots[0];s.narration='Copy.';s.establish=.1;s.action=.3;s.hold=.1;s.narration_seconds=0
        with tempfile.TemporaryDirectory()as tmp:
            path=Path(tmp)/'clip.mp4';meta=export_mp4(m,path,r,fps=5)
            probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-show_streams','-of','json',str(path)]))['streams'][0]
            self.assertEqual((probe['codec_name'],probe['pix_fmt']),('h264','yuv420p'))
            self.assertEqual(int(probe['nb_frames']),meta['frames'])
            self.assertTrue(path.with_suffix('.manifest.json').exists())
            with self.assertRaisesRegex(RuntimeError,'FFmpeg export failed'):
                export_mp4(m,Path(tmp)/'bad.mp4',r,fps=5,ffmpeg='false')
            self.assertFalse((Path(tmp)/'bad.mp4').exists())


if __name__=='__main__':
    suite=unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(MolecularTests))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(ReactionTests))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(PresentationTests))
    if '--native'in sys.argv:suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(NativeTests))
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(not result.wasSuccessful())
