"""Decode native MP4 review frames and verify them against current definitions.

Creates an audit record pending human visual inspection, never declares its own
pictures good. Does not fetch sources or render video using a browser.
"""
from pathlib import Path
import argparse,hashlib,json,math,subprocess
from .examples import EXAMPLES
from .batch import probe
from .render import Renderer,Style
from .model import Movie


def run(directory,names):
    root=Path(directory);records={}
    for name in names:
        candidates=list(root.rglob(name+'.mp4'))
        if len(candidates)!=1:raise ValueError(f'{name}: expected one full MP4, got {len(candidates)}')
        video=candidates[0];manifest=json.loads(video.with_suffix('.manifest.json').read_text());export=manifest.pop('export')
        manifest=json.loads(json.dumps(Movie.from_dict(manifest).to_dict()))
        m=EXAMPLES[name]();current=json.loads(json.dumps(m.to_dict()))
        if manifest!=current:raise ValueError(name+': exported scene differs from current definition')
        if export['style']['timeline']:raise ValueError('presentation footer should be opt-in')
        info=probe(video,export['style']['width'],export['style']['height'],export['fps'],export['frames'])
        reviewdir=video.parent/(name+'-review');reviewdir.mkdir(exist_ok=True)
        frames=[];offset=0;speed=export['speed'];renderer=Renderer(Style(1920,1080))
        for shot in m.shots:
            samples=[('action35',shot.establish+shot.action*.35),('action85',shot.establish+shot.action*.85)]
            if shot.transition:samples.insert(0,('transition',shot.transition*.5))
            for label,local in samples:
                timestamp=(offset+local)/speed;path=reviewdir/(shot.id+'-'+label+'.png')
                subprocess.run(['ffmpeg','-v','error','-y','-ss',str(timestamp),'-i',str(video),'-frames:v','1','-vf','scale=960:540',str(path)],check=True)
                frames.append(dict(shot=shot.id,phase=label,seconds=timestamp,file=path.name))
            offset+=shot.timing()[2]
        cairo=renderer.cairo;columns=3;width=640;height=360
        sheet=cairo.ImageSurface(cairo.FORMAT_RGB24,columns*width,math.ceil(len(frames)/columns)*height);c=cairo.Context(sheet)
        for i,item in enumerate(frames):
            surface=cairo.ImageSurface.create_from_png(str(reviewdir/item['file']));c.save();c.translate((i%columns)*width,(i//columns)*height);c.scale(width/960,height/540);c.set_source_surface(surface);c.paint();c.restore()
        sheet.write_to_png(str(reviewdir/'decoded-contact.png'))
        # Separate exportable vector figure; this is not proof of MP4 content.
        renderer.svg(reviewdir/'final.svg',m,max(0,m.duration()-1))
        files=list(Path(__file__).parent.glob('*.py'))
        if m.provenance.get('protocol_module'):files.append(Path(m.provenance['protocol_module']))
        record=dict(example=name,video=str(video),movie_sha256=hashlib.sha256(video.read_bytes()).hexdigest(),source_scene_match=True,export=export,video_info=info,shots=len(m.shots),frames=frames,limits=m.provenance.get('limits',[]),visual_review='pending',code_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in files})
        (reviewdir/'audit.json').write_text(json.dumps(record,indent=2)+'\n');records[name]=record
        print(name,'decoded',len(frames),'frames',flush=True)
    (root/'audit.json').write_text(json.dumps(records,indent=2)+'\n')
    return records


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,required=True);p.add_argument('--examples',nargs='+',choices=EXAMPLES,required=True);a=p.parse_args();run(a.directory,a.examples)

if __name__=='__main__':main()
