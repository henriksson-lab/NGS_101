"""Render and verify the agreed presentation example suite without a browser."""
from __future__ import annotations
import argparse
from dataclasses import replace,asdict
import json
import math
from pathlib import Path
import shutil
import subprocess
import time
from .examples import EXAMPLES
from .export import export_mp4
from .render import Renderer,Style
from .model import Movie, DEFAULT_SPEED

FOCUS={'reference':'copy','smartseq':'rt','crispr-mip':'gap','boloramis':'rca',
       'atac':'tagment','smallseq':'lig3','gibson':'resect','hic':'contact',
       'dropseq':'capture','scrrbs':'conversion','celseq':'ivt','smartseq2':'rt','rnaseh':'cleave','golden-gate':'seal'}


def probe(path,width,height,fps,frames=None):
    result=json.loads(subprocess.check_output(['ffprobe','-v','error','-show_streams','-show_format','-of','json',str(path)]))
    stream=result['streams'][0]
    a,b=map(int,stream['r_frame_rate'].split('/'))
    if (stream['codec_name'],stream['pix_fmt'],stream['width'],stream['height'],a/b)!=('h264','yuv420p',width,height,fps):
        raise ValueError('export codec/dimensions/rate mismatch')
    if frames is not None and int(stream['nb_frames'])!=frames:raise ValueError('export frame count mismatch')
    return {'codec':stream['codec_name'],'width':width,'height':height,'fps':fps,'frames':int(stream['nb_frames']),
            'duration':float(result['format']['duration']),'bytes':Path(path).stat().st_size}


def contact_sheet(movie,renderer,path):
    cairo=renderer.cairo
    width=1440;row_height=270
    sheet=cairo.ImageSurface(cairo.FORMAT_RGB24,width,row_height*len(movie.shots));c=cairo.Context(sheet)
    offset=0;overflows=[]
    for row,shot in enumerate(movie.shots):
        for col,f in enumerate([0,.5,1]):
            surface=renderer.image(movie,offset+shot.establish+shot.action*f, speed=1)
            overflows.extend(renderer.overflows)
            c.save();c.translate(col*480,row*row_height)
            c.scale(480/renderer.style.width,270/renderer.style.height)
            c.set_source_surface(surface);c.paint();c.restore()
        offset+=shot.timing()[2]
    sheet.write_to_png(str(path))
    if overflows:raise ValueError('Text outside layout bounds: '+str(sorted(set(overflows))))


def run(directory,names,size='1920x1080',fps=30,clips=True,resume=False,speed=DEFAULT_SPEED):
    if not shutil.which('ffprobe'):raise RuntimeError('Suite verification needs ffprobe')
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    w,h=map(int,size.split('x'));renderer=Renderer(Style(w,h))
    record=directory/'verification.json'
    results=json.loads(record.read_text())if record.exists()else{}
    for name in names:
        started=time.monotonic();movie=EXAMPLES[name]();movie.validate()
        output=directory/(name+'.mp4');serialized=movie.to_dict()
        (directory/(name+'.json')).write_text(json.dumps(serialized,indent=2)+'\n')
        contact_sheet(movie,renderer,directory/(name+'-contact.png'))
        reused=resume and output.exists()
        if reused:
            manifest=json.loads(output.with_suffix('.manifest.json').read_text());meta=manifest.pop('export')
            if json.loads(json.dumps(Movie.from_dict(manifest).to_dict()))!=json.loads(json.dumps(serialized)):raise ValueError('Resume scene differs from existing export: '+name)
            if meta['style']!=asdict(renderer.style) or meta['fps']!=fps or meta['speed']!=speed:raise ValueError('Resume settings mismatch')
        else:meta=export_mp4(movie,output,renderer,fps=fps,speed=speed)
        video=probe(output,w,h,fps,meta['frames'])
        offset=0
        for shot in movie.shots:
            if shot.id==FOCUS[name]:break
            offset+=shot.timing()[2]
        still_t=(offset+shot.establish+shot.action*.55)/speed
        renderer.png(directory/(name+'.png'),movie,still_t,speed=speed)
        renderer.svg(directory/(name+'.svg'),movie,still_t,speed=speed)
        # Decode real exported video as well as inspecting native stills.
        decoded=directory/(name+'-decoded.png')
        subprocess.run(['ffmpeg','-hide_banner','-loglevel','error','-y','-ss',str(still_t),'-i',str(output),
                        '-frames:v','1',str(decoded)],check=True)
        entry={'video':video,'operation_seconds':time.monotonic()-started,'reused_video':reused,'shots':len(movie.shots),
               'selected_still_time':still_t,'scope':movie.provenance.get('representation')}
        if clips:
            focus=directory/(name+'-focus.mp4');duration=shot.timing()[2]/speed
            subprocess.run(['ffmpeg','-hide_banner','-loglevel','error','-y','-ss',str(offset/speed),'-i',str(output),
                            '-t',str(duration),'-an','-c:v','libx264','-crf','18','-preset','fast',
                            '-pix_fmt','yuv420p','-movflags','+faststart',str(focus)],check=True)
            entry['focus']=probe(focus,w,h,fps)
            entry['focus']['shot']=shot.id
            focus.with_suffix('.manifest.json').write_text(json.dumps({'parent':output.name,'parent_source_manifest':output.with_suffix('.manifest.json').name,
                'source_interval':[offset/speed,offset/speed+duration],'shot':shot.id,'narration':shot.narration,
                'note':'Frame-accurate clip from the native-rendered parent MP4; no browser capture.'},indent=2)+'\n')
        if not clips and name in results and 'focus'in results[name]:entry['focus']=results[name]['focus']
        results[name]=entry
        (directory/'verification.json').write_text(json.dumps(results,indent=2)+'\n')
        print(name,json.dumps(entry),flush=True)
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--examples',nargs='+',choices=EXAMPLES,default=list(EXAMPLES))
    p.add_argument('--size',default='1920x1080');p.add_argument('--fps',type=int,default=30)
    p.add_argument('--no-clips',action='store_true');p.add_argument('--resume',action='store_true')
    p.add_argument('--speed',type=float,default=DEFAULT_SPEED,help='Playback multiplier; default 1.5')
    a=p.parse_args();run(a.directory,a.examples,a.size,a.fps,not a.no_clips,a.resume,a.speed)


if __name__=='__main__':main()
