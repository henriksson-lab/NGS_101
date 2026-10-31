"""python3 -m animation --help"""
from __future__ import annotations
import argparse
from dataclasses import replace
import json
from pathlib import Path
from .model import Movie, DEFAULT_SPEED, finite
from .examples import EXAMPLES


def main():
    p=argparse.ArgumentParser(description='Native chemistry presentation animation; no browser.')
    p.add_argument('example',choices=[*EXAMPLES,'json','list'])
    p.add_argument('--input',type=Path,help='Versioned movie JSON when example=json')
    p.add_argument('--output',type=Path,help='.mp4, .svg, .png or .json')
    p.add_argument('--list-shots',action='store_true')
    p.add_argument('--timings',type=Path,help='JSON mapping shot ids to timing/narration overrides')
    p.add_argument('--font',default='DejaVu Sans')
    p.add_argument('--shot',help='Export one shot by id')
    p.add_argument('--time',type=float,default=4.,help='Timestamp for SVG/PNG stills')
    p.add_argument('--size',default='1920x1080')
    p.add_argument('--fps',type=int,default=30)
    p.add_argument('--wpm',type=float,default=125,help='Narration word-rate estimate')
    p.add_argument('--speed',type=float,default=DEFAULT_SPEED,help='Playback multiplier (default: 1.5; use 1 for original timing), including narration allowance')
    p.add_argument('--duration',type=float,help='Override selected shot duration; rejects truncated narration/holds')
    p.add_argument('--narration-seconds',type=float,help='Known narration length for selected shot')
    p.add_argument('--bases',action='store_true',default=True,help='Nucleotide display (default)')
    p.add_argument('--directions',action='store_true',default=True,help='Strand arrows (default)')
    for flag in ['evidence','coordinates','reagents','timeline']:
        p.add_argument('--'+flag,action='store_true')
    a=p.parse_args()
    finite(a.speed,'speed')
    if a.speed<=0:raise ValueError('speed must be positive')
    if a.example=='list':
        for key,fn in EXAMPLES.items():
            m=fn();print(f'{key:14s} {m.title:24s} {m.duration(a.wpm,a.speed):6.1f}s   '+', '.join(s.id for s in m.shots))
        return
    if a.example=='json':
        if a.input is None:p.error('--input required for json')
        movie=Movie.from_dict(json.loads(a.input.read_text()))
    else:movie=EXAMPLES[a.example]()
    if a.timings:
        edits=json.loads(a.timings.read_text());allowed={'duration','narration_seconds','establish','action','hold','narration','cue'}
        if set(edits)-{s.id for s in movie.shots}:p.error('timing file names unknown shot')
        for edit in edits.values():
            if set(edit)-allowed:p.error('unsupported timing override field')
        movie=replace(movie,shots=[replace(s,**edits.get(s.id,{}))for s in movie.shots])
    if a.list_shots:
        for s in movie.shots:print(s.id,round(s.timing(a.wpm)[2]/a.speed,2),s.narration)
        return
    if a.output is None:p.error('--output required for export')
    if a.shot:
        shots=[s for s in movie.shots if s.id==a.shot]
        if not shots:p.error('unknown shot id')
        movie=replace(movie,shots=shots)
    if a.duration is not None or a.narration_seconds is not None:
        if len(movie.shots)!=1:p.error('timing overrides require --shot or a one-shot movie')
        shot=movie.shots[0]
        movie=replace(movie,shots=[replace(shot,duration=a.duration,narration_seconds=a.narration_seconds)])
    movie.validate();movie.duration(a.wpm,a.speed)
    output=a.output
    if output.exists():p.error('output already exists; choose a new path or remove it explicitly')
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.suffix=='.json':output.write_text(json.dumps(movie.to_dict(),indent=2)+'\n');return
    from .render import Renderer,Style
    try:w,h=map(int,a.size.split('x'))
    except ValueError:p.error('--size must be WIDTHxHEIGHT')
    renderer=Renderer(Style(w,h,font=a.font,**{k:getattr(a,k)for k in ['bases','directions','evidence','coordinates','reagents','timeline']}))
    if output.suffix=='.png':renderer.png(output,movie,a.time,a.wpm,a.speed)
    elif output.suffix=='.svg':renderer.svg(output,movie,a.time,a.wpm,a.speed)
    elif output.suffix=='.mp4':
        from .export import export_mp4
        result=export_mp4(movie,output,renderer,fps=a.fps,wpm=a.wpm,speed=a.speed)
        print(json.dumps(result))
    else:p.error('output must be .mp4/.png/.svg/.json')


if __name__=='__main__':
    try:main()
    except (ValueError,RuntimeError,OSError,TypeError,KeyError) as exc:
        raise SystemExit('Animation error: '+str(exc)) from None
