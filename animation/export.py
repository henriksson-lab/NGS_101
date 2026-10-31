"""Native frames streamed to FFmpeg; no browser and no temporary frame directory."""
from __future__ import annotations
from dataclasses import asdict
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from .model import DEFAULT_SPEED


def export_mp4(movie, output, renderer, *, fps=30, wpm=125, speed=DEFAULT_SPEED, ffmpeg='ffmpeg'):
    movie.validate()
    if not isinstance(fps,int) or fps<=0: raise ValueError('fps must be a positive integer')
    if renderer.style.width%2 or renderer.style.height%2:
        raise ValueError('H.264 yuv420p requires even frame dimensions')
    binary=shutil.which(ffmpeg)
    if not binary: raise RuntimeError('MP4 export requires FFmpeg on PATH')
    output=Path(output)
    if output.suffix.lower()!='.mp4': raise ValueError('video output must end in .mp4')
    if output.exists() or output.with_suffix('.manifest.json').exists():raise FileExistsError('output or manifest already exists; choose a new path')
    output.parent.mkdir(parents=True,exist_ok=True)
    duration=movie.duration(wpm,speed); frames=math.ceil(duration*fps)
    # Cairo RGB24 uses native-endian32bit pixels; stream without a Python pixel loop.
    pixel_format='bgr0'if sys.byteorder=='little'else'0rgb'
    with tempfile.TemporaryDirectory(prefix='.animation-',dir=output.parent)as tmp:
        partial=Path(tmp)/'movie.mp4'
        args=[binary,'-hide_banner','-loglevel','error','-y','-f','rawvideo','-pixel_format',pixel_format,
              '-video_size',f'{renderer.style.width}x{renderer.style.height}','-framerate',str(fps),
              '-i','pipe:0','-an','-c:v','libx264','-preset','fast','-crf','18',
              '-pix_fmt','yuv420p','-movflags','+faststart',str(partial)]
        with tempfile.TemporaryFile()as err:
            proc=subprocess.Popen(args,stdin=subprocess.PIPE,stderr=err)
            try:
                for i in range(frames):
                    surface=renderer.image(movie,i/fps,wpm,speed)
                    if surface.get_stride()!=renderer.style.width*4: raise RuntimeError('unsupported frame stride')
                    proc.stdin.write(surface.get_data())
                proc.stdin.close(); code=proc.wait()
            except BaseException as exc:
                try:proc.stdin.close()
                except OSError:pass
                try:proc.wait(timeout=15)
                except subprocess.TimeoutExpired:proc.kill();proc.wait()
                err.seek(0)
                message=err.read().decode(errors='replace')
                if proc.returncode and isinstance(exc,Exception): raise RuntimeError('FFmpeg export failed: '+message)from exc
                raise
            if code:
                err.seek(0);raise RuntimeError('FFmpeg export failed: '+err.read().decode(errors='replace'))
        partial.replace(output)
    manifest=movie.to_dict()
    manifest['export']={'fps':fps,'frames':frames,'duration_seconds':frames/fps,'wpm':wpm,'speed':speed,
                        'style':asdict(renderer.style),'codec':'H.264','pixel_format':'yuv420p',
                        'timing_note':'Illustrative movement; not experimental kinetics. Word counts estimate speech, not forced audio alignment.'}
    output.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest['export']
