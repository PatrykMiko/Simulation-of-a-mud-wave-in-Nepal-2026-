import csv
import subprocess
import numpy as np
from PIL import Image,ImageDraw
from ravaflow_tools import metrics as M,viz


def render(folder,data):
    z,meta=M.read_raster(data/'elevation.tif')
    ref,rm=M.read_raster(data/'reference_extent.tif')
    rel,_=M.read_raster(data/'release_height.tif')
    h,hm=M.read_raster(folder/'hflow_max.tif')
    obs,valid=M.reference_masks(ref,rm)
    viz.plot_extent_comparison(z,M.model_extent(h,hm,.1),obs,valid&(rel==0),meta,
                               folder/'zasieg.png','Maksymalny zasięg ≥ 0,1 m • siatka 30 m')
    viz.plot_hmax(z,h,hm,folder/'glebokosc.png',title='Maksymalna głębokość [m]')
    base=(100+120*viz.hillshade(z,30,meta['nodata'])).astype('uint8')
    # One depth scale and an undistorted native aspect ratio across all frames.
    from matplotlib import colormaps
    palette = (255 * colormaps['viridis'](np.linspace(0, 1, 256))[:,:3]).astype('uint8')
    maximum = max(1., float(np.max(h)))
    width = 864
    size = (width, 2 * round(width * z.shape[0] / (2*z.shape[1])))
    def depth_indices(depth):
        return np.clip(255*np.log10(depth/.1)/np.log10(maximum/.1),0,255).astype('uint8')

    movie=folder/'fala.tmp.mp4'
    cmd=['ffmpeg','-loglevel','error','-y','-f','rawvideo','-pixel_format','rgb24',
         '-video_size',f'{size[0]}x{size[1]}','-framerate','5','-i','-','-c:v','libx264','-threads','1','-pix_fmt','yuv420p',str(movie)]
    dtype=np.dtype([('index','<i4'),('depth','<f4'),('solid_v','<f4'),('water_v','<f4')])
    with subprocess.Popen(cmd,stdin=subprocess.PIPE) as proc:
        with (folder/'frames/timeline.csv').open() as f:
            for frame,t,*_ in csv.reader(f):
                records=np.fromfile(folder/'frames'/f'{int(frame):05}.bin',dtype=dtype)
                indices=records['index']
                if (indices<0).any() or (indices>=z.size).any():raise ValueError('Poza siatką')
                rgb=np.repeat(base[:,:,None],3,axis=2)
                rgb.reshape(-1,3)[indices]=palette[depth_indices(records['depth'])]
                image=Image.fromarray(rgb).resize(size)
                draw=ImageDraw.Draw(image)
                draw.rectangle((4,4,342,52),fill='white')
                draw.text((10,6),f'r.avaflow  |  t = {float(t):.0f} s  |  depth [m]',fill='black')
                for j in range(256):
                    draw.line((10+j,24,10+j,31),fill=tuple(palette[j]))
                draw.text((10,33),'0.1',fill='black')
                draw.text((215,33),f'{maximum:.1f} m',fill='black')
                draw.text((277,20),'log10',fill='black')
                proc.stdin.write(image.tobytes())
        proc.stdin.close()
        if proc.wait():raise RuntimeError('Błąd eksportu MP4')
    movie.replace(folder/'fala.mp4')
