import math
from ovito.vis import *
import pandas as pd
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np
from ovito.data import *
import os


path1='/home/heruhe/Desktop/Ga2o3/cascade/5type/1500ev'
os.chdir(path1) 
pipeline = import_file("ws0")
pipeline.add_to_scene()


vp = Viewport()
vp.type = Viewport.Type.Perspective
vp.camera_pos = (-100, -150, 150)
vp.camera_dir = (2, 3, -3)
vp.fov = math.radians(60.0)
vp.render_image(filename='0.png', size=(320, 240),renderer=TachyonRenderer())
# Trajectory rendering loop over all frames of the loaded animation:
#for frame in range(pipeline.source.num_frames):
    # Render the current frame.
#    print("Rendering frame", frame)
#    img = vp.render_image(size = (400,400), frame=frame)
    # Compute alpha value for current frame.
#    alpha = 1.0 - frame / pipeline.source.num_frames
    # Replace alpha channel of rendered image.
#    alpha_img = QImage(img.size(), QImage.Format_Alpha8)
#    alpha_img.fill(QColor.fromRgbF(1,1,1,alpha))
#    img.setAlphaChannel(alpha_img)
    # Write rendered image to an output file. 
#    img.save('frame%i.png' % frame)