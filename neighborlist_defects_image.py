
##export xyz file for ovito visualization
import pandas as pd
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np
from ovito.data import *
from collections import Counter
import os
from ovito.vis import ParticlesVis
from ovito.modifiers import ComputePropertyModifier
from ovito.vis import ColorLegendOverlay, Viewport
from ovito.vis import OSPRayRenderer, TachyonRenderer, OpenGLRenderer

def setup_particle_types(frame, data):
    types = data.particles_.particle_types_
    #types.type_by_id_(1).name = "Ga"
    types.type_by_id_(1).radius = 0.6
    #types.type_by_id_(2).name = "Ga"
    types.type_by_id_(2).radius = 0.6
    #types.type_by_id_(3).name = "O"
    types.type_by_id_(3).radius = 0.3
    #types.type_by_id_(4).name = "O"
    types.type_by_id_(4).radius = 0.3
    #types.type_by_id_(5).name = "O"
    types.type_by_id_(5).radius = 0.3


#path1='/home/heruhe/Desktop/Ga2o3/cascade/5type_overlapping/out'
path1='/home/heruhe/Desktop/Ga2o3/cascade/5type/1500ev/anneal'
#natom=81920
natom=160000
os.chdir(path1)

for i in range(1):

    #I_atomsdf=pd.read_csv('frame{}/defects{}.csv'.format(i,i))
    I_atomsdf=pd.read_csv('anneal_npt_defects8.csv')
    for a in range(len(I_atomsdf['Particle Identifier'])):
    #for a in [0]:
        #pipeline = import_file("data.lastframe-{}".format(i))
        pipeline = import_file("data.lastframe_anneal_npt_8")
        pipeline.add_to_scene()
        index=I_atomsdf['Particle Identifier'][a]
        print(index)
        #obtain interstitial atom's positon
        x=I_atomsdf[I_atomsdf['Particle Identifier']==index].values[0][3]
        y=I_atomsdf[I_atomsdf['Particle Identifier']==index].values[0][4]
        z=I_atomsdf[I_atomsdf['Particle Identifier']==index].values[0][5]

        hx=x+6
        lx=x-6
        hy=y+6
        ly=y-6
        hz=z+6
        lz=z-6
        
        mod=ExpressionSelectionModifier(expression = 'Position.X<{}&&Position.X>{}&&Position.Y<{}&&Position.Y>{}&&Position.Z<{}&&Position.Z>{}'.format(hx,lx,hy,ly,hz,lz))
        pipeline.modifiers.append(mod)
        pipeline.modifiers.append(setup_particle_types)
        pipeline.modifiers.append(InvertSelectionModifier())
        pipeline.modifiers.append(DeleteSelectedModifier())
        pipeline.modifiers.append(CreateBondsModifier(cutoff=2.2))
        pipeline.modifiers.append(ExpressionSelectionModifier(expression ='ParticleIdentifier=={}'.format(index) ))
        pipeline.modifiers.append(AssignColorModifier(color=(1.15, 1.15, 1.15)))
    
        export_file(pipeline, "npt_anneal_frame/defects{}.xyz".format(index), "xyz", columns =["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])
        cell_vis = pipeline.compute().cell.vis
        cell_vis.enabled = False
        vis_element = pipeline.compute().particles.vis
        vis_element.shape = ParticlesVis.Shape.Sphere
        #vis_element.radius = 0.6
        bonds_vis = pipeline.compute().particles.bonds.vis
        bonds_vis.width = 0.2
        vp = Viewport(type=Viewport.Type.Front)


        vp.zoom_all()
        #vp.render_image(filename='frame{}/{}.png'.format(i,index),size=(600, 600),renderer=OSPRayRenderer())
        vp.render_image(filename='npt_anneal_frame/{}.png'.format(index),size=(600, 600),renderer=OSPRayRenderer())
        pipeline.remove_from_scene()


