

from ovito.io import *
from ovito.data import *
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np
import os
path='/home/heruhe/Desktop/Ga2o3/cascade/5type/1500ev/dump/pka2'
path='/home/ad/home/h/heruhe/ga2o3/cascade/1500ev/dump'
os.chdir(path) 
for i in range(80,90):
    fn='dump.casc{}'.format(i)
    if os.path.exists(fn):
        pipeline = import_file(fn)

# Perform Wigner-Seitz analysis:
        ws = WignerSeitzAnalysisModifier(
            per_type_occupancies = True,
            affine_mapping = ReferenceConfigurationModifier.AffineMapping.ToReference)
        pipeline.modifiers.append(ws)

# Define a modifier function that selects sites of type A=1 which
# are occupied by exactly one atom of type B=2.
        def modify(frame, data):

    # Retrieve the two-dimensional array with the site occupancy numbers.
    # Use [...] to cast it to a Numpy array.
            occupancies = data.particles['Occupancy']

    # Get the site types as additional input:
            site_type = data.particles['Particle Type']

    # Calculate total occupancy of every site:
            total_occupancy = np.sum(occupancies, axis=1)

    # Set up a particle selection by creating the Selection property:
            selection = data.particles_.create_property('Selection')

    # Select O vacnacies.
            selection[...] = (site_type >2 )  & (total_occupancy == 0)

    # Additionally output the total number of O & Ga vacancies as a global attribute:
            data.attributes['O_vancancies'] = np.count_nonzero(selection)
            data.attributes['Ga_vancancies'] = data.attributes['WignerSeitz.vacancy_count']-np.count_nonzero(selection)
# Insert Python modifier into the data pipeline.
        pipeline.modifiers.append(modify)

# Let OVITO do the computation and export the number of identified
# antisites as a function of simulation time to a text file:
        export_file(pipeline, "vacancies{}.txt".format(i), "txt/attr",
            columns = ['Timestep', 'O_vancancies','Ga_vancancies'],
            multiple_frames = True)
        os.remove(fn)

