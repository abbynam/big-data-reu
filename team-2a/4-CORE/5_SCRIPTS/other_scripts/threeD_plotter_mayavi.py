#!/usr/bin/python
"""
Description:
  Python script that reads in DAT file from a core reconstruction (typically named: output.dat)
    and can: plot a 3D image of the reconstructed data (using mayavi library and Intensity3D.py script)
             create a series of rotated 3D images which can be used to create an animation

  There are 3 plotting options (uncomment the one you want to use, also they can be combined):
    1. Create 3D contour plot (contours sets the number of drawn contours)
        - mlab.contour3d(normData, contours = 10, transparent = True)
    2. Create 3D volume rendering (vmax and vmin set the upper and lower percentiles of visible data)
        - mlab.pipeline.volume(mlab.pipeline.scalar_field(normData), vmin = 0.2, vmax = 0.8)
    3. Create 3D volume with three 2D planes through the location of the maximum value
        - mlab.volume_slice(normData, plane_orientation = 'x_axes', slice_index = max_index[0])
        - mlab.volume_slice(normData, plane_orientation = 'y_axes', slice_index = max_index[1])
        - mlab.volume_slice(normData, plane_orientation = 'z_axes', slice_index = max_index[2])

  To Run: python3 threeD_plotter_mayavi.py output.dat (plots diretory)
    Settings are handled through Global Variable under SETTINGS
      - ANIMATION                   = True (create series of images for animation) / False (open interactive plot)
      - OUTPUT_FILE                 = name of output file, will save in [output path] directory
      - IMAGE_SIZE                  = size of image in pixels, will be square, so assigning both height and width / default is: 400 px
      - ANIMATION_DELAY             = time step of animation in milliseconds

NOTES: 
- written using python v3.6.2 / has not been tested using python2
- Intensity3D.py required to run this code
- Mayavi: 3D scientific data visualization python library requied to run this code
  - link: https://docs.enthought.com/mayavi/mayavi/
- Unfortunately, animation images do not save in [output path].  Look under Documents/mayavi_movies
- I used GIMP to create an animated GIF from the series of images
  - you can find instructions here: http://elearnhub.org/how-to-make-an-animated-gif-in-gimp/

Intensity3D code written by Dennis Mackin <dsmackin@mdanderson.org>
"""
__author__ = "Steve Peterson <steve.peterson@uct.ac.za>"
__date__ = "September 06, 2018"
__version__ = "$Revision: 1.0.0$"

#------------------------------------------------------------------
# PYTHON IMPORT STATEMENTS
#------------------------------------------------------------------

import sys, os
import numpy
from mayavi import mlab

import Intensity3D


#------------------------------------------------------------------
# SETTINGS / DECLARATIONS
#------------------------------------------------------------------

ANIMATION = True
OUTPUT_FILE = 'temp'
IMAGE_SIZE = 400
ANIMATION_DELAY = 100


USAGE = '''
USAGE: %s [Path to 3D dose file] [output path]

The following example generates plots using the image intensity information from the file output.dat.

   python %s output.dat outputfolder
'''



#------------------------------------------------------------------
# FUNCTION/CLASS DEFINITIONS
#------------------------------------------------------------------

@mlab.animate(delay = ANIMATION_DELAY)
def anim(deg = 10):
    f = mlab.gcf()
    print ('\n  Looping endlessly through 3D dose intensity')
    print ('   - rotation: {} degrees in {} ms increments'.format(deg, ANIMATION_DELAY))
    print ('   - use Animation GUI to control animation')
    print ('   - NOTE: no images are being saved!')

    while 1:
        f.scene.camera.azimuth(10)
        f.scene.render()
        yield


@mlab.animate(delay = ANIMATION_DELAY)
def anim_step(N = 72, deg = 5):
    f = mlab.gcf()
    f.scene.movie_maker.record = True
    #f.scene.foreground = (0, 0, 0)
    #f.scene.background = (1, 1, 1)
    print ('\n  Creating series of 3D dose intensity images (format: PNG)')
    print ('   - rotation: {} steps of {} degrees in {} ms increments'.format(N, deg, ANIMATION_DELAY))
    print ('   - NOTE: images are being saved!  But not to {} directory!'.format(output_folder_path))
    
    # Looping through steps to creation animation
    i = 0
    while i in range(N):
        f.scene.camera.azimuth(deg)
        f.scene.render()
        i += 1
        yield



#------------------------------------------------------------------
# MAIN PROGRAM
#------------------------------------------------------------------
if __name__ == "__main__":

    def print_usage():
        print (USAGE % (sys.argv[0], sys.argv[0]))
        sys.exit(100)

    if len(sys.argv) == 3:

        data_file_path = sys.argv[1]
        output_folder_path = sys.argv[2]
    else:
        print_usage()

    # Storing data into Intensity3D
    my3DObj = Intensity3D.Intensity3D(data_file_path, output_folder_path)
    print ('\n  Plotting 3D dose intensity from file: {}'.format(data_file_path))
    print ('   - size: {}'.format(numpy.shape(my3DObj.values)))
    
    # Determining the coordinate limits from the data
    X = my3DObj.bin_centers[0]
    xmin = my3DObj.edges[0][0]
    xmax = my3DObj.edges[0][numpy.shape(my3DObj.edges)[1] - 1]
    Y = my3DObj.bin_centers[1]
    ymin = my3DObj.edges[1][0]
    ymax = my3DObj.edges[1][numpy.shape(my3DObj.edges)[1] - 1]
    Z = my3DObj.bin_centers[2]
    zmin = my3DObj.edges[2][0]
    zmax = my3DObj.edges[2][numpy.shape(my3DObj.edges)[1] - 1]
    print ('   - coordinate limits (mm) -> Xmin: {}, Xmax: {}, Ymin: {}, Ymax: {}, Zmin: {}, Zmax: {}'.format(xmin, xmax, ymin, ymax, zmin, zmax))

    # Normalizing the data
    rawData = my3DObj.values
    maxValue = numpy.amax(rawData)
    normData = rawData / maxValue
    print ('   - normalizing data to 1 / Maximum value: {}'.format(maxValue))
    
    # Finding position of maximum valus
    max_pos = my3DObj.get_maximum()
    print ('   - position of maximum point in reconstruction (mm): {}'.format(max_pos))
    max_value = my3DObj.f3D(max_pos)
    print ('   - value of maximum point in reconstruction (mm): {}'.format(max_value))
    max_index = numpy.unravel_index(my3DObj.values.argmax(), my3DObj.values.shape)
    print ('   - index of maximum point in reconstruction (mm): {}'.format(max_index))
    
    # Create figure with white background and black foreground with size: IMAGE_SIZE X IMAGE_SIZE
    mlab.figure(bgcolor = (1., 1., 1.), fgcolor = (0., 0., 0.), size = (IMAGE_SIZE, IMAGE_SIZE))


    ### THREE PLOTTING OPTIONS (Uncomment appropriate lines, possible to combine)
    # 1. Create 3D contour plot (contours sets the number of drawn contours)
    mlab.contour3d(normData, contours = 10, transparent = True)
    # 2. Create 3D volume rendering (vmax and vmin set the upper and lower percentiles of visible data)
    #mlab.pipeline.volume(mlab.pipeline.scalar_field(normData), vmin = 0.2, vmax = 0.8)
    # 3. Create 3D volume with three 2D planes through the location of the maximum value
    #mlab.volume_slice(normData, plane_orientation = 'x_axes', slice_index = max_index[0])
    #mlab.volume_slice(normData, plane_orientation = 'y_axes', slice_index = max_index[1])
    #mlab.volume_slice(normData, plane_orientation = 'z_axes', slice_index = max_index[2])


    # Add axes to plot
    mlab.axes(ranges=[xmin, xmax, ymin, ymax, zmin, zmax])
    mlab.outline()
    
    # Create animation of the rotating plot
    #  NOTE: There are two options -> anim() & anim_step()
    #   anim() -> rotates image endlessly, controlled by GUI, no images are saved
    #   anim_step() -> rotates image fixed number of steps and saves each step as a PNG image
    if ANIMATION:
        #   format: anim (rotation_step_in_deg) / default values are (10)
        anim(5)             # Loops through animation indefinitely
        #   format: anim_step (number_of_steps, rotation_step_in_deg) / default values are (72, 5)
        #anim_step(10, 36)  # Loops animation through set number of images

    # Save last image
    img_path = output_folder_path+"/"+OUTPUT_FILE+".png"
    mlab.savefig(img_path)
    print ('\n  Writing image to {} -> Image size: {} x {} px'.format(img_path, IMAGE_SIZE, IMAGE_SIZE))
    
    mlab.show()

