import os
import time

import numpy
from matplotlib import pyplot 

from amuse.units import nbody_system
from amuse.units import units,constants,trigo

from amuse.datamodel import Particles

from amuse.community.fi.interface import FiViewer,FiMap

from amuse.io import read_set_from_file

from ubvi.ubvi import rgb_frame

from interactive import go_interactive

#~ import logging
#~ logging.basicConfig(level=logging.DEBUG)

def all_parts(parts):
    result=Particles()
    result.add_particles(parts)
    for p in parts:
        if p.subsystem:
            result.remove_particle(p)
            sub           = all_parts(p.subsystem)
            sub.position += p.position
            sub.velocity += p.velocity
            result.add_particles(sub)
    return result

def _all_parts(parts,i=0):
    result=Particles()
    for p in parts:
        if p.subsystem:
            sub           = All_parts(p.subsystem)
            sub.position += P.position
            sub.velocity += P.velocity
            result.add_particles(sub)
        else:          
            result.add_particle(p)
    return result


def column_density_map(mapper, part):  
    if callable(mapper):
        mapper  = mapper()
  
    p           = mapper.particles.add_particles(part)
    p.weight    = part.mass.value_in(units.amu)
    projected   = mapper.image.pixel_value 
    #~ print part.weight.sum()
    #~ print ">>>", projected.min(),projected.max(),projected.mean(),projected.sum()
    mapper.stop()
    im          = numpy.transpose(projected)  
    return im

def make_rgb(
        mapper, 
        stars, 
        gas,
        vmax    = None,
        ):  
    if len(stars)==0: 
        print "no stars to show!"
        
    a   = numpy.where(stars.temperature.number >0)[0]
    if len(a)==0: 
        print "no stars with temperature > 0...!"

    print "min,max radius:", stars.radius.min(),stars.radius.max()
    print "min,max temp:", stars.temperature.min(),stars.temperature.max()
        
    vmax, image = rgb_frame(
            stars, 
            dryrun          = False,
            vmax            = vmax,
            percentile      = 0.9995,
            multi_psf       = False,
            sourcebands     = "ubvri",
            mapper_factory  = mapper, 
            gas             = gas,
            #image_width     = image_width,
            #image_size      = image_size,
            )
    
    shape       = image.shape
    image       = numpy.transpose(image)

    return image
    
def makemap(
        mode, 
        converter, 
        gas, 
        stars, 
        image_size,
        image_width,
        #image_target, 
        #viewpoint, 
        #horiz_angle, 
        #ratio, 
        #projection,
        vmax    = None,
        ):
    #image_width   = 2*(image_target-viewpoint).length()*trigo.tan(horiz_angle/2)
  
    def mapper():
        mapper    = FiMap(
                converter, 
                use_gl      = False,
                )#,mode="openmp",redirection="none")
        
        #mapper.parameters.minimum_distance      = 1. | units.AU
        mapper.parameters.image_size            = image_size
        #mapper.parameters.image_target          = image_target
        mapper.parameters.image_width           = image_width
        #mapper.parameters.projection_direction  = (image_target-viewpoint)/(image_target-viewpoint).length()
        #mapper.parameters.projection_mode       = projection
        #mapper.parameters.image_angle           = horiz_angle
        #mapper.parameters.viewpoint             = viewpoint
        mapper.parameters.projection_direction  = [0,0,1]
        mapper.parameters.upvector              = [0,-1,0]
        mapper.parameters.extinction_flag       = True if mode=="stars+gas" else False
        return mapper
  
    if mode=="gas":
        image   = column_density_map(
                mapper, 
                gas,
                )
    else:
        image   = make_rgb(
                mapper, 
                stars, 
                gas if mode=="stars+gas" else None,
                vmax  = vmax,
                )
    
    return image

def calculate_effective_temperature(luminosity, radius):
    return ((luminosity/(constants.four_pi_stefan_boltzmann*radius**2))**.25).in_(units.K)


class FRESCO(object):
    helpstring  = """
    """

    def __init__(
            self, 
            stars_file_name     = "stars.hdf5",
            gas_file_name       = "gas.hdf5",
            file_type           = "amuse",
            converter_mass      = 1000. | units.MSun,
            converter_length    = 1. | units.parsec, 
            ):
        self.stars_file_name    = stars_file_name
        self.gas_file_name      = gas_file_name
        self.file_type          = file_type
        self.length_unit        = units.parsec
        if not os.path.isfile(self.stars_file_name):
            raise Exception("%s does not exist"%self.file_name)

        self.converter  = nbody_system.nbody_to_si(
                converter_mass,
                converter_length,
                )

        self.vmax       = None
        
    def load_files(
            self,
            ):

        try:
            self.stars  = read_set_from_file(self.stars_file_name,self.file_type)
            self.mode   = "stars"
        except:
            self.stars  = None
        try:
            self.gas    = read_set_from_file(self.gas_file_name,self.file_type)
            if self.mode == "stars":
                self.mode   = "stars+gas"
            else:
                self.mode   = "gas"
        except:
            self.gas    = None

    def create_image(
            self,
            age = 0|units.yr,
            #projection  = "perspective", 
            ):
        
        converter   = self.converter
        gas         = self.gas
        stars       = self.stars
        mode        = self.mode

        try:
            print self.stars[0].temperature
        except:
            try:
                self.stars.temperature = calculate_effective_temperature(
                        self.stars.luminosity,
                        self.stars.radius,
                        )
            except:
                print "Calculating luminosity/temperature for %s old stars..."%(age)
                from amuse.community.sse.interface import SSE
                se = SSE()
                se.particles.add_particles(self.stars)
                if age > 0|units.Myr:
                    print "Evolve start"
                    se.evolve_model(age)
                    print "Evolve done"
                self.stars.luminosity = se.particles.luminosity
                self.stars.radius = se.particles.radius
                self.stars.temperature = calculate_effective_temperature(
                        self.stars.luminosity,
                        self.stars.radius,
                        )
                se.stop()


        res_x           = 1024
        res_y           = 1024
        ratio           = res_x/res_y
        #image_target    = self.viz.parameters.image_target #[0, 0, 0] length
        #viewpoint       = self.viz.parameters.viewpoint    #[0.0, 1.0, 0.0] length
        #angle           = self.viz.parameters.image_angle  #45.0 deg
        #horiz_angle     = 2*trigo.arctan(ratio*trigo.tan(angle/2))
        image_width     = 10.|self.length_unit#2*(image_target-viewpoint).length()*trigo.tan(horiz_angle/2)
        image_size      = [res_x,res_y]
        
        image           = makemap(
                mode, 
                converter, 
                gas, 
                stars, 
                image_size,
                image_width,
                #image_target, 
                #viewpoint,
                #horiz_angle, 
                #ratio, 
                #projection,
                vmax    = self.vmax,
                )

        iw              = image_width.value_in(self.length_unit)
        self.extent     = [-0.5*iw,0.5*iw,-0.5*iw/ratio,0.5*iw/ratio]

        if mode=="gas":
            mean_density    = numpy.mean(image)
            image           = numpy.log10(image+mean_density/1.e6)
  
            vmax    = numpy.log10(mean_density)+2
            vmin    = numpy.log10(mean_density)-2

        self.image = image


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-i", 
        dest="file_name", default = "", type="str", help="input filename [%default]")
    return result

if __name__ in ('__main__','__plot__'):
    o, arguments    = new_option_parser().parse_args()

    fig = pyplot.figure(
            figsize=(8,8),
            )
    ax  = fig.add_subplot(1,1,1)
    
    fresco  = FRESCO(
            stars_file_name = o.file_name,
            gas_file_name   = None,
            )
    fresco.load_files()
    fresco.create_image()
    
    ax.imshow(
            fresco.image, 
            origin  = "lower",
            extent  = fresco.extent,
            )
    ax.set_xlim(fresco.extent[:2])
    ax.set_ylim(fresco.extent[2:])
    ax.set_xlabel("[%s]"%(fresco.length_unit))
    ax.set_ylabel("[%s]"%(fresco.length_unit))

    pyplot.show()
