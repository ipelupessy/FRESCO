# -*- coding: utf-8 -*-
from __future__ import (
        print_function,
        division
        )

from amuse.units import units, constants
from amuse.io import read_set_from_file
from ubvinew import rgb_frame

# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse


def new_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '-i',
            dest='filename',
            default='test.hdf5',
            help='file containing the stars [test.hdf5]',
            )
    parser.add_argument(
            '-o',
            dest='imagefilename',
            default='test.png',
            help='write image to this file [test.png]',
            )
    parser.add_argument(
            '-b',
            dest='sourcebands',
            default='ubvri',
            help='colour bands to use [ubvri]',
            )
    parser.add_argument(
            '-x',
            dest='plot_axes',
            action='store_true',
            default=False,
            help='plot axes [False]',
            )
    return parser.parse_args()


def calculate_effective_temperature(luminosity, radius):
    return (
            (
                luminosity / (constants.four_pi_stefan_boltzmann*radius**2)
                )**.25
            ).in_(units.K)


def image_from_stars(
        stars,
        image_width=10. | units.parsec,
        image_size=[1024, 1024],
        percentile=0.9995,
        calc_temperature=True,
        age=0. | units.Myr,
        sourcebands="ubvri",
        ):
    if calc_temperature:
        # calculates the temperature of the stars from their total luminosity
        # and radius, calculates those first if needed
        try:
            stars.temperature = calculate_effective_temperature(
                    stars.luminosity,
                    stars.radius,
                    )
        except:
            print(
                    "Calculating luminosity/temperature for %s old stars..."
                    % (age)
                    )
            from amuse.community.sse.interface import SSE
            se = SSE()
            se.particles.add_particles(stars)
            if age > 0 | units.Myr:
                se.evolve_model(age)
            stars.luminosity = se.particles.luminosity
            stars.radius = se.particles.radius
            stars.temperature = calculate_effective_temperature(
                    stars.luminosity,
                    stars.radius,
                    )
            se.stop()

    vmax, rgb = rgb_frame(
            stars,
            dryrun=False,
            image_width=image_width,
            multi_psf=False,  # True,
            image_size=image_size,
            percentile=percentile,
            sourcebands=sourcebands,
            )
    return rgb['pixels']


if __name__ == "__main__":
    args = new_argument_parser()
    filename = args.filename
    imagefilename = args.imagefilename

    plot_axes = args.plot_axes
    sourcebands = args.sourcebands

    length_unit = units.parsec
    dpi = 600
    image_width_arcsec = 160

    # FIXME this should depend on the distance!
    # size = fixed at nr of arcmin
    image_width = 10. | units.parsec

    # FIXME the psf is fixed pixel size, so the pixels in the image here
    # reflects how much pixels will be spread!  In principle, this should be a
    # fixed number, equal to the number of pixels of the ccd for which the psf
    # was made.  Also, the image should reflect the distance to the "observed"
    # cluster.
    image_size = [2048, 2048]
    if plot_axes:
        left = 0.15
        bottom = 0.15
    else:
        left = 0.
        bottom = 0.
    right = 1.0
    top = 1.0
    figwidth = image_size[0] / dpi / (right - left)
    figheight = image_size[1] / dpi / (top - bottom)
    figsize = (figwidth, figheight)

    age = 500. | units.Myr
    percentile = 0.9995  # for determining vmax

    xmin = -0.5 * image_width.value_in(length_unit)
    xmax = 0.5 * image_width.value_in(length_unit)
    ymin = -0.5 * image_width.value_in(length_unit)
    ymax = 0.5 * image_width.value_in(length_unit)

    stars = read_set_from_file(
            filename,
            "amuse",
            close_file=True,
            )

    # FIXME: add these features
    # - Rotate so that xy = observed x/y axes of figure
    # - Scale positions to desired ra/dec (script Alison)
    # - calculate vmax based on nr of photons/exposure time

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, dpi=dpi)
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_xlabel("[%s]" % (length_unit))
    ax.set_ylabel("[%s]" % (length_unit))
    ax.set_aspect(1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    image = image_from_stars(
            stars,
            image_width=image_width,
            image_size=image_size,
            percentile=percentile,
            calc_temperature=True,
            age=age,
            sourcebands=sourcebands,
            )

    plt.imshow(
            image,
            origin='lower',
            extent=[
                xmin,
                xmax,
                ymin,
                ymax,
                ],
            )

    plt.savefig(
            imagefilename,
            dpi=dpi,
            )
