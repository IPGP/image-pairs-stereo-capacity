#!/usr/bin/python3
# -*- coding:utf-8 -*-

import sys
import argparse
import textwrap
from math import radians, degrees, sin, cos, acos, tan, atan2, pi
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import combinations

# howstereo.py - Computes the B to H ratio of pairs of Pleiades or SPOT6|7
# images
# Copyright (C) 2022 Arthur Delorme
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>.

# (Contact: delorme@ipgp.fr)

class Image(object):
    """
    Regroup the parameters attached to an image
    
    A is a point on the look direction with OA = 1
    
    Input:
        scan        [float] incidence angle in the Scan axis direction
        ortho       [float] incidence angle in the OrthoScan axis direction
        az          [float] azimuth of the Scan axis
        s_comp      [float] coordinate of A on the Scan axis
        o_comp      [float] coordinate of A on the OrthoScan axis
        z_comp      [float] coordinate of A on the z axis
        geo_n_comp  [float] coordinate of A on the North axis
        geo_w_comp  [float] coordinate of A on the West axis
        geo_z_comp  [float] coordinate of A on the z axis
    """
    
    def __init__(self, name, scan, ortho, az,
            s_comp=None, o_comp=None, z_comp=None,
            geo_n_comp=None, geo_w_comp=None, geo_z_comp=None):
        self._name = name
        self._scan = scan
        self._ortho = ortho
        self._az = az
        self._s_comp = s_comp
        self._o_comp = o_comp
        self._z_comp = z_comp
        self._geo_n_comp = geo_n_comp
        self._geo_w_comp = geo_w_comp
        self._geo_z_comp = geo_z_comp
    
    @property
    def name(self):
        return self._name
    
    @property
    def scan(self):
        return self._scan
    
    @property
    def ortho(self):
        return self._ortho
    
    @property
    def az(self):
        return self._az
    
    @property
    def s_comp(self):
        return self._s_comp
    @s_comp.setter
    def s_comp(self, value):
        self._s_comp = value
    
    @property
    def o_comp(self):
        return self._o_comp
    @o_comp.setter
    def o_comp(self, value):
        self._o_comp = value
    
    @property
    def z_comp(self):
        return self._z_comp
    @z_comp.setter
    def z_comp(self, value):
        self._z_comp = value
    
    @property
    def geo_n_comp(self):
        if not self._geo_n_comp:
            self.compute_geo_coord()
        return self._geo_n_comp
    @geo_n_comp.setter
    def geo_n_comp(self, value):
        self._geo_n_comp = value
    
    @property
    def geo_w_comp(self):
        if not self._geo_n_comp:
            self.compute_geo_coord()
        return self._geo_w_comp
    @geo_w_comp.setter
    def geo_w_comp(self, value):
        self._geo_w_comp = value
    
    @property
    def geo_z_comp(self):
        if not self._geo_n_comp:
            self.compute_geo_coord()
        return self._geo_z_comp
    @geo_z_comp.setter
    def geo_z_comp(self, value):
        self._geo_z_comp = value
    
    def compute_geo_coord(self):
        """
        Computes the geographic coordinates of point A, if the required
        variables (i.e. self._s_comp, self._o_comp, self._z_comp) have been set
        """
        
        if self._s_comp is not None and self._o_comp is not None and \
                self._z_comp is not None:
            # We express A coordinates in the geographic reference through a
            # rotation of self._az around the z axis
            rot_z = [[cos(self._az), -sin(self._az),       0],
                     [sin(self._az),  cos(self._az),       0],
                     [            0,              0,       1]]
            
            sym_z = [[1,  0,  0], # We also need to apply a symmetry with
                     [0, -1,  0], # respect to the (N,Z) plane, because the
                     [0,  0,  1]] # (O, N, E, Z) coordinate system is not direct
            
            sym_rot_z = np.dot(rot_z, sym_z)
            
            self._geo_n_comp, self._geo_w_comp, self._geo_z_comp = np.dot(
                    sym_rot_z,
                    [self._s_comp, self._o_comp, self._z_comp]
                )

def compute_b_to_h(im1, im2):
    """
    Computes the B/H ratio of a couple of images
    
    See the Figures 45 and 47 in the Pléiades Imagery User Guide.
    See also the Figures 35, 36 and 37 in the SPOT 6 & SPOT 7 Imagery User
    Guide.
    The coordinate system is (O, Scan axis, OrthoScan axis, perpendicular to
    the ground) -> (O, Scan, Ortho, z).
    P is a point on the look direction of the first acquisition with OP = 1.
    
    s_comp          coordinate of P (or Q or P') on the Scan axis
    o_comp          coordinate of P (or Q or P') on the OrthoScan axis
    z_comp          coordinate of P (or Q or P') on the z axis
    
    Input:
        im1         Image instance for the first image
        im2         Image instance for the second image
    
    Returns:
        stereo_angle (in degrees)
        b_to_h
    """
    
    # 1) We express P coordinates in the reference of the first acquisition
    
    im1.s_comp = cos(im1.ortho) * sin(im1.scan)
    im1.o_comp = cos(im1.scan) * sin(im1.ortho)
    im1.z_comp = cos(im1.ortho) * cos(im1.scan)
    
    # 2) We express P coordinates in the reference of the second acquisition
    #    with a rotation around the z axis
    
    teta = im2.az - im1.az # Difference of azimuth between the two acquisitions
    
    p = [im1.s_comp, im1.o_comp, im1.z_comp]
    rot_z = [[cos(teta), -sin(teta),         0], # Rotation around the z axis; a
             [sin(teta),  cos(teta),         0], # counterclockwise angle is
             [        0,          0,         1]] # positive
    
    p_prime = np.dot(rot_z, p)
    
    # Q is a point on the look direction of the second acquisition
    
    # 3) We compute the scalar product of vectors OP' and OQ to get the stereo
    #    angle and the B/H
    
    im2.s_comp = cos(im2.ortho) * sin(im2.scan)
    im2.o_comp = cos(im2.scan) * sin(im2.ortho)
    im2.z_comp = cos(im2.ortho) * cos(im2.scan)
    
    q = [im2.s_comp, im2.o_comp, im2.z_comp]
    
    # Scalar product:
    # (1) vect(OP').vect(OQ) = dist(OP') * dist(OQ) * cos(alpha) = cos(alpha)
    # (2) vect(OP').vect(OQ) = P'_s * Q_s + P'_o * Q_o + P'_z * Q_z
    # (1)&(2) -> alpha = acos(P'_s * Q_s + P'_o * Q_o + P'_z * Q_z)
    stereo_angle = acos(p_prime[0] * q[0]+ p_prime[1] * q[1] + p_prime[2] * \
        q[2])
    half_angle = stereo_angle / 2.
    b_to_h = tan(half_angle) * 2
    
    return [degrees(stereo_angle), b_to_h]

def repeatForEach(elements, times):
    """
    "a more user-friendly way to specify color for each arrow in 3D quiver?"
    https://github.com/matplotlib/matplotlib/issues/8484
    """
    return [e for e in elements for _ in range(times)]

if __name__ == "__main__":
    print('howstereo.py Copyright (C) 2020 Arthur Delorme\n')

    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=("Computes the B to H ratio of pairs of Pleiades or "
                "SPOT6|7 images"),
            epilog=textwrap.dedent('''
Note on azimuth angle:
    TL;DR: if the angle is from a SPOT6|7 DIMAP file, select "target" for --az-mode. Otherwise,
    select "scan" (default).

    The B/H is calculated using the azimuth of the scan axis (i.e. the angle between geographic
    north and the image line axis on the ground). In the Geostore, this angle corresponds to
    the Orientation angle. In the DIMAP file, for Pléiades, AZIMUTH_ANGLE also corresponds to
    this angle, but for SPOT6|7, AZIMUTH_ANGLE corresponds to the target azimuth. If you
    provide the target azimuth, you must select "target" for --az_mode, and the program will
    perform the conversion.

This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute
it under certain conditions.
See the GNU General Public License for more details.''')
        )
    parser.add_argument('--inc1', type=float, nargs=2,
        metavar=('ALONG', 'ACROSS'),
        help="incidence angles for image 1 (in degrees)")
    parser.add_argument('--az1', type=float, metavar='AZIMUTH',
        help="azimuth angle for image 1 (in degrees)")
    parser.add_argument('--inc2', type=float, nargs=2,
        metavar=('ALONG', 'ACROSS'),
        help="incidence angles for image 2 (in degrees)")
    parser.add_argument('--az2', type=float, metavar='AZIMUTH',
        help="azimuth angle for image 2 (in degrees)")
    parser.add_argument('--az_mode', type=str, choices=['scan', 'target'],
            default='scan', help=("type of azimuth angle (see the note below; "
                "default: %(default)s)")
        )
    parser.add_argument('--input_file', metavar='FILE',
            help=("input from a file instead (in csv format: "
                  "inc_along,inc_across,az)")
        )
    parser.add_argument('--show_plot', action='store_true',
        help="show a 3D, interactive plot")
    if len(sys.argv) == 1: # https://stackoverflow.com/a/4042861/13433994
        parser.print_help(sys.stderr)
        sys.exit()

    args = parser.parse_args()


    # The convention used for the angles is: positive counterclockwise. Thus,
    # some of the angles need to be converted from CNES convention to this one.


    # Read the input

    images = [] # A list of all the images
    if args.input_file:
        with open(args.input_file) as f:
            for i, l in enumerate(f):
                scan, ortho, az = l.split(',')
                scan = radians(float(scan))
                ortho = radians(float(ortho))
                az = radians(float(az))

                if args.az_mode == 'target':
                    az = (az + atan2(tan(ortho), tan(scan))) % (2 * pi)

                # Minus sign: conversion from the CNES convention to our
                ortho = -ortho
                az = -az

                images.append(Image('im{}'.format(i+1), scan, ortho, az))
    else:
        if not (args.inc1 and args.az1 and args.inc2 and args.az2):
            sys.exit("Error: missing arguments")
        scan1, ortho1 = args.inc1
        scan1 = radians(scan1)
        ortho1 = radians(ortho1)
        az1 = radians(args.az1)

        scan2, ortho2 = args.inc2
        scan2 = radians(scan2)
        ortho2 = radians(ortho2)
        az2 = radians(args.az2)

        if args.az_mode == 'target':
            az1 = (az1 + atan2(tan(ortho1), tan(scan1))) % (2 * pi)
            az2 = (az2 + atan2(tan(ortho2), tan(scan2))) % (2 * pi)

        # Minus sign: conversion from the CNES convention to our
        ortho1 = -ortho1 
        az1 = -az1
        ortho2 = -ortho2
        az2 = -az2

        images.extend([Image('im1', scan1, ortho1, az1),
                       Image('im2', scan2, ortho2, az2)])


    # Compute the B/H

    pairs = [] # A list of all possible pairs of images
    for p in list(combinations(images, 2)):
        im1, im2 = p
        stereo_angle, b_to_h = compute_b_to_h(im1, im2)
        pairs.append({
                'name': '{}-{}'.format(im1.name, im2.name),
                'stereo_angle': stereo_angle,
                'b_to_h': b_to_h
            })

    pairs = sorted(pairs, key=lambda k: k['b_to_h'])

    print('pair\tb/h\tangle')
    for p in pairs:
        print('{}\t{:.2f}\t{:.1f}°'.format(
                p['name'], p['b_to_h'], p['stereo_angle']
            ))


    # 3D plot

    if not args.show_plot:
        sys.exit()

    rays = np.zeros(shape=(len(images),6)) # Initialization; the rays to draw,
                                           # one per image
    for i, im in enumerate(images):
        rays[i] = [ 2*im.geo_n_comp,  2*im.geo_w_comp,  2*im.geo_z_comp,
                   -4*im.geo_n_comp, -4*im.geo_w_comp, -4*im.geo_z_comp]

    north = np.array([[-1, 0, 0, 2, 0, 0]])

    # Colors (https://github.com/matplotlib/matplotlib/issues/8484)
    cmap = plt.get_cmap('rainbow')
    c = list(np.linspace(0, 1, num=len(images)))
    c = c + repeatForEach(c, 2)

    x_r, y_r, z_r, u_r, v_r, w_r = zip(*rays)
    x_n, y_n, z_n, u_n, v_n, w_n = zip(*north)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.quiver(x_r, y_r, z_r, u_r, v_r, w_r, arrow_length_ratio=0.1,
        color=cmap(c), linewidths=3)
    ax.quiver(x_n, y_n, z_n, u_n, v_n, w_n, colors=[0,0,0], linewidths=2)

    xx, yy = np.meshgrid(range(-2,3), range(-2,3))
    zz = np.array(5*[5*[0]])
    ax.plot_surface(xx, yy, zz, alpha=0.1)

    ax.text(1, 0, 0, 'N', fontsize=30)
    for i, im in enumerate(images):
        ax.text(2*im.geo_n_comp, 2*im.geo_w_comp, 2*im.geo_z_comp,
            'i{}'.format(i+1), fontsize=20)

    lim = [-2.5, 2.5]
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_zlim(lim)
    ax.set_xlabel('North')
    ax.set_ylabel('West')
    ax.set_zlabel('Vertical')

    ax.view_init(elev=30., azim=-160)

    plt.tight_layout()
    plt.show()
