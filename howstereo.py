#!/usr/bin/python3
# -*- coding:utf-8 -*-

import sys
import argparse
from math import radians, degrees, sin, cos, acos, tan, atan
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# howstereo.py - Compute the B to H ratio of a couple of Pleiades or SPOT6|7
# images
# Copyright (C) 2020 Arthur Delorme
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

print('''howstereo.py Copyright (C) 2020 Arthur Delorme
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute
it under certain conditions.
See the GNU General Public License for more details.''')

parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Compute the B to H ratio of a couple of Pleiades or SPOT6|7 images'
    )
parser.add_argument('--incid1', type=float, nargs=2,
    metavar=('scan', 'ortho'),
    help="incidence angles for image 1 (in degrees)")
parser.add_argument('--az1', type=float, metavar='azimuth',
    help="scan azimuth for image 1 (in degrees)")
parser.add_argument('--incid2', type=float, nargs=2,
    metavar=('scan', 'ortho'),
    help="incidence angles for image 2 (in degrees)")
parser.add_argument('--az2', type=float, metavar='azimuth',
    help="scan azimuth for image 2 (in degrees)")
parser.add_argument('--show_plot', action='store_true',
    help="show a 3D, interactive plot")

args = parser.parse_args()

# The convention used for the angles is: positive counterclockwise. Thus, some
# of the angles need to be converted from CNES convention to this one.

scan1, ortho1 = args.incid1
scan1 = radians(scan1)
ortho1 = -radians(ortho1) # Minus sign: conversion from the CNES convention to our
scan2, ortho2 = args.incid2
scan2 = radians(scan2)
ortho2 = -radians(ortho2) # Minus sign: conversion from the CNES convention to our
az1 = -radians(args.az1) # Minus sign: conversion from the CNES convention to our
az2 = -radians(args.az2) # Minus sign: conversion from the CNES convention to our

# See the Figures 45 and 47 in the Pléiades Imagery User Guide
# See also the Figures 35, 36 and 37 in the SPOT 6 & SPOT 7 Imagery User Guide
# The coordinate system is (O, Scan axis, OrthoScan axis, Normale to the ground)
# -> (O, Scan, Ortho, z)
# We define a point P on the look direction of the first acquisition with OP = 1
# beta: incidence angle
# i_s : incidence angle in the Scan axis direction
# i_o : incidence angle in the Ortho axis direction
# p_s : coordinate of P on the Scan axis
# p_o : coordinate of P on the Ortho axis
# p_z : coordinate of P on the z axis
# 1) We express P's coordinates in the reference of the first acquisition

p_s = cos(ortho1) * sin(scan1)
p_o = cos(scan1) * sin(ortho1)
p_z = cos(ortho1) * cos(scan1)

# 2) We express P's coordinates in the reference of the second acquisition with
#    a rotation around the z axis

teta = az2 - az1 # Azimuth difference between the two acquisitions
print('\nRotation: {:.1f}°\n'.format(-degrees(teta))) # Minus sign: conversion
                                                      # from our convention to
                                                      # the CNES convention

p = [p_s, p_o, p_z]
rot_z = [[cos(teta), -sin(teta), 0], [sin(teta), cos(teta), 0], [0, 0, 1]]
# rot_z: rotation around the z axis; a counterclockwise angle is positive

p_prime = np.dot(rot_z, p)

# We define a point Q on the look direction of the second acquisition
# 3) We compute the scalar product of vectors OP' and OQ to get the angle and
#    B/H

q_s = cos(ortho2) * sin(scan2)
q_o = cos(scan2) * sin(ortho2)
q_z = cos(ortho2) * cos(scan2)

q = [q_s, q_o, q_z]

# Scalar product:
# (1) vect(OP').vect(OQ) = dist(OP')*dist(OQ)*cos(alpha) = cos(alpha)
# (2) vect(OP').vect(OQ) = P'x*Qx + P'y*Qy + P'z*Qz
# (1) and (2) -> alpha = acos(P'x*Qx + P'y*Qy + P'z*Qz)
stereo_angle = acos(q[0] * p_prime[0] + q[1] * p_prime[1] + q[2] * p_prime[2])
half_angle = stereo_angle / 2.
b_to_h = tan(half_angle) * 2

scan1_prime = atan(p_prime[0]/p_prime[2])
ortho1_prime = atan(p_prime[1]/p_prime[2])
print('     Scan   OrthoScan')
print('P : {:5.1f}° {:5.1f}°'.format(degrees(scan1), -degrees(ortho1)))
print('P\': {:5.1f}° {:5.1f}°'.format(degrees(scan1_prime),
    -degrees(ortho1_prime)))
print('Q : {:5.1f}° {:5.1f}°'.format(degrees(scan2), -degrees(ortho2)))
print('\nStereo angle: {:.1f}°'.format(degrees(stereo_angle)))
print('-> B to H ratio: {:.3f}'.format(b_to_h))
# (Minus signs for the OrthoScan angles: conversions from our convention to the
# CNES convention)

# 3D plot

if not args.show_plot:
    sys.exit()

# We express Q's and P''s coordinates in the geographic reference through a
# rotation around the z axis, using the azimuth of the second acquisition as
# the angle (with a minus sign because the transformation goes in the opposite
# direction with respect to the azimuth)

rot_z = [[cos(-az2), -sin(-az2), 0], [sin(-az2), cos(-az2), 0], [0, 0, 1]]
q_north = np.dot(rot_z, q)
p_prime_north = np.dot(rot_z, p_prime)

q_north_n = q_north[0]
q_north_w = q_north[1]
q_north_z = q_north[2]
p_prime_north_n = p_prime_north[0]
p_prime_north_w = p_prime_north[1]
p_prime_north_z = p_prime_north[2]

rays = np.array([
        [2*p_prime_north_n, 2*p_prime_north_w, 2*p_prime_north_z,
            -4*p_prime_north_n, -4*p_prime_north_w, -4*p_prime_north_z],
        [2*q_north_n, 2*q_north_w, 2*q_north_z,
            -4*q_north_n, -4*q_north_w,-4*q_north_z]
    ])
north = np.array([[-1, 0, 0, 2, 0, 0]])

color1 = (0,0,1)
color2 = (1,0,0)
color3 = (0,0,0)
x_r, y_r, z_r, u_r, v_r, w_r = zip(*rays)
x_n, y_n, z_n, u_n, v_n, w_n = zip(*north)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver3D(x_r, y_r, z_r, u_r, v_r, w_r, arrow_length_ratio=0.1,
    colors=[color1, color2, color1, color1, color2, color2], linewidths=2)
ax.quiver3D(x_n, y_n, z_n, u_n, v_n, w_n, colors=[color3], linewidths=1)
xx, yy = np.meshgrid(range(-2,3), range(-2,3))
zz = np.array(5*[5*[0]])
ax.plot_surface(xx, yy, zz, alpha=0.1)
ax.text(1, 0, 0, 'N', 'x')
ax.text(2*p_prime_north_n, 2*p_prime_north_w, 2*p_prime_north_z, 'I1', 'x')
ax.text(2*q_north_n, 2*q_north_w, 2*q_north_z, 'I2', 'x')
ax.set_xlim([-4, 4])
ax.set_ylim([-4, 4])
ax.set_zlim([-4, 4])
ax.set_xlabel('North')
ax.set_ylabel('West')
ax.set_zlabel('Vertical')
plt.show()
