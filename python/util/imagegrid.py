# This file is in the public domain.

import numpy as np

def gen_sectorial_grid_xz(min_radius, max_radius, radius_step,
                          min_angle, max_angle, angle_step,
                          origin_x, origin_z):
    """Generate a sectorial grid and return the X-Z values.

    Parameters
    ----------
    min_radius : numeric
        Minimum radius (m).
    max_radius : numeric
        Maximum radius (m).
    radius_step : numeric
        Radius step (m).
    min_angle : numeric
        Minimum angle (degree).
    max_angle : numeric
        Maximum angle (degree).
    angle_step : numeric
        Angle step (degree).
    origin_x : numeric
        X-coordinate of the origin (m).
    origin_z : numeric
        Z-coordinate of the origin (m).

    Returns
    -------
    x : ndarray
        2D array with x values.
    z : ndarray
        2D array with z values.
    """

    radius_list = np.arange(min_radius, max_radius + 1.0e-6, radius_step)
    angle_list = np.deg2rad(np.arange(min_angle, max_angle + 1.0e-3, angle_step))
    angle, radius = np.meshgrid(angle_list, radius_list)
    x = radius * np.sin(angle) + origin_x
    z = radius * np.cos(angle) + origin_z
    return x, z

def gen_rectangular_grid():

    pass
