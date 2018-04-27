# This file is in the public domain.

def propagation_speed_marczak1997(temp_celsius):
    return 1.402385e3 +\
           5.038813 * temp_celsius - 5.799136e-2 * temp_celsius**2 +\
           3.287156e-4 * temp_celsius**3 - 1.398845e-6 * temp_celsius**4 +\
           2.787860e-9 * temp_celsius**5
