import os
import numpy as np
import pandas as pd
import ggmcalc
from pyshtools.shio import read_icgem_gfc


class GGM:

    def __init__(self, filename):

        self.filename = filename
        params = self.read_header()
        self.earth_gravity_constant = params['earth_gravity_constant']
        self.radius = params['radius']
        self.max_degree = params['max_degree']
        self.errors = params['errors']
        self.product_type = params['product_type']
        self.modelname = params['modelname']
        self.errors = params['errors']
        self.norm = params['norm']
        self.tide_system = params['tide_system']
        coeffs_and_size = self.read_coeffs()
        self.coeffs = coeffs_and_size[0]
        self.size = coeffs_and_size[1]

    def read_header(self):

        params = {}

        params['norm'] = 'unknown'
        params['tide_system'] = 'unknown'

        with open(self.filename, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                items = line.split()
                key, value = items[0:2]
                match key:
                    case 'product_type':
                        params[key] = value.strip()
                    case 'modelname':
                        params[key] = value.strip()
                    case 'earth_gravity_constant':
                        params[key] = float(value)
                    case 'radius':
                        params[key] = float(value)
                    case 'max_degree':
                        params[key] = int(value)
                    case 'errors':
                        params[key] = value.strip()
                    case 'norm':
                        print(value)
                        params[key] = value.strip()
                    case 'tide_system':
                        params[key] = value.strip()
                
                if 'end_of_head' in line:
                    break
        return params

    def read_coeffs(self):

        clim, gm, r0 = read_icgem_gfc(self.filename)
        nmax = self.max_degree

        C = clim[0]
        S = clim[1]

        size = nmax * (nmax + 2) + 1
        coeffs = np.zeros(size, dtype=np.float64)

        for n in range(nmax + 1):
            for m in range(n + 1):
                if m == 0:
                    idx = (nmax - m) * (nmax - m + 1) + nmax - n + 1
                    coeffs[idx - 1] = C[n, m]
                else:
                    idx = (nmax - m) * (nmax - m + 1) + 2 * (nmax - n) + 1
                    coeffs[idx - 1] = C[n, m]
                    coeffs[idx]     = S[n, m]

        coeffs = np.asfortranarray(coeffs, dtype=np.float64)
        return coeffs, len(coeffs) #, gm, r0

coef_file = os.path.join('dat', 'coeff.dat')
ellipsoid_config = os.path.join('dat', 'ellipsoid.dat')
input_file = os.path.join('dat', 'input.dat')

points = pd.read_csv(
    input_file,
    sep='\t',
    header=None,
    names=['latitude', 'longitude', 'elevation'],
)

with open(ellipsoid_config, 'r', encoding='utf8') as conf:

    gm = float(conf.readline().strip())
    a = float(conf.readline().strip())
    b = float(conf.readline().strip())
    f = float(conf.readline().strip())
    omega = float(conf.readline().strip())
    u0 = float(conf.readline().strip())
    ellps = conf.readline().strip()

ggm = GGM(coef_file)

ggmcalc.ggmcalc_mod.ggmcalc_main(
    gm,
    a,
    b,
    f,
    omega,
    u0,
    ellps,
    input_file,
    True,
    True,
    0,
    0.001,
    0.001,
    'test',
    ggm.coeffs,
    ggm.product_type,
    ggm.modelname,
    ggm.earth_gravity_constant,
    ggm.radius,
    ggm.max_degree,
    ggm.errors,
    ggm.norm,
    ggm.tide_system,
    points.latitude.astype(float),
    points.longitude.astype(float),
    points.elevation.astype(float),
)

print(ggmcalc.result.height_anomaly_p)