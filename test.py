import os
import ggmcalc as ggm

coef_file = os.path.join('dat', 'coeff.dat')
ellipsoid_config = os.path.join('dat', 'ellipsoid.dat')
input_file = os.path.join('dat', 'input.dat')

ggm.ggmcalc_mod.ggmcalc_main(
    coef_file,
    ellipsoid_config,
    input_file,
    True,
    True,
    0,
    0.001,
    0.001,
    'test'
)

# print(ggm.inout_mod.height_anomaly_p)