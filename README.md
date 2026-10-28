# cassini-parameterization

Pashkevich's Cassini-oval nuclear shape parameterization: a Fortran 2018
library with a C API and Python bindings. It takes `[alpha, alpha1, …, alphaN]`
(a Cassini oval in lemniscate coordinates plus a Legendre series, Nucl. Phys.
A 169 (1971) 275) and returns the axial surface: R(theta) and dR/dtheta at
caller thetas, the parametric rho(z) profile, the poles, the neck, the
resolved geometry, eps and the quadrupole moment, on the two-tier shape
parameterization contract shared with beta- and fos-parameterization.
