"""Python bindings for Pashkevich's Cassini-oval nuclear-shape parameterization.

``params = [alpha, alpha1, ..., alphaN]`` (up to 21 values): a Cassini oval in
lemniscate coordinates plus a Legendre series; ``[alpha]`` alone is the bare
oval and ``[0]`` the sphere. Two tiers: one-shot module-level functions, and
:class:`Cache`, a read-only cache built once and shared freely between
threads. Validation failures come back as a :class:`Status` on the result
object; only usage errors raise :class:`CassiniParamError`.
"""
from ._cdefs import MAX_PARAMS, N_POINTS_FLOOR, N_X_GRID_FLOOR
from ._libloader import load_library
from .api import (
    N_X_GRID_DEFAULT,
    Cache,
    CassiniParamError,
    CassiniShape,
    Geometry,
    NeckResult,
    QuadrupoleResult,
    RadiusDerivativeResult,
    RadiusGridResult,
    RhoProfileResult,
    Status,
    epsilon,
    geometry,
    neck,
    quadrupole,
    radius_and_derivative,
    radius_grid,
    rho_profile,
    rho_profile_unchecked,
    shape,
    status_message,
    theta_grid,
)

__all__ = [
    "Cache", "CassiniParamError", "CassiniShape", "Geometry", "NeckResult", "QuadrupoleResult",
    "RadiusDerivativeResult", "RadiusGridResult", "RhoProfileResult", "Status",
    "epsilon", "geometry", "neck", "quadrupole", "radius_and_derivative", "radius_grid",
    "rho_profile", "rho_profile_unchecked", "shape", "status_message", "theta_grid",
    "MAX_PARAMS", "N_POINTS_FLOOR", "N_X_GRID_FLOOR", "N_X_GRID_DEFAULT",
    "load_library",
]
