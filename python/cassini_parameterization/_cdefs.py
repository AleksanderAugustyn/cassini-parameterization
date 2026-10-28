"""ctypes signatures for the cassini-parameterization 0.1.0 C API.

One entry per exported symbol in ``include/cassini_parameterization.h``. Every
integer the C API takes or returns is ``int32_t``; every buffer is ``double*``.
"""
from __future__ import annotations

import ctypes

c_dbl = ctypes.c_double
c_dbl_p = ctypes.POINTER(ctypes.c_double)
c_i32 = ctypes.c_int32
c_i32_p = ctypes.POINTER(ctypes.c_int32)
c_void = ctypes.c_void_p

#: Largest parameter vector, ``[alpha, alpha1, ..., alpha20]`` (``CASSINI_PARAM_MAX_PARAMS``).
MAX_PARAMS = 21
#: Fewest scan nodes (``CASSINI_PARAM_N_X_GRID_FLOOR``).
N_X_GRID_FLOOR = 65
#: Fewest profile nodes (``CASSINI_PARAM_N_POINTS_FLOOR``).
N_POINTS_FLOOR = 4


class GeometryC(ctypes.Structure):
    """``cassini_param_geometry_t``, field for field in header order."""

    _fields_ = [
        ("eps", ctypes.c_double),
        ("c", ctypes.c_double),
        ("z_shift", ctypes.c_double),
        ("r_north", ctypes.c_double),
        ("r_south", ctypes.c_double),
        ("gap", ctypes.c_double),
        ("d_min", ctypes.c_double),
        ("is_separated", ctypes.c_int32),
        ("is_pinched", ctypes.c_int32),
        ("is_folded", ctypes.c_int32),
    ]


c_geom_p = ctypes.POINTER(GeometryC)

# name -> (argtypes, restype)
_SIGNATURES = {
    # --- diagnostics and eps ---
    "cassini_param_status_message": ([c_i32], ctypes.c_char_p),
    "cassini_param_epsilon": ([c_dbl_p, c_i32], c_dbl),

    # --- cache lifecycle; the trailing int32_t* receives the create status ---
    "cassini_param_cache_create": ([c_i32, c_i32, c_dbl_p, c_i32, c_i32_p], c_void),
    "cassini_param_cache_destroy": ([c_void], None),

    # --- cached computes, status is the return value ---
    "cassini_param_cache_radius_grid": ([c_void, c_dbl_p, c_i32, c_dbl_p, c_i32], c_i32),
    "cassini_param_cache_radius_and_derivative": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32], c_i32),
    "cassini_param_cache_radius_and_derivative_at_thetas": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_i32, c_dbl_p, c_dbl_p], c_i32),
    "cassini_param_cache_shape": ([c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_dbl_p], c_i32),
    "cassini_param_cache_rho_profile": (
        [c_void, c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_dbl_p], c_i32),
    "cassini_param_cache_rho_profile_unchecked": (
        [c_void, c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_dbl_p], c_i32),
    "cassini_param_cache_neck": ([c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32_p], c_i32),
    "cassini_param_cache_geometry": ([c_void, c_dbl_p, c_i32, c_geom_p], c_i32),
    "cassini_param_cache_quadrupole": ([c_void, c_dbl_p, c_i32, c_dbl_p], c_i32),

    # --- one-shot computes, trailing nullable status out-parameter ---
    "cassini_param_radius_grid": ([c_dbl_p, c_i32, c_i32, c_dbl_p, c_i32, c_dbl_p, c_i32_p], None),
    "cassini_param_radius_and_derivative": (
        [c_dbl_p, c_i32, c_i32, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32_p], None),
    "cassini_param_shape": ([c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_dbl_p, c_i32_p], None),
    "cassini_param_rho_profile": ([c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_dbl_p, c_i32_p], None),
    "cassini_param_rho_profile_unchecked": (
        [c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_dbl_p, c_i32_p], None),
    "cassini_param_neck": ([c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_i32_p, c_i32_p], None),
    "cassini_param_geometry": ([c_dbl_p, c_i32, c_i32, c_geom_p, c_i32_p], None),
    "cassini_param_quadrupole": ([c_dbl_p, c_i32, c_dbl_p, c_i32_p], None),
}


def configure(lib: ctypes.CDLL) -> ctypes.CDLL:
    """Set argtypes/restype on every exported symbol (idempotent).

    Parameters
    ----------
    lib : ctypes.CDLL
        The loaded shared library.

    Returns
    -------
    ctypes.CDLL
        The same object, with every signature attached.

    Raises
    ------
    AttributeError
        If the library lacks a 0.1.0 symbol.
    """
    for name, (argtypes, restype) in _SIGNATURES.items():
        func = getattr(lib, name)
        func.argtypes = argtypes
        func.restype = restype
    return lib
