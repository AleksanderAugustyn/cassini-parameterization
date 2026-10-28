"""High-level API: Status, result objects, Cache, one-shot functions.

Two tiers, matching the C API. The module-level functions are one-shot: they
compute a shape and keep nothing. :class:`Cache` is the read-only cache: build
it once, then request any number of shapes against it. It holds nothing
derived from a shape parameter, so no call influences a later one. Both tiers
return the same status and bitwise-identical outputs (with the same
``n_x_grid``).

A :class:`Cache` is immutable after construction and may be shared between
threads: ctypes releases the GIL during a call, and the library computes on a
cache concurrently. :meth:`Cache.close` must not race with a compute.

Shape-validation failures come back as result objects carrying a
:class:`Status` with every numeric output zero-filled; :class:`CassiniParamError`
is raised only for usage errors -- a failed create, a closed handle, non-1-D
input, or more than ``MAX_PARAMS`` values given to :func:`epsilon`.

Units and frame: R0 = 1 (volume 4 pi / 3), centre of mass at the origin.

Precondition -- every ``params`` and ``thetas`` array must be finite. Release
builds assume no NaN, so non-finite input is undefined behavior.
"""
from __future__ import annotations

import ctypes
from dataclasses import dataclass
from enum import IntEnum
from typing import Optional

import numpy as np
import numpy.typing as npt

from ._cdefs import MAX_PARAMS, N_POINTS_FLOOR, N_X_GRID_FLOOR, GeometryC, c_dbl_p, configure
from ._libloader import load_library

#: Scan resolution of the one-shot functions. The star gate, the fold flag,
#: the neck search and the R(theta) brackets are sampled on it; at 129 nodes
#: the sampled gate passes a few lattice shapes whose dense margin is below
#: the threshold, at 2049 none.
N_X_GRID_DEFAULT = 2049

_lib: Optional[ctypes.CDLL] = None


def _get_lib() -> ctypes.CDLL:
    """Load and configure the shared library once per process."""
    global _lib
    if _lib is None:
        _lib = configure(load_library())
    return _lib


class CassiniParamError(RuntimeError):
    """Raised for usage errors: failed create, closed handle, non-1-D input.

    Attributes
    ----------
    status : Status or None
        The library's status code when the error came from a rejected cache
        construction; ``None`` for errors raised on the Python side.
    """

    def __init__(self, message: str, status: "Optional[Status]" = None) -> None:
        super().__init__(message)
        self.status = status


class Status(IntEnum):
    """Status codes of the C API (``CASSINI_*`` in the header).

    0-5 are the shared shape-parameterization contract codes, identical in
    every library of the family; 100+ are cassini's own.
    """

    valid = 0
    too_many_params = 1
    cache_not_initialized = 2
    invalid_grid = 3
    wrong_param_count = 4
    invalid_init = 5
    r_nonpositive = 100
    toroidal = 101
    not_star_convex = 102
    separated = 103


def status_message(status: int) -> str:
    """Fixed description of a status code.

    Parameters
    ----------
    status : int
        A :class:`Status` value or its integer code.

    Returns
    -------
    str
        The library's static message; ``"unknown status code"`` for codes it
        does not know.
    """
    return _get_lib().cassini_param_status_message(int(status)).decode()


# --- result objects --------------------------------------------------------

class _StatusMixin:
    """``ok`` and ``message`` for every result object."""

    status: Status

    @property
    def ok(self) -> bool:
        """True when the evaluation succeeded."""
        return self.status == Status.valid

    @property
    def message(self) -> str:
        """Description of :attr:`status`."""
        return status_message(self.status)


@dataclass(frozen=True)
class RadiusGridResult(_StatusMixin):
    """R(theta) on a theta set; zero-filled when ``status`` is nonzero."""

    radii: npt.NDArray[np.float64]
    status: Status


@dataclass(frozen=True)
class RadiusDerivativeResult(_StatusMixin):
    """R(theta) and dR/dtheta; both zero-filled on failure."""

    radii: npt.NDArray[np.float64]
    dr_dtheta: npt.NDArray[np.float64]
    status: Status


@dataclass(frozen=True)
class CassiniShape(_StatusMixin):
    """z_shift and the pole distances, with the R(theta) verdict.

    A pinched shape or one below the star margin is ``not_star_convex`` and a
    separated one ``separated``, with zeros: :func:`geometry` has the poles
    and flags of every shape that resolves.
    """

    z_shift: float
    r_north: float
    r_south: float
    status: Status


@dataclass(frozen=True)
class RhoProfileResult(_StatusMixin):
    """The parametric profile (z, rho, drho/dz), north pole first.

    The nodes follow the curve (cosine-spaced in the lemniscate coordinate),
    not a uniform z grid: folded (dimpled) shapes give a multi-valued rho(z).
    rho and drho/dz are 0 at both poles and, for pinched and separated
    shapes, at the two tips, which are consecutive nodes.
    """

    z: npt.NDArray[np.float64]
    rho: npt.NDArray[np.float64]
    drho_dz: npt.NDArray[np.float64]
    status: Status


@dataclass(frozen=True)
class NeckResult(_StatusMixin):
    """The narrowest interior minimum of rho, refined.

    Pinched and separated shapes are ``found`` with ``rho_neck = 0`` at the
    split plane. A valid shape without a neck is ``ok`` with ``found`` False
    and zeros.
    """

    z_neck: float
    rho_neck: float
    found: bool
    status: Status


@dataclass(frozen=True)
class Geometry(_StatusMixin):
    """The resolved shape.

    ``eps`` is the squared focal distance s, ``c`` the volume scale,
    ``z_shift`` the raw centre of mass over c, ``gap`` the distance between
    the two tips of a separated shape, ``d_min`` the star margin (the smallest
    distance from the centre of mass to a tangent line; 0 for separated
    shapes). Valid for pinched and separated shapes too.
    """

    eps: float
    c: float
    z_shift: float
    r_north: float
    r_south: float
    gap: float
    d_min: float
    is_separated: bool
    is_pinched: bool
    is_folded: bool
    status: Status


@dataclass(frozen=True)
class QuadrupoleResult(_StatusMixin):
    """q2 = <r^2 P2(cos theta)> over the volume in R0 = 1 units.

    Pashkevich's Q is ``16 q2 / 3``; a charge quadrupole moment is
    ``2 Z R0^2 q2``.
    """

    q2: float
    status: Status


# --- helpers ---------------------------------------------------------------

def theta_grid(n: int) -> npt.NDArray[np.float64]:
    """The open uniform theta grid ``theta_i = i*pi/(n+1)``, ``i = 1..n``.

    Open: neither pole is a node, and no node can round past pi. Pole radii
    come from :func:`geometry`.

    Parameters
    ----------
    n : int
        Number of interior points.

    Returns
    -------
    numpy.ndarray
        ``n`` thetas in radians, strictly inside ``(0, pi)``.
    """
    return np.arange(1, int(n) + 1, dtype=np.float64) * np.pi / (int(n) + 1.0)


def _as_1d(values: npt.ArrayLike, name: str) -> npt.NDArray[np.float64]:
    """Contiguous 1-D float64 copy or view of ``values``; CassiniParamError otherwise."""
    arr = np.ascontiguousarray(values, dtype=np.float64)
    if arr.ndim != 1:
        raise CassiniParamError(f"{name} must be 1-D, got shape {arr.shape}")
    return arr


def _ptr(arr: npt.NDArray[np.float64]) -> ctypes.POINTER(ctypes.c_double):  # type: ignore[valid-type]
    """Raw double pointer to a contiguous float64 array."""
    return arr.ctypes.data_as(c_dbl_p)


def _zeros(n: int) -> npt.NDArray[np.float64]:
    """An output buffer; a negative size gives an empty one and the library judges it."""
    return np.zeros(max(int(n), 0), dtype=np.float64)


def _geometry(g: GeometryC, status: int) -> Geometry:
    return Geometry(eps=g.eps, c=g.c, z_shift=g.z_shift, r_north=g.r_north, r_south=g.r_south,
                    gap=g.gap, d_min=g.d_min, is_separated=bool(g.is_separated),
                    is_pinched=bool(g.is_pinched), is_folded=bool(g.is_folded),
                    status=Status(status))


def epsilon(params: npt.ArrayLike) -> float:
    """eps (the squared focal distance s) from params, by eq. 10.

    Parameters
    ----------
    params : array_like
        ``[alpha, alpha1, ...]``, 0 .. ``MAX_PARAMS`` values; none is the sphere.

    Returns
    -------
    float
        eps; equal to alpha for the bare oval and to R(0)^2 at alpha = 1.

    Raises
    ------
    CassiniParamError
        For more than ``MAX_PARAMS`` values (never silently truncated).
    """
    p = _as_1d(params, "params")
    if p.size > MAX_PARAMS:
        raise CassiniParamError(f"epsilon takes at most {MAX_PARAMS} params, got {p.size}",
                                status=Status.too_many_params)
    return float(_get_lib().cassini_param_epsilon(_ptr(p), p.size))


# --- one-shot functions ----------------------------------------------------

def radius_grid(params: npt.ArrayLike, thetas: npt.ArrayLike, *,
                n_x_grid: int = N_X_GRID_DEFAULT) -> RadiusGridResult:
    """R(theta) for one shape.

    Parameters
    ----------
    params : array_like
        ``[alpha, alpha1, ...]``, 1 .. ``MAX_PARAMS`` values.
    thetas : array_like
        Evaluation angles in radians, at least one, all inside ``[0, pi]``.
    n_x_grid : int, keyword-only
        Scan resolution, at least ``N_X_GRID_FLOOR``.

    Returns
    -------
    RadiusGridResult
        Radii from the centre of mass and the status; zero-filled on failure.
    """
    p = _as_1d(params, "params")
    t = _as_1d(thetas, "thetas")
    radii = _zeros(t.size)
    status = ctypes.c_int32(0)
    _get_lib().cassini_param_radius_grid(_ptr(p), p.size, int(n_x_grid), _ptr(t), t.size,
                                         _ptr(radii), ctypes.byref(status))
    return RadiusGridResult(radii=radii, status=Status(status.value))


def radius_and_derivative(params: npt.ArrayLike, thetas: npt.ArrayLike, *,
                          n_x_grid: int = N_X_GRID_DEFAULT) -> RadiusDerivativeResult:
    """R(theta) and dR/dtheta for one shape (arguments as :func:`radius_grid`)."""
    p = _as_1d(params, "params")
    t = _as_1d(thetas, "thetas")
    radii = _zeros(t.size)
    dr_dtheta = _zeros(t.size)
    status = ctypes.c_int32(0)
    _get_lib().cassini_param_radius_and_derivative(_ptr(p), p.size, int(n_x_grid), _ptr(t), t.size,
                                                   _ptr(radii), _ptr(dr_dtheta), ctypes.byref(status))
    return RadiusDerivativeResult(radii=radii, dr_dtheta=dr_dtheta, status=Status(status.value))


def shape(params: npt.ArrayLike, *, n_x_grid: int = N_X_GRID_DEFAULT) -> CassiniShape:
    """z_shift and the poles with the R(theta) verdict (see :class:`CassiniShape`)."""
    p = _as_1d(params, "params")
    out = [ctypes.c_double(0.0) for _ in range(3)]
    status = ctypes.c_int32(0)
    _get_lib().cassini_param_shape(_ptr(p), p.size, int(n_x_grid), *(ctypes.byref(v) for v in out),
                                   ctypes.byref(status))
    z_shift, r_north, r_south = (v.value for v in out)
    return CassiniShape(z_shift=z_shift, r_north=r_north, r_south=r_south,
                        status=Status(status.value))


def _profile(params: npt.ArrayLike, n_points: int, checked: bool) -> RhoProfileResult:
    p = _as_1d(params, "params")
    z, rho, drho_dz = _zeros(n_points), _zeros(n_points), _zeros(n_points)
    status = ctypes.c_int32(0)
    lib = _get_lib()
    entry = lib.cassini_param_rho_profile if checked else lib.cassini_param_rho_profile_unchecked
    entry(_ptr(p), p.size, int(n_points), _ptr(z), _ptr(rho), _ptr(drho_dz), ctypes.byref(status))
    return RhoProfileResult(z=z, rho=rho, drho_dz=drho_dz, status=Status(status.value))


def rho_profile(params: npt.ArrayLike, n_points: int) -> RhoProfileResult:
    """The profile on ``n_points >= N_POINTS_FLOOR`` nodes; ``separated`` for two fragments."""
    return _profile(params, n_points, checked=True)


def rho_profile_unchecked(params: npt.ArrayLike, n_points: int) -> RhoProfileResult:
    """As :func:`rho_profile`, never ``separated``: draws both lobes."""
    return _profile(params, n_points, checked=False)


def neck(params: npt.ArrayLike, *, n_x_grid: int = N_X_GRID_DEFAULT) -> NeckResult:
    """The neck (see :class:`NeckResult`)."""
    p = _as_1d(params, "params")
    z_neck = ctypes.c_double(0.0)
    rho_neck = ctypes.c_double(0.0)
    found = ctypes.c_int32(0)
    status = ctypes.c_int32(0)
    _get_lib().cassini_param_neck(_ptr(p), p.size, int(n_x_grid), ctypes.byref(z_neck),
                                  ctypes.byref(rho_neck), ctypes.byref(found), ctypes.byref(status))
    return NeckResult(z_neck=z_neck.value, rho_neck=rho_neck.value, found=bool(found.value),
                      status=Status(status.value))


def geometry(params: npt.ArrayLike, *, n_x_grid: int = N_X_GRID_DEFAULT) -> Geometry:
    """The resolved shape (see :class:`Geometry`)."""
    p = _as_1d(params, "params")
    g = GeometryC()
    status = ctypes.c_int32(0)
    _get_lib().cassini_param_geometry(_ptr(p), p.size, int(n_x_grid), ctypes.byref(g),
                                      ctypes.byref(status))
    return _geometry(g, status.value)


def quadrupole(params: npt.ArrayLike) -> QuadrupoleResult:
    """The quadrupole moment (see :class:`QuadrupoleResult`)."""
    p = _as_1d(params, "params")
    q2 = ctypes.c_double(0.0)
    status = ctypes.c_int32(0)
    _get_lib().cassini_param_quadrupole(_ptr(p), p.size, ctypes.byref(q2), ctypes.byref(status))
    return QuadrupoleResult(q2=q2.value, status=Status(status.value))


# --- read-only cache -------------------------------------------------------

class Cache:
    """Read-only cache: Legendre tables, the scan grid and a primary theta set.

    Build once, evaluate any number of shapes. Nothing derived from a shape
    parameter is stored, so no call influences a later one. Immutable after
    construction and safe to share between threads; :meth:`close` must not
    race with a compute. The handle is released by :meth:`close`, by leaving
    a ``with`` block, or at garbage collection.

    Parameters
    ----------
    max_params : int
        Longest parameter vector the cache accepts, 1 .. ``MAX_PARAMS``.
    n_x_grid : int
        Scan resolution, at least ``N_X_GRID_FLOOR``.
    thetas : array_like
        Primary evaluation angles in radians, at least one, all inside
        ``[0, pi]``. See :func:`theta_grid`. The cache keeps its own copy.

    Attributes
    ----------
    max_params, n_x_grid, n_thetas : int
        The construction arguments.

    Raises
    ------
    CassiniParamError
        If the library refuses the arguments. Its ``status`` carries the code:
        ``too_many_params`` (``max_params > 21``), ``invalid_init``
        (``max_params < 1``) or ``invalid_grid`` (the scan grid or the theta
        set).
    """

    def __init__(self, max_params: int, n_x_grid: int, thetas: npt.ArrayLike) -> None:
        lib = _get_lib()
        t = _as_1d(thetas, "thetas")
        status = ctypes.c_int32(0)
        handle = lib.cassini_param_cache_create(int(max_params), int(n_x_grid), _ptr(t), t.size,
                                                ctypes.byref(status))
        if not handle:
            code = Status(status.value)
            raise CassiniParamError(
                f"Cache creation failed (max_params={max_params}, n_x_grid={n_x_grid}, "
                f"n_thetas={t.size}): {status_message(code)} [{code.name}]; need "
                f"1 <= max_params <= {MAX_PARAMS}, n_x_grid >= {N_X_GRID_FLOOR} and at least "
                f"one theta in [0, pi]", status=code)
        self._handle: Optional[ctypes.c_void_p] = ctypes.c_void_p(handle)
        # Bound here so __del__ never needs module globals during shutdown.
        self._destroy = lib.cassini_param_cache_destroy
        self.max_params = int(max_params)
        self.n_x_grid = int(n_x_grid)
        self.n_thetas = int(t.size)

    def close(self) -> None:
        """Destroy the underlying handle. Idempotent."""
        handle = getattr(self, "_handle", None)
        if handle is not None:
            self._handle = None
            self._destroy(handle)

    def __enter__(self) -> "Cache":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()

    def _require_handle(self) -> ctypes.c_void_p:
        handle = getattr(self, "_handle", None)
        if handle is None:
            raise CassiniParamError("Cache is closed")
        return handle

    def radius_grid(self, params: npt.ArrayLike) -> RadiusGridResult:
        """R(theta) at the cache's own thetas."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        radii = _zeros(self.n_thetas)
        status = _get_lib().cassini_param_cache_radius_grid(handle, _ptr(p), p.size, _ptr(radii),
                                                            radii.size)
        return RadiusGridResult(radii=radii, status=Status(status))

    def radius_and_derivative(self, params: npt.ArrayLike) -> RadiusDerivativeResult:
        """R(theta) and dR/dtheta at the cache's own thetas."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        radii = _zeros(self.n_thetas)
        dr_dtheta = _zeros(self.n_thetas)
        status = _get_lib().cassini_param_cache_radius_and_derivative(
            handle, _ptr(p), p.size, _ptr(radii), _ptr(dr_dtheta), radii.size)
        return RadiusDerivativeResult(radii=radii, dr_dtheta=dr_dtheta, status=Status(status))

    def radius_and_derivative_at_thetas(self, params: npt.ArrayLike,
                                        thetas: npt.ArrayLike) -> RadiusDerivativeResult:
        """R and dR/dtheta at caller thetas (at least one, all inside ``[0, pi]``).

        Given the cache's own thetas this returns the same bits as
        :meth:`radius_and_derivative`.
        """
        handle = self._require_handle()
        p = _as_1d(params, "params")
        t = _as_1d(thetas, "thetas")
        radii = _zeros(t.size)
        dr_dtheta = _zeros(t.size)
        status = _get_lib().cassini_param_cache_radius_and_derivative_at_thetas(
            handle, _ptr(p), p.size, _ptr(t), t.size, _ptr(radii), _ptr(dr_dtheta))
        return RadiusDerivativeResult(radii=radii, dr_dtheta=dr_dtheta, status=Status(status))

    def shape(self, params: npt.ArrayLike) -> CassiniShape:
        """z_shift and the poles with the R(theta) verdict."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        out = [ctypes.c_double(0.0) for _ in range(3)]
        status = _get_lib().cassini_param_cache_shape(handle, _ptr(p), p.size,
                                                      *(ctypes.byref(v) for v in out))
        z_shift, r_north, r_south = (v.value for v in out)
        return CassiniShape(z_shift=z_shift, r_north=r_north, r_south=r_south, status=Status(status))

    def _profile(self, params: npt.ArrayLike, n_points: int, checked: bool) -> RhoProfileResult:
        handle = self._require_handle()
        p = _as_1d(params, "params")
        z, rho, drho_dz = _zeros(n_points), _zeros(n_points), _zeros(n_points)
        lib = _get_lib()
        entry = (lib.cassini_param_cache_rho_profile if checked
                 else lib.cassini_param_cache_rho_profile_unchecked)
        status = entry(handle, _ptr(p), p.size, int(n_points), _ptr(z), _ptr(rho), _ptr(drho_dz))
        return RhoProfileResult(z=z, rho=rho, drho_dz=drho_dz, status=Status(status))

    def rho_profile(self, params: npt.ArrayLike, n_points: int) -> RhoProfileResult:
        """The checked profile on ``n_points`` nodes."""
        return self._profile(params, n_points, checked=True)

    def rho_profile_unchecked(self, params: npt.ArrayLike, n_points: int) -> RhoProfileResult:
        """The unchecked profile on ``n_points`` nodes."""
        return self._profile(params, n_points, checked=False)

    def neck(self, params: npt.ArrayLike) -> NeckResult:
        """The neck."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        z_neck = ctypes.c_double(0.0)
        rho_neck = ctypes.c_double(0.0)
        found = ctypes.c_int32(0)
        status = _get_lib().cassini_param_cache_neck(handle, _ptr(p), p.size, ctypes.byref(z_neck),
                                                     ctypes.byref(rho_neck), ctypes.byref(found))
        return NeckResult(z_neck=z_neck.value, rho_neck=rho_neck.value, found=bool(found.value),
                          status=Status(status))

    def geometry(self, params: npt.ArrayLike) -> Geometry:
        """The resolved shape."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        g = GeometryC()
        status = _get_lib().cassini_param_cache_geometry(handle, _ptr(p), p.size, ctypes.byref(g))
        return _geometry(g, status)

    def quadrupole(self, params: npt.ArrayLike) -> QuadrupoleResult:
        """The quadrupole moment."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        q2 = ctypes.c_double(0.0)
        status = _get_lib().cassini_param_cache_quadrupole(handle, _ptr(p), p.size, ctypes.byref(q2))
        return QuadrupoleResult(q2=q2.value, status=Status(status))


__all__ = [
    "Cache", "CassiniParamError", "CassiniShape", "Geometry", "NeckResult", "QuadrupoleResult",
    "RadiusDerivativeResult", "RadiusGridResult", "RhoProfileResult", "Status",
    "MAX_PARAMS", "N_POINTS_FLOOR", "N_X_GRID_FLOOR", "N_X_GRID_DEFAULT",
    "epsilon", "geometry", "neck", "quadrupole", "radius_and_derivative", "radius_grid",
    "rho_profile", "rho_profile_unchecked", "shape", "status_message", "theta_grid",
]
