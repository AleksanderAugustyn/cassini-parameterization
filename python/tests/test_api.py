"""Contract tests for the 0.1.0 Python surface.

Shape-validation failures are results, not exceptions: the buffers come back
zero-filled with a nonzero :class:`Status`. Only usage errors -- a failed
create, a closed handle, non-1-D input, an over-long eps vector -- raise
``CassiniParamError``.

The golden block is captured by ``tests/golden_capture.f08`` through the same
C entries these tests call. The recapture procedure lives in that file's
header; never loosen ``GOLDEN_ATOL`` to make a recapture fit.
"""
from __future__ import annotations

import threading

import numpy as np
import pytest

import cassini_parameterization as cp

# --- fixtures (params = [alpha, alpha1, alpha2, alpha3, alpha4]) ------------

SPHERE = [0.0]
GS = [0.25, 0.0, 0.0, 0.0, 0.075]
SADDLE1 = [0.35, 0.0, 0.0, 0.0, -0.075]
SD_MIN = [0.525, 0.0, 0.0, 0.0, 0.025]
SADDLE2 = [0.65, 0.2, 0.0, 0.025, 0.05]
ASYM = [0.8, 0.15, -0.05, 0.1, -0.05]
DIMPLED = [0.5, 0.0, 0.0, 0.0, -0.2]
FOLD_BOX = [0.125, 0.25, 0.0, 0.2, -0.175]
STAR_FAIL = [0.95, -0.2, 0.0, -0.05, -0.05]
PINCHED = [1.0]
PINCHED_ASYM = [1.0, 0.2, 0.0, 0.025, 0.05]
SEPARATED = [1.05]
R_NEG = [0.0, -1.1]
TOROIDAL = [-1.1]

N_X = 129
N_Z = 181
GOLDEN_IDX = [0, 45, 90, 135, 180]
# pi/8, pi/2, 7pi/8 as pre-rounded literals, bit-identical to golden_capture.f08.
DERIV_THETAS = [0.39269908169872414, 1.5707963267948966, 2.748893571891069]
# A Debug capture against Release and wheel builds: a few ulp of values of
# order unity, two orders of magnitude inside this.
GOLDEN_ATOL = 1e-13


def _same_bits(a, b) -> bool:
    """True when ``a`` and ``b`` hold the same IEEE754 bit patterns.

    ``np.array_equal`` calls ``+0.0`` and ``-0.0`` equal; the contract's
    bitwise claims need the bytes.
    """
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    return a.shape == b.shape and a.tobytes() == b.tobytes()


def _inclusive_grid(n: int) -> np.ndarray:
    """``n`` uniform nodes on the CLOSED ``[0, pi]``, both endpoints pinned."""
    grid = np.arange(n, dtype=np.float64) * np.pi / (n - 1.0)
    grid[0] = 0.0
    grid[-1] = np.pi
    return grid


# --- golden block (tests/golden_capture.f08, Debug) ------------------------
GS_GEOMETRY = [
     2.2729125976562503E-001,
     9.5781323167470889E-001,
     0.0000000000000000E+000,
     1.2277702594586477E+000,
     1.2277702594586477E+000,
     0.0000000000000000E+000,
     9.5102675357042044E-001,
]
GS_FLAGS = (False, False, False)
GS_NECK_FOUND = False
GS_NECK = [
     0.0000000000000000E+000,
     0.0000000000000000E+000,
]
GS_Q2 = [
     1.0000052816422994E-001,
]
GS_DERIV_R = [
     1.1299479938682273E+000,
     9.5102675357042044E-001,
     1.1299479938682273E+000,
]
GS_DERIV_DR = [
    -4.0875546775773469E-001,
    -0.0000000000000000E+000,
     4.0875546775773453E-001,
]
GS_Z = [
     1.2277702594586477E+000,
     8.0681483985226243E-001,
     0.0000000000000000E+000,
    -8.0681483985226243E-001,
    -1.2277702594586477E+000,
]
GS_RHO = [
     0.0000000000000000E+000,
     6.3497713326402205E-001,
     9.5102675357042044E-001,
     6.3497713326402205E-001,
     0.0000000000000000E+000,
]
SADDLE1_GEOMETRY = [
     3.5948706054687496E-001,
     9.2476837149657332E-001,
     0.0000000000000000E+000,
     1.1919970976005940E+000,
     1.1919970976005940E+000,
     0.0000000000000000E+000,
     8.2591409898816615E-001,
]
SADDLE1_FLAGS = (False, False, False)
SADDLE1_NECK_FOUND = True
SADDLE1_NECK = [
     0.0000000000000000E+000,
     8.2711327651249000E-001,
]
SADDLE1_Q2 = [
     1.9101842039893721E-001,
]
SADDLE1_DERIV_R = [
     1.1929333230336638E+000,
     8.2711327651249000E-001,
     1.1929333230336638E+000,
]
SADDLE1_DERIV_DR = [
    -9.1518322606482905E-002,
    -0.0000000000000000E+000,
     9.1518322606482905E-002,
]
SADDLE1_Z = [
     1.1919970976005940E+000,
     9.3038797959933661E-001,
     0.0000000000000000E+000,
    -9.3038797959933661E-001,
    -1.1919970976005940E+000,
]
SADDLE1_RHO = [
     0.0000000000000000E+000,
     6.6728233788824476E-001,
     8.2711327651249000E-001,
     6.6728233788824476E-001,
     0.0000000000000000E+000,
]
SD_MIN_GEOMETRY = [
     5.2734045410156227E-001,
     8.8026121291944814E-001,
     0.0000000000000000E+000,
     1.4270433945713412E+000,
     1.4270433945713412E+000,
     0.0000000000000000E+000,
     7.9643266260391687E-001,
]
SD_MIN_FLAGS = (False, False, False)
SD_MIN_NECK_FOUND = False
SD_MIN_NECK = [
     0.0000000000000000E+000,
     0.0000000000000000E+000,
]
SD_MIN_Q2 = [
     3.0749773364673766E-001,
]
SD_MIN_DERIV_R = [
     1.2924894592567362E+000,
     7.9643266260391687E-001,
     1.2924894592567415E+000,
]
SD_MIN_DERIV_DR = [
    -6.0577524313106768E-001,
    -0.0000000000000000E+000,
     6.0577524313106079E-001,
]
SD_MIN_Z = [
     1.4270433945713412E+000,
     1.0287250408071629E+000,
     0.0000000000000000E+000,
    -1.0287250408071629E+000,
    -1.4270433945713412E+000,
]
SD_MIN_RHO = [
     0.0000000000000000E+000,
     6.1458343668056781E-001,
     7.9643266260391687E-001,
     6.1458343668056781E-001,
     0.0000000000000000E+000,
]
SADDLE2_GEOMETRY = [
     6.5443066406249994E-001,
     8.5822913759649999E-001,
     3.1615265662871961E-001,
     1.4432671218507112E+000,
     1.6624668729108560E+000,
     0.0000000000000000E+000,
     7.4851232695938064E-001,
]
SADDLE2_FLAGS = (False, False, False)
SADDLE2_NECK_FOUND = False
SADDLE2_NECK = [
     0.0000000000000000E+000,
     0.0000000000000000E+000,
]
SADDLE2_Q2 = [
     3.9468761615721543E-001,
]
SADDLE2_DERIV_R = [
     1.3166621763197803E+000,
     7.7445086634339744E-001,
     1.3630585594390094E+000,
]
SADDLE2_DERIV_DR = [
    -5.8340265164910665E-001,
    -1.1789847988949095E-001,
     1.1578528206652472E+000,
]
SADDLE2_Z = [
     1.4432671218507112E+000,
     8.6716626556398069E-001,
    -3.1615265662871961E-001,
    -1.3648612346373206E+000,
    -1.6624668729108560E+000,
]
SADDLE2_RHO = [
     0.0000000000000000E+000,
     7.1536304713677912E-001,
     7.2149717529294322E-001,
     4.5966324540727166E-001,
     0.0000000000000000E+000,
]
ASYM_GEOMETRY = [
     8.2403515624999979E-001,
     7.8828491503764142E-001,
     2.4870281667044572E-001,
     1.6098965597763146E+000,
     1.6650478043224941E+000,
     0.0000000000000000E+000,
     4.2730024593304999E-001,
]
ASYM_FLAGS = (False, False, True)
ASYM_NECK_FOUND = True
ASYM_NECK = [
    -2.4870281667044569E-001,
     5.5077847766065346E-001,
]
ASYM_Q2 = [
     6.8763992262183427E-001,
]
ASYM_DERIV_R = [
     1.4731262369543894E+000,
     5.8524174480757318E-001,
     1.5620325145406855E+000,
]
ASYM_DERIV_DR = [
    -7.0975252530674637E-001,
    -1.4808426429869420E-001,
     1.0298330024387086E+000,
]
ASYM_Z = [
     1.6098965597763146E+000,
     1.1060906512630073E+000,
    -2.4870281667044572E-001,
    -1.5187244983713699E+000,
    -1.6650478043224941E+000,
]
ASYM_RHO = [
     0.0000000000000000E+000,
     7.1369186389161021E-001,
     5.5077847766065346E-001,
     5.3557987996224465E-001,
     0.0000000000000000E+000,
]
DIMPLED_GEOMETRY = [
     4.8171875000000003E-001,
     8.9005337018471975E-001,
     0.0000000000000000E+000,
     1.1899424030675905E+000,
     1.1899424030675905E+000,
     0.0000000000000000E+000,
     6.2733974554103478E-001,
]
DIMPLED_FLAGS = (False, False, True)
DIMPLED_NECK_FOUND = True
DIMPLED_NECK = [
     0.0000000000000000E+000,
     6.8701356673122360E-001,
]
DIMPLED_Q2 = [
     3.1808468010521934E-001,
]
DIMPLED_DERIV_R = [
     1.2855754844654665E+000,
     6.8701356673122360E-001,
     1.2855754844654665E+000,
]
DIMPLED_DERIV_DR = [
     1.0694723470036263E-001,
    -0.0000000000000000E+000,
    -1.0694723470036263E-001,
]
DIMPLED_Z = [
     1.1899424030675905E+000,
     1.0498153660903118E+000,
     0.0000000000000000E+000,
    -1.0498153660903118E+000,
    -1.1899424030675905E+000,
]
DIMPLED_RHO = [
     0.0000000000000000E+000,
     7.0287358480296014E-001,
     6.8701356673122360E-001,
     7.0287358480296014E-001,
     0.0000000000000000E+000,
]
# --- end golden block --------------------------------------------------------

GOLDEN_SHAPES = {
    "GS": GS, "SADDLE1": SADDLE1, "SD_MIN": SD_MIN, "SADDLE2": SADDLE2,
    "ASYM": ASYM, "DIMPLED": DIMPLED,
}


def _golden(name: str, suffix: str):
    return globals()[f"{name}{suffix}"]


# --- theta grid --------------------------------------------------------------

def test_theta_grid_is_open_and_uniform() -> None:
    grid = cp.theta_grid(16)
    assert grid.size == 16 and np.all(grid > 0.0) and np.all(grid < np.pi)
    assert np.allclose(np.diff(grid), np.pi / 17.0, rtol=0.0, atol=1e-15)
    with cp.Cache(5, N_X, grid) as cache:
        assert cache.n_thetas == 16 and cache.n_x_grid == N_X and cache.max_params == 5


# --- physics -----------------------------------------------------------------

def test_sphere_is_unit() -> None:
    grid = cp.theta_grid(32)
    radii = cp.radius_grid(SPHERE, grid, n_x_grid=N_X)
    assert radii.ok and np.allclose(radii.radii, 1.0, rtol=0.0, atol=1e-14)
    shape = cp.shape(SPHERE, n_x_grid=N_X)
    assert shape.ok and abs(shape.r_north - 1.0) < 1e-14 and abs(shape.r_south - 1.0) < 1e-14
    assert abs(shape.z_shift) < 1e-15
    assert not cp.neck(SPHERE, n_x_grid=N_X).found
    assert abs(cp.geometry(SPHERE, n_x_grid=N_X).c - 1.0) < 1e-14
    assert cp.epsilon([0.0]) == 0.0
    padded = cp.radius_grid([0.0, 0.0, 0.0], grid, n_x_grid=N_X)
    assert _same_bits(radii.radii, padded.radii)


@pytest.mark.parametrize("name", list(GOLDEN_SHAPES))
def test_goldens_match_fortran(name: str) -> None:
    params = GOLDEN_SHAPES[name]
    with cp.Cache(5, N_X, DERIV_THETAS) as cache:
        g = cache.geometry(params)
        assert g.ok
        geo = [g.eps, g.c, g.z_shift, g.r_north, g.r_south, g.gap, g.d_min]
        np.testing.assert_allclose(geo, _golden(name, "_GEOMETRY"), rtol=0.0, atol=GOLDEN_ATOL)
        assert (g.is_separated, g.is_pinched, g.is_folded) == _golden(name, "_FLAGS")
        nk = cache.neck(params)
        assert nk.ok and nk.found == _golden(name, "_NECK_FOUND")
        np.testing.assert_allclose([nk.z_neck, nk.rho_neck], _golden(name, "_NECK"), rtol=0.0,
                                   atol=GOLDEN_ATOL)
        q = cache.quadrupole(params)
        assert q.ok
        np.testing.assert_allclose([q.q2], _golden(name, "_Q2"), rtol=0.0, atol=GOLDEN_ATOL)
        rd = cache.radius_and_derivative(params)
        assert rd.ok
        np.testing.assert_allclose(rd.radii, _golden(name, "_DERIV_R"), rtol=0.0, atol=GOLDEN_ATOL)
        np.testing.assert_allclose(rd.dr_dtheta, _golden(name, "_DERIV_DR"), rtol=0.0,
                                   atol=GOLDEN_ATOL)
        prof = cache.rho_profile(params, N_Z)
        assert prof.ok
        np.testing.assert_allclose(prof.z[GOLDEN_IDX], _golden(name, "_Z"), rtol=0.0, atol=GOLDEN_ATOL)
        np.testing.assert_allclose(prof.rho[GOLDEN_IDX], _golden(name, "_RHO"), rtol=0.0,
                                   atol=GOLDEN_ATOL)


@pytest.mark.parametrize("params", [SADDLE2, DIMPLED])
def test_derivative_matches_finite_difference(params) -> None:
    thetas = np.linspace(0.05, np.pi - 0.05, 41)
    h = 1e-6
    with cp.Cache(5, N_X, thetas) as cache:
        rd = cache.radius_and_derivative(params)
        up = cache.radius_and_derivative_at_thetas(params, thetas + h)
        down = cache.radius_and_derivative_at_thetas(params, thetas - h)
    assert rd.ok and up.ok and down.ok
    fd = (up.radii - down.radii) / (2.0 * h)
    assert np.all(np.abs(fd - rd.dr_dtheta) <= 1e-6 * (1.0 + np.abs(rd.dr_dtheta)))


def test_profile_tips_are_exact_zero() -> None:
    prof = cp.rho_profile(GS, 101)
    assert prof.ok and prof.rho[0] == 0.0 and prof.rho[-1] == 0.0
    assert prof.drho_dz[0] == 0.0 and prof.drho_dz[-1] == 0.0
    sep = cp.rho_profile_unchecked(SEPARATED, 100)
    assert sep.ok and sep.rho[49] == 0.0 and sep.rho[50] == 0.0
    assert sep.z[49] > 0.0 > sep.z[50]


def test_unchecked_profile_draws_a_separated_shape() -> None:
    checked = cp.rho_profile(SEPARATED, 101)
    assert checked.status == cp.Status.separated and not checked.ok
    assert np.all(checked.rho == 0.0) and np.all(checked.z == 0.0)
    unchecked = cp.rho_profile_unchecked(SEPARATED, 101)
    assert unchecked.ok
    zeros = np.flatnonzero(unchecked.rho[1:-1] == 0.0) + 1
    assert list(zeros) == [50, 51]
    assert unchecked.rho[49] > 0.0 and unchecked.rho[52] > 0.0


def test_unchecked_equals_checked_on_connected() -> None:
    a = cp.rho_profile(FOLD_BOX, 257)
    b = cp.rho_profile_unchecked(FOLD_BOX, 257)
    assert a.ok and b.ok
    assert _same_bits(a.z, b.z) and _same_bits(a.rho, b.rho) and _same_bits(a.drho_dz, b.drho_dz)


def test_geometry_flags() -> None:
    sep = cp.geometry(SEPARATED, n_x_grid=N_X)
    assert sep.ok and sep.is_separated and sep.gap > 0.0 and sep.d_min == 0.0
    pin = cp.geometry(PINCHED, n_x_grid=N_X)
    assert pin.ok and pin.is_pinched and not pin.is_separated
    dim = cp.geometry(DIMPLED, n_x_grid=N_X)
    assert dim.ok and dim.is_folded
    assert cp.radius_grid(DIMPLED, cp.theta_grid(8), n_x_grid=N_X).ok
    star = cp.geometry(STAR_FAIL, n_x_grid=N_X)
    assert star.ok and star.d_min < 0.0
    assert cp.radius_grid(STAR_FAIL, cp.theta_grid(8), n_x_grid=N_X).status == cp.Status.not_star_convex
    bad = cp.geometry(R_NEG, n_x_grid=N_X)
    assert bad.status == cp.Status.r_nonpositive and not bad.ok
    assert all(getattr(bad, f) == 0.0 for f in ("eps", "c", "z_shift", "r_north", "r_south", "gap", "d_min"))
    assert not (bad.is_separated or bad.is_pinched or bad.is_folded)


def test_epsilon() -> None:
    assert cp.epsilon([0.5]) == 0.5
    g = cp.geometry(PINCHED_ASYM, n_x_grid=N_X)
    assert g.is_pinched
    r0 = 1.0 + 0.05 * 0.375
    assert abs(cp.epsilon(PINCHED_ASYM) - r0**2) < 1e-15
    assert cp.epsilon([]) == 0.0
    with pytest.raises(cp.CassiniParamError):
        cp.epsilon(np.zeros(22))


@pytest.mark.parametrize("params,code", [(R_NEG, cp.Status.r_nonpositive),
                                         (TOROIDAL, cp.Status.toroidal)])
def test_r_nonpositive_and_toroidal_are_statuses(params, code) -> None:
    for res in _all_one_shot(params):
        assert res.status == code and not res.ok


def test_wrong_count() -> None:
    with cp.Cache(5, N_X, cp.theta_grid(8)) as cache:
        assert cache.radius_grid(np.zeros(6)).status == cp.Status.wrong_param_count
    assert cp.radius_grid(np.zeros(22), cp.theta_grid(8)).status == cp.Status.too_many_params
    assert cp.radius_grid([], cp.theta_grid(8)).status == cp.Status.wrong_param_count


@pytest.mark.parametrize("args,code", [
    ((5, N_X, []), cp.Status.invalid_grid),
    ((0, N_X, [1.0]), cp.Status.invalid_init),
    ((22, N_X, [1.0]), cp.Status.too_many_params),
    ((5, 64, [1.0]), cp.Status.invalid_grid),
])
def test_cache_rejections_raise(args, code) -> None:
    with pytest.raises(cp.CassiniParamError) as info:
        cp.Cache(*args)
    assert info.value.status == code


def test_closed_cache_raises() -> None:
    cache = cp.Cache(5, N_X, cp.theta_grid(4))
    cache.close()
    cache.close()
    with pytest.raises(cp.CassiniParamError):
        cache.radius_grid(GS)


def test_inputs_are_converted_to_contiguous_float64() -> None:
    thetas = cp.theta_grid(8)
    with cp.Cache(5, N_X, thetas) as cache:
        ref = cache.radius_grid(GS).radii
        as_list = cache.radius_grid([0.25, 0, 0, 0, 0.075]).radii
        as_f32 = cache.radius_grid(np.array([0.25, 0, 0, 0, 0.075], dtype=np.float32))
        strided = np.zeros(10)
        strided[::2] = GS
        as_view = cache.radius_grid(strided[::2]).radii
        as_tuple = cache.radius_grid(tuple(GS)).radii
    assert _same_bits(ref, as_list) and _same_bits(ref, as_view) and _same_bits(ref, as_tuple)
    assert as_f32.ok
    with pytest.raises(cp.CassiniParamError):
        cp.radius_grid(np.zeros((2, 2)), thetas)


def _all_one_shot(params, thetas=None, n_x_grid=N_X, n_points=N_Z):
    thetas = cp.theta_grid(8) if thetas is None else thetas
    return [
        cp.radius_grid(params, thetas, n_x_grid=n_x_grid),
        cp.radius_and_derivative(params, thetas, n_x_grid=n_x_grid),
        cp.shape(params, n_x_grid=n_x_grid),
        cp.rho_profile(params, n_points),
        cp.rho_profile_unchecked(params, n_points),
        cp.neck(params, n_x_grid=n_x_grid),
        cp.geometry(params, n_x_grid=n_x_grid),
        cp.quadrupole(params),
    ]


def _all_cached(cache, params, thetas, n_points=N_Z):
    return [
        cache.radius_grid(params),
        cache.radius_and_derivative(params),
        cache.shape(params),
        cache.rho_profile(params, n_points),
        cache.rho_profile_unchecked(params, n_points),
        cache.neck(params),
        cache.geometry(params),
        cache.quadrupole(params),
        cache.radius_and_derivative_at_thetas(params, thetas),
    ]


def _numbers(res):
    out = []
    for value in vars(res).values():
        if isinstance(value, np.ndarray):
            out.extend(value.tolist())
        elif isinstance(value, float):
            out.append(value)
    return out


@pytest.mark.parametrize("params", [R_NEG, TOROIDAL, SEPARATED, PINCHED, STAR_FAIL])
def test_failures_zero_fill_every_output(params) -> None:
    for res in _all_one_shot(params):
        if res.ok:
            continue
        assert all(v == 0.0 for v in _numbers(res)), type(res).__name__
        for flag in ("found", "is_separated", "is_pinched", "is_folded"):
            assert not getattr(res, flag, False)


def test_profile_floor() -> None:
    assert cp.N_POINTS_FLOOR == 4
    prof = cp.rho_profile(GS, 4)
    assert prof.ok and prof.rho[1] > 0.0 and prof.rho[2] > 0.0
    assert cp.rho_profile(GS, 3).status == cp.Status.invalid_grid
    with cp.Cache(5, N_X, cp.theta_grid(8)) as cache:
        assert cache.rho_profile(GS, 3).status == cp.Status.invalid_grid


def _same_result(a, b) -> bool:
    if type(a) is not type(b) or a.status != b.status:
        return False
    for name, value in vars(a).items():
        other = getattr(b, name)
        if isinstance(value, (np.ndarray, float)):
            if not _same_bits(value, other):
                return False
        elif value != other:
            return False
    return True


@pytest.mark.parametrize("params", [SADDLE2, DIMPLED, SEPARATED, R_NEG])
def test_one_shot_equals_cached_bits(params) -> None:
    thetas = cp.theta_grid(8)
    with cp.Cache(len(params), N_X, thetas) as cache:
        cached = _all_cached(cache, params, thetas)
    for a, b in zip(_all_one_shot(params, thetas), cached):
        assert _same_result(a, b), type(a).__name__


@pytest.mark.parametrize("params", [SADDLE2, DIMPLED])
def test_at_thetas_equals_grid_bits(params) -> None:
    thetas = cp.theta_grid(16)
    with cp.Cache(5, N_X, thetas) as cache:
        grid = cache.radius_and_derivative(params)
        at = cache.radius_and_derivative_at_thetas(params, thetas)
    assert _same_bits(grid.radii, at.radii) and _same_bits(grid.dr_dtheta, at.dr_dtheta)


@pytest.mark.parametrize("short,padded", [([0.5], [0.5, 0.0, 0.0]), (GS, GS + [0.0, 0.0])])
def test_short_equals_padded_bits(short, padded) -> None:
    thetas = cp.theta_grid(8)
    for a, b in zip(_all_one_shot(short, thetas), _all_one_shot(padded, thetas)):
        assert _same_result(a, b), type(a).__name__
    with cp.Cache(7, N_X, thetas) as cache:
        for a, b in zip(_all_cached(cache, short, thetas), _all_cached(cache, padded, thetas)):
            assert _same_result(a, b), type(a).__name__


def _walk_shape(t: int, s: int):
    return [0.3 + 0.006 * s + 0.01 * t, 0.002 * t, 0.0, 0.001 * s - 0.05, 0.05 - 0.001 * t]


def test_cache_is_shared_between_threads() -> None:
    thetas = _inclusive_grid(32)
    n_threads, n_steps = 8, 100
    with cp.Cache(5, N_X, thetas) as cache:
        serial = [[cache.radius_and_derivative(_walk_shape(t, s)) for s in range(n_steps)]
                  for t in range(n_threads)]
        threaded = [None] * n_threads

        def work(t: int) -> None:
            threaded[t] = [cache.radius_and_derivative(_walk_shape(t, s)) for s in range(n_steps)]

        workers = [threading.Thread(target=work, args=(t,)) for t in range(n_threads)]
        for w in workers:
            w.start()
        for w in workers:
            w.join()
    n_valid = 0
    for t in range(n_threads):
        for a, b in zip(serial[t], threaded[t]):
            assert _same_result(a, b)
            n_valid += a.ok
    assert n_valid > 0


def test_status_message_and_enum_coverage() -> None:
    for code in cp.Status:
        assert cp.status_message(code)
    assert cp.status_message(999) == "unknown status code"
    assert cp.MAX_PARAMS == 21 and cp.N_X_GRID_FLOOR == 65 and cp.N_X_GRID_DEFAULT == 2049
