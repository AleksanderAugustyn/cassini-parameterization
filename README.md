# cassini-parameterization

Pashkevich's Cassini-oval nuclear shape parameterization: a Fortran 2018
library with a C API and Python bindings. It takes `[alpha, alpha1, …, alphaN]`
(a Cassini oval in lemniscate coordinates plus a Legendre series, Nucl. Phys.
A 169 (1971) 275) and returns the axial surface: R(theta) and dR/dtheta at
caller thetas, the parametric rho(z) profile, the poles, the neck, the
resolved geometry, eps and the quadrupole moment.

The library follows the two-tier shape parameterization contract shared with
[beta-parameterization](https://github.com/AleksanderAugustyn/beta-parameterization),
[fos-parameterization](https://github.com/AleksanderAugustyn/fos-parameterization)
and [tqs-parameterization](https://github.com/AleksanderAugustyn/tqs-parameterization):
one-shot functions, plus a read-only cache that is built once and shared
across threads. Like beta and fos it is an expansion library: the consumer
chooses how many Legendre orders to carry.

## The shape

In Pashkevich's lemniscate coordinates (R, x) the lines R = const are Cassini
ovals with foci at z = ±sqrt(s), and the nuclear surface is

    R(x) = 1 + sum_{n=1}^{N} alpha_n P_n(x),    x in [-1, 1]

rotated about the z axis. With u = 2x² − 1 and p = sqrt(R⁴ + 2 s R² u + s²),

    rho(x) = sqrt((p − R² u − s) / 2),    z(x) = sign(x) sqrt((p + R² u + s) / 2).

The squared focal distance s = eps follows from alpha (Pashkevich eq. 10,
Jachimowicz et al. 2024 eq. 6), evaluated in its collapsed form

    eps = (alpha − 1)/4 [R(1)² + R(−1)²] + (alpha + 1)/2 R(0)²,

so the bare oval (all alpha_n = 0) has eps = alpha, and alpha = 1 gives
eps = R(0)² bit for bit: the pinch. `[0]` is the sphere. The volume is
normalized to 4π/3 (R0 = 1) and the centre of mass is the origin; the library
has no notion of A, Z or N.

**Parameters.** `params = [alpha, alpha1, …, alphaN]`, at most 21 values
(Legendre orders to 20). A cache takes `max_params` in 1 … 21 and accepts any
vector of 1 … max_params values; missing trailing values are zero, and a short
vector gives bitwise the same outputs as its zero-padded form. The validated
production setting is `max_params = 5`, the paper's {alpha, alpha1 … alpha4}.

**What is accepted and flagged.**

- Separated shapes (alpha > 1 on the bare oval, s > R(0)²) resolve: the
  geometry has `is_separated` and the gap, the unchecked profile draws both
  lobes with the two tips as consecutive nodes on the axis. The R(theta)
  outputs, `shape` and the checked profile report 103.
- Pinched shapes (s = R(0)² exactly) resolve with `is_pinched`; the R(theta)
  outputs and `shape` report 102 by rule.
- Folded shapes resolve and are served on every path when star-convex: a pole
  dents inward (a polar dimple, common for negative alpha4 in the paper's
  box), rho(z) becomes multi-valued, and `geometry` reports `is_folded`. The
  profile is the parametric curve, so it draws a dimple without inversion.
- Code 100: R(x) <= 0 somewhere on a fixed 1025-node check grid. Code 101: a
  pole leaves the axis (R(±1)² + s <= 0, a toroidal figure).

**The star gate.** R(theta) needs the surface to be star-convex from the
centre of mass. The gate quantity is the signed distance from the origin to
each tangent line of the profile, d = (rho z_x − z rho_x)/|(rho_x, z_x)|, and
the R(theta) outputs require min d >= `STAR_CONVEXITY_MARGIN` = 1e-3 (R0
units). The value is provisional and inclusive by design; `geometry` exports
`d_min` on every shape that resolves, so a later analysis can pick the value
from data. The scan grid (`n_x_grid`) locates the minima of d, and each
sampled local minimum is then refined by a golden-section search between its
scan neighbours: between nodes d can dip well below its samples (by 1.7e-2
on the 65-node floor scan), and a purely sampled gate passes shapes that are
not star-convex. With the refinement the gate verdict agrees with a
100,001-node scan on every point of the paper's lattice at alpha2 = 0 and
±0.1, at 65 and at 129 nodes, and d_min does not depend on `n_x_grid`.

**Numerics.** Volume and centre of mass are Pashkevich's x-integrals by
Gauss–Legendre on each half of [−1, 1] under x = ±sin²(πt/2), 512 nodes per
half. The substitution is quadratic at both ends of each half: at x = 0 it
smooths the exact pinch (where the integrand goes like sqrt|x|) and the
near-pinch, and at x = ±1 it does the same for a pole close to lifting off.
The relative volume error against independent references stays below 2e-14
through the pinch and the gap, below 6e-11 toward the toroidal boundary on
the paper's box, and below 3e-9 for the bare oval at alpha = −1 + delta,
where the pole radius and its conditioning both vanish. The lemniscate map
uses three conditioning rules — p as a hypot, the two radicands from their
product identity, and the derivative brackets chosen by a sign test so that
every bracket that appears is a same-sign sum — so rho is exact at the poles,
the pinch and in the void, and z is exact at a connected equator. R(theta) is
a safeguarded Newton root-find of theta(x) along the curve, bracketed on the
scan, with the analytic dtheta/dx; dR/dtheta is analytic. No finite
differences anywhere.

## Two tiers

1. **One-shot.** One call computes a shape and keeps nothing.
2. **Read-only cache.** Build a cache once (`max_params`, the scan resolution
   `n_x_grid`, a primary theta set), then request any number of shapes
   against it. Nothing derived from a shape parameter is stored, so no call
   influences a later one, and any number of threads may compute on one
   cache at once.

The two tiers return the same status and bitwise-identical outputs.

| Output | Cached | One-shot |
|---|---|---|
| R(theta) at the primary thetas | `cache_radius_grid_s` | `compute_radius_grid_standalone_s` |
| R(theta) and dR/dtheta | `cache_radius_and_derivative_s` | `compute_radius_and_derivative_standalone_s` |
| R(theta) and dR/dtheta at caller thetas | `cache_radius_and_derivative_at_thetas_s` | (the one-shot forms take caller thetas) |
| z_shift and pole distances, with the R(theta) verdict | `cache_shape_s` | `compute_shape_standalone_s` |
| parametric profile (z, rho, drho/dz) | `cache_rho_profile_s` | `compute_rho_profile_standalone_s` |
| profile, separated shapes included | `cache_rho_profile_unchecked_s` | `compute_rho_profile_unchecked_standalone_s` |
| neck position and radius | `cache_neck_s` | `compute_neck_standalone_s` |
| resolved geometry: eps, c, z_shift, poles, gap, d_min, flags | `cache_geometry_s` | `compute_geometry_standalone_s` |
| quadrupole moment q2 | `cache_quadrupole_s` | `compute_quadrupole_standalone_s` |
| eps from params | — | `compute_epsilon_f` |

The C functions are `cassini_param_cache_<output>` and `cassini_param_<output>`
(and `cassini_param_epsilon`); the Python names are the methods of `Cache` and
the module-level functions. `q2` is the dimensionless <r² P2(cos theta)> over
the volume; Pashkevich's Q is 16 q2 / 3.

### Fortran

```fortran
use cassini_parameterization_mod, only: cache_t, cache_init_s, cache_radius_grid_s, &
        compute_geometry_standalone_s, cassini_geometry_t

type(cache_t) :: cache
type(cassini_geometry_t) :: geometry
real(kind = rk) :: thetas(64), radii(64)
integer(kind = ik) :: status

! One-shot
call compute_geometry_standalone_s([0.65_rk, 0.2_rk, 0.0_rk, 0.025_rk, 0.05_rk], 2049_ik, &
        geometry, status)

! Cached
call cache_init_s(cache, 5_ik, 2049_ik, thetas, status)
call cache_radius_grid_s(cache, [0.25_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.075_rk], radii, status)
```

### C

```c
#include "cassini_parameterization.h"

const double params[5] = {0.65, 0.2, 0.0, 0.025, 0.05};
int32_t status;
cassini_param_cache_t *cache = cassini_param_cache_create(5, 2049, thetas, n_thetas, &status);
status = cassini_param_cache_radius_grid(cache, params, 5, radii, n_thetas);
cassini_param_geometry_t geometry;
cassini_param_geometry(params, 5, 2049, &geometry, &status);   /* one-shot */
cassini_param_cache_destroy(cache);
```

### Python

```python
import cassini_parameterization as cp

params = [0.65, 0.2, 0.0, 0.025, 0.05]
print(cp.geometry(params))                          # one-shot, n_x_grid = 2049

with cp.Cache(5, 2049, cp.theta_grid(64)) as cache:  # cached
    res = cache.radius_and_derivative(params)
    if res.ok:
        print(res.radii, res.dr_dtheta)
```

## Status codes

Every fallible routine reports an integer status; `0` is success. On any
nonzero status every output is zero-filled.

| Code | Name | Meaning |
|---|---|---|
| 1 | `too_many_params` | `max_params > 21`, or a one-shot call with more than 21 parameters |
| 2 | `cache_not_initialized` | compute on a cache that was not built |
| 3 | `invalid_grid` | no thetas, a theta outside `[0, pi]`, `n_x_grid < 65`, a profile under 4 nodes, or an output array of the wrong size |
| 4 | `wrong_param_count` | no parameters, or more than the cache's `max_params` |
| 5 | `invalid_init` | `max_params < 1` (and, in C, a negative size argument to a one-shot call) |
| 100 | `r_nonpositive` | R(x) <= 0 on the check grid |
| 101 | `toroidal` | a pole has left the symmetry axis |
| 102 | `not_star_convex` | R(theta) outputs and `shape` only: pinched, or star margin below 1e-3 |
| 103 | `separated` | R(theta) outputs, `shape` and the checked profile: two fragments |

Precedence on every path: usage codes, then 100, 101, 103, 102. `neck`,
`geometry`, `quadrupole` and the unchecked profile fail only with usage
codes, 100 and 101. Codes 1-5 are shared with the other shape libraries.
Inputs must be finite: the release build cannot detect NaN.

## Build and install

**Python.** 0.1.0 is not on PyPI yet. Build the wheel from a checkout:

    python -m pip install build
    python -m build --wheel --outdir dist_raw
    python -m pip install dist_raw/cassini_parameterization-*.whl

**CMake.** Fetch the library into a Fortran or C/C++ project:

```cmake
include(FetchContent)
FetchContent_Declare(
        cassini-parameterization
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/cassini-parameterization.git
        GIT_TAG 0.1.0
)
set(CASSINI_PARAM_BUILD_TESTS OFF)
FetchContent_MakeAvailable(cassini-parameterization)

target_link_libraries(my_target PRIVATE CassiniParameterization::cassini_parameterization)
```

Targets: `CassiniParameterization::cassini_parameterization` (static, Fortran
modules), `CassiniParameterization::cassini_parameterization_shared` (shared
library) and `CassiniParameterization::cassini_parameterization_cxx` (the C
header plus the shared library).

**From a checkout.**

    cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
    cmake --build build/release -j
    ctest --test-dir build/release --output-on-failure -E geometry_sweep
    CASSINI_PARAM_LIB=$PWD/build/release/libcassini_parameterization.so PYTHONPATH=$PWD/python \
        python -m pytest python/tests -q

The geometry sweep over the paper's lattice is its own ctest
(`ctest -R geometry_sweep`, label `sweep`); `--full` and `--a2 <value>` on
the executable run the spec's finer and alpha2 passes. The Python loader
prefers `CASSINI_PARAM_LIB`, then the library inside an installed wheel, then
`build/release`: rebuild `build/release` after pulling, or a stale build can
shadow the source you are testing.

## Version pins

| Dependency | Version |
|---|---|
| [gcc-compiler-options](https://github.com/AleksanderAugustyn/gcc-compiler-options) | 2.0.0 |
| [fortran-foundations](https://github.com/AleksanderAugustyn/fortran-foundations) | 3.0.0 |

Both are fetched by CMake. The library needs gfortran 10 or newer; it is
developed with gfortran 13. Python 3.9 or newer, NumPy 1.21 or newer.

## Changelog and license

See the
[changelog](https://github.com/AleksanderAugustyn/cassini-parameterization/blob/master/CHANGELOG.md).
MIT license.
