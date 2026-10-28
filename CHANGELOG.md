# Changelog

All notable changes to this project are documented here. Versions follow
the git tags.

## 0.1.0 — unreleased

First release: Pashkevich's Cassini ovals with a Legendre series,
`params = [alpha, alpha1, …, alphaN]`, on the two-tier
shape-parameterization contract as an expansion library (hard cap 21,
validated at 5).

- Nine outputs in both tiers, bitwise identical between them: R(theta),
  R(theta) with dR/dtheta (at the cache's thetas and at caller thetas),
  z_shift and the pole distances, the checked and unchecked parametric
  profile, the neck, the resolved geometry (eps, c, z_shift, poles, gap,
  d_min and the separated, pinched and folded flags), and the quadrupole
  moment q2. Plus `compute_epsilon_f` for eps alone.
- eps by eq. 10 in its collapsed form; alpha = 1 gives eps = R(0)² bit for
  bit, so the pinch test is exact.
- Volume and centre of mass by Gauss–Legendre under x = ±sin²(πt/2),
  512 nodes per half: relative volume error below 2e-14 through the pinch
  and the gap, below 6e-11 toward the toroidal boundary on the paper's box.
- The lemniscate map with three conditioning rules (hypot p, radicands from
  their product identity, derivative brackets by a sign test) and regular
  forms of the derivatives at the connected equator and the poles.
- R(theta) by a bracketed, safeguarded Newton root-find of theta(x) along the
  curve; analytic dR/dtheta. No finite differences anywhere.
- C API (`cassini_param_`) with a `cassini_param_geometry_t` record; Python
  ctypes wheel (`cassini_parameterization`) with `Cache` and one-shot
  functions.

**Gate policy.** The R(theta) outputs require the signed distance from the
centre of mass to every tangent line of the profile to be at least
`STAR_CONVEXITY_MARGIN` = 1e-3 (R0 units), provisional and inclusive: folded
(dimpled) shapes are accepted when star-convex, and `d_min` is exported so
the value can be tightened from data. The gate is sampled on the scan grid;
at 129 nodes 6 of the 26,620 points of the paper's step-0.05 lattice pass
with a dense d_min below the margin (2 not star-convex), at 2049 none.
Pinched shapes are 102 by rule; separated shapes 103.

**Builder decisions** (the owner can veto):

1. A buffer whose extent disagrees with the cache or with its sibling
   arrays is code 3; a negative output size from C is the same code.
2. Floors: scan grid 65 nodes, profile 4 nodes, R-check grid 1025 nodes;
   the Python one-shot default `n_x_grid` is 2049. The quadrature has 512
   nodes per half (raised from the planned 64: see the numerics note above).
3. `geometry` is a record on every surface: `cassini_geometry_t`,
   `cassini_param_geometry_t`, `Geometry`.
4. The neck is the bracketed root of the analytic d rho/dx between the scan
   neighbours of the narrowest interior minimum (regula falsi, Illinois);
   when a minimum and a maximum share one scan interval, d rho/dx is sampled
   on 16 sub-intervals first. Ties within 1e-12 relative go north.
5. R(theta): safeguarded Newton on theta(x) - theta to 1e-14 in x, at most
   60 iterations; theta = 0 and pi return the analytic poles.
6. The one-shot forms that do not scan build their local cache on the floor
   grid and a one-node theta set, neither of which they consult.
7. Goldens live in `python/tests/test_api.py`, captured through the C API
   by `tests/golden_capture.f08`; `GOLDEN_ATOL = 1e-13`.
8. No standalone at-thetas form: the one-shot radius forms take thetas.
9. No tolerance on the toroidal test.
10. `compute_epsilon_f` accepts 0 … 21 parameters.
11. Legendre tables are `p(0:order, n_nodes)` per grid.
12. The quadrature substitution is x = ±sin²(πt/2), not the spec's x = ±t²:
    the latter treats only x = 0 and misses by up to 1.2e-6 near a pole
    close to lifting off (5.6e-5 on the bare oval).
