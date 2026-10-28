/**
 * @file cassini_parameterization.h
 * @brief C API for the Fortran Cassini-oval nuclear shape library (v0.1.0).
 *
 * The shape is Pashkevich's (Nucl. Phys. A 169 (1971) 275): in lemniscate
 * coordinates (R, x) the surface is R(x) = 1 + sum_n alpha_n P_n(x) on a
 * Cassini oval of squared focal distance eps, with eps from alpha by eq. 10.
 * params = [alpha, alpha1, ..., alphaN]; [alpha] alone is the bare oval and
 * [0] the sphere. All lengths are in R0 = 1 units with the centre of mass at
 * the origin.
 *
 * Two tiers of computation:
 *   - One-shot (`cassini_param_<compute>`): one call computes a shape and
 *     keeps nothing.
 *   - Read-only cache (`cassini_param_cache_*`): build a cache once with
 *     cassini_param_cache_create(), then request any number of shapes against
 *     it. The cache holds Legendre tables on fixed grids, the scan grid and
 *     the primary theta set; nothing derived from a shape parameter is
 *     stored, so no call influences a later one.
 *
 *   The two tiers return the same status and bitwise-identical outputs.
 *
 * Parameters:
 *   An expansion library: up to CASSINI_PARAM_MAX_PARAMS (21) parameters,
 *   alpha plus Legendre orders 1 .. 20. max_params at creation is 1 .. 21
 *   (larger -> 1, smaller -> 5). A cached call needs 1 <= n_params <=
 *   max_params, anything else is 4. A one-shot call with n_params > 21 is 1
 *   and with none is 4. Missing trailing parameters are zero: a short vector
 *   and its zero-padded form give bitwise-identical outputs.
 *
 * Thread safety:
 *   A cache is IMMUTABLE after creation. Every compute takes it `const`, and
 *   any number of threads may compute on one handle concurrently.
 *   cassini_param_cache_create() and cassini_param_cache_destroy() on a
 *   handle must not race with any other call on it. The one-shot calls hold
 *   no state.
 *
 * Diagnostics:
 *   cassini_param_cache_create() returns NULL on failure and reports the
 *   rejecting code through its nullable `status`. The cached computes return
 *   the status directly; the one-shot calls report through a trailing
 *   nullable `int32_t* status`. CASSINI_VALID (0) means success.
 *   cassini_param_status_message() maps a code to a static, null-terminated
 *   string owned by the library and safe to read from any thread.
 *
 * Order of the checks (the first failing one sets the status):
 *   create        max_params > 21 -> 1; < 1 -> 5; n_x_grid below
 *                 CASSINI_PARAM_N_X_GRID_FLOOR (negative included) -> 3; then
 *                 the theta set -> 3 (none, or one outside [0, pi])
 *   cached call   NULL handle -> 2; n_params outside 1 .. max_params -> 4
 *                 (a negative n_params is an empty vector); output size -> 3
 *                 (a negative size included); at-thetas only: the theta set
 *                 -> 3; then the value codes
 *   one-shot      a negative n_params, n_thetas or n_z -> 5; n_params > 21
 *                 -> 1; 0 -> 4; output size, n_x_grid and theta set -> 3;
 *                 then the value codes
 *   Value codes, in order: 100 (R(x) <= 0 on the check grid), 101 (a pole
 *   off the axis: toroidal), then 103 (separated) and 102 (pinched, or the
 *   star margin d_min below 1e-3) where an output gates them: the R(theta)
 *   outputs and shape apply both, the checked profile applies 103, the
 *   unchecked profile, neck, geometry and quadrupole neither.
 *
 * Failure behavior:
 *   On any nonzero status every output is zero-filled over its stated
 *   extent, `found` is 0, and a geometry record is all zero.
 *
 * Buffer sizes:
 *   Every compute takes the size of its output buffers explicitly: n_radii
 *   must equal the cache's theta count, n_thetas sizes both the caller's
 *   theta set and the outputs of the at-thetas form, and n_z (at least
 *   CASSINI_PARAM_N_POINTS_FLOOR) sizes all three profile arrays. A stated
 *   size must be the actual extent of the caller's buffer.
 *
 * Grids:
 *   n_x_grid (at least CASSINI_PARAM_N_X_GRID_FLOOR) is the scan resolution:
 *   the fold flag, the neck search and the R(theta) brackets are sampled on
 *   it, and it locates the minima of the star margin d, each of which is
 *   then refined between its scan neighbours, so the gate does not depend on
 *   n_x_grid. The neck is found from the scan: a neck narrower than the scan
 *   spacing can be missed at the floor.
 */
#ifndef CASSINI_PARAMETERIZATION_H
#define CASSINI_PARAMETERIZATION_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --- Limits and floors --- */
#define CASSINI_PARAM_MAX_PARAMS 21
#define CASSINI_PARAM_N_X_GRID_FLOOR 65
#define CASSINI_PARAM_N_POINTS_FLOOR 4

/* --- Status codes --- */
#define CASSINI_VALID 0
#define CASSINI_ERROR_TOO_MANY_PARAMS 1
#define CASSINI_ERROR_CACHE_NOT_INITIALIZED 2
#define CASSINI_ERROR_INVALID_GRID 3
#define CASSINI_ERROR_WRONG_PARAM_COUNT 4
#define CASSINI_ERROR_INVALID_INIT 5
#define CASSINI_ERROR_R_NONPOSITIVE 100
#define CASSINI_ERROR_TOROIDAL 101
#define CASSINI_ERROR_NOT_STAR_CONVEX 102
#define CASSINI_ERROR_SEPARATED 103

/** Opaque read-only cache handle. */
typedef struct cassini_param_cache cassini_param_cache_t;

/**
 * The resolved shape. eps (= s), the volume scale c, z_shift (the raw centre
 * of mass over c), the pole distances, the gap (separated shapes only), the
 * star margin d_min (0 for separated shapes) and three 0/1 flags.
 */
typedef struct {
    double eps;
    double c;
    double z_shift;
    double r_north;
    double r_south;
    double gap;
    double d_min;
    int32_t is_separated;
    int32_t is_pinched;
    int32_t is_folded;
} cassini_param_geometry_t;

/** Static, null-terminated description of a status code; never NULL. */
const char *cassini_param_status_message(int32_t status);

/** eps from params (0 .. 21 of them; none is the sphere). n_params < 0 or > 21 returns 0.0. */
double cassini_param_epsilon(const double *params, int32_t n_params);

/* --- Cache lifecycle --- */
cassini_param_cache_t *cassini_param_cache_create(int32_t max_params, int32_t n_x_grid,
                                                  const double *thetas, int32_t n_thetas,
                                                  int32_t *status);
void cassini_param_cache_destroy(cassini_param_cache_t *cache);

/* --- Cached computes (return the status) --- */
int32_t cassini_param_cache_radius_grid(const cassini_param_cache_t *cache, const double *params,
                                        int32_t n_params, double *radii, int32_t n_radii);
int32_t cassini_param_cache_radius_and_derivative(const cassini_param_cache_t *cache,
                                                  const double *params, int32_t n_params,
                                                  double *radii, double *dr_dtheta, int32_t n_radii);
int32_t cassini_param_cache_radius_and_derivative_at_thetas(const cassini_param_cache_t *cache,
                                                            const double *params, int32_t n_params,
                                                            const double *thetas, int32_t n_thetas,
                                                            double *radii, double *dr_dtheta);
int32_t cassini_param_cache_shape(const cassini_param_cache_t *cache, const double *params,
                                  int32_t n_params, double *z_shift, double *r_north, double *r_south);
int32_t cassini_param_cache_rho_profile(const cassini_param_cache_t *cache, const double *params,
                                        int32_t n_params, int32_t n_z, double *z, double *rho,
                                        double *drho_dz);
int32_t cassini_param_cache_rho_profile_unchecked(const cassini_param_cache_t *cache,
                                                  const double *params, int32_t n_params,
                                                  int32_t n_z, double *z, double *rho,
                                                  double *drho_dz);
int32_t cassini_param_cache_neck(const cassini_param_cache_t *cache, const double *params,
                                 int32_t n_params, double *z_neck, double *rho_neck, int32_t *found);
int32_t cassini_param_cache_geometry(const cassini_param_cache_t *cache, const double *params,
                                     int32_t n_params, cassini_param_geometry_t *geometry);
int32_t cassini_param_cache_quadrupole(const cassini_param_cache_t *cache, const double *params,
                                       int32_t n_params, double *q2);

/* --- One-shot computes (nullable trailing status) --- */
void cassini_param_radius_grid(const double *params, int32_t n_params, int32_t n_x_grid,
                               const double *thetas, int32_t n_thetas, double *radii,
                               int32_t *status);
void cassini_param_radius_and_derivative(const double *params, int32_t n_params, int32_t n_x_grid,
                                         const double *thetas, int32_t n_thetas, double *radii,
                                         double *dr_dtheta, int32_t *status);
void cassini_param_shape(const double *params, int32_t n_params, int32_t n_x_grid, double *z_shift,
                         double *r_north, double *r_south, int32_t *status);
void cassini_param_rho_profile(const double *params, int32_t n_params, int32_t n_z, double *z,
                               double *rho, double *drho_dz, int32_t *status);
void cassini_param_rho_profile_unchecked(const double *params, int32_t n_params, int32_t n_z,
                                         double *z, double *rho, double *drho_dz, int32_t *status);
void cassini_param_neck(const double *params, int32_t n_params, int32_t n_x_grid, double *z_neck,
                        double *rho_neck, int32_t *found, int32_t *status);
void cassini_param_geometry(const double *params, int32_t n_params, int32_t n_x_grid,
                            cassini_param_geometry_t *geometry, int32_t *status);
void cassini_param_quadrupole(const double *params, int32_t n_params, double *q2, int32_t *status);

#ifdef __cplusplus
}
#endif

#endif /* CASSINI_PARAMETERIZATION_H */
