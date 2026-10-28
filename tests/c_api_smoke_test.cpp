// C API smoke test: argument validation, the status matrix on a few
// fixtures, NULL handling, one-shot entries with and without a status
// pointer, eps and the static messages. Exit code = number of failures.
#include "cassini_parameterization.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <numbers>
#include <vector>

namespace {
int failures = 0;

void check(const bool ok, const char* label) {
    if (!ok) {
        std::printf("FAIL: %s\n", label);
        ++failures;
    }
}

constexpr int n_theta = 64;
constexpr int n_x = 65;
constexpr int n_z = 101;

const double GS[5] = {0.25, 0.0, 0.0, 0.0, 0.075};
const double DIMPLED[5] = {0.5, 0.0, 0.0, 0.0, -0.2};
const double SEPARATED[1] = {1.05};

std::vector<double> pinned_thetas() {
    std::vector<double> t(n_theta);
    for (int i = 0; i < n_theta; ++i) t[i] = i * std::numbers::pi / (n_theta - 1);
    t[0] = 0.0;
    t[n_theta - 1] = std::numbers::pi;
    return t;
}

bool all_zero(const std::vector<double>& v) {
    for (const double x : v) if (x != 0.0) return false;
    return true;
}

bool finite_positive(const std::vector<double>& v) {
    for (const double x : v) if (!(std::isfinite(x) && x > 0.0)) return false;
    return true;
}

// Every cached compute on one parameter vector; returns the nine statuses.
struct Cached {
    int32_t st[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
    std::vector<double> radii = std::vector<double>(n_theta), radii_d = std::vector<double>(n_theta);
    std::vector<double> dr = std::vector<double>(n_theta), at_r = std::vector<double>(3);
    std::vector<double> at_dr = std::vector<double>(3);
    std::vector<double> z = std::vector<double>(n_z), rho = std::vector<double>(n_z);
    std::vector<double> drho = std::vector<double>(n_z), zu = std::vector<double>(n_z);
    std::vector<double> rhou = std::vector<double>(n_z), drhou = std::vector<double>(n_z);
    double z_shift = -1.0, r_north = -1.0, r_south = -1.0, z_neck = -1.0, rho_neck = -1.0, q2 = -1.0;
    int32_t found = -1;
    cassini_param_geometry_t geom{};
};

Cached run_cached(const cassini_param_cache_t* cache, const double* p, const int32_t n) {
    Cached c;
    const double at[3] = {0.4, 1.5, 2.7};
    c.st[0] = cassini_param_cache_radius_grid(cache, p, n, c.radii.data(), n_theta);
    c.st[1] = cassini_param_cache_radius_and_derivative(cache, p, n, c.radii_d.data(), c.dr.data(), n_theta);
    c.st[2] = cassini_param_cache_radius_and_derivative_at_thetas(cache, p, n, at, 3, c.at_r.data(),
                                                                  c.at_dr.data());
    c.st[3] = cassini_param_cache_shape(cache, p, n, &c.z_shift, &c.r_north, &c.r_south);
    c.st[4] = cassini_param_cache_rho_profile(cache, p, n, n_z, c.z.data(), c.rho.data(), c.drho.data());
    c.st[5] = cassini_param_cache_rho_profile_unchecked(cache, p, n, n_z, c.zu.data(), c.rhou.data(),
                                                        c.drhou.data());
    c.st[6] = cassini_param_cache_neck(cache, p, n, &c.z_neck, &c.rho_neck, &c.found);
    c.st[7] = cassini_param_cache_geometry(cache, p, n, &c.geom);
    c.st[8] = cassini_param_cache_quadrupole(cache, p, n, &c.q2);
    return c;
}

bool all_status(const Cached& c, const int32_t want) {
    for (const int32_t s : c.st) if (s != want) return false;
    return true;
}
}  // namespace

int main() {
    const std::vector<double> thetas = pinned_thetas();
    int32_t status = -1;

    // create
    check(cassini_param_cache_create(0, n_x, thetas.data(), n_theta, &status) == nullptr && status == 5,
          "create max_params 0 -> NULL, 5");
    check(cassini_param_cache_create(22, n_x, thetas.data(), n_theta, &status) == nullptr && status == 1,
          "create max_params 22 -> 1");
    check(cassini_param_cache_create(5, 64, thetas.data(), n_theta, &status) == nullptr && status == 3,
          "create n_x_grid 64 -> 3");
    check(cassini_param_cache_create(5, -1, thetas.data(), n_theta, &status) == nullptr && status == 3,
          "create negative n_x_grid -> 3");
    check(cassini_param_cache_create(5, n_x, nullptr, 0, &status) == nullptr && status == 3,
          "create no thetas -> 3");
    check(cassini_param_cache_create(5, n_x, thetas.data(), -2, &status) == nullptr && status == 3,
          "create negative n_thetas -> 3");
    {
        cassini_param_cache_t* tmp = cassini_param_cache_create(5, n_x, thetas.data(), n_theta, nullptr);
        check(tmp != nullptr, "create with a NULL status pointer runs");
        cassini_param_cache_destroy(tmp);
    }
    cassini_param_cache_t* cache = cassini_param_cache_create(5, n_x, thetas.data(), n_theta, &status);
    check(cache != nullptr && status == 0, "create 5 -> handle, 0");
    cassini_param_cache_t* cache21 = cassini_param_cache_create(21, n_x, thetas.data(), n_theta, &status);
    check(cache21 != nullptr && status == 0, "create 21 -> handle, 0");
    if (cache == nullptr || cache21 == nullptr) return failures + 1;

    // NULL handle and argument validation on the cached tier
    {
        const Cached c = run_cached(nullptr, GS, 5);
        check(all_status(c, 2), "NULL handle -> 2 on every cached compute");
        check(all_zero(c.radii) && c.r_north == 0.0 && c.found == 0 && c.geom.c == 0.0,
              "NULL handle zero-fills");
        std::vector<double> radii(n_theta, 7.0);
        const double six[6] = {0.25, 0, 0, 0, 0.075, 0.1};
        check(cassini_param_cache_radius_grid(cache, six, 6, radii.data(), n_theta) == 4, "n_params 6 -> 4");
        check(cassini_param_cache_radius_grid(cache, GS, 5, radii.data(), 63) == 3, "n_radii 63 -> 3");
        check(cassini_param_cache_radius_grid(cache, GS, 5, radii.data(), -1) == 3, "negative n_radii -> 3");
        check(cassini_param_cache_radius_grid(cache, nullptr, 0, radii.data(), n_theta) == 4,
              "NULL params, n_params 0 -> 4");
        check(cassini_param_cache_radius_grid(cache, GS, -3, radii.data(), n_theta) == 4,
              "negative n_params (cached) -> 4");
    }

    // GS
    {
        const Cached c = run_cached(cache, GS, 5);
        check(all_status(c, 0), "GS every cached compute 0");
        check(finite_positive(c.radii), "GS radii finite positive");
        check(c.geom.is_folded == 0 && c.found == 0, "GS not folded, no neck");
    }
    // DIMPLED
    {
        const Cached c = run_cached(cache, DIMPLED, 5);
        check(all_status(c, 0), "DIMPLED every cached compute 0");
        check(c.geom.is_folded == 1, "DIMPLED is_folded");
        check(c.found == 1 && c.rho_neck > 0.0, "DIMPLED neck found");
    }
    // SEPARATED
    {
        const Cached c = run_cached(cache, SEPARATED, 1);
        check(c.st[0] == 103 && all_zero(c.radii), "SEPARATED radius_grid 103, zero-filled");
        check(c.st[4] == 103, "SEPARATED rho_profile 103");
        check(c.st[5] == 0 && c.rhou[50] == 0.0 && c.rhou[51] == 0.0 && c.rhou[49] > 0.0 && c.rhou[52] > 0.0,
              "SEPARATED unchecked 0 with two zero interior nodes");
        check(c.st[7] == 0 && c.geom.is_separated == 1 && c.geom.gap > 0.0, "SEPARATED geometry");
        check(c.st[3] == 103 && c.r_north == 0.0, "SEPARATED shape 103 with r_north 0");
    }

    // one-shot entries
    {
        std::vector<double> radii(n_theta), dr(n_theta), z(n_z), rho(n_z), drho(n_z);
        double a = 0, b = 0, c = 0;
        int32_t found = -1;
        cassini_param_geometry_t g{};
        cassini_param_radius_grid(GS, 5, n_x, thetas.data(), n_theta, radii.data(), nullptr);
        cassini_param_radius_and_derivative(GS, 5, n_x, thetas.data(), n_theta, radii.data(), dr.data(), nullptr);
        cassini_param_shape(GS, 5, n_x, &a, &b, &c, nullptr);
        cassini_param_rho_profile(GS, 5, n_z, z.data(), rho.data(), drho.data(), nullptr);
        cassini_param_rho_profile_unchecked(GS, 5, n_z, z.data(), rho.data(), drho.data(), nullptr);
        cassini_param_neck(GS, 5, n_x, &a, &b, &found, nullptr);
        cassini_param_geometry(GS, 5, n_x, &g, nullptr);
        cassini_param_quadrupole(GS, 5, &a, nullptr);
        check(g.c > 0.0, "one-shot entries run with a NULL status pointer");

        cassini_param_radius_grid(SEPARATED, 1, n_x, thetas.data(), n_theta, radii.data(), &status);
        check(status == 103 && all_zero(radii), "one-shot SEPARATED radius_grid 103");
        cassini_param_shape(DIMPLED, 5, n_x, &a, &b, &c, &status);
        check(status == 0 && b > 0.0, "one-shot DIMPLED shape 0");
        cassini_param_radius_grid(GS, -1, n_x, thetas.data(), n_theta, radii.data(), &status);
        check(status == 5, "one-shot negative n_params -> 5");
        cassini_param_radius_grid(GS, 5, n_x, thetas.data(), -1, radii.data(), &status);
        check(status == 5, "one-shot negative n_thetas -> 5");
        cassini_param_rho_profile(GS, 5, -1, z.data(), rho.data(), drho.data(), &status);
        check(status == 5, "one-shot negative n_z -> 5");
        cassini_param_radius_grid(GS, 5, -1, thetas.data(), n_theta, radii.data(), &status);
        check(status == 3, "one-shot negative n_x_grid -> 3");
        cassini_param_geometry(GS, 5, -7, &g, &status);
        check(status == 3 && g.c == 0.0, "one-shot geometry negative n_x_grid -> 3, zeroed");
        std::vector<double> p22(22, 0.0);
        p22[0] = 0.25;
        cassini_param_quadrupole(p22.data(), 22, &a, &status);
        check(status == 1, "one-shot n_params 22 -> 1");
        cassini_param_quadrupole(nullptr, 0, &a, &status);
        check(status == 4 && a == 0.0, "one-shot n_params 0 -> 4");
    }

    // epsilon
    {
        const double half[1] = {0.5};
        check(cassini_param_epsilon(half, 1) == 0.5, "epsilon [0.5] -> 0.5");
        check(cassini_param_epsilon(nullptr, 0) == 0.0, "epsilon [] -> 0");
        std::vector<double> p22(22, 0.3);
        check(cassini_param_epsilon(p22.data(), 22) == 0.0, "epsilon 22 params -> 0");
        check(cassini_param_epsilon(half, -1) == 0.0, "epsilon negative n -> 0");
    }

    // messages
    {
        const int32_t codes[] = {0, 3, 100, 101, 102, 103, 999};
        for (const int32_t code : codes) {
            const char* m = cassini_param_status_message(code);
            check(m != nullptr && std::strlen(m) > 0, "status_message non-empty");
        }
        check(cassini_param_status_message(101) == cassini_param_status_message(101),
              "status_message pointer stable");
    }

    cassini_param_cache_destroy(cache);
    cassini_param_cache_destroy(cache21);
    cassini_param_cache_destroy(nullptr);
    std::printf("C API smoke test: %d failure(s)\n", failures);
    return failures;
}
