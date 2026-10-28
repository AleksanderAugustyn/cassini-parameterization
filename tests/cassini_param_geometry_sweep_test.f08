!> Geometry sweep over the paper's lattice (J eq. 7): alpha in [0, 0.95],
!! alpha1, alpha3, alpha4 in [-0.25, 0.25], alpha2 fixed.
!!
!! Tier A: counts of the R(theta)-path status {0, 100, 101, 102, 103} on the
!! 129-node scan, of folded shapes among the accepted ones, of shapes with a
!! neck, and of sampled-gate misses: accepted shapes whose d_min on a
!! 2049-node scan is below STAR_CONVEXITY_MARGIN (the gate is sampled, and
!! between 129-scan nodes d can dip by up to ~5e-3; two lattice points are
!! not star-convex at all). On the default grid they are asserted against
!! the goldens below (locked at STAR_CONVEXITY_MARGIN = 1e-3; re-lock when
!! the margin changes).
!! Tier B: for every accepted shape the 2049-node scan confirms, the volume and centre of mass from
!! R(theta) alone (|V/(4 pi/3) - 1|, |zcm| <= 1e-10), on composite
!! Gauss-Legendre panels in theta whose breakpoints are theta(x_k) of the
!! shape's 129-node profile (8 nodes per panel): a uniform theta rule misses
!! by 1.3e-3 at 200 thetas and 1.2e-5 at 4000 near the margin, where R(theta)
!! turns within ~d_min. And the round trip R(atan2(rho, z)) = hypot(z, rho)
!! at 21 profile nodes, x = 0 excluded (<= 1e-12). The maxima are printed.
!!
!! Default: step 0.05 at alpha2 = 0 (26,620 points). `--full`: step 0.025
!! (361,179 points). `--a2 <value>`: alpha2. Goldens apply to the default only.
!! OpenMP is optional; serial builds just run longer.
program cassini_param_geometry_sweep_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use cassini_parameterization_mod, only: cache_t, cassini_geometry_t, cache_init_s, &
            cache_radius_grid_s, cache_geometry_s, cache_neck_s, cache_rho_profile_s, &
            cache_radius_and_derivative_at_thetas_s, STAR_CONVEXITY_MARGIN, CASSINI_ERROR_R_NONPOSITIVE, &
            CASSINI_ERROR_TOROIDAL, CASSINI_ERROR_NOT_STAR_CONVEX, CASSINI_ERROR_SEPARATED
    use cassini_test_reference_mod, only: gl_thetas_s, panel_thetas_s, reference_theta_volume_s
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary
    implicit none

    real(kind = rk), parameter :: FOUR_THIRDS_PI = 4.0_rk * PI_C / 3.0_rk
    integer(kind = ik), parameter :: N_THETA = 200_ik, N_PROFILE = 21_ik
    integer(kind = ik), parameter :: N_PANEL_SCAN = 257_ik, N_PER_PANEL = 8_ik
    integer(kind = ik), parameter :: N_B = (N_PANEL_SCAN - 1_ik) * N_PER_PANEL
    !> Tier A goldens of the default grid: R(theta) statuses 0, 100, 101, 102,
    !! 103, folded among accepted, neck found, sampled-gate misses. Captured on
    !! the first green run.
    integer(kind = ik), parameter :: GOLDEN(8) = [25926_ik, 0_ik, 0_ik, 694_ik, 0_ik, 11497_ik, 12204_ik, 6_ik]
    integer(kind = ik), parameter :: N_X_DENSE = 2049_ik

    type(cache_t) :: cache, dense
    real(kind = rk) :: gl_t(N_THETA), gl_w(N_THETA), step, a2
    real(kind = rk) :: worst_v, worst_z, worst_rt
    integer(kind = ik) :: counts(8), n_alpha, n_box, ia, ib, status, n_total
    logical :: default_grid

    call parse_args(step, a2, default_grid)
    n_alpha = nint(0.95_rk / step, ik) + 1_ik
    n_box = nint(0.5_rk / step, ik) + 1_ik
    n_total = n_alpha * n_box**3
    call gl_thetas_s(N_THETA, gl_t, gl_w)
    call cache_init_s(cache, 5_ik, 129_ik, gl_t, status)
    call assert_int_eq(status, 0_ik, 'sweep cache')
    call cache_init_s(dense, 5_ik, N_X_DENSE, gl_t(1:1), status)
    call assert_int_eq(status, 0_ik, 'dense-scan cache')

    counts = 0_ik
    worst_v = 0.0_rk
    worst_z = 0.0_rk
    worst_rt = 0.0_rk
    !$omp parallel do collapse(2) schedule(dynamic) reduction(+:counts) &
    !$omp reduction(max:worst_v, worst_z, worst_rt)
    do ia = 0_ik, n_alpha - 1_ik
        do ib = 0_ik, n_box**3 - 1_ik
            call evaluate_shape(ia, ib, counts, worst_v, worst_z, worst_rt)
        end do
    end do
    !$omp end parallel do

    write(*, '(A,F6.3,A,F6.3,A,I0,A)') 'sweep step ', step, ', alpha2 = ', a2, ', ', n_total, ' points'
    write(*, '(A,5I8)') 'Tier A R(theta) status 0 / 100 / 101 / 102 / 103: ', counts(1:5)
    write(*, '(A,F7.3,A,F7.3,A)') 'Tier A accepted ', 100.0_rk * real(counts(1), rk) / real(n_total, rk), &
            '%, 102 ', 100.0_rk * real(counts(4), rk) / real(n_total, rk), '%'
    write(*, '(A,I8,A,F7.3,A)') 'Tier A folded among accepted: ', counts(6), ' (', &
            100.0_rk * real(counts(6), rk) / real(max(counts(1), 1_ik), rk), '%)'
    write(*, '(A,I8)') 'Tier A neck found: ', counts(7)
    write(*, '(A,I8)') 'Tier A sampled-gate misses (2049-scan d_min < margin): ', counts(8)
    write(*, '(A,3ES10.2)') 'Tier B worst |V/(4pi/3) - 1|, |zcm|, round trip: ', worst_v, worst_z, worst_rt
    write(*, '(A,8(I0,:,", "))') 'golden block: ', counts

    call assert_true(worst_v <= 1.0e-10_rk, 'Tier B volume')
    call assert_true(worst_z <= 1.0e-10_rk, 'Tier B centre of mass')
    call assert_true(worst_rt <= 1.0e-12_rk, 'Tier B round trip')
    if (default_grid) call assert_true(all(counts == GOLDEN), 'Tier A counts match the goldens')
    call test_summary()

contains

    !> One lattice point: Tier A tallies and Tier B maxima (the round trip is
    !! checked before the dense-scan confirmation, on every accepted shape).
    subroutine evaluate_shape(ia, ib, counts, worst_v, worst_z, worst_rt)
        integer(kind = ik), intent(in) :: ia, ib
        integer(kind = ik), intent(inout) :: counts(8)
        real(kind = rk), intent(inout) :: worst_v, worst_z, worst_rt
        type(cassini_geometry_t) :: geometry, geometry_dense
        real(kind = rk) :: params(5), radii(N_THETA), vol, zcm, z_neck, rho_neck
        real(kind = rk) :: zs(N_PANEL_SCAN), rhos(N_PANEL_SCAN), drhos(N_PANEL_SCAN), breaks(N_PANEL_SCAN)
        real(kind = rk) :: tb(N_B), wb(N_B), rb(N_B), drb(N_B)
        real(kind = rk) :: z(N_PROFILE), rho(N_PROFILE), drho(N_PROFILE)
        real(kind = rk) :: th(N_PROFILE - 1_ik), r_rt(N_PROFILE - 1_ik), dr_rt(N_PROFILE - 1_ik)
        integer(kind = ik) :: st, st_r, j, m
        logical :: found
        params(1) = real(ia, rk) * step
        params(2) = -0.25_rk + real(mod(ib, n_box), rk) * step
        params(3) = a2
        params(4) = -0.25_rk + real(mod(ib / n_box, n_box), rk) * step
        params(5) = -0.25_rk + real(ib / (n_box * n_box), rk) * step
        call cache_radius_grid_s(cache, params, radii, st_r)
        select case (st_r)
        case (0_ik)
            counts(1) = counts(1) + 1_ik
        case (CASSINI_ERROR_R_NONPOSITIVE)
            counts(2) = counts(2) + 1_ik
        case (CASSINI_ERROR_TOROIDAL)
            counts(3) = counts(3) + 1_ik
        case (CASSINI_ERROR_NOT_STAR_CONVEX)
            counts(4) = counts(4) + 1_ik
        case (CASSINI_ERROR_SEPARATED)
            counts(5) = counts(5) + 1_ik
        end select
        call cache_neck_s(cache, params, z_neck, rho_neck, found, st)
        if (found) counts(7) = counts(7) + 1_ik
        if (st_r /= 0_ik) return
        call cache_geometry_s(cache, params, geometry, st)
        if (geometry%is_folded) counts(6) = counts(6) + 1_ik
        call cache_rho_profile_s(cache, params, z, rho, drho, st)
        m = 0_ik
        do j = 1_ik, N_PROFILE
            if (j == (N_PROFILE + 1_ik) / 2_ik) cycle
            m = m + 1_ik
            th(m) = atan2(rho(j), z(j))
        end do
        call cache_radius_and_derivative_at_thetas_s(cache, params, th, r_rt, dr_rt, st)
        m = 0_ik
        do j = 1_ik, N_PROFILE
            if (j == (N_PROFILE + 1_ik) / 2_ik) cycle
            m = m + 1_ik
            worst_rt = max(worst_rt, abs(r_rt(m) - hypot(z(j), rho(j))))
        end do
        call cache_geometry_s(dense, params, geometry_dense, st)
        if (geometry_dense%d_min < STAR_CONVEXITY_MARGIN) then
            counts(8) = counts(8) + 1_ik
            return
        end if
        call cache_rho_profile_s(cache, params, zs, rhos, drhos, st)
        do j = 1_ik, N_PANEL_SCAN
            breaks(j) = atan2(rhos(j), zs(j))
        end do
        call panel_thetas_s(breaks, N_PER_PANEL, tb, wb)
        call cache_radius_and_derivative_at_thetas_s(cache, params, tb, rb, drb, st)
        call reference_theta_volume_s(tb, wb, rb, vol, zcm)
        worst_v = max(worst_v, abs(vol / FOUR_THIRDS_PI - 1.0_rk))
        worst_z = max(worst_z, abs(zcm))
    end subroutine evaluate_shape

    subroutine parse_args(step, a2, default_grid)
        real(kind = rk), intent(out) :: step, a2
        logical, intent(out) :: default_grid
        character(len = 64) :: arg
        integer :: i, n
        step = 0.05_rk
        a2 = 0.0_rk
        default_grid = .true.
        n = command_argument_count()
        i = 1
        do while (i <= n)
            call get_command_argument(i, arg)
            select case (trim(arg))
            case ('--full')
                step = 0.025_rk
                default_grid = .false.
            case ('--a2')
                i = i + 1
                call get_command_argument(i, arg)
                read(arg, *) a2
                default_grid = .false.
            case default
                write(*, '(A)') 'Usage: cassini_param_geometry_sweep_test [--full] [--a2 <value>]'
                stop 2
            end select
            i = i + 1
        end do
    end subroutine parse_args

end program cassini_param_geometry_sweep_test
