!> Parity with the legacy ShapePlotterCassini.py for the four 236U shapes of
!! J fig. 6: a transcription cross-check, not a precision reference.
!!
!! Literals captured 2026-10-28 by a scratch script (not in the repo) that
!! built each shape with the legacy CassiniParameters(92, 144, alpha,
!! alpha_params = [alpha1, alpha2, alpha3, alpha4, 0]) and printed, in R0 =
!! r0 A^(1/3) units: calculate_epsilon(); c = (integrate_volume / (4 pi/3
!! R0^3))^(1/3) and calculate_zcm / R0, both on its 2000-point uniform
!! x_points (a midpoint rule, which sets the 1e-4 tolerance).
program cassini_param_legacy_test
    use precision_utilities_mod, only: ik, rk
    use cassini_parameterization_mod, only: cassini_geometry_t, compute_epsilon_f, &
            compute_geometry_standalone_s
    use cassini_test_fixtures_mod, only: GS, SADDLE1, SD_MIN, SADDLE2
    use test_utils_mod, only: assert_int_eq, assert_abs_close, test_summary
    implicit none

    call check(GS, 0.22729125976562492_rk, 0.9578121445945146_rk, 0.0_rk, 'GS')
    call check(SADDLE1, 0.359487060546875_rk, 0.9247680169321069_rk, 0.0_rk, 'SADDLE1')
    call check(SD_MIN, 0.5273404541015623_rk, 0.8802604703504848_rk, 0.0_rk, 'SD_MIN')
    call check(SADDLE2, 0.6544306640625_rk, 0.8582280542337426_rk, 0.2713283633363508_rk, 'SADDLE2')
    call test_summary()

contains

    !> G1: eps; G2: c and the raw centre of mass z_shift * c.
    subroutine check(params, eps, c, zcm, name)
        real(kind = rk), intent(in) :: params(:), eps, c, zcm
        character(len = *), intent(in) :: name
        type(cassini_geometry_t) :: geometry
        integer(kind = ik) :: status
        call assert_abs_close(compute_epsilon_f(params), eps, 1.0e-12_rk, 'G1 eps ' // name)
        call compute_geometry_standalone_s(params, 129_ik, geometry, status)
        call assert_int_eq(status, 0_ik, 'G2 status ' // name)
        call assert_abs_close(geometry%c, c, 1.0e-4_rk, 'G2 c ' // name)
        call assert_abs_close(geometry%z_shift * geometry%c, zcm, 1.0e-4_rk, 'G2 zm_raw ' // name)
    end subroutine check

end program cassini_param_legacy_test
