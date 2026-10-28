!> Status-code tests: cassini codes are >= 100, shared shape_core codes
!! re-export, every code maps to a non-empty message, and the fixed constants
!! hold.
program cassini_param_status_test
    use precision_utilities_mod, only: ik, rk
    use cassini_parameterization_mod, only: CASSINI_VALID, CASSINI_ERROR_R_NONPOSITIVE, &
            CASSINI_ERROR_TOROIDAL, CASSINI_ERROR_NOT_STAR_CONVEX, CASSINI_ERROR_SEPARATED, &
            CASSINI_MAX_PARAMS, CASSINI_MAX_ORDER, CASSINI_N_GL_HALF, CASSINI_N_R_CHECK, &
            CASSINI_N_X_GRID_FLOOR, CASSINI_N_POINTS_FLOOR, STAR_CONVEXITY_MARGIN, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, status_message
    use test_utils_mod, only: assert_int_eq, assert_true, assert_abs_close, test_summary
    implicit none

    call assert_int_eq(int(CASSINI_VALID), int(SHAPE_VALID), 'CASSINI_VALID aliases SHAPE_VALID')
    call assert_int_eq(int(CASSINI_ERROR_R_NONPOSITIVE), 100, 'r_nonpositive = 100')
    call assert_int_eq(int(CASSINI_ERROR_TOROIDAL), 101, 'toroidal = 101')
    call assert_int_eq(int(CASSINI_ERROR_NOT_STAR_CONVEX), 102, 'not_star_convex = 102')
    call assert_int_eq(int(CASSINI_ERROR_SEPARATED), 103, 'separated = 103')
    call assert_int_eq(int(CASSINI_MAX_PARAMS), 21, 'hard cap 21')
    call assert_int_eq(int(CASSINI_MAX_ORDER), 20, 'Legendre orders to 20')
    call assert_int_eq(int(CASSINI_N_GL_HALF), 64, 'GL order per half')
    call assert_int_eq(int(CASSINI_N_R_CHECK), 1025, 'R-check grid')
    call assert_int_eq(int(CASSINI_N_X_GRID_FLOOR), 65, 'scan floor')
    call assert_int_eq(int(CASSINI_N_POINTS_FLOOR), 4, 'profile floor')
    call assert_abs_close(STAR_CONVEXITY_MARGIN, 1.0e-3_rk, 0.0_rk, 'inclusive margin, provisional')
    call assert_true(len_trim(status_message(101_ik)) > 0, 'message for 101 nonempty')
    call assert_true(len_trim(status_message(SHAPE_ERROR_TOO_MANY_PARAMS)) > 0, 'shared code has message')
    call assert_true(len_trim(status_message(999_ik)) > 0, 'unknown code has fallback message')
    call test_summary()
end program cassini_param_status_test
