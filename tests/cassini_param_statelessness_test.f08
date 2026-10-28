!> Contract family 2, statelessness: repeat calls, calls after rejections,
!! and independent caches all give the same bits; the cache owns its thetas.
program cassini_param_statelessness_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use cassini_parameterization_mod, only: cache_t
    use cassini_test_fixtures_mod, only: ASYM, R_NEG, TOROIDAL, SEPARATED, STAR_FAIL
    use cassini_test_tiers_mod, only: result_t, run_cached, build_cache, results_bits_eq
    use test_utils_mod, only: assert_true, test_summary
    implicit none

    integer(kind = ik), parameter :: N_X = 129_ik, N_Z = 257_ik
    real(kind = rk) :: thetas(32)
    integer(kind = ik) :: k

    do k = 1_ik, 32_ik
        thetas(k) = PI_C * (real(k, rk) - 0.5_rk) / 32.0_rk
    end do

    call test_repeat()
    call test_after_rejection(R_NEG, 'R_NEG')
    call test_after_rejection(TOROIDAL, 'TOROIDAL')
    call test_after_rejection(SEPARATED, 'SEPARATED')
    call test_after_rejection(STAR_FAIL, 'STAR_FAIL')
    call test_after_wrong_count()
    call test_two_caches_and_owned_thetas()
    call test_summary()

contains

    !> S1: the same call twice on one cache.
    subroutine test_repeat()
        type(cache_t) :: cache
        type(result_t) :: a, b
        integer(kind = ik) :: status
        call build_cache(5_ik, N_X, thetas, cache, status)
        call run_cached(cache, ASYM, thetas, N_Z, a)
        call run_cached(cache, ASYM, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'S1 repeat == first call')
    end subroutine test_repeat

    !> S2: a rejection leaves nothing behind.
    subroutine test_after_rejection(bad, name)
        real(kind = rk), intent(in) :: bad(:)
        character(len = *), intent(in) :: name
        type(cache_t) :: used, fresh
        type(result_t) :: a, b, junk
        integer(kind = ik) :: status
        call build_cache(5_ik, N_X, thetas, used, status)
        call build_cache(5_ik, N_X, thetas, fresh, status)
        call run_cached(used, bad, thetas, N_Z, junk)
        call run_cached(used, ASYM, thetas, N_Z, a)
        call run_cached(fresh, ASYM, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'S2 ASYM after ' // name // ' == fresh cache')
    end subroutine test_after_rejection

    !> S2: a wrong-count vector in between.
    subroutine test_after_wrong_count()
        type(cache_t) :: used, fresh
        type(result_t) :: a, b, junk
        real(kind = rk) :: long(6)
        integer(kind = ik) :: status
        long = 0.1_rk
        call build_cache(5_ik, N_X, thetas, used, status)
        call build_cache(5_ik, N_X, thetas, fresh, status)
        call run_cached(used, long, thetas, N_Z, junk)
        call run_cached(used, ASYM, thetas, N_Z, a)
        call run_cached(fresh, ASYM, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'S2 ASYM after a wrong count == fresh cache')
    end subroutine test_after_wrong_count

    !> S3: identical caches agree; overwriting the caller's thetas changes nothing.
    subroutine test_two_caches_and_owned_thetas()
        type(cache_t) :: c1, c2
        type(result_t) :: a, b, c
        real(kind = rk) :: mine(32)
        integer(kind = ik) :: status
        mine = thetas
        call build_cache(5_ik, N_X, mine, c1, status)
        call build_cache(5_ik, N_X, thetas, c2, status)
        call run_cached(c1, ASYM, thetas, N_Z, a)
        call run_cached(c2, ASYM, thetas, N_Z, b)
        call assert_true(results_bits_eq(a, b), 'S3 two caches agree')
        mine = 0.3_rk
        call run_cached(c1, ASYM, thetas, N_Z, c)
        call assert_true(results_bits_eq(a, c), 'S3 the cache owns a copy of its thetas')
    end subroutine test_two_caches_and_owned_thetas

end program cassini_param_statelessness_test
