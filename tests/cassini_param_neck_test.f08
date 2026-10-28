!> The neck: absent on convex and slightly deformed ovals, the bare oval's
!! onset at eps = 1/2, exact symmetric anchors, asymmetric necks against a
!! dense scan, the tie rule, pinched and separated shapes, the refinement.
program cassini_param_neck_test
    use precision_utilities_mod, only: ik, rk
    use cassini_parameterization_workers_mod, only: resolve_t, tables_t, neck_s, cosine_nodes_s, &
            build_tables_s, curve_on_tables_s, series_at_s, lemniscate_point_s, CASSINI_MAX_ORDER
    use cassini_test_fixtures_mod, only: SPHERE, BARE_05, BARE_095, OBLATE, GS, SADDLE1, SD_MIN, &
            SADDLE2, ASYM, DIMPLED, FOLD_BOX, TWO_NECKS, STAR_FAIL, PINCHED, PINCHED_ASYM, SEPARATED
    use cassini_test_reference_mod, only: test_resolve_s
    use test_utils_mod, only: assert_true, assert_abs_close, assert_bits_eq, is_zero_f, test_summary
    implicit none

    call test_no_neck(SPHERE, 'SPHERE')
    call test_no_neck(BARE_05, 'BARE_05')
    call test_no_neck(OBLATE, 'OBLATE')
    call test_no_neck(GS, 'GS')
    call test_no_neck(SD_MIN, 'SD_MIN')
    call test_no_neck(SADDLE2, 'SADDLE2')
    call test_bare_onset()
    call test_symmetric_anchors()
    call test_asymmetric(STAR_FAIL, 0.561617_rk, 0.295547_rk, 'STAR_FAIL')
    call test_asymmetric(FOLD_BOX, -0.154052_rk, 0.865867_rk, 'FOLD_BOX')
    call test_two_necks()
    call test_rule_necks(PINCHED, 'PINCHED')
    call test_rule_necks(PINCHED_ASYM, 'PINCHED_ASYM')
    call test_rule_necks(SEPARATED, 'SEPARATED')
    call test_refinement(STAR_FAIL, 'STAR_FAIL')
    call test_refinement(FOLD_BOX, 'FOLD_BOX')
    call test_refinement(TWO_NECKS, 'TWO_NECKS')
    call test_ranked_after_refinement()
    call test_summary()

contains

    !> The neck on an n-node scan; also the scan's smallest interior local minimum.
    subroutine neck_on_scan_s(params, n, res, found, z_neck, rho_neck, x_neck, scan_min)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n
        type(resolve_t), intent(out) :: res
        logical, intent(out) :: found
        real(kind = rk), intent(out) :: z_neck, rho_neck, x_neck, scan_min
        type(tables_t) :: tab
        real(kind = rk) :: x(n), rho(n), z(n), drho_dx(n), dz_dx(n)
        integer(kind = ik) :: status, k
        call test_resolve_s(params, res, status)
        call cosine_nodes_s(n, x)
        call build_tables_s(x, CASSINI_MAX_ORDER, tab)
        call curve_on_tables_s(params, res, tab, rho, z, drho_dx, dz_dx)
        call neck_s(params, res, tab, rho, found, z_neck, rho_neck, x_neck)
        scan_min = huge(1.0_rk)
        do k = 2_ik, n - 1_ik
            if (rho(k) < rho(k - 1_ik) .and. rho(k) <= rho(k + 1_ik)) scan_min = min(scan_min, rho(k))
        end do
    end subroutine neck_on_scan_s

    !> N1: convex ovals and small deformations have no neck; outputs are 0.
    subroutine test_no_neck(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        logical :: found
        real(kind = rk) :: z_neck, rho_neck, x_neck, scan_min
        call neck_on_scan_s(params, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        call assert_true(.not. found, 'N1 no neck ' // name)
        call assert_true(is_zero_f(z_neck) .and. is_zero_f(rho_neck), 'N1 outputs zero ' // name)
    end subroutine test_no_neck

    !> N2: the bare oval acquires its neck at eps = 1/2.
    subroutine test_bare_onset()
        type(resolve_t) :: res
        logical :: found
        real(kind = rk) :: z_neck, rho_neck, x_neck, scan_min, bare(1)
        bare(1) = 0.49_rk
        call neck_on_scan_s(bare, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        call assert_true(.not. found, 'N2 [0.49] no neck')
        bare(1) = 0.51_rk
        call neck_on_scan_s(bare, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        call assert_true(found, 'N2 [0.51] neck')
        bare(1) = 0.5_rk
        call neck_on_scan_s(bare, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        write(*, '(A,L1)') 'N2 [0.5] found (recorded, not asserted): ', found
    end subroutine test_bare_onset

    !> N3: x = 0 is the neck of symmetric necked shapes, exactly.
    subroutine test_symmetric_anchors()
        type(resolve_t) :: res
        logical :: found
        real(kind = rk) :: z_neck, rho_neck, x_neck, scan_min
        call neck_on_scan_s(BARE_095, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        call assert_true(found, 'N3 BARE_095 found')
        call assert_bits_eq(z_neck, 0.0_rk, 'N3 BARE_095 z_neck = 0 bits')
        call assert_abs_close(rho_neck, sqrt(0.05_rk) / res%c, 1.0e-14_rk, 'N3 BARE_095 rho_neck')
        call neck_on_scan_s(DIMPLED, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        call assert_true(found, 'N3 DIMPLED found')
        call assert_abs_close(rho_neck, 0.687014_rk, 1.0e-5_rk, 'N3 DIMPLED rho_neck')
        call neck_on_scan_s(SADDLE1, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        call assert_true(found, 'N3 SADDLE1 found')
        call assert_abs_close(rho_neck, 0.827113_rk, 1.0e-5_rk, 'N3 SADDLE1 rho_neck')
        call neck_on_scan_s(ASYM, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        call assert_true(found, 'N3 ASYM found')
        call assert_abs_close(z_neck, -res%z_shift, 1.0e-14_rk, 'N3 ASYM z_neck = -z_shift')
        call assert_abs_close(rho_neck, 0.550778_rk, 1.0e-5_rk, 'N3 ASYM rho_neck')
    end subroutine test_symmetric_anchors

    !> N4: asymmetric necks, and against a 20001-node scan minimum.
    subroutine test_asymmetric(params, z_want, rho_want, name)
        real(kind = rk), intent(in) :: params(:), z_want, rho_want
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        logical :: found
        real(kind = rk) :: z_neck, rho_neck, x_neck, scan_min, z_d, rho_d, x_d, dense_min
        call neck_on_scan_s(params, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        call assert_true(found, 'N4 found ' // name)
        call assert_abs_close(z_neck, z_want, 1.0e-5_rk, 'N4 z_neck ' // name)
        call assert_abs_close(rho_neck, rho_want, 1.0e-5_rk, 'N4 rho_neck ' // name)
        call neck_on_scan_s(params, 20001_ik, res, found, z_d, rho_d, x_d, dense_min)
        call assert_abs_close(rho_neck, dense_min, 1.0e-6_rk, 'N4 rho_neck vs 20001-node minimum ' // name)
    end subroutine test_asymmetric

    !> N5: a central maximum between two equal necks; the tie goes north.
    subroutine test_two_necks()
        type(resolve_t) :: res
        logical :: found
        real(kind = rk) :: z_neck, rho_neck, x_neck, scan_min
        call neck_on_scan_s(TWO_NECKS, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        call assert_true(found, 'N5 TWO_NECKS found')
        call assert_abs_close(z_neck, 1.032412_rk, 1.0e-5_rk, 'N5 TWO_NECKS z_neck (north)')
        call assert_abs_close(rho_neck, 0.242723_rk, 1.0e-5_rk, 'N5 TWO_NECKS rho_neck')
    end subroutine test_two_necks

    !> N6: pinched and separated necks are set by rule.
    subroutine test_rule_necks(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        logical :: found
        real(kind = rk) :: z_neck, rho_neck, x_neck, scan_min
        call neck_on_scan_s(params, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        call assert_true(found, 'N6 found ' // name)
        call assert_true(is_zero_f(rho_neck), 'N6 rho_neck = 0 ' // name)
        call assert_bits_eq(z_neck, -res%z_shift, 'N6 z_neck = -z_shift ' // name)
    end subroutine test_rule_necks

    !> N7: the refined neck is no wider than the scan's and sits on drho/dx = 0.
    subroutine test_refinement(params, name)
        real(kind = rk), intent(in) :: params(:)
        character(len = *), intent(in) :: name
        type(resolve_t) :: res
        logical :: found
        real(kind = rk) :: z_neck, rho_neck, x_neck, scan_min, r, dr, rb, zb, drb, dzb
        call neck_on_scan_s(params, 129_ik, res, found, z_neck, rho_neck, x_neck, scan_min)
        call assert_true(rho_neck <= scan_min, 'N7 refined rho_neck <= scan minimum ' // name)
        call series_at_s(params, res%n_active, x_neck, r, dr)
        call lemniscate_point_s(x_neck, r, dr, res%eps, rb, zb, drb, dzb)
        call assert_true(abs(drb) <= 1.0e-12_rk, 'N7 drho_bar/dx = 0 at x_neck ' // name)
    end subroutine test_refinement

    !> N8: candidates are ranked by their REFINED rho. Two necks whose sampled
    !! rho differ by 2e-5 on the 65-node scan: the northern samples narrower,
    !! the southern refines narrower (0.386437 vs 0.386827). The neck is the
    !! southern one at every scan size (cross-family compliance review).
    subroutine test_ranked_after_refinement()
        real(kind = rk), parameter :: DOUBLE(5) = [0.657054_rk, 0.014275_rk, 0.009004_rk, &
                0.049565_rk, 0.498834_rk]
        integer(kind = ik), parameter :: GRIDS(3) = [65_ik, 129_ik, 2049_ik]
        type(resolve_t) :: res
        logical :: found
        real(kind = rk) :: z_neck, rho_neck, x_neck, scan_min
        integer(kind = ik) :: j
        do j = 1_ik, 3_ik
            call neck_on_scan_s(DOUBLE, GRIDS(j), res, found, z_neck, rho_neck, x_neck, scan_min)
            call assert_true(found, 'N8 found')
            call assert_abs_close(z_neck, -1.1220019_rk, 1.0e-6_rk, 'N8 the narrower (southern) neck')
            call assert_abs_close(rho_neck, 0.3864367_rk, 1.0e-6_rk, 'N8 rho_neck')
        end do
    end subroutine test_ranked_after_refinement

end program cassini_param_neck_test
