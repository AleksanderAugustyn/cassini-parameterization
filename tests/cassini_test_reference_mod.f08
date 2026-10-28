!> Independent reference computations for the cassini test suites: the
!! printed formulas evaluated literally, and integrations along paths the
!! library does not take.
module cassini_test_reference_mod

    use precision_utilities_mod, only: ik, rk

    implicit none

    private

    public :: reference_epsilon_f

contains

    !> J eq. 6 as printed: the pole brackets from the plain sums of alpha_n and
    !! (-1)^n alpha_n, the equator bracket from (2m - 1)!! / (2^m m!) built by a
    !! loop. params = [alpha, alpha1, ..., alphaN], any length >= 1.
    pure function reference_epsilon_f(params) result(eps)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk) :: eps
        real(kind = rk) :: alpha, sum_plus, sum_minus, sum_zero, coeff, sign_n
        integer(kind = ik) :: n, m, k, n_max

        alpha = params(1)
        n_max = size(params, kind = ik) - 1_ik
        sum_plus = 1.0_rk
        sum_minus = 1.0_rk
        sign_n = 1.0_rk
        do n = 1_ik, n_max
            sign_n = -sign_n
            sum_plus = sum_plus + params(n + 1_ik)
            sum_minus = sum_minus + sign_n * params(n + 1_ik)
        end do

        sum_zero = 1.0_rk
        do m = 1_ik, n_max / 2_ik
            coeff = 1.0_rk
            do k = 1_ik, m
                ! (2k - 1) / (2 k) builds (2m - 1)!! / (2^m m!)
                coeff = coeff * real(2_ik * k - 1_ik, rk) / real(2_ik * k, rk)
            end do
            if (mod(m, 2_ik) == 1_ik) coeff = -coeff
            sum_zero = sum_zero + coeff * params(2_ik * m + 1_ik)
        end do

        eps = (alpha - 1.0_rk) / 4.0_rk * (sum_plus**2 + sum_minus**2) &
                + (alpha + 1.0_rk) / 2.0_rk * sum_zero**2
    end function reference_epsilon_f

end module cassini_test_reference_mod
