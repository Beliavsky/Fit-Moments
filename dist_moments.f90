module dist_moments_mod
    use kind_mod
    implicit none
    private
    public :: nig_symmetric_moments, vg_symmetric_moments, gh_symmetric_moments, &
              normal_symmetric_moments, student_t_symmetric_moments, &
              logistic_symmetric_moments, sech_symmetric_moments, &
              laplace_symmetric_moments, ged_symmetric_moments, &
              cauchy_symmetric_moments

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: catalan = 0.91596559417721901505_dp

contains

    subroutine logistic_symmetric_moments(s, m2, m4, m1_abs)
        real(dp), intent(in)  :: s
        real(dp), intent(out) :: m2, m4
        real(dp), intent(out), optional :: m1_abs
        m2 = (pi * s)**2 / 3.0_dp
        m4 = 4.2_dp * m2**2
        if (present(m1_abs)) m1_abs = 2.0_dp * s * log(2.0_dp)
    end subroutine logistic_symmetric_moments

    subroutine sech_symmetric_moments(s, m2, m4, m1_abs)
        real(dp), intent(in)  :: s
        real(dp), intent(out) :: m2, m4
        real(dp), intent(out), optional :: m1_abs
        m2 = s**2
        m4 = 5.0_dp * m2**2
        if (present(m1_abs)) m1_abs = 4.0_dp * s * catalan / pi
    end subroutine sech_symmetric_moments

    subroutine laplace_symmetric_moments(b, m2, m4, m1_abs)
        real(dp), intent(in)  :: b
        real(dp), intent(out) :: m2, m4
        real(dp), intent(out), optional :: m1_abs
        m2 = 2.0_dp * b**2
        m4 = 6.0_dp * m2**2
        if (present(m1_abs)) m1_abs = b
    end subroutine laplace_symmetric_moments

    subroutine ged_symmetric_moments(beta, s, m2, m4, m1_abs)
        real(dp), intent(in)  :: beta, s
        real(dp), intent(out) :: m2, m4
        real(dp), intent(out), optional :: m1_abs
        real(dp) :: g1, g2, g3, g5
        g1 = exp(log_gamma(1.0_dp / beta))
        g3 = exp(log_gamma(3.0_dp / beta))
        g5 = exp(log_gamma(5.0_dp / beta))
        m2 = s**2 * (g3 / g1)
        m4 = s**4 * (g5 / g1)
        if (present(m1_abs)) then
            g2 = exp(log_gamma(2.0_dp / beta))
            m1_abs = s * (g2 / g1)
        end if
    end subroutine ged_symmetric_moments

    subroutine cauchy_symmetric_moments(gamma, m2, m4, m1_abs)
        real(dp), intent(in)  :: gamma
        real(dp), intent(out) :: m2, m4
        real(dp), intent(out), optional :: m1_abs
        m2 = -1.0_dp
        m4 = -1.0_dp
        if (present(m1_abs)) m1_abs = -1.0_dp
    end subroutine cauchy_symmetric_moments

    subroutine normal_symmetric_moments(mu, sigma, m2, m4, m1_abs)
        real(dp), intent(in)  :: mu, sigma
        real(dp), intent(out) :: m2, m4
        real(dp), intent(out), optional :: m1_abs
        m2 = sigma**2
        m4 = 3.0_dp * sigma**4
        if (present(m1_abs)) m1_abs = sigma * sqrt(2.0_dp / pi)
    end subroutine normal_symmetric_moments

    subroutine student_t_symmetric_moments(df, m2, m4, s, m1_abs)
        real(dp), intent(in)  :: df
        real(dp), intent(out) :: m2, m4
        real(dp), intent(in), optional :: s
        real(dp), intent(out), optional :: m1_abs
        real(dp) :: ls
        ls = 1.0_dp
        if (present(s)) ls = s
        m2 = ls**2 * df / (df - 2.0_dp)
        m4 = (3.0_dp + 6.0_dp / (df - 4.0_dp)) * m2**2
        if (present(m1_abs)) then
            if (df > 1.0_dp) then
                m1_abs = ls * (2.0_dp * sqrt(df) * exp(log_gamma((df + 1.0_dp) / 2.0_dp))) / &
                         ((df - 1.0_dp) * sqrt(pi) * exp(log_gamma(df / 2.0_dp)))
            else
                m1_abs = -1.0_dp
            end if
        end if
    end subroutine student_t_symmetric_moments

    subroutine nig_symmetric_moments(alpha, delta, m2, m4, m1_abs)
        real(dp), intent(in)  :: alpha, delta
        real(dp), intent(out) :: m2, m4
        real(dp), intent(out), optional :: m1_abs
        m2 = delta / alpha
        m4 = 3.0_dp * (delta / alpha**3 + (delta / alpha)**2)
        if (present(m1_abs)) then
            m1_abs = (2.0_dp * delta * exp(alpha * delta) * kv(0.0_dp, alpha * delta)) / pi
        end if
    end subroutine nig_symmetric_moments

    subroutine vg_symmetric_moments(lambda, alpha, m2, m4, m1_abs)
        real(dp), intent(in)  :: lambda, alpha
        real(dp), intent(out) :: m2, m4
        real(dp), intent(out), optional :: m1_abs
        m2 = 2.0_dp * lambda / alpha**2
        m4 = 12.0_dp * (lambda**2 + lambda) / alpha**4
        if (present(m1_abs)) then
            m1_abs = (2.0_dp * exp(log_gamma(lambda + 0.5_dp))) / &
                     (alpha * sqrt(pi) * exp(log_gamma(lambda)))
        end if
    end subroutine vg_symmetric_moments

    subroutine gh_symmetric_moments(lambda, alpha, delta, m2, m4, m1_abs)
        real(dp), intent(in)  :: lambda, alpha, delta
        real(dp), intent(out) :: m2, m4
        real(dp), intent(out), optional :: m1_abs
        real(dp) :: k_lam, k_lam1, k_lam2, ad

        ad = alpha * delta
        k_lam = kv(lambda, ad)
        k_lam1 = kv(lambda + 1.0_dp, ad)
        k_lam2 = kv(lambda + 2.0_dp, ad)
        
        m2 = (delta / alpha) * (k_lam1 / k_lam)
        m4 = 3.0_dp * (delta / alpha)**2 * (k_lam2 / k_lam)
        if (present(m1_abs)) then
            m1_abs = sqrt(2.0_dp / pi) * sqrt(delta / alpha) * kv(lambda + 0.5_dp, ad) / k_lam
        end if
    end subroutine gh_symmetric_moments

    function kv(v, x) result(res)
        real(dp), intent(in) :: v, x
        real(dp) :: res
        integer, parameter :: m = 2000
        real(dp) :: t_max, dt, t, sum_val
        integer :: i
        
        t_max = acosh( (50.0_dp + abs(v) * log(10.0_dp)) / x + 1.0_dp )
        if (t_max < 10.0_dp) t_max = 10.0_dp
        if (t_max > 20.0_dp) t_max = 20.0_dp
        
        dt = t_max / dble(m)
        sum_val = 0.5_dp * (exp(-x * cosh(0.0_dp)) * cosh(v * 0.0_dp) + &
                            exp(-x * cosh(t_max)) * cosh(v * t_max))
        do i = 1, m - 1
            t = dble(i) * dt
            sum_val = sum_val + exp(-x * cosh(t)) * cosh(v * t)
        end do
        res = sum_val * dt
    end function kv

end module dist_moments_mod
