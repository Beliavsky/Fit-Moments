program xdist_moments
    use kind_mod
    use dist_moments_mod
    implicit none
    real(dp) :: alpha, delta, lambda, mu, sigma, df, m2, m4

    print *, "Testing Normal symmetric moments (mu=0):"
    print *
    mu = 0.0_dp
    sigma = 1.0_dp
    call normal_symmetric_moments(mu, sigma, m2, m4)
    print *, "mu = ", mu, " sigma = ", sigma
    print *, "  m2 (variance)   = ", m2
    print *, "  m4              = ", m4
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp
    print *

    print *, "Testing Student's t symmetric moments (mu=0):"
    print *
    df = 10.0_dp
    call student_t_symmetric_moments(df, m2, m4)
    print *, "df = ", df
    print *, "  m2 (variance)   = ", m2
    print *, "  m4              = ", m4
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp
    print *

    print *, "Testing NIG symmetric moments (mu=0, beta=0):"
    print *
    
    alpha = 1.0_dp
    delta = 1.0_dp
    call nig_symmetric_moments(alpha, delta, m2, m4)
    print *, "alpha = ", alpha, " delta = ", delta
    print *, "  m2 (variance)   = ", m2
    print *, "  m4              = ", m4
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp
    print *

    alpha = 2.0_dp
    delta = 0.5_dp
    call nig_symmetric_moments(alpha, delta, m2, m4)
    print *, "alpha = ", alpha, " delta = ", delta
    print *, "  m2 (variance)   = ", m2
    print *, "  m4              = ", m4
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp
    print *

    print *, "Testing VG symmetric moments (mu=0, beta=0):"
    print *
    
    alpha = 1.0_dp
    lambda = 1.0_dp
    call vg_symmetric_moments(lambda, alpha, m2, m4)
    print *, "alpha = ", alpha, " lambda = ", lambda
    print *, "  m2 (variance)   = ", m2
    print *, "  m4              = ", m4
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp
    print *

    alpha = 2.0_dp
    lambda = 0.5_dp
    call vg_symmetric_moments(lambda, alpha, m2, m4)
    print *, "alpha = ", alpha, " lambda = ", lambda
    print *, "  m2 (variance)   = ", m2
    print *, "  m4              = ", m4
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp
    print *

    print *, "Testing GH symmetric moments (mu=0, beta=0):"
    print *

    alpha = 1.0_dp
    delta = 1.0_dp
    lambda = 1.0_dp
    call gh_symmetric_moments(lambda, alpha, delta, m2, m4)
    print *, "alpha = ", alpha, " delta = ", delta, " lambda = ", lambda
    print *, "  m2 (variance)   = ", m2
    print *, "  m4              = ", m4
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp
    print *

    alpha = 2.0_dp
    delta = 0.5_dp
    lambda = 2.0_dp
    call gh_symmetric_moments(lambda, alpha, delta, m2, m4)
    print *, "alpha = ", alpha, " delta = ", delta, " lambda = ", lambda
    print *, "  m2 (variance)   = ", m2
    print *, "  m4              = ", m4
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp
    print *

    print *, "Testing Logistic symmetric moments:"
    print *
    call logistic_symmetric_moments(1.0_dp, m2, m4)
    print *, "s = 1.0"
    print *, "  m2 (variance)   = ", m2
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp
    print *

    print *, "Testing Sech symmetric moments:"
    print *
    call sech_symmetric_moments(1.0_dp, m2, m4)
    print *, "s = 1.0"
    print *, "  m2 (variance)   = ", m2
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp
    print *

    print *, "Testing Laplace symmetric moments:"
    print *
    call laplace_symmetric_moments(1.0_dp, m2, m4)
    print *, "b = 1.0"
    print *, "  m2 (variance)   = ", m2
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp
    print *

    print *, "Testing GED symmetric moments:"
    print *
    call ged_symmetric_moments(1.5_dp, 1.0_dp, m2, m4)
    print *, "beta = 1.5, s = 1.0"
    print *, "  m2 (variance)   = ", m2
    print *, "  Excess Kurtosis = ", m4/(m2**2) - 3.0_dp

end program xdist_moments
