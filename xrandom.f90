program xrandom
    use random_mod
    implicit none
    integer, parameter :: n = 1000000
    real(dp), dimension(n) :: v
    real(dp) :: alpha, delta, lambda

    call random_seed()

    ! --- Normal ---
    v = rnorm(n, 0.0_dp, 1.0_dp)
    call print_moments("Normal (0,1)", v)

    ! --- Student t ---
    v = rt(n, 10.0_dp)
    call print_moments("Student t (df=10)", v)

    ! --- Symmetric NIG ---
    alpha = 1.0_dp
    delta = 1.0_dp
    v = rnig(n, alpha, delta)
    call print_moments("Symmetric NIG (alpha=1, delta=1)", v)

    ! --- Symmetric VG ---
    lambda = 1.0_dp
    alpha = 1.0_dp
    v = rvg(n, lambda, alpha)
    call print_moments("Symmetric VG (lambda=1, alpha=1)", v)

    ! --- Symmetric GH ---
    lambda = 1.0_dp
    alpha = 1.0_dp
    delta = 1.0_dp
    v = rgh(n, lambda, alpha, delta)
    call print_moments("Symmetric GH (lambda=1, alpha=1, delta=1)", v)

    ! --- Logistic ---
    v = rlogistic(n, 0.0_dp, 1.0_dp)
    call print_moments("Logistic (mu=0, s=1)", v)

    ! --- Sech ---
    v = rsech(n, 0.0_dp, 1.0_dp)
    call print_moments("Sech (mu=0, s=1)", v)

    ! --- Laplace ---
    v = rlaplace(n, 0.0_dp, 1.0_dp)
    call print_moments("Laplace (mu=0, b=1)", v)

    ! --- GED ---
    v = rged(n, 0.0_dp, 1.0_dp, 1.5_dp)
    call print_moments("GED (mu=0, s=1, beta=1.5)", v)

    ! --- Cauchy ---
    v = rcauchy(n, 0.0_dp, 1.0_dp)
    call print_moments("Cauchy (x0=0, gamma=1)", v)

contains

    subroutine print_moments(label, x)
        character(len=*), intent(in) :: label
        real(dp), dimension(:), intent(in) :: x
        real(dp) :: m, m2, m4, k
        integer :: sz
        sz = size(x)
        m = sum(x) / sz
        m2 = sum((x - m)**2) / sz
        m4 = sum((x - m)**4) / sz
        k = (m4 / (m2**2)) - 3.0_dp
        print *, label
        print *, "  Mean:            ", m
        print *, "  Variance:        ", m2
        print *, "  Excess Kurtosis: ", k
        print *
    end subroutine print_moments

end program xrandom
