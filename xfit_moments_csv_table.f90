program xfit_moments_csv_table
    use kind_mod
    use dist_moments_mod
    implicit none

    real(dp), parameter :: missing_param = -999.0_dp
    real(dp), parameter :: bad_dev = 1.0e99_dp
    character(len=*), parameter :: missing_name = '-'
    integer, parameter :: n_dist = 10
    character(len=*), parameter :: csv_file = 'temp.csv'
    character(len=*), parameter :: csv_header = &
        'rank,distribution,npar,scale_name,scale,shape1_name,shape1,shape2_name,shape2,' // &
        'm2_fit,m4_fit,dm2,dm4,rel_dev,note,m2_target,m4_target,kurt_target'

    type :: fit_row
        character(len=24) :: name = ''
        integer :: npar = 0
        character(len=12) :: scale_name = missing_name
        character(len=12) :: shape1_name = missing_name
        character(len=12) :: shape2_name = missing_name
        real(dp) :: scale = missing_param
        real(dp) :: shape1 = missing_param
        real(dp) :: shape2 = missing_param
        real(dp) :: m2_fit = 0.0_dp
        real(dp) :: m4_fit = 0.0_dp
        real(dp) :: dm2 = 0.0_dp
        real(dp) :: dm4 = 0.0_dp
        real(dp) :: dev_score = 0.0_dp
        logical :: ok = .false.
        character(len=80) :: note = ''
    end type fit_row

    type(fit_row) :: rows(n_dist)
    real(dp) :: m2_target, m4_target, kurt_target, gh_lambda

    print *, 'enter target m2 and m4:'
    read(*, *) m2_target, m4_target

    if (m2_target <= 0.0_dp) then
        error stop 'm2 must be positive'
    end if
    if (m4_target <= 0.0_dp) then
        error stop 'm4 must be positive'
    end if

    kurt_target = m4_target / m2_target**2
    gh_lambda = 1.0_dp

    call add_normal(rows(1), m2_target, m4_target)
    call add_logistic(rows(2), m2_target, m4_target)
    call add_sech(rows(3), m2_target, m4_target)
    call add_laplace(rows(4), m2_target, m4_target)
    call add_student_t(rows(5), m2_target, m4_target)
    call add_nig(rows(6), m2_target, m4_target)
    call add_vg(rows(7), m2_target, m4_target)
    call add_ged(rows(8), m2_target, m4_target)
    call add_gh(rows(9), m2_target, m4_target, gh_lambda)
    call add_cauchy(rows(10))

    call sort_rows(rows)
    call write_csv(rows, m2_target, m4_target, kurt_target, csv_file)

    print '(a,es16.8)', 'target m2   = ', m2_target
    print '(a,es16.8)', 'target m4   = ', m4_target
    print '(a,es16.8)', 'target kurt = ', kurt_target
    print '(a)', 'wrote csv table to temp.csv'

contains

    subroutine add_normal(row, m2_target, m4_target)
        type(fit_row), intent(out) :: row
        real(dp), intent(in) :: m2_target, m4_target
        real(dp) :: sigma, m2_fit, m4_fit

        sigma = sqrt(m2_target)
        call normal_symmetric_moments(0.0_dp, sigma, m2_fit, m4_fit)
        call set_row(row, 'normal', 1, 'sigma', sigma, missing_name, missing_param, missing_name, missing_param, &
            m2_target, m4_target, m2_fit, m4_fit, .true., 'fit to m2 only')
    end subroutine add_normal

    subroutine add_logistic(row, m2_target, m4_target)
        type(fit_row), intent(out) :: row
        real(dp), intent(in) :: m2_target, m4_target
        real(dp) :: s, m2_fit, m4_fit
        real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp

        s = sqrt(3.0_dp * m2_target) / pi
        call logistic_symmetric_moments(s, m2_fit, m4_fit)
        call set_row(row, 'logistic', 1, 's', s, missing_name, missing_param, missing_name, missing_param, &
            m2_target, m4_target, m2_fit, m4_fit, .true., 'fit to m2 only')
    end subroutine add_logistic

    subroutine add_sech(row, m2_target, m4_target)
        type(fit_row), intent(out) :: row
        real(dp), intent(in) :: m2_target, m4_target
        real(dp) :: s, m2_fit, m4_fit

        s = sqrt(m2_target)
        call sech_symmetric_moments(s, m2_fit, m4_fit)
        call set_row(row, 'sech', 1, 's', s, missing_name, missing_param, missing_name, missing_param, &
            m2_target, m4_target, m2_fit, m4_fit, .true., 'fit to m2 only')
    end subroutine add_sech

    subroutine add_laplace(row, m2_target, m4_target)
        type(fit_row), intent(out) :: row
        real(dp), intent(in) :: m2_target, m4_target
        real(dp) :: b, m2_fit, m4_fit

        b = sqrt(0.5_dp * m2_target)
        call laplace_symmetric_moments(b, m2_fit, m4_fit)
        call set_row(row, 'laplace', 1, 'b', b, missing_name, missing_param, missing_name, missing_param, &
            m2_target, m4_target, m2_fit, m4_fit, .true., 'fit to m2 only')
    end subroutine add_laplace

    subroutine add_student_t(row, m2_target, m4_target)
        type(fit_row), intent(out) :: row
        real(dp), intent(in) :: m2_target, m4_target
        real(dp) :: kurt, df, scale, m2_std, m4_std, m2_fit, m4_fit
        logical :: ok
        character(len=80) :: note

        kurt = m4_target / m2_target**2
        ok = .true.
        note = 'exact m2,m4 fit'
        df = missing_param
        scale = missing_param
        m2_fit = 0.0_dp
        m4_fit = 0.0_dp

        if (kurt <= 3.0_dp) then
            ok = .false.
            note = 'requires kurtosis > 3'
        else
            df = 4.0_dp + 6.0_dp / (kurt - 3.0_dp)
            call student_t_symmetric_moments(df, m2_std, m4_std)
            scale = sqrt(m2_target / m2_std)
            m2_fit = scale**2 * m2_std
            m4_fit = scale**4 * m4_std
        end if

        call set_row(row, 'student_t', 2, 'scale', scale, 'df', df, missing_name, missing_param, &
            m2_target, m4_target, m2_fit, m4_fit, ok, note)
    end subroutine add_student_t

    subroutine add_nig(row, m2_target, m4_target)
        type(fit_row), intent(out) :: row
        real(dp), intent(in) :: m2_target, m4_target
        real(dp) :: kurt, p, alpha, delta, m2_fit, m4_fit
        logical :: ok
        character(len=80) :: note

        kurt = m4_target / m2_target**2
        ok = .true.
        note = 'exact m2,m4 fit'
        alpha = missing_param
        delta = missing_param
        m2_fit = 0.0_dp
        m4_fit = 0.0_dp

        if (kurt <= 3.0_dp) then
            ok = .false.
            note = 'requires kurtosis > 3'
        else
            p = 3.0_dp / (kurt - 3.0_dp)
            alpha = sqrt(p / m2_target)
            delta = sqrt(p * m2_target)
            call nig_symmetric_moments(alpha, delta, m2_fit, m4_fit)
        end if

        call set_row(row, 'nig', 2, 'delta', delta, 'alpha', alpha, missing_name, missing_param, &
            m2_target, m4_target, m2_fit, m4_fit, ok, note)
    end subroutine add_nig

    subroutine add_vg(row, m2_target, m4_target)
        type(fit_row), intent(out) :: row
        real(dp), intent(in) :: m2_target, m4_target
        real(dp) :: kurt, lambda, alpha, m2_fit, m4_fit
        logical :: ok
        character(len=80) :: note

        kurt = m4_target / m2_target**2
        ok = .true.
        note = 'exact m2,m4 fit'
        lambda = missing_param
        alpha = missing_param
        m2_fit = 0.0_dp
        m4_fit = 0.0_dp

        if (kurt <= 3.0_dp) then
            ok = .false.
            note = 'requires kurtosis > 3'
        else
            lambda = 3.0_dp / (kurt - 3.0_dp)
            alpha = sqrt(2.0_dp * lambda / m2_target)
            call vg_symmetric_moments(lambda, alpha, m2_fit, m4_fit)
        end if

        call set_row(row, 'vg', 2, 'alpha', alpha, 'lambda', lambda, missing_name, missing_param, &
            m2_target, m4_target, m2_fit, m4_fit, ok, note)
    end subroutine add_vg

    subroutine add_ged(row, m2_target, m4_target)
        type(fit_row), intent(out) :: row
        real(dp), intent(in) :: m2_target, m4_target
        real(dp) :: kurt, beta, s, g1, g3, m2_fit, m4_fit
        logical :: ok
        character(len=80) :: note

        kurt = m4_target / m2_target**2
        beta = missing_param
        s = missing_param
        m2_fit = 0.0_dp
        m4_fit = 0.0_dp
        call solve_ged_beta(kurt, beta, ok, note)

        if (ok) then
            g1 = exp(log_gamma(1.0_dp / beta))
            g3 = exp(log_gamma(3.0_dp / beta))
            s = sqrt(m2_target * g1 / g3)
            call ged_symmetric_moments(beta, s, m2_fit, m4_fit)
            note = 'exact m2,m4 fit'
        else
            beta = missing_param
            s = missing_param
        end if

        call set_row(row, 'ged', 2, 's', s, 'beta', beta, missing_name, missing_param, &
            m2_target, m4_target, m2_fit, m4_fit, ok, note)
    end subroutine add_ged

    subroutine add_gh(row, m2_target, m4_target, lambda)
        type(fit_row), intent(out) :: row
        real(dp), intent(in) :: m2_target, m4_target, lambda
        real(dp) :: kurt, z, q, r1, alpha, delta, m2_fit, m4_fit
        logical :: ok
        character(len=80) :: note

        kurt = m4_target / m2_target**2
        alpha = missing_param
        delta = missing_param
        m2_fit = 0.0_dp
        m4_fit = 0.0_dp
        call solve_gh_z(lambda, kurt, z, ok, note)

        if (ok) then
            r1 = kv_local(lambda + 1.0_dp, z) / kv_local(lambda, z)
            q = m2_target / r1
            alpha = sqrt(z / q)
            delta = sqrt(z * q)
            call gh_symmetric_moments(lambda, alpha, delta, m2_fit, m4_fit)
            note = 'exact m2,m4 fit'
        else
            alpha = missing_param
            delta = missing_param
        end if

        call set_row(row, 'gh_fixed_lambda', 3, 'delta', delta, 'lambda', lambda, 'alpha', alpha, &
            m2_target, m4_target, m2_fit, m4_fit, ok, note)
    end subroutine add_gh

    subroutine add_cauchy(row)
        type(fit_row), intent(out) :: row

        row%name = 'cauchy'
        row%npar = 1
        row%scale_name = 'gamma'
        row%scale = missing_param
        row%shape1_name = missing_name
        row%shape1 = missing_param
        row%shape2_name = missing_name
        row%shape2 = missing_param
        row%m2_fit = 0.0_dp
        row%m4_fit = 0.0_dp
        row%dm2 = bad_dev
        row%dm4 = bad_dev
        row%dev_score = bad_dev
        row%ok = .false.
        row%note = 'm2 and m4 do not exist'
    end subroutine add_cauchy

    subroutine set_row(row, name, npar, scale_name, scale, shape1_name, shape1, shape2_name, shape2, &
        m2_target, m4_target, m2_fit, m4_fit, ok, note)
        type(fit_row), intent(out) :: row
        character(len=*), intent(in) :: name, scale_name, shape1_name, shape2_name, note
        integer, intent(in) :: npar
        real(dp), intent(in) :: scale, shape1, shape2, m2_target, m4_target, m2_fit, m4_fit
        logical, intent(in) :: ok
        real(dp) :: scale2, scale4

        row%name = name
        row%npar = npar
        row%scale_name = scale_name
        row%shape1_name = shape1_name
        row%shape2_name = shape2_name
        row%scale = scale
        row%shape1 = shape1
        row%shape2 = shape2
        row%m2_fit = m2_fit
        row%m4_fit = m4_fit
        row%ok = ok
        row%note = note

        if (ok) then
            row%dm2 = m2_fit - m2_target
            row%dm4 = m4_fit - m4_target
            scale2 = max(1.0_dp, abs(m2_target))
            scale4 = max(1.0_dp, abs(m4_target))
            row%dev_score = abs(row%dm2) / scale2 + abs(row%dm4) / scale4
        else
            row%dm2 = bad_dev
            row%dm4 = bad_dev
            row%dev_score = bad_dev
        end if
    end subroutine set_row

    subroutine sort_rows(rows)
        type(fit_row), intent(inout) :: rows(:)
        type(fit_row) :: tmp
        integer :: i, j, k

        do i = 1, size(rows) - 1
            k = i
            do j = i + 1, size(rows)
                if (rows(j)%dev_score < rows(k)%dev_score) k = j
            end do
            if (k /= i) then
                tmp = rows(i)
                rows(i) = rows(k)
                rows(k) = tmp
            end if
        end do
    end subroutine sort_rows

    subroutine write_csv(rows, m2_target, m4_target, kurt_target, filename)
        type(fit_row), intent(in) :: rows(:)
        real(dp), intent(in) :: m2_target, m4_target, kurt_target
        character(len=*), intent(in) :: filename
        integer :: i, unit

        open(newunit=unit, file=filename, status='replace', action='write')
        write(unit, '(a)') csv_header
        do i = 1, size(rows)
            call write_csv_row(unit, i, rows(i), m2_target, m4_target, kurt_target)
        end do
        close(unit)
    end subroutine write_csv

    subroutine write_csv_row(unit, rank, row, m2_target, m4_target, kurt_target)
        integer, intent(in) :: unit, rank
        type(fit_row), intent(in) :: row
        real(dp), intent(in) :: m2_target, m4_target, kurt_target
        character(len=32) :: srank, snpar
        character(len=:), allocatable :: f_distribution, f_scale_name, f_shape1_name, f_shape2_name
        character(len=:), allocatable :: f_scale, f_shape1, f_shape2
        character(len=:), allocatable :: f_m2f, f_m4f, f_dm2, f_dm4, f_dev
        character(len=:), allocatable :: f_note, f_m2t, f_m4t, f_kurt

        write(srank, '(i0)') rank
        write(snpar, '(i0)') row%npar

        f_distribution = csv_text_field(row%name)
        f_scale_name = csv_text_field(row%scale_name)
        f_shape1_name = csv_text_field(row%shape1_name)
        f_shape2_name = csv_text_field(row%shape2_name)
        f_scale = csv_real_field(row%scale)
        f_shape1 = csv_real_field(row%shape1)
        f_shape2 = csv_real_field(row%shape2)
        f_m2f = csv_real_field(row%m2_fit)
        f_m4f = csv_real_field(row%m4_fit)
        f_dm2 = csv_real_field(row%dm2)
        f_dm4 = csv_real_field(row%dm4)
        f_dev = csv_real_field(row%dev_score)
        f_note = csv_note_field(row%note)
        f_m2t = csv_real_field(m2_target)
        f_m4t = csv_real_field(m4_target)
        f_kurt = csv_real_field(kurt_target)

        write(unit, '(a)') trim(srank) // ',' // f_distribution // ',' // trim(snpar) // ',' // &
            f_scale_name // ',' // f_scale // ',' // &
            f_shape1_name // ',' // f_shape1 // ',' // &
            f_shape2_name // ',' // f_shape2 // ',' // &
            f_m2f // ',' // f_m4f // ',' // f_dm2 // ',' // f_dm4 // ',' // f_dev // ',' // &
            f_note // ',' // f_m2t // ',' // f_m4t // ',' // f_kurt
    end subroutine write_csv_row

    function csv_quote(s) result(out)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: out

        out = '"' // trim(s) // '"'
    end function csv_quote

    function csv_text_field(s) result(out)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: out

        if (trim(s) == '' .or. trim(s) == missing_name) then
            out = ''
        else
            out = csv_quote(trim(s))
        end if
    end function csv_text_field

    function csv_note_field(s) result(out)
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: out

        if (trim(s) == '') then
            out = ''
        else
            out = csv_quote(trim(s))
        end if
    end function csv_note_field

    function csv_real_field(x) result(out)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: out
        character(len=64) :: sx

        if (abs(x - missing_param) <= 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(missing_param))) then
            out = ''
        else
            write(sx, '(es24.16)') x
            out = trim(adjustl(sx))
        end if
    end function csv_real_field

    subroutine solve_ged_beta(kurt_target, beta, ok, note)
        real(dp), intent(in) :: kurt_target
        real(dp), intent(out) :: beta
        logical, intent(out) :: ok
        character(len=*), intent(out) :: note
        real(dp) :: lo, hi, mid, flo, fhi, fmid
        integer :: iter

        lo = 0.10_dp
        hi = 200.0_dp
        flo = ged_kurtosis(lo)
        fhi = ged_kurtosis(hi)

        if (kurt_target > flo) then
            ok = .false.
            beta = missing_param
            note = 'target kurtosis above GED bracket'
            return
        end if

        if (kurt_target < fhi) then
            ok = .false.
            beta = missing_param
            note = 'target kurtosis below GED range'
            return
        end if

        do iter = 1, 200
            mid = 0.5_dp * (lo + hi)
            fmid = ged_kurtosis(mid)
            if (abs(fmid - kurt_target) <= 1.0e-12_dp * max(1.0_dp, kurt_target)) exit
            if (fmid > kurt_target) then
                lo = mid
            else
                hi = mid
            end if
        end do

        beta = 0.5_dp * (lo + hi)
        ok = .true.
        note = ''
    end subroutine solve_ged_beta

    function ged_kurtosis(beta) result(kurt)
        real(dp), intent(in) :: beta
        real(dp) :: kurt
        real(dp) :: lg1, lg3, lg5

        lg1 = log_gamma(1.0_dp / beta)
        lg3 = log_gamma(3.0_dp / beta)
        lg5 = log_gamma(5.0_dp / beta)
        kurt = exp(lg5 + lg1 - 2.0_dp * lg3)
    end function ged_kurtosis

    subroutine solve_gh_z(lambda, kurt_target, z, ok, note)
        real(dp), intent(in) :: lambda, kurt_target
        real(dp), intent(out) :: z
        logical, intent(out) :: ok
        character(len=*), intent(out) :: note
        real(dp) :: lo, hi, mid, flo, fhi, fmid
        integer :: iter

        lo = 1.0e-4_dp
        hi = 200.0_dp
        flo = gh_kurtosis_fixed_lambda(lambda, lo)
        fhi = gh_kurtosis_fixed_lambda(lambda, hi)

        if (kurt_target > flo) then
            ok = .false.
            z = 0.0_dp
            note = 'target kurtosis too large for fixed lambda'
            return
        end if

        if (kurt_target < fhi) then
            ok = .false.
            z = 0.0_dp
            note = 'target kurtosis too close to 3'
            return
        end if

        do iter = 1, 200
            mid = 0.5_dp * (lo + hi)
            fmid = gh_kurtosis_fixed_lambda(lambda, mid)
            if (abs(fmid - kurt_target) <= 1.0e-11_dp * max(1.0_dp, kurt_target)) exit
            if (fmid > kurt_target) then
                lo = mid
            else
                hi = mid
            end if
        end do

        z = 0.5_dp * (lo + hi)
        ok = .true.
        note = ''
    end subroutine solve_gh_z

    function gh_kurtosis_fixed_lambda(lambda, z) result(kurt)
        real(dp), intent(in) :: lambda, z
        real(dp) :: kurt
        real(dp) :: k0, k1, k2

        k0 = kv_local(lambda, z)
        k1 = kv_local(lambda + 1.0_dp, z)
        k2 = kv_local(lambda + 2.0_dp, z)
        kurt = 3.0_dp * k2 * k0 / (k1**2)
    end function gh_kurtosis_fixed_lambda

    function kv_local(v, x) result(res)
        real(dp), intent(in) :: v, x
        real(dp) :: res
        integer, parameter :: m = 2000
        real(dp) :: t_max, dt, t, sum_val
        integer :: i

        t_max = acosh((50.0_dp + abs(v) * log(10.0_dp)) / x + 1.0_dp)
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
    end function kv_local

end program xfit_moments_csv_table
