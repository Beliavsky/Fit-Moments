module fit_moments_families_mod
use, intrinsic :: iso_fortran_env, only: dp => real64
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
implicit none
private

real(kind=dp), parameter :: tiny = 1.0e-12_dp
real(kind=dp), parameter :: ged_kurt_min = 1.8_dp
real(kind=dp), parameter :: huge_obj = 1.0e300_dp

public :: moment_summary_t, fit_result_t
public :: validate_target_moments, print_fit_result
public :: student_t_fit_from_kurtosis, fit_jf, fit_skew_t_fs, fit_azzalini_skew_t
public :: fit_ged_from_kurtosis, fit_skew_ged
public :: fit_symmetric_nig_from_kurtosis, fit_nig
public :: fit_symmetric_vg_from_kurtosis, fit_vg
public :: summary_from_mean_var_skew_kurt, central_stats_from_raw

abstract interface
   subroutine residual_2d_interface(x, r, ok)
      import :: dp
      real(kind=dp), intent(in) :: x(2)
      real(kind=dp), intent(out) :: r(2)
      logical, intent(out) :: ok
   end subroutine residual_2d_interface
end interface

type :: moment_summary_t
   real(kind=dp) :: m1 = 0.0_dp
   real(kind=dp) :: m2 = 0.0_dp
   real(kind=dp) :: m3 = 0.0_dp
   real(kind=dp) :: m4 = 0.0_dp
   real(kind=dp) :: mean = 0.0_dp
   real(kind=dp) :: var = 0.0_dp
   real(kind=dp) :: skew = 0.0_dp
   real(kind=dp) :: kurt = 0.0_dp
   real(kind=dp) :: exkurt = 0.0_dp
end type moment_summary_t

type :: fit_result_t
   character(len=64) :: name = ''
   integer :: nparam = 0
   character(len=24) :: param_names(4) = ''
   real(kind=dp) :: params(4) = 0.0_dp
   type(moment_summary_t) :: moments
   real(kind=dp) :: objective = huge_obj
   logical :: success = .false.
   character(len=160) :: message = ''
   integer :: nfev = 0
end type fit_result_t

real(kind=dp), save :: g_m3_target = 0.0_dp
real(kind=dp), save :: g_m4_target = 3.0_dp

contains

function make_nan() result(x)
real(kind=dp) :: x
x = ieee_value(0.0_dp, ieee_quiet_nan)
end function make_nan

subroutine validate_target_moments(m3_target, m4_target)
real(kind=dp), intent(in) :: m3_target, m4_target
if (m4_target < m3_target*m3_target + 1.0_dp - 1.0e-12_dp) then
   error stop 'invalid target: need m4 >= m3^2 + 1'
end if
end subroutine validate_target_moments

function central_stats_from_raw(m1, m2, m3, m4, ok) result(stats)
real(kind=dp), intent(in) :: m1, m2, m3, m4
logical, intent(out) :: ok
type(moment_summary_t) :: stats
real(kind=dp) :: mean, var, c3, c4

mean = m1
var = m2 - mean*mean
if (var <= 0.0_dp) then
   ok = .false.
   stats = bad_moment_summary()
   return
end if

c3 = m3 - 3.0_dp*mean*m2 + 2.0_dp*mean**3
c4 = m4 - 4.0_dp*mean*m3 + 6.0_dp*mean*mean*m2 - 3.0_dp*mean**4
stats%m1 = m1
stats%m2 = m2
stats%m3 = m3
stats%m4 = m4
stats%mean = mean
stats%var = var
stats%skew = c3 / var**1.5_dp
stats%kurt = c4 / (var*var)
stats%exkurt = stats%kurt - 3.0_dp
ok = .true.
end function central_stats_from_raw

function summary_from_mean_var_skew_kurt(mean, var, skew, kurt) result(stats)
real(kind=dp), intent(in) :: mean, var, skew, kurt
type(moment_summary_t) :: stats
real(kind=dp) :: c3, c4
c3 = skew * var**1.5_dp
c4 = kurt * var * var
stats%mean = mean
stats%var = var
stats%skew = skew
stats%kurt = kurt
stats%exkurt = kurt - 3.0_dp
stats%m1 = mean
stats%m2 = var + mean*mean
stats%m3 = c3 + 3.0_dp*mean*var + mean**3
stats%m4 = c4 + 4.0_dp*mean*stats%m3 - 6.0_dp*mean*mean*stats%m2 + 3.0_dp*mean**4
end function summary_from_mean_var_skew_kurt

function bad_moment_summary() result(stats)
type(moment_summary_t) :: stats
stats%m1 = make_nan()
stats%m2 = make_nan()
stats%m3 = make_nan()
stats%m4 = make_nan()
stats%mean = make_nan()
stats%var = make_nan()
stats%skew = make_nan()
stats%kurt = make_nan()
stats%exkurt = make_nan()
end function bad_moment_summary

subroutine print_fit_result(result)
type(fit_result_t), intent(in) :: result
integer :: i
print '(a)', trim(result%name)
print '(a)', repeat('-', len_trim(result%name))
do i = 1, result%nparam
   if (result%params(i) > 0.5_dp*huge(1.0_dp)) then
      write(*,'(a16,a)') trim(result%param_names(i)), ' = inf'
   else if (result%params(i) /= result%params(i)) then
      write(*,'(a16,a)') trim(result%param_names(i)), ' = nan'
   else
      write(*,'(a16,a,es20.12)') trim(result%param_names(i)), ' = ', result%params(i)
   end if
end do
if (ieee_is_nan(result%moments%mean)) then
   print '(a)', 'fit mean         = nan'
   print '(a)', 'fit variance     = nan'
   print '(a)', 'fit skew         = nan'
   print '(a)', 'fit kurtosis     = nan'
   print '(a)', 'fit exkurt       = nan'
else
   write(*,'(a,es20.12)') 'fit mean         = ', result%moments%mean
   write(*,'(a,es20.12)') 'fit variance     = ', result%moments%var
   write(*,'(a,es20.12)') 'fit skew         = ', result%moments%skew
   write(*,'(a,es20.12)') 'fit kurtosis     = ', result%moments%kurt
   write(*,'(a,es20.12)') 'fit exkurt       = ', result%moments%exkurt
end if
write(*,'(a,es14.6)') 'objective        = ', result%objective
write(*,'(a,l1)') 'success          = ', result%success
if (result%nfev > 0) then
   write(*,'(a,i0)') 'nfev             = ', result%nfev
end if
write(*,'(a,a)') 'message          = ', trim(result%message)
print *
end subroutine print_fit_result

pure real(kind=dp) function log_beta_fn(x, y) result(val)
real(kind=dp), intent(in) :: x, y
val = log_gamma(x) + log_gamma(y) - log_gamma(x + y)
end function log_beta_fn

pure real(kind=dp) function beta_fn(x, y) result(val)
real(kind=dp), intent(in) :: x, y
val = exp(log_beta_fn(x, y))
end function beta_fn

pure real(kind=dp) function binom_int(n, k) result(val)
integer, intent(in) :: n, k
integer :: j, kk
if (k < 0 .or. k > n) then
   val = 0.0_dp
   return
end if
if (k == 0 .or. k == n) then
   val = 1.0_dp
   return
end if
kk = min(k, n-k)
val = 1.0_dp
do j = 1, kk
   val = val * real(n-kk+j, dp) / real(j, dp)
end do
end function binom_int

function student_t_fit_from_kurtosis(m4_target) result(res)
real(kind=dp), intent(in) :: m4_target
type(fit_result_t) :: res
real(kind=dp) :: exkurt, nu
res%name = 'symmetric_student_t'
res%nparam = 1
res%param_names(1) = 'nu'
exkurt = m4_target - 3.0_dp
if (exkurt < -tiny) then
   res%params(1) = make_nan()
   res%moments = bad_moment_summary()
   res%objective = huge_obj
   res%success = .false.
   res%message = 'target kurtosis is below 3; no symmetric Student t fit exists'
   return
end if
if (abs(exkurt) <= tiny) then
   res%params(1) = huge(1.0_dp)
   res%moments = summary_from_mean_var_skew_kurt(0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp)
   res%objective = 0.0_dp
   res%success = .true.
   res%message = 'normal limit'
   return
end if
nu = 4.0_dp + 6.0_dp / exkurt
res%params(1) = nu
res%moments = summary_from_mean_var_skew_kurt(0.0_dp, 1.0_dp, 0.0_dp, m4_target)
res%objective = 0.0_dp
res%success = .true.
res%message = 'matched kurtosis exactly'
end function student_t_fit_from_kurtosis

pure real(kind=dp) function t_log_abs_moment(r, nu) result(val)
integer, intent(in) :: r
real(kind=dp), intent(in) :: nu
val = 0.5_dp*real(r, dp)*log(nu - 2.0_dp) + log_gamma(0.5_dp*(real(r, dp) + 1.0_dp)) + &
      log_gamma(0.5_dp*(nu - real(r, dp))) - 0.5_dp*log(acos(-1.0_dp)) - log_gamma(0.5_dp*nu)
end function t_log_abs_moment

pure real(kind=dp) function t_abs_moment(r, nu) result(val)
integer, intent(in) :: r
real(kind=dp), intent(in) :: nu
val = exp(t_log_abs_moment(r, nu))
end function t_abs_moment

pure real(kind=dp) function jf_raw_moment(n, a, b) result(mu)
integer, intent(in) :: n
real(kind=dp), intent(in) :: a, b
integer :: i
real(kind=dp) :: hn, sum_terms, sign_i, log_beta_ab
hn = 0.5_dp*real(n, dp)
log_beta_ab = log_beta_fn(a, b)
sum_terms = 0.0_dp
do i = 0, n
   if (mod(i,2) == 0) then
      sign_i = 1.0_dp
   else
      sign_i = -1.0_dp
   end if
   sum_terms = sum_terms + sign_i * binom_int(n, i) * exp(log_beta_fn(a + hn - real(i, dp), b - hn + real(i, dp)) - log_beta_ab)
end do
mu = (a + b)**hn * sum_terms / (2.0_dp**n)
end function jf_raw_moment

function jf_stats(a, b, ok) result(stats)
real(kind=dp), intent(in) :: a, b
logical, intent(out) :: ok
type(moment_summary_t) :: stats
stats = central_stats_from_raw(jf_raw_moment(1, a, b), jf_raw_moment(2, a, b), jf_raw_moment(3, a, b), jf_raw_moment(4, a, b), ok)
end function jf_stats

subroutine fs_raw_from_abs_moments(m_abs, xi, raw)
real(kind=dp), intent(in) :: m_abs(4), xi
real(kind=dp), intent(out) :: raw(4)
real(kind=dp) :: denom
integer :: n
denom = xi + 1.0_dp/xi
do n = 1, 4
   raw(n) = m_abs(n) * (xi**(n+1) + (-1.0_dp)**n * xi**(-(n+1))) / denom
end do
end subroutine fs_raw_from_abs_moments

pure real(kind=dp) function ged_log_abs_moment(r, p) result(val)
integer, intent(in) :: r
real(kind=dp), intent(in) :: p
val = log_gamma((real(r, dp) + 1.0_dp)/p) - log_gamma(1.0_dp/p) + 0.5_dp*real(r, dp)*(log_gamma(1.0_dp/p) - log_gamma(3.0_dp/p))
end function ged_log_abs_moment

pure real(kind=dp) function ged_abs_moment(r, p) result(val)
integer, intent(in) :: r
real(kind=dp), intent(in) :: p
val = exp(ged_log_abs_moment(r, p))
end function ged_abs_moment

pure real(kind=dp) function ged_kurtosis(p) result(val)
real(kind=dp), intent(in) :: p
val = exp(log_gamma(5.0_dp/p) + log_gamma(1.0_dp/p) - 2.0_dp*log_gamma(3.0_dp/p))
end function ged_kurtosis

function fit_ged_from_kurtosis(m4_target) result(res)
real(kind=dp), intent(in) :: m4_target
type(fit_result_t) :: res
real(kind=dp) :: p, kurt
logical :: conv
integer :: nfev
res%name = 'symmetric_ged'
res%nparam = 1
res%param_names(1) = 'p'
if (m4_target < ged_kurt_min - 1.0e-12_dp) then
   res%params(1) = make_nan()
   res%moments = bad_moment_summary()
   res%objective = huge_obj
   res%success = .false.
   res%message = 'target kurtosis is below 1.8; no symmetric GED fit exists'
   return
end if
if (abs(m4_target - 3.0_dp) <= 1.0e-12_dp) then
   res%params(1) = 2.0_dp
   res%moments = summary_from_mean_var_skew_kurt(0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp)
   res%objective = 0.0_dp
   res%success = .true.
   res%message = 'normal case'
   return
end if
if (abs(m4_target - ged_kurt_min) <= 1.0e-10_dp) then
   res%params(1) = huge(1.0_dp)
   res%moments = summary_from_mean_var_skew_kurt(0.0_dp, 1.0_dp, 0.0_dp, ged_kurt_min)
   res%objective = 0.0_dp
   res%success = .true.
   res%message = 'uniform limit p = infinity'
   return
end if
call solve_ged_shape(m4_target, p, conv, nfev)
kurt = ged_kurtosis(p)
res%params(1) = p
res%moments = summary_from_mean_var_skew_kurt(0.0_dp, 1.0_dp, 0.0_dp, kurt)
res%objective = abs(kurt - m4_target)
res%success = conv
res%nfev = nfev
res%message = 'matched kurtosis by bisection'
end function fit_ged_from_kurtosis

subroutine solve_ged_shape(m4_target, p, converged, nfev)
real(kind=dp), intent(in) :: m4_target
real(kind=dp), intent(out) :: p
logical, intent(out) :: converged
integer, intent(out) :: nfev
real(kind=dp) :: lo, hi, mid, flo, fhi, fm
integer :: iter
lo = exp(-8.0_dp)
hi = exp(8.0_dp)
flo = ged_kurtosis(lo) - m4_target
fhi = ged_kurtosis(hi) - m4_target
nfev = 2
if (flo*fhi > 0.0_dp) then
   p = 2.0_dp
   converged = .false.
   return
end if
do iter = 1, 200
   mid = sqrt(lo*hi)
   fm = ged_kurtosis(mid) - m4_target
   nfev = nfev + 1
   if (abs(fm) <= 1.0e-13_dp) then
      p = mid
      converged = .true.
      return
   end if
   if (flo*fm <= 0.0_dp) then
      hi = mid
      fhi = fm
   else
      lo = mid
      flo = fm
   end if
   if (abs(log(hi/lo)) <= 1.0e-12_dp) then
      p = sqrt(lo*hi)
      converged = .true.
      return
   end if
end do
p = sqrt(lo*hi)
converged = .false.
end subroutine solve_ged_shape

function skew_ged_stats(xi, p, ok) result(stats)
real(kind=dp), intent(in) :: xi, p
logical, intent(out) :: ok
type(moment_summary_t) :: stats
real(kind=dp) :: m_abs(4), raw(4)
integer :: r
do r = 1, 4
   m_abs(r) = ged_abs_moment(r, p)
end do
call fs_raw_from_abs_moments(m_abs, xi, raw)
stats = central_stats_from_raw(raw(1), raw(2), raw(3), raw(4), ok)
end function skew_ged_stats

function skew_t_fs_stats(xi, nu, ok) result(stats)
real(kind=dp), intent(in) :: xi, nu
logical, intent(out) :: ok
type(moment_summary_t) :: stats
real(kind=dp) :: m_abs(4), raw(4)
integer :: r
do r = 1, 4
   m_abs(r) = t_abs_moment(r, nu)
end do
call fs_raw_from_abs_moments(m_abs, xi, raw)
stats = central_stats_from_raw(raw(1), raw(2), raw(3), raw(4), ok)
end function skew_t_fs_stats

subroutine azzalini_skew_t_raw_moments(alpha, nu, raw)
real(kind=dp), intent(in) :: alpha, nu
real(kind=dp), intent(out) :: raw(4)
real(kind=dp) :: delta, bconst

delta = alpha / sqrt(1.0_dp + alpha*alpha)
bconst = sqrt(2.0_dp / acos(-1.0_dp))
raw(1) = bconst * delta * mix_coeff(1, nu)
raw(2) = mix_coeff(2, nu)
raw(3) = bconst * delta * (3.0_dp - delta*delta) * mix_coeff(3, nu)
raw(4) = 3.0_dp * mix_coeff(4, nu)
contains
   pure real(kind=dp) function mix_coeff(r, nu) result(val)
   integer, intent(in) :: r
   real(kind=dp), intent(in) :: nu
   val = exp(0.5_dp*real(r, dp)*log(nu/2.0_dp) + log_gamma(0.5_dp*(nu - real(r, dp))) - log_gamma(0.5_dp*nu))
   end function mix_coeff
end subroutine azzalini_skew_t_raw_moments

function azzalini_skew_t_stats(alpha, nu, ok) result(stats)
real(kind=dp), intent(in) :: alpha, nu
logical, intent(out) :: ok
type(moment_summary_t) :: stats
real(kind=dp) :: raw(4)
call azzalini_skew_t_raw_moments(alpha, nu, raw)
stats = central_stats_from_raw(raw(1), raw(2), raw(3), raw(4), ok)
end function azzalini_skew_t_stats

function nig_stats(alpha, beta, delta, mu, ok) result(stats)
real(kind=dp), intent(in) :: alpha, beta, delta, mu
logical, intent(out) :: ok
type(moment_summary_t) :: stats
real(kind=dp) :: gamma, mean, var, skew, exkurt, kurt
if (alpha <= 0.0_dp .or. delta <= 0.0_dp .or. abs(beta) >= alpha) then
   ok = .false.
   stats = bad_moment_summary()
   return
end if
gamma = sqrt(alpha*alpha - beta*beta)
mean = mu + delta*beta/gamma
var = delta*alpha*alpha / gamma**3
skew = 3.0_dp*beta / (alpha*sqrt(delta*gamma))
exkurt = 3.0_dp*(1.0_dp + 4.0_dp*(beta/alpha)**2) / (delta*gamma)
kurt = 3.0_dp + exkurt
stats = summary_from_mean_var_skew_kurt(mean, var, skew, kurt)
ok = .true.
end function nig_stats

function vg_stats(r, theta, sigma, mu, ok) result(stats)
real(kind=dp), intent(in) :: r, theta, sigma, mu
logical, intent(out) :: ok
type(moment_summary_t) :: stats
real(kind=dp) :: mean, var, c3, c4, skew, kurt
if (r <= 0.0_dp .or. sigma <= 0.0_dp) then
   ok = .false.
   stats = bad_moment_summary()
   return
end if
mean = mu + r*theta
var = r*(sigma*sigma + 2.0_dp*theta*theta)
c3 = 2.0_dp*r*theta*(3.0_dp*sigma*sigma + 4.0_dp*theta*theta)
c4 = 3.0_dp*r*((r + 2.0_dp)*sigma**4 + (4.0_dp*r + 16.0_dp)*theta*theta*sigma*sigma + (4.0_dp*r + 16.0_dp)*theta**4)
skew = c3 / var**1.5_dp
kurt = c4 / (var*var)
stats = summary_from_mean_var_skew_kurt(mean, var, skew, kurt)
ok = .true.
end function vg_stats

function fit_symmetric_nig_from_kurtosis(m4_target) result(res)
real(kind=dp), intent(in) :: m4_target
type(fit_result_t) :: res
real(kind=dp) :: exkurt, alpha
logical :: ok
res%name = 'symmetric_normal_inverse_gaussian'
res%nparam = 4
res%param_names(1:4) = [character(len=24) :: 'alpha','beta','delta','mu']
exkurt = m4_target - 3.0_dp
if (exkurt < -tiny) then
   res%params(1) = make_nan()
   res%params(2) = 0.0_dp
   res%params(3) = make_nan()
   res%params(4) = 0.0_dp
   res%moments = bad_moment_summary()
   res%objective = huge_obj
   res%success = .false.
   res%message = 'target kurtosis is below 3; no symmetric NIG fit exists'
   return
end if
if (abs(exkurt) <= tiny) then
   res%params = [huge(1.0_dp), 0.0_dp, 1.0_dp, 0.0_dp]
   res%moments = summary_from_mean_var_skew_kurt(0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp)
   res%objective = 0.0_dp
   res%success = .true.
   res%message = 'normal limit'
   return
end if
alpha = 3.0_dp / exkurt
res%params = [alpha, 0.0_dp, 1.0_dp, 0.0_dp]
res%moments = nig_stats(alpha, 0.0_dp, 1.0_dp, 0.0_dp, ok)
res%objective = 0.0_dp
res%success = .true.
res%message = 'matched kurtosis exactly'
end function fit_symmetric_nig_from_kurtosis

function fit_nig(m3_target, m4_target) result(res)
real(kind=dp), intent(in) :: m3_target, m4_target
type(fit_result_t) :: res
real(kind=dp) :: exkurt, s, denom, q, rho, eta, alpha, beta
logical :: ok
res%name = 'normal_inverse_gaussian'
res%nparam = 4
res%param_names(1:4) = [character(len=24) :: 'alpha','beta','delta','mu']
exkurt = m4_target - 3.0_dp
s = m3_target
if (exkurt < -tiny) then
   call set_failed(res, 'target kurtosis is below 3; no NIG fit exists')
   return
end if
if (abs(exkurt) <= tiny) then
   if (abs(s) <= tiny) then
      res%params = [huge(1.0_dp), 0.0_dp, 1.0_dp, 0.0_dp]
      res%moments = summary_from_mean_var_skew_kurt(0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp)
      res%objective = 0.0_dp
      res%success = .true.
      res%message = 'normal limit'
   else
      call set_failed(res, 'nonzero skew with zero excess kurtosis cannot be matched by NIG')
   end if
   return
end if

denom = 3.0_dp*exkurt - 4.0_dp*s*s
if (denom <= tiny) then
   call set_failed(res, 'target moments violate NIG region: need excess kurtosis > 4*skew^2/3')
   return
end if
q = s*s / denom
if (q >= 1.0_dp) then
   call set_failed(res, 'target moments fall outside NIG region')
   return
end if
if (abs(s) > tiny) then
   rho = sign(sqrt(q), s)
else
   rho = 0.0_dp
end if
eta = 3.0_dp*(1.0_dp + 4.0_dp*q) / (exkurt*sqrt(1.0_dp - q))
alpha = eta
beta = rho*eta
res%params = [alpha, beta, 1.0_dp, 0.0_dp]
res%moments = nig_stats(alpha, beta, 1.0_dp, 0.0_dp, ok)
res%objective = (res%moments%skew - s)**2 + (res%moments%kurt - m4_target)**2
res%success = .true.
res%message = 'matched skewness and kurtosis by closed form'
contains
   subroutine set_failed(res, msg)
   type(fit_result_t), intent(inout) :: res
   character(len=*), intent(in) :: msg
   res%params = [make_nan(), make_nan(), make_nan(), 0.0_dp]
   res%moments = bad_moment_summary()
   res%objective = huge_obj
   res%success = .false.
   res%message = msg
   end subroutine set_failed
end function fit_nig

function fit_symmetric_vg_from_kurtosis(m4_target) result(res)
real(kind=dp), intent(in) :: m4_target
type(fit_result_t) :: res
real(kind=dp) :: exkurt, r
logical :: ok
res%name = 'symmetric_variance_gamma'
res%nparam = 4
res%param_names(1:4) = [character(len=24) :: 'r','theta','sigma','mu']
exkurt = m4_target - 3.0_dp
if (exkurt < -tiny) then
   res%params = [make_nan(), 0.0_dp, 1.0_dp, 0.0_dp]
   res%moments = bad_moment_summary()
   res%objective = huge_obj
   res%success = .false.
   res%message = 'target kurtosis is below 3; no symmetric VG fit exists'
   return
end if
if (abs(exkurt) <= tiny) then
   res%params = [huge(1.0_dp), 0.0_dp, 1.0_dp, 0.0_dp]
   res%moments = summary_from_mean_var_skew_kurt(0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp)
   res%objective = 0.0_dp
   res%success = .true.
   res%message = 'normal limit'
   return
end if
r = 6.0_dp / exkurt
res%params = [r, 0.0_dp, 1.0_dp, 0.0_dp]
res%moments = vg_stats(r, 0.0_dp, 1.0_dp, 0.0_dp, ok)
res%objective = 0.0_dp
res%success = .true.
res%message = 'matched kurtosis exactly'
end function fit_symmetric_vg_from_kurtosis

function fit_vg(m3_target, m4_target) result(res)
real(kind=dp), intent(in) :: m3_target, m4_target
type(fit_result_t) :: res
real(kind=dp) :: exkurt, s, ratio_target, y_hi, y, theta, r
logical :: ok, conv
integer :: nfev
res%name = 'variance_gamma'
res%nparam = 4
res%param_names(1:4) = [character(len=24) :: 'r','theta','sigma','mu']
exkurt = m4_target - 3.0_dp
s = m3_target
if (exkurt < -tiny) then
   call set_failed('target kurtosis is below 3; no VG fit exists')
   return
end if
if (abs(exkurt) <= tiny) then
   if (abs(s) <= tiny) then
      res%params = [huge(1.0_dp), 0.0_dp, 1.0_dp, 0.0_dp]
      res%moments = summary_from_mean_var_skew_kurt(0.0_dp, 1.0_dp, 0.0_dp, 3.0_dp)
      res%objective = 0.0_dp
      res%success = .true.
      res%message = 'normal limit'
   else
      call set_failed('nonzero skew with zero excess kurtosis cannot be matched by VG')
   end if
   return
end if
ratio_target = (s*s) / exkurt
if (ratio_target >= 2.0_dp/3.0_dp - 1.0e-12_dp) then
   call set_failed('target moments violate VG region: need skew^2 / excess kurtosis < 2/3')
   return
end if
if (abs(s) <= tiny) then
   r = 6.0_dp / exkurt
   res%params = [r, 0.0_dp, 1.0_dp, 0.0_dp]
   res%moments = vg_stats(r, 0.0_dp, 1.0_dp, 0.0_dp, ok)
   res%objective = (res%moments%kurt - m4_target)**2
   res%success = .true.
   res%message = 'symmetric solution'
   return
end if
y_hi = 1.0_dp
do while (vg_ratio_y(y_hi) - ratio_target < 0.0_dp .and. y_hi < 1.0e8_dp)
   y_hi = 2.0_dp*y_hi
end do
if (vg_ratio_y(y_hi) - ratio_target < 0.0_dp) then
   call set_failed('failed to bracket VG shape root')
   return
end if
call solve_vg_y(ratio_target, 0.0_dp, y_hi, y, conv, nfev)
theta = sign(sqrt(y), s)
r = 6.0_dp*(1.0_dp + 8.0_dp*y + 8.0_dp*y*y) / (exkurt*(1.0_dp + 2.0_dp*y)**2)
res%params = [r, theta, 1.0_dp, 0.0_dp]
res%moments = vg_stats(r, theta, 1.0_dp, 0.0_dp, ok)
res%objective = (res%moments%skew - s)**2 + (res%moments%kurt - m4_target)**2
res%success = conv
res%nfev = nfev
res%message = 'matched skewness and kurtosis by 1D root solve'
contains
   subroutine set_failed(msg)
   character(len=*), intent(in) :: msg
   res%params = [make_nan(), make_nan(), 1.0_dp, 0.0_dp]
   res%moments = bad_moment_summary()
   res%objective = huge_obj
   res%success = .false.
   res%message = msg
   end subroutine set_failed
end function fit_vg

pure real(kind=dp) function vg_ratio_y(y) result(val)
real(kind=dp), intent(in) :: y
val = 4.0_dp*y*(3.0_dp + 4.0_dp*y)**2 / (6.0_dp*(1.0_dp + 8.0_dp*y + 8.0_dp*y*y)*(1.0_dp + 2.0_dp*y))
end function vg_ratio_y

subroutine solve_vg_y(target, lo0, hi0, y, converged, nfev)
real(kind=dp), intent(in) :: target, lo0, hi0
real(kind=dp), intent(out) :: y
logical, intent(out) :: converged
integer, intent(out) :: nfev
real(kind=dp) :: lo, hi, mid, flo, fhi, fm
integer :: iter
lo = lo0
hi = hi0
flo = vg_ratio_y(lo) - target
fhi = vg_ratio_y(hi) - target
nfev = 2
do iter = 1, 200
   mid = 0.5_dp*(lo + hi)
   fm = vg_ratio_y(mid) - target
   nfev = nfev + 1
   if (abs(fm) <= 1.0e-13_dp .or. abs(hi - lo) <= 1.0e-13_dp*max(1.0_dp, abs(mid))) then
      y = mid
      converged = .true.
      return
   end if
   if (flo*fm <= 0.0_dp) then
      hi = mid
      fhi = fm
   else
      lo = mid
      flo = fm
   end if
end do
y = 0.5_dp*(lo + hi)
converged = .false.
end subroutine solve_vg_y

function fit_jf(m3_target, m4_target) result(res)
real(kind=dp), intent(in) :: m3_target, m4_target
type(fit_result_t) :: res
real(kind=dp) :: starts(2,8), lower(2), upper(2), xbest(2), best_obj, a, b
logical :: success, ok
integer :: nfev
character(len=160) :: message
res%name = 'jones_faddy_skew_t'
res%nparam = 2
res%param_names(1:2) = [character(len=24) :: 'a','b']
g_m3_target = m3_target
g_m4_target = m4_target
starts(:,1) = [0.0_dp, 0.0_dp]
starts(:,2) = [log(2.0_dp), 0.0_dp]
starts(:,3) = [0.0_dp, log(2.0_dp)]
starts(:,4) = [log(4.0_dp), 0.0_dp]
starts(:,5) = [0.0_dp, log(4.0_dp)]
starts(:,6) = [log(8.0_dp), log(2.0_dp)]
starts(:,7) = [log(2.0_dp), log(8.0_dp)]
starts(:,8) = [log(20.0_dp), log(20.0_dp)]
lower = [-3.0_dp, -3.0_dp]
upper = [8.0_dp, 8.0_dp]
call solve_least_squares_2d(jf_residuals, starts, lower, upper, xbest, best_obj, success, nfev, message)
a = 2.0_dp + exp(xbest(1))
b = 2.0_dp + exp(xbest(2))
res%params(1:2) = [a, b]
res%moments = jf_stats(a, b, ok)
res%objective = best_obj
res%success = success
res%nfev = nfev
res%message = message
end function fit_jf

subroutine jf_residuals(x, r, ok)
real(kind=dp), intent(in) :: x(2)
real(kind=dp), intent(out) :: r(2)
logical, intent(out) :: ok
type(moment_summary_t) :: stats
real(kind=dp) :: a, b
a = 2.0_dp + exp(x(1))
b = 2.0_dp + exp(x(2))
stats = jf_stats(a, b, ok)
if (ok) then
   r = [stats%skew - g_m3_target, stats%kurt - g_m4_target]
else
   r = 0.0_dp
end if
end subroutine jf_residuals

function fit_skew_ged(m3_target, m4_target) result(res)
real(kind=dp), intent(in) :: m3_target, m4_target
type(fit_result_t) :: res
real(kind=dp) :: starts(2,7), lower(2), upper(2), xbest(2), best_obj, xi, p
logical :: success, ok
integer :: nfev
character(len=160) :: message
res%name = 'skew_ged_fernandez_steel'
res%nparam = 2
res%param_names(1:2) = [character(len=24) :: 'xi','p']
g_m3_target = m3_target
g_m4_target = m4_target
starts(:,1) = [0.0_dp, log(2.0_dp)]
starts(:,2) = [log(1.5_dp), log(2.0_dp)]
starts(:,3) = [-log(1.5_dp), log(2.0_dp)]
starts(:,4) = [log(2.0_dp), log(1.0_dp)]
starts(:,5) = [-log(2.0_dp), log(1.0_dp)]
starts(:,6) = [log(3.0_dp), log(0.7_dp)]
starts(:,7) = [-log(3.0_dp), log(0.7_dp)]
lower = [-3.5_dp, -3.5_dp]
upper = [3.5_dp, 3.5_dp]
call solve_least_squares_2d(skew_ged_residuals, starts, lower, upper, xbest, best_obj, success, nfev, message)
xi = exp(xbest(1))
p = exp(xbest(2))
res%params(1:2) = [xi, p]
res%moments = skew_ged_stats(xi, p, ok)
res%objective = best_obj
res%success = success
res%nfev = nfev
res%message = message
end function fit_skew_ged

subroutine skew_ged_residuals(x, r, ok)
real(kind=dp), intent(in) :: x(2)
real(kind=dp), intent(out) :: r(2)
logical, intent(out) :: ok
type(moment_summary_t) :: stats
stats = skew_ged_stats(exp(x(1)), exp(x(2)), ok)
if (ok) then
   r = [stats%skew - g_m3_target, stats%kurt - g_m4_target]
else
   r = 0.0_dp
end if
end subroutine skew_ged_residuals

function fit_skew_t_fs(m3_target, m4_target) result(res)
real(kind=dp), intent(in) :: m3_target, m4_target
type(fit_result_t) :: res
real(kind=dp) :: starts(2,7), lower(2), upper(2), xbest(2), best_obj, xi, nu
logical :: success, ok
integer :: nfev
character(len=160) :: message
res%name = 'skew_student_t_fernandez_steel'
res%nparam = 2
res%param_names(1:2) = [character(len=24) :: 'xi','nu']
g_m3_target = m3_target
g_m4_target = m4_target
starts(:,1) = [0.0_dp, 0.0_dp]
starts(:,2) = [log(1.5_dp), log(2.0_dp)]
starts(:,3) = [-log(1.5_dp), log(2.0_dp)]
starts(:,4) = [log(2.0_dp), log(4.0_dp)]
starts(:,5) = [-log(2.0_dp), log(4.0_dp)]
starts(:,6) = [log(3.0_dp), log(8.0_dp)]
starts(:,7) = [-log(3.0_dp), log(8.0_dp)]
lower = [-3.5_dp, -3.0_dp]
upper = [3.5_dp, 6.0_dp]
call solve_least_squares_2d(skew_t_fs_residuals, starts, lower, upper, xbest, best_obj, success, nfev, message)
xi = exp(xbest(1))
nu = 4.0_dp + exp(xbest(2))
res%params(1:2) = [xi, nu]
res%moments = skew_t_fs_stats(xi, nu, ok)
res%objective = best_obj
res%success = success
res%nfev = nfev
res%message = message
end function fit_skew_t_fs

subroutine skew_t_fs_residuals(x, r, ok)
real(kind=dp), intent(in) :: x(2)
real(kind=dp), intent(out) :: r(2)
logical, intent(out) :: ok
type(moment_summary_t) :: stats
stats = skew_t_fs_stats(exp(x(1)), 4.0_dp + exp(x(2)), ok)
if (ok) then
   r = [stats%skew - g_m3_target, stats%kurt - g_m4_target]
else
   r = 0.0_dp
end if
end subroutine skew_t_fs_residuals

function fit_azzalini_skew_t(m3_target, m4_target) result(res)
real(kind=dp), intent(in) :: m3_target, m4_target
type(fit_result_t) :: res
real(kind=dp) :: starts(2,7), lower(2), upper(2), xbest(2), best_obj, alpha, nu
logical :: success, ok
integer :: nfev
character(len=160) :: message
res%name = 'azzalini_skew_student_t'
res%nparam = 2
res%param_names(1:2) = [character(len=24) :: 'alpha','nu']
g_m3_target = m3_target
g_m4_target = m4_target
starts(:,1) = [0.0_dp, 0.0_dp]
starts(:,2) = [1.0_dp, log(2.0_dp)]
starts(:,3) = [-1.0_dp, log(2.0_dp)]
starts(:,4) = [3.0_dp, log(4.0_dp)]
starts(:,5) = [-3.0_dp, log(4.0_dp)]
starts(:,6) = [8.0_dp, log(8.0_dp)]
starts(:,7) = [-8.0_dp, log(8.0_dp)]
lower = [-30.0_dp, -3.0_dp]
upper = [30.0_dp, 6.0_dp]
call solve_least_squares_2d(azzalini_residuals, starts, lower, upper, xbest, best_obj, success, nfev, message)
alpha = xbest(1)
nu = 4.0_dp + exp(xbest(2))
res%params(1:2) = [alpha, nu]
res%moments = azzalini_skew_t_stats(alpha, nu, ok)
res%objective = best_obj
res%success = success
res%nfev = nfev
res%message = message
end function fit_azzalini_skew_t

subroutine azzalini_residuals(x, r, ok)
real(kind=dp), intent(in) :: x(2)
real(kind=dp), intent(out) :: r(2)
logical, intent(out) :: ok
type(moment_summary_t) :: stats
stats = azzalini_skew_t_stats(x(1), 4.0_dp + exp(x(2)), ok)
if (ok) then
   r = [stats%skew - g_m3_target, stats%kurt - g_m4_target]
else
   r = 0.0_dp
end if
end subroutine azzalini_residuals

subroutine solve_least_squares_2d(func, starts, lower, upper, xbest, best_obj, success, nfev_total, message)
procedure(residual_2d_interface) :: func
real(kind=dp), intent(in) :: starts(:,:), lower(2), upper(2)
real(kind=dp), intent(out) :: xbest(2), best_obj
logical, intent(out) :: success
integer, intent(out) :: nfev_total
character(len=*), intent(out) :: message
integer :: nstart, i, nfev_local
real(kind=dp) :: x(2), obj
logical :: succ_local
character(len=160) :: msg_local
nstart = size(starts, 2)
best_obj = huge_obj
xbest = starts(:,1)
success = .false.
message = 'no start attempted'
nfev_total = 0
do i = 1, nstart
   x = starts(:,i)
   call project_bounds(x, lower, upper)
   call solve_one_start(func, x, lower, upper, obj, succ_local, nfev_local, msg_local)
   nfev_total = nfev_total + nfev_local
   if (obj < best_obj) then
      best_obj = obj
      xbest = x
      success = succ_local
      message = msg_local
   end if
end do
if (best_obj > 1.0e-8_dp) then
   success = .false.
   if (len_trim(message) == 0) message = 'best residual remains nonzero'
else
   success = .true.
end if
end subroutine solve_least_squares_2d

subroutine solve_one_start(func, x, lower, upper, obj, success, nfev, message)
procedure(residual_2d_interface) :: func
real(kind=dp), intent(inout) :: x(2)
real(kind=dp), intent(in) :: lower(2), upper(2)
real(kind=dp), intent(out) :: obj
logical, intent(out) :: success
integer, intent(out) :: nfev
character(len=*), intent(out) :: message
real(kind=dp) :: r(2), rnew(2), j(2,2), step(2), xnew(2), jtjr(2), a11, a12, a22, det
real(kind=dp) :: lambda, obj_new, h, xt(2), rt(2), step_norm
logical :: ok, ok2
integer :: iter, k, inner

call func(x, r, ok)
nfev = 1
if (.not. ok) then
   obj = huge_obj
   success = .false.
   message = 'initial evaluation failed'
   return
end if
obj = dot_product(r, r)
lambda = 1.0e-3_dp

do iter = 1, 120
   do k = 1, 2
      xt = x
      h = 1.0e-6_dp*max(1.0_dp, abs(x(k)))
      xt(k) = min(upper(k), max(lower(k), xt(k) + h))
      if (xt(k) == x(k)) xt(k) = min(upper(k), x(k) - h)
      call func(xt, rt, ok2)
      nfev = nfev + 1
      if (.not. ok2) then
         j(:,k) = 0.0_dp
      else
         j(:,k) = (rt - r) / (xt(k) - x(k))
      end if
   end do
   jtjr = matmul(transpose(j), r)
   a11 = dot_product(j(:,1), j(:,1)) + lambda
   a12 = dot_product(j(:,1), j(:,2))
   a22 = dot_product(j(:,2), j(:,2)) + lambda
   det = a11*a22 - a12*a12
   if (abs(det) <= 1.0e-30_dp) then
      step = -0.1_dp*jtjr
   else
      step(1) = (-a22*jtjr(1) + a12*jtjr(2)) / det
      step(2) = ( a12*jtjr(1) - a11*jtjr(2)) / det
   end if
   step_norm = max(abs(step(1)), abs(step(2)))
   if (step_norm <= 1.0e-12_dp) then
      success = (obj <= 1.0e-10_dp)
      message = 'step tolerance reached'
      return
   end if
   ok2 = .false.
   do inner = 1, 20
      xnew = x + step
      call project_bounds(xnew, lower, upper)
      call func(xnew, rnew, ok)
      nfev = nfev + 1
      if (ok) then
         obj_new = dot_product(rnew, rnew)
      else
         obj_new = huge_obj
      end if
      if (obj_new < obj) then
         x = xnew
         r = rnew
         obj = obj_new
         lambda = max(lambda/5.0_dp, 1.0e-12_dp)
         ok2 = .true.
         exit
      else
         lambda = min(lambda*10.0_dp, 1.0e12_dp)
         step = 0.5_dp*step
      end if
   end do
   if (.not. ok2) then
      success = (obj <= 1.0e-10_dp)
      message = 'no improving step found'
      return
   end if
   if (obj <= 1.0e-20_dp) then
      success = .true.
      message = 'residual tolerance reached'
      return
   end if
end do
success = (obj <= 1.0e-10_dp)
message = 'maximum iterations reached'
end subroutine solve_one_start

subroutine project_bounds(x, lower, upper)
real(kind=dp), intent(inout) :: x(2)
real(kind=dp), intent(in) :: lower(2), upper(2)
x = min(upper, max(lower, x))
end subroutine project_bounds

end module fit_moments_families_mod

module fit_moments_data_mod
use, intrinsic :: iso_fortran_env, only: dp => real64
use fit_moments_families_mod
implicit none
private
public :: read_data_file, sample_raw_moments, standardized_moments, print_data_summary, print_original_fit

contains

subroutine read_data_file(filename, x)
character(len=*), intent(in) :: filename
real(kind=dp), allocatable, intent(out) :: x(:)
integer :: n, i, ios, unit
character (len=1000) :: text

open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
if (ios /= 0) error stop "could not open data file " // trim(filename)
do
   read (unit, "(a)", iostat=ios) text
   if (text(1:1) == "#" .or. text(1:1) == "!") then
      print "(a)", trim(text)
   else
      read(text, *, iostat=ios) n
      if (ios /= 0 .or. n <= 0) error stop "could not read positive number of observations from first line"
      exit
   end if
end do
allocate(x(n))
do i = 1, n
   read(unit, *, iostat=ios) x(i)
   if (ios /= 0) error stop "could not read observation from data file"
end do
close(unit)
end subroutine read_data_file

function sample_raw_moments(x) result(stats)
real(kind=dp), intent(in) :: x(:)
type(moment_summary_t) :: stats
logical :: ok
real(kind=dp) :: ninv, m1, m2, m3, m4

ninv = 1.0_dp / real(size(x), dp)
m1 = ninv * sum(x)
m2 = ninv * sum(x**2)
m3 = ninv * sum(x**3)
m4 = ninv * sum(x**4)
stats = central_stats_from_raw(m1, m2, m3, m4, ok)
if (.not. ok) error stop 'sample variance is not positive'
end function sample_raw_moments

function standardized_moments(x) result(stats)
real(kind=dp), intent(in) :: x(:)
type(moment_summary_t) :: stats
real(kind=dp), allocatable :: z(:)
real(kind=dp) :: mean_x, sd_x

mean_x = sum(x) / real(size(x), dp)
sd_x = sqrt(sum((x - mean_x)**2) / real(size(x), dp))
if (sd_x <= 0.0_dp) error stop 'sample standard deviation is not positive'
allocate(z(size(x)))
z = (x - mean_x) / sd_x
stats = sample_raw_moments(z)
deallocate(z)
end function standardized_moments

subroutine print_data_summary(n, raw_stats, std_stats)
integer, intent(in) :: n
type(moment_summary_t), intent(in) :: raw_stats, std_stats
write(*,'(a,i0)') 'n               = ', n
write(*,'(a,es20.12)') 'sample mean     = ', raw_stats%mean
write(*,'(a,es20.12)') 'sample variance = ', raw_stats%var
write(*,'(a,es20.12)') 'sample std dev  = ', sqrt(raw_stats%var)
write(*,'(a,es20.12)') 'sample skew     = ', raw_stats%skew
write(*,'(a,es20.12)') 'sample kurtosis = ', raw_stats%kurt
write(*,'(a,es20.12)') 'sample exkurt   = ', raw_stats%exkurt
print *
write(*,'(a,es20.12)') 'std m1          = ', std_stats%m1
write(*,'(a,es20.12)') 'std m2          = ', std_stats%m2
write(*,'(a,es20.12)') 'std m3          = ', std_stats%m3
write(*,'(a,es20.12)') 'std m4          = ', std_stats%m4
write(*,'(a,es20.12)') 'std skew        = ', std_stats%skew
write(*,'(a,es20.12)') 'std kurtosis    = ', std_stats%kurt
write(*,'(a,es20.12)') 'std exkurt      = ', std_stats%exkurt
print *
end subroutine print_data_summary

subroutine print_original_fit(res, data_mean, data_var)
type(fit_result_t), intent(in) :: res
real(kind=dp), intent(in) :: data_mean, data_var
real(kind=dp) :: data_sd, loc_shift, scale_fac, mean_u, var_u, mu_native
real(kind=dp) :: a, b, xi, nu, alpha, p, beta, delta, mu, r, theta, sigma
real(kind=dp) :: tscale

print '(a)', trim(res%name)
print '(a)', repeat('-', len_trim(res%name))
write(*,'(a,l1)') 'success               = ', res%success
write(*,'(a,es14.6)') 'objective             = ', res%objective
if (res%nfev > 0) write(*,'(a,i0)') 'nfev                  = ', res%nfev
write(*,'(a,a)') 'message               = ', trim(res%message)
if (.not. res%success) then
   print *
   return
end if

data_sd = sqrt(data_var)
mean_u = res%moments%mean
var_u = res%moments%var
if (var_u <= 0.0_dp) then
   print '(a)', 'invalid fitted variance'
   print *
   return
end if
scale_fac = data_sd / sqrt(var_u)
loc_shift = data_mean - scale_fac*mean_u

select case (trim(res%name))
case ('symmetric_student_t')
   nu = res%params(1)
   if (nu > 0.5_dp*huge(1.0_dp)) then
      tscale = data_sd
   else
      tscale = data_sd * sqrt((nu - 2.0_dp) / nu)
   end if
   write(*,'(a,es20.12)') 'nu                    = ', nu
   write(*,'(a,es20.12)') 'loc                   = ', data_mean
   write(*,'(a,es20.12)') 'scale                 = ', tscale
case ('jones_faddy_skew_t')
   a = res%params(1)
   b = res%params(2)
   write(*,'(a,es20.12)') 'a                     = ', a
   write(*,'(a,es20.12)') 'b                     = ', b
   write(*,'(a,es20.12)') 'loc                   = ', loc_shift
   write(*,'(a,es20.12)') 'scale                 = ', scale_fac
case ('skew_student_t_fernandez_steel')
   xi = res%params(1)
   nu = res%params(2)
   write(*,'(a,es20.12)') 'xi                    = ', xi
   write(*,'(a,es20.12)') 'nu                    = ', nu
   write(*,'(a,es20.12)') 'loc                   = ', loc_shift
   write(*,'(a,es20.12)') 'scale                 = ', scale_fac
case ('azzalini_skew_student_t')
   alpha = res%params(1)
   nu = res%params(2)
   write(*,'(a,es20.12)') 'alpha                 = ', alpha
   write(*,'(a,es20.12)') 'nu                    = ', nu
   write(*,'(a,es20.12)') 'loc                   = ', loc_shift
   write(*,'(a,es20.12)') 'scale                 = ', scale_fac
case ('symmetric_ged')
   p = res%params(1)
   write(*,'(a,es20.12)') 'p                     = ', p
   write(*,'(a,es20.12)') 'loc                   = ', data_mean
   write(*,'(a,es20.12)') 'scale                 = ', data_sd
case ('skew_ged_fernandez_steel')
   xi = res%params(1)
   p = res%params(2)
   write(*,'(a,es20.12)') 'xi                    = ', xi
   write(*,'(a,es20.12)') 'p                     = ', p
   write(*,'(a,es20.12)') 'loc                   = ', loc_shift
   write(*,'(a,es20.12)') 'scale                 = ', scale_fac
case ('symmetric_normal_inverse_gaussian', 'normal_inverse_gaussian')
   alpha = res%params(1)
   beta = res%params(2)
   delta = res%params(3)
   mu = res%params(4)
   write(*,'(a,es20.12)') 'alpha                 = ', alpha / scale_fac
   write(*,'(a,es20.12)') 'beta                  = ', beta / scale_fac
   write(*,'(a,es20.12)') 'delta                 = ', delta * scale_fac
   mu_native = loc_shift + scale_fac*mu
   write(*,'(a,es20.12)') 'mu                    = ', mu_native
case ('symmetric_variance_gamma', 'variance_gamma')
   r = res%params(1)
   theta = res%params(2)
   sigma = res%params(3)
   mu = res%params(4)
   write(*,'(a,es20.12)') 'r                     = ', r
   write(*,'(a,es20.12)') 'theta                 = ', theta * scale_fac
   write(*,'(a,es20.12)') 'sigma                 = ', sigma * scale_fac
   mu_native = loc_shift + scale_fac*mu
   write(*,'(a,es20.12)') 'mu                    = ', mu_native
case default
   write(*,'(a)', advance='no') ''
end select

write(*,'(a,es20.12)') 'fit mean              = ', data_mean
write(*,'(a,es20.12)') 'fit variance          = ', data_var
write(*,'(a,es20.12)') 'fit skew              = ', res%moments%skew
write(*,'(a,es20.12)') 'fit kurtosis          = ', res%moments%kurt
write(*,'(a,es20.12)') 'fit exkurt            = ', res%moments%exkurt
print *
end subroutine print_original_fit

end module fit_moments_data_mod

module fit_dist_mle_mod
use, intrinsic :: iso_fortran_env, only: dp => real64
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
use fit_moments_families_mod
implicit none
private

real(kind=dp), parameter :: pi_dp = acos(-1.0_dp)
real(kind=dp), parameter :: log2_dp = log(2.0_dp)
real(kind=dp), parameter :: euler_gamma_dp = 0.5772156649015328606_dp
real(kind=dp), parameter :: tiny_pdf = 1.0e-300_dp
real(kind=dp), parameter :: huge_nll = 1.0e300_dp
integer, parameter :: fam_student_t = 1
integer, parameter :: fam_jf = 2
integer, parameter :: fam_fs_t = 3
integer, parameter :: fam_azzalini = 4
integer, parameter :: fam_ged = 5
integer, parameter :: fam_skew_ged = 6
integer, parameter :: fam_sym_nig = 7
integer, parameter :: fam_nig = 8
integer, parameter :: fam_sym_vg = 9
integer, parameter :: fam_vg = 10

public :: mle_result_t, fit_family_mle, print_mle_result

type :: mle_result_t
   character(len=64) :: name = ''
   integer :: nparam = 0
   character(len=24) :: param_names(4) = ''
   real(kind=dp) :: start_params(4) = 0.0_dp
   real(kind=dp) :: params(4) = 0.0_dp
   real(kind=dp) :: negloglik = huge_nll
   logical :: success = .false.
   character(len=160) :: message = ''
   integer :: nfev = 0
end type mle_result_t

real(kind=dp), allocatable, save :: g_x(:)
integer, save :: g_family = 0
real(kind=dp), save :: g_data_scale = 1.0_dp

contains

pure real(kind=dp) function log_beta_local(x, y) result(val)
real(kind=dp), intent(in) :: x, y
val = log_gamma(x) + log_gamma(y) - log_gamma(x + y)
end function log_beta_local

pure real(kind=dp) function pos_from_log(x) result(val)
real(kind=dp), intent(in) :: x
val = max(exp(max(min(x, 700.0_dp), -700.0_dp)), 1.0e-12_dp)
end function pos_from_log

subroutine set_mle_data(x)
real(kind=dp), intent(in) :: x(:)
if (allocated(g_x)) deallocate(g_x)
allocate(g_x(size(x)))
g_x = x
g_data_scale = sqrt(sum((x - sum(x)/real(size(x), dp))**2) / real(size(x), dp))
if (g_data_scale <= 0.0_dp) g_data_scale = 1.0_dp
end subroutine set_mle_data

subroutine fit_family_mle(ifam, start_fit, data_mean, data_var, x, res)
integer, intent(in) :: ifam
real(kind=dp), intent(in) :: data_mean, data_var
real(kind=dp), intent(in) :: x(:)
type(fit_result_t), intent(in) :: start_fit
type(mle_result_t), intent(out) :: res
real(kind=dp) :: start_native(4), x0(4), step(4), xbest(4), fbest
integer :: n
logical :: ok, success
character(len=160) :: message
integer :: nfev

call set_mle_data(x)
res%name = trim(start_fit%name)
call native_start_from_moments(ifam, start_fit, data_mean, data_var, start_native, n, ok)
res%nparam = n
call set_param_names(ifam, res%param_names, res%nparam)
if (.not. ok) then
   res%success = .false.
   res%message = 'moment fit failed; mle skipped'
   res%negloglik = huge_nll
   return
end if
res%start_params(1:n) = start_native(1:n)
call pack_unconstrained(ifam, start_native, n, x0, ok)
if (.not. ok) then
   res%success = .false.
   res%message = 'could not transform starting values'
   res%negloglik = huge_nll
   return
end if
call default_steps(ifam, n, x0, step)
g_family = ifam
call nelder_mead(n, x0, step, xbest, fbest, success, nfev, message)
res%success = success
res%nfev = nfev
res%message = message
res%negloglik = fbest
if (success) then
   call unpack_unconstrained(ifam, xbest, n, res%params, ok)
   if (.not. ok) then
      res%success = .false.
      res%message = 'optimizer finished at invalid parameter vector'
      res%negloglik = huge_nll
   end if
else
   call unpack_unconstrained(ifam, xbest, n, res%params, ok)
end if
end subroutine fit_family_mle

subroutine set_param_names(ifam, names, n)
integer, intent(in) :: ifam, n
character(len=24), intent(out) :: names(4)
names = ''
select case (ifam)
case (fam_student_t)
   names(1:3) = [character(len=24) :: 'nu', 'loc', 'scale']
case (fam_jf)
   names(1:4) = [character(len=24) :: 'a', 'b', 'loc', 'scale']
case (fam_fs_t)
   names(1:4) = [character(len=24) :: 'xi', 'nu', 'loc', 'scale']
case (fam_azzalini)
   names(1:4) = [character(len=24) :: 'alpha', 'nu', 'loc', 'scale']
case (fam_ged)
   names(1:3) = [character(len=24) :: 'p', 'loc', 'scale']
case (fam_skew_ged)
   names(1:4) = [character(len=24) :: 'xi', 'p', 'loc', 'scale']
case (fam_sym_nig)
   names(1:3) = [character(len=24) :: 'alpha', 'delta', 'mu']
case (fam_nig)
   names(1:4) = [character(len=24) :: 'alpha', 'beta', 'delta', 'mu']
case (fam_sym_vg)
   names(1:3) = [character(len=24) :: 'r', 'sigma', 'mu']
case (fam_vg)
   names(1:4) = [character(len=24) :: 'r', 'theta', 'sigma', 'mu']
end select
end subroutine set_param_names

subroutine native_start_from_moments(ifam, start_fit, data_mean, data_var, p, n, ok)
integer, intent(in) :: ifam
real(kind=dp), intent(in) :: data_mean, data_var
type(fit_result_t), intent(in) :: start_fit
real(kind=dp), intent(out) :: p(4)
integer, intent(out) :: n
logical, intent(out) :: ok
real(kind=dp) :: data_sd, loc_shift, scale_fac, mean_u, var_u, nu, alpha, beta, delta, mu
p = 0.0_dp
ok = start_fit%success
if (.not. ok) then
   n = 0
   return
end if
if (ieee_is_nan(start_fit%moments%var)) then
   ok = .false.
   n = 0
   return
end if

data_sd = sqrt(data_var)
mean_u = start_fit%moments%mean
var_u = start_fit%moments%var
if (var_u <= 0.0_dp) then
   ok = .false.
   n = 0
   return
end if
scale_fac = data_sd / sqrt(var_u)
loc_shift = data_mean - scale_fac*mean_u

select case (ifam)
case (fam_student_t)
   nu = start_fit%params(1)
   p(1) = nu
   p(2) = data_mean
   if (nu > 0.5_dp*huge(1.0_dp)) then
      p(3) = max(data_sd, 1.0e-6_dp)
   else
      p(3) = max(data_sd * sqrt((nu - 2.0_dp) / nu), 1.0e-6_dp)
   end if
   n = 3
case (fam_jf)
   p(1:2) = start_fit%params(1:2)
   p(3) = loc_shift
   p(4) = max(scale_fac, 1.0e-6_dp)
   n = 4
case (fam_fs_t)
   p(1:2) = start_fit%params(1:2)
   p(3) = loc_shift
   p(4) = max(scale_fac, 1.0e-6_dp)
   n = 4
case (fam_azzalini)
   p(1:2) = start_fit%params(1:2)
   p(3) = loc_shift
   p(4) = max(scale_fac, 1.0e-6_dp)
   n = 4
case (fam_ged)
   p(1) = start_fit%params(1)
   p(2) = data_mean
   p(3) = max(data_sd, 1.0e-6_dp)
   n = 3
case (fam_skew_ged)
   p(1:2) = start_fit%params(1:2)
   p(3) = loc_shift
   p(4) = max(scale_fac, 1.0e-6_dp)
   n = 4
case (fam_sym_nig)
   alpha = start_fit%params(1)
   delta = start_fit%params(3)
   mu = start_fit%params(4)
   p(1) = alpha / scale_fac
   p(2) = max(delta * scale_fac, 1.0e-6_dp)
   p(3) = loc_shift + scale_fac * mu
   n = 3
case (fam_nig)
   alpha = start_fit%params(1)
   beta = start_fit%params(2)
   delta = start_fit%params(3)
   mu = start_fit%params(4)
   p(1) = alpha / scale_fac
   p(2) = beta / scale_fac
   p(3) = max(delta * scale_fac, 1.0e-6_dp)
   p(4) = loc_shift + scale_fac * mu
   n = 4
case (fam_sym_vg)
   p(1) = start_fit%params(1)
   p(2) = max(start_fit%params(3) * scale_fac, 1.0e-6_dp)
   p(3) = loc_shift + scale_fac * start_fit%params(4)
   n = 3
case (fam_vg)
   p(1) = start_fit%params(1)
   p(2) = start_fit%params(2) * scale_fac
   p(3) = max(start_fit%params(3) * scale_fac, 1.0e-6_dp)
   p(4) = loc_shift + scale_fac * start_fit%params(4)
   n = 4
case default
   ok = .false.
   n = 0
end select
end subroutine native_start_from_moments

subroutine pack_unconstrained(ifam, p, n, x, ok)
integer, intent(in) :: ifam, n
real(kind=dp), intent(in) :: p(4)
real(kind=dp), intent(out) :: x(4)
logical, intent(out) :: ok
real(kind=dp) :: gamma
x = 0.0_dp
ok = .true.
select case (ifam)
case (fam_student_t)
   ok = p(1) > 2.0_dp .and. p(3) > 0.0_dp
   if (ok) x(1:3) = [log(p(1) - 2.0_dp), p(2), log(p(3))]
case (fam_jf)
   ok = minval(p(1:2)) > 0.0_dp .and. p(4) > 0.0_dp
   if (ok) x(1:4) = [log(p(1)), log(p(2)), p(3), log(p(4))]
case (fam_fs_t)
   ok = p(1) > 0.0_dp .and. p(2) > 2.0_dp .and. p(4) > 0.0_dp
   if (ok) x(1:4) = [log(p(1)), log(p(2) - 2.0_dp), p(3), log(p(4))]
case (fam_azzalini)
   ok = p(2) > 2.0_dp .and. p(4) > 0.0_dp
   if (ok) x(1:4) = [p(1), log(p(2) - 2.0_dp), p(3), log(p(4))]
case (fam_ged)
   ok = p(1) > 0.0_dp .and. p(3) > 0.0_dp
   if (ok) x(1:3) = [log(p(1)), p(2), log(p(3))]
case (fam_skew_ged)
   ok = p(1) > 0.0_dp .and. p(2) > 0.0_dp .and. p(4) > 0.0_dp
   if (ok) x(1:4) = [log(p(1)), log(p(2)), p(3), log(p(4))]
case (fam_sym_nig)
   ok = p(1) > 0.0_dp .and. p(2) > 0.0_dp
   if (ok) x(1:3) = [log(p(1)), log(p(2)), p(3)]
case (fam_nig)
   ok = p(1) > abs(p(2)) .and. p(3) > 0.0_dp
   if (ok) then
      gamma = sqrt(p(1)*p(1) - p(2)*p(2))
      x(1:4) = [log(gamma), asinh(p(2) / gamma), log(p(3)), p(4)]
   end if
case (fam_sym_vg)
   ok = p(1) > 0.0_dp .and. p(2) > 0.0_dp
   if (ok) x(1:3) = [log(p(1)), log(p(2)), p(3)]
case (fam_vg)
   ok = p(1) > 0.0_dp .and. p(3) > 0.0_dp
   if (ok) x(1:4) = [log(p(1)), p(2), log(p(3)), p(4)]
case default
   ok = .false.
end select
end subroutine pack_unconstrained

subroutine unpack_unconstrained(ifam, x, n, p, ok)
integer, intent(in) :: ifam, n
real(kind=dp), intent(in) :: x(4)
real(kind=dp), intent(out) :: p(4)
logical, intent(out) :: ok
real(kind=dp) :: gamma, v
p = 0.0_dp
ok = .true.
select case (ifam)
case (fam_student_t)
   p(1:3) = [2.0_dp + pos_from_log(x(1)), x(2), pos_from_log(x(3))]
case (fam_jf)
   p(1:4) = [pos_from_log(x(1)), pos_from_log(x(2)), x(3), pos_from_log(x(4))]
case (fam_fs_t)
   p(1:4) = [pos_from_log(x(1)), 2.0_dp + pos_from_log(x(2)), x(3), pos_from_log(x(4))]
case (fam_azzalini)
   p(1:4) = [x(1), 2.0_dp + pos_from_log(x(2)), x(3), pos_from_log(x(4))]
case (fam_ged)
   p(1:3) = [pos_from_log(x(1)), x(2), pos_from_log(x(3))]
case (fam_skew_ged)
   p(1:4) = [pos_from_log(x(1)), pos_from_log(x(2)), x(3), pos_from_log(x(4))]
case (fam_sym_nig)
   p(1:3) = [pos_from_log(x(1)), pos_from_log(x(2)), x(3)]
case (fam_nig)
   gamma = pos_from_log(x(1))
   v = x(2)
   p(1) = gamma * cosh(v)
   p(2) = gamma * sinh(v)
   p(3) = pos_from_log(x(3))
   p(4) = x(4)
case (fam_sym_vg)
   p(1:3) = [pos_from_log(x(1)), pos_from_log(x(2)), x(3)]
case (fam_vg)
   p(1:4) = [pos_from_log(x(1)), x(2), pos_from_log(x(3)), x(4)]
case default
   ok = .false.
end select
end subroutine unpack_unconstrained

subroutine default_steps(ifam, n, x0, step)
integer, intent(in) :: ifam, n
real(kind=dp), intent(in) :: x0(4)
real(kind=dp), intent(out) :: step(4)
real(kind=dp) :: s
step = 0.2_dp
s = max(g_data_scale, 1.0e-3_dp)
select case (ifam)
case (fam_student_t)
   step(1:3) = [0.3_dp, 0.2_dp*s, 0.2_dp]
case (fam_jf)
   step(1:4) = [0.3_dp, 0.3_dp, 0.2_dp*s, 0.2_dp]
case (fam_fs_t)
   step(1:4) = [0.2_dp, 0.3_dp, 0.2_dp*s, 0.2_dp]
case (fam_azzalini)
   step(1:4) = [0.5_dp, 0.3_dp, 0.2_dp*s, 0.2_dp]
case (fam_ged)
   step(1:3) = [0.3_dp, 0.2_dp*s, 0.2_dp]
case (fam_skew_ged)
   step(1:4) = [0.2_dp, 0.3_dp, 0.2_dp*s, 0.2_dp]
case (fam_sym_nig)
   step(1:3) = [0.3_dp, 0.2_dp, 0.2_dp*s]
case (fam_nig)
   step(1:4) = [0.3_dp, 0.3_dp, 0.2_dp, 0.2_dp*s]
case (fam_sym_vg)
   step(1:3) = [0.3_dp, 0.2_dp, 0.2_dp*s]
case (fam_vg)
   step(1:4) = [0.3_dp, 0.2_dp*s, 0.2_dp, 0.2_dp*s]
end select
end subroutine default_steps

subroutine print_mle_result(res)
type(mle_result_t), intent(in) :: res
integer :: i
print '(a)', trim(res%name) // ' mle'
print '(a)', repeat('-', len_trim(res%name) + 4)
write(*,'(a,l1)') 'success               = ', res%success
write(*,'(a,es14.6)') 'negloglik             = ', res%negloglik
if (res%nfev > 0) write(*,'(a,i0)') 'nfev                  = ', res%nfev
write(*,'(a,a)') 'message               = ', trim(res%message)
do i = 1, res%nparam
   write(*,'(a16,a,es20.12)') trim(res%param_names(i)) // '_start', ' = ', res%start_params(i)
end do
do i = 1, res%nparam
   if (ieee_is_nan(res%params(i))) then
      write(*,'(a16,a)') trim(res%param_names(i)), ' = nan'
   else
      write(*,'(a16,a,es20.12)') trim(res%param_names(i)), ' = ', res%params(i)
   end if
end do
print *
end subroutine print_mle_result

subroutine nelder_mead(n, x0, step, xbest, fbest, success, nfev, message)
integer, intent(in) :: n
real(kind=dp), intent(in) :: x0(4), step(4)
real(kind=dp), intent(out) :: xbest(4), fbest
logical, intent(out) :: success
integer, intent(out) :: nfev
character(len=*), intent(out) :: message
real(kind=dp) :: simp(4,5), f(5), centroid(4), xr(4), xe(4), xc(4)
real(kind=dp) :: fr, fe, fc, tolx, tolf, scale0
integer :: i, iter
logical :: ok
real(kind=dp), parameter :: alpha = 1.0_dp, gamma = 2.0_dp, rho = 0.5_dp, sigma = 0.5_dp

simp = 0.0_dp
simp(:,1) = x0
do i = 1, n
   simp(:,i+1) = x0
   simp(i,i+1) = x0(i) + merge(step(i), 0.1_dp, abs(step(i)) > 0.0_dp)
end do
nfev = 0
do i = 1, n + 1
   f(i) = objective_current(n, simp(:,i), ok)
   nfev = nfev + 1
end do

tolx = 1.0e-6_dp
tolf = 1.0e-8_dp

do iter = 1, 400
   call sort_simplex(n, simp, f)
   scale0 = max(1.0_dp, maxval(abs(simp(1:n,1))))
   if (maxval(abs(simp(1:n,2:n+1) - spread(simp(1:n,1), 2, n))) <= tolx*scale0 .and. &
       maxval(abs(f(2:n+1) - f(1))) <= tolf*(1.0_dp + abs(f(1)))) then
      exit
   end if
   centroid(1:n) = sum(simp(1:n,1:n), dim=2) / real(n, dp)
   xr(1:n) = centroid(1:n) + alpha * (centroid(1:n) - simp(1:n,n+1))
   fr = objective_current(n, xr, ok)
   nfev = nfev + 1
   if (fr < f(1)) then
      xe(1:n) = centroid(1:n) + gamma * (xr(1:n) - centroid(1:n))
      fe = objective_current(n, xe, ok)
      nfev = nfev + 1
      if (fe < fr) then
         simp(1:n,n+1) = xe(1:n)
         f(n+1) = fe
      else
         simp(1:n,n+1) = xr(1:n)
         f(n+1) = fr
      end if
   else if (fr < f(n)) then
      simp(1:n,n+1) = xr(1:n)
      f(n+1) = fr
   else
      if (fr < f(n+1)) then
         xc(1:n) = centroid(1:n) + rho * (xr(1:n) - centroid(1:n))
      else
         xc(1:n) = centroid(1:n) - rho * (centroid(1:n) - simp(1:n,n+1))
      end if
      fc = objective_current(n, xc, ok)
      nfev = nfev + 1
      if (fc < min(fr, f(n+1))) then
         simp(1:n,n+1) = xc(1:n)
         f(n+1) = fc
      else
         do i = 2, n + 1
            simp(1:n,i) = simp(1:n,1) + sigma * (simp(1:n,i) - simp(1:n,1))
            f(i) = objective_current(n, simp(:,i), ok)
            nfev = nfev + 1
         end do
      end if
   end if
end do
call sort_simplex(n, simp, f)
xbest = 0.0_dp
xbest(1:n) = simp(1:n,1)
fbest = f(1)
success = fbest < huge_nll / 10.0_dp
if (iter <= 400) then
   message = 'nelder-mead completed'
else
   message = 'maximum iterations reached'
end if
end subroutine nelder_mead

subroutine sort_simplex(n, simp, f)
integer, intent(in) :: n
real(kind=dp), intent(inout) :: simp(4,5), f(5)
integer :: i, j, m
real(kind=dp) :: ft, xt(4)
m = n + 1
do i = 1, m - 1
   do j = i + 1, m
      if (f(j) < f(i)) then
         ft = f(i)
         f(i) = f(j)
         f(j) = ft
         xt = simp(:,i)
         simp(:,i) = simp(:,j)
         simp(:,j) = xt
      end if
   end do
end do
end subroutine sort_simplex

real(kind=dp) function objective_current(n, x, ok) result(f)
integer, intent(in) :: n
real(kind=dp), intent(in) :: x(4)
logical, intent(out) :: ok
real(kind=dp) :: p(4)
logical :: ok2
call unpack_unconstrained(g_family, x, n, p, ok2)
if (.not. ok2) then
   ok = .false.
   f = huge_nll
   return
end if
select case (g_family)
case (fam_student_t)
   f = nll_student_t(p(1), p(2), p(3), ok)
case (fam_jf)
   f = nll_jf(p(1), p(2), p(3), p(4), ok)
case (fam_fs_t)
   f = nll_fs_t(p(1), p(2), p(3), p(4), ok)
case (fam_azzalini)
   f = nll_azzalini(p(1), p(2), p(3), p(4), ok)
case (fam_ged)
   f = nll_ged(p(1), p(2), p(3), ok)
case (fam_skew_ged)
   f = nll_skew_ged(p(1), p(2), p(3), p(4), ok)
case (fam_sym_nig)
   f = nll_sym_nig(p(1), p(2), p(3), ok)
case (fam_nig)
   f = nll_nig(p(1), p(2), p(3), p(4), ok)
case (fam_sym_vg)
   f = nll_sym_vg(p(1), p(2), p(3), ok)
case (fam_vg)
   f = nll_vg(p(1), p(2), p(3), p(4), ok)
case default
   ok = .false.
   f = huge_nll
end select
if (.not. ok) f = huge_nll
end function objective_current

real(kind=dp) function nll_student_t(nu, mu, scale, ok) result(f)
real(kind=dp), intent(in) :: nu, mu, scale
logical, intent(out) :: ok
integer :: i
real(kind=dp) :: lp
ok = nu > 2.0_dp .and. scale > 0.0_dp
if (.not. ok) then
   f = huge_nll
   return
end if
f = 0.0_dp
do i = 1, size(g_x)
   lp = student_t_logpdf(g_x(i), nu, mu, scale)
   if (lp <= log(tiny_pdf)) then
      ok = .false.
      f = huge_nll
      return
   end if
   f = f - lp
end do
end function nll_student_t

real(kind=dp) function nll_jf(a, b, mu, scale, ok) result(f)
real(kind=dp), intent(in) :: a, b, mu, scale
logical, intent(out) :: ok
integer :: i
real(kind=dp) :: lp
ok = a > 0.0_dp .and. b > 0.0_dp .and. scale > 0.0_dp
if (.not. ok) then
   f = huge_nll
   return
end if
f = 0.0_dp
do i = 1, size(g_x)
   lp = jf_logpdf(g_x(i), a, b, mu, scale)
   if (lp <= log(tiny_pdf)) then
      ok = .false.
      f = huge_nll
      return
   end if
   f = f - lp
end do
end function nll_jf

real(kind=dp) function nll_fs_t(xi, nu, mu, scale, ok) result(f)
real(kind=dp), intent(in) :: xi, nu, mu, scale
logical, intent(out) :: ok
integer :: i
real(kind=dp) :: lp
ok = xi > 0.0_dp .and. nu > 2.0_dp .and. scale > 0.0_dp
if (.not. ok) then
   f = huge_nll
   return
end if
f = 0.0_dp
do i = 1, size(g_x)
   lp = fs_skew_t_logpdf(g_x(i), xi, nu, mu, scale)
   if (lp <= log(tiny_pdf)) then
      ok = .false.
      f = huge_nll
      return
   end if
   f = f - lp
end do
end function nll_fs_t

real(kind=dp) function nll_azzalini(alpha, nu, mu, scale, ok) result(f)
real(kind=dp), intent(in) :: alpha, nu, mu, scale
logical, intent(out) :: ok
integer :: i
real(kind=dp) :: lp
ok = nu > 2.0_dp .and. scale > 0.0_dp
if (.not. ok) then
   f = huge_nll
   return
end if
f = 0.0_dp
do i = 1, size(g_x)
   lp = azzalini_skew_t_logpdf(g_x(i), alpha, nu, mu, scale)
   if (lp <= log(tiny_pdf)) then
      ok = .false.
      f = huge_nll
      return
   end if
   f = f - lp
end do
end function nll_azzalini

real(kind=dp) function nll_ged(p, mu, scale, ok) result(f)
real(kind=dp), intent(in) :: p, mu, scale
logical, intent(out) :: ok
integer :: i
real(kind=dp) :: lp
ok = p > 0.0_dp .and. scale > 0.0_dp
if (.not. ok) then
   f = huge_nll
   return
end if
f = 0.0_dp
do i = 1, size(g_x)
   lp = ged_logpdf(g_x(i), p, mu, scale)
   if (lp <= log(tiny_pdf)) then
      ok = .false.
      f = huge_nll
      return
   end if
   f = f - lp
end do
end function nll_ged

real(kind=dp) function nll_skew_ged(xi, p, mu, scale, ok) result(f)
real(kind=dp), intent(in) :: xi, p, mu, scale
logical, intent(out) :: ok
integer :: i
real(kind=dp) :: lp
ok = xi > 0.0_dp .and. p > 0.0_dp .and. scale > 0.0_dp
if (.not. ok) then
   f = huge_nll
   return
end if
f = 0.0_dp
do i = 1, size(g_x)
   lp = skew_ged_logpdf(g_x(i), xi, p, mu, scale)
   if (lp <= log(tiny_pdf)) then
      ok = .false.
      f = huge_nll
      return
   end if
   f = f - lp
end do
end function nll_skew_ged

real(kind=dp) function nll_sym_nig(alpha, delta, mu, ok) result(f)
real(kind=dp), intent(in) :: alpha, delta, mu
logical, intent(out) :: ok
f = nll_nig(alpha, 0.0_dp, delta, mu, ok)
end function nll_sym_nig

real(kind=dp) function nll_nig(alpha, beta, delta, mu, ok) result(f)
real(kind=dp), intent(in) :: alpha, beta, delta, mu
logical, intent(out) :: ok
integer :: i
real(kind=dp) :: lp
ok = alpha > abs(beta) .and. delta > 0.0_dp
if (.not. ok) then
   f = huge_nll
   return
end if
f = 0.0_dp
do i = 1, size(g_x)
   lp = nig_logpdf(g_x(i), alpha, beta, delta, mu)
   if (lp <= log(tiny_pdf)) then
      ok = .false.
      f = huge_nll
      return
   end if
   f = f - lp
end do
end function nll_nig

real(kind=dp) function nll_sym_vg(r, sigma, mu, ok) result(f)
real(kind=dp), intent(in) :: r, sigma, mu
logical, intent(out) :: ok
f = nll_vg(r, 0.0_dp, sigma, mu, ok)
end function nll_sym_vg

real(kind=dp) function nll_vg(r, theta, sigma, mu, ok) result(f)
real(kind=dp), intent(in) :: r, theta, sigma, mu
logical, intent(out) :: ok
integer :: i
real(kind=dp) :: lp
ok = r > 0.0_dp .and. sigma > 0.0_dp
if (.not. ok) then
   f = huge_nll
   return
end if
f = 0.0_dp
do i = 1, size(g_x)
   lp = vg_logpdf(g_x(i), r, theta, sigma, mu)
   if (lp <= log(tiny_pdf)) then
      ok = .false.
      f = huge_nll
      return
   end if
   f = f - lp
end do
end function nll_vg

pure real(kind=dp) function student_t_logpdf(x, nu, mu, scale) result(lp)
real(kind=dp), intent(in) :: x, nu, mu, scale
real(kind=dp) :: z
if (nu <= 0.0_dp .or. scale <= 0.0_dp) then
   lp = log(tiny_pdf)
   return
end if
z = (x - mu) / scale
lp = log_gamma(0.5_dp*(nu + 1.0_dp)) - log_gamma(0.5_dp*nu) - 0.5_dp*(log(pi_dp) + log(nu)) - log(scale) - &
     0.5_dp*(nu + 1.0_dp)*log(1.0_dp + z*z/nu)
end function student_t_logpdf

pure real(kind=dp) function student_t_std_logpdf(x, nu) result(lp)
real(kind=dp), intent(in) :: x, nu
real(kind=dp) :: s
if (nu <= 2.0_dp) then
   lp = log(tiny_pdf)
   return
end if
s = sqrt((nu - 2.0_dp) / nu)
lp = student_t_logpdf(x, nu, 0.0_dp, s)
end function student_t_std_logpdf

pure real(kind=dp) function ged_std_logpdf(x, p) result(lp)
real(kind=dp), intent(in) :: x, p
real(kind=dp) :: lam, axp
if (p <= 0.0_dp) then
   lp = log(tiny_pdf)
   return
end if
lam = sqrt(exp(log_gamma(1.0_dp/p) - log_gamma(3.0_dp/p)))
axp = (abs(x) / lam)**p
lp = log(p) - log(2.0_dp) - log(lam) - log_gamma(1.0_dp/p) - axp
end function ged_std_logpdf

pure real(kind=dp) function jf_logpdf_raw(x, a, b) result(lp)
real(kind=dp), intent(in) :: x, a, b
real(kind=dp) :: s, u, v
if (a <= 0.0_dp .or. b <= 0.0_dp) then
   lp = log(tiny_pdf)
   return
end if
s = sqrt(a + b + x*x)
u = 1.0_dp + x / s
v = 1.0_dp - x / s
if (u <= 0.0_dp .or. v <= 0.0_dp) then
   lp = log(tiny_pdf)
   return
end if
lp = -(a + b - 1.0_dp)*log2_dp - log_beta_local(a, b) - 0.5_dp*log(a + b) + (a + 0.5_dp)*log(u) + (b + 0.5_dp)*log(v)
end function jf_logpdf_raw

pure real(kind=dp) function jf_logpdf(y, a, b, mu, scale) result(lp)
real(kind=dp), intent(in) :: y, a, b, mu, scale
lp = jf_logpdf_raw((y - mu)/scale, a, b) - log(scale)
end function jf_logpdf

pure real(kind=dp) function fs_skew_t_logpdf(y, xi, nu, mu, scale) result(lp)
real(kind=dp), intent(in) :: y, xi, nu, mu, scale
real(kind=dp) :: z, c
if (xi <= 0.0_dp .or. nu <= 2.0_dp .or. scale <= 0.0_dp) then
   lp = log(tiny_pdf)
   return
end if
z = (y - mu) / scale
c = log(2.0_dp) - log(xi + 1.0_dp/xi)
if (z < 0.0_dp) then
   lp = c + student_t_std_logpdf(z*xi, nu) - log(scale)
else
   lp = c + student_t_std_logpdf(z/xi, nu) - log(scale)
end if
end function fs_skew_t_logpdf

pure real(kind=dp) function skew_ged_logpdf(y, xi, p, mu, scale) result(lp)
real(kind=dp), intent(in) :: y, xi, p, mu, scale
real(kind=dp) :: z, c
if (xi <= 0.0_dp .or. p <= 0.0_dp .or. scale <= 0.0_dp) then
   lp = log(tiny_pdf)
   return
end if
z = (y - mu) / scale
c = log(2.0_dp) - log(xi + 1.0_dp/xi)
if (z < 0.0_dp) then
   lp = c + ged_std_logpdf(z*xi, p) - log(scale)
else
   lp = c + ged_std_logpdf(z/xi, p) - log(scale)
end if
end function skew_ged_logpdf

pure real(kind=dp) function ged_logpdf(y, p, mu, scale) result(lp)
real(kind=dp), intent(in) :: y, p, mu, scale
if (scale <= 0.0_dp .or. p <= 0.0_dp) then
   lp = log(tiny_pdf)
   return
end if
lp = ged_std_logpdf((y - mu)/scale, p) - log(scale)
end function ged_logpdf

real(kind=dp) function azzalini_skew_t_logpdf(y, alpha, nu, mu, scale) result(lp)
real(kind=dp), intent(in) :: y, alpha, nu, mu, scale
real(kind=dp) :: z, arg, cdfv
if (nu <= 0.0_dp .or. scale <= 0.0_dp) then
   lp = log(tiny_pdf)
   return
end if
z = (y - mu) / scale
arg = alpha * z * sqrt((nu + 1.0_dp) / (nu + z*z))
cdfv = student_t_cdf(arg, nu + 1.0_dp)
if (cdfv <= 0.0_dp) then
   lp = log(tiny_pdf)
   return
end if
lp = log(2.0_dp) + student_t_logpdf(z, nu, 0.0_dp, 1.0_dp) + log(max(cdfv, tiny_pdf)) - log(scale)
end function azzalini_skew_t_logpdf

real(kind=dp) function student_t_cdf(x, nu) result(cdf)
real(kind=dp), intent(in) :: x, nu
real(kind=dp) :: z, ib
if (nu <= 0.0_dp) then
   cdf = 0.5_dp
   return
end if
if (abs(x) <= 1.0e-15_dp) then
   cdf = 0.5_dp
   return
end if
z = nu / (nu + x*x)
ib = betai(0.5_dp*nu, 0.5_dp, z)
if (x > 0.0_dp) then
   cdf = 1.0_dp - 0.5_dp*ib
else
   cdf = 0.5_dp*ib
end if
end function student_t_cdf

real(kind=dp) function betai(a, b, x) result(bt)
real(kind=dp), intent(in) :: a, b, x
if (x <= 0.0_dp) then
   bt = 0.0_dp
   return
end if
if (x >= 1.0_dp) then
   bt = 1.0_dp
   return
end if
if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
   bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a*log(x) + b*log(1.0_dp - x)) * betacf(a, b, x) / a
else
   bt = 1.0_dp - exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + b*log(1.0_dp - x) + a*log(x)) * betacf(b, a, 1.0_dp - x) / b
end if
end function betai

real(kind=dp) function betacf(a, b, x) result(cf)
real(kind=dp), intent(in) :: a, b, x
integer, parameter :: maxit = 200
real(kind=dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
integer :: m, m2
real(kind=dp) :: aa, c, d, del, qab, qam, qap
qab = a + b
qap = a + 1.0_dp
qam = a - 1.0_dp
c = 1.0_dp
d = 1.0_dp - qab*x/qap
if (abs(d) < fpmin) d = fpmin
d = 1.0_dp / d
cf = d
do m = 1, maxit
   m2 = 2*m
   aa = m*(b - m)*x / ((qam + real(m2, dp))*(a + real(m2, dp)))
   d = 1.0_dp + aa*d
   if (abs(d) < fpmin) d = fpmin
   c = 1.0_dp + aa/c
   if (abs(c) < fpmin) c = fpmin
   d = 1.0_dp / d
   cf = cf * d*c
   aa = -(a + real(m, dp))*(qab + real(m, dp))*x / ((a + real(m2, dp))*(qap + real(m2, dp)))
   d = 1.0_dp + aa*d
   if (abs(d) < fpmin) d = fpmin
   c = 1.0_dp + aa/c
   if (abs(c) < fpmin) c = fpmin
   d = 1.0_dp / d
   del = d*c
   cf = cf * del
   if (abs(del - 1.0_dp) <= eps) exit
end do
end function betacf

real(kind=dp) function nig_logpdf(x, alpha, beta, delta, mu) result(lp)
real(kind=dp), intent(in) :: x, alpha, beta, delta, mu
real(kind=dp) :: d, gamma, q, logk
if (alpha <= abs(beta) .or. delta <= 0.0_dp) then
   lp = log(tiny_pdf)
   return
end if
gamma = sqrt(alpha*alpha - beta*beta)
d = x - mu
q = sqrt(delta*delta + d*d)
logk = log_bessel_k(1.0_dp, alpha*q)
lp = log(alpha) + log(delta) - log(pi_dp) - log(q) + logk + delta*gamma + beta*d
end function nig_logpdf

real(kind=dp) function vg_logpdf(x, r, theta, sigma, mu) result(lp)
real(kind=dp), intent(in) :: x, r, theta, sigma, mu
real(kind=dp) :: k, nu, d, absd, arg, logk, c2
if (r <= 0.0_dp .or. sigma <= 0.0_dp) then
   lp = log(tiny_pdf)
   return
end if
k = 0.5_dp*r
nu = k - 0.5_dp
d = x - mu
absd = max(abs(d), 1.0e-12_dp*sigma)
c2 = sigma*sigma + theta*theta
arg = absd * sqrt(c2) / (sigma*sigma)
logk = log_bessel_k(abs(nu), arg)
lp = theta*d/(sigma*sigma) + (1.0_dp - k)*log(2.0_dp) - log_gamma(k) - 0.5_dp*log(2.0_dp*pi_dp) - log(sigma) + &
     nu*(log(absd) - 0.5_dp*log(c2)) + logk
end function vg_logpdf

real(kind=dp) function log_bessel_k(v, z) result(lv)
real(kind=dp), intent(in) :: v, z
real(kind=dp) :: vv, kval, corr
vv = abs(v)
if (z <= 0.0_dp) then
   lv = log(huge(1.0_dp))
   return
end if
if (z < 1.0e-5_dp) then
   if (vv < 1.0e-8_dp) then
      kval = -log(0.5_dp*z) - euler_gamma_dp
      lv = log(max(kval, tiny_pdf))
   else
      lv = log(0.5_dp) + log_gamma(vv) + vv*log(2.0_dp/z)
   end if
   return
end if
if (z > 50.0_dp) then
   corr = 1.0_dp + (4.0_dp*vv*vv - 1.0_dp) / (8.0_dp*z)
   lv = 0.5_dp*(log(pi_dp) - log(2.0_dp*z)) - z + log(max(corr, tiny_pdf))
   return
end if
kval = bessel_k_integral(vv, z)
lv = log(max(kval, tiny_pdf))
end function log_bessel_k

real(kind=dp) function bessel_k_integral(v, z) result(val)
real(kind=dp), intent(in) :: v, z
real(kind=dp), parameter :: h = 0.05_dp
integer, parameter :: nsub = 10, maxchunk = 800
integer :: chunk, j
real(kind=dp) :: t0, t, sumc, whole
whole = 0.0_dp
t0 = 0.0_dp
do chunk = 1, maxchunk
   sumc = scaled_k_integrand(v, z, t0) + scaled_k_integrand(v, z, t0 + real(nsub, dp)*h)
   do j = 1, nsub - 1
      t = t0 + real(j, dp)*h
      if (mod(j, 2) == 0) then
         sumc = sumc + 2.0_dp*scaled_k_integrand(v, z, t)
      else
         sumc = sumc + 4.0_dp*scaled_k_integrand(v, z, t)
      end if
   end do
   sumc = sumc * h / 3.0_dp
   whole = whole + sumc
   if (chunk > 10 .and. sumc < 1.0e-12_dp*whole) exit
   t0 = t0 + real(nsub, dp)*h
end do
val = exp(-z) * whole
end function bessel_k_integral

pure real(kind=dp) function scaled_k_integrand(v, z, t) result(val)
real(kind=dp), intent(in) :: v, z, t
val = exp(-z*(cosh(t) - 1.0_dp)) * cosh(v*t)
end function scaled_k_integrand

end module fit_dist_mle_mod

program xfit_dist
use, intrinsic :: iso_fortran_env, only: dp => real64
use fit_moments_families_mod
use fit_moments_data_mod
use fit_dist_mle_mod
implicit none
real(kind=dp), allocatable :: x(:)
character(len=256) :: filename, arg
integer :: nargs, i
logical :: file_set
real(kind=dp) :: m3_target, m4_target
type(moment_summary_t) :: raw_stats, std_stats
type(fit_result_t) :: fits(10)
type(mle_result_t) :: mlefits(10)

filename = "fs_skew_t_data.txt"
file_set = .false.
nargs = command_argument_count()
i = 1
do while (i <= nargs)
   call get_command_argument(i, arg)
   select case (trim(arg))
   case ('--file')
      if (i == nargs) error stop 'missing value after --file'
      call get_command_argument(i+1, filename)
      file_set = .true.
      i = i + 2
   case default
      if (index(arg, '--file=') == 1) then
         filename = arg(8:)
         file_set = .true.
         i = i + 1
      else if (.not. file_set .and. arg(1:2) /= '--') then
         filename = trim(arg)
         file_set = .true.
         i = i + 1
      else
         error stop 'unknown argument'
      end if
   end select
end do

call read_data_file(filename, x)
raw_stats = sample_raw_moments(x)
std_stats = standardized_moments(x)
m3_target = std_stats%m3
m4_target = std_stats%m4
call validate_target_moments(m3_target, m4_target)

write(*,'(a,a)') 'file            = ', trim(filename)
call print_data_summary(size(x), raw_stats, std_stats)

fits(1) = student_t_fit_from_kurtosis(m4_target)
fits(2) = fit_jf(m3_target, m4_target)
fits(3) = fit_skew_t_fs(m3_target, m4_target)
fits(4) = fit_azzalini_skew_t(m3_target, m4_target)
fits(5) = fit_ged_from_kurtosis(m4_target)
fits(6) = fit_skew_ged(m3_target, m4_target)
fits(7) = fit_symmetric_nig_from_kurtosis(m4_target)
fits(8) = fit_nig(m3_target, m4_target)
fits(9) = fit_symmetric_vg_from_kurtosis(m4_target)
fits(10) = fit_vg(m3_target, m4_target)

print '(a)', 'moment fits on original data scale'
print '(a)', repeat('=', 34)
do i = 1, size(fits)
   call print_original_fit(fits(i), raw_stats%mean, raw_stats%var)
end do

print '(a)', 'maximum-likelihood fits started from moment estimates'
print '(a)', repeat('=', 52)
do i = 1, size(fits)
   call fit_family_mle(i, fits(i), raw_stats%mean, raw_stats%var, x, mlefits(i))
   call print_mle_result(mlefits(i))
end do

end program xfit_dist
