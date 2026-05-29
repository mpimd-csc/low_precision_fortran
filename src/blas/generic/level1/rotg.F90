!> \brief \b SROTG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SROTG constructs a plane rotation
!>    [  c  s ] [ a ] = [ r ]
!>    [ -s  c ] [ b ]   [ 0 ]
!> satisfying c**2 + s**2 = 1.
!>
!> The computation uses the formulas
!>    sigma = sgn(a)    if |a| >  |b|
!>          = sgn(b)    if |b| >= |a|
!>    r = sigma*sqrt( a**2 + b**2 )
!>    c = 1; s = 0      if r = 0
!>    c = a/r; s = b/r  if r != 0
!> The subroutine also computes
!>    z = s    if |a| > |b|,
!>      = 1/c  if |b| >= |a| and c != 0
!>      = 1    if c = 0
!> This allows c and s to be reconstructed from z as follows:
!>    If z = 1, set c = 0, s = 1.
!>    If |z| < 1, set c = sqrt(1 - z**2) and s = z.
!>    If |z| > 1, set c = 1/z and s = sqrt( 1 - c**2).
!>
!> \endverbatim
!>
!> @see lartg, @see lartgp
!
!  Arguments:
!  ==========
!
!> \param[in,out] A
!> \verbatim
!>          A is REAL
!>          On entry, the scalar a.
!>          On exit, the scalar r.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL
!>          On entry, the scalar b.
!>          On exit, the scalar z.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL
!>          The scalar c.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL
!>          The scalar s.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!
!> \ingroup rotg
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!> \endverbatim
!
!  =====================================================================
#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_rotg_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_rotg_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine rotg_impl(a,b,c,s)
        implicit none
        type(DT), intent(inout) :: a, b
        type(DT), intent(out) :: c, s

        type(DT) :: anorm, bnorm, scl, sigma, r, z
        type(DT) :: zero, one, safmin, safmax

        zero = 0.0
        one = 1.0
        safmin = DT(radix(zero))**max(minexponent(zero)-1, 1-maxexponent(zero))
        safmax = DT(radix(zero))**max(1-minexponent(zero), maxexponent(zero)-1)

        anorm = abs(a)
        bnorm = abs(b)
        if( bnorm == zero ) then
            c = one
            s = zero
            b = zero
        else if( anorm == zero ) then
            c = zero
            s = one
            a = b
            b = one
        else
            scl = min( safmax, max( safmin, anorm, bnorm ) )
            if( anorm > bnorm ) then
                sigma = sign(one, a)
            else
                sigma = sign(one, b)
            end if
            r = sigma*( scl*sqrt((a/scl)**2 + (b/scl)**2) )
            c = a/r
            s = b/r

            if( anorm > bnorm ) then
                z = s
            else if( c /= zero ) then
                z = one/c
            else
                z = one
            end if
            a = r
            b = z
        end if
    end subroutine

    module subroutine rotg(sa,sb,sc,ss)
        implicit none
        type(DT), intent(inout) :: sa, sb
        type(DT), intent(out) :: sc, ss

        call rotg_impl(sa, sb, sc, ss)
    end subroutine

end submodule
