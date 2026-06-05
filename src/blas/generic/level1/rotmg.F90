! SPDX-License-Identifier: LGPL-3.0-or-later
!> \brief \b SROTMG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SROTMG(SD1,SD2,SX1,SY1,SPARAM)
!
!       .. Scalar Arguments ..
!       REAL SD1,SD2,SX1,SY1
!       ..
!       .. Array Arguments ..
!       type(DT) SPARAM(5)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
!>    THE SECOND COMPONENT OF THE 2-VECTOR
!>         (SQRT(SD1)*SX1,SQRT(SD2)*SY2)**T
!>    WITH SPARAM(1)=SFLAG.
!>
!>    H HAS ONE OF THE FOLLOWING FORMS:
!>
!>      SFLAG=-1.E0     SFLAG=0.E0      SFLAG=1.E0      SFLAG=-2.E0
!>
!>      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
!>    H=(          )    (          )    (          )    (          )
!>      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
!>
!>    LOCATIONS 2-4 OF SPARAM CONTAIN SH11, SH21, SH12, AND SH22
!>    RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE
!>    VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.)
!>
!>    THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
!>    INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
!>    OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] SD1
!> \verbatim
!>          SD1 is REAL
!> \endverbatim
!>
!> \param[in,out] SD2
!> \verbatim
!>          SD2 is REAL
!> \endverbatim
!>
!> \param[in,out] SX1
!> \verbatim
!>          SX1 is REAL
!> \endverbatim
!>
!> \param[in] SY1
!> \verbatim
!>          SY1 is REAL
!> \endverbatim
!>
!> \param[out] SPARAM
!> \verbatim
!>          SPARAM is type(DT) array, dimension (5)
!>     SPARAM(1)=SFLAG
!>     SPARAM(2)=SH11
!>     SPARAM(3)=SH21
!>     SPARAM(4)=SH12
!>     SPARAM(5)=SH22
!> \endverbatim
!  =====================================================================
#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_rotmg_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_rotmg_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure subroutine rotmg_impl(sd1,sd2,sx1,sy1,sparams)
        implicit none
        type(DT), intent(inout) :: sd1, sd2, sx1
        type(DT), intent(in) :: sy1
        type(DT), intent(out) :: sparams(*)

        type(DT) :: gam, gamsq, one, rgamsq, sflag, sh11, sh12, sh21, sh22, sp1, sp2, sq1, sq2, stemp, su, two, zero

        zero = 0.0
        one = 1.0
        two = 2.0
        gam = 16.0
        gamsq = 256.0
        rgamsq = 3.90625E-3

        if (sd1.lt.zero) then
            sflag = -one
            sh11 = zero
            sh12 = zero
            sh21 = zero
            sh22 = zero
            sd1 = zero
            sd2 = zero
            sx1 = zero
        else
            sp2 = sd2*sy1
            if (sp2.eq.zero) then
                sflag = -two
                sparams(1) = sflag
                return
            end if
            sp1 = sd1*sx1
            sq2 = sp2*sy1
            sq1 = sp1*sx1

            if (abs(sq1).gt.abs(sq2)) then
                sh21 = -sy1/sx1
                sh12 = sp2/sp1
                su = one - sh12*sh21
                if (su.gt.zero) then
                    sflag = zero
                    sd1 = sd1/su
                    sd2 = sd2/su
                    sx1 = sx1*su
                else
                    sflag = -one
                    sh11 = zero
                    sh12 = zero
                    sh21 = zero
                    sh22 = zero
                    sd1 = zero
                    sd2 = zero
                    sx1 = zero
                end if
            else
                if (sq2.lt.zero) then
                    sflag = -one
                    sh11 = zero
                    sh12 = zero
                    sh21 = zero
                    sh22 = zero
                    sd1 = zero
                    sd2 = zero
                    sx1 = zero
                else
                    sflag = one
                    sh11 = sp1/sp2
                    sh22 = sx1/sy1
                    su = one + sh11*sh22
                    stemp = sd2/su
                    sd2 = sd1/su
                    sd1 = stemp
                    sx1 = sy1*su
                end if
            end if

            if (sd1.ne.zero) then
                do while ((sd1.le.rgamsq) .or. (sd1.ge.gamsq))
                    if (sflag.eq.zero) then
                        sh11 = one
                        sh22 = one
                        sflag = -one
                    else
                        sh21 = -one
                        sh12 = one
                        sflag = -one
                    end if
                    if (sd1.le.rgamsq) then
                        sd1 = sd1*gam**2
                        sx1 = sx1/gam
                        sh11 = sh11/gam
                        sh12 = sh12/gam
                    else
                        sd1 = sd1/gam**2
                        sx1 = sx1*gam
                        sh11 = sh11*gam
                        sh12 = sh12*gam
                    end if
                end do
            end if

            if (sd2.ne.zero) then
                do while ( (abs(sd2).le.rgamsq) .or. (abs(sd2).ge.gamsq) )
                    if (sflag.eq.zero) then
                        sh11 = one
                        sh22 = one
                        sflag = -one
                    else
                        sh21 = -one
                        sh12 = one
                        sflag = -one
                    end if
                    if (abs(sd2).le.rgamsq) then
                        sd2 = sd2*gam**2
                        sh21 = sh21/gam
                        sh22 = sh22/gam
                    else
                        sd2 = sd2/gam**2
                        sh21 = sh21*gam
                        sh22 = sh22*gam
                    end if
                end do
            end if
        end if

        if (sflag.lt.zero) then
            sparams(2) = sh11
            sparams(3) = sh21
            sparams(4) = sh12
            sparams(5) = sh22
        else if (sflag.eq.zero) then
            sparams(3) = sh21
            sparams(4) = sh12
        else
            sparams(2) = sh11
            sparams(5) = sh22
        end if
        sparams(1) = sflag
    end subroutine

    module subroutine rotmg(sa,sb,sc,sd,sparam)
        implicit none
        type(DT), intent(inout) :: sa, sb, sc
        type(DT), intent(inout) :: sd
        type(DT), target, intent(out) :: sparam(..)

        type(DT), CONTIGUOUS, pointer :: psparams(:)
        integer(int32) :: lenp
        integer(int32) :: lp(1)
        type(c_ptr) :: ptr

        lenp = 5
        lp(1) = lenp
        ptr = c_loc(sparam)
        call c_f_pointer(ptr, psparams, lp)

        call rotmg_impl(sa, sb, sc, sd, psparams)
    end subroutine

end submodule
