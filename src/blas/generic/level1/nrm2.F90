! SPDX-License-Identifier: LGPL-3.0-or-later

#ifdef LPF_FP8_E5M2
submodule (lpf_blas_fp8_e5m2) lpf_blas_nrm2_fp8_e5m2
    use lpf_fp8_e5m2
#endif
#ifdef LPF_FP8_E4M3
submodule (lpf_blas_fp8_e4m3) lpf_blas_nrm2_fp8_e4m3
    use lpf_fp8_e4m3
#endif
    use iso_fortran_env
    use iso_c_binding

contains

    pure function nrm2_impl(n,x,incx) result(out)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), intent(in) :: x(*)
        type(DT) :: out

        integer(int64) :: ix
        type(DT) :: norm , scale, ssq, absxi, one, zero

        one  = 1.0
        zero = 0.0

        if (n.lt.1 .or. incx.lt.1) then
            norm = 0.0
        else if (n.eq.1) then
            norm = abs(x(1))
        else
            scale = 0.0
            ssq = 1.0
            !        the following loop is equivalent to this call to the lapack
            !        auxiliary routine:
            !        call slassq( n, x, incx, scale, ssq )
            !
            do ix = 1,1 + (n-1)*incx,incx
                if (x(ix).ne.zero) then
                    absxi = abs(x(ix))
                    if (scale.lt.absxi) then
                        ssq = one + ssq* (scale/absxi)**2
                        scale = absxi
                    else
                        ssq = ssq + (absxi/scale)**2
                    end if
                end if
            end do
            norm = scale*sqrt(ssq)
        end if

        out = norm
        return
    end function

    module function nrm2_64(n,sx,incx) result(out)
        implicit none
        integer(int64), intent(in) :: incx, n
        type(DT), target, intent(in) :: sx(..)
        type(DT) :: out

        type(DT), CONTIGUOUS, pointer :: px(:)
        integer(int64) :: len
        integer(int64) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, px, lx)

        out = nrm2_impl(n, px, incx)
    end function

    module function nrm2_32(n,sx,incx) result(out)
        implicit none
        integer(int32), intent(in) :: incx, n
        type(DT), target, intent(in) :: sx(..)
        type(DT) :: out

        type(DT), CONTIGUOUS, pointer :: px(:)
        integer(int32) :: len
        integer(int32) :: lx(1)
        type(c_ptr) :: ptr

        len = 1 + ( (n - 1) * abs(incx))
        lx(1) = len
        ptr = c_loc(sx)
        call c_f_pointer(ptr, px, lx)

        out = nrm2_impl(int(n, int64), px, int(incx, int64))
    end function

end submodule
