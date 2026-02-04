submodule(lpf_lapack_fp16) lpf_lapack_gecholqr_shift__fp16

contains
        !> \brief \b DGECHOLQR3
    !
    !  Definition:
    !  ===========
    !
    !       SUBROUTINE DGECHOLQR3( M, N, A, LDA, R, LDR, WORK, INFO )
    !
    !       .. Scalar Arguments ..
    !       INTEGER INFO, M, N
    !       ..
    !       .. Array Arguments ..
    !       DOUBLE PRECISION A( LDA, * ), R( LDR, * ), WORK( * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> DGECHOLQR3 computes a QR factorization of a real M-by-N matrix A
    !> using the cholQR3 algorithm from
    !>
    !>    Shifted CholeskyQR for computing the QR factorization of ill-conditioned matrices
    !>    Fukaya, Takeshi; Kannan, Ramaseshan; Nakatsukasa, Yuji; Yamamoto, Yusaku; Yanagisawa, Yuka
    !>    2018, http://arxiv.org/pdf/1809.11085v1
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] M
    !> \verbatim
    !>          M is INTEGER
    !>          The number of rows of the matrix A.  M >= 0.
    !> \endverbatim
    !>
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The number of columns of the matrix A. M >= N >= 0.
    !> \endverbatim
    !>
    !> \param[in,out] A
    !> \verbatim
    !>          A is DOUBLE PRECISION array, dimension (LDA,N)
    !>          On entry, the M-by-N matrix A.
    !>          On exit, it contains the matrix Q such that A = Q * R.
    !> \endverbatim
    !>
    !> \param[in] LDA
    !> \verbatim
    !>          LDA is INTEGER
    !>          The leading dimension of the array A.  LDA >= max(1,M).
    !> \endverbatim
    !>
    !> \param[out] R
    !> \verbatim
    !>          R is DOUBLE PRECISION array, dimension (LDR,N)
    !>          On exit, R contains the upper triangular matrix R of A = Q * R.
    !> \endverbatim
    !>
    !> \param[in] LDR
    !> \verbatim
    !>          LDR is INTEGER
    !>          The leading dimension of the array R.  LDR >= max(1,N).
    !> \endverbatim
    !>
    !> \param[in,out] WORK
    !> \verbatim
    !>          WORK is DOUBLE PRECISION array, dimension (M*N)
    !>          WORKSPACE
    !> \endverbatim
    !>
    !> \param[out] INFO
    !> \verbatim
    !>          INFO is INTEGER
    !>          = 0:  successful exit
    !>          < 0:  if INFO = -i, the i-th argument had an illegal value
    !>          > 0:  if INFO = i, the Cholesky decomposition failed in iteration i
    !> \endverbatim
    !
    !  Authors:
    !  ========
    !
    !> \author Martin Koehler
    !> \date December 2019
    !
    !> \ingroup cholqr
    !
    module subroutine gecholqr_shift(m, n, a, lda, r, ldr, work, info)
        implicit none
        ! arguments
        integer(lpf_default_int_kind), intent(in)   :: m, n, lda, ldr
        integer(lpf_default_int_kind), intent(inout) :: info
        type(fp16), intent(inout) :: a(lda, *), r(ldr, *), work(*)

        ! locals
        type(fp16) :: shift, eps
        integer(lpf_default_int_kind) :: it, iinfo, k
        type(fp16) :: one, zero

        ! check arguments
        one = 1.0
        zero = 0.0
        info = 0

        if ( m .lt. 0 ) then
            info = -1
        else if ( n .lt. 0 .or. n .gt. m) then
            info = -2
        else if ( lda .lt. max(1, m)) then
            info = -4
        else if ( ldr .lt. max(1, n)) then
            info = -6
        end if

        if ( info .ne. 0 ) then
            call lpf_blas_xerbla("hgecholqr_shift", -info)
            return
        end if

        ! quick return
        if ( m.eq.0 .or. n.eq.0 ) then
            return
        endif


        ! compute the shift
        eps = lamch("eps")
        shift = lange("f", m, n, a, lda, work)

        !
        !  fp16(m) * fp16(n) overflows easily !
        !
        shift = fp16(11.0)* (fp16(m)*fp16(n) + fp16(n*(n+1))) * eps * (shift**2)


        r(1:n, 1:n) = 0.0
        do k = 1, n
            r(k,k) = 1.0
        end do

        do it = 1, 3
            ! call hsyrk( "upper", "trans", n, m, one, a, lda, zero, work, n)

            call hgemm("T", "N", n, n, m, one, a, lda, a, lda, zero, work, n)

            if ( it .le. 1 ) then
                do k = 1, n
                    work(k+(k-1)*n) = work(k+(k-1)*n) + shift
                end do
            end if

            call potrf("upper", n, work, n, iinfo)
            if ( iinfo .ne. 0 ) then
                info = it
                return
            end if

            call htrsm("right", "upper", "notrans", "nounit", m, n, one, work, n, a, lda )
            if (it .eq. 1) then
                call hlacpy("u", n, n, work, n, r, ldr)
            else
                call htrmm("left", "upper", "notrans", "nounit", n, n, one, work, n, r, ldr)
            endif

        end do

        return


    end subroutine
end submodule

