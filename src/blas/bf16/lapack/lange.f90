! SPDX-License-Identifier: LGPL-3.0-or-later
submodule (lpf_lapack_bf16) lpf_lapack_lange_bf16
    use lpf_bf16
    use lpf_types
    implicit none

contains

    !> \brief \b lange returns the value of the 1-norm, Frobenius norm, inf
    !
    !
    !  Definition:
    !  ===========
    !
    !       type(bf16)             FUNCTION lange( NORM, M, N, A, LDA, WORK )
    !
    !       .. Scalar Arguments ..
    !       CHARACTER          NORM
    !       INTEGER            LDA, M, N
    !       ..
    !       .. Array Arguments ..
    !       type(bf16)               A( LDA, * ), WORK( * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> lange  returns the value of the one norm,  or the Frobenius norm, or
    !> the  infinity norm,  or the  element of  largest absolute value  of a
    !> real matrix A.
    !> \endverbatim
    !>
    !> \return lange
    !> \verbatim
    !>
    !>    lange = ( max(abs(A(i,j))), NORM = 'M' or 'm'
    !>             (
    !>             ( norm1(A),         NORM = '1', 'O' or 'o'
    !>             (
    !>             ( normI(A),         NORM = 'I' or 'i'
    !>             (
    !>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
    !>
    !> where  norm1  denotes the  one norm of a matrix (maximum column sum),
    !> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
    !> normF  denotes the Frobenius norm of a matrix (square root of sum of
    !> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix no
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] NORM
    !> \verbatim
    !>          NORM is CHARACTER*1
    !>          Specifies the value to be returned in lange as described
    !>          above.
    !> \endverbatim
    !>
    !> \param[in] M
    !> \verbatim
    !>          M is INTEGER
    !>          The number of rows of the matrix A.  M >= 0.  When M = 0,
    !>          lange is set to zero.
    !> \endverbatim
    !>
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The number of columns of the matrix A.  N >= 0.  When N = 0,
    !>          lange is set to zero.
    !> \endverbatim
    !>
    !> \param[in] A
    !> \verbatim
    !>          A is type(bf16) array, dimension (LDA,N)
    !>          The m by n matrix A.
    !> \endverbatim
    !>
    !> \param[in] LDA
    !> \verbatim
    !>          LDA is INTEGER
    !>          The leading dimension of the array A.  LDA >= max(M,1).
    !> \endverbatim
    !>
    !> \param[out] WORK
    !> \verbatim
    !>          WORK is type(bf16) array, dimension (MAX(1,LWORK)),
    !>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
    !>          referenced.
    !> \endverbatim
    !
    module function lange( norm, m, n, a, lda, work ) result(nrm)
        character, intent(in) :: norm
        integer(lpf_default_int_kind), intent(in)   :: lda, m, n
        type(bf16), intent(in)      :: a( lda, * )
        type(bf16), intent(inout)   :: work( * )
        type(bf16) :: nrm
        !     ..
        !
        !     .. parameters ..
        type(bf16) ::            one, zero
        !     .. local scalars ..
        integer(lpf_default_int_kind)            :: i, j
        integer(lpf_default_int_kind), parameter :: ione = 1

        real(real32) :: sca, sum
        type(bf16) :: value, temp
        !     ..

        one  = 1.0e+0
        zero = 0.0e+0
        if( min( m, n ).eq.0 ) then
            value = zero
        else if( lsame( norm, 'm' ) ) then
            !
            !        find max(abs(a(i,j))).
            !
            value = zero
            do j = 1, n
                do i = 1, m
                    temp = abs( a( i, j ) )
                    if( value.lt.temp .or. isnan( temp ) ) value = temp
                end do
            end do
        else if( ( lsame( norm, 'o' ) ) .or. ( norm.eq.'1' ) ) then
            !
            !        find norm1(a).
            !
            value = zero
            do j = 1, n
                sum = zero
                do i = 1, m
                    sum = sum + abs( a( i, j ) )
                end do
                if( value.lt.sum .or. isnan( sum ) ) value = sum
            end do
        else if( lsame( norm, 'i' ) ) then
            !
            !        find normi(a).
            !
            do  i = 1, m
                work( i ) = zero
            end do
            do j = 1, n
                do i = 1, m
                    work( i ) = work( i ) + abs( a( i, j ) )
                end do
            end do

            value = zero
            do i = 1, m
                temp = work( i )
                if( value.lt.temp .or. isnan( temp ) ) value = temp
            end do
        else if( ( lsame( norm, 'f' ) ) .or. ( lsame( norm, 'e' ) ) ) then
            !
            !        find normf(a).
            !
            sca = 0.0
            sum = 1.0
            do  j = 1, n
                call lassq( m, a( 1, j ), ione, sca, sum )
            end do
            value = bf16(sca*sqrt( sum ))
        end if
        !
        nrm = value
        return
    end function

    !> \brief \b lassq updates a sum of squares represented in scaled form.
    !
    !
    !  Definition:
    !  ===========
    !
    !       SUBROUTINE lassq( N, X, INCX, SCALE, SUMSQ )
    !
    !       .. Scalar Arguments ..
    !       INTEGER            INCX, N
    !       REAL               SCALE, SUMSQ
    !       ..
    !       .. Array Arguments ..
    !       type(bf16)               X( * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> lassq  returns the values  scl  and  smsq  such that
    !>
    !>    ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
    !>
    !> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
    !> assumed to be non-negative and  scl  returns the value
    !>
    !>    scl = max( scale, abs( x( i ) ) ).
    !>
    !> scale and sumsq must be supplied in SCALE and SUMSQ and
    !> scl and smsq are overwritten on SCALE and SUMSQ respectively.
    !>
    !> The routine makes only one pass through the vector x.
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The number of elements to be used from the vector X.
    !> \endverbatim
    !>
    !> \param[in] X
    !> \verbatim
    !>          X is type(bf16) array, dimension (N)
    !>          The vector for which a scaled sum of squares is computed.
    !>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
    !> \endverbatim
    !>
    !> \param[in] INCX
    !> \verbatim
    !>          INCX is INTEGER
    !>          The increment between successive values of the vector X.
    !>          INCX > 0.
    !> \endverbatim
    !>
    !> \param[in,out] SCALE
    !> \verbatim
    !>          SCALE is REAL
    !>          On entry, the value  scale  in the equation above.
    !>          On exit, SCALE is overwritten with  scl , the scaling factor
    !>          for the sum of squares.
    !> \endverbatim
    !>
    !> \param[in,out] SUMSQ
    !> \verbatim
    !>          SUMSQ is REAL
    !>          On entry, the value  sumsq  in the equation above.
    !>          On exit, SUMSQ is overwritten with  smsq , the basic sum of
    !>          squares from which  scl  has been factored out.
    !> \endverbatim
    !
    module subroutine lassq( n, x, incx, sca, sumsq )
        integer(lpf_default_int_kind), intent(in) ::   incx, n
        real(real32), intent(inout) :: sca, sumsq
        !     ..
        !     .. array arguments ..
        type(bf16), intent(in)  :: x( * )
        !     ..
        !
        ! =====================================================================
        !
        !     .. parameters ..
        real(real32) ::    zero
        !     ..
        !     .. local scalars ..
        integer(lpf_default_int_kind) :: ix
        real(real32) :: absxi
        !     ..

        zero = 0.0
        if( n.gt.0 ) then
            do ix = 1, 1 + ( n-1 )*incx, incx
                absxi = abs( real(x( ix ) ))
                if( absxi.gt.zero.or.isnan( absxi ) ) then
                    if( sca.lt.absxi ) then
                        sumsq = 1.0 + sumsq*( sca / absxi )**2
                        sca = absxi
                    else
                        sumsq = sumsq + ( absxi / sca )**2
                    end if
                end if
            end do
        end if

        return
    end subroutine


end submodule
