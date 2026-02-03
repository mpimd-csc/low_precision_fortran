submodule(lpf_lapack_fp16) lpf_lapack_geqrt2v_fp16
contains
    !> \brief \b HGEQRT2 computes a QR factorization of a general real or co
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> Download HGEQRT2 + dependencies
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&fil
    !> [TGZ]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&fil
    !> [ZIP]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&fil
    !> [TXT]</a>
    !
    !  Definition:
    !  ===========
    !
    !       SUBROUTINE HGEQRT2( M, N, A, LDA, T, LDT, INFO )
    !
    !       .. Scalar Arguments ..
    !       INTEGER   INFO, LDA, LDT, M, N
    !       ..
    !       .. Array Arguments ..
    !       TYPE(FP16)   A( LDA, * ), T( LDT, * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> HGEQRT2 computes a QR factorization of a real M-by-N matrix A,
    !> using the compact WY representation of Q.
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] M
    !> \verbatim
    !>          M is INTEGER
    !>          The number of rows of the matrix A.  M >= N.
    !> \endverbatim
    !>
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The number of columns of the matrix A.  N >= 0.
    !> \endverbatim
    !>
    !> \param[in,out] A
    !> \verbatim
    !>          A is TYPE(FP16) array, dimension (LDA,N)
    !>          On entry, the real M-by-N matrix A.  On exit, the elements o
    !>          above the diagonal contain the N-by-N upper triangular matri
    !>          elements below the diagonal are the columns of V.  See below
    !>          further details.
    !> \endverbatim
    !>
    !> \param[in] LDA
    !> \verbatim
    !>          LDA is INTEGER
    !>          The leading dimension of the array A.  LDA >= max(1,M).
    !> \endverbatim
    !>
    !> \param[out] T
    !> \verbatim
    !>          T is TYPE(FP16) array, dimension (LDT,N)
    !>          The N-by-N upper triangular factor of the block reflector.
    !>          The elements on and above the diagonal contain the block
    !>          reflector T; the elements below the diagonal are not used.
    !>          See below for further details.
    !> \endverbatim
    !>
    !> \param[in] LDT
    !> \verbatim
    !>          LDT is INTEGER
    !>          The leading dimension of the array T.  LDT >= max(1,N).
    !> \endverbatim
    !>
    !> \param[out] INFO
    !> \verbatim
    !>          INFO is INTEGER
    !>          = 0: successful exit
    !>          < 0: if INFO = -i, the i-th argument had an illegal value
    !> \endverbatim
    !
    !  Authors:
    !  ========
    !
    !> \author Univ. of Tennessee
    !> \author Univ. of California Berkeley
    !> \author Univ. of Colorado Denver
    !> \author NAG Ltd.
    !
    !> \ingroup geqrt2
    !
    !> \par Further Details:
    !  =====================
    !>
    !> \verbatim
    !>
    !>  The matrix V stores the elementary reflectors H(i) in the i-th colum
    !>  below the diagonal. For example, if M=5 and N=3, the matrix V is
    !>
    !>               V = (  1       )
    !>                   ( v1  1    )
    !>                   ( v1 v2  1 )
    !>                   ( v1 v2 v3 )
    !>                   ( v1 v2 v3 )
    !>
    !>  where the vi's represent the vectors which define H(i), which are re
    !>  in the matrix A.  The 1's along the diagonal of V are not stored in
    !>  block reflector H is then given by
    !>
    !>               H = I - V * T * V**T
    !>
    !>  where V**T is the transpose of V.
    !> \endverbatim
    !>
    !  =====================================================================
    module subroutine geqrt2v( norm, m, n, a, lda, diagr, t, ldt, info )
        integer(lpf_default_int_kind), intent(in) :: lda, ldt, m, n
        integer(lpf_default_int_kind), intent(inout) :: info
        type(fp16), intent(inout) :: a( lda, * ), diagr(*),  t( ldt, * )
        character(len=*), intent(in)  :: norm

        type(fp16)  :: one, zero
        integer(lpf_default_int_kind) ::  i, k
        type(fp16)  ::  aii, alpha
        !     .. external subroutines ..

        one = 1.0
        zero = 0.0
        info = 0
        if( n.lt.0 ) then
            info = -2
        else if( m.lt.n ) then
            info = -1
        else if( lda.lt.max( 1, m ) ) then
            info = -4
        else if( ldt.lt.max( 1, n ) ) then
            info = -6
        end if
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'geqrt2v', -info )
            return
        end if
        !
        k = min( m, n )
        !
        do i = 1, k
            !
            !        generate elem. refl. h(i) to annihilate a(i+1:m,i), tau(i) -> t
            !
            diagr(i) = a(i,i)
            call larfgv1( norm, m-i+1, diagr(i), a( i, i ), a( min( i+1, m ), i ), 1_lpf_default_int_kind,       &
                &                t( i, 1 ) )

            if( i.lt.n ) then
                !
                !           apply h(i) to a(i:m,i+1:n) from the left
                !
                aii = a( i, i )
                a( i, i ) = one
                !
                !           w(1:n-i) := a(i:m,i+1:n)^h * a(i:m,i) [w = t(:,n)]
                !
                call hgemv( 't',m-i+1, n-i, one, a( i, i+1 ), lda,          &
                    &                  a( i, i ), 1_lpf_default_int_kind, zero, t( 1, n ), 1_lpf_default_int_kind )
                !
                !           a(i:m,i+1:n) = a(i:m,i+1:n) + alpha*a(i:m,i)*w(1:n-1)^h
                !
                alpha = -(t( i, 1 ))
                call hger( m-i+1, n-i, alpha, a( i, i ), 1_lpf_default_int_kind,                 &
                    &           t( 1, n ), 1_lpf_default_int_kind, a( i, i+1 ), lda )
                a( i, i ) = aii
            end if
        end do
        !
        do i = 2, n
            aii = a( i, i )
            a( i, i ) = one
            !
            !        t(1:i-1,i) := alpha * a(i:m,1:i-1)**t * a(i:m,i)
            !
            alpha = -t( i, 1 )
            call hgemv( 't', m-i+1, i-1, alpha, a( i, 1 ), lda,            &
                &               a( i, i ), 1_lpf_default_int_kind, zero, t( 1, i ), 1_lpf_default_int_kind)
            a( i, i ) = aii
            !
            !        t(1:i-1,i) := t(1:i-1,1:i-1) * t(1:i-1,i)
            !
            call htrmv( 'u', 'n', 'n', i-1, t, ldt, t( 1, i ), 1_lpf_default_int_kind)
            !
            !           t(i,i) = tau(i)
            !
            t( i, i ) = t( i, 1 )
            t( i, 1) = zero
        end do

        !
        !     end of geqrt2
        !
    end subroutine

end submodule
