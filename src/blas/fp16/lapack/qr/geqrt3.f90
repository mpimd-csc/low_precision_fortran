submodule(lpf_lapack_fp16) lpf_lapack_geqrt3_fp16
contains
    !> \brief \b HGEQRT3 recursively computes a QR factorization of a genera
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> Download HGEQRT3 + dependencies
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
    !       RECURSIVE SUBROUTINE HGEQRT3( M, N, A, LDA, T, LDT, INFO )
    !
    !       .. Scalar Arguments ..
    !       INTEGER   INFO, LDA, M, N, LDT
    !       ..
    !       .. Array Arguments ..
    !       REAL   A( LDA, * ), T( LDT, * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> HGEQRT3 recursively computes a QR factorization of a real M-by-N
    !> matrix A, using the compact WY representation of Q.
    !>
    !> Based on the algorithm of Elmroth and Gustavson,
    !> IBM J. Res. Develop. Vol 44 No. 4 July 2000.
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
    !>          A is REAL array, dimension (LDA,N)
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
    !>          T is REAL array, dimension (LDT,N)
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
    !> \ingroup geqrt3
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
    !>
    !>  For details of the algorithm, see Elmroth and Gustavson (cited above
    !> \endverbatim
    !>
    !  =====================================================================
    module recursive subroutine geqrt3( norm, m, n, a, lda, t, ldt, info )
        integer(lpf_default_int_kind), intent(in) :: lda, m, n, ldt
        integer(lpf_default_int_kind), intent(out) :: info

      character(len=*), intent(in)  :: norm
        !     ..
        !     .. array arguments ..
        type(fp16), intent(inout) ::  a( lda, * ), t( ldt, * )
        !     ..
        !
        !  =====================================================================
        !
        !     ..h parameters ..
        type(fp16) ::   one
        !     ..
        !     .. local scalars ..
        integer(lpf_default_int_kind) ::  i, i1, j, j1, n1, n2, iinfo
        !     ..

        one = 1.0
        info = 0
        if( n .lt. 0 ) then
            info = -2
        else if( m .lt. n ) then
            info = -1
        else if( lda .lt. max( 1, m ) ) then
            info = -4
        else if( ldt .lt. max( 1, n ) ) then
            info = -6
        end if
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'geqrt3', -info )
            return
        end if
        !
        if( n.eq.1 ) then
            !
            !        compute householder transform when n=1
            !
            call larfg( norm, m, a(1,1), a( min( 2, m ), 1 ), 1, t(1,1) )
            !
        else
            !
            !        otherwise, split a into blocks...
            !
            n1 = n/2
            n2 = n-n1
            j1 = min( n1+1, n )
            i1 = min( n+1, m )
            !
            !        compute a(1:m,1:n1) <- (y1,r1,t1), where q1 = i - y1 t1 y1^h
            !
            call geqrt3( norm, m, n1, a, lda, t, ldt, iinfo )
            !
            !        compute a(1:m,j1:n) = q1^h a(1:m,j1:n) [workspace: t(1:n1,j1:n)
            !
            do j=1,n2
                do i=1,n1
                    t( i, j+n1 ) = a( i, j+n1 )
                end do
            end do
            call htrmm( 'l', 'l', 't', 'u', n1, n2, one,                   &
                &               a, lda, t( 1, j1 ), ldt )
            !
            call hgemm( 't', 'n', n1, n2, m-n1, one, a( j1, 1 ), lda,      &
                &               a( j1, j1 ), lda, one, t( 1, j1 ), ldt)
            !
            call htrmm( 'l', 'u', 't', 'n', n1, n2, one,                   &
                &               t, ldt, t( 1, j1 ), ldt )
            !
            call hgemm( 'n', 'n', m-n1, n2, n1, -one, a( j1, 1 ), lda,     &
                &               t( 1, j1 ), ldt, one, a( j1, j1 ), lda )
            !
            call htrmm( 'l', 'l', 'n', 'u', n1, n2, one,                   &
                &               a, lda, t( 1, j1 ), ldt )
            !
            do j=1,n2
                do i=1,n1
                    a( i, j+n1 ) = a( i, j+n1 ) - t( i, j+n1 )
                end do
            end do
            !
            !        compute a(j1:m,j1:n) <- (y2,r2,t2) where q2 = i - y2 t2 y2^h
            !
            call geqrt3( norm, m-n1, n2, a( j1, j1 ), lda,                      &
                &                t( j1, j1 ), ldt, iinfo )
            !
            !        compute t3 = t(1:n1,j1:n) = -t1 y1^h y2 t2
            !
            do i=1,n1
                do j=1,n2
                    t( i, j+n1 ) = (a( j+n1, i ))
                end do
            end do
            !
            call htrmm( 'r', 'l', 'n', 'u', n1, n2, one,                   &
                &               a( j1, j1 ), lda, t( 1, j1 ), ldt )
            !
            call hgemm( 't', 'n', n1, n2, m-n, one, a( i1, 1 ), lda,       &
                &               a( i1, j1 ), lda, one, t( 1, j1 ), ldt )
            !
            call htrmm( 'l', 'u', 'n', 'n', n1, n2, -one, t, ldt,          &
                &               t( 1, j1 ), ldt )
            !
            call htrmm( 'r', 'u', 'n', 'n', n1, n2, one,                   &
                &               t( j1, j1 ), ldt, t( 1, j1 ), ldt )
            !
            !        y = (y1,y2); r = [ r1  a(1:n1,j1:n) ];  t = [t1 t3]
            !                         [  0        r2     ]       [ 0 t2]
            !
        end if
        !
        return
        !
        !     end of geqrt3
        !
    end subroutine

end submodule
