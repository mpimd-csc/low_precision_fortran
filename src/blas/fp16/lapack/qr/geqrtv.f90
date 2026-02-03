submodule(lpf_lapack_fp16) lpf_lapack_geqrtv_fp16

contains
    !> \brief \b HGEQRT
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> Download HGEQRT + dependencies
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
    !       SUBROUTINE HGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
    !
    !       .. Scalar Arguments ..
    !       INTEGER INFO, LDA, LDT, M, N, NB
    !       ..
    !       .. Array Arguments ..
    !       TYPE(FP16) A( LDA, * ), T( LDT, * ), WORK( * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> HGEQRT computes a blocked QR factorization of a type(fp16) M-by-N matrix A
    !> using the compact WY representation of Q.
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
    !>          The number of columns of the matrix A.  N >= 0.
    !> \endverbatim
    !>
    !> \param[in] NB
    !> \verbatim
    !>          NB is INTEGER
    !>          The block size to be used in the blocked QR.  MIN(M,N) >= NB
    !> \endverbatim
    !>
    !> \param[in,out] A
    !> \verbatim
    !>          A is TYPE(FP16) array, dimension (LDA,N)
    !>          On entry, the M-by-N matrix A.
    !>          On exit, the elements on and above the diagonal of the array
    !>          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
    !>          upper triangular if M >= N); the elements below the diagonal
    !>          are the columns of V.
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
    !>          T is TYPE(FP16) array, dimension (LDT,MIN(M,N))
    !>          The upper triangular block reflectors stored in compact form
    !>          as a sequence of upper triangular blocks.  See below
    !>          for further details.
    !> \endverbatim
    !>
    !> \param[in] LDT
    !> \verbatim
    !>          LDT is INTEGER
    !>          The leading dimension of the array T.  LDT >= NB.
    !> \endverbatim
    !>
    !> \param[out] WORK
    !> \verbatim
    !>          WORK is TYPE(FP16) array, dimension (NB*N)
    !> \endverbatim
    !>
    !> \param[out] INFO
    !> \verbatim
    !>          INFO is INTEGER
    !>          = 0:  successful exit
    !>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
    !> \ingroup geqrt
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
    !>
    !>  Let K=MIN(M,N).  The number of blocks is B = ceiling(K/NB), where ea
    !>  block is of order NB except for the last block, which is of order
    !>  IB = K - (B-1)*NB.  For each of the B blocks, a upper triangular blo
    !>  reflector factor is computed: T1, T2, ..., TB.  The NB-by-NB (and IB
    !>  for the last block) T's are stored in the NB-by-K matrix T as
    !>
    !>               T = (T1 T2 ... TB).
    !> \endverbatim
    !>
    !  =====================================================================
    module subroutine geqrtv( norm, m, n, nb, a, lda, diagr, t, ldt, work, info )
        integer(lpf_default_int_kind), intent(in)  :: lda, ldt, m, n, nb
        integer(lpf_default_int_kind), intent(out) :: info
        type(fp16), intent(inout) :: a( lda, * ), diagr(*), t( ldt, * ), work( * )
        character(len=*), intent(in)  :: norm

        !     .. local scalars ..
        integer(lpf_default_int_kind)    :: i, ib, iinfo, k
        logical, parameter    :: use_recursive_qr = .true.
        !     .. external subroutines ..
        !     ..
        !
        !     test the input arguments
        !

        info = 0
        if( m.lt.0 ) then
            info = -1
        else if( n.lt.0 ) then
            info = -2
        else if( nb.lt.1 .or. ( nb.gt.min(m,n) .and. min(m,n).gt.0 ) )then
            info = -3
        else if( lda.lt.max( 1, m ) ) then
            info = -5
        else if( ldt.lt.nb ) then
            info = -7
        end if
        if( info.ne.0 ) then
            call lpf_blas_xerbla( 'geqrtv', -info )
            return
        end if
        !
        !     quick return if possible
        !
        k = min( m, n )
        if( k.eq.0 ) return
        !
        !     blocked loop of length k
        !
        do i = 1, k,  nb
            ib = min( k-i+1, nb )
            !
            !     compute the qr factorization of the current block a(i:m,i:i+ib-1)
            !
            if( use_recursive_qr ) then
                call geqrt3v( norm, m-i+1, ib, a(i,i), lda, diagr(i), t(1,i), ldt,          &
                    &                    iinfo )
            else
                call geqrt2v( norm, m-i+1, ib, a(i,i), lda, diagr(i), t(1,i), ldt,          &
                    &                    iinfo )
            end if
            if( i+ib.le.n ) then
                !
                !     update by applying h**t to a(i:m,i+ib:n) from the left
                !
                call larfbv( 'l', 't', 'f', 'c', m-i+1, n-i-ib+1, ib,       &
                    &                   a( i, i ), lda, t( 1, i ), ldt,                &
                    &                   a( i, i+ib ), lda, work , n-i-ib+1 )
            end if
        end do
        return
        !
        !     end of geqrt
        !
    end subroutine

end submodule
