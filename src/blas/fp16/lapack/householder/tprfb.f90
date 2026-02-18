submodule(lpf_lapack_fp16) lpf_lapack_tprfb_fp16

contains



    !> \brief \b STPRFB applies a type(fp16) "triangular-pentagonal" block reflect
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> Download STPRFB + dependencies
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
    !       SUBROUTINE STPRFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, L,
    !                          V, LDV, T, LDT, A, LDA, B, LDB, WORK, LDWORK
    !
    !       .. Scalar Arguments ..
    !       CHARACTER DIRECT, SIDE, STOREV, TRANS
    !       INTEGER   K, L, LDA, LDB, LDT, LDV, LDWORK, M, N
    !       ..
    !       .. Array Arguments ..
    !       TYPE(FP16)   A( LDA, * ), B( LDB, * ), T( LDT, * ),
    !      $          V( LDV, * ), WORK( LDWORK, * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> STPRFB applies a type(fp16) "triangular-pentagonal" block reflector H or it
    !> transpose H**T to a type(fp16) matrix C, which is composed of two
    !> blocks A and B, either from the left or right.
    !>
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] SIDE
    !> \verbatim
    !>          SIDE is CHARACTER*1
    !>          = 'L': apply H or H**T from the Left
    !>          = 'R': apply H or H**T from the Right
    !> \endverbatim
    !>
    !> \param[in] TRANS
    !> \verbatim
    !>          TRANS is CHARACTER*1
    !>          = 'N': apply H (No transpose)
    !>          = 'T': apply H**T (Transpose)
    !> \endverbatim
    !>
    !> \param[in] DIRECT
    !> \verbatim
    !>          DIRECT is CHARACTER*1
    !>          Indicates how H is formed from a product of elementary
    !>          reflectors
    !>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
    !>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
    !> \endverbatim
    !>
    !> \param[in] STOREV
    !> \verbatim
    !>          STOREV is CHARACTER*1
    !>          Indicates how the vectors which define the elementary
    !>          reflectors are stored:
    !>          = 'C': Columns
    !>          = 'R': Rows
    !> \endverbatim
    !>
    !> \param[in] M
    !> \verbatim
    !>          M is INTEGER
    !>          The number of rows of the matrix B.
    !>          M >= 0.
    !> \endverbatim
    !>
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The number of columns of the matrix B.
    !>          N >= 0.
    !> \endverbatim
    !>
    !> \param[in] K
    !> \verbatim
    !>          K is INTEGER
    !>          The order of the matrix T, i.e. the number of elementary
    !>          reflectors whose product defines the block reflector.
    !>          K >= 0.
    !> \endverbatim
    !>
    !> \param[in] L
    !> \verbatim
    !>          L is INTEGER
    !>          The order of the trapezoidal part of V.
    !>          K >= L >= 0.  See Further Details.
    !> \endverbatim
    !>
    !> \param[in] V
    !> \verbatim
    !>          V is TYPE(FP16) array, dimension
    !>                                (LDV,K) if STOREV = 'C'
    !>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
    !>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
    !>          The pentagonal matrix V, which contains the elementary refle
    !>          H(1), H(2), ..., H(K).  See Further Details.
    !> \endverbatim
    !>
    !> \param[in] LDV
    !> \verbatim
    !>          LDV is INTEGER
    !>          The leading dimension of the array V.
    !>          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
    !>          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
    !>          if STOREV = 'R', LDV >= K.
    !> \endverbatim
    !>
    !> \param[in] T
    !> \verbatim
    !>          T is TYPE(FP16) array, dimension (LDT,K)
    !>          The triangular K-by-K matrix T in the representation of the
    !>          block reflector.
    !> \endverbatim
    !>
    !> \param[in] LDT
    !> \verbatim
    !>          LDT is INTEGER
    !>          The leading dimension of the array T.
    !>          LDT >= K.
    !> \endverbatim
    !>
    !> \param[in,out] A
    !> \verbatim
    !>          A is TYPE(FP16) array, dimension
    !>          (LDA,N) if SIDE = 'L' or (LDA,K) if SIDE = 'R'
    !>          On entry, the K-by-N or M-by-K matrix A.
    !>          On exit, A is overwritten by the corresponding block of
    !>          H*C or H**T*C or C*H or C*H**T.  See Further Details.
    !> \endverbatim
    !>
    !> \param[in] LDA
    !> \verbatim
    !>          LDA is INTEGER
    !>          The leading dimension of the array A.
    !>          If SIDE = 'L', LDA >= max(1,K);
    !>          If SIDE = 'R', LDA >= max(1,M).
    !> \endverbatim
    !>
    !> \param[in,out] B
    !> \verbatim
    !>          B is TYPE(FP16) array, dimension (LDB,N)
    !>          On entry, the M-by-N matrix B.
    !>          On exit, B is overwritten by the corresponding block of
    !>          H*C or H**T*C or C*H or C*H**T.  See Further Details.
    !> \endverbatim
    !>
    !> \param[in] LDB
    !> \verbatim
    !>          LDB is INTEGER
    !>          The leading dimension of the array B.
    !>          LDB >= max(1,M).
    !> \endverbatim
    !>
    !> \param[out] WORK
    !> \verbatim
    !>          WORK is TYPE(FP16) array, dimension
    !>          (LDWORK,N) if SIDE = 'L',
    !>          (LDWORK,K) if SIDE = 'R'.
    !> \endverbatim
    !>
    !> \param[in] LDWORK
    !> \verbatim
    !>          LDWORK is INTEGER
    !>          The leading dimension of the array WORK.
    !>          If SIDE = 'L', LDWORK >= K;
    !>          if SIDE = 'R', LDWORK >= M.
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
    !> \ingroup tprfb
    !
    !> \par Further Details:
    !  =====================
    !>
    !> \verbatim
    !>
    !>  The matrix C is a composite matrix formed from blocks A and B.
    !>  The block B is of size M-by-N; if SIDE = 'R', A is of size M-by-K,
    !>  and if SIDE = 'L', A is of size K-by-N.
    !>
    !>  If SIDE = 'R' and DIRECT = 'F', C = [A B].
    !>
    !>  If SIDE = 'L' and DIRECT = 'F', C = [A]
    !>                                      [B].
    !>
    !>  If SIDE = 'R' and DIRECT = 'B', C = [B A].
    !>
    !>  If SIDE = 'L' and DIRECT = 'B', C = [B]
    !>                                      [A].
    !>
    !>  The pentagonal matrix V is composed of a rectangular block V1 and a
    !>  trapezoidal block V2.  The size of the trapezoidal block is determin
    !>  the parameter L, where 0<=L<=K.  If L=K, the V2 block of V is triang
    !>  if L=0, there is no trapezoidal block, thus V = V1 is rectangular.
    !>
    !>  If DIRECT = 'F' and STOREV = 'C':  V = [V1]
    !>                                         [V2]
    !>     - V2 is upper trapezoidal (first L rows of K-by-K upper triangula
    !>
    !>  If DIRECT = 'F' and STOREV = 'R':  V = [V1 V2]
    !>
    !>     - V2 is lower trapezoidal (first L columns of K-by-K lower triang
    !>
    !>  If DIRECT = 'B' and STOREV = 'C':  V = [V2]
    !>                                         [V1]
    !>     - V2 is lower trapezoidal (last L rows of K-by-K lower triangular
    !>
    !>  If DIRECT = 'B' and STOREV = 'R':  V = [V2 V1]
    !>
    !>     - V2 is upper trapezoidal (last L columns of K-by-K upper triangu
    !>
    !>  If STOREV = 'C' and SIDE = 'L', V is M-by-K with V2 L-by-K.
    !>
    !>  If STOREV = 'C' and SIDE = 'R', V is N-by-K with V2 L-by-K.
    !>
    !>  If STOREV = 'R' and SIDE = 'L', V is K-by-M with V2 K-by-L.
    !>
    !>  If STOREV = 'R' and SIDE = 'R', V is K-by-N with V2 K-by-L.
    !> \endverbatim
    !>
    !  =====================================================================
    module subroutine tprfb( side, trans, direct, storev, m, n, k, l,       &
            &                   v, ldv, t, ldt, a, lda, b, ldb, work, ldwork )
        implicit none
        !
        !  -- lapack auxiliary routine --
        !  -- lapack is a software package provided by univ. of tennessee,    --
        !  -- univ. of california berkeley, univ. of colorado denver and nag ltd
        !
        !     .. scalar arguments ..
        character, intent(in) ::  direct, side, storev, trans
        integer(lpf_default_int_kind), intent(in) ::  k, l, lda, ldb, ldt, ldv, ldwork, m, n
        !     ..
        !     .. array arguments ..
        type(fp16), intent(inout)  ::  a( lda, * ), b( ldb, * ), t( ldt, * ), v( ldv, * ), work( ldwork, * )
        !     ..
        !
        !  =====================================================================
        !
        !     .. parameters ..
        type(fp16) ::  one, zero
        !     ..
        !     .. local scalars ..
        integer(lpf_default_int_kind) ::  i, j, mp, np, kp
        logical   left, forward, column, right, backward, row
        !     ..
        !     ..
        !     .. executable statements ..
        !
        !     quick return if possible
        !

        one = fp16(1.0)
        zero = fp16( 0.0 )

        if( m.le.0 .or. n.le.0 .or. k.le.0 .or. l.lt.0 ) return
        !
        if( lsame( storev, 'c' ) ) then
            column = .true.
            row = .false.
        else if ( lsame( storev, 'r' ) ) then
            column = .false.
            row = .true.
        else
            column = .false.
            row = .false.
        end if
        !
        if( lsame( side, 'l' ) ) then
            left = .true.
            right = .false.
        else if( lsame( side, 'r' ) ) then
            left = .false.
            right = .true.
        else
            left = .false.
            right = .false.
        end if
        !
        if( lsame( direct, 'f' ) ) then
            forward = .true.
            backward = .false.
        else if( lsame( direct, 'b' ) ) then
            forward = .false.
            backward = .true.
        else
            forward = .false.
            backward = .false.
        end if
        !
        ! ----------------------------------------------------------------------
        !
        if( column .and. forward .and. left  ) then
            !
            ! ----------------------------------------------------------------------
            !
            !        let  w =  [ i ]    (k-by-k)
            !                  [ v ]    (m-by-k)
            !
            !        form  h c  or  h**t c  where  c = [ a ]  (k-by-n)
            !                                          [ b ]  (m-by-n)
            !
            !        h = i - w t w**t          or  h**t = i - w t**t w**t
            !
            !        a = a -   t (a + v**t b)  or  a = a -   t**t (a + v**t b)
            !        b = b - v t (a + v**t b)  or  b = b - v t**t (a + v**t b)
            !
            ! ----------------------------------------------------------------------
            !
            mp = min( m-l+1, m )
            kp = min( l+1, k )
            !
            do j = 1, n
                do i = 1, l
                    work( i, j ) = b( m-l+i, j )
                end do
            end do
            call htrmm( 'l', 'u', 't', 'n', l, n, one, v( mp, 1 ), ldv,    &
                &               work, ldwork )
            call gemm( 't', 'n', l, n, m-l, one, v, ldv, b, ldb,          &
                &               one, work, ldwork )
            call gemm( 't', 'n', k-l, n, m, one, v( 1, kp ), ldv,         &
                &               b, ldb, zero, work( kp, 1 ), ldwork )
            !
            do j = 1, n
                do i = 1, k
                    work( i, j ) = work( i, j ) + a( i, j )
                end do
            end do
            !
            call htrmm( 'l', 'u', trans, 'n', k, n, one, t, ldt,           &
                &               work, ldwork )
            !
            do j = 1, n
                do i = 1, k
                    a( i, j ) = a( i, j ) - work( i, j )
                end do
            end do
            !
            call gemm( 'n', 'n', m-l, n, k, -one, v, ldv, work, ldwork,   &
                &               one, b, ldb )
            call gemm( 'n', 'n', l, n, k-l, -one, v( mp, kp ), ldv,       &
                &               work( kp, 1 ), ldwork, one, b( mp, 1 ),  ldb )
            call htrmm( 'l', 'u', 'n', 'n', l, n, one, v( mp, 1 ), ldv,    &
                &               work, ldwork )
            do j = 1, n
                do i = 1, l
                    b( m-l+i, j ) = b( m-l+i, j ) - work( i, j )
                end do
            end do
            !
            ! ----------------------------------------------------------------------
            !
        else if( column .and. forward .and. right ) then
            !
            ! ----------------------------------------------------------------------
            !
            !        let  w =  [ i ]    (k-by-k)
            !                  [ v ]    (n-by-k)
            !
            !        form  c h or  c h**t  where  c = [ a b ] (a is m-by-k, b is m-b
            !
            !        h = i - w t w**t          or  h**t = i - w t**t w**t
            !
            !        a = a - (a + b v) t       or  a = a - (a + b v) t**t
            !        b = b - (a + b v) t v**t  or  b = b - (a + b v) t**t v**t
            !
            ! ----------------------------------------------------------------------
            !
            np = min( n-l+1, n )
            kp = min( l+1, k )
            !
            do j = 1, l
                do i = 1, m
                    work( i, j ) = b( i, n-l+j )
                end do
            end do
            call htrmm( 'r', 'u', 'n', 'n', m, l, one, v( np, 1 ), ldv,    &
                &               work, ldwork )
            call gemm( 'n', 'n', m, l, n-l, one, b, ldb,                  &
                &               v, ldv, one, work, ldwork )
            call gemm( 'n', 'n', m, k-l, n, one, b, ldb,                  &
                &               v( 1, kp ), ldv, zero, work( 1, kp ), ldwork )
            !
            do j = 1, k
                do i = 1, m
                    work( i, j ) = work( i, j ) + a( i, j )
                end do
            end do
            !
            call htrmm( 'r', 'u', trans, 'n', m, k, one, t, ldt,           &
                &               work, ldwork )
            !
            do j = 1, k
                do i = 1, m
                    a( i, j ) = a( i, j ) - work( i, j )
                end do
            end do
            !
            call gemm( 'n', 't', m, n-l, k, -one, work, ldwork,           &
                &               v, ldv, one, b, ldb )
            call gemm( 'n', 't', m, l, k-l, -one, work( 1, kp ),          &
                &               ldwork,                                            &
                &               v( np, kp ), ldv, one, b( 1, np ), ldb )
            call htrmm( 'r', 'u', 't', 'n', m, l, one, v( np, 1 ), ldv,    &
                &               work, ldwork )
            do j = 1, l
                do i = 1, m
                    b( i, n-l+j ) = b( i, n-l+j ) - work( i, j )
                end do
            end do
            !
            ! ----------------------------------------------------------------------
            !
        else if( column .and. backward .and. left ) then
            !
            ! ----------------------------------------------------------------------
            !
            !        let  w =  [ v ]    (m-by-k)
            !                  [ i ]    (k-by-k)
            !
            !        form  h c  or  h**t c  where  c = [ b ]  (m-by-n)
            !                                          [ a ]  (k-by-n)
            !
            !        h = i - w t w**t         or  h**t = i - w t**t w**t
            !
            !        a = a -   t (a + v**t b)  or  a = a -   t**t (a + v**t b)
            !        b = b - v t (a + v**t b)  or  b = b - v t**t (a + v**t b)
            !
            ! ----------------------------------------------------------------------
            !
            mp = min( l+1, m )
            kp = min( k-l+1, k )
            !
            do j = 1, n
                do i = 1, l
                    work( k-l+i, j ) = b( i, j )
                end do
            end do
            !
            call htrmm( 'l', 'l', 't', 'n', l, n, one, v( 1, kp ), ldv,    &
                &               work( kp, 1 ), ldwork )
            call gemm( 't', 'n', l, n, m-l, one, v( mp, kp ), ldv,        &
                &               b( mp, 1 ), ldb, one, work( kp, 1 ), ldwork )
            call gemm( 't', 'n', k-l, n, m, one, v, ldv,                  &
                &               b, ldb, zero, work, ldwork )
            !
            do j = 1, n
                do i = 1, k
                    work( i, j ) = work( i, j ) + a( i, j )
                end do
            end do
            !
            call htrmm( 'l', 'l', trans, 'n', k, n, one, t, ldt,           &
                &               work, ldwork )
            !
            do j = 1, n
                do i = 1, k
                    a( i, j ) = a( i, j ) - work( i, j )
                end do
            end do
            !
            call gemm( 'n', 'n', m-l, n, k, -one, v( mp, 1 ), ldv,        &
                &               work, ldwork, one, b( mp, 1 ), ldb )
            call gemm( 'n', 'n', l, n, k-l, -one, v, ldv,                 &
                &               work, ldwork, one, b,  ldb )
            call htrmm( 'l', 'l', 'n', 'n', l, n, one, v( 1, kp ), ldv,    &
                &               work( kp, 1 ), ldwork )
            do j = 1, n
                do i = 1, l
                    b( i, j ) = b( i, j ) - work( k-l+i, j )
                end do
            end do
            !
            ! ----------------------------------------------------------------------
            !
        else if( column .and. backward .and. right ) then
            !
            ! ----------------------------------------------------------------------
            !
            !        let  w =  [ v ]    (n-by-k)
            !                  [ i ]    (k-by-k)
            !
            !        form  c h  or  c h**t  where  c = [ b a ] (b is m-by-n, a is m-
            !
            !        h = i - w t w**t          or  h**t = i - w t**t w**t
            !
            !        a = a - (a + b v) t       or  a = a - (a + b v) t**t
            !        b = b - (a + b v) t v**t  or  b = b - (a + b v) t**t v**t
            !
            ! ----------------------------------------------------------------------
            !
            np = min( l+1, n )
            kp = min( k-l+1, k )
            !
            do j = 1, l
                do i = 1, m
                    work( i, k-l+j ) = b( i, j )
                end do
            end do
            call htrmm( 'r', 'l', 'n', 'n', m, l, one, v( 1, kp ), ldv,    &
                &               work( 1, kp ), ldwork )
            call gemm( 'n', 'n', m, l, n-l, one, b( 1, np ), ldb,         &
                &               v( np, kp ), ldv, one, work( 1, kp ), ldwork )
            call gemm( 'n', 'n', m, k-l, n, one, b, ldb,                  &
                &               v, ldv, zero, work, ldwork )
            !
            do j = 1, k
                do i = 1, m
                    work( i, j ) = work( i, j ) + a( i, j )
                end do
            end do
            !
            call htrmm( 'r', 'l', trans, 'n', m, k, one, t, ldt,           &
                &               work, ldwork )
            !
            do j = 1, k
                do i = 1, m
                    a( i, j ) = a( i, j ) - work( i, j )
                end do
            end do
            !
            call gemm( 'n', 't', m, n-l, k, -one, work, ldwork,           &
                &               v( np, 1 ), ldv, one, b( 1, np ), ldb )
            call gemm( 'n', 't', m, l, k-l, -one, work, ldwork,           &
                &               v, ldv, one, b, ldb )
            call htrmm( 'r', 'l', 't', 'n', m, l, one, v( 1, kp ), ldv,    &
                &               work( 1, kp ), ldwork )
            do j = 1, l
                do i = 1, m
                    b( i, j ) = b( i, j ) - work( i, k-l+j )
                end do
            end do
            !
            ! ----------------------------------------------------------------------
            !
        else if( row .and. forward .and. left ) then
            !
            ! ----------------------------------------------------------------------
            !
            !        let  w =  [ i v ] ( i is k-by-k, v is k-by-m )
            !
            !        form  h c  or  h**t c  where  c = [ a ]  (k-by-n)
            !                                          [ b ]  (m-by-n)
            !
            !        h = i - w**t t w          or  h**t = i - w**t t**t w
            !
            !        a = a -      t (a + v b)  or  a = a -      t**t (a + v b)
            !        b = b - v**t t (a + v b)  or  b = b - v**t t**t (a + v b)
            !
            ! ----------------------------------------------------------------------
            !
            mp = min( m-l+1, m )
            kp = min( l+1, k )
            !
            do j = 1, n
                do i = 1, l
                    work( i, j ) = b( m-l+i, j )
                end do
            end do
            call htrmm( 'l', 'l', 'n', 'n', l, n, one, v( 1, mp ), ldv,    &
                &               work, ldb )
            call gemm( 'n', 'n', l, n, m-l, one, v, ldv,b, ldb,           &
                &               one, work, ldwork )
            call gemm( 'n', 'n', k-l, n, m, one, v( kp, 1 ), ldv,         &
                &               b, ldb, zero, work( kp, 1 ), ldwork )
            !
            do j = 1, n
                do i = 1, k
                    work( i, j ) = work( i, j ) + a( i, j )
                end do
            end do
            !
            call htrmm( 'l', 'u', trans, 'n', k, n, one, t, ldt,           &
                &               work, ldwork )
            !
            do j = 1, n
                do i = 1, k
                    a( i, j ) = a( i, j ) - work( i, j )
                end do
            end do
            !
            call gemm( 't', 'n', m-l, n, k, -one, v, ldv, work, ldwork,   &
                &               one, b, ldb )
            call gemm( 't', 'n', l, n, k-l, -one, v( kp, mp ), ldv,       &
                &               work( kp, 1 ), ldwork, one, b( mp, 1 ), ldb )
            call htrmm( 'l', 'l', 't', 'n', l, n, one, v( 1, mp ), ldv,    &
                &               work, ldwork )
            do j = 1, n
                do i = 1, l
                    b( m-l+i, j ) = b( m-l+i, j ) - work( i, j )
                end do
            end do
            !
            ! ----------------------------------------------------------------------
            !
        else if( row .and. forward .and. right ) then
            !
            ! ----------------------------------------------------------------------
            !
            !        let  w =  [ i v ] ( i is k-by-k, v is k-by-n )
            !
            !        form  c h  or  c h**t  where  c = [ a b ] (a is m-by-k, b is m-
            !
            !        h = i - w**t t w            or  h**t = i - w**t t**t w
            !
            !        a = a - (a + b v**t) t      or  a = a - (a + b v**t) t**t
            !        b = b - (a + b v**t) t v    or  b = b - (a + b v**t) t**t v
            !
            ! ----------------------------------------------------------------------
            !
            np = min( n-l+1, n )
            kp = min( l+1, k )
            !
            do j = 1, l
                do i = 1, m
                    work( i, j ) = b( i, n-l+j )
                end do
            end do
            call htrmm( 'r', 'l', 't', 'n', m, l, one, v( 1, np ), ldv,    &
                &               work, ldwork )
            call gemm( 'n', 't', m, l, n-l, one, b, ldb, v, ldv,          &
                &               one, work, ldwork )
            call gemm( 'n', 't', m, k-l, n, one, b, ldb,                  &
                &               v( kp, 1 ), ldv, zero, work( 1, kp ), ldwork )
            !
            do j = 1, k
                do i = 1, m
                    work( i, j ) = work( i, j ) + a( i, j )
                end do
            end do
            !
            call htrmm( 'r', 'u', trans, 'n', m, k, one, t, ldt,           &
                &               work, ldwork )
            !
            do j = 1, k
                do i = 1, m
                    a( i, j ) = a( i, j ) - work( i, j )
                end do
            end do
            !
            call gemm( 'n', 'n', m, n-l, k, -one, work, ldwork,           &
                &               v, ldv, one, b, ldb )
            call gemm( 'n', 'n', m, l, k-l, -one, work( 1, kp ),          &
                &               ldwork,                                            &
                &               v( kp, np ), ldv, one, b( 1, np ), ldb )
            call htrmm( 'r', 'l', 'n', 'n', m, l, one, v( 1, np ), ldv,    &
                &               work, ldwork )
            do j = 1, l
                do i = 1, m
                    b( i, n-l+j ) = b( i, n-l+j ) - work( i, j )
                end do
            end do
            !
            ! ----------------------------------------------------------------------
            !
        else if( row .and. backward .and. left ) then
            !
            ! ----------------------------------------------------------------------
            !
            !        let  w =  [ v i ] ( i is k-by-k, v is k-by-m )
            !
            !        form  h c  or  h**t c  where  c = [ b ]  (m-by-n)
            !                                          [ a ]  (k-by-n)
            !
            !        h = i - w**t t w          or  h**t = i - w**t t**t w
            !
            !        a = a -      t (a + v b)  or  a = a -      t**t (a + v b)
            !        b = b - v**t t (a + v b)  or  b = b - v**t t**t (a + v b)
            !
            ! ----------------------------------------------------------------------
            !
            mp = min( l+1, m )
            kp = min( k-l+1, k )
            !
            do j = 1, n
                do i = 1, l
                    work( k-l+i, j ) = b( i, j )
                end do
            end do
            call htrmm( 'l', 'u', 'n', 'n', l, n, one, v( kp, 1 ), ldv,    &
                &               work( kp, 1 ), ldwork )
            call gemm( 'n', 'n', l, n, m-l, one, v( kp, mp ), ldv,        &
                &               b( mp, 1 ), ldb, one, work( kp, 1 ), ldwork )
            call gemm( 'n', 'n', k-l, n, m, one, v, ldv, b, ldb,          &
                &               zero, work, ldwork )
            !
            do j = 1, n
                do i = 1, k
                    work( i, j ) = work( i, j ) + a( i, j )
                end do
            end do
            !
            call htrmm( 'l', 'l ', trans, 'n', k, n, one, t, ldt,          &
                &               work, ldwork )
            !
            do j = 1, n
                do i = 1, k
                    a( i, j ) = a( i, j ) - work( i, j )
                end do
            end do
            !
            call gemm( 't', 'n', m-l, n, k, -one, v( 1, mp ), ldv,        &
                &               work, ldwork, one, b( mp, 1 ), ldb )
            call gemm( 't', 'n', l, n, k-l, -one, v, ldv,                 &
                &               work, ldwork, one, b, ldb )
            call htrmm( 'l', 'u', 't', 'n', l, n, one, v( kp, 1 ), ldv,    &
                &               work( kp, 1 ), ldwork )
            do j = 1, n
                do i = 1, l
                    b( i, j ) = b( i, j ) - work( k-l+i, j )
                end do
            end do
            !
            ! ----------------------------------------------------------------------
            !
        else if( row .and. backward .and. right ) then
            !
            ! ----------------------------------------------------------------------
            !
            !        let  w =  [ v i ] ( i is k-by-k, v is k-by-n )
            !
            !        form  c h  or  c h**t  where  c = [ b a ] (a is m-by-k, b is m-
            !
            !        h = i - w**t t w            or  h**t = i - w**t t**t w
            !
            !        a = a - (a + b v**t) t      or  a = a - (a + b v**t) t**t
            !        b = b - (a + b v**t) t v    or  b = b - (a + b v**t) t**t v
            !
            ! ----------------------------------------------------------------------
            !
            np = min( l+1, n )
            kp = min( k-l+1, k )
            !
            do j = 1, l
                do i = 1, m
                    work( i, k-l+j ) = b( i, j )
                end do
            end do
            call htrmm( 'r', 'u', 't', 'n', m, l, one, v( kp, 1 ), ldv,    &
                &               work( 1, kp ), ldwork )
            call gemm( 'n', 't', m, l, n-l, one, b( 1, np ), ldb,         &
                &               v( kp, np ), ldv, one, work( 1, kp ), ldwork )
            call gemm( 'n', 't', m, k-l, n, one, b, ldb, v, ldv,          &
                &               zero, work, ldwork )
            !
            do j = 1, k
                do i = 1, m
                    work( i, j ) = work( i, j ) + a( i, j )
                end do
            end do
            !
            call htrmm( 'r', 'l', trans, 'n', m, k, one, t, ldt,           &
                &               work, ldwork )
            !
            do j = 1, k
                do i = 1, m
                    a( i, j ) = a( i, j ) - work( i, j )
                end do
            end do
            !
            call gemm( 'n', 'n', m, n-l, k, -one, work, ldwork,           &
                &               v( 1, np ), ldv, one, b( 1, np ), ldb )
            call gemm( 'n', 'n', m, l, k-l , -one, work, ldwork,          &
                &               v, ldv, one, b, ldb )
            call htrmm( 'r', 'u', 'n', 'n', m, l, one, v( kp, 1 ), ldv,    &
                &               work( 1, kp ), ldwork )
            do j = 1, l
                do i = 1, m
                    b( i, j ) = b( i, j ) - work( i, k-l+j )
                end do
            end do
            !
        end if
        !
        return
        !
        !     end of stprfb
        !
    end subroutine
end submodule
