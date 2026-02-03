submodule (lpf_lapack_fp16) lpf_lapack_larfbv_fp16
contains

    !> \brief \b HLARFB applies a block reflector or its transpose to a general rectangular matrix.
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !        http://www.netlib.org/lapack/explore-html/
    !
    !> Download HLARFB + dependencies
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfb.f">
    !> [TGZ]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfb.f">
    !> [ZIP]</a>
    !> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfb.f">
    !> [TXT]</a>
    !
    !  Definition:
    !  ===========
    !
    !   SUBROUTINE HLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
    !                      T, LDT, C, LDC, WORK, LDWORK )
    !
    !   .. Scalar Arguments ..
    !   CHARACTER          DIRECT, SIDE, STOREV, TRANS
    !   INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
    !   ..
    !   .. Array Arguments ..
    !   REAL               C( LDC, * ), T( LDT, * ), V( LDV, * ),
    !  $                   WORK( LDWORK, * )
    !   ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> HLARFB applies a real block reflector H or its transpose H**T to a
    !> real m by n matrix C, from either the left or the right.
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
    !>          = 'C': Columnwise
    !>          = 'R': Rowwise
    !> \endverbatim
    !>
    !> \param[in] M
    !> \verbatim
    !>          M is INTEGER
    !>          The number of rows of the matrix C.
    !> \endverbatim
    !>
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The number of columns of the matrix C.
    !> \endverbatim
    !>
    !> \param[in] K
    !> \verbatim
    !>          K is INTEGER
    !>          The order of the matrix T (= the number of elementary
    !>          reflectors whose product defines the block reflector).
    !>          If SIDE = 'L', M >= K >= 0;
    !>          if SIDE = 'R', N >= K >= 0.
    !> \endverbatim
    !>
    !> \param[in] V
    !> \verbatim
    !>          V is REAL array, dimension
    !>                                (LDV,K) if STOREV = 'C'
    !>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
    !>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
    !>          The matrix V. See Further Details.
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
    !>          T is REAL array, dimension (LDT,K)
    !>          The triangular k by k matrix T in the representation of the
    !>          block reflector.
    !> \endverbatim
    !>
    !> \param[in] LDT
    !> \verbatim
    !>          LDT is INTEGER
    !>          The leading dimension of the array T. LDT >= K.
    !> \endverbatim
    !>
    !> \param[in,out] C
    !> \verbatim
    !>          C is REAL array, dimension (LDC,N)
    !>          On entry, the m by n matrix C.
    !>          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.
    !> \endverbatim
    !>
    !> \param[in] LDC
    !> \verbatim
    !>          LDC is INTEGER
    !>          The leading dimension of the array C. LDC >= max(1,M).
    !> \endverbatim
    !>
    !> \param[out] WORK
    !> \verbatim
    !>          WORK is REAL array, dimension (LDWORK,K)
    !> \endverbatim
    !>
    !> \param[in] LDWORK
    !> \verbatim
    !>          LDWORK is INTEGER
    !>          The leading dimension of the array WORK.
    !>          If SIDE = 'L', LDWORK >= max(1,N);
    !>          if SIDE = 'R', LDWORK >= max(1,M).
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
    !> \ingroup larfb
    !
    !> \par Further Details:
    !  =====================
    !>
    !> \verbatim
    !>
    !>  The shape of the matrix V and the storage of the vectors which define
    !>  the H(i) is best illustrated by the following example with n = 5 and
    !>  k = 3. The triangular part of V (including its diagonal) is not
    !>  referenced.
    !>
    !>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
    !>
    !>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
    !>                   ( v1  1    )                     (     1 v2 v2 v2 )
    !>                   ( v1 v2  1 )                     (        1 v3 v3 )
    !>                   ( v1 v2 v3 )
    !>                   ( v1 v2 v3 )
    !>
    !>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
    !>
    !>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
    !>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
    !>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
    !>                   (     1 v3 )
    !>                   (        1 )
    !> \endverbatim
    !>
    !  =====================================================================
    module subroutine larfbv( side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork ) bind(c, name="hlarfbv_")
        character, intent(in) :: direct, side, storev, trans
        integer(lpf_default_int_kind), intent(in) :: k, ldc, ldt, ldv, ldwork, m, n
        type(fp16), intent(in) :: t( ldt, * ), v( ldv, * )
        type(fp16), intent(inout) :: c( ldc, * ), work( ldwork, * )

        ! .. parameters ..
        type(fp16) :: one
        ! ..
        ! .. local scalars ..
        character :: transt
        integer(lpf_default_int_kind) :: i, j
        ! ..
        !
        ! quick return if possible
        !

        one = fp16(1.0)
        if( m.le.0 .or. n.le.0 ) &
            return
        !
        if( lsame( trans, 'n' ) ) then
            transt = 't'
        else
            transt = 'n'
        end if
        !
        if( lsame( storev, 'c' ) ) then
            !
            if( lsame( direct, 'f' ) ) then
                !
                !       let  v =  ( v1 )    (first k rows)
                !                 ( v2 )
                !       where  v1  is unit lower triangular.
                !
                if( lsame( side, 'l' ) ) then
                    !
                    !          form  h * c  or  h**t * c  where  c = ( c1 )
                    !                                                ( c2 )
                    !
                    !          w := c**t * v  =  (c1**t * v1 + c2**t * v2)  (stored in work)
                    !
                    !          w := c1**t
                    !
                    do j = 1, k
                        call hcopy( n, c( j, 1 ), ldc, work( 1, j ), 1_lpf_default_int_kind )
                    end do
                    !
                    !          w := w * v1
                    !
                    call htrmm( 'right', 'lower', 'no transpose', 'non-unit', n, k, one, v, ldv, work, ldwork )
                    if( m.gt.k ) then
                        !
                        !             w := w + c2**t * v2
                        !
                        call hgemm( 'transpose','no transpose', n, k, m-k, one, c( k+1, 1 ),ldc, v( k+1, 1 ),ldv, one, work,ldwork)
                    end if
                    !
                    !          w := w * t**t  or  w * t
                    !
                    call htrmm( 'right', 'upper', transt, 'non-unit', n, k, one, t, ldt, work, ldwork )
                    !
                    !          c := c - v * w**t
                    !
                    if( m.gt.k ) then
                        !
                        !             c2 := c2 - v2 * w**t
                        !
                        call hgemm( 'no transpose', 'transpose', m-k, n, k, -one, v( k+1, 1 ), ldv, work, ldwork, one, &
                            c( k+1, 1 ), ldc )
                    end if
                    !
                    !          w := w * v1**t
                    !
                    call htrmm( 'right', 'lower', 'transpose', 'non-unit', n, k, one, v, ldv, work, ldwork )
                    !
                    !          c1 := c1 - w**t
                    !
                    do j = 1, k
                        do i = 1, n
                            c( j, i ) = c( j, i ) - work( i, j )
                        end do
                    end do
                    !
                else if( lsame( side, 'r' ) ) then
                    !
                    !          form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    !
                    !          w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    !
                    !          w := c1
                    !
                    do j = 1, k
                        call hcopy( m, c( 1, j ), 1_lpf_default_int_kind, work( 1, j ), 1_lpf_default_int_kind )
                    end do
                    !
                    !          w := w * v1
                    !
                    call htrmm( 'right', 'lower', 'no transpose', 'non-unit', m, k, one, v, ldv, work, ldwork )
                    if( n.gt.k ) then
                        !
                        !             w := w + c2 * v2
                        !
                        call hgemm( 'no transpose', 'no transpose', m, k, n-k, one, c( 1, k+1 ), ldc, v( k+1, 1 ), ldv, &
                            one, work, ldwork )
                    end if
                    !
                    !          w := w * t  or  w * t**t
                    !
                    call htrmm( 'right', 'upper', trans, 'non-unit', m, k, one, t, ldt, work, ldwork )
                    !
                    !          c := c - w * v**t
                    !
                    if( n.gt.k ) then
                        !
                        !             c2 := c2 - w * v2**t
                        !
                        call hgemm( 'no transpose', 'transpose', m, n-k, k, -one, work, ldwork, v( k+1, 1 ), ldv, one, &
                            c( 1, k+1 ), ldc )
                    end if
                    !
                    !          w := w * v1**t
                    !
                    call htrmm( 'right', 'lower', 'transpose', 'non-unit', m, k, one, v, ldv, work, ldwork )
                    !
                    !          c1 := c1 - w
                    !
                    do j = 1, k
                        do i = 1, m
                            c( i, j ) = c( i, j ) - work( i, j )
                        end do
                    end do
                end if
                !
            else
                !
                !       let  v =  ( v1 )
                !                 ( v2 )    (last k rows)
                !       where  v2  is unit upper triangular.
                !
                if( lsame( side, 'l' ) ) then
                    !
                    !          form  h * c  or  h**t * c  where  c = ( c1 )
                    !                                                ( c2 )
                    !
                    !          w := c**t * v  =  (c1**t * v1 + c2**t * v2)  (stored in work)
                    !
                    !          w := c2**t
                    !
                    do j = 1, k
                        call hcopy( n, c( m-k+j, 1 ), ldc, work( 1, j ), 1_lpf_default_int_kind )
                    end do
                    !
                    !          w := w * v2
                    !
                    call htrmm( 'right', 'upper', 'no transpose', 'non-unit', n, k, one, v( m-k+1, 1 ), ldv, work, ldwork )
                    if( m.gt.k ) then
                        !
                        !             w := w + c1**t * v1
                        !
                        call hgemm( 'transpose', 'no transpose', n, k, m-k, one, c, ldc, v, ldv, one, work, ldwork )
                    end if
                    !
                    !          w := w * t**t  or  w * t
                    !
                    call htrmm( 'right', 'lower', transt, 'non-unit', n, k, one, t, ldt, work, ldwork )
                    !
                    !          c := c - v * w**t
                    !
                    if( m.gt.k ) then
                        !
                        !             c1 := c1 - v1 * w**t
                        !
                        call hgemm( 'no transpose', 'transpose', m-k, n, k, -one, v, ldv, work, ldwork, one, c, ldc )
                    end if
                    !
                    !          w := w * v2**t
                    !
                    call htrmm( 'right', 'upper', 'transpose', 'non-unit', n, k, one, v( m-k+1, 1 ), ldv, work, ldwork )
                    !
                    !          c2 := c2 - w**t
                    !
                    do j = 1, k
                        do i = 1, n
                            c( m-k+j, i ) = c( m-k+j, i ) - work( i, j )
                        end do
                    end do
                    !
                else if( lsame( side, 'r' ) ) then
                    !
                    !          form  c * h  or  c * h'  where  c = ( c1  c2 )
                    !
                    !          w := c * v  =  (c1*v1 + c2*v2)  (stored in work)
                    !
                    !          w := c2
                    !
                    do j = 1, k
                        call hcopy( m, c( 1, n-k+j ), 1_lpf_default_int_kind, work( 1, j ), 1_lpf_default_int_kind )
                    end do
                    !
                    !          w := w * v2
                    !
                    call htrmm( 'right', 'upper', 'no transpose', 'non-unit', m, k, one, v( n-k+1, 1 ), ldv, work, ldwork )
                    if( n.gt.k ) then
                        !
                        !             w := w + c1 * v1
                        !
                        call hgemm( 'no transpose', 'no transpose', m, k, n-k, one, c, ldc, v, ldv, one, work, ldwork )
                    end if
                    !
                    !          w := w * t  or  w * t**t
                    !
                    call htrmm( 'right', 'lower', trans, 'non-unit', m, k, one, t, ldt, work, ldwork )
                    !
                    !          c := c - w * v**t
                    !
                    if( n.gt.k ) then
                        !
                        !             c1 := c1 - w * v1**t
                        !
                        call hgemm( 'no transpose', 'transpose', m, n-k, k, -one, work, ldwork, v, ldv, one, c, ldc )
                    end if
                    !
                    !          w := w * v2**t
                    !
                    call htrmm( 'right', 'upper', 'transpose', 'non-unit', m, k, one, v( n-k+1, 1 ), ldv, work, ldwork )
                    !
                    !          c2 := c2 - w
                    !
                    do j = 1, k
                        do i = 1, m
                            c( i, n-k+j ) = c( i, n-k+j ) - work( i, j )
                        end do
                    end do
                end if
            end if
            !
        else if( lsame( storev, 'r' ) ) then
            !
            if( lsame( direct, 'f' ) ) then
                !
                !       let  v =  ( v1  v2 )    (v1: first k columns)
                !       where  v1  is unit upper triangular.
                !
                if( lsame( side, 'l' ) ) then
                    !
                    !          form  h * c  or  h**t * c  where  c = ( c1 )
                    !                                                ( c2 )
                    !
                    !          w := c**t * v**t  =  (c1**t * v1**t + c2**t * v2**t) (stored in work)
                    !
                    !          w := c1**t
                    !
                    do j = 1, k
                        call hcopy( n, c( j, 1 ), ldc, work( 1, j ), 1_lpf_default_int_kind )
                    end do
                    !
                    !          w := w * v1**t
                    !
                    call htrmm( 'right', 'upper', 'transpose', 'non-unit', n, k, one, v, ldv, work, ldwork )
                    if( m.gt.k ) then
                        !
                        !             w := w + c2**t * v2**t
                        !
                        call hgemm( 'transpose', 'transpose', n, k, m-k, one,  c( k+1, 1 ), ldc, v( 1, k+1 ), ldv, one, &
                            work, ldwork )
                    end if
                    !
                    !          w := w * t**t  or  w * t
                    !
                    call htrmm( 'right', 'upper', transt, 'non-unit', n, k, one, t, ldt, work, ldwork )
                    !
                    !          c := c - v**t * w**t
                    !
                    if( m.gt.k ) then
                        !
                        !             c2 := c2 - v2**t * w**t
                        !
                        call hgemm( 'transpose', 'transpose', m-k, n, k, -one,  v( 1, k+1 ), ldv, work, ldwork, one, &
                            c( k+1, 1 ), ldc )
                    end if
                    !
                    !          w := w * v1
                    !
                    call htrmm( 'right', 'upper', 'no transpose', 'non-unit', n, k, one, v, ldv, work, ldwork )
                    !
                    !          c1 := c1 - w**t
                    !
                    do j = 1, k
                        do i = 1, n
                            c( j, i ) = c( j, i ) - work( i, j )
                        end do
                    end do
                    !
                else if( lsame( side, 'r' ) ) then
                    !
                    !          form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    !
                    !          w := c * v**t  =  (c1*v1**t + c2*v2**t)  (stored in work)
                    !
                    !          w := c1
                    !
                    do j = 1, k
                        call hcopy( m, c( 1, j ), 1_lpf_default_int_kind, work( 1, j ), 1_lpf_default_int_kind )
                    end do
                    !
                    !          w := w * v1**t
                    !
                    call htrmm( 'right', 'upper', 'transpose', 'non-unit', m, k, one, v, ldv, work, ldwork )
                    if( n.gt.k ) then
                        !
                        !             w := w + c2 * v2**t
                        !
                        call hgemm( 'no transpose', 'transpose', m, k, n-k, one, c( 1, k+1 ), ldc, v( 1, k+1 ), ldv, &
                            one, work, ldwork )
                    end if
                    !
                    !          w := w * t  or  w * t**t
                    !
                    call htrmm( 'right', 'upper', trans, 'non-unit', m, k,  one, t, ldt, work, ldwork )
                    !
                    !          c := c - w * v
                    !
                    if( n.gt.k ) then
                        !
                        !             c2 := c2 - w * v2
                        !
                        call hgemm( 'no transpose', 'no transpose', m, n-k, k, -one, work, ldwork, v( 1, k+1 ), ldv, one, &
                            c( 1, k+1 ), ldc )
                    end if
                    !
                    !          w := w * v1
                    !
                    call htrmm( 'right', 'upper', 'no transpose', 'non-unit', m, k, one, v, ldv, work, ldwork )
                    !
                    !          c1 := c1 - w
                    !
                    do j = 1, k
                        do i = 1, m
                            c( i, j ) = c( i, j ) - work( i, j )
                        end do
                    end do
                    !
                end if
                !
            else
                !
                !       let  v =  ( v1  v2 )    (v2: last k columns)
                !       where  v2  is unit lower triangular.
                !
                if( lsame( side, 'l' ) ) then
                    !
                    !          form  h * c  or  h**t * c  where  c = ( c1 )
                    !                                                ( c2 )
                    !
                    !          w := c**t * v**t  =  (c1**t * v1**t + c2**t * v2**t) (stored in work)
                    !
                    !          w := c2**t
                    !
                    do j = 1, k
                        call hcopy( n, c( m-k+j, 1 ), ldc, work( 1, j ), 1_lpf_default_int_kind )
                    end do
                    !
                    !          w := w * v2**t
                    !
                    call htrmm( 'right', 'lower', 'transpose', 'non-unit', n, k, one, v( 1, m-k+1 ), ldv, work, ldwork )
                    if( m.gt.k ) then
                        !
                        !             w := w + c1**t * v1**t
                        !
                        call hgemm( 'transpose', 'transpose', n, k, m-k, one, c, ldc, v, ldv, one, work, ldwork )
                    end if
                    !
                    !          w := w * t**t  or  w * t
                    !
                    call htrmm( 'right', 'lower', transt, 'non-unit', n, k, one, t, ldt, work, ldwork )
                    !
                    !          c := c - v**t * w**t
                    !
                    if( m.gt.k ) then
                        !
                        !             c1 := c1 - v1**t * w**t
                        !
                        call hgemm( 'transpose', 'transpose', m-k, n, k, -one, v, ldv, work, ldwork, one, c, ldc )
                    end if
                    !
                    !          w := w * v2
                    !
                    call htrmm( 'right', 'lower', 'no transpose', 'non-unit', n, k, one, v( 1, m-k+1 ), ldv, work, ldwork )
                    !
                    !          c2 := c2 - w**t
                    !
                    do j = 1, k
                        do i = 1, n
                            c( m-k+j, i ) = c( m-k+j, i ) - work( i, j )
                        end do
                    end do
                    !
                else if( lsame( side, 'r' ) ) then
                    !
                    !          form  c * h  or  c * h**t  where  c = ( c1  c2 )
                    !
                    !          w := c * v**t  =  (c1*v1**t + c2*v2**t)  (stored in work)
                    !
                    !          w := c2
                    !
                    do j = 1, k
                        call hcopy( m, c( 1, n-k+j ), 1_lpf_default_int_kind, work( 1, j ), 1_lpf_default_int_kind )
                    end do
                    !
                    !          w := w * v2**t
                    !
                    call htrmm( 'right', 'lower', 'transpose', 'non-unit', m, k, one, v( 1, n-k+1 ), ldv, work, ldwork )
                    if( n.gt.k ) then
                        !
                        !             w := w + c1 * v1**t
                        !
                        call hgemm( 'no transpose', 'transpose', m, k, n-k, one, c, ldc, v, ldv, one, work, ldwork )
                    end if
                    !
                    !          w := w * t  or  w * t**t
                    !
                    call htrmm( 'right', 'lower', trans, 'non-unit', m, k, one, t, ldt, work, ldwork )
                    !
                    !          c := c - w * v
                    !
                    if( n.gt.k ) then
                        !
                        !             c1 := c1 - w * v1
                        !
                        call hgemm( 'no transpose', 'no transpose', m, n-k, k, -one, work, ldwork, v, ldv, one, c, ldc )
                    end if
                    !
                    !          w := w * v2
                    !
                    call htrmm( 'right', 'lower', 'no transpose', 'non-unit', m, k, one, v( 1, n-k+1 ), ldv, work, ldwork )
                    !
                    !          c1 := c1 - w
                    !
                    do j = 1, k
                        do i = 1, m
                            c( i, n-k+j ) = c( i, n-k+j ) - work( i, j )
                        end do
                    end do
                    !
                end if
                !
            end if
        end if
        !
        return
        !
        ! end of hlarfb
        !
    end subroutine larfbv
end submodule
