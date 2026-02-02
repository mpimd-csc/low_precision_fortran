submodule (lpf_lapack_fp16) lpf_lapack_larft_fp16
contains
    !> \brief \b HLARFT forms the triangular factor T of a block reflector H
    !
    !  =========== DOCUMENTATION ===========
    !
    ! Online html documentation available at
    !            http://www.netlib.org/lapack/explore-html/
    !
    !> Download HLARFT + dependencies
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
    !       RECURSIVE SUBROUTINE HLARFT( DIRECT, STOREV, N, K, V, LDV, TAU,
    !
    !       .. Scalar Arguments ..
    !       CHARACTER          DIRECT, STOREV
    !       INTEGER            K, LDT, LDV, N
    !       ..
    !       .. Array Arguments ..
    !       REAL               T( LDT, * ), TAU( * ), V( LDV, * )
    !       ..
    !
    !
    !> \par Purpose:
    !  =============
    !>
    !> \verbatim
    !>
    !> HLARFT forms the triangular factor T of a real block reflector H
    !> of order n, which is defined as a product of k elementary reflectors.
    !>
    !> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
    !>
    !> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
    !>
    !> If STOREV = 'C', the vector which defines the elementary reflector
    !> H(i) is stored in the i-th column of the array V, and
    !>
    !>    H  =  I - V * T * V**T
    !>
    !> If STOREV = 'R', the vector which defines the elementary reflector
    !> H(i) is stored in the i-th row of the array V, and
    !>
    !>    H  =  I - V**T * T * V
    !> \endverbatim
    !
    !  Arguments:
    !  ==========
    !
    !> \param[in] DIRECT
    !> \verbatim
    !>          DIRECT is CHARACTER*1
    !>          Specifies the order in which the elementary reflectors are
    !>          multiplied to form the block reflector:
    !>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
    !>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
    !> \endverbatim
    !>
    !> \param[in] STOREV
    !> \verbatim
    !>          STOREV is CHARACTER*1
    !>          Specifies how the vectors which define the elementary
    !>          reflectors are stored (see also Further Details):
    !>          = 'C': columnwise
    !>          = 'R': rowwise
    !> \endverbatim
    !>
    !> \param[in] N
    !> \verbatim
    !>          N is INTEGER
    !>          The order of the block reflector H. N >= 0.
    !> \endverbatim
    !>
    !> \param[in] K
    !> \verbatim
    !>          K is INTEGER
    !>          The order of the triangular factor T (= the number of
    !>          elementary reflectors). K >= 1.
    !> \endverbatim
    !>
    !> \param[in] V
    !> \verbatim
    !>          V is REAL array, dimension
    !>                               (LDV,K) if STOREV = 'C'
    !>                               (LDV,N) if STOREV = 'R'
    !>          The matrix V. See further details.
    !> \endverbatim
    !>
    !> \param[in] LDV
    !> \verbatim
    !>          LDV is INTEGER
    !>          The leading dimension of the array V.
    !>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
    !> \endverbatim
    !>
    !> \param[in] TAU
    !> \verbatim
    !>          TAU is REAL array, dimension (K)
    !>          TAU(i) must contain the scalar factor of the elementary
    !>          reflector H(i).
    !> \endverbatim
    !>
    !> \param[out] T
    !> \verbatim
    !>          T is REAL array, dimension (LDT,K)
    !>          The k by k triangular factor T of the block reflector.
    !>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T i
    !>          lower triangular. The rest of the array is not used.
    !> \endverbatim
    !>
    !> \param[in] LDT
    !> \verbatim
    !>          LDT is INTEGER
    !>          The leading dimension of the array T. LDT >= K.
    !> \endverbatim
    !
    !  Authors:
    !  ========
    !
    !> \author Univ. of Tennessee
    !> \author Univ. of California Berkeley
    !> \author Johnathan Rhyne, Univ. of Colorado Denver (original author, 2
    !> \author NAG Ltd.
    !
    !> \ingroup larft
    !
    !> \par Further Details:
    !  =====================
    !>
    !> \verbatim
    !>
    !>  The shape of the matrix V and the storage of the vectors which defin
    !>  the H(i) is best illustrated by the following example with n = 5 and
    !>  k = 3. The elements equal to 1 are not stored.
    !>
    !>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R'
    !>
    !>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
    !>                   ( v1  1    )                     (     1 v2 v2 v2 )
    !>                   ( v1 v2  1 )                     (        1 v3 v3 )
    !>                   ( v1 v2 v3 )
    !>                   ( v1 v2 v3 )
    !>
    !>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R'
    !>
    !>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
    !>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
    !>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
    !>                   (     1 v3 )
    !>                   (        1 )
    !> \endverbatim
    !>
    !  =====================================================================
    module recursive subroutine larft( direct, storev, n, k, v, ldv, tau, t, ldt ) bind(c, name="hlarft_")
        character, intent(in) ::         direct, storev
        integer(lpf_default_int_kind), intent(in)   ::     k, ldt, ldv, n
        !     ..
        !     .. array arguments ..
        !
        type(fp16), intent(out) ::   t( ldt, * )
        type(fp16), intent(in)  ::   tau( * ), v( ldv, * )
        !     ..
        !
        !     .. parameters ..
        !
        type(fp16)               :: one, neg_one, zero
        !
        !     .. local scalars ..
        !
        integer(lpf_default_int_kind)            :: i,j,l
        logical            :: qr,lq,ql,dirf,colv
        !
        !     the general scheme used is inspired by the approach inside dgeqrt3
        !     which was (at the time of writing this code):
        !     based on the algorithm of elmroth and gustavson,
        !     ibm j. res. develop. vol 44 no. 4 july 2000.
        !     ..
        !     .. executable statements ..
        !
        !     quick return if possible
        !

        one=1.0
        zero = 0.0
        neg_one=-1.0

        if(n.eq.0.or.k.eq.0) then
            return
        end if
        !
        !     base case
        !
        if(n.eq.1.or.k.eq.1) then
            t(1,1) = tau(1)
            return
        end if
        !
        !     beginning of executable statements
        !
        l = k / 2
        !
        !     determine what kind of q we need to compute
        !     we assume that if the user doesn't provide 'f' for direct,
        !     then they meant to provide 'b' and if they don't provide
        !     'c' for storev, then they meant to provide 'r'
        !
        dirf = lsame(direct,'f')
        colv = lsame(storev,'c')
        !
        !     qr happens when we have forward direction in column storage
        !
        qr = dirf.and.colv
        !
        !     lq happens when we have forward direction in row storage
        !
        lq = dirf.and.(.not.colv)
        !
        !     ql happens when we have backward direction in column storage
        !
        ql = (.not.dirf).and.colv
        !
        !     the last case is rq. due to how we structured this, if the
        !     above 3 are false, then rq must be true, so we never store
        !     this
        !     rq happens when we have backward direction in row storage
        !     rq = (.not.dirf).and.(.not.colv)
        !
        if(qr) then
            !
            !        break v apart into 6 components
            !
            !        v = |---------------|
            !            |v_{1,1} 0      |
            !            |v_{2,1} v_{2,2}|
            !            |v_{3,1} v_{3,2}|
            !            |---------------|
            !
            !        v_{1,1}\in\r^{l,l}      unit lower triangular
            !        v_{2,1}\in\r^{k-l,l}    rectangular
            !        v_{3,1}\in\r^{n-k,l}    rectangular
            !
            !        v_{2,2}\in\r^{k-l,k-l}  unit lower triangular
            !        v_{3,2}\in\r^{n-k,k-l}  rectangular
            !
            !        we will construct the t matrix
            !        t = |---------------|
            !            |t_{1,1} t_{1,2}|
            !            |0       t_{2,2}|
            !            |---------------|
            !
            !        t is the triangular factor obtained from block reflectors.
            !        to motivate the structure, assume we have already computed t_{1
            !        and t_{2,2}. then collect the associated reflectors in v_1 and
            !
            !        t_{1,1}\in\r^{l, l}     upper triangular
            !        t_{2,2}\in\r^{k-l, k-l} upper triangular
            !        t_{1,2}\in\r^{l, k-l}   rectangular
            !
            !        where l = floor(k/2)
            !
            !        then, consider the product:
            !
            !        (i - v_1*t_{1,1}*v_1')*(i - v_2*t_{2,2}*v_2')
            !        = i - v_1*t_{1,1}*v_1' - v_2*t_{2,2}*v_2' + v_1*t_{1,1}*v_1'*v_
            !
            !        define t_{1,2} = -t_{1,1}*v_1'*v_2*t_{2,2}
            !
            !        then, we can define the matrix v as
            !        v = |-------|
            !            |v_1 v_2|
            !            |-------|
            !
            !        so, our product is equivalent to the matrix product
            !        i - v*t*v'
            !        this means, we can compute t_{1,1} and t_{2,2}, then use this i
            !        to compute t_{1,2}
            !
            !        compute t_{1,1} recursively
            !
            call larft(direct, storev, n, l, v, ldv, tau, t, ldt)
            !
            !        compute t_{2,2} recursively
            !
            call larft(direct, storev, n-l, k-l, v(l+1, l+1), ldv, tau(l+1), t(l+1, l+1), ldt)
            !
            !        compute t_{1,2}
            !        t_{1,2} = v_{2,1}'
            !
            do j = 1, l
                do i = 1, k-l
                    t(j, l+i) = v(l+i, j)
                end do
            end do
            !
            !        t_{1,2} = t_{1,2}*v_{2,2}
            !
            call htrmm('right', 'lower', 'no transpose', 'unit', l, k-l, one, v(l+1, l+1), ldv, t(1, l+1), ldt)

            !
            !        t_{1,2} = v_{3,1}'*v_{3,2} + t_{1,2}
            !        note: we assume k <= n, and gemm will do nothing if n=k
            !
            call hgemm('transpose', 'no transpose', l, k-l, n-k, one, v(k+1, 1), ldv, v(k+1, l+1), ldv, one, t(1, l+1), ldt)
            !
            !        at this point, we have that t_{1,2} = v_1'*v_2
            !        all that is left is to pre and post multiply by -t_{1,1} and t_
            !        respectively.
            !
            !        t_{1,2} = -t_{1,1}*t_{1,2}
            !
            call htrmm('left', 'upper', 'no transpose', 'non-unit', l, k-l, neg_one, t, ldt, t(1, l+1), ldt)
            !
            !        t_{1,2} = t_{1,2}*t_{2,2}
            !
            call htrmm('right', 'upper', 'no transpose', 'non-unit', l, k-l, one, t(l+1, l+1), ldt, t(1, l+1), ldt)

        else if(lq) then
            !
            !        break v apart into 6 components
            !
            !        v = |----------------------|
            !            |v_{1,1} v_{1,2} v{1,3}|
            !            |0       v_{2,2} v{2,3}|
            !            |----------------------|
            !
            !        v_{1,1}\in\r^{l,l}      unit upper triangular
            !        v_{1,2}\in\r^{l,k-l}    rectangular
            !        v_{1,3}\in\r^{l,n-k}    rectangular
            !
            !        v_{2,2}\in\r^{k-l,k-l}  unit upper triangular
            !        v_{2,3}\in\r^{k-l,n-k}  rectangular
            !
            !        where l = floor(k/2)
            !
            !        we will construct the t matrix
            !        t = |---------------|
            !            |t_{1,1} t_{1,2}|
            !            |0       t_{2,2}|
            !            |---------------|
            !
            !        t is the triangular factor obtained from block reflectors.
            !        to motivate the structure, assume we have already computed t_{1
            !        and t_{2,2}. then collect the associated reflectors in v_1 and
            !
            !        t_{1,1}\in\r^{l, l}     upper triangular
            !        t_{2,2}\in\r^{k-l, k-l} upper triangular
            !        t_{1,2}\in\r^{l, k-l}   rectangular
            !
            !        then, consider the product:
            !
            !        (i - v_1'*t_{1,1}*v_1)*(i - v_2'*t_{2,2}*v_2)
            !        = i - v_1'*t_{1,1}*v_1 - v_2'*t_{2,2}*v_2 + v_1'*t_{1,1}*v_1*v_
            !
            !        define t_{1,2} = -t_{1,1}*v_1*v_2'*t_{2,2}
            !
            !        then, we can define the matrix v as
            !        v = |---|
            !            |v_1|
            !            |v_2|
            !            |---|
            !
            !        so, our product is equivalent to the matrix product
            !        i - v'*t*v
            !        this means, we can compute t_{1,1} and t_{2,2}, then use this i
            !        to compute t_{1,2}
            !
            !        compute t_{1,1} recursively
            !
            call larft(direct, storev, n, l, v, ldv, tau, t, ldt)
            !
            !        compute t_{2,2} recursively
            !
            call larft(direct, storev, n-l, k-l, v(l+1, l+1), ldv, tau(l+1), t(l+1, l+1), ldt)

            !
            !        compute t_{1,2}
            !        t_{1,2} = v_{1,2}
            !
            call hlacpy('all', l, k-l, v(1, l+1), ldv, t(1, l+1), ldt)
            !
            !        t_{1,2} = t_{1,2}*v_{2,2}'
            !
            call htrmm('right', 'upper', 'transpose', 'unit', l, k-l, one, v(l+1, l+1), ldv, t(1, l+1), ldt)

            !
            !        t_{1,2} = v_{1,3}*v_{2,3}' + t_{1,2}
            !        note: we assume k <= n, and gemm will do nothing if n=k
            !
            call hgemm('no transpose', 'transpose', l, k-l, n-k, one, v(1, k+1), ldv, v(l+1, k+1), ldv, one, t(1, l+1), ldt)
            !
            !        at this point, we have that t_{1,2} = v_1*v_2'
            !        all that is left is to pre and post multiply by -t_{1,1} and t_
            !        respectively.
            !
            !        t_{1,2} = -t_{1,1}*t_{1,2}
            !
            call htrmm('left', 'upper', 'no transpose', 'non-unit', l, k-l, neg_one, t, ldt, t(1, l+1), ldt)

            !
            !        t_{1,2} = t_{1,2}*t_{2,2}
            !
            call htrmm('right', 'upper', 'no transpose', 'non-unit', l, k-l, one, t(l+1, l+1), ldt, t(1, l+1), ldt)
        else if(ql) then
            !
            !        break v apart into 6 components
            !
            !        v = |---------------|
            !            |v_{1,1} v_{1,2}|
            !            |v_{2,1} v_{2,2}|
            !            |0       v_{3,2}|
            !            |---------------|
            !
            !        v_{1,1}\in\r^{n-k,k-l}  rectangular
            !        v_{2,1}\in\r^{k-l,k-l}  unit upper triangular
            !
            !        v_{1,2}\in\r^{n-k,l}    rectangular
            !        v_{2,2}\in\r^{k-l,l}    rectangular
            !        v_{3,2}\in\r^{l,l}      unit upper triangular
            !
            !        we will construct the t matrix
            !        t = |---------------|
            !            |t_{1,1} 0      |
            !            |t_{2,1} t_{2,2}|
            !            |---------------|
            !
            !        t is the triangular factor obtained from block reflectors.
            !        to motivate the structure, assume we have already computed t_{1
            !        and t_{2,2}. then collect the associated reflectors in v_1 and
            !
            !        t_{1,1}\in\r^{k-l, k-l} non-unit lower triangular
            !        t_{2,2}\in\r^{l, l}     non-unit lower triangular
            !        t_{2,1}\in\r^{k-l, l}   rectangular
            !
            !        where l = floor(k/2)
            !
            !        then, consider the product:
            !
            !        (i - v_2*t_{2,2}*v_2')*(i - v_1*t_{1,1}*v_1')
            !        = i - v_2*t_{2,2}*v_2' - v_1*t_{1,1}*v_1' + v_2*t_{2,2}*v_2'*v_
            !
            !        define t_{2,1} = -t_{2,2}*v_2'*v_1*t_{1,1}
            !
            !        then, we can define the matrix v as
            !        v = |-------|
            !            |v_1 v_2|
            !            |-------|
            !
            !        so, our product is equivalent to the matrix product
            !        i - v*t*v'
            !        this means, we can compute t_{1,1} and t_{2,2}, then use this i
            !        to compute t_{2,1}
            !
            !        compute t_{1,1} recursively
            !
            call larft(direct, storev, n-l, k-l, v, ldv, tau, t, ldt)
            !
            !        compute t_{2,2} recursively
            !
            call larft(direct, storev, n, l, v(1, k-l+1), ldv, tau(k-l+1), t(k-l+1, k-l+1), ldt)
            !
            !        compute t_{2,1}
            !        t_{2,1} = v_{2,2}'
            !
            do j = 1, k-l
                do i = 1, l
                    t(k-l+i, j) = v(n-k+j, k-l+i)
                end do
            end do
            !
            !        t_{2,1} = t_{2,1}*v_{2,1}
            !
            call htrmm('right', 'upper', 'no transpose', 'unit', l, k-l, one, v(n-k+1, 1), ldv, t(k-l+1, 1), ldt)

            !
            !        t_{2,1} = v_{2,2}'*v_{2,1} + t_{2,1}
            !        note: we assume k <= n, and gemm will do nothing if n=k
            !
            call hgemm('transpose', 'no transpose', l, k-l, n-k, one, v(1, k-l+1), ldv, v, ldv, one, t(k-l+1, 1), ldt)
            !
            !        at this point, we have that t_{2,1} = v_2'*v_1
            !        all that is left is to pre and post multiply by -t_{2,2} and t_
            !        respectively.
            !
            !        t_{2,1} = -t_{2,2}*t_{2,1}
            !
            call htrmm('left', 'lower', 'no transpose', 'non-unit', l,     &
                &               k-l, neg_one, t(k-l+1, k-l+1), ldt,                &
                &               t(k-l+1, 1), ldt)
            !
            !        t_{2,1} = t_{2,1}*t_{1,1}
            !
            call htrmm('right', 'lower', 'no transpose', 'non-unit', l,    &
                &               k-l, one, t, ldt, t(k-l+1, 1), ldt)
        else
            !
            !        else means rq case
            !
            !        break v apart into 6 components
            !
            !        v = |-----------------------|
            !            |v_{1,1} v_{1,2} 0      |
            !            |v_{2,1} v_{2,2} v_{2,3}|
            !            |-----------------------|
            !
            !        v_{1,1}\in\r^{k-l,n-k}  rectangular
            !        v_{1,2}\in\r^{k-l,k-l}  unit lower triangular
            !
            !        v_{2,1}\in\r^{l,n-k}    rectangular
            !        v_{2,2}\in\r^{l,k-l}    rectangular
            !        v_{2,3}\in\r^{l,l}      unit lower triangular
            !
            !        we will construct the t matrix
            !        t = |---------------|
            !            |t_{1,1} 0      |
            !            |t_{2,1} t_{2,2}|
            !            |---------------|
            !
            !        t is the triangular factor obtained from block reflectors.
            !        to motivate the structure, assume we have already computed t_{1
            !        and t_{2,2}. then collect the associated reflectors in v_1 and
            !
            !        t_{1,1}\in\r^{k-l, k-l} non-unit lower triangular
            !        t_{2,2}\in\r^{l, l}     non-unit lower triangular
            !        t_{2,1}\in\r^{k-l, l}   rectangular
            !
            !        where l = floor(k/2)
            !
            !        then, consider the product:
            !
            !        (i - v_2'*t_{2,2}*v_2)*(i - v_1'*t_{1,1}*v_1)
            !        = i - v_2'*t_{2,2}*v_2 - v_1'*t_{1,1}*v_1 + v_2'*t_{2,2}*v_2*v_
            !
            !        define t_{2,1} = -t_{2,2}*v_2*v_1'*t_{1,1}
            !
            !        then, we can define the matrix v as
            !        v = |---|
            !            |v_1|
            !            |v_2|
            !            |---|
            !
            !        so, our product is equivalent to the matrix product
            !        i - v'tv
            !        this means, we can compute t_{1,1} and t_{2,2}, then use this i
            !        to compute t_{2,1}
            !
            !        compute t_{1,1} recursively
            !
            call larft(direct, storev, n-l, k-l, v, ldv, tau, t, ldt)
            !
            !        compute t_{2,2} recursively
            !
            call larft(direct, storev, n, l, v(k-l+1, 1), ldv,            &
                &               tau(k-l+1), t(k-l+1, k-l+1), ldt)
            !
            !        compute t_{2,1}
            !        t_{2,1} = v_{2,2}
            !
            call hlacpy('all', l, k-l, v(k-l+1, n-k+1), ldv,               &
                &      t(k-l+1, 1), ldt)

            !
            !        t_{2,1} = t_{2,1}*v_{1,2}'
            !
            call htrmm('right', 'lower', 'transpose', 'unit', l, k-l,      &
                &               one, v(1, n-k+1), ldv, t(k-l+1, 1), ldt)

            !
            !        t_{2,1} = v_{2,1}*v_{1,1}' + t_{2,1}
            !        note: we assume k <= n, and gemm will do nothing if n=k
            !
            call hgemm('no transpose', 'transpose', l, k-l, n-k, one,      &
                &               v(k-l+1, 1), ldv, v, ldv, one, t(k-l+1, 1),        &
                &               ldt)

            !
            !        at this point, we have that t_{2,1} = v_2*v_1'
            !        all that is left is to pre and post multiply by -t_{2,2} and t_
            !        respectively.
            !
            !        t_{2,1} = -t_{2,2}*t_{2,1}
            !
            call htrmm('left', 'lower', 'no tranpose', 'non-unit', l,      &
                &               k-l, neg_one, t(k-l+1, k-l+1), ldt,                &
                &               t(k-l+1, 1), ldt)

            !
            !        t_{2,1} = t_{2,1}*t_{1,1}
            !
            call htrmm('right', 'lower', 'no tranpose', 'non-unit', l,     &
                &               k-l, one, t, ldt, t(k-l+1, 1), ldt)
        end if
    end subroutine
end submodule
