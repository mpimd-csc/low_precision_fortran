subroutine qrf2(m, n, house_norm, srcA, lda, st)
    use qr_stat
    use lpf_fp16
    use lpf_lapack_fp16
    use lpf_types
    use iso_fortran_env, only: real32, real64
    implicit none

    integer(lpf_default_int_kind) :: m, n, lda
    type(qr_stats) :: st
    real(real32) :: srcA(lda, *)
    character :: house_norm

    external slange
    real(real32) :: slange

    type(fp16), allocatable :: A(:,:)
    type(fp16), allocatable :: R(:,:)
    type(fp16), allocatable :: Q(:,:)
    type(fp16), allocatable :: QQ(:,:)
    type(fp16), allocatable :: tau(:)
    type(fp16), allocatable :: work(:)

    integer(lpf_default_int_kind) :: k, l, info
    real(real64) :: tic, toc
    real(real32) :: nrm_A

    k = min(m,n)

    allocate(A(m,n))
    allocate(R(k,n))
    allocate(Q(m,m))
    allocate(QQ(m,m))
    allocate(tau(k))
    allocate(work(max(m,n)))

    st = qr_stats_init(m,n, house_norm, "geqr2")

    A(1:m, 1:n) = srcA(1:m,1:n)

    tic = lpf_get_wtime()
    call geqr2(house_norm, m, n, A, m, tau, work, info)
    toc = lpf_get_wtime()



    st % walltime = toc - tic
    st % info = info

    ! Get R
    R = 0.0
    call hlacpy("u", k, n, A, m , R, k)

    ! Get Q
    Q = 0.0
    QQ = 0.0
    do l = 1, m
        Q(l,l) = 1.0
        QQ(l,l) = 1.0
    end do
    call orm2r("L", "N", m, m, k, A, m, tau, Q, m, work, info )

    call hgemm("T", "N", m, m, m, fp16(-1.0), Q, m, Q, m, fp16(1.0), QQ, m)

    st % orth_err = lange("F", m, m, QQ, m, work) / sqrt(real(m))


    A(1:m, 1:n) = srcA(1:m,1:n)

    call hgemm("N", "N", m, n, k, fp16(-1.0), Q, m, R, k, fp16(1.0), A, m )
    nrm_A = slange("F", m, n, srcA, lda, work)
    st % reconstruct_err = lange ("F", m, n, A, m, work) / fp16(nrm_A)

    deallocate(A, R, Q, QQ, tau, work)

end subroutine

program qr_benchmark
    use qr_stat
    implicit none

    integer :: m, n, runs

    real (real32), allocatable :: A(:,:)
    character(len = 1024) :: arg
    type(qr_stats) :: st

    m = 100
    n = 100
    runs = 5

    if ( command_argument_count() .ge. 2) then
        call get_command_argument(1,  arg )
        read(arg, *) m
        call get_command_argument(2,  arg )
        read(arg, *)
    end if
    if ( command_argument_count() .eq. 3) then
        call get_command_argument(3,  arg )
        read(arg, *) runs
    endif

    write(*,'("# m    = ", I8)') m
    write(*,'("# n    = ", I8)') n
    write(*,'("# runs = ", I8)') runs

    allocate(A(m,n))

    call RANDOM_NUMBER(A)

    call print_qr_stats_header()

    call qrf2(m, n, 'L',A, m, st)

    call print_qr_stats(st)
    call qrf2(m, n, 'H',A, m, st)

    call print_qr_stats(st)

    deallocate(A)


    !     use fp16_support
!     use bf16_support
!     use hpblas_fp16
!     use hpblas_bf16
!     use hpblas_lapack_fp16
!     use hpblas_lapack_bf16
!     use hpblas_bf16_aux
!     use hpblas_fp16_aux
!     use hpblas_scale_diag
!     use iso_fortran_env
!     use printmatrix
!
!     implicit none
!     integer, parameter :: n = 3, m = 3, lda = m, ldq = m, ldt = m, ldr = m
!     integer :: info
!     integer :: k
!
!     type(fp16), dimension(lda, n) :: A_fp16, R_fp16, VV_fp16
!     type(fp16), dimension(lda, n) :: AH32_fp16, RH32_fp16, VVH32_fp16
!     type(fp16), dimension(ldq, n) :: Q_fp16, QQ_fp16
!     type(fp16), dimension(ldq, n) :: QH32_fp16, QQH32_fp16
!     type(fp16), dimension(ldt, n) :: T_fp16, TH32_fp16
!
!     type(fp16), dimension(m*n) :: work_fp16
!     type(fp16), dimension(m)   :: s_fp16, v_fp16
!     type(bf16), dimension(lda, n) :: A_bf16, R_bf16, VV_bf16
!     type(bf16), dimension(lda, n) :: AH32_bf16, RH32_bf16, VVH32_bf16
!     type(bf16), dimension(ldt, n) :: T_bf16
!     type(bf16), dimension(ldt, n) :: TH32_bf16
!     type(bf16), dimension(ldq, n) :: Q_bf16, QQ_bf16
!     type(bf16), dimension(ldq, n) :: QH32_bf16, QQH32_bf16
!     type(bf16), dimension(m*n) :: work_bf16
!     type(bf16), dimension(m) :: s_bf16, v_bf16
!
!     real(real32), dimension(lda, n) :: A_fp32, R_fp32, tmp, VV_fp32
!     real(real32), dimension(ldq, n) :: Q_fp32, QQ_fp32
!     real(real32), dimension(ldt, n) :: T_fp32
!     real(real32), dimension(n*m) :: work_fp32, s_fp32, v_fp32, tau_fp32
!
!     real(real32) :: nrm_r, nrm_h32_fp16, nrm_fp16, nrm_bf16, nrm_h32_bf16, nrm_qq_fp32
!     type(fp16) :: nrm_qq_fp16, nrm_qq_h32_fp16
!     type(bf16) :: nrm_qq_bf16, nrm_qq_h32_bf16
!
!     external :: sgeqrt2, sgeqrt, sgemqrt, slacpy, slange, sgemm, sormqr
!     real(real32) :: slange
!
!     A_fp32 = reshape([1.0, 2.0, 3.0, &
!         -2.0, -1.0, -3.0, &
!         -1.5, 2.0, 3.0], shape(A_fp32))
!
!     WRITE(*,*) "Original Matrix (FP32) :"
!     call print_matrix(A_fp32)
!     A_fp16 = A_fp32
!     A_bf16 = A_fp32
!     WRITE(*,*) "Original Matrix (FP16) :"
!     call print_matrix(A_fp16)
!     WRITE(*,*) "Original Matrix (BF16) :"
!     call print_matrix(A_bf16)
!
!     ! Scale
!     call scale_diag(m ,n, A_fp32, lda, s_fp32, v_fp32, info)
!     call scale_diag(m ,n, A_fp16, lda, s_fp16, v_fp16, info)
!     call scale_diag(m ,n, A_bf16, lda, s_bf16, v_bf16, info)
!
!     WRITE(*,*) ""
!     WRITE(*,*) "================== Scaling =================="
!     WRITE(*,*) ""
!     write(*,*) "Scaled Matrix (FP32) :"
!     call print_matrix(A_fp32)
!     write(*,*) "Scaled Matrix (FP16) :"
!     call print_matrix(A_fp16)
!     write(*,*) "Scaled Matrix (BF16) :"
!     call print_matrix(A_bf16)
!
!     AH32_fp16 = A_fp16
!     AH32_bf16 = A_bf16
!
!     ! QR Decomposition
!     T_fp32 = 0
!     WRITE(*,*) ""
!     WRITE(*,*) "================== QR Decomposition ============="
!     call sgeqrt2(m, n, A_fp32, lda, T_fp32, ldt, info)
!     call hgeqrt2(m, n, A_fp16, lda, T_fp16, ldt, info)
!     ! call hgeqrt2_v(m, n, A_fp16, lda, VV_fp16, lda, T_fp16, ldt, info)
!
!     call bgeqrt2(m, n, A_bf16, lda, T_bf16, ldt, info)
!
!     ! call hgeqrt2(m, n, AH32_fp16, lda, TH32_fp16, ldt, info )
!     call hgeqrt2_h32(m, n, AH32_fp16, lda, TH32_fp16, ldt, work_fp32, info )
!     call bgeqrt2_h32(m, n, AH32_bf16, lda, TH32_bf16, ldt, info )
!
!     R_fp32 = 0
!     R_fp16 = 0
!     RH32_fp16 = 0
!     R_bf16 = 0
!     RH32_bf16= 0
!
!     call slacpy("U", m, n, A_fp32, lda, R_fp32, ldr)
!     call hlacpy("U", m, n, A_fp16, lda, R_fp16, ldr)
!     call hlacpy("U", m, n, AH32_fp16, lda, RH32_fp16, ldr)
!     call blacpy("U", m, n, A_bf16, lda, R_bf16, ldr)
!     call blacpy("U", m, n, AH32_bf16, lda, RH32_bf16, ldr)
!
!     nrm_r = slange("F", m , n, R_fp32, ldr, work_fp32)
!     TMP = R_fp16
!     TMP = TMP - R_fp32
!     nrm_fp16 = slange("F", m, n, TMP, ldr, work_fp32) / nrm_r
!
!     TMP = RH32_fp16
!     TMP = TMP - R_fp32
!     nrm_h32_fp16 = slange("F", m, n, TMP, ldr, work_fp32) / nrm_r
!
!
!     TMP = R_bf16
!     TMP = TMP - R_fp32
!     nrm_bf16 = slange("F", m, n, TMP, ldr, work_fp32) / nrm_r
!
!     TMP = RH32_bf16
!     TMP = TMP - R_fp32
!     nrm_h32_bf16 = slange("F", m, n, TMP, ldr, work_fp32) / nrm_r
!
!
!
!     write(*,'(A)') "R (FP32) :"
!     call print_matrix(R_fp32)
!     write(*,'(A, E16.9)') "R (FP16) || R_fp32 - R_fp16 ||_F / || R_fp32 ||_F = ", nrm_fp16
!     call print_matrix(R_fp16)
!     write(*,'(A, E16.9)') "R (FP16/H32) || R_fp32 - R_fp16 ||_F / || R_fp32 ||_F = ", nrm_h32_fp16
!     call print_matrix(RH32_fp16)
!     write(*,'(A, E16.9)') "R (BF16) || R_fp32 - R_bf16 ||_F / || R_bf16 ||_F = ", nrm_bf16
!     call print_matrix(R_bf16)
!     write(*,'(A, E16.9)') "R (BF16/H32) || R_fp32 - R_bf16 ||_F / || R_bf16 ||_F = ", nrm_h32_bf16
!     call print_matrix(RH32_bf16)
!
!
!     ! VV_fp16 = 0
!     VV_fp32 = 0
!     VVH32_fp16 = 0
!     VV_bf16 = 0
!     VVH32_bf16 = 0
!
!     call slacpy("L", m, n, A_fp32, lda, VV_fp32, ldr)
!     ! call hlacpy("L", m, n, A_fp16, lda, VV_fp16, ldr)
!     call hlacpy("L", m, n, AH32_fp16, lda, VVH32_fp16, ldr)
!     call blacpy("L", m, n, A_bf16, lda, VV_bf16, ldr)
!     call blacpy("L", m, n, AH32_bf16, lda, VVH32_bf16, ldr)
!     do k = 1, m
!         VV_fp32(k,k) = 1.0
!         ! VV_fp16(k,k) = 1.0
!         VV_bf16(k,k) = 1.0
!         VVH32_fp16(k,k) = 1.0
!         VVH32_bf16(k,k) = 1.0
!     end do
!
!     write(*,*) ""
!     write(*,*) "============= V ================="
!     write(*,*) ""
!
!     write(*,'(A)') "V (FP32) :"
!     call print_matrix(VV_fp32)
!     write(*,'(A)') "V (FP16) :"
!     call print_matrix(VV_fp16)
!     write(*,'(A)') "V (FP16/H32) :"
!     call print_matrix(VVH32_fp16)
!     write(*,'(A)') "V (BF16) :"
!     call print_matrix(VV_bf16)
!     write(*,'(A)') "V (BF16/H32) :"
!     call print_matrix(VVH32_bf16)
!
!     write(*,*) ""
!     write(*,*) "============= T ================="
!     write(*,*) ""
!
!     write(*,'(A)') "T (FP32) :"
!     call print_matrix(T_fp32)
!     write(*,'(A)') "T (FP16) :"
!     call print_matrix(T_fp16)
!     write(*,'(A)') "T (FP16/H32) :"
!     call print_matrix(TH32_fp16)
!     write(*,'(A)') "T (BF16) :"
!     call print_matrix(T_bf16)
!     write(*,'(A)') "T (BF16/H32) :"
!     call print_matrix(TH32_bf16)
!
!
!     ! Setup Q
!     Q_fp32 = 0.0
!     Q_fp16 = 0.0
!     QH32_fp16 = 0.0
!     Q_bf16 = 0.0
!     QH32_bf16 = 0.0
!     do k = 1, m
!         Q_fp32(k,k) = 1.0
!         Q_fp16(k,k) = 1.0
!         Q_bf16(k,k) = 1.0
!         QH32_fp16(k,k) = 1.0
!         QH32_bf16(k,k) = 1.0
!         tau_fp32(k)  = T_fp32(k,k)
!     end do
!
!     call sgemqrt('L', 'N', m, n, m, m, A_fp32, lda, T_fp32, ldt, Q_fp32, ldq, work_fp32, info)
!     call hgemqrt('L', 'N', m, n, m, m, A_fp16, lda, T_fp16, ldt, Q_fp16, ldq, work_fp16, info)
!     call hgemqrt('L', 'N', m, n, m, m, AH32_fp16, lda, TH32_fp16, ldt, QH32_fp16, ldq, work_fp16, info)
!     call bgemqrt('L', 'N', m, n, m, m, A_bf16, lda, T_bf16, ldt, Q_bf16, ldq, work_bf16, info)
!     call bgemqrt('L', 'N', m, n, m, m, AH32_bf16, lda, TH32_bf16, ldt, QH32_bf16, ldq, work_bf16, info)
!
!     write(*,*) ""
!     write(*,*) "=======================  Q  =============="
!     write(*,*) ""
!     write(*,*) "Q_fp32"
!     call print_matrix(Q_fp32)
!     write(*,*) "Q_fp16"
!     call print_matrix(Q_fp16)
!     write(*,*) "Q_fp16 / H32"
!     call print_matrix(QH32_fp16)
!     write(*,*) "Q_bf16"
!     call print_matrix(Q_bf16)
!     write(*,*) "Q_bf16 / H32"
!     call print_matrix(QH32_bf16)
!
!
!
!     QQ_fp32 = 0.0
!     QQ_fp16 = 0.0
!     QQH32_fp16  = 0.0
!     QQ_bf16 = 0.0
!     QQH32_bf16 = 0.0
!     do k = 1, m
!         QQ_fp32(k,k) = 1.0
!         QQ_fp16(k,k) = 1.0
!         QQ_bf16(k,k) = 1.0
!         QQH32_fp16 (k,k) = 1.0
!         QQH32_bf16 (k,k) = 1.0
!     end do
!
!     call sgemm("T", "N", m, m, m, -1.0, Q_fp32, ldq, Q_fp32, ldq, 1.0, QQ_fp32, ldq)
!     call hgemm("T", "N", m, m, m, -fp16(1.0), Q_fp16, ldq, Q_fp16, ldq, fp16(1.0), QQ_fp16, ldq)
!     call hgemm("T", "N", m, m, m, -fp16(1.0), QH32_fp16, ldq, QH32_fp16, ldq, fp16(1.0), QQH32_fp16, ldq)
!     call bgemm("T", "N", m, m, m, -bf16(1.0), Q_bf16, ldq, Q_bf16, ldq, bf16(1.0), QQ_bf16, ldq)
!     call bgemm("T", "N", m, m, m, -bf16(1.0), QH32_bf16, ldq, QH32_bf16, ldq, bf16(1.0), QQH32_bf16, ldq)
!
!     nrm_qq_fp32 = slange("F", m, m, QQ_fp32, ldq, work_fp32) / sqrt(real(m))
!     nrm_qq_fp16 = hlange("F", m, m, QQ_fp16, ldq, work_fp16) / sqrt(fp16(m))
!     nrm_qq_h32_fp16 = hlange("F", m, m, QQH32_fp16, ldq, work_fp16) / sqrt(fp16(m))
!     nrm_qq_bf16 = blange("F", m, m, QQ_bf16, ldq, work_bf16) / sqrt(bf16(m))
!     nrm_qq_h32_bf16 = blange("F", m, m, QQH32_bf16, ldq, work_bf16) / sqrt(bf16(m))
!
!     WRITE(*,*) ""
!     WRITE(*,*) "=================== QQ = Q**T Q - I ========================"
!     WRITE(*,*) ""
!
!     WRITE(*,'(A, E16.9)') "QQ_FP32, nrm = ", nrm_qq_fp32
!     call print_matrix(QQ_fp32)
!     WRITE(*,'(A, DT"E"(16,9))') "QQ_FP16, nrm = ", nrm_qq_fp16
!     call print_matrix(QQ_fp16)
!     WRITE(*,'(A, DT"E"(16,9))') "QQ_FP16 / H32, nrm = ", nrm_qq_h32_fp16
!     call print_matrix(QQH32_fp16)
!
!     WRITE(*,'(A, DT"E"(16,9))') "QQ_BF16, nrm = ", nrm_qq_bf16
!     call print_matrix(QQ_bf16)
!     WRITE(*,'(A, DT"E"(16,9))') "QQ_BF16 / H32, nrm = ", nrm_qq_h32_bf16
!     call print_matrix(QQH32_bf16)
!
!
!
!
!
!
end program qr_benchmark
