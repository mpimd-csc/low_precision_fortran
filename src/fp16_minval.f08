submodule (fp16_support) fp16_minval
    use iso_c_binding
    use iso_fortran_env
    implicit none

contains
    ! Overall min in 1D
    pure function minval_fp16_1d(array) result(min_value)
        type(fp16), dimension(:), intent(in) :: array
        type(fp16) :: min_value
        integer :: i
        type(fp16) :: hval;


        if (size(array) .eq. 0) then
            hval = fp16(1.0)
            min_value = huge(hval)
            return
        end if

        min_value = array(1)

        do i = 2, size(array)
            if (array(i) .lt. min_value) then
                min_value = array(i)
            end if
        end do
    end function minval_fp16_1d

    ! Overall min in 2D
    pure function minval_fp16_2d(array) result(min_value)
        type(fp16), dimension(:,:), intent(in) :: array
        type(fp16) :: min_value
        integer :: i, j

        min_value = array(1, 1)

        do i = 1, size(array, 1)
            do j = 1, size(array, 2)
                if (array(i, j) .lt. min_value) then
                    min_value = array(i, j)
                end if
            end do
        end do
    end function minval_fp16_2d

    ! Overall min in 3D
    pure function minval_fp16_3d(array) result(min_value)
        type(fp16), dimension(:,:,:), intent(in) :: array
        type(fp16) :: min_value
        integer :: i, j, k

        min_value = array(1, 1, 1)

        do i = 1, size(array, 1)
            do j = 1, size(array, 2)
                do k = 1, size(array, 3)
                    if (array(i, j, k) .lt. min_value) then
                        min_value = array(i, j, k)
                    end if
                end do
            end do
        end do
    end function minval_fp16_3d

    ! Overall min in 4D
    pure function minval_fp16_4d(array) result(min_value)
        type(fp16), dimension(:,:,:,:), intent(in) :: array
        type(fp16) :: min_value
        integer :: i, j, k, l

        min_value = array(1, 1, 1,1)

        do i = 1, size(array, 1)
            do j = 1, size(array, 2)
                do k = 1, size(array, 3)
                    do l = 1, size(array, 4)
                        if (array(i, j, k,l) .lt. min_value) then
                            min_value = array(i, j, k, l)
                        end if
                    end do
                end do
            end do
        end do
    end function minval_fp16_4d


    ! Overall Max in 1D but with dim argument.
    pure function minval_fp16_1d_dim(array, dim) result(min_value)
        type(fp16), dimension(:), intent(in) :: array
        integer, intent(in) :: dim
        type(fp16) :: min_value

        if (dim .ne. 1) then
            error stop 'Invalid dimension for 1D array'
        end if

        min_value = minval_fp16_1d(array)
    end function minval_fp16_1d_dim

    ! Max along dim in 2D
    pure function minval_fp16_2d_dim(array, dim) result(min_value)
        type(fp16), dimension(:,:), intent(in) :: array
        integer, intent(in) :: dim
        type(fp16), dimension(size(array, merge(2, 1, dim == 1))) :: min_value
        integer :: i, j

        if (dim < 1 .or. dim > 2) then
            error stop 'Invalid dimension for 2D array'
        end if

        ! Initialisierung mit den Elementen des Arrays
        if (dim == 1) then
            do j = 1, size(array, 2)
                min_value(j) = array(1, j)
                do i = 2, size(array, 1)
                    if (array(i, j) .lt. min_value(j)) then
                        min_value(j) = array(i, j)
                    end if
                end do
            end do
        else
            do i = 1, size(array, 1)
                min_value(i) = array(i,1)
                do j = 2, size(array, 2)
                    if (array(i, j) .lt. min_value(i)) then
                        min_value(i) = array(i, j)
                    end if
                end do
            end do
        end if
    end function minval_fp16_2d_dim

    ! Max along dim in 3D
    pure function minval_fp16_3d_dim(array, dim) result(min_value)
        type(fp16), dimension(:,:,:), intent(in) :: array
        integer, intent(in) :: dim
        !
        ! The input is a m x n x k array.
        ! if dim == 1, the output is n x k, we need dimension 2 and 3
        ! if dim == 2, the output is m x k, we need dimension 1 and 3
        ! if dim == 3, the output is m x n, we need dimension 1 and 2
        !
        type(fp16), dimension( size(array, merge(2, 1, dim == 1)), &
            & size(array, merge(2, 3, dim == 3))) :: min_value
        integer :: i, j, k

        if (dim < 1 .or. dim > 3) then
            error stop 'Invalid dimension for 3D array'
        end if

        if (dim == 1) then
            ! dim == 1, for each (j,k) iterate over i
            do j = 1, size(array, 2)
                do k = 1, size(array, 3)
                    min_value(j,k) = array(1, j, k)
                    do i = 2, size(array, 1)
                        if (array(i, j, k) .lt. min_value(j,k)) then
                            min_value(j,k) = array(i, j, k)
                        end if
                    end do
                end do
            end do
        elseif (dim == 2) then
            ! dim == 2, for each (i,k) iterate over j
            do i = 1, size(array, 1)
                do k = 1, size(array, 3)
                    min_value(i,k) = array(i, 1, k)
                    do j = 2, size(array, 2)
                        if (array(i, j, k) .lt. min_value(i,k)) then
                            min_value(i,k) = array(i, j, k)
                        end if
                    end do
                end do
            end do
        else
            ! dim == 3, for each (i,j) iterate over k
            do i = 1, size(array, 1)
                do j = 1, size(array, 2)
                    min_value(i,j) = array(i,j,1)
                    do k = 2, size(array, 3)
                        if (array(i, j, k) .lt. min_value(i,j)) then
                            min_value(i,j) = array(i, j, k)
                        end if
                    end do
                end do
            end do
        end if

    end function minval_fp16_3d_dim

    ! Max along dim in 4D
    pure function minval_fp16_4d_dim(array, dim) result(min_value)
        type(fp16), dimension(:,:,:,:), intent(in) :: array
        integer, intent(in) :: dim
        !
        ! The input is a m x n x k x l array.
        ! if dim == 1, the output is n x k x l, we need dimension 2, 3, 4
        ! if dim == 2, the output is m x k x l, we need dimension 1, 3, 4
        ! if dim == 3, the output is m x n x l, we need dimension 1, 2, 4
        ! if dim == 4, the output is m x n x k, we need dimension 1, 2, 3
        !
        type(fp16), dimension( size(array, merge(2, 1, dim == 1)), &
            & size(array, merge(3, 2, dim < 3)), &
            & size(array, merge(4, 3, dim == 4))) :: min_value
        integer :: i, j, k, l

        if (dim < 1 .or. dim > 4) then
            error stop 'Invalid dimension for 3D array'
        end if

        if (dim == 1) then
            ! dim == 1, for each (j,k,l) iterate over i
            do j = 1, size(array, 2)
                do k = 1, size(array, 3)
                    do l = 1, size(array, 4)
                        min_value(j,k,l) = array(1, j, k, l)
                        do i = 2, size(array, 1)
                            if (array(i, j, k, l) .lt. min_value(j,k,l)) then
                                min_value(j,k,l) = array(i, j, k,l)
                            end if
                        end do
                    end do
                end do
            end do
        elseif (dim == 2) then
            ! dim == 2, for each (i,k,l) iterate over j
            do i = 1, size(array, 1)
                do k = 1, size(array, 3)
                    do l =1, size(array, 4)
                        min_value(i,k,l) = array(i, 1, k,l)
                        do j = 2, size(array, 2)
                            if (array(i, j, k,l) .lt. min_value(i,k,l)) then
                                min_value(i,k, l) = array(i, j, k, l)
                            end if
                        end do
                    end do
                end do
            end do
        elseif ( dim == 3) then
            ! dim == 3, for each (i,j, l) iterate over k
            do i = 1, size(array, 1)
                do j = 1, size(array, 2)
                    do l = 1, size(array, 4)
                        min_value(i,j, l) = array(i,j,1, l)
                        do k = 2, size(array, 3)
                            if (array(i, j, k, l) .lt. min_value(i,j, l)) then
                                min_value(i,j,l) = array(i, j, k, l)
                            end if
                        end do
                    end do
                end do
            end do
        else
            ! dim == 4, for each (i,j, k) iterate over l
            do i = 1, size(array, 1)
                do j = 1, size(array, 2)
                    do k = 2, size(array, 3)
                        min_value(i,j, k) = array(i,j,k,1)
                        do l = 1, size(array, 4)
                            if (array(i, j, k, l) .lt. min_value(i,j, k)) then
                                min_value(i,j,k) = array(i, j, k, l)
                            end if
                        end do
                    end do
                end do
            end do

        end if

    end function minval_fp16_4d_dim



end submodule fp16_minval
