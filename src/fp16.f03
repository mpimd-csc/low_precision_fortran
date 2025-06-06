MODULE FP16_SUPPORT
    USE iso_c_binding
    USE iso_fortran_env
    IMPLICIT NONE

    PRIVATE 

    PUBLIC :: FP16, ASSIGNMENT(=), write(formatted), operator(+), operator(-), operator(*), operator(/)
    PUBLIC :: operator(**), abs


    TYPE, BIND(C) :: FP16
        INTEGER(c_int16_t) :: value
    END TYPE

    interface FP16 
        module procedure :: construct_real
        module procedure :: construct_double
        module procedure :: construct_int
    end interface

    INTERFACE ASSIGNMENT(=)
        MODULE PROCEDURE assign_int, assign_real, assign_double, assign_fp16, assign_to_real
    END INTERFACE

    INTERFACE write(formatted)
        MODULE PROCEDURE  write_formatted
    END INTERFACE

    ! Interface für den + Operator
    interface operator(+)
        module procedure add_fp16_fp16, add_fp16_real, add_real_fp16
    end interface operator(+)

    ! Interface für den - Operator
    interface operator(-)
        module procedure subtract_fp16_fp16, subtract_fp16_real, subtract_real_fp16
    end interface operator(-)

    ! Interface für den * Operator
    interface operator(*)
        module procedure multiply_fp16_fp16, multiply_fp16_real, multiply_real_fp16
    end interface operator(*)

    ! Interface für den / Operator
    interface operator(/)
        module procedure divide_fp16_fp16, divide_fp16_real, divide_real_fp16
    end interface operator(/)

    ! Interface für den ** Operator
    interface operator(**) 
        module procedure power_fp16_fp16, power_fp16_real, power_fp16_int
    end interface operator(**)

    interface abs
        module procedure abs_fp16
    end interface

    ! C Interfaces 
    INTERFACE 
        SUBROUTINE SET_FP16_FROM_INT(OUT, IN) BIND(C, name = "__fp16_helper_set_from_int")
            USE, INTRINSIC:: iso_c_binding
            INTEGER(c_int16_t), INTENT(OUT) :: OUT 
            INTEGER(c_int), VALUE, INTENT(IN) :: IN 
        END SUBROUTINE 
        SUBROUTINE SET_FP16_FROM_REAL(OUT, IN) BIND(C, name = "__fp16_helper_set_from_float")
            USE, INTRINSIC:: iso_c_binding
            INTEGER(c_int16_t), INTENT(OUT) :: OUT 
            REAL(c_float), VALUE, INTENT(IN) :: IN 
        END SUBROUTINE 
        SUBROUTINE SET_FP16_FROM_DOUBLE(OUT, IN) BIND(C, name = "__fp16_helper_set_from_double")
            USE, INTRINSIC:: iso_c_binding
            INTEGER(c_int16_t), INTENT(OUT) :: OUT 
            REAL(c_double), VALUE, INTENT(IN) :: IN 
        END SUBROUTINE 
        REAL(c_float) FUNCTION GET_FP16(IN) BIND(C, name = "__fp16_helper_get_float")
            USE, INTRINSIC :: iso_c_binding
            INTEGER(c_int16_t), INTENT(IN), VALUE :: IN 
        END FUNCTION 
    
        subroutine helper_add_fp16_fp16(out, a, b) bind(c, name = "__fp16_helper_add_fp16_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine 
        subroutine helper_add_fp16_real(out, a, b) bind(c, name = "__fp16_helper_add_fp16_real")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: a
            real(c_float), intent(in), value :: b
        end subroutine 

        subroutine helper_mul_fp16_fp16(out, a, b) bind(c, name = "__fp16_helper_mul_fp16_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine 
        subroutine helper_mul_fp16_real(out, a, b) bind(c, name = "__fp16_helper_mul_fp16_real")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: a
            real(c_float), intent(in), value :: b
        end subroutine 

        subroutine helper_sub_fp16_fp16(out, a, b) bind(c, name = "__fp16_helper_sub_fp16_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine 
        subroutine helper_sub_fp16_real(out, a, b) bind(c, name = "__fp16_helper_sub_fp16_real")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: a
            real(c_float), intent(in), value :: b
        end subroutine 
        subroutine helper_sub_real_fp16(out, a, b) bind(c, name = "__fp16_helper_sub_real_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            real(c_float), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine 

        subroutine helper_div_fp16_fp16(out, a, b) bind(c, name = "__fp16_helper_div_fp16_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine 
        subroutine helper_div_fp16_real(out, a, b) bind(c, name = "__fp16_helper_div_fp16_real")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: a
            real(c_float), intent(in), value :: b
        end subroutine 
        subroutine helper_div_real_fp16(out, a, b) bind(c, name = "__fp16_helper_div_real_fp16")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            real(c_float), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine 
    
        subroutine helper_power_fp16_fp16(out, a, b) bind(c, name = "__fp16_helper_pow_fp16_fp16")
             use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: a
            integer(c_int16_t), intent(in), value :: b
        end subroutine 
        subroutine helper_power_fp16_real(out, a, b) bind(c, name = "__fp16_helper_pow_fp16_real")
             use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: a
            real(c_float), intent(in), value :: b
        end subroutine 
        subroutine helper_power_fp16_int(out, a, b) bind(c, name = "__fp16_helper_pow_fp16_int")
             use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: a
            integer(c_int), intent(in), value :: b
        end subroutine 

        subroutine helper_abs_fp16(out, in) bind(c, name = "__fp16_helper_abs")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out 
            integer(c_int16_t), intent(in), value :: in 
        end subroutine 

    END INTERFACE



CONTAINS
    ! Contructors 
    function construct_real(in) result(r)
        real(real32), intent(in) :: in
        type(fp16) :: r
        call SET_FP16_FROM_REAL(r%value, in)
    end function 

    function construct_double(in) result(r)
        real(real64), intent(in) :: in
        type(fp16) :: r
        call SET_FP16_FROM_DOUBLE(r%value, in)
    end function 

    function construct_int(in) result(r)
        integer, intent(in) :: in
        type(fp16) :: r
        call SET_FP16_FROM_INT(r%value, in)
    end function 


    !
    ! Assignment Operator 
    !
    SUBROUTINE assign_int(this, that)
        TYPE(FP16), INTENT(OUT) :: this
        INTEGER, INTENT(IN) :: that
        CALL SET_FP16_FROM_INT(this%value, that)
    END SUBROUTINE

    SUBROUTINE assign_real(this, that)
        TYPE(FP16), INTENT(OUT) :: this
        REAL(real32), INTENT(IN) :: that
        CALL SET_FP16_FROM_REAL(this%value, that) 
    END SUBROUTINE

    SUBROUTINE assign_double(this, that)
        TYPE(FP16), INTENT(OUT) :: this
        REAL(real64), INTENT(IN) :: that
        CALL SET_FP16_FROM_DOUBLE(this%value, that)
    END SUBROUTINE

    SUBROUTINE assign_fp16(this, that)
        TYPE(FP16), INTENT(OUT) :: this
        TYPE(FP16), INTENT(IN) :: that
        this%value = that%value
    END SUBROUTINE

    SUBROUTINE assign_to_real(this, that)
        REAL(real32), INTENT(OUT) :: this
        TYPE(FP16), INTENT(IN) :: that
        this = GET_FP16(that%value)
    END SUBROUTINE

    !
    ! Operator(+)
    !
    function add_fp16_fp16(this, that) result(sum)
        type(fp16), intent(in) :: this
        type(fp16), intent(in) :: that 
        type(fp16) :: sum 
        
        call helper_add_fp16_fp16(sum%value, this%value, that%value)
    end function 

    function add_fp16_real(this, that) result(sum)
        type(fp16), intent(in) :: this
        real(real32), intent(in) :: that 
        type(fp16) :: sum 
        call helper_add_fp16_real(sum%value, this%value, that)
    end function 

    function add_real_fp16(this, that) result(sum)
        real(real32), intent(in) :: this
        type(fp16), intent(in) :: that 
        type(fp16) :: sum 
        call helper_add_fp16_real(sum%value, that%value, this)
    end function 

    !
    ! Operator (-)
    !
    function subtract_fp16_fp16(this, that) result(diff)
        type(fp16), intent(in) :: this
        type(fp16), intent(in) :: that 
        type(fp16) :: diff 

        call helper_sub_fp16_fp16(diff%value, this%value, that%value)
    end function 

    function subtract_fp16_real(this, that) result(diff)
        type(fp16), intent(in) :: this
        real(real32), intent(in) :: that 
        type(fp16) :: diff 

        call helper_sub_fp16_real(diff%value, this%value, that)
    end function 

    function subtract_real_fp16(this, that) result(diff)
        real(real32), intent(in) :: this
        type(fp16), intent(in) :: that 
        type(fp16) :: diff 

        call helper_sub_real_fp16(diff%value, this, that%value)
    end function 

    !
    ! Operator(*)
    !
    function multiply_fp16_fp16(this, that) result(prod)
        type(fp16), intent(in) :: this
        type(fp16), intent(in) :: that 
        type(fp16) :: prod 

        call helper_mul_fp16_fp16(prod%value, this%value, that%value)
    end function 

    function multiply_fp16_real(this, that) result(prod)
        type(fp16), intent(in) :: this
        real(real32), intent(in) :: that 
        type(fp16) :: prod 

        call helper_mul_fp16_real(prod%value, this%value, that)
    end function 

    function multiply_real_fp16(this, that) result(prod)
        real(real32), intent(in) :: this
        type(fp16), intent(in) :: that 
        type(fp16) :: prod 

        call helper_mul_fp16_real(prod%value, that%value, this)
    end function 

    !
    ! Operator(/)
    !
    function divide_fp16_fp16(this, that) result(quot)
        type(fp16), intent(in) :: this
        type(fp16), intent(in) :: that 
        type(fp16) :: quot 

        call helper_div_fp16_fp16(quot%value, this%value, that%value)
    end function 

    function divide_fp16_real(this, that) result(quot)
        type(fp16), intent(in) :: this
        real(real32), intent(in) :: that 
        type(fp16) :: quot 

        call helper_div_fp16_real(quot%value, this%value, that)
    end function 

    function divide_real_fp16(this, that) result(quot)
        real(real32), intent(in) :: this
        type(fp16), intent(in) :: that 
        type(fp16) :: quot 

        call helper_div_real_fp16(quot%value, this, that%value)
    end function 

    !
    ! Operator(**)
    !
    function power_fp16_fp16(this, that) result(power)
        type(fp16), intent(in) :: this
        type(fp16), intent(in) :: that 
        type(fp16) :: power 

        call helper_power_fp16_fp16(power%value, this%value, that%value)
    end function 

    function power_fp16_real(this, that) result(power)
        type(fp16), intent(in) :: this
        real(real32), intent(in) :: that 
        type(fp16) :: power 

        call helper_power_fp16_real(power%value, this%value, that)
    end function 

    function power_fp16_int(this, that) result(power)
        type(fp16), intent(in) :: this
        integer, intent(in) :: that
        type(fp16) :: power 

        call helper_power_fp16_int(power%value, this%value, that)
    end function 

    !
    ! Functions 
    !

    ! ABS 
    function abs_fp16(in) result (r)
        type(fp16), intent(in) :: in 
        type(fp16) :: r

        call helper_abs_fp16(r%value, in %value)
    end function 
   



    ! Formated output
    subroutine write_formatted(dtv, unit, iotype, v_list, iostat, iomsg)                                                       
        type(FP16), intent(in) :: dtv                                                                                                
        integer, intent(in) :: unit                                                                                                      
        character(*), intent(in) :: iotype                                                                                               
        integer, intent(in) :: v_list(:)                                                                                                 
        integer, intent(out) :: iostat                                                                                                   
        character(*), intent(inout), optional :: iomsg                                                                                   
        character(10) :: pfmt
        
        REAL(real32) :: x
        x = GET_FP16(dtv%value)
        
        if (iotype == 'LISTDIRECTED') then                                                                                               
            write(unit, '(F0.4)', iostat=iostat, iomsg=iomsg) x 
        else if (iotype == 'DT') then     
            if (size(v_list) == 0 ) then 
                write(unit, '(F0.4)', iostat=iostat, iomsg=iomsg) x               
            else if (size(v_list) == 1 ) then
                write(pfmt, '(a,i2,a)') '(F0.', v_list(1),')'
                write(unit, pfmt, iostat = iostat, iomsg = iomsg) x 
            else if (size(v_list) == 2 ) then 
                write(pfmt, '(a,i0,a,i0,a)') '(F', v_list(1),'.',v_list(2),')'
                write(unit, pfmt, iostat = iostat, iomsg = iomsg) x 
            else 
                iostat = -1 
                if (present(iomsg)) iomsg = 'Too many options in DT setting'                                                                               
            endif 
        else                                                                                                                             
            iostat = -1                                                                                                                    
            if (present(iomsg)) iomsg = 'Unsupported iotype'                                                                               
        end if                                                                                                                           
    end subroutine write_formatted

END MODULE
