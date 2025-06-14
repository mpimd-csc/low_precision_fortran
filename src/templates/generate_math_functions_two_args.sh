#!/usr/bin/env sh

OUTPUTC=../fp16_math_two_args.c
OUTPUTF=../fp16_math_two_args.f08

rm -f "${OUTPUTC}"

cat header.tmpl.c > "${OUTPUTC}"

while read line
do
    FNAME=$(echo $line | cut -f1 -d, )
    CNAME=$(echo $line | cut -f2 -d, )
    sed -e "s/@FNAME@/$FNAME/g" -e "s/@CNAME@/$CNAME/g" function_two_arg.tmpl.c >> "${OUTPUTC}"
done < math_functions_two_args.lst

cat > ${OUTPUTF} <<EOF
submodule (fp16_support) fp16_math_two_args
    use iso_c_binding
    use iso_fortran_env
    implicit none

    interface
EOF
while read line
do
    FNAME=$(echo $line | cut -f1 -d, )
    CNAME=$(echo $line | cut -f2 -d, )

    cat >> ${OUTPUTF} <<EOF
        pure subroutine helper_fp16_${FNAME}(out, in1, in2) bind(c, name="__fp16_helper_${FNAME}")
            use, intrinsic :: iso_c_binding
            integer(c_int16_t), intent(out) :: out
            integer(c_int16_t), intent(in), value :: in1
            integer(c_int16_t), intent(in), value :: in2
        end subroutine
EOF

done < math_functions_two_args.lst
cat >> ${OUTPUTF} <<EOF
    end interface
contains
EOF

while read line
do
    FNAME=$(echo $line | cut -f1 -d, )
    CNAME=$(echo $line | cut -f2 -d, )

    cat >> ${OUTPUTF} << EOF
    module elemental function ${FNAME}_fp16(in1, in2) result(out)
            type(FP16) :: out
            type(FP16), intent(in) :: in1
            type(FP16), intent(in) :: in2

            call helper_fp16_${FNAME}(out%value, in1%value, in2%value)
    end function
EOF
    cat << EOF
    interface ${FNAME}
        module elemental function ${FNAME}_fp16(in1, in2) result(out)
            type(FP16), intent(in) :: in1
            type(FP16), intent(in) :: in2
            type(FP16) :: out
        end function ${FNAME}_fp16
    end interface
EOF


done < math_functions_two_args.lst
cat >> ${OUTPUTF} <<EOF
end submodule fp16_math_two_args
EOF

while read line
do
    FNAME=$(echo $line | cut -f1 -d, )
    CNAME=$(echo $line | cut -f2 -d, )
    echo "    PUBLIC :: ${FNAME}"
done < math_functions_two_args.lst

