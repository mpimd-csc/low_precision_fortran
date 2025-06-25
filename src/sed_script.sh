#!/bin/sh
# A small helper to convert FP16 to BF16.
# (C) Martin Koehler
# LGPLv3+
for i in fp16*
do
    x=bf16${i#fp16}
    echo "$i >> $x"
    sed -e 's,f16,bf16,g' \
        -e 's,FP16,BF16,g' \
        -e 's,_Float16,__bf16,g' \
        -e 's,fp16,bf16,g' \
        -e 's,FLT16,BFLT16,g' \
        $i > $x
done

