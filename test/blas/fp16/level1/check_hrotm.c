//    SPDX-License-Identifier: LGPL-3.0-or-later
/*
   This file is part of HPBLAS, a BLAS Implementation for Half-Precision
   Copyright (C) 2025 Martin Koehler

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 3 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
   */

#include <stdio.h>
#include <stdbool.h>
#include "lpf_internal.h"
#include "fp16_helper.h"



void LPF_GLOBAL(hrotm,HROTM)(lpf_blas_int_t *, lpf_float16_t *, lpf_blas_int_t *, lpf_float16_t *, lpf_blas_int_t *, lpf_float16_t *);

void check_hrotm( bool *ok) {
#define NTEST 8
    lpf_blas_int_t n[NTEST]    = { 1, 1, 2, 4, 1, 1, 3, 7};
    lpf_blas_int_t incx[NTEST] = { 1, 1, 1, 2, 2, 1, 1, 1};
    lpf_blas_int_t incy[NTEST] = { 1, 1, 1, 2, 1, 2, 1, 1};
    lpf_float16_t dx [10] = { 0.6e0, 0.1e0, -0.5e0, 0.8e0, 0.9e0, -0.3e0,-0.4e0, 0, 0, 0};
    lpf_float16_t dy [10] = { 0.5e0, -0.9e0, 0.3e0, 0.7e0, -0.6e0, 0.2e0, 0.8e0, 0, 0 ,0};
    lpf_float16_t dpar[4][5] = {
        {-2.E0,  0.E0,  0.E0,  0.E0,  0.E0},
        {-1.E0,  2.E0, -3.E0, -4.E0,  5.E0},
        { 0.E0,  0.E0,  2.E0, -3.E0,  0.E0},
        { 1.E0,  5.E0,  2.E0,  0.E0, -4.E0}
    };

    lpf_float16_t xresult[NTEST][10] = {
        {  0.600098,   0.099976,  -0.500000,   0.799805,   0.899902,  -0.300049,  -0.399902,   0.000000,   0.000000,   0.000000  },
        { -0.799805,   0.099976,  -0.500000,   0.799805,   0.899902,  -0.300049,  -0.399902,   0.000000,   0.000000,   0.000000  },
        { -0.899902,   2.798828,  -0.500000,   0.799805,   0.899902,  -0.300049,  -0.399902,   0.000000,   0.000000,   0.000000  },
        {  3.500000,   0.099976,  -2.199219,   0.799805,   3.898438,  -0.300049,  -1.199219,   0.000000,   0.000000,   0.000000  },
        {  0.600098,   0.099976,  -0.500000,   0.799805,   0.899902,  -0.300049,  -0.399902,   0.000000,   0.000000,   0.000000  },
        { -0.799805,   0.099976,  -0.500000,   0.799805,   0.899902,  -0.300049,  -0.399902,   0.000000,   0.000000,   0.000000  },
        { -0.899902,   2.798828,  -1.400391,   0.799805,   0.899902,  -0.300049,  -0.399902,   0.000000,   0.000000,   0.000000  },
        {  3.500000,  -0.399902,  -2.199219,   4.699219,   3.898438,  -1.300781,  -1.199219,   0.000000,   0.000000,   0.000000  }
    };

    lpf_float16_t yresult[NTEST][10] = {
        {  0.500000,  -0.899902,   0.300049,   0.700195,  -0.600098,   0.199951,   0.799805,   0.000000,   0.000000,   0.000000  },
        {  0.699707,  -0.899902,   0.300049,   0.700195,  -0.600098,   0.199951,   0.799805,   0.000000,   0.000000,   0.000000  },
        {  1.700195,  -0.700195,   0.300049,   0.700195,  -0.600098,   0.199951,   0.799805,   0.000000,   0.000000,   0.000000  },
        { -2.599609,  -0.899902,  -0.700195,   0.700195,   1.500000,   0.199951,  -2.798828,   0.000000,   0.000000,   0.000000  },
        {  0.500000,  -0.899902,   0.300049,   0.700195,  -0.600098,   0.199951,   0.799805,   0.000000,   0.000000,   0.000000  },
        {  0.699707,  -0.899902,   0.300049,   0.700195,  -0.600098,   0.199951,   0.799805,   0.000000,   0.000000,   0.000000  },
        {  1.700195,  -0.700195,  -0.700195,   0.700195,  -0.600098,   0.199951,   0.799805,   0.000000,   0.000000,   0.000000  },
        { -2.599609,   3.500000,  -0.700195,  -3.601562,   1.500000,  -0.499756,  -2.798828,   0.000000,   0.000000,   0.000000  },
    };
    lpf_float16_t ytemp[10];
    lpf_float16_t xtemp[10];
    lpf_float16_t partemp[5];
    lpf_blas_int_t k;



    for ( int i = 0; i < NTEST; i++) {
        bool lok = true;
        for (k = 0; k < 10; ++k) {
            ytemp[k] = dy[k];
            xtemp[k] = dx[k];
        }
        for (k = 0; k < 5; ++k) {
            partemp[k] = dpar[i%4][k];
        }
        LPF_GLOBAL(hrotm, HROTM)(&n[i], xtemp, &incx[i], ytemp, &incy[i], partemp );

        /** printf("{"); */
        /** for (k = 0; k < 10; ++k) { */
        /**     printf("%10f%c ", (float) ytemp[k], (k == 9) ? ' ':',' ); */
        /** } */
        /** printf("},\n"); */
        /**  */

        for ( k = 0;  k < n[i]; k++) {
            lpf_float16_t mag1 = abs_fp16(xresult[i][k]);
            lpf_float16_t mag2 = abs_fp16(yresult[i][k]);

            if ( diff_fp16( xtemp[k], xresult[i][k]) > __FLT16_EPSILON__ * mag1 * 10.0 ) {
                printf("HROTM   -- FAIL -- Testcase %2d: N = %3d, INCX = %3d, INCY = %3d, K = %2d, X, computed = %10f, expected = %10f\n",
                        i+1, n[i], incx[i], incy[i], k, (float) xtemp[k], (float) xresult[i][k]);
                *ok = false;
                lok = false;
            }
            if ( diff_fp16( ytemp[k], yresult[i][k]) > __FLT16_EPSILON__ * mag2 * 10.0 ) {
                printf("HROTM   -- FAIL -- Testcase %2d: N = %3d, INCX = %3d, INCY = %3d, K = %2d, Y, computed = %10f, expected = %10f\n",
                        i+1, n[i], incx[i], incy[i], k, (float) ytemp[k], (float) yresult[i][k]);
                *ok = false;
                lok = false;
            }

        }
        if (lok) {
            printf("HROTM   -- PASS -- Testcase %2d: N = %3d, INCX = %3d, INCY = %3d\n", i+1, n[i], incx[i], incy[i]);
        }

    }
}

int main(int argc, char *argv[])
{

    (void) argc;
    (void) argv;

    bool ok = true;
    check_hrotm(&ok);
    return (ok)?0:1;
}
