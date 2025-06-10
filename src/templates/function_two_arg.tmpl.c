/*
 * Function: @FNAME@ ( C: @CNAME@ )
 */
HIDDEN void __fp16_helper_@FNAME@(int16_t *out, int16_t in1, int16_t in2)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in1};
    fp16_handler_t _b = { .i16 = in2};
    r->f16 = (_Float16) @CNAME@ ((float) _a.f16, (float) _b.f16);
}

