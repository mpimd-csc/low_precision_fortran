/*
 * Function: @FNAME@ ( C: @CNAME@ )
 */
HIDDEN void __fp16_helper_@FNAME@(int16_t *out, int16_t in)
{
    fp16_handler_t *r = (fp16_handler_t * ) out;
    fp16_handler_t _a = { .i16 = in};
    r->f16 = (_Float16) @CNAME@ ((float) _a.f16);
}

