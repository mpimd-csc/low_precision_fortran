#include <stdio.h>
#include <stdlib.h>
#include <nanobench.h>

extern "C" void test_linspace_();

int main(int argc, char *argv[])
{
    ankerl::nanobench::Bench b;
    b.performanceCounters(true);
    b.clockResolutionMultiple(10000000);
    b.run("fp16_linspace", [](){ test_linspace_();});

    return 0;
}
