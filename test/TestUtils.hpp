#include <stdio.h>
#include <stdlib.h>
#include <cmath>

namespace TestUtils {

    bool compare(double x, double y, double relativeError = 0.01) {
        return fabs(x - y) / x < relativeError;
    }

    void fail() {
        printf("A coevol automated test has failed!\n");
        exit(1);
    }

    void assert(double x, double y) {
        if (!compare(x, y)) fail();
    }
}
