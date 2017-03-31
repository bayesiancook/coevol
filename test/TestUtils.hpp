#include <stdio.h>
#include <stdlib.h>
#include <cmath>

namespace TestUtils {

    bool compare(double x, double y, double relativeError = 0.01) {
        return fabs(x - y) / x < relativeError;
    }

    void fail() {
        printf("-- A coevol automated test has failed!\n");
        exit(1);
    }

    void cassert(double x, double y, double relativeError = 0.01) {
        if (!compare(x, y, relativeError)) {
            printf("-- Expected valued %f (+/- %.2f%%) but got %f instead!\n", x,
                   relativeError * 100, y);
            fail();
        }
    }
}
