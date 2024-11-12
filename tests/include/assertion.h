#ifndef MPM_ASSERTION_H_
#define MPM_ASSERTION_H_

#include <cmath>

inline bool almost_equal(double a, double b, double tolerance) {
    return (a == b) || (std::abs(a - b) < tolerance) ||
           (std::abs(a) < tolerance && std::abs(b) < tolerance);
}
#endif
