//
// Created by Sahand Adibnia on 3/30/24.
//

#include "Utils.h"

double factorial(int n) {
    if (n == 1 or n == 0) {
        return 1;
    } else {
        return n*factorial(n-1);
    }
}

double double_factorial(int n) {
    if (n <= 1) {
        return 1;
    } else if (n == 2) {
        return 2;
    } else {
        return n*double_factorial(n-2);
    }
}

double nChooseR(int n, int r) {
    return factorial(n)/(factorial(n-r)*factorial(r));
}