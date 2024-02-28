#include <vector>
#include <iostream>
#include "util.cpp"

#include "complex"

#include "qcd.cpp"

using std::complex;
using std::vector;

int main()
{
    // std::complex<double> a(1, 1);
    Vector xs = {1.0, 2.0, 3.0, 4.0, 5.0};
    Vector sliced(xs.begin() + 0, xs.begin() + 3);

    for (auto i = 0; i < 10; i++)
    {
        std::cout << i << std::endl;
    }

    // std::cout << sliced << std::endl;

    return 0;
}
