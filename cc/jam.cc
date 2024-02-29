#include "process.cc"
#include <iostream>

int main()
{
    Process process {
        { Plus, { 1, 2, 3, 4 } },
        { Minus, { 4, 5, 6, 7 } },
    };

    std::cout << process;
}
