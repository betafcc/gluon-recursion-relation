#include <complex>
#include <vector>
#include <iostream>

using Number = std::complex<double>;
using Vector = std::vector<Number>;

// Overload to perform mathematical vector addition
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    std::vector<T> acc;
    acc.reserve(a.size());

    for (size_t i = 0; i < a.size(); ++i)
        acc.push_back(a[i] + b[i]);

    return acc;
}

// Overload to perform mathematical vector subtraction
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
    std::vector<T> acc;
    acc.reserve(a.size());

    for (size_t i = 0; i < a.size(); ++i)
        acc.push_back(a[i] - b[i]);

    return acc;
}

// Overload to perform mathematical scalar multiplication to a vector
Vector operator*(const Number& a, const Vector& b)
{
    Vector acc;
    acc.reserve(b.size());

    for (size_t i = 0; i < b.size(); ++i)
        acc.push_back(a * b[i]);

    return acc;
}

Vector operator*(const Vector& b, const Number& a)
{
    Vector acc;
    acc.reserve(b.size());

    for (size_t i = 0; i < b.size(); ++i)
        acc.push_back(a * b[i]);

    return acc;
}

Vector conjugate(const Vector& v)
{
    Vector result;
    result.reserve(v.size());

    for (const auto& element : v)
        result.push_back(std::conj(element));

    return result;
}