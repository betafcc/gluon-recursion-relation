#include <complex>
#include <stdexcept>
#include "overloads.cpp"

typedef std::complex<double> Number;
typedef std::vector<Number> Vector;
typedef std::vector<std::vector<Number>> Matrix;

auto dot(Vector xs, Vector ys) -> Number
{
    Number acc = 0;

    for (auto i = 0; i < xs.size(); i++)
        acc += xs[i] * ys[i];

    return acc;
}

Vector operator+(Vector xs, Vector ys)
{
    Vector acc;
    acc.reserve(xs.size());

    for (auto i = 0; i < xs.size(); i++)
        acc[i] = xs[i] + ys[i];

    return acc;
}

Vector operator*(Number c, Vector xs)
{
    Vector acc;
    acc.reserve(xs.size());

    for (auto x : xs)
        acc.push_back(c * x);

    return acc;
}

Vector operator*(double c, Vector xs)
{
    Vector acc;
    acc.reserve(xs.size());

    for (auto x : xs)
        acc.push_back(c * x);

    return acc;
}

Vector operator*(Vector xs, Number c)
{
    Vector acc;
    acc.reserve(xs.size());

    for (auto x : xs)
        acc.push_back(c * x);

    return acc;
}

Vector operator*(Vector xs, double c)
{
    Vector acc;
    acc.reserve(xs.size());

    for (auto x : xs)
        acc.push_back(c * x);

    return acc;
}

Vector operator-(Vector xs, Vector ys)
{
    return xs + (-1 * ys);
}

Vector operator-(Vector xs)
{
    return -1 * xs;
}

std::ostream &operator<<(std::ostream &os, const Matrix &m)
{

    for (size_t i = 0; i < m.size(); ++i)
    {
        os << m[i];
        os << std::endl;
    }

    return os;
}

auto norm(Vector xs) -> Number
{
    Number acc = 0;

    for (auto x : xs)
        acc += x * x;

    return sqrt(acc);
}

auto is_linear_dependent(Vector xs, Vector ys) -> bool
{
    return abs(abs(dot(xs, ys)) - norm(xs) * norm(ys)) < 10e-6;
}
