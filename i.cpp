#include <complex>

template <typename T>
auto i(T z) -> std::complex<T>
{
    return std::complex<T>(0, z);
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::complex<T> &c)
{

    if (c.real() == 0 && c.imag() == 0)
        os << "0";
    else if (c.real() == 0 && c.imag() == 1)
        os << "i";
    else if (c.imag() == 0)
    {
        os << c.real();
    }
    else if (c.real() == 0)
    {
        os << c.imag() << "i";
    }
    else
    {
        os << c.real();
        if (c.imag() < 0)
            os << " - ";
        else
            os << " + ";
        os << abs(c.imag());
        os << "i";
    }

    return os;
}
