#include <complex>
#include <iostream>
#include <vector>

typedef std::complex<double> Number;
typedef std::vector<Number> Vector;
typedef std::vector<Vector> Matrix;

enum class Helicity {
    Plus,
    Minus
};

struct Gluon {
    Helicity helicity;
    Vector momentum;
};

const Helicity Plus = Helicity::Plus;
const Helicity Minus = Helicity::Minus;

// Dirac Matrices
const std::vector<Matrix> Gamma {
    {
        { 1, 0, 0, 0 },
        { 0, 1, 0, 0 },
        { 0, 0, -1, 0 },
        { 0, 0, 0, -1 },
    },
    {
        { 0, 0, 0, 1 },
        { 0, 0, 1, 0 },
        { 0, -1, 0, 0 },
        { -1, 0, 0, 0 },
    },
    {
        { 0, 0, 0, Number(0, -1) },
        { 0, 0, Number(0, 1), 0 },
        { 0, Number(0, 1), 0, 0 },
        { Number(0, -1), 0, 0, 0 },
    },
    {
        { 0, 0, 1, 0 },
        { 0, 0, 0, -1 },
        { -1, 0, 0, 0 },
        { 0, 1, 0, 0 },
    }
};

// Overload to pretty print complex numbers
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::complex<T>& c)
{

    if (c.real() == 0 && c.imag() == 0)
        os << "0";
    else if (c.real() == 0 && c.imag() == 1)
        os << "i";
    else if (c.imag() == 0) {
        os << c.real();
    } else if (c.real() == 0) {
        os << c.imag() << "i";
    } else {
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

// Minkowski Inner Product
Number mp(Vector a, Vector b)
{
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
}

// Same as mp(a, a)
Number mp2(Vector a)
{
    return mp(a, a);
}

// [q, gx, k⟩
auto fish(Vector q, Matrix gx, Vector k) -> Number
{
    throw std::invalid_argument("Not Implemented");
}

// ⟨q | k⟩
auto braket(Vector q, Vector k) -> Number
{
    throw std::invalid_argument("Not Implemented");
}

// [q | k]
auto box(Vector q, Vector k) -> Number
{
    throw std::invalid_argument("Not Implemented");
}

class Process {
    std::vector<Gluon> gluons;

public:
    Process(std::initializer_list<Gluon> g) :
        gluons(g) { }

    // κ (2.10)
    // but with CS conventions (index from 0, excluding end)
    // Add all vectors `gluons[i].vector` for i ∈ [start, end)
    Vector kappa(std::size_t start, std::size_t end) const
    {
        // TODO improve performance by creating single empty vector of 0s and adding in-place
        Vector acc = gluons[start].momentum;

        for (auto i = start; i < end; ++i)
            acc = acc + gluons[i].momentum;

        return acc;
    }

    // Grab momentum of another gluon on the process
    Vector auxiliar(std::size_t i)
    {
        return gluons[(i + 1) % gluons.size()].momentum;
    }

    Number polarization_component(std::size_t component, std::size_t i)
    {
        auto gc = Gamma[component];
        auto q = auxiliar(i);

        switch (gluons[i].helicity) {
        case Helicity::Plus:
            return (1 / sqrt(2)) * (fish(gluons[i].momentum, gc, q) / braket(q, gluons[i].momentum));
        case Helicity::Minus:
            return (-1 / sqrt(2)) * (fish(q, gc, gluons[i].momentum) / box(q, gluons[i].momentum));
        }
    }

    // J_ξ(1) -> current_component(ξ, 0)
    // Number current_component(std::size_t component, std::size_t i)
    // {
    // }

    // J(1) -> current(0)
    // Vector current(std::size_t i)
    // {
    //     return gpv(i, gluons[i].helicity, gluons[i].vector);
    // }

    friend std::ostream& operator<<(std::ostream& os, const Process& p)
    {
        auto indent = "    ";

        os << "Process: \n";
        os << indent << "gluon helicity  momentum" << '\n';
        for (auto i = 0; i < p.gluons.size(); ++i) {
            os << indent << i << "     ";
            if (p.gluons[i].helicity == Helicity::Plus)
                os << "+";
            else
                os << "-";

            os << "         ⟨ ";
            for (const auto& component : p.gluons[i].momentum) {
                os << component << " ";
            }
            os << "⟩\n";
        }
        return os;
    }
};

int main()
{
    Process process {
        { Plus, { 1, 2, 3, 4 } },
        { Minus, { 5, 6, 7, 8 } },
    };

    std::cout << process;

    // std::cout << Number(0, 1) / Number(1, 2);
    // std::cout << x;
}