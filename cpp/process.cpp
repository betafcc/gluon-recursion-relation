#include <complex>
#include <vector>

using Number = std::complex<double>;
using Vector = std::vector<Number>;
using Matrix = std::vector<Vector>;

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

// Overload to pretty print vectors
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
    os << "⟨ ";
    for (size_t i = 0; i < vec.size(); ++i)
        os << vec[i] << " ";
    os << "⟩";
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

Vector conjugate(const Vector& v)
{
    Vector result;
    result.reserve(v.size());

    for (const auto& element : v)
        result.push_back(std::conj(element));

    return result;
}

// Inner Product
Number dot(Vector a, Vector b)
{
    Number acc = 0;

    for (size_t i = 0; i < a.size(); ++i)
        acc = acc + a[i] * b[i];

    return acc;
}

Vector dot(Matrix a, Vector b)
{
    Vector acc;
    acc.reserve(b.size());

    for (size_t i = 0; i < a.size(); ++i)
        acc.push_back(dot(a[i], b));

    return acc;
}

Vector dot(Vector a, Matrix b)
{
    Vector acc;
    acc.reserve(a.size());

    for (size_t i = 0; i < b.size(); ++i)
        acc.push_back(dot(a, b[i]));

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

// ⟨q | k⟩ - Srednicki Quantum Field Theory page 356
Number braket(Vector q, Vector k)
{
    return dot(conjugate(q), k);
}

// [q | k] - Srednicki Quantum Field Theory page 356
Number box(Vector q, Vector k)
{
    return dot(q, conjugate(k));
}

// [q | g | k⟩ - Srednicki Quantum Field Theory page 357
Number fish(Vector q, Matrix g, Vector k)
{
    return dot(q, dot(g, k));
}

// ⟨q | g | k] - Srednicki Quantum Field Theory page 357
Number backfish(Vector q, Matrix g, Vector k)
{
    return dot(conjugate(q), dot(g, k));
}

class Process {
    std::vector<Gluon> gluons;

public:
    Process(std::initializer_list<Gluon> gs) :
        gluons(gs) { }

    Process(std::vector<Gluon> gs) :
        gluons(gs) { }

    auto current(std::vector<std::size_t> gis)
    {
        if (gis.size() == 1)
            return polarization(gis[0]);
    }

    // κ (2.10) from berends
    // but with CS conventions (index from 0, excluding end)
    // Add all vectors `gluons[i].vector` for i ∈ [start, end)
    Vector kappa(std::size_t start, std::size_t end) const
    {
        // TODO: improve performance by creating single empty vector of 0s and adding in-place
        Vector acc = gluons[start].momentum;

        for (auto i = start; i < end; ++i)
            acc = acc + gluons[i].momentum;

        return acc;
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
        // os << indent << "g\tλ  p\u20D7" << '\n';
        os << indent << "g\tλ  p" << '\n';
        for (auto i = 0; i < p.gluons.size(); ++i) {
            os << indent << i << "\t";

            if (p.gluons[i].helicity == Helicity::Plus)
                os << "+";
            else
                os << "-";

            os << "  " << p.gluons[i].momentum << '\n';
        }
        return os;
    }

private:
    Number current_component(std::vector<std::size_t> gis, std::size_t component)
    {
        if (gis.size() == 1)
            return polarization_component(gis[0], component);
        // else if (gis.size() == 2) {
        //     return square_brackets({ 0 }, { 1 }, component) / mp2(gluons[0].momentum + gluons[1].momentum);
        // }
    }

    // FIXME: this should be a public Gluon method
    Vector polarization(std::size_t gi)
    {
        Vector acc;
        acc.reserve(4);

        for (auto j = 0; j < 4; ++j)
            acc.push_back(polarization_component(gi, j));

        return acc;
    }

    // (60.7) and (60.8) - Srednicki Quantum Field Theory page 357
    Number polarization_component(std::size_t gi, std::size_t component)
    {
        auto gm = Gamma[component];
        auto q = auxiliar(gi); // this prevents polarization from being method of Gluon
        auto k = gluons[gi].momentum;

        switch (gluons[gi].helicity) {
        case Helicity::Plus: // (60.7)
            return (-1 / sqrt(2)) * (backfish(q, gm, k) / braket(q, k));
        case Helicity::Minus: // (60.8)
            return (-1 / sqrt(2)) * (fish(q, gm, k) / box(q, k));
        }
    }

    // Grab momentum of another gluon on the process
    // FIXME: this should not exist, it's only being used to satisfy "arbitrary massless reference momentum"
    //   - based on Srednicki Quantum Field Theory page 357
    //   it could be any non parallel vector, and a simpler one
    //   eg `is_linear_dependent(vector, {1, 1, 0, 0}) ? {1, 1, 0, 0} : {1, 0, 0, 1}`
    Vector auxiliar(std::size_t gi)
    {
        return gluons[(gi + 1) % gluons.size()].momentum;
    }
};