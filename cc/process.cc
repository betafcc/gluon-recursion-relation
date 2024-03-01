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
        return os << '0';
    if (c.real() != 0)
        os << c.real();
    if (c.imag() != 0) {
        if (c.imag() > 0 && c.real() != 0)
            os << '+';
        if (c.imag() == 1)
            os << 'i';
        else if (c.imag() == -1)
            os << "-i";
        else
            os << c.imag() << 'i';
    }
    return os;
}

// Overload to pretty print vectors
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
    os << "⟨ ";
    for (size_t i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i != vec.size() - 1)
            os << " ";
    }
    os << " ⟩";
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

    Number current(std::vector<std::size_t> gis, std::size_t xi)
    {
        if (gis.size() == 1)
            return polarization(gis[0], xi);
        else if (gis.size() == 2)
            // berends (2.7)
            return (1.0 / mp2(gluons[0].momentum + gluons[1].momentum)) * sb({ 0 }, { 1 }, xi);
        else if (gis.size() == 3)
            // berends (2.11)
            return (1.0 / mp2(kappa(gis))) * (sb({ 1 }, { 2, 3 }, xi) + sb({ 1, 2 }, { 3 }, xi));
        else
            throw std::runtime_error("Not implemented");
    }

    Vector current(std::vector<std::size_t> gis)
    {
        if (gis.size() == 1)
            return polarization(gis[0]);
        else if (gis.size() == 2)
            return {
                current(gis, 0),
                current(gis, 1),
                current(gis, 2),
                current(gis, 3),
            };
        else if (gis.size() == 3)
            // berends (2.14)
            return (1.0 / mp2(kappa(gis))) * (sb({ 0 }, { 1, 2 }) + sb({ 0, 1 }, { 2 }) + cb({ 0 }, { 1 }, { 2 }));
        else {
            // berends (2.19)
            auto factor = 1.0 / mp2(kappa(gis));
            auto n = gis.size();

            // FIXME: we can use only one acc
            auto sb_acc = Vector(4, 0);
            for (auto m = 0; m < (n - 1); ++m)
                sb_acc = sb_acc + sb( // clang-format off
                    std::vector<std::size_t>(gis.begin(), gis.begin() + m + 1),
                    std::vector<std::size_t>(gis.begin() + m + 1, gis.begin() + n)
                ); // clang-format on

            auto cb_acc = Vector(4, 0);
            for (auto m = 0; m < (n - 2); ++m)
                for (auto k = m + 1; k < (n - 1); ++k)
                    cb_acc = cb_acc + cb( // clang-format off
                        std::vector<std::size_t>(gis.begin(), gis.begin() + m + 1),
                        std::vector<std::size_t>(gis.begin() + m + 1, gis.begin() + k + 1),
                        std::vector<std::size_t>(gis.begin() + k + 1, gis.begin() + n)
                    ); // clang-format on

            return factor * (sb_acc + cb_acc);
        }
    }

    // κ (2.10) from berends
    // but with CS conventions (index from 0, excluding end)
    // Add all vectors `gluons[i].momentum` for i ∈ [start, end)
    Vector kappa(std::size_t start, std::size_t end) const
    {
        // TODO: improve performance by creating single empty vector of 0s and adding in-place
        Vector acc = gluons[start].momentum;

        for (auto i = start; i < end; ++i)
            acc = acc + gluons[i].momentum;

        return acc;
    }

    // this version takes a vector of indexes,
    // since it's not clear how kappa should behave when indexes out of order
    // or not continuous
    // eg κ({1, 4, 3}) != κ(1, 3)
    Vector kappa(std::vector<std::size_t> gis) const
    {
        // TODO: improve performance by creating single empty vector of 0s and adding in-place
        Vector acc = gluons[gis[0]].momentum;

        for (auto i = 1; i < gis.size(); ++i)
            acc = acc + gluons[gis[i]].momentum;

        return acc;
    }

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
    /*
    berends (2.7) = (2.4)
    [J(1), J(2)]_ξ
        = 2 K₂ . J(1) J_ξ(2)
        - 2 K₁ . J(2) J_ξ(1)
        + (K₁ - K₂)_ξ J(1) . J(2)

    =>

    [J(1 ... m), J(m + 1 ... n)]_ξ
        = 2 κ(m + 1, n) . J(1 ... m) J_ξ(m + 1 ... n)
        - 2 κ(1, m)     . J(m + 1 ... n) J_ξ(1 ... m)
        + (κ(1, m) - κ(m + 1, n))_ξ J(1 ... m) . J(m + 1 ... n)

    =>

    [J(xs), J(ys)][xi]
        = 2 κ(ys) . J(xs) J(ys)[xi]
        - 2 κ(xs) . J(ys) J(xs)[xi]
        + (κ(xs) - κ(ys))[xi] J(xs) . J(ys)
    */
    Number sb(
        std::vector<std::size_t> xs,
        std::vector<std::size_t> ys,
        std::size_t xi
    )
    {
        return dot(2 * kappa(ys), current(xs) * current(ys, xi))
            - dot(2 * kappa(xs), current(ys) * current(xs, xi))
            + (kappa(xs) - kappa(ys))[xi] * dot(current(xs), current(ys));
    }

    Vector sb(
        std::vector<std::size_t> xs,
        std::vector<std::size_t> ys
    )
    {
        return {
            sb(xs, ys, 0),
            sb(xs, ys, 1),
            sb(xs, ys, 2),
            sb(xs, ys, 3),
        };
    }

    /*
    berends (2.13)

    {J(1), J(2), J(3)}_ξ
      = J(1) . (J(3) J_ξ(2) - J(2) J_ξ(3))
      - J(3) . (J(2) J_ξ(1) - J(1) J_ξ(2))

    =>

    {J(a), J(b), J(c)}[xi]
      = J(a) . (J(c) * J(b)[xi] - J(b) * J(c)[xi])
      - J(c) . (J(b) * J(a)[xi] - J(a) * J(b)[xi])
    */
    Number cb(
        std::vector<std::size_t> a,
        std::vector<std::size_t> b,
        std::vector<std::size_t> c,
        std::size_t xi
    )
    {

        return dot(current(a), current(c) * current(b, xi) - current(b) * current(c, xi))
            - dot(current(c), current(b) * current(a, xi) - current(a) * current(b, xi));
    }

    Vector cb(
        std::vector<std::size_t> a,
        std::vector<std::size_t> b,
        std::vector<std::size_t> c
    )
    {
        return {
            cb(a, b, c, 0),
            cb(a, b, c, 1),
            cb(a, b, c, 2),
            cb(a, b, c, 3),
        };
    }

    // FIXME: this should be a Gluon method
    // (60.7) and (60.8) - Srednicki Quantum Field Theory page 357
    Number polarization(std::size_t gi, std::size_t xi)
    {
        auto gm = Gamma[xi];
        auto q = auxiliar(gi); // this prevents polarization from being method of Gluon
        auto k = gluons[gi].momentum;

        switch (gluons[gi].helicity) {
        case Helicity::Plus: // (60.7)
            return (-1 / sqrt(2)) * (backfish(q, gm, k) / braket(q, k));
        case Helicity::Minus: // (60.8)
            return (-1 / sqrt(2)) * (fish(q, gm, k) / box(q, k));
        }
    }

    Vector polarization(std::size_t gi)
    {
        return {
            polarization(gi, 0),
            polarization(gi, 1),
            polarization(gi, 2),
            polarization(gi, 3),
        };
    }

    // Grab momentum of another gluon on the process
    // FIXME: this should not exist, it's only being used to satisfy "arbitrary massless reference momentum"
    //   - based on Srednicki Quantum Field Theory page 357
    //   it could be any non-parallel vector, and a simpler one
    //   eg `is_linear_dependent(vector, {1, 1, 0, 0}) ? {1, 1, 0, 0} : {1, 0, 0, 1}`
    Vector auxiliar(std::size_t gi)
    {
        return gluons[(gi + 1) % gluons.size()].momentum;
    }
};
