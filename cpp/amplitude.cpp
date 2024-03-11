#include "vec.cpp"
#include <complex>
#include <iostream>
#include <vector>

using Number = std::complex<double>;
using Vector = std::vector<Number>;
using Matrix = std::vector<Vector>;

enum Helicity {
    Minus = -1,
    Plus = 1
};

struct Gluon {
    Helicity helicity;
    Vector momentum;
};

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

// Overload to print complex numbers
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

// Overload to print vectors
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

Matrix dot(Matrix A, Matrix B)
{
    int m = A.size(); // Number of rows in A
    int n = B[0].size(); // Number of columns in B
    int p = B.size(); // Number of rows in B

    // Resultant matrix of size m x n
    Matrix result(m, Vector(n, 0));

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < p; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

// Minkowski Inner Product
Number mp(Vector a, Vector b)
{
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
}

Number mp2(Vector a)
{
    return mp(a, a);
}

// ⟨k | - definition
Matrix langle(Vector k)
{
    double num1 = sqrt((abs(k[0].real() + k[3].real())) / (k[0].real() + k[3].real()));
    Matrix kmod1 = { { Number(0, 0), Number(0, 0), Number(-(num1 * k[1].real()), -(num1 * k[2].real())), Number(num1 * (k[0].real() + k[3].real()), 0) } };
    return kmod1;
}

// | k⟩ - definition
Matrix rangle(Vector k)
{
    double num2 = sqrt((abs(k[0].real() + k[3].real())) / (k[0].real() + k[3].real()));
    Matrix kmod2 = { { Number(0, 0) },
                     { Number(0, 0) },
                     { Number(num2 * (k[0].real() + k[3].real()), 0) },
                     { Number(num2 * k[1].real(), num2 * k[2].real()) } };
    return kmod2;
}

// [k | - definition
Matrix lbox(Vector k)
{
    double num3 = sqrt(1 / (abs(k[0].real() + k[3].real())));
    Matrix kmod3 = { { Number(num3 * (k[0].real() + k[3].real()), 0), Number(num3 * k[1].real(), -(num3 * k[2].real())), Number(0, 0), Number(0, 0) } };
    return kmod3;
}

// | k] - definition
Matrix rbox(Vector k)
{
    double num4 = sqrt(1 / (abs(k[0].real() + k[3].real())));
    Matrix kmod4 = { { Number(-num4 * k[1].real(), num4 * k[2].real()) },
                     { Number(num4 * (k[0].real() + k[3].real()), 0) },
                     { Number(0, 0) },
                     { Number(0, 0) } };
    return kmod4;
}

// ⟨q | k⟩ - definition needs to be modified
Matrix braket(Vector q, Vector k)
{
    return dot(langle(q), rangle(k));
}

// [q | k] - definition needs to be modified
Matrix box(Vector q, Vector k)
{
    return dot(lbox(q), rbox(k));
}

// [q | g | k⟩ - definition needs to be modified
Matrix numerator1(Vector q, Matrix g, Vector k)
{
    return dot(lbox(q), dot(g, rangle(k)));
}

// [k | g | q⟩ - definition needs to be modified
Matrix numerator2(Vector q, Matrix g, Vector k)
{
    return dot(lbox(k), dot(g, rangle(q)));
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

            auto sb_acc = Vector(4, 0);
            for (auto m = 0; m < (n - 1); ++m)
                sb_acc = sb_acc + sb(std::vector<std::size_t>(gis.begin(), gis.begin() + m + 1), std::vector<std::size_t>(gis.begin() + m + 1, gis.begin() + n));

            auto cb_acc = Vector(4, 0);
            for (auto m = 0; m < (n - 2); ++m)
                for (auto k = m + 1; k < (n - 1); ++k)
                    cb_acc = cb_acc + cb(std::vector<std::size_t>(gis.begin(), gis.begin() + m + 1), std::vector<std::size_t>(gis.begin() + m + 1, gis.begin() + k + 1), std::vector<std::size_t>(gis.begin() + k + 1, gis.begin() + n));

            return factor * (sb_acc + cb_acc);
        }
    }

    // berends 2.10
    // All vectors `gluons[i].momentum` for i ∈ [start, end]	//Query: should the end be open or closed????
    Vector kappa(std::size_t start, std::size_t end) const
    {
        Vector acc = gluons[start].momentum;

        for (auto i = start; i <= end; ++i)
            acc = acc + gluons[i].momentum;

        return acc;
    } /// this part is most probably not needed, since we are using gis as input for kappa everywhere!! but need to confirm.
    // takes a vector of indexes,
    // it's not clear how kappa should behave if indices out of order!!!

    Vector kappa(std::vector<std::size_t> gis) const
    {
        Vector acc = gluons[gis[0]].momentum;

        for (auto i = 1; i <= gis.size(); ++i)
            acc = acc + gluons[gis[i]].momentum;

        return acc;
    }

    friend std::ostream& operator<<(std::ostream& os, const Process& p)
    {
        auto indent = "    ";

        os << "Process: \n";

        os << indent << "g \t λ  p" << '\n';
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
    Number sb(
        std::vector<std::size_t> xs,
        std::vector<std::size_t> ys,
        std::size_t xi
    )
    {
        return dot(2 * kappa(ys), current(xs) * current(ys, xi))
            - dot(2 * kappa(xs), current(ys) * current(xs, xi))
            + (kappa(xs) - kappa(ys)).at(xi)
            * dot(current(xs), current(ys));
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

    Number cb(
        std::vector<std::size_t> xs,
        std::vector<std::size_t> ys,
        std::vector<std::size_t> zs,
        std::size_t xi
    )
    {

        return dot(current(xs), current(zs) * current(ys, xi) - current(ys) * current(zs, xi))
            - dot(current(zs), current(ys) * current(xs, xi) - current(xs) * current(ys, xi));
    }

    Vector cb(
        std::vector<std::size_t> xs,
        std::vector<std::size_t> ys,
        std::vector<std::size_t> zs
    )
    {
        return {
            cb(xs, ys, zs, 0),
            cb(xs, ys, zs, 1),
            cb(xs, ys, zs, 2),
            cb(xs, ys, zs, 3),
        };
    }

    Number polarization(std::size_t gi, std::size_t xi)
    {
        auto gm = Gamma[xi];
        auto q = auxiliary(gi);
        auto k = gluons[gi].momentum;

        switch (gluons[gi].helicity) {
        case Helicity::Plus:
            return (1 / sqrt(2)) * (numerator2(q, gm, k) / braket(q, k));
        case Helicity::Minus:
            return (-1 / sqrt(2)) * (numerator1(q, gm, k) / box(q, k));
        }
        return 0;
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

    Vector auxiliary(std::size_t gi)
    {
        return gluons[(gi + 1) % gluons.size()].momentum;
    }
};
