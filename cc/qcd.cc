#include "i.cpp"
#include <stdexcept>

#include "Vector.cpp"

typedef struct
{
    int helicity;
    Vector vector;
} VecInfo;

static const Matrix GammaU0 {
    { 1, 0, 0, 0 },
    { 0, 1, 0, 0 },
    { 0, 0, -1, 0 },
    { 0, 0, 0, -1 },
};

static const Matrix GammaU1 {
    { 0, 0, 0, 1 },
    { 0, 0, 1, 0 },
    { 0, -1, 0, 0 },
    { -1, 0, 0, 0 },
};

static const Matrix GammaU2 {
    { 0, 0, 0, i<double>(-1) },
    { 0, 0, i<double>(1), 0 },
    { 0, i<double>(1), 0, 0 },
    { i<double>(-1), 0, 0, 0 },
};

static const Matrix GammaU3 {
    { 0, 0, 1, 0 },
    { 0, 0, 0, -1 },
    { -1, 0, 0, 0 },
    { 0, 1, 0, 0 },
};

static const std::vector<Matrix> Gamma { GammaU0, GammaU1, GammaU2, GammaU3 };

// [q, gx, k⟩
auto fish(Vector q, Matrix gx, Vector k) -> Vector
{
    throw std::invalid_argument("Not Implemented");
}

// ⟨q | k⟩
auto braket(Vector q, Vector k) -> double
{
    throw std::invalid_argument("Not Implemented");
}

// [q | k]
auto box(Vector q, Vector k) -> double
{
    throw std::invalid_argument("Not Implemented");
}

auto gen_aux_vec(Vector xs) -> Vector
{
    if (is_linear_dependent(xs, { 1, 1, 0, 0 }))
        return { 1, 0, 1, 0 };
    return { 1, 1, 0, 0 };
}

// Gluon Polarization Vector
auto gpv(int i, int helicity, Vector k) -> Vector
{
    auto q = gen_aux_vec(k);
    auto gx = Gamma[i];

    if (helicity > 0)
        return (1 / sqrt(2)) * (fish(k, gx, q) * (1 / braket(q, k)));
    else
        return (-1 / sqrt(2)) * (fish(q, gx, k) * (1 / box(q, k)));
}

auto kappa(std::vector<VecInfo> vs) -> Vector
{
    auto acc = vs.front().vector;

    for (auto it = next(vs.begin()); it != vs.end(); ++it)
        acc = acc + it->vector;

    return acc;
}

auto current(int i, std::vector<VecInfo> vs) -> Vector;

// Minkowski Inner Product
auto mp(Vector a, Vector b) -> std::complex<double>
{
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
}

auto mp2(Vector a) -> std::complex<double>
{
    return mp(a, a);
}

std::complex<double> square_brackets(
    int i,
    std::vector<VecInfo> vsl,
    std::vector<VecInfo> vsr
)
{
    // clang-format off
    return 2.0 * mp(kappa(vsr), current(i, vsl)) * current(i, vsr)[i]
        - 2.0 * mp(kappa(vsl), current(i, vsr)) * current(i, vsl)[i]
        + (kappa(vsl) - kappa(vsr))[i] * mp(current(i, vsl), current(i, vsr));
    // clang-format on
}

std::complex<double> curly_brackets(
    int i, std::vector<VecInfo> vsl,
    std::vector<VecInfo> vsc,
    std::vector<VecInfo> vsr
)
{
    // clang-format off
    return mp(current(i, vsl), (current(i, vsr) * current(i, vsc)[i] - current(i, vsc) * current(i, vsr)[i])) -
           mp(current(i, vsr), (current(i, vsc) * current(i, vsl)[i] - current(i, vsl) * current(i, vsc)[i]));
    // clang-format on
}

auto current(int i, std::vector<VecInfo> vs) -> Vector
{
    if (vs.size() == 1)
        return gpv(i, vs[0].helicity, vs[0].vector);
    else if (vs.size() == 2) {
        Vector result;
        result.reserve(vs[0].vector.size());

        for (auto j = 0; j < 4; j++)
            result[j] = square_brackets(j, { vs[0] }, { vs[1] }) * (1.0 / mp2(vs[0].vector + vs[1].vector));

        return result;
    } else {
        Vector result;
        result.reserve(vs[0].vector.size());

        for (auto j = 0; j < 4; j++) {

            auto factor = 1.0 / mp2(kappa(vs));
            auto n = vs.size();

            std::complex<double> sb_acc = 0;
            for (auto m = 0; m < (n - 1); m++) {
                std::vector<VecInfo> left(vs.begin() + 0, vs.begin() + m);
                std::vector<VecInfo> right(vs.begin() + m + 1, vs.begin() + n);
                sb_acc = sb_acc + square_brackets(j, left, right);
            }

            std::complex<double> cb_acc = 0;
            for (auto m = 0; m < (n - 2); m++)
                for (auto k = m + 1; k < (n - 1); k++) {
                    std::vector<VecInfo> left(vs.begin() + 0, vs.begin() + m);
                    std::vector<VecInfo> center(vs.begin() + m + 1, vs.begin() + k);
                    std::vector<VecInfo> right(vs.begin() + k + 1, vs.begin() + n);
                    cb_acc = cb_acc + curly_brackets(j, left, center, right);
                }

            result[j] = factor * (sb_acc + cb_acc);
        }

        return result;
    }
}
