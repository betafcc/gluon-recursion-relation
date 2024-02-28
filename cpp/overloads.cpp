#include <functional>
#include <complex>

using std::function;

template <typename A, typename B, typename C>
function<C(A)> operator<<(function<B(A)> f, function<C(B)> g);

template <typename A, typename B, typename C>
function<C(A)> operator>>(function<C(B)> g, function<B(A)> f);

template <typename A, typename B, typename C>
function<C(A)> operator<<(function<B(A)> f, function<C(B)> g)
{
    return [f, g](A a) -> C
    { return g(f(a)); };
}

template <typename A, typename B, typename C>
function<C(A)> operator>>(function<C(B)> g, function<B(A)> f)
{
    return [f, g](A a) -> C
    { return g(f(a)); };
}

// Overload << operator to print std::vector<int>
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec)
{
    os << "⟨ ";
    for (size_t i = 0; i < vec.size(); ++i)
    {
        os << vec[i];
        if (i < vec.size() - 1)
            os << " ";
    }
    os << " ⟩";
    return os;
}
