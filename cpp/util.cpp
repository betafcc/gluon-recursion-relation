#include <vector>
#include <functional>
#include <stdexcept>

using std::function;
using std::invalid_argument;
using std::min;
using std::next;
using std::vector;

auto range(int start, int end, int step = 1) -> vector<int>
{
    vector<int> sequence;

    for (int i = start; (step > 0) ? i < end : i > end; i += step)
        sequence.push_back(i);

    return sequence;
}

template <typename A, typename B>
auto map(function<B(A)> f, const vector<A> &vec) -> vector<B>
{
    vector<B> result;
    result.reserve(vec.size());

    for (const A &item : vec)
        result.push_back(f(item));

    return result;
}

template <typename A, typename B>
auto map(function<B(A)> f) -> function<vector<B>(const vector<A> &)>
{
    return [f](const vector<A> &vec) -> vector<B>
    {
        return map(f, vec);
    };
}

template <typename A, typename B, typename C>
auto zip_with(function<C(A, B)> f, const vector<A> &va, const vector<B> &vb) -> vector<C>
{
    vector<C> result;
    size_t minSize = min(va.size(), vb.size());
    result.reserve(minSize);

    for (size_t i = 0; i < minSize; ++i)
        result.push_back(f(va[i], vb[i]));

    return result;
}

template <typename A, typename B, typename C>
auto zip_with(function<C(A, B)> f, const vector<A> &va) -> function<vector<C>(const vector<B> &)>
{
    return [&f, &va](const vector<B> &vb) -> vector<C>
    {
        return zip_with(f, va, vb);
    };
}

template <typename A, typename B, typename C>
auto zip_with(function<C(A, B)> f) -> function<function<vector<C>(const vector<B> &)>(const vector<A> &)>
{
    return [&f](const vector<A> &va) -> function<vector<C>(const vector<B> &)>
    {
        return zip_with<A, B, C>(f, va);
    };
}

template <typename A, typename B>
auto reduce(function<A(A, B)> f, const vector<B> &vec, A initial) -> A
{
    A acc = initial;
    for (const B &item : vec)
        acc = f(acc, item);

    return acc;
}

template <typename A, typename B>
auto reduce(function<A(A, B)> f, const vector<B> &vec) -> A
{
    if (vec.empty())
        return A{};
    // throw invalid_argument("Vector empty");

    A acc = vec.front();
    for (auto it = next(vec.begin()); it != vec.end(); ++it)
        acc = f(acc, *it);

    return acc;
}

template <typename A, typename B>
class ReduceFunctor
{
public:
    function<A(A, B)> f;

    ReduceFunctor(function<A(A, B)> f) : f(std::move(f)) {}

    A operator()(const vector<B> &vec, A initial) const
    {
        return reduce(f, vec, initial);
    }

    A operator()(const vector<B> &vec) const
    {
        return reduce(f, vec);
    }
};

template <typename A, typename B>
auto reduce(function<A(A, B)> f) -> ReduceFunctor<A, B>
{
    return ReduceFunctor<A, B>(std::move(f));
}