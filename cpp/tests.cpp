#include "process.cpp"
#include <iostream>
#include <random>
#include <vector>

class Random {
private:
    std::mt19937 gen;

public:
    explicit Random(std::mt19937::result_type seed = std::random_device {}()) :
        gen(seed) { }

    template <typename T>
    T randint(T min, T max)
    {
        std::uniform_int_distribution<T> distr(min, max);
        return distr(gen);
    }

    double random()
    {
        std::uniform_real_distribution<> distr(0, 1);
        return distr(gen);
    }

    template <typename T>
    T random(T min, T max)
    {
        std::uniform_real_distribution<T> distr(min, max);
        return distr(gen);
    }

    template <typename T>
    T choice(const std::vector<T>& vec)
    {
        std::uniform_int_distribution<std::size_t> distr(0, vec.size() - 1);
        return vec[distr(gen)];
    }

    template <typename T>
    std::vector<T> vector(std::size_t size)
    {
        std::vector<T> vec(size);
        for (std::size_t i = 0; i < size; ++i)
            vec[i] = random();
        return vec;
    }

    Helicity helicity()
    {
        return choice(std::vector<Helicity>({ Helicity::Minus, Helicity::Plus }));
    }

    Vector momentum()
    {
        return vector<Number>(4);
    }

    // Vector momentum(energy, mass)
    // {

    //     double cos_theta = rand_cosine();
    //     double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    //     double phi = 2 * M_PI * rand_uniform();

    //     double pmod = std::sqrt(energy * energy - mass * mass);
    //     double px = pmod * sin_theta * std::cos(phi);
    //     double py = pmod * sin_theta * std::sin(phi);
    //     double pz = pmod * cos_theta;

    //     return { energy, px, py, pz };
    // }

    Gluon gluon()
    {
        return { helicity(), momentum() };
    }

    Process process()
    {
        std::vector<Gluon> gluons;
        for (std::size_t i = 0; i < 4; ++i)
            gluons.push_back(gluon());

        return Process(gluons);
    }
};

void test(const char* message, bool condition)
{
    if (condition)
        std::cout << message << " ✅\n";
    else
        std::cerr << message << " ❌\n";
}

int main()
{
    std::vector<Gluon> gluons = {
        { Plus, { 1, 2, 3, 4 } },
        { Minus, { 5, 6, 7, 8 } },
        { Plus, { 2, 3, 4, 5 } },

    };

    Process process(gluons);

    auto xs = Vector({ 1, 2, 3, 4 });

    std::cout << xs << '\n';

    std::cout << Vector(xs.begin() + 0, xs.begin() + 1) << '\n';

    std::cout << process << '\n';
    std::cout << process.current({ 1, 0 }) << '\n';
    std::cout << process.current({ 0, 1, 2 });
    std::cout << '\n';

    // TODO:
    // berends (2.5)
    // (K₁ + K₂) · J(1, 2) = 0

    // TODO:
    // berends (2.6)
    // J_ξ(1, 2) = -J_ξ(2, 1)

    // test(
    //     "berends (2.15) -> J(3, 2, 1) = J(1, 2, 3) ",
    //     process.current(2, 1, 0) == process.current(0, 1, 2)
    // );

    // test(
    //     "berends (2.16) -> J(1, 2, 3) + J(2, 3, 1)  + J(3, 1, 2) = 0",
    //     process.current(0, 1, 2) + process.current(1, 2, 0) == process.current(2, 0, 1)
    // );

    // test(
    //     "berends (2.17) -> (K₁ + K₂ + K₃) · J(1, 2, 3) = 0",
    //     dot(gluons[0].momentum + gluons[1].momentum + gluons[2].momentum, process.current(0, 1, 2))
    //         == 0.0
    // );

    // std::cout << Number(0, -1);

    std::cout << Number(0, -1) << " == -i\n";
    std::cout << Number(0, 1) << " == i\n";
    std::cout << Number(0, -2) << " == -2i\n";
    std::cout << Number(0, 2) << " == 2i\n";

    std::cout << Number(0, 0) << " == 0\n";
    std::cout << Number(1, 0) << " == 1\n";
    std::cout << Number(-1, 0) << " == -1\n";
    std::cout << Number(-1, -1) << " == -1-i\n";

    std::cout
        << "\nDone 💅\n";
}
