#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <ctime>

class FourVector {
 
private:
    double t, x, y, z;

public:
    // Constructors
    FourVector() : t(0), x(0), y(0), z(0) {}
    FourVector(double t_, double x_, double y_, double z_) : t(t_), x(x_), y(y_), z(z_) {}
    FourVector(const FourVector& other) : t(other.t), x(other.x), y(other.y), z(other.z) {}

    // Get and Set functions
    double getT() const { return t; }
    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }

 
    // Vector operations
    FourVector operator+(const FourVector& other) const { return FourVector(t + other.t, x + other.x, y + other.y, z + other.z);}

    FourVector operator-(const FourVector& other) const { return FourVector(t - other.t, x - other.x, y - other.y, z - other.z);}

    FourVector operator*(double scalar) const { return FourVector(t * scalar, x * scalar, y * scalar, z * scalar);}

    double operator*(const FourVector& other) const { return t * other.t - x * other.x - y * other.y - z * other.z;}

    double norm() const { return std::sqrt(t * t - x * x - y * y - z * z);}
    
    void print() const { std::cout << "(" << t << ", " << x << ", " << y << ", " << z << ")";}
};

class ThreeVector {
 
private:
    double x, y, z;

public:
    // Constructors
    ThreeVector() : x(0), y(0), z(0) {}
    ThreeVector(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    ThreeVector(const ThreeVector& other) : x(other.x), y(other.y), z(other.z) {}

    // Get and Set functions
    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }
    
    double modulus() const { return std::sqrt( x * x + y * y + z * z);}
    
    void print3() const { std::cout << "(" << x << ", " << y << ", " << z << ")";}
    
};

// Define a function to generate a random double in the range [0, 1)
double rand_uniform() {
    static std::mt19937_64 generator(std::random_device{}());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

// Function for a random cosine in the range [-1, 1]
double rand_cosine() {
    return 2.0 * rand_uniform() - 1.0;
}

// Function to calculate the invariant mass in a 2-body collision
std::vector<double> invariant_mass(const double centerOfMassEnergy, const std::vector<double>& masses) {
    int n = masses.size();
    std::vector<double> invariant_m(n, 0.0);
    double sum = 0.0;
    const double mass_total = std::accumulate(masses.begin(), masses.end(), 0.0);
    
    // Vector of random numbers and sort it
    std::vector<double> random_numbers;
    for (size_t i = 0; i < masses.size() - 1; ++i) {
        random_numbers.push_back(rand_uniform());
    }
    std::sort(random_numbers.begin(), random_numbers.end());
    
    for (size_t i = 0; i < masses.size() - 1; ++i) {
        sum += masses[i];
        invariant_m[i] = sum + random_numbers[i] * (std::sqrt(centerOfMassEnergy) - mass_total);
    }
    invariant_m[n - 1] = std::sqrt(centerOfMassEnergy);
    return invariant_m;
}


// Function to generate random particle momenta
FourVector generate4momenta(double mass, double energy) {
    double cos_theta = rand_cosine();
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = 2 * M_PI * rand_uniform();

    double pmod = std::sqrt(energy * energy - mass * mass);
    double px = pmod * sin_theta * std::cos(phi);
    double py = pmod * sin_theta * std::sin(phi);
    double pz = pmod * cos_theta;

    return FourVector(energy, px, py, pz);
}

// Function to perform Lorentz transformations
FourVector lorentzTransform(const FourVector& vector, const ThreeVector& boostVector) {
    double beta = boostVector.modulus();
    double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    double factor1 = (gamma - 1.0) / (beta * beta);
    double factor2 = vector.getX() * boostVector.getX() + vector.getY() * boostVector.getY() + vector.getZ() * boostVector.getZ();

    FourVector transformed(
        gamma * (vector.getT() + factor2),
        vector.getX() + boostVector.getX() * factor1 * factor2 + gamma * boostVector.getX() * vector.getT(),
        vector.getY() + boostVector.getY() * factor1 * factor2 + gamma * boostVector.getY() * vector.getT(),
        vector.getZ() + boostVector.getZ() * factor1 * factor2 + gamma * boostVector.getZ() * vector.getT()
    );

    return transformed;
}

void generatePhasespace(const double centerOfMassEnergy, const std::vector<double>& masses, const std::vector<double>& invariant_mass_values) {

    // Initialize particles in the center-of-mass frame
    FourVector pA(std::sqrt(centerOfMassEnergy) / 2.0, 0.0, 0.0, std::sqrt(centerOfMassEnergy) / 2.0);
    FourVector pB(std::sqrt(centerOfMassEnergy) / 2.0, 0.0, 0.0, -std::sqrt(centerOfMassEnergy) / 2.0);
    
    //"FinalStateParticles" stores p1, p2,...., pn
    int n = masses.size();
    std::vector<FourVector> FinalStateParticles;

    // Calculating 4-momenta of particle 1 and 2 in CM12

    double E2 = (invariant_mass_values[1] * invariant_mass_values[1] + masses[1] * masses[1] - masses[0] * masses[0]) / (2 * invariant_mass_values[1]);
    double E1 = invariant_mass_values[1] - E2;

    FourVector p2 = generate4momenta(masses[1], E2);
    FourVector p1(E1, -p2.getX(), -p2.getY(), -p2.getZ());
    
    FinalStateParticles.push_back(p1);
    FinalStateParticles.push_back(p2);
    

    //moving to more than 2 particles
    
    for(size_t i = 3; i < n + 1; ++i){
        
        for(size_t j = 0; j < i - 1; ++j){
            FourVector k(0, 0, 0, 0);
            k = k + FinalStateParticles[j];
	}
        
        double Ek = (invariant_mass_values[i - 1] * invariant_mass_values[i - 1] + invariant_mass_values[i - 2] * invariant_mass_values[i - 2] - masses[i - 1] * masses[i - 1]) / (2 * invariant_mass_values[i - 1]);
        double E = invariant_mass_values[i - 1] - Ek;
        
        FourVector p = generate4momenta(masses[i - 1], E);
        FinalStateParticles.push_back(p);
        
        ThreeVector boostVector = {-p.getX() / Ek, -p.getY() / Ek, -p.getZ() / Ek};
            
        // Perform Lorentz transformations
        for(size_t j = 0; j < i - 1; ++j ){
            FinalStateParticles[j] = lorentzTransform(FinalStateParticles[j], boostVector);
	}
    }
    
    // Print the final state 4 momenta
    for (int i = 0; i < FinalStateParticles.size(); ++i) {
        std::cout << "Particle " << i + 1 << ": ";
        FinalStateParticles[i].print();
        std::cout<<std::endl;
        //std::cout << FinalStateParticles[i].norm();
        std::cout << std::endl;
    }

    // Check energy-momentum conservations
    FourVector TotalInitialMomentum = pA + pB ;
    FourVector TotalFinalMomentum(0.0, 0.0, 0.0, 0.0);

    for (const FourVector& particle : FinalStateParticles) {
        TotalFinalMomentum = TotalFinalMomentum + particle;
    }

    FourVector TotalMomentum = TotalInitialMomentum - TotalFinalMomentum;
    std::cout << "Total momentum: ";
    TotalMomentum.print();
    std::cout << std::endl;
    
}


int main() {
    
    double centerOfMassEnergy = 5000.0;  
    
    int num_events = 10;
    
    std::vector<double> masses = {2.0, 5.0, 4.5, 7.5, 4.0};
    if(masses.size() > 15){
    std::cout << "Too many final states." << std::endl;
    }
    
    std::vector<double> invariant_mass_values = invariant_mass(centerOfMassEnergy, masses);

    for(int q = 0; q < num_events; q++){
    std::cout << "Event number:" << q + 1 <<std::endl;
    generatePhasespace(centerOfMassEnergy, masses, invariant_mass_values);
    }
    
    

    return 0;
}


