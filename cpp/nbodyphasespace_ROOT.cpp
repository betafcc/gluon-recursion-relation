#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <ctime>
#include "TVectorD.h"
#include "TLorentzVector.h"

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
TLorentzVector generate4momenta(double mass, double energy) {
    double cos_theta = rand_cosine();
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = 2 * M_PI * rand_uniform();

    double pmod = std::sqrt(energy * energy - mass * mass);
    double px = pmod * sin_theta * std::cos(phi);
    double py = pmod * sin_theta * std::sin(phi);
    double pz = pmod * cos_theta;

    return TLorentzVector(px, py, pz, energy);
}


void generatePhasespace(const double centerOfMassEnergy, const std::vector<double>& masses, const std::vector<double>& invariant_mass_values) {

    // Initialize particles in the center-of-mass frame
    TLorentzVector pA( 0.0, 0.0, std::sqrt(centerOfMassEnergy) / 2.0, std::sqrt(centerOfMassEnergy) / 2.0);
    TLorentzVector pB( 0.0, 0.0, -std::sqrt(centerOfMassEnergy) / 2.0, std::sqrt(centerOfMassEnergy) / 2.0);
    
    //"FinalStateParticles" stores p1, p2,...., pn
    int n = masses.size();
    std::vector<TLorentzVector> FinalStateParticles;

    // Calculating 4-momenta of particle 1 and 2 in CM12

    double E2 = (invariant_mass_values[1] * invariant_mass_values[1] + masses[1] * masses[1] - masses[0] * masses[0]) / (2 * invariant_mass_values[1]);
    double E1 = invariant_mass_values[1] - E2;

    TLorentzVector p2 = generate4momenta(masses[1], E2);
    TLorentzVector p1(-p2.Px(), -p2.Py(), -p2.Pz(), E1);
    
    FinalStateParticles.push_back(p1);
    FinalStateParticles.push_back(p2);
    

    //moving to more than 2 particles
    
    for(size_t i = 3; i < n + 1; ++i){
        
        for(size_t j = 0; j < i - 1; ++j){
            TLorentzVector k(0, 0, 0, 0);
            k = k + FinalStateParticles[j];
	}
        
        double Ek = (invariant_mass_values[i - 1] * invariant_mass_values[i - 1] + invariant_mass_values[i - 2] * invariant_mass_values[i - 2] - masses[i - 1] * masses[i - 1]) / (2 * invariant_mass_values[i - 1]);
        double E = invariant_mass_values[i - 1] - Ek;
        
        TLorentzVector p = generate4momenta(masses[i - 1], E);
        FinalStateParticles.push_back(p);
        
        TVector3 boostVector(-p.Px() / Ek, -p.Py() / Ek, -p.Pz() / Ek);
            
        // Perform Lorentz transformations
        for(size_t j = 0; j < i - 1; ++j ){
            FinalStateParticles[j].Boost(boostVector);
	}
    }
    std::cout << "Initial states:" << std::endl;
    std::cout << "(" << pA.E() << "," << pA.Px() << "," << pA.Py() << "," << pA.Pz() << ")" << std::endl;
    std::cout << "(" << pB.E() << "," << pB.Px() << "," << pB.Py() << "," << pB.Pz() << ")" << std::endl;
    
    // Print the final state 4 momenta
    for (int i = 0; i < FinalStateParticles.size(); ++i) {
        std::cout << "Particle " << i + 1 << ": ";
        std::cout<< "(" << FinalStateParticles[i].E() << "," << FinalStateParticles[i].Px() << "," << FinalStateParticles[i].Py() << "," << FinalStateParticles[i].Pz() << ")" << std::endl;
        //std::cout<<std::endl;
        std::cout << FinalStateParticles[i].M();
        std::cout << std::endl;
    }

    // Check energy-momentum conservations
    TLorentzVector TotalInitialMomentum = pA + pB ;
    TLorentzVector TotalFinalMomentum(0.0, 0.0, 0.0, 0.0);

    for (const TLorentzVector& particle : FinalStateParticles) {
        TotalFinalMomentum = TotalFinalMomentum + particle;
    }

    TLorentzVector TotalMomentum = TotalInitialMomentum - TotalFinalMomentum;
    std::cout << "Total momentum: ";
    std::cout<< "(" << TotalMomentum.E() << "," << TotalMomentum.Px() << "," << TotalMomentum.Py() << "," << TotalMomentum.Pz() << "," << ")" << std::endl;
    std::cout << std::endl;
    
}


int main() {
    
    double centerOfMassEnergy = 1000.0;  
    
    std::vector<double> masses = {0.0, 0.0, 0.0};
    //std::vector<double> masses = {0, 0};
    
    std::vector<double> invariant_mass_values = invariant_mass(centerOfMassEnergy, masses);

    
    generatePhasespace(centerOfMassEnergy, masses, invariant_mass_values);
    

    return 0;
}


