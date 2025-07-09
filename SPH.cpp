#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <algorithm>


int main() 
{
    //initialization particles
    std::vector<Particle> particles = initParticles();
    
    for (int step = 0; step < 1000; step++) 
    {
        computeDensityPressure(particles);
        computeForces(particles);
        integrate(particles);
        
        if (step % 10 == 0) 
        {
            saveToFile(particles, step/10);
            std::cout << "Step " << step << " completed." << std::endl;
        }
    }
    
    return 0;
}