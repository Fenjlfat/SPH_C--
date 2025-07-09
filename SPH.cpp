#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <algorithm>



// Ядро сглаживания (кубический сплайн в 3D)
float W(float r, float h) 
{
    float q = r / h;
    if (q <= 1.0f) 
    {
        return (1.0f / (M_PI * h * h * h)) * (1.0f - 1.5f*q*q + 0.75f*q*q*q);
    } 
    else if (q <= 2.0f) 
    {
        return (1.0f / (4.0f * M_PI * h * h * h)) * (2.0f - q)*(2.0f - q)*(2.0f - q);
    } 
    else 
    {
        return 0.0f;
    }
}

// Градиент ядра сглаживания
vec3 gradW(const vec3& r_vec, float h) 
{
    float r = r_vec.length();
    float q = r / h;
    
    if (r < 1e-6f) return vec3(0, 0, 0);
    
    vec3 result;
    if (q <= 1.0f) 
    {
        float factor = (1.0f / (M_PI * h * h * h * h)) * (-3.0f*q + 2.25f*q*q);
        result = r_vec * (factor / r);
    } 
    else if (q <= 2.0f) 
    {
        float factor = (1.0f / (4.0f * M_PI * h * h * h * h)) * (-0.75f * (2.0f - q)*(2.0f - q));
        result = r_vec * (factor / r);
    } 
    else 
    {
        result = vec3(0, 0, 0);
    }
    
    return result;
}

// Инициализация частиц
std::vector<Particle> initParticles() 
{
    std::vector<Particle> particles(NUM_PARTICLES);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(0.0f, 1.0f);
    
    // Металлический блок в центре
    int i = 0;
    for (float y = HEIGHT/2 - 1.0f; y < HEIGHT/2 + 1.0f && i < NUM_PARTICLES; y += 0.15f) 
    {
        for (float x = WIDTH/2 - 1.0f; x < WIDTH/2 + 1.0f && i < NUM_PARTICLES; x += 0.15f) 
        {
            for (float z = DEPTH/2 - 1.0f; z < DEPTH/2 + 1.0f && i < NUM_PARTICLES; z += 0.15f) 
            {
                particles[i].position = vec3(
                    x + dis(gen)*0.02f,
                    y + dis(gen)*0.02f,
                    z + dis(gen)*0.02f
                );
                particles[i].velocity = vec3(0, 0, 0);
                particles[i].force = vec3(0, 0, 0);
                particles[i].density = 0.0f;
                particles[i].pressure = 0.0f;
                particles[i].boundary = false;
                particles[i].plastic_strain = vec3(0, 0, 0);
                i++;
            }
        }
    }
    
    // Граничные частицы (стены и пол)
    for (; i < NUM_PARTICLES; i++) 
    {
        if (i % 4 == 0) 
        {
            particles[i].position = vec3(dis(gen) * WIDTH, 0.0f, dis(gen) * DEPTH);
        } 
        else if (i % 4 == 1) 
        {
            particles[i].position = vec3(dis(gen) * WIDTH, HEIGHT, dis(gen) * DEPTH);
        } 
        else if (i % 4 == 2) 
        {
            particles[i].position = vec3(0.0f, dis(gen) * HEIGHT, dis(gen) * DEPTH);
        } 
        else 
        {
            particles[i].position = vec3(WIDTH, dis(gen) * HEIGHT, dis(gen) * DEPTH);
        }
        
        particles[i].velocity = vec3(0, 0, 0);
        particles[i].force = vec3(0, 0, 0);
        particles[i].density = 0.0f;
        particles[i].pressure = 0.0f;
        particles[i].boundary = true;
        particles[i].plastic_strain = vec3(0, 0, 0);
    }
    
    return particles;
}

// Расчет плотности и давления
void computeDensityPressure(std::vector<Particle>& particles) 
{
    for (auto& pi : particles) 
    {
        pi.density = 0.0f;
        for (const auto& pj : particles) 
        {
            vec3 r = pi.position - pj.position;
            float dist = r.length();
            
            if (dist < H) 
            {
                pi.density += PARTICLE_MASS * W(dist, H);
            }
        }
        
        // Уравнение состояния (Tait equation для металлов)
        pi.pressure = K * (pi.density - REST_DENSITY);
    }
}

// Модель пластичности (von Mises yield criterion)
void applyPlasticity(Particle& pi, const vec3& strain_rate, float dt) 
{
    float strain_rate_mag = strain_rate.length();
    if (strain_rate_mag > 0) 
    {
        vec3 deviatoric = strain_rate - strain_rate * (1.0f/3.0f) * strain_rate_mag;
        float deviatoric_mag = deviatoric.length();
        
        float yield_factor = std::min(1.0f, YIELD_STRESS / (2.0f * MU * deviatoric_mag));
        pi.plastic_strain = pi.plastic_strain + (1.0f - yield_factor) * deviatoric * dt;
    }
}




int main() 
{
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