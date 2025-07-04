#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <algorithm>

// Параметры модели
const int NUM_PARTICLES = 2000;
const float WIDTH = 5.0f;
const float HEIGHT = 5.0f;
const float DEPTH = 5.0f;
const float PARTICLE_MASS = 0.05f;
const float K = 5000.0f;       // Модуль упругости
const float MU = 200.0f;       // Вязкость
const float YIELD_STRESS = 50.0f; // Предел текучести
const float REST_DENSITY = 1.0f;
const float H = 0.5f;          // Радиус сглаживания
const float DT = 0.001f;       // Шаг по времени
const float G = -9.81f;       // Гравитация
const vec3 GRAVITY(0.0f, G, 0.0f);

// Вектор в 3D пространстве
struct vec3 
{
    float x, y, z;
    
    vec3() : x(0), y(0), z(0) {}
    vec3(float x, float y, float z) : x(x), y(y), z(z) {}
    
    vec3 operator+(const vec3& other) const 
    {
        return vec3(x + other.x, y + other.y, z + other.z);
    }
    
    vec3 operator-(const vec3& other) const 
    {
        return vec3(x - other.x, y - other.y, z - other.z);
    }
    
    vec3 operator*(float scalar) const 
    {
        return vec3(x * scalar, y * scalar, z * scalar);
    }
    
    vec3& operator+=(const vec3& other) 
    {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
    
    float length() const 
    {
        return sqrt(x*x + y*y + z*z);
    }
    
    vec3 normalized() const 
    {
        float len = length();
        if (len > 0) return *this * (1.0f / len);
        return *this;
    }
};

// Структура для частицы
struct Particle 
{
    vec3 position;
    vec3 velocity;
    vec3 force;
    float density;
    float pressure;
    bool boundary;
    vec3 plastic_strain; // Для модели пластичности
};

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

// Расчет сил
void computeForces(std::vector<Particle>& particles) 
{
    for (auto& pi : particles) {
        if (pi.boundary) continue;
        
        pi.force = GRAVITY * PARTICLE_MASS; // Гравитация
        vec3 strain_rate(0, 0, 0);
        
        for (auto& pj : particles) {
            if (&pi == &pj) continue;
            
            vec3 r = pi.position - pj.position;
            float dist = r.length();
            
            if (dist < H && dist > 1e-6f) 
            {
                vec3 grad = gradW(r, H);
                
                // Сила давления
                float pressure_force = -PARTICLE_MASS * PARTICLE_MASS * 
                                      (pi.pressure / (pi.density * pi.density) + 
                                      (pj.pressure / (pj.density * pj.density));
                pi.force += grad * pressure_force;
                
                // Сила вязкости
                vec3 v_ij = pi.velocity - pj.velocity;
                pi.force += v_ij * (MU * PARTICLE_MASS * PARTICLE_MASS / 
                                   (pi.density * pj.density) * 
                                   (-45.0f / (M_PI * H * H * H * H * H)) * (H - dist));
                
                // Для модели пластичности
                strain_rate += (v_ij * PARTICLE_MASS / pj.density) * grad;
            }
        }
        
        // Применяем модель пластичности
        applyPlasticity(pi, strain_rate, DT);
    }
}

// Интегрирование (метод Верле)
void integrate(std::vector<Particle>& particles) 
{
    for (auto& p : particles) 
    {
        if (p.boundary) continue;
        
        // Сохраняем предыдущую скорость
        vec3 old_velocity = p.velocity;
        
        // Обновляем скорость
        p.velocity += p.force * (DT / PARTICLE_MASS);
        
        // Обновляем позицию
        p.position += p.velocity * DT;
        
        // Граничные условия
        if (p.position.x < 0.0f) {
            p.position.x = 0.0f;
            p.velocity.x *= -0.5f;
        }
        if (p.position.x > WIDTH) {
            p.position.x = WIDTH;
            p.velocity.x *= -0.5f;
        }
        if (p.position.y < 0.0f) {
            p.position.y = 0.0f;
            p.velocity.y *= -0.5f;
        }
        if (p.position.y > HEIGHT) {
            p.position.y = HEIGHT;
            p.velocity.y *= -0.5f;
        }
        if (p.position.z < 0.0f) {
            p.position.z = 0.0f;
            p.velocity.z *= -0.5f;
        }
        if (p.position.z > DEPTH) {
            p.position.z = DEPTH;
            p.velocity.z *= -0.5f;
        }
    }
}

// Сохранение в файл для визуализации
void saveToFile(const std::vector<Particle>& particles, int step) 
{
    std::ofstream file("frame_" + std::to_string(step) + ".xyz");
    for (const auto& p : particles) 
    {
        file << p.position.x << " " << p.position.y << " " << p.position.z << "\n";
    }
    file.close();
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