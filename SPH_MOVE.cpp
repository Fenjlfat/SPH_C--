class SPH_MOVE
{
private:
    /* data */
public:
    SPH_MOVE(/* args */);
    ~SPH_MOVE();
    void MOVE();
};

SPH_MOVE::SPH_MOVE(/* args */)
{
}

SPH_MOVE::~SPH_MOVE()
{
}


void SPH_MOVE::MOVE()
{
    
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
        if (p.position.x < 0.0f) 
        {
            p.position.x = 0.0f;
            p.velocity.x *= -0.5f;
        }
        if (p.position.x > WIDTH) 
        {
            p.position.x = WIDTH;
            p.velocity.x *= -0.5f;
        }
        if (p.position.y < 0.0f) 
        {
            p.position.y = 0.0f;
            p.velocity.y *= -0.5f;
        }
        if (p.position.y > HEIGHT) 
        {
            p.position.y = HEIGHT;
            p.velocity.y *= -0.5f;
        }
        if (p.position.z < 0.0f) 
        {
            p.position.z = 0.0f;
            p.velocity.z *= -0.5f;
        }
        if (p.position.z > DEPTH) 
        {
            p.position.z = DEPTH;
            p.velocity.z *= -0.5f;
        }
    }
}


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
