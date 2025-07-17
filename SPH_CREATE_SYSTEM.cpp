class SPH_CREATE_SYSTEM
{
private:
    /* data */
public:
    SPH_CREATE_SYSTEM(/* args */);
    ~SPH_CREATE_SYSTEM();
};

SPH_CREATE_SYSTEM::SPH_CREATE_SYSTEM(/* args */)     //konstruktor
{
}

SPH_CREATE_SYSTEM::~SPH_CREATE_SYSTEM()      //destruktor
{
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
                particles[i].position = vec3(x + dis(gen)*0.02f, y + dis(gen)*0.02f, z + dis(gen)*0.02f);
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