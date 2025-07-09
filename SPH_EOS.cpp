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