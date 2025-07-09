
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