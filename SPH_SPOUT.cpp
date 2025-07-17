// Сохранение в файл для визуализации

void saveToFile(const std::vector<Particle> &particles, int step) 
{
    std::ofstream file("frame_" + std::to_string(step) + ".xyz");
    for (const auto& p : particles) 
    {
        file << p.position.x << " " << p.position.y << " " << p.position.z << "\n";
    }
    file.close();
}
