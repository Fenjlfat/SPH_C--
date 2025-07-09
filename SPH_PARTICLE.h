//параметры частицы
struct SPH_Particle
{
    //kordinate
    double X = 0.0;
    double Y = 0.0;
    double Z = 0.0;
    //velocity
    double VX = 0.0;
    double VY = 0.0;
    double VZ = 0.0;
    //stress
    double SXX = 0.0;
    double SYY = 0.0;
    double SZZ = 0.0;
    double SXY = 0.0;
    double SXZ = 0.0;
    double SYZ = 0.0;

    
    /* data */
};
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