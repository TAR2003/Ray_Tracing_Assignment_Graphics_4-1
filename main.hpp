#ifndef RAY_TRACING_MAIN_H
#define RAY_TRACING_MAIN_H

#include <vector>
#include <cmath>
#include <iostream>
#include "bitmap_image.hpp"

// OpenGL / GLUT Headers
#ifdef __APPLE__
#include <GLUT/glut.h> // Use GLUT framework on macOS
#else
#include <GL/glut.h> // Use standard GLUT location on Linux/Windows
#endif

// Forward declarations
class Object;
class PointLight;
class SpotLight;

// Vector3D class for 3D points and directions
class Vector3D {
public:
    double x, y, z;
    
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}
    
    Vector3D operator+(const Vector3D& v) const {
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }
    
    Vector3D operator-(const Vector3D& v) const {
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }
    
    Vector3D operator*(double scalar) const {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }
    
    Vector3D operator/(double scalar) const {
        return Vector3D(x / scalar, y / scalar, z / scalar);
    }
    
    double dot(const Vector3D& v) const {
        return x * v.x + y * v.y + z * v.z;
    }
    
    Vector3D cross(const Vector3D& v) const {
        return Vector3D(y * v.z - z * v.y,
                       z * v.x - x * v.z,
                       x * v.y - y * v.x);
    }
    
    double length() const {
        return sqrt(x * x + y * y + z * z);
    }
    
    Vector3D normalize() const {
        double len = length();
        if (len == 0) return *this;
        return *this / len;
    }
};

// Color class for RGB colors
class Color {
public:
    double r, g, b;
    
    Color() : r(0), g(0), b(0) {}
    Color(double r, double g, double b) : r(r), g(g), b(b) {}
    
    Color operator+(const Color& c) const {
        return Color(r + c.r, g + c.g, b + c.b);
    }
    
    Color operator*(double scalar) const {
        return Color(r * scalar, g * scalar, b * scalar);
    }
    
    Color operator*(const Color& c) const {
        return Color(r * c.r, g * c.g, b * c.b);
    }
    
    void clamp() {
        r = std::max(0.0, std::min(1.0, r));
        g = std::max(0.0, std::min(1.0, g));
        b = std::max(0.0, std::min(1.0, b));
    }
};

// Ray class for ray tracing
class Ray {
public:
    Vector3D start;
    Vector3D dir;
    
    Ray() {}
    Ray(const Vector3D& start, const Vector3D& dir) : start(start), dir(dir.normalize()) {}
};

// Base Object class
class Object {
public:
    Vector3D reference_point;
    double height, width, length;
    Color color;
    double coEfficients[4]; // ambient, diffuse, specular, reflection
    int shine;
    
    Object() : shine(0) {
        coEfficients[0] = coEfficients[1] = coEfficients[2] = coEfficients[3] = 0;
    }
    
    virtual void draw() = 0;
    virtual double intersect(Ray* r, Color* color, int level) = 0;
    virtual Vector3D getNormal(const Vector3D& point) const = 0;
    virtual Color getColorAt(const Vector3D& point) const = 0;
    
    void setColor(double r, double g, double b) {
        color = Color(r, g, b);
    }
    
    void setShine(int shine) {
        this->shine = shine;
    }
    
    void setCoEfficients(double ambient, double diffuse, double specular, double reflection) {
        coEfficients[0] = ambient;
        coEfficients[1] = diffuse;
        coEfficients[2] = specular;
        coEfficients[3] = reflection;
    }
};

// Sphere class
class Sphere : public Object {
public:
    Sphere(const Vector3D& center, double radius) {
        reference_point = center;
        length = radius;
    }
    
    void draw() override;
    double intersect(Ray* r, Color* color, int level) override;
    Vector3D getNormal(const Vector3D& point) const override;
    Color getColorAt(const Vector3D& point) const override;
};

// Triangle class
class Triangle : public Object {
public:
    Vector3D a, b, c;
    
    Triangle(const Vector3D& a, const Vector3D& b, const Vector3D& c) : a(a), b(b), c(c) {}
    
    void draw() override;
    double intersect(Ray* r, Color* color, int level) override;
    Vector3D getNormal(const Vector3D& point) const override;
    Color getColorAt(const Vector3D& point) const override;
};

// Floor class
class Floor : public Object {
public:
    double floorWidth, tileWidth;
    
    Floor(double floorWidth, double tileWidth);
    
    void draw() override;
    double intersect(Ray* r, Color* color, int level) override;
    Vector3D getNormal(const Vector3D& point) const override;
    Color getColorAt(const Vector3D& point) const override;
};

// General Quadric Surface class
class GeneralQuadric : public Object {
public:
    double coefficients[10]; // A-J coefficients
    
    GeneralQuadric(double A, double B, double C, double D, double E,
                   double F, double G, double H, double I, double J);
    
    void draw() override;
    double intersect(Ray* r, Color* color, int level) override;
    Vector3D getNormal(const Vector3D& point) const override;
    Color getColorAt(const Vector3D& point) const override;
    
private:
    bool isWithinReferenceCube(const Vector3D& point) const;
};

// Light classes
class PointLight {
public:
    Vector3D light_pos;
    Color color;
    
    PointLight(const Vector3D& pos, const Color& color) : light_pos(pos), color(color) {}
    
    void draw();
};

class SpotLight {
public:
    PointLight pointLight;
    Vector3D light_dir;
    double cutoff_angle;
    
    SpotLight(const Vector3D& pos, const Color& color, const Vector3D& dir, double angle)
        : pointLight(pos, color), light_dir(dir.normalize()), cutoff_angle(angle) {}
    
    void draw();
};

// Global variables
extern std::vector<Object*> objects;
extern std::vector<PointLight> pointLights;
extern std::vector<SpotLight> spotLights;
extern int recursionLevel;
extern int imageWidth, imageHeight;
extern int captureCount;

// Function declarations
void initGL();
void display();
void reshapeListener(GLsizei width, GLsizei height);
void keyboardListener(unsigned char key, int x, int y);
void specialKeyListener(int key, int x, int y);
void drawAxes();
void loadData();
void capture();
Color calculateLighting(Ray* ray, const Vector3D& point, const Vector3D& normal, 
                       const Color& objectColor, const Object* obj, int level);
#endif // RAY_TRACING_MAIN_H