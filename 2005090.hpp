#ifndef RAYTRACER_CLASSES_H
#define RAYTRACER_CLASSES_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <bits/stdc++.h>
#include <GL/glew.h>
using namespace std;

// OpenGL / GLUT Headers
#ifdef __APPLE__
#include <GLUT/glut.h> // Use GLUT framework on macOS
#else
#include <GL/glut.h> // Use standard GLUT location on Linux/Windows
#endif

// Vector3D class for 3D points and directions
class Vector3D
{
public:
    double x, y, z;

    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

    Vector3D operator+(const Vector3D &v) const
    {
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }

    Vector3D operator-(const Vector3D &v) const
    {
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }

    Vector3D operator*(double scalar) const
    {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }

    Vector3D operator/(double scalar) const
    {
        return Vector3D(x / scalar, y / scalar, z / scalar);
    }

    double dot(const Vector3D &v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }

    Vector3D cross(const Vector3D &v) const
    {
        return Vector3D(y * v.z - z * v.y,
                        z * v.x - x * v.z,
                        x * v.y - y * v.x);
    }

    double length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    Vector3D normalize() const
    {
        double len = length();
        if (len == 0)
            return *this;
        return *this / len;
    }
};

// Color class for RGB colors
class Color
{
public:
    double r, g, b;

    Color() : r(0), g(0), b(0) {}
    Color(double r, double g, double b) : r(r), g(g), b(b) {}

    Color operator+(const Color &c) const
    {
        return Color(r + c.r, g + c.g, b + c.b);
    }

    Color operator*(double scalar) const
    {
        return Color(r * scalar, g * scalar, b * scalar);
    }

    Color operator*(const Color &c) const
    {
        return Color(r * c.r, g * c.g, b * c.b);
    }

    void clamp()
    {
        r = std::max(0.0, std::min(1.0, r));
        g = std::max(0.0, std::min(1.0, g));
        b = std::max(0.0, std::min(1.0, b));
    }
};

// Ray class for ray tracing
class Ray
{
public:
    Vector3D start;
    Vector3D dir;

    Ray() {}
    Ray(const Vector3D &start, const Vector3D &dir) : start(start), dir(dir.normalize()) {}
};

// Base Object class
class Object
{
public:
    Vector3D reference_point;
    double height, width, length;
    Color color;
    double coEfficients[4]; // ambient, diffuse, specular, reflection
    int shine;

    Object() : shine(0)
    {
        coEfficients[0] = coEfficients[1] = coEfficients[2] = coEfficients[3] = 0;
    }

    virtual void draw() = 0;
    virtual double intersect(Ray *r, Color *color, int level) = 0;

    void setColor(double r, double g, double b)
    {
        color = Color(r, g, b);
    }

    void setShine(int shine)
    {
        this->shine = shine;
    }

    void setCoEfficients(double ambient, double diffuse, double specular, double reflection)
    {
        coEfficients[0] = ambient;
        coEfficients[1] = diffuse;
        coEfficients[2] = specular;
        coEfficients[3] = reflection;
    }

    virtual Vector3D getNormal(const Vector3D &point) const = 0;
    virtual Color getColorAt(const Vector3D &point) const = 0;
};

// Sphere class
class Sphere : public Object
{
public:
    Sphere(const Vector3D &center, double radius)
    {
        reference_point = center;
        length = radius;
    }

    void draw() override
    {
        // OpenGL code to draw sphere
        glPushMatrix();
        glTranslatef(reference_point.x, reference_point.y, reference_point.z);
        glutSolidSphere(length, 100, 100);
        glPopMatrix();
    }

    double intersect(Ray *r, Color *color, int level) override
    {
        // Ray-sphere intersection implementation
        Vector3D oc = r->start - reference_point;
        double a = r->dir.dot(r->dir);
        double b = 2.0 * oc.dot(r->dir);
        double c = oc.dot(oc) - length * length;
        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0)
            return -1;

        double t = (-b - sqrt(discriminant)) / (2.0 * a);
        if (t < 0)
        {
            t = (-b + sqrt(discriminant)) / (2.0 * a);
            if (t < 0)
                return -1;
        }

        if (level == 0)
            return t;

        // Calculate lighting and color
        Vector3D intersectionPoint = r->start + r->dir * t;
        Vector3D normal = getNormal(intersectionPoint);
        *color = calculateLighting(r, intersectionPoint, normal, level);
        return t;
    }

    Vector3D getNormal(const Vector3D &point) const override
    {
        return (point - reference_point).normalize();
    }

    Color getColorAt(const Vector3D &point) const override
    {
        return color;
    }

private:
    Color calculateLighting(Ray *ray, const Vector3D &point, const Vector3D &normal, int level);
};

// Triangle class
class Triangle : public Object
{
public:
    Vector3D a, b, c;

    Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c) : a(a), b(b), c(c) {}

    void draw() override
    {
        // OpenGL code to draw triangle
        glBegin(GL_TRIANGLES);
        glVertex3f(a.x, a.y, a.z);
        glVertex3f(b.x, b.y, b.z);
        glVertex3f(c.x, c.y, c.z);
        glEnd();
    }

    double intersect(Ray *r, Color *color, int level) override
    {
        // Ray-triangle intersection implementation
        Vector3D edge1 = b - a;
        Vector3D edge2 = c - a;
        Vector3D h = r->dir.cross(edge2);
        double det = edge1.dot(h);

        if (fabs(det) < 1e-6)
            return -1; // Ray parallel to triangle

        double invDet = 1.0 / det;
        Vector3D s = r->start - a;
        double u = invDet * s.dot(h);

        if (u < 0.0 || u > 1.0)
            return -1;

        Vector3D q = s.cross(edge1);
        double v = invDet * r->dir.dot(q);

        if (v < 0.0 || u + v > 1.0)
            return -1;

        double t = invDet * edge2.dot(q);
        if (t < 0)
            return -1;

        if (level == 0)
            return t;

        // Calculate lighting and color
        Vector3D intersectionPoint = r->start + r->dir * t;
        Vector3D normal = getNormal(intersectionPoint);
        *color = calculateLighting(r, intersectionPoint, normal, level);
        return t;
    }

    Vector3D getNormal(const Vector3D &point) const override
    {
        return (b - a).cross(c - a).normalize();
    }

    Color getColorAt(const Vector3D &point) const override
    {
        return color;
    }

private:
    Color calculateLighting(Ray *ray, const Vector3D &point, const Vector3D &normal, int level);
};

// Floor class
class Floor : public Object
{
public:
    double floorWidth, tileWidth;

    Floor(double floorWidth, double tileWidth) : floorWidth(floorWidth), tileWidth(tileWidth)
    {
        reference_point = Vector3D(-floorWidth / 2, -floorWidth / 2, 0);
        length = tileWidth;
    }

    void draw() override
    {
        // OpenGL code to draw checkerboard floor
        bool colorToggle = false;
        for (double x = -floorWidth / 2; x < floorWidth / 2; x += tileWidth)
        {
            for (double y = -floorWidth / 2; y < floorWidth / 2; y += tileWidth)
            {
                if (colorToggle)
                    glColor3f(1, 1, 1); // White
                else
                    glColor3f(0, 0, 0); // Black

                glBegin(GL_QUADS);
                glVertex3f(x, y, 0);
                glVertex3f(x + tileWidth, y, 0);
                glVertex3f(x + tileWidth, y + tileWidth, 0);
                glVertex3f(x, y + tileWidth, 0);
                glEnd();

                colorToggle = !colorToggle;
            }
            if ((int)(floorWidth / tileWidth) % 2 == 0)
                colorToggle = !colorToggle;
        }
    }

    double intersect(Ray *r, Color *color, int level) override
    {
        // Ray-floor intersection implementation
        if (fabs(r->dir.z) < 1e-6)
            return -1;

        double t = -(r->start.z) / r->dir.z;
        if (t < 0)
            return -1;

        Vector3D intersectionPoint = r->start + r->dir * t;
        if (fabs(intersectionPoint.x) > floorWidth / 2 || fabs(intersectionPoint.y) > floorWidth / 2)
            return -1;

        if (level == 0)
            return t;

        // Calculate lighting and color
        Vector3D normal = getNormal(intersectionPoint);
        *color = calculateLighting(r, intersectionPoint, normal, level);
        return t;
    }

    Vector3D getNormal(const Vector3D &point) const override
    {
        return Vector3D(0, 0, 1);
    }

    Color getColorAt(const Vector3D &point) const override
    {
        // Checkerboard pattern
        int x = (int)((point.x - reference_point.x) / tileWidth);
        int y = (int)((point.y - reference_point.y) / tileWidth);

        if ((x + y) % 2 == 0)
            return Color(1, 1, 1); // White
        else
            return Color(0, 0, 0); // Black
    }

private:
    Color calculateLighting(Ray *ray, const Vector3D &point, const Vector3D &normal, int level);
};

// General Quadric Surface class
class GeneralQuadric : public Object
{
public:
    double coefficients[10]; // A-J coefficients

    GeneralQuadric(double A, double B, double C, double D, double E,
                   double F, double G, double H, double I, double J)
    {
        coefficients[0] = A;
        coefficients[1] = B;
        coefficients[2] = C;
        coefficients[3] = D;
        coefficients[4] = E;
        coefficients[5] = F;
        coefficients[6] = G;
        coefficients[7] = H;
        coefficients[8] = I;
        coefficients[9] = J;
    }

    void draw() override
    {
        // Not required to draw for general quadric surfaces
    }

    double intersect(Ray *r, Color *color, int level) override
    {
        // Ray-quadric intersection implementation
        double Aq = coefficients[0] * r->dir.x * r->dir.x +
                    coefficients[1] * r->dir.y * r->dir.y +
                    coefficients[2] * r->dir.z * r->dir.z +
                    coefficients[3] * r->dir.x * r->dir.y +
                    coefficients[4] * r->dir.x * r->dir.z +
                    coefficients[5] * r->dir.y * r->dir.z;

        double Bq = 2 * coefficients[0] * r->start.x * r->dir.x +
                    2 * coefficients[1] * r->start.y * r->dir.y +
                    2 * coefficients[2] * r->start.z * r->dir.z +
                    coefficients[3] * (r->start.x * r->dir.y + r->start.y * r->dir.x) +
                    coefficients[4] * (r->start.x * r->dir.z + r->start.z * r->dir.x) +
                    coefficients[5] * (r->start.y * r->dir.z + r->start.z * r->dir.y) +
                    coefficients[6] * r->dir.x +
                    coefficients[7] * r->dir.y +
                    coefficients[8] * r->dir.z;

        double Cq = coefficients[0] * r->start.x * r->start.x +
                    coefficients[1] * r->start.y * r->start.y +
                    coefficients[2] * r->start.z * r->start.z +
                    coefficients[3] * r->start.x * r->start.y +
                    coefficients[4] * r->start.x * r->start.z +
                    coefficients[5] * r->start.y * r->start.z +
                    coefficients[6] * r->start.x +
                    coefficients[7] * r->start.y +
                    coefficients[8] * r->start.z +
                    coefficients[9];

        double discriminant = Bq * Bq - 4 * Aq * Cq;
        if (discriminant < 0)
            return -1;

        double t1 = (-Bq - sqrt(discriminant)) / (2 * Aq);
        double t2 = (-Bq + sqrt(discriminant)) / (2 * Aq);

        double t = -1;
        if (t1 > 0)
            t = t1;
        if (t2 > 0 && (t2 < t1 || t < 0))
            t = t2;
        if (t < 0)
            return -1;

        // Check if intersection is within reference cube
        Vector3D intersectionPoint = r->start + r->dir * t;
        if (!isWithinReferenceCube(intersectionPoint))
            return -1;

        if (level == 0)
            return t;

        // Calculate lighting and color
        Vector3D normal = getNormal(intersectionPoint);
        *color = calculateLighting(r, intersectionPoint, normal, level);
        return t;
    }

    Vector3D getNormal(const Vector3D &point) const override
    {
        // Normal is gradient of the quadric equation
        double nx = 2 * coefficients[0] * point.x +
                    coefficients[3] * point.y +
                    coefficients[4] * point.z +
                    coefficients[6];

        double ny = 2 * coefficients[1] * point.y +
                    coefficients[3] * point.x +
                    coefficients[5] * point.z +
                    coefficients[7];

        double nz = 2 * coefficients[2] * point.z +
                    coefficients[4] * point.x +
                    coefficients[5] * point.y +
                    coefficients[8];

        return Vector3D(nx, ny, nz).normalize();
    }

    Color getColorAt(const Vector3D &point) const override
    {
        return color;
    }

private:
    bool isWithinReferenceCube(const Vector3D &point) const
    {
        if (width > 0 && (point.x < reference_point.x || point.x > reference_point.x + width))
            return false;
        if (height > 0 && (point.y < reference_point.y || point.y > reference_point.y + height))
            return false;
        if (length > 0 && (point.z < reference_point.z || point.z > reference_point.z + length))
            return false;
        return true;
    }

    Color calculateLighting(Ray *ray, const Vector3D &point, const Vector3D &normal, int level);
};

// Light classes
class PointLight
{
public:
    Vector3D light_pos;
    Color color;

    PointLight(const Vector3D &pos, const Color &color) : light_pos(pos), color(color) {}

    void draw()
    {
        // Draw point light as a small sphere
        glPushMatrix();
        glTranslatef(light_pos.x, light_pos.y, light_pos.z);
        glColor3f(color.r, color.g, color.b);
        glutSolidSphere(0.2, 10, 10);
        glPopMatrix();
    }
};

class SpotLight
{
public:
    PointLight pointLight;
    Vector3D light_dir;
    double cutoff_angle;

    SpotLight(const Vector3D &pos, const Color &color, const Vector3D &dir, double angle)
        : pointLight(pos, color), light_dir(dir.normalize()), cutoff_angle(angle) {}

    void draw()
    {
        pointLight.draw();
        // Additional visualization for spotlight direction
    }
};

// Texture class for floor
class Texture
{
private:
    unsigned char *textureData;
    int textureWidth, textureHeight, textureChannels;

public:
    Texture() : textureData(nullptr), textureWidth(0), textureHeight(0), textureChannels(0) {}

    bool loadTexture(const char *filename)
    {
        // Implementation to load texture from file
        // You can use libraries like stb_image.h for this
        return true;
    }

    Color sampleTexture(double u, double v) const
    {
        if (!textureData || textureWidth <= 0 || textureHeight <= 0)
        {
            return Color(0.5, 0.5, 0.5); // Gray fallback
        }

        // Clamp u and v to [0,1]
        u = std::max(0.0, std::min(1.0, u));
        v = std::max(0.0, std::min(1.0, v));

        // Normalized -> pixel coords
        int pixel_x = (int)(u * (textureWidth - 1));
        int pixel_y = (int)((1.0 - v) * (textureHeight - 1)); // Flip Y

        // Safety clamp
        pixel_x = std::max(0, std::min(textureWidth - 1, pixel_x));
        pixel_y = std::max(0, std::min(textureHeight - 1, pixel_y));

        // Compute array index
        int index = (pixel_y * textureWidth + pixel_x) * textureChannels;
        int max_index = textureWidth * textureHeight * textureChannels;
        if (index < 0 || index + 2 >= max_index)
        {
            return Color(1.0, 0.0, 1.0); // Magenta = error
        }

        Color color;
        color.r = textureData[index] / 255.0;

        if (textureChannels >= 2)
        {
            color.g = textureData[index + 1] / 255.0;
        }
        else
        {
            color.g = color.r; // Grayscale
        }

        if (textureChannels >= 3)
        {
            color.b = textureData[index + 2] / 255.0;
        }
        else
        {
            color.b = color.r; // Grayscale
        }

        return color;
    }
};

// Global variables (extern declarations)
extern std::vector<Object *> objects;
extern std::vector<PointLight> pointLights;
extern std::vector<SpotLight> spotLights;
extern Texture floorTexture;
extern bool useTexture;

#endif // RAYTRACER_CLASSES_H