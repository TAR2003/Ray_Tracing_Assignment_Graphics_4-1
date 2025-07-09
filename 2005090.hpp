#ifndef SCENE_OBJECTS_HPP
#define SCENE_OBJECTS_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>

// Include stb_image for texture loading
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// OpenGL / GLUT Headers
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

// Forward declarations
class Object;
class PointLight;
class SpotLight;

class Vector3D
{
public:
    double x, y, z;

    Vector3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

    Vector3D operator+(const Vector3D &v) const
    {
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }

    Vector3D operator-(const Vector3D &v) const
    {
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }

    Vector3D operator*(double s) const
    {
        return Vector3D(x * s, y * s, z * s);
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

    void normalize()
    {
        double len = length();
        if (len > 0)
        {
            x /= len;
            y /= len;
            z /= len;
        }
    }
};

class Color
{
public:
    double r, g, b;
    
    Color(double r = 0.0, double g = 0.0, double b = 0.0) : r(r), g(g), b(b) {}
    
    void toArray(double *array) const
    {
        array[0] = r;
        array[1] = g;
        array[2] = b;
    }
};

// Texture loading and sampling utilities
class TextureManager
{
public:
    static unsigned char* textureData;
    static int textureWidth;
    static int textureHeight;
    static int textureChannels;
    static bool textureLoaded;
    
    static bool loadTexture(const char* filename)
    {
        if (textureData) {
            stbi_image_free(textureData);
        }
        
        textureData = stbi_load(filename, &textureWidth, &textureHeight, &textureChannels, 0);
        textureLoaded = (textureData != nullptr);
        
        if (textureLoaded) {
            std::cout << "Texture loaded: " << filename << " (" << textureWidth << "x" << textureHeight << ", " << textureChannels << " channels)" << std::endl;
        } else {
            std::cout << "Failed to load texture: " << filename << std::endl;
        }
        
        return textureLoaded;
    }
    
    static Color sampleTexture(double u, double v)
    {
        if (!textureData || textureWidth <= 0 || textureHeight <= 0) {
            return Color(0.5, 0.5, 0.5); // Gray fallback
        }
        
        // Clamp u and v to [0,1]
        u = std::max(0.0, std::min(1.0, u));
        v = std::max(0.0, std::min(1.0, v));
        
        // Normalized-> pixel coords
        int pixel_x = (int)(u * (textureWidth - 1));
        int pixel_y = (int)((1.0 - v) * (textureHeight - 1)); // Flip Y
        
        // Safety clamp
        pixel_x = std::max(0, std::min(textureWidth - 1, pixel_x));
        pixel_y = std::max(0, std::min(textureHeight - 1, pixel_y));
        
        // Compute array index
        int index = (pixel_y * textureWidth + pixel_x) * textureChannels;
        int max_index = textureWidth * textureHeight * textureChannels;
        if (index < 0 || index + 2 >= max_index) {
            return Color(1.0, 0.0, 1.0); // Magenta = error
        }
        
        Color color;
        color.r = textureData[index] / 255.0;
        
        if (textureChannels >= 2) {
            color.g = textureData[index + 1] / 255.0;
        } else {
            color.g = color.r; // Grayscale
        }
        
        if (textureChannels >= 3) {
            color.b = textureData[index + 2] / 255.0;
        } else {
            color.b = color.r; // Grayscale
        }
        
        return color;
    }
    
    static void cleanup()
    {
        if (textureData) {
            stbi_image_free(textureData);
            textureData = nullptr;
            textureLoaded = false;
        }
    }
};

class Ray
{
public:
    Vector3D start;
    Vector3D dir;

    Ray(const Vector3D &start, const Vector3D &dir) : start(start), dir(dir)
    {
        this->dir.normalize();
    }
};

class PointLight
{
public:
    Vector3D position;
    double color[3];

    PointLight(const Vector3D &pos, double r, double g, double b) : position(pos)
    {
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void draw()
    {
        glPushMatrix();
        glTranslatef(position.x, position.y, position.z);
        glColor3f(color[0], color[1], color[2]);
        glutSolidSphere(0.2, 10, 10);
        glPopMatrix();
    }
};

class SpotLight : public PointLight
{
public:
    Vector3D direction;
    double cutoffAngle;

    SpotLight(const Vector3D &pos, const Vector3D &dir, double angle,
              double r, double g, double b)
        : PointLight(pos, r, g, b), direction(dir), cutoffAngle(angle)
    {
        direction.normalize();
    }

    void draw()
    {
        PointLight::draw();
        // Additional visualization for spotlight direction
    }
};

// External references to light arrays (defined in main cpp file)
extern std::vector<PointLight *> pointLights;
extern std::vector<SpotLight *> spotLights;
extern std::vector<Object *> objects;
extern int maxRecursionLevel;

class Object
{
public:
    Vector3D reference_point;
    double height, width, length;
    double color[3];
    double coefficients[4]; // ambient, diffuse, specular, reflection
    int shine;

    Object() : height(0), width(0), length(0), shine(0)
    {
        color[0] = color[1] = color[2] = 0.0;
        coefficients[0] = coefficients[1] = coefficients[2] = coefficients[3] = 0.0;
    }

    virtual ~Object() {}
    virtual void draw() = 0;
    virtual double intersect(const Ray &r, double *color, int level) const = 0;
    virtual void getNormal(const Vector3D &point, Vector3D &normal) const = 0;

    void setColor(double r, double g, double b)
    {
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }

    void setShine(int sh)
    {
        shine = sh;
    }

    void setCoefficients(double ambient, double diffuse, double specular, double reflection)
    {
        coefficients[0] = ambient;
        coefficients[1] = diffuse;
        coefficients[2] = specular;
        coefficients[3] = reflection;
    }

protected:
    // Helper method to get color at intersection point (can be overridden for textured objects)
    virtual void getColorAt(const Vector3D &point, double *colorAt) const
    {
        colorAt[0] = color[0];
        colorAt[1] = color[1];
        colorAt[2] = color[2];
    }
};

class Sphere : public Object
{
public:
    double radius;

    Sphere(const Vector3D &center, double r) : radius(r)
    {
        reference_point = center;
        height = width = length = r * 2;
    }

    void draw() override
    {
        glPushMatrix();
        glTranslatef(reference_point.x, reference_point.y, reference_point.z);
        glColor3f(color[0], color[1], color[2]);
        glutSolidSphere(radius, 50, 50);
        glPopMatrix();
    }

    double intersect(const Ray &r, double *color, int level) const override
    {
        Vector3D oc = r.start - reference_point;
        double a = r.dir.dot(r.dir);
        double b = 2.0 * oc.dot(r.dir);
        double c = oc.dot(oc) - radius * radius;
        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0)
            return -1.0;

        double sqrt_discriminant = sqrt(discriminant);
        double t1 = (-b - sqrt_discriminant) / (2.0 * a);
        double t2 = (-b + sqrt_discriminant) / (2.0 * a);

        double t = (t1 > 0) ? t1 : t2;
        if (t <= 0)
            return -1.0;

        // If level is 0, only return intersection distance
        if (level == 0)
            return t;

        // Calculate intersection point and lighting
        Vector3D intersectionPoint = r.start + r.dir * t;
        Vector3D normal;
        getNormal(intersectionPoint, normal);

        // Get color at intersection point
        double intersectionPointColor[3];
        getColorAt(intersectionPoint, intersectionPointColor);

        // Initialize with ambient component
        color[0] = intersectionPointColor[0] * coefficients[0]; // ambient
        color[1] = intersectionPointColor[1] * coefficients[0];
        color[2] = intersectionPointColor[2] * coefficients[0];

        // Process each point light
        for (PointLight *pl : pointLights)
        {
            // Cast ray from light to intersection point
            Vector3D lightDir = intersectionPoint - pl->position;
            double lightDistance = lightDir.length();
            lightDir.normalize();
            Ray lightRay(pl->position, lightDir);

            // Check if intersection point is in shadow
            bool inShadow = false;
            for (Object *obj : objects)
            {
                if (obj != this) // Don't check self-intersection
                {
                    double shadowT = obj->intersect(lightRay, nullptr, 0);
                    if (shadowT > 0.001 && shadowT < lightDistance - 0.001) // Small epsilon for numerical stability
                    {
                        inShadow = true;
                        break;
                    }
                }
            }

            if (!inShadow)
            {
                // Calculate diffuse component (Lambert)
                Vector3D toLight = pl->position - intersectionPoint;
                toLight.normalize();
                double lambertValue = std::max(0.0, normal.dot(toLight));

                // Calculate specular component (Phong)
                Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                reflected.normalize();
                Vector3D toCamera = r.start - intersectionPoint;
                toCamera.normalize();
                double phongValue = std::max(0.0, reflected.dot(toCamera));
                phongValue = pow(phongValue, shine);

                // Add diffuse and specular contributions
                for (int i = 0; i < 3; i++)
                {
                    color[i] += pl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i]; // diffuse
                    color[i] += pl->color[i] * coefficients[2] * phongValue; // specular
                }
            }
        }

        // Process each spot light
        for (SpotLight *sl : spotLights)
        {
            Vector3D lightToPoint = intersectionPoint - sl->position;
            double lightDistance = lightToPoint.length();
            lightToPoint.normalize();

            // Check if point is within spotlight cone
            double angle = acos(lightToPoint.dot(sl->direction));
            if (angle <= sl->cutoffAngle)
            {
                Ray lightRay(sl->position, lightToPoint);

                // Check shadow
                bool inShadow = false;
                for (Object *obj : objects)
                {
                    if (obj != this) // Don't check self-intersection
                    {
                        double shadowT = obj->intersect(lightRay, nullptr, 0);
                        if (shadowT > 0.001 && shadowT < lightDistance - 0.001)
                        {
                            inShadow = true;
                            break;
                        }
                    }
                }

                if (!inShadow)
                {
                    Vector3D toLight = sl->position - intersectionPoint;
                    toLight.normalize();
                    double lambertValue = std::max(0.0, normal.dot(toLight));

                    Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                    reflected.normalize();
                    Vector3D toCamera = r.start - intersectionPoint;
                    toCamera.normalize();
                    double phongValue = std::max(0.0, reflected.dot(toCamera));
                    phongValue = pow(phongValue, shine);

                    for (int i = 0; i < 3; i++)
                    {
                        color[i] += sl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i];
                        color[i] += sl->color[i] * coefficients[2] * phongValue;
                    }
                }
            }
        }

        // Handle reflection
        if (level < maxRecursionLevel && coefficients[3] > 0) // coefficients[3] is reflection coefficient
        {
            // Calculate reflected ray direction
            Vector3D incident = r.dir;
            Vector3D reflected = incident - normal * (2.0 * incident.dot(normal));
            reflected.normalize();
            
            // Create reflected ray starting slightly forward from intersection point to avoid self-intersection
            Vector3D reflectedStart = intersectionPoint + reflected * 0.001;
            Ray reflectedRay(reflectedStart, reflected);
            
            // Find nearest intersecting object for the reflected ray
            double minT = std::numeric_limits<double>::max();
            Object *nearestObject = nullptr;
            
            for (Object *obj : objects)
            {
                double t = obj->intersect(reflectedRay, nullptr, 0);
                if (t > 0 && t < minT)
                {
                    minT = t;
                    nearestObject = obj;
                }
            }
            
            // If reflection ray hits an object, calculate reflected color
            if (nearestObject != nullptr)
            {
                double reflectedColor[3] = {0, 0, 0};
                nearestObject->intersect(reflectedRay, reflectedColor, level + 1);
                
                // Add reflected color contribution
                for (int i = 0; i < 3; i++)
                {
                    color[i] += reflectedColor[i] * coefficients[3]; // coefficients[3] is reflection coefficient
                }
            }
        }

        return t;
    }

    void getNormal(const Vector3D &point, Vector3D &normal) const override
    {
        normal = point - reference_point;
        normal.normalize();
    }
};

class Floor : public Object
{
public:
    double floorWidth;
    double tileWidth;
    bool useTexture; // Flag to switch between checkerboard and texture

    Floor(double width, double tileSize, bool texture = false) : floorWidth(width), tileWidth(tileSize), useTexture(texture)
    {
        reference_point = Vector3D(-width / 2, -width / 2, 0);
        length = tileSize;
        setColor(0.7, 0.7, 0.7);
        setCoefficients(0.4, 0.3, 0.2, 0.1);
        setShine(5);
    }
    
    void setTextureMode(bool texture)
    {
        useTexture = texture;
    }

protected:
    // Override getColorAt to handle both checkerboard and texture
    void getColorAt(const Vector3D &point, double *colorAt) const override
    {
        if (useTexture && TextureManager::textureLoaded)
        {
            // Calculate texture coordinates (u, v) from world coordinates
            // Map floor coordinates to [0,1] texture space
            double u = (point.x + floorWidth / 2) / floorWidth;
            double v = (point.y + floorWidth / 2) / floorWidth;
            
            // Sample texture
            Color texColor = TextureManager::sampleTexture(u, v);
            colorAt[0] = texColor.r;
            colorAt[1] = texColor.g;
            colorAt[2] = texColor.b;
        }
        else
        {
            // Use checkerboard pattern
            int tx = static_cast<int>((point.x + floorWidth / 2) / tileWidth);
            int ty = static_cast<int>((point.y + floorWidth / 2) / tileWidth);
            
            if ((tx + ty) % 2 == 0)
            {
                colorAt[0] = colorAt[1] = colorAt[2] = 0.9; // White tile
            }
            else
            {
                colorAt[0] = colorAt[1] = colorAt[2] = 0.1; // Black tile
            }
        }
    }

    void draw() override
    {
        int tileCount = floorWidth / tileWidth;
        glPushMatrix();
        glBegin(GL_QUADS);

        for (int i = 0; i < tileCount; i++)
        {
            for (int j = 0; j < tileCount; j++)
            {
                if ((i + j) % 2 == 0)
                {
                    glColor3f(0.9, 0.9, 0.9);
                }
                else
                {
                    glColor3f(0.1, 0.1, 0.1);
                }

                double x1 = reference_point.x + i * tileWidth;
                double y1 = reference_point.y + j * tileWidth;
                double x2 = x1 + tileWidth;
                double y2 = y1 + tileWidth;

                glVertex3f(x1, y1, 0);
                glVertex3f(x2, y1, 0);
                glVertex3f(x2, y2, 0);
                glVertex3f(x1, y2, 0);
            }
        }

        glEnd();
        glPopMatrix();
    }

    double intersect(const Ray &r, double *color, int level) const override
    {
        if (fabs(r.dir.z) < 1e-6)
            return -1.0;

        double t = -r.start.z / r.dir.z;
        if (t <= 0)
            return -1.0;

        Vector3D intersection = r.start + r.dir * t;
        if (fabs(intersection.x) > floorWidth / 2 || fabs(intersection.y) > floorWidth / 2)
        {
            return -1.0;
        }

        if (level == 0)
            return t;

        // Get color at intersection point (handles both texture and checkerboard)
        double intersectionPointColor[3];
        getColorAt(intersection, intersectionPointColor);

        // Initialize with ambient component
        color[0] = intersectionPointColor[0] * coefficients[0];
        color[1] = intersectionPointColor[1] * coefficients[0];
        color[2] = intersectionPointColor[2] * coefficients[0];

        Vector3D normal = Vector3D(0, 0, 1); // Floor normal always points up

        // Process each point light
        for (PointLight *pl : pointLights)
        {
            Vector3D lightDir = intersection - pl->position;
            double lightDistance = lightDir.length();
            lightDir.normalize();
            Ray lightRay(pl->position, lightDir);

            // Check if intersection point is in shadow
            bool inShadow = false;
            for (Object *obj : objects)
            {
                if (obj != this)
                {
                    double shadowT = obj->intersect(lightRay, nullptr, 0);
                    if (shadowT > 0.001 && shadowT < lightDistance - 0.001)
                    {
                        inShadow = true;
                        break;
                    }
                }
            }

            if (!inShadow)
            {
                Vector3D toLight = pl->position - intersection;
                toLight.normalize();
                double lambertValue = std::max(0.0, normal.dot(toLight));

                Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                reflected.normalize();
                Vector3D toCamera = r.start - intersection;
                toCamera.normalize();
                double phongValue = std::max(0.0, reflected.dot(toCamera));
                phongValue = pow(phongValue, shine);

                for (int i = 0; i < 3; i++)
                {
                    color[i] += pl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i];
                    color[i] += pl->color[i] * coefficients[2] * phongValue;
                }
            }
        }

        // Process each spot light
        for (SpotLight *sl : spotLights)
        {
            Vector3D lightToPoint = intersection - sl->position;
            double lightDistance = lightToPoint.length();
            lightToPoint.normalize();

            double angle = acos(lightToPoint.dot(sl->direction));
            if (angle <= sl->cutoffAngle)
            {
                Ray lightRay(sl->position, lightToPoint);

                bool inShadow = false;
                for (Object *obj : objects)
                {
                    if (obj != this)
                    {
                        double shadowT = obj->intersect(lightRay, nullptr, 0);
                        if (shadowT > 0.001 && shadowT < lightDistance - 0.001)
                        {
                            inShadow = true;
                            break;
                        }
                    }
                }

                if (!inShadow)
                {
                    Vector3D toLight = sl->position - intersection;
                    toLight.normalize();
                    double lambertValue = std::max(0.0, normal.dot(toLight));

                    Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                    reflected.normalize();
                    Vector3D toCamera = r.start - intersection;
                    toCamera.normalize();
                    double phongValue = std::max(0.0, reflected.dot(toCamera));
                    phongValue = pow(phongValue, shine);

                    for (int i = 0; i < 3; i++)
                    {
                        color[i] += sl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i];
                        color[i] += sl->color[i] * coefficients[2] * phongValue;
                    }
                }
            }
        }

        // Handle reflection
        if (level < maxRecursionLevel && coefficients[3] > 0) // coefficients[3] is reflection coefficient
        {
            // Calculate reflected ray direction
            Vector3D incident = r.dir;
            Vector3D reflected = incident - normal * (2.0 * incident.dot(normal));
            reflected.normalize();
            
            // Create reflected ray starting slightly forward from intersection point to avoid self-intersection
            Vector3D reflectedStart = intersection + reflected * 0.001;
            Ray reflectedRay(reflectedStart, reflected);
            
            // Find nearest intersecting object for the reflected ray
            double minT = std::numeric_limits<double>::max();
            Object *nearestObject = nullptr;
            
            for (Object *obj : objects)
            {
                double t = obj->intersect(reflectedRay, nullptr, 0);
                if (t > 0 && t < minT)
                {
                    minT = t;
                    nearestObject = obj;
                }
            }
            
            // If reflection ray hits an object, calculate reflected color
            if (nearestObject != nullptr)
            {
                double reflectedColor[3] = {0, 0, 0};
                nearestObject->intersect(reflectedRay, reflectedColor, level + 1);
                
                // Add reflected color contribution
                for (int i = 0; i < 3; i++)
                {
                    color[i] += reflectedColor[i] * coefficients[3]; // coefficients[3] is reflection coefficient
                }
            }
        }

        return t;
    }

    void getNormal(const Vector3D &point, Vector3D &normal) const override
    {
        normal = Vector3D(0, 0, 1);
    }
};

class Triangle : public Object
{
public:
    Vector3D a, b, c; // Three vertices

    Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c) : a(a), b(b), c(c)
    {
        // Calculate reference point as centroid
        reference_point = Vector3D((a.x + b.x + c.x) / 3.0,
                                   (a.y + b.y + c.y) / 3.0,
                                   (a.z + b.z + c.z) / 3.0);
    }

    void draw() override
    {
        glPushMatrix();
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
        glVertex3f(a.x, a.y, a.z);
        glVertex3f(b.x, b.y, b.z);
        glVertex3f(c.x, c.y, c.z);
        glEnd();
        glPopMatrix();
    }

    double intersect(const Ray &r, double *color, int level) const override
    {
        // Ray-triangle intersection using Möller-Trumbore algorithm
        Vector3D edge1 = b - a;
        Vector3D edge2 = c - a;
        Vector3D h = r.dir.cross(edge2);
        double a_det = edge1.dot(h);

        if (a_det > -1e-6 && a_det < 1e-6)
            return -1.0; // Ray is parallel to triangle

        double f = 1.0 / a_det;
        Vector3D s = r.start - a;
        double u = f * s.dot(h);

        if (u < 0.0 || u > 1.0)
            return -1.0;

        Vector3D q = s.cross(edge1);
        double v = f * r.dir.dot(q);

        if (v < 0.0 || u + v > 1.0)
            return -1.0;

        double t = f * edge2.dot(q);

        if (t > 1e-6) // Ray intersection
        {
            if (level == 0)
                return t;

            // Calculate intersection point and lighting
            Vector3D intersectionPoint = r.start + r.dir * t;
            Vector3D normal;
            getNormal(intersectionPoint, normal);

            // Get color at intersection point
            double intersectionPointColor[3];
            getColorAt(intersectionPoint, intersectionPointColor);

            // Initialize with ambient component
            color[0] = intersectionPointColor[0] * coefficients[0]; // ambient
            color[1] = intersectionPointColor[1] * coefficients[0];
            color[2] = intersectionPointColor[2] * coefficients[0];

            // Process each point light
            for (PointLight *pl : pointLights)
            {
                // Cast ray from light to intersection point
                Vector3D lightDir = intersectionPoint - pl->position;
                double lightDistance = lightDir.length();
                lightDir.normalize();
                Ray lightRay(pl->position, lightDir);

                // Check if intersection point is in shadow
                bool inShadow = false;
                for (Object *obj : objects)
                {
                    if (obj != this) // Don't check self-intersection
                    {
                        double shadowT = obj->intersect(lightRay, nullptr, 0);
                        if (shadowT > 0.001 && shadowT < lightDistance - 0.001) // Small epsilon for numerical stability
                        {
                            inShadow = true;
                            break;
                        }
                    }
                }

                if (!inShadow)
                {
                    // Calculate diffuse component (Lambert)
                    Vector3D toLight = pl->position - intersectionPoint;
                    toLight.normalize();
                    double lambertValue = std::max(0.0, normal.dot(toLight));

                    // Calculate specular component (Phong)
                    Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                    reflected.normalize();
                    Vector3D toCamera = r.start - intersectionPoint;
                    toCamera.normalize();
                    double phongValue = std::max(0.0, reflected.dot(toCamera));
                    phongValue = pow(phongValue, shine);

                    // Add diffuse and specular contributions
                    for (int i = 0; i < 3; i++)
                    {
                        color[i] += pl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i]; // diffuse
                        color[i] += pl->color[i] * coefficients[2] * phongValue; // specular
                    }
                }
            }

            // Process each spot light
            for (SpotLight *sl : spotLights)
            {
                Vector3D lightToPoint = intersectionPoint - sl->position;
                double lightDistance = lightToPoint.length();
                lightToPoint.normalize();

                // Check if point is within spotlight cone
                double angle = acos(lightToPoint.dot(sl->direction));
                if (angle <= sl->cutoffAngle)
                {
                    Ray lightRay(sl->position, lightToPoint);

                    // Check shadow
                    bool inShadow = false;
                    for (Object *obj : objects)
                    {
                        if (obj != this) // Don't check self-intersection
                        {
                            double shadowT = obj->intersect(lightRay, nullptr, 0);
                            if (shadowT > 0.001 && shadowT < lightDistance - 0.001)
                            {
                                inShadow = true;
                                break;
                            }
                        }
                    }

                    if (!inShadow)
                    {
                        Vector3D toLight = sl->position - intersectionPoint;
                        toLight.normalize();
                        double lambertValue = std::max(0.0, normal.dot(toLight));

                        Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                        reflected.normalize();
                        Vector3D toCamera = r.start - intersectionPoint;
                        toCamera.normalize();
                        double phongValue = std::max(0.0, reflected.dot(toCamera));
                        phongValue = pow(phongValue, shine);

                        for (int i = 0; i < 3; i++)
                        {
                            color[i] += sl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i];
                            color[i] += sl->color[i] * coefficients[2] * phongValue;
                        }
                    }
                }
            }

            // Handle reflection
            if (level < maxRecursionLevel && coefficients[3] > 0) // coefficients[3] is reflection coefficient
            {
                // Calculate reflected ray direction
                Vector3D incident = r.dir;
                Vector3D reflected = incident - normal * (2.0 * incident.dot(normal));
                reflected.normalize();
                
                // Create reflected ray starting slightly forward from intersection point to avoid self-intersection
                Vector3D reflectedStart = intersectionPoint + reflected * 0.001;
                Ray reflectedRay(reflectedStart, reflected);
                
                // Find nearest intersecting object for the reflected ray
                double minT = std::numeric_limits<double>::max();
                Object *nearestObject = nullptr;
                
                for (Object *obj : objects)
                {
                    double t = obj->intersect(reflectedRay, nullptr, 0);
                    if (t > 0 && t < minT)
                    {
                        minT = t;
                        nearestObject = obj;
                    }
                }
                
                // If reflection ray hits an object, calculate reflected color
                if (nearestObject != nullptr)
                {
                    double reflectedColor[3] = {0, 0, 0};
                    nearestObject->intersect(reflectedRay, reflectedColor, level + 1);
                    
                    // Add reflected color contribution
                    for (int i = 0; i < 3; i++)
                    {
                        color[i] += reflectedColor[i] * coefficients[3]; // coefficients[3] is reflection coefficient
                    }
                }
            }

            return t;
        }

        return -1.0; // Line intersection but not ray intersection
    }

    void getNormal(const Vector3D &point, Vector3D &normal) const override
    {
        Vector3D edge1 = b - a;
        Vector3D edge2 = c - a;
        normal = edge1.cross(edge2);
        normal.normalize();
    }
};

class GeneralQuadric : public Object
{
public:
    double coeffs[10]; // A B C D E F G H I J coefficients

    GeneralQuadric(double A, double B, double C, double D, double E,
                   double F, double G, double H, double I, double J)
    {
        coeffs[0] = A; coeffs[1] = B; coeffs[2] = C; coeffs[3] = D; coeffs[4] = E;
        coeffs[5] = F; coeffs[6] = G; coeffs[7] = H; coeffs[8] = I; coeffs[9] = J;
    }

    void draw() override
    {
        // Simple visualization - draw a wireframe representation
        glPushMatrix();
        glColor3f(color[0], color[1], color[2]);
        glTranslatef(reference_point.x, reference_point.y, reference_point.z);
        glutWireCube(2.0);
        glPopMatrix();
    }

    double intersect(const Ray &r, double *color, int level) const override
    {
        // Quadric surface equation: Ax² + By² + Cz² + Dxy + Exz + Fyz + Gx + Hy + Iz + J = 0
        // For ray: P = P0 + t*d, substitute and solve quadratic equation

        Vector3D o = r.start;
        Vector3D d = r.dir;

        double a = coeffs[0] * d.x * d.x + coeffs[1] * d.y * d.y + coeffs[2] * d.z * d.z +
                   coeffs[3] * d.x * d.y + coeffs[4] * d.x * d.z + coeffs[5] * d.y * d.z;

        double b = 2 * coeffs[0] * o.x * d.x + 2 * coeffs[1] * o.y * d.y + 2 * coeffs[2] * o.z * d.z +
                   coeffs[3] * (o.x * d.y + o.y * d.x) + coeffs[4] * (o.x * d.z + o.z * d.x) +
                   coeffs[5] * (o.y * d.z + o.z * d.y) + coeffs[6] * d.x + coeffs[7] * d.y + coeffs[8] * d.z;

        double c = coeffs[0] * o.x * o.x + coeffs[1] * o.y * o.y + coeffs[2] * o.z * o.z +
                   coeffs[3] * o.x * o.y + coeffs[4] * o.x * o.z + coeffs[5] * o.y * o.z +
                   coeffs[6] * o.x + coeffs[7] * o.y + coeffs[8] * o.z + coeffs[9];

        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0)
            return -1.0;

        double sqrt_discriminant = sqrt(discriminant);
        double t1 = (-b - sqrt_discriminant) / (2.0 * a);
        double t2 = (-b + sqrt_discriminant) / (2.0 * a);

        double t = (t1 > 1e-6) ? t1 : t2;
        if (t <= 1e-6)
            return -1.0;

        // Check if intersection point is within reference cube bounds
        Vector3D intersection = r.start + r.dir * t;
        if (!isWithinReferenceCube(intersection))
            return -1.0;

        if (level == 0)
            return t;

        // Calculate intersection point and lighting
        Vector3D intersectionPoint = r.start + r.dir * t;
        Vector3D normal;
        getNormal(intersectionPoint, normal);

        // Get color at intersection point
        double intersectionPointColor[3];
        getColorAt(intersectionPoint, intersectionPointColor);

        // Initialize with ambient component
        color[0] = intersectionPointColor[0] * coefficients[0]; // ambient
        color[1] = intersectionPointColor[1] * coefficients[0];
        color[2] = intersectionPointColor[2] * coefficients[0];

        // Process each point light
        for (PointLight *pl : pointLights)
        {
            // Cast ray from light to intersection point
            Vector3D lightDir = intersectionPoint - pl->position;
            double lightDistance = lightDir.length();
            lightDir.normalize();
            Ray lightRay(pl->position, lightDir);

            // Check if intersection point is in shadow
            bool inShadow = false;
            for (Object *obj : objects)
            {
                if (obj != this) // Don't check self-intersection
                {
                    double shadowT = obj->intersect(lightRay, nullptr, 0);
                    if (shadowT > 0.001 && shadowT < lightDistance - 0.001) // Small epsilon for numerical stability
                    {
                        inShadow = true;
                        break;
                    }
                }
            }

            if (!inShadow)
            {
                // Calculate diffuse component (Lambert)
                Vector3D toLight = pl->position - intersectionPoint;
                toLight.normalize();
                double lambertValue = std::max(0.0, normal.dot(toLight));

                // Calculate specular component (Phong)
                Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                reflected.normalize();
                Vector3D toCamera = r.start - intersectionPoint;
                toCamera.normalize();
                double phongValue = std::max(0.0, reflected.dot(toCamera));
                phongValue = pow(phongValue, shine);

                // Add diffuse and specular contributions
                for (int i = 0; i < 3; i++)
                {
                    color[i] += pl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i]; // diffuse
                    color[i] += pl->color[i] * coefficients[2] * phongValue; // specular
                }
            }
        }

        // Process each spot light
        for (SpotLight *sl : spotLights)
        {
            Vector3D lightToPoint = intersectionPoint - sl->position;
            double lightDistance = lightToPoint.length();
            lightToPoint.normalize();

            // Check if point is within spotlight cone
            double angle = acos(lightToPoint.dot(sl->direction));
            if (angle <= sl->cutoffAngle)
            {
                Ray lightRay(sl->position, lightToPoint);

                // Check shadow
                bool inShadow = false;
                for (Object *obj : objects)
                {
                    if (obj != this) // Don't check self-intersection
                    {
                        double shadowT = obj->intersect(lightRay, nullptr, 0);
                        if (shadowT > 0.001 && shadowT < lightDistance - 0.001)
                        {
                            inShadow = true;
                            break;
                        }
                    }
                }

                if (!inShadow)
                {
                    Vector3D toLight = sl->position - intersectionPoint;
                    toLight.normalize();
                    double lambertValue = std::max(0.0, normal.dot(toLight));

                    Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                    reflected.normalize();
                    Vector3D toCamera = r.start - intersectionPoint;
                    toCamera.normalize();
                    double phongValue = std::max(0.0, reflected.dot(toCamera));
                    phongValue = pow(phongValue, shine);

                    for (int i = 0; i < 3; i++)
                    {
                        color[i] += sl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i];
                        color[i] += sl->color[i] * coefficients[2] * phongValue;
                    }
                }
            }
        }

        // Handle reflection
        if (level < maxRecursionLevel && coefficients[3] > 0) // coefficients[3] is reflection coefficient
        {
            // Calculate reflected ray direction
            Vector3D incident = r.dir;
            Vector3D reflected = incident - normal * (2.0 * incident.dot(normal));
            reflected.normalize();
            
            // Create reflected ray starting slightly forward from intersection point to avoid self-intersection
            Vector3D reflectedStart = intersectionPoint + reflected * 0.001;
            Ray reflectedRay(reflectedStart, reflected);
            
            // Find nearest intersecting object for the reflected ray
            double minT = std::numeric_limits<double>::max();
            Object *nearestObject = nullptr;
            
            for (Object *obj : objects)
            {
                double t = obj->intersect(reflectedRay, nullptr, 0);
                if (t > 0 && t < minT)
                {
                    minT = t;
                    nearestObject = obj;
                }
            }
            
            // If reflection ray hits an object, calculate reflected color
            if (nearestObject != nullptr)
            {
                double reflectedColor[3] = {0, 0, 0};
                nearestObject->intersect(reflectedRay, reflectedColor, level + 1);
                
                // Add reflected color contribution
                for (int i = 0; i < 3; i++)
                {
                    color[i] += reflectedColor[i] * coefficients[3]; // coefficients[3] is reflection coefficient
                }
            }
        }

        return t;
    }

    void getNormal(const Vector3D &point, Vector3D &normal) const override
    {
        // Gradient of quadric surface gives normal
        normal.x = 2 * coeffs[0] * point.x + coeffs[3] * point.y + coeffs[4] * point.z + coeffs[6];
        normal.y = 2 * coeffs[1] * point.y + coeffs[3] * point.x + coeffs[5] * point.z + coeffs[7];
        normal.z = 2 * coeffs[2] * point.z + coeffs[4] * point.x + coeffs[5] * point.y + coeffs[8];
        normal.normalize();
    }

private:
    bool isWithinReferenceCube(const Vector3D &point) const
    {
        // Check if point is within the reference cube bounds
        if (length > 0 && fabs(point.x - reference_point.x) > length / 2) return false;
        if (width > 0 && fabs(point.y - reference_point.y) > width / 2) return false;
        if (height > 0 && fabs(point.z - reference_point.z) > height / 2) return false;
        return true;
    }
};

#endif // SCENE_OBJECTS_HPP