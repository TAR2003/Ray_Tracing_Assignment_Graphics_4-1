#ifndef SCENE_OBJECTS_HPP
#define SCENE_OBJECTS_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
using namespace std;

// including the image file for loading the texture
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// OpenGL GLUT headers
#ifndef __APPLE__
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

// Forward Declarations of different objects
class Object;
class PointLight;
class SpotLight;

// vector3D class for 3D vector operations
class Vector3D
{
public:
    double x, y, z; // coordinates of the vector

    // the constructor initializes the vector with given coordinates
    Vector3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

    // the operator+ function adds two vectors
    Vector3D operator+(const Vector3D &v) const
    {
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }
    // the operator- function subtracts two vectors
    Vector3D operator-(const Vector3D &v) const
    {
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }
    // the operator* function multiplies a vector by a scalar
    Vector3D operator*(double s) const
    {
        return Vector3D(x * s, y * s, z * s);
    }
    // the dot function calculates the dot product of two vectors
    double dot(const Vector3D &v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }
    // the cross function calculates the cross product of two vectors
    Vector3D cross(const Vector3D &v) const
    {
        return Vector3D(y * v.z - z * v.y,
                        z * v.x - x * v.z,
                        x * v.y - y * v.x);
    }
    // the length function calculates the length of the vector
    double length() const
    {
        return sqrt(x * x + y * y + z * z);
    }
    // the normalize function normalizes the vector to unit length
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
    // the print function prints the vector coordinates
    void print() const
    {
        cout << "Vector3D(" << x << ", " << y << ", " << z << ")" << endl;
    }
};

class Color
{
public:
    double r, g, b; // color components red, green, blue
    // constructor initializes the color with given components
    Color(double r = 0, double g = 0, double b = 0) : r(r), g(g), b(b) {}

    /// @brief Convert the color to an array
    /// @param array The array to store the color components
    void toArray(double *array) const
    {
        array[0] = r;
        array[1] = g;
        array[2] = b;
    }
    // function to print the color components
    void print() const
    {
        cout << "Color(" << r << ", " << g << ", " << b << ")" << endl;
    }
};

// Texture loading and sampling utilities for the scene
class TextureManager
{
public:
    static unsigned char *textureData;                       // Pointer to texture data
    static int textureWidth, textureHeight, textureChannels; // Texture dimensions and channels
    static bool textureLoaded;                               // Flag to check if texture is loaded

    /// @brief Load a texture from file
    /// @param filename The name of the texture file
    /// @return True if the texture was loaded successfully, false otherwise
    static bool loadTexture(const char *filename)
    {
        cleanup();                                                                             // Clean up previous texture if any remained after the previous stages or loading times
        textureData = stbi_load(filename, &textureWidth, &textureHeight, &textureChannels, 0); // Load the texture using stb_image
        if (!textureData)                                                                      // Check if texture loading failed
        {
            cerr << "Failed to load texture: " << filename << endl; // Print error message
            return false;
        }
        textureLoaded = true; // Set the texture loaded flag to true
        return true;
    }

    /// @brief Sample the texture at given (u, v) coordinates
    /// @param u The horizontal texture coordinate
    /// @param v The vertical texture coordinate
    /// @return The color sampled from the texture
    static Color sampleTexture(double u, double v)
    {

        if (!textureLoaded)
            return Color(1.0, 0.0, 1.0); // Return magenta if no texture loaded
        // Ensure u and v are in the range [0, 1]
        u = max(0.0, min(1.0, u)); // Clamp u to [0, 1]
        v = max(0.0, min(1.0, v)); // Clamp v to [0, 1]

        // Check if texture data is valid
        // and dimensions are positive
        if (!textureData || textureWidth <= 0 || textureHeight <= 0)
        {
            return Color(0.5, 0.5, 0.5); // Gray fallback
        }
        int pixel_x = (int)(u * (textureWidth - 1)); // Convert u to pixel x coordinate
        // Flip v to match OpenGL texture coordinate system (0 at bottom)
        int pixel_y = (int)((1.0 - v) * (textureHeight - 1)); // Flipping y

        // Safety clamp to ensure pixel coordinates are within bounds
        pixel_x = max(0, min(textureWidth - 1, pixel_x));
        pixel_y = max(0, min(textureHeight - 1, pixel_y));

        // Compute array index for the pixel in the texture data
        // Each pixel has textureChannels components (e.g., RGB or RGBA)
        int index = (pixel_y * textureWidth + pixel_x) * textureChannels; // Calculate index in the texture data array
        int max_index = textureWidth * textureHeight * textureChannels;   // Maximum index in the texture data array
        if (index < 0 || index + 2 >= max_index)                          // Check if index is out of bounds
            return Color(1.0, 0.0, 1.0);                                  // Magenta = error

        Color color;                          // Create a Color object to store the sampled color
        color.r = textureData[index] / 255.0; // Normalize the red component (0-255 to 0.0-1.0)

        if (textureChannels >= 2)
            color.g = textureData[index + 1] / 255.0; // Normalize the green component
        else
            color.g = color.r; // Grayscale

        if (textureChannels >= 3)
            color.b = textureData[index + 2] / 255.0; // Normalize the blue component
        else
            color.b = color.r; // Grayscale
        // every time grayscale is being used as a fallback, it is set to the red component

        return color;
    }

    /// @brief Clean up the texture data
    /// Frees the texture data if it was loaded and sets the pointer to nullptr
    static void cleanup()
    {
        if (textureData) // Check if texture data is loaded
        {
            stbi_image_free(textureData); // Free the texture data using stb_image
            textureData = nullptr;        // Set pointer to nullptr to avoid dangling pointer
            textureLoaded = false;        // Set texture loaded flag to false
        }
    }
};

class Ray
{
public:
    Vector3D start; // Starting point of the ray
    Vector3D dir;   // Direction of the ray

    /// @brief Constructor to initialize a ray with a starting point and direction
    /// @param start The starting point of the ray
    /// @param dir The direction of the ray
    Ray(const Vector3D &start, const Vector3D &dir) : start(start), dir(dir)
    {
        this->dir.normalize(); // Normalize the direction vector
    }
};

class PointLight
{
public:
    Vector3D position; // Position of the light in 3D space
    double color[3];   // Color components of the light (RGB)

    /// @brief Constructor to initialize a point light
    /// @param pos The position of the light
    /// @param r The red component of the light color
    /// @param g The green component of the light color
    /// @param b The blue component of the light color
    PointLight(const Vector3D &pos, double r, double g, double b) : position(pos)
    {
        color[0] = r; // Set red component
        color[1] = g; // Set green component
        color[2] = b; // Set blue component
    }

    /// @brief Draw the point light as a small sphere in the scene at the position variable
    /// This function uses OpenGL to visualize the light in the scene.
    void draw()
    {
        glPushMatrix();                                   // Save the current matrix state
        glTranslatef(position.x, position.y, position.z); // Move to the light position
        glColor3f(color[0], color[1], color[2]);          // Set the light color using the array color
        glutSolidSphere(0.2, 10, 10);                     // Draw a small sphere to represent the light as kind a light bulb
        glPopMatrix();                                    // Restore the previous matrix state before we changed the the current state to draw the light
    }
};

class SpotLight : public PointLight
{
public:
    Vector3D direction; // Direction of the spotlight beam
    double cutoffAngle; // Cutoff angle for the spotlight

    /// @brief Constructor to initialize a spotlight
    /// @param pos The position of the spotlight
    /// @param dir The direction of the spotlight beam
    /// @param angle The cutoff angle for the spotlight
    /// @param r The red component of the light color
    /// @param g The green component of the light color
    /// @param b The blue component of the light color
    SpotLight(const Vector3D &pos, const Vector3D &dir, double angle,
              double r, double g, double b)
        : PointLight(pos, r, g, b), direction(dir), cutoffAngle(angle)
    {
        direction.normalize(); // Normalize the direction vector
    }

    /// @brief Draw the spotlight as a point light and visualize its direction
    void draw()
    {
        PointLight::draw(); // Draw the point light representation
    }
};

// External References to the object and light array
extern vector<PointLight *> pointLights; // Array of point lights in the scene
extern vector<SpotLight *> spotLights;   // Array of spotlights in the scene
extern vector<Object *> objects;         // Array of objects in the scene
extern int maxRecursionLevel;            // Maximum recursion level for reflections

class Object
{
public:
    Vector3D reference_point;     // Reference point of the object in 3D space
    double height, width, length; // Dimensions of the object
    double color[3];              // Color components of the object (RGB)
    double coefficients[4];       // ambient, diffuse, specular, reflection coefficients
    int shine;                    // Shininess factor for specular highlights

    Object()
    {
        height = width = length = 0.0;                                               // Default dimensions of the object
        color[0] = color[1] = color[2] = 0.0;                                        // Default color is black
        coefficients[0] = coefficients[1] = coefficients[2] = coefficients[3] = 0.0; // Default coefficients
        shine = 0;                                                                   // Default shininess factor
    }

    virtual ~Object() {}                                                        // virtual destructor to allow proper cleanup of derived classes
    virtual void draw() = 0;                                                    // Pure virtual function to draw the object
    virtual double intersect(const Ray &r, double *color, int level) const = 0; // Pure virtual function to check intersection with a ray
    virtual void getNormal(const Vector3D &point, Vector3D &normal) const = 0;  // Pure virtual function to get the normal at a point on the object

    void setColor(double r, double g, double b)
    {
        color[0] = r; // Set red component
        color[1] = g; // Set green component
        color[2] = b; // Set blue component
    }

    void setShine(int sh)
    {
        shine = sh; // Set the shininess factor
    }

    /// @brief Set the coefficients for ambient, diffuse, specular, and reflection
    void setCoefficients(double ambient, double diffuse, double specular, double reflection)
    {
        coefficients[0] = ambient;
        coefficients[1] = diffuse;
        coefficients[2] = specular;
        coefficients[3] = reflection;
    }

    /// @brief Get the color of the object at a specific point
    /// @param point The point at which to get the color
    /// @param colorAt The array to store the color components (RGB)
    virtual void getColorAt(const Vector3D &point, double *colorAt) const
    {
        // Default implementation returns the object's color
        colorAt[0] = color[0];
        colorAt[1] = color[1];
        colorAt[2] = color[2];
    }
};

class Sphere : public Object
{
public:
    double radius; // Radius of the sphere

    Sphere(const Vector3D &center, double r) : radius(r)
    {
        reference_point = center;        // Set the reference point to the center of the sphere
        height = width = length = r * 2; // Set dimensions based on radius
    }
    void draw() override
    {
        glPushMatrix();                                                        // Save the current matrix state
        glTranslatef(reference_point.x, reference_point.y, reference_point.z); // Move to the sphere's center
        glColor3f(color[0], color[1], color[2]);                               // Set the sphere's color
        glutSolidSphere(radius, 50, 50);                                       // Draw the sphere with a solid appearance
        glPopMatrix();                                                         // Restore the previous matrix state
    }

    // @brief get the intersection point of the ray with the sphere
    double intersect(const Ray &r, double *color, int level) const override
    {
        Vector3D oc = r.start - reference_point; // Vector from ray start to sphere center
        double a = r.dir.dot(r.dir);             // Dot product of ray direction with itself
        double b = 2.0 * oc.dot(r.dir);          // Dot product of oc with ray direction, multiplied by 2
        double c = oc.dot(oc) - radius * radius; // Dot product of oc with itself minus radius squared
        double discriminant = b * b - 4 * a * c; // Discriminant for quadratic equation
        if (discriminant < 0)
            return -1.0;

        double sqrt_discriminant = sqrt(discriminant);    // Square root of the discriminant
        double t1 = (-b - sqrt_discriminant) / (2.0 * a); // First root of the quadratic equation
        double t2 = (-b + sqrt_discriminant) / (2.0 * a); // Second root of the quadratic equation

        // we will take the value of the smallest t possible
        double t = (t1 > 0) ? t1 : t2;
        if (t <= 0)
            return -1.0; // if no t root is possible we will return -1.0

        // If level is 0, only return intersection distance
        if (level == 0)
            return t;

        // Calculate intersection point and lighting
        Vector3D intersectionPoint = r.start + r.dir * t;
        Vector3D normal;
        getNormal(intersectionPoint, normal); // Get the normal at the intersection point

        // Get color at intersection point
        double intersectionPointColor[3];
        getColorAt(intersectionPoint, intersectionPointColor);

        // Initialize with ambient component of the object
        color[0] = intersectionPointColor[0] * coefficients[0];
        color[1] = intersectionPointColor[1] * coefficients[0];
        color[2] = intersectionPointColor[2] * coefficients[0];

        // Process each point light
        for (PointLight *pl : pointLights)
        {
            // Cast ray from light to intersection point
            Vector3D lightDir = intersectionPoint - pl->position; // the direction from the light position to the intersection point, used vector difference operator overloading
            double lightDistance = lightDir.length();
            lightDir.normalize();
            Ray lightRay(pl->position, lightDir); // Create a ray from the light to the intersection point

            // Check if intersection point is in shadow
            bool inShadow = false;
            for (Object *obj : objects)
            {
                // in this step we are checking any of the objects in the scene stand middle of the light and the intersection point and on the ray line
                // if so, then the ray can not reach the intersection point, thus it'll be shadowed
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
                // Calculate diffuse component Lambert
                // represents how much light hits the surface based on its angle to the light source
                Vector3D toLight = pl->position - intersectionPoint;
                toLight.normalize();
                double lambertValue = max(0.0, normal.dot(toLight));

                // Calculate specular component Phong
                // represents the shiny highlight on the surface based on the angle to the light source and viewer
                Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                reflected.normalize();
                Vector3D toCamera = r.start - intersectionPoint;
                toCamera.normalize();
                double phongValue = max(0.0, reflected.dot(toCamera));
                phongValue = pow(phongValue, shine);

                // Add diffuse and specular contributions
                // coefficients[1] is diffuse coefficient, coefficients[2] is specular coefficient
                for (int i = 0; i < 3; i++)
                {
                    color[i] += pl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i]; // diffuse
                    color[i] += pl->color[i] * coefficients[2] * phongValue;                               // specular
                }
            }
        }

        // Process each spot light
        for (SpotLight *sl : spotLights)
        {
            // Spotlight cone check
            Vector3D lightToPoint = intersectionPoint - sl->position; // Direction from spotlight position to intersection point
            double lightDistance = lightToPoint.length();
            lightToPoint.normalize(); // Normalize the direction vector

            // Check if point is within spotlight cone
            double angle = acos(lightToPoint.dot(sl->direction));
            if (angle <= sl->cutoffAngle)
            {
                // Create ray from spotlight to intersection point
                // This ray will be used to check if the intersection point is in shadow
                Ray lightRay(sl->position, lightToPoint);

                // Check shadow
                bool inShadow = false;
                for (Object *obj : objects)
                {
                    // to check if there is any objects between the spotlight and the intersection point
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
                    Vector3D toLight = sl->position - intersectionPoint; // direction from intersection point to spotlight position
                    toLight.normalize();                                 // Normalize the direction vector
                    double lambertValue = max(0.0, normal.dot(toLight)); // Lambert's cosine law for diffuse reflection

                    // Calculate specular component using Phong reflection model
                    Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                    reflected.normalize();
                    Vector3D toCamera = r.start - intersectionPoint;
                    toCamera.normalize();
                    double phongValue = max(0.0, reflected.dot(toCamera));
                    phongValue = pow(phongValue, shine);

                    for (int i = 0; i < 3; i++)
                    {
                        color[i] += sl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i];
                        color[i] += sl->color[i] * coefficients[2] * phongValue;
                    }
                }
            }
        }

        // Handle reflection of the ray if recursion level allows
        if (level < maxRecursionLevel && coefficients[3] > 0) // coefficients[3] is reflection coefficient
        {
            // Calculate reflected ray direction
            Vector3D incident = r.dir;                                             // the incoming ray direction
            Vector3D reflected = incident - normal * (2.0 * incident.dot(normal)); // computed using the reflection formula
            reflected.normalize();

            // Create reflected ray starting slightly forward from intersection point to avoid self-intersection
            Vector3D reflectedStart = intersectionPoint + reflected * 0.001;
            Ray reflectedRay(reflectedStart, reflected);

            // tracks the smallest intersection distance t along the reflected ray
            double minT = numeric_limits<double>::max();
            Object *nearestObject = nullptr;

            for (Object *obj : objects)
            {
                // tries to check for all objects which one is the closest object in the ray
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
                double reflectedColor[3] = {0, 0, 0};                              // init color array for reflected color
                nearestObject->intersect(reflectedRay, reflectedColor, level + 1); // traces the reflected ray firther with incrementing the level

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
    double floorWidth; // Width of the floor
    double tileWidth;  // Width of each tile in the checkerboard pattern
    bool useTexture;   // Flag to switch between checkerboard and texture

    /// @brief Constructor for the Floor class
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
            // there is no texture loaded or we are in checkerboard mode
            // Use checkerboard pattern
            int tx = static_cast<int>((point.x + floorWidth / 2) / tileWidth);
            int ty = static_cast<int>((point.y + floorWidth / 2) / tileWidth);

            // straight forward checkerboard pattern alternating black and white tiles
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

    // drawing the floor as a checkerboard or texture
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

    /// @brief Intersect a ray with the floor
    /// @param r The ray to intersect
    /// @param color The color at the intersection point
    /// @param level The recursion level
    /// @return The distance to the intersection point, or -1 if no intersection
    double intersect(const Ray &r, double *color, int level) const override
    {
        // NOw we will first test if the ray is parallel to the floor
        // if so, we will return -1.0 to indicate no intersection
        if (fabs(r.dir.z) < 1e-6)
            return -1.0;

        double t = -r.start.z / r.dir.z;
        if (t <= 0)
            return -1.0;

        // checks with the intersectiojn distance t if the intersection point is within the floor bounds
        // rejects if the intersection point is outside the floor width
        Vector3D intersection = r.start + r.dir * t;
        if (fabs(intersection.x) > floorWidth / 2 || fabs(intersection.y) > floorWidth / 2)
        {
            return -1.0;
        }

        if (level == 0)
            return t;

        // Get color at intersection point for both checkerboard and texture
        double intersectionPointColor[3];
        getColorAt(intersection, intersectionPointColor);

        // Initialize with ambient component
        color[0] = intersectionPointColor[0] * coefficients[0];
        color[1] = intersectionPointColor[1] * coefficients[0];
        color[2] = intersectionPointColor[2] * coefficients[0];

        Vector3D normal = Vector3D(0, 0, 1); // Floor normal always points up, this is a must

        // Process each point light
        for (PointLight *pl : pointLights)
        {
            // for every point light, we will cast a ray from the light to the intersection point
            Vector3D lightDir = intersection - pl->position;
            double lightDistance = lightDir.length();
            lightDir.normalize();
            Ray lightRay(pl->position, lightDir);

            // Check if intersection point is in shadow
            bool inShadow = false;
            // Now shadow checking, if there is any object between the light and the intersection point
            // if so, then the intersection point is shadowed
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
                // Calculate diffuse component Lambert
                // represents how much light hits the surface based on its angle to the light source
                Vector3D toLight = pl->position - intersection;
                toLight.normalize();
                double lambertValue = max(0.0, normal.dot(toLight));

                // Calculate specular component Phong
                // represents the shiny highlight on the surface based on the angle to the light source and viewer
                Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                reflected.normalize();
                Vector3D toCamera = r.start - intersection;
                toCamera.normalize();
                double phongValue = max(0.0, reflected.dot(toCamera));
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
            // Check if intersection point is within spotlight cone
            // Cast ray from spotlight to intersection point
            Vector3D lightToPoint = intersection - sl->position;
            double lightDistance = lightToPoint.length();
            lightToPoint.normalize();

            double angle = acos(lightToPoint.dot(sl->direction)); // Calculate angle between light direction and vector to intersection point
            // Check if the angle is within the spotlight's cutoff angle

            if (angle <= sl->cutoffAngle)
            {
                // the angle is less or equal  to the cutoff angle, so we will proceed with the spotlight
                // Create ray from spotlight to intersection point
                // This ray will be used to check if the intersection point is in shadow
                Ray lightRay(sl->position, lightToPoint);

                bool inShadow = false;
                // we will check if there is any object between the spotlight and the intersection point
                // if so, then the intersection point is shadowed
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
                    // Calculate diffuse component Lambert
                    // represents how much light hits the surface based on its angle to the light source
                    Vector3D toLight = sl->position - intersection;
                    toLight.normalize();
                    double lambertValue = max(0.0, normal.dot(toLight));

                    // Calculate specular component using Phong reflection model
                    // represents the shiny highlight on the surface based on the angle to the light source and viewer
                    Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                    reflected.normalize();
                    Vector3D toCamera = r.start - intersection;
                    toCamera.normalize();
                    double phongValue = max(0.0, reflected.dot(toCamera));
                    phongValue = pow(phongValue, shine);

                    for (int i = 0; i < 3; i++)
                    {
                        color[i] += sl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i];
                        color[i] += sl->color[i] * coefficients[2] * phongValue;
                    }
                }
            }
        }

        // Handle reflection of the ray if recursion level allows
        if (level < maxRecursionLevel && coefficients[3] > 0) // coefficients[3] is reflection coefficient
        {
            // Calculate reflected ray direction
            Vector3D incident = r.dir; // the incoming ray direction
            Vector3D reflected = incident - normal * (2.0 * incident.dot(normal)); // computed using the reflection formula 
            reflected.normalize();

            // Create reflected ray starting slightly forward from intersection point to avoid self-intersection
            Vector3D reflectedStart = intersection + reflected * 0.001;
            Ray reflectedRay(reflectedStart, reflected); // now we  created a new Ray object which represents the reflected ray

            // Find nearest intersecting object for the reflected ray
            double minT = numeric_limits<double>::max(); // Initialize minimum intersection distance
            Object *nearestObject = nullptr;

            // Now we will iterate over all objects in the scene 
            // and check which one is the closest object in the reflected ray
            // and we will store the minimum intersection distance t
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
                // nearestObject is the closest object in the reflected ray
                double reflectedColor[3] = {0, 0, 0};
                nearestObject->intersect(reflectedRay, reflectedColor, level + 1); // we will increase the recursion level by 1

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

    /// @brief draw the triangle using OpenGL
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

    /// @brief Intersect a ray with the triangle
    /// @param r The ray to intersect
    /// @param color The color at the intersection point
    /// @param level The recursion level
    /// @return The distance to the intersection point, or -1 if no intersection
    double intersect(const Ray &r, double *color, int level) const override
    {
        // Ray-triangle intersection 
        Vector3D edge1 = b - a;
        Vector3D edge2 = c - a;
        Vector3D h = r.dir.cross(edge2);
        double a_det = edge1.dot(h);

        // we will check if the ray is parallel to the triangle
        // if so, we will return -1.0 to indicate no intersection
        if (a_det > -1e-6 && a_det < 1e-6)
            return -1.0; // Ray is parallel to triangle

        // Barycentric coordinates calculation
        // If the ray is not parallel, we will proceed with the intersection calculation
        // use early rejection technique to avoid unnecessary calculations
        double f = 1.0 / a_det;
        Vector3D s = r.start - a;
        double u = f * s.dot(h);


        if (u < 0.0 || u > 1.0)
            return -1.0; // if the barycentric coordinate u is outside the triangle

        Vector3D q = s.cross(edge1);
        double v = f * r.dir.dot(q);

        if (v < 0.0 || u + v > 1.0)
            return -1.0; // if the barycentric coordinate v is outside the triangle or u + v exceeds 1

        double t = f * edge2.dot(q);

        // Check if the intersection point is in front of the ray start
        // if t is less than or equal to 0, then the intersection point is behind, so we dont want it
        if (t > 1e-6) // Ray intersection , here 1e-6 is used to avoid double value problems
        {
            if (level == 0)
                return t; // level is 0, so there is no need to calculate color or lighting

            // Calculate intersection point and lighting
            Vector3D intersectionPoint = r.start + r.dir * t; // we compute the intersection point using the ray start and direction
            // Get the normal at the intersection point
            Vector3D normal;
            getNormal(intersectionPoint, normal);

            // Get color at intersection point
            double intersectionPointColor[3];
            getColorAt(intersectionPoint, intersectionPointColor);

            // scales the base color intersection point color by ambient coefficient
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
                    // Check if any object is between the light and the intersection point
                    // if so, then the intersection point is shadowed
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
                    // represents how much light hits the surface based on its angle to the light source
                    Vector3D toLight = pl->position - intersectionPoint;
                    toLight.normalize();
                    double lambertValue = max(0.0, normal.dot(toLight));

                    // Calculate specular component (Phong)
                    // represents the shiny highlight on the surface based on the angle to the light source and viewer
                    // reflected vector is calculated using the normal and the light direction
                    Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                    reflected.normalize();
                    Vector3D toCamera = r.start - intersectionPoint;
                    toCamera.normalize();
                    double phongValue = max(0.0, reflected.dot(toCamera));
                    phongValue = pow(phongValue, shine);

                    // Add diffuse and specular contributions
                    for (int i = 0; i < 3; i++)
                    {
                        color[i] += pl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i]; // diffuse
                        color[i] += pl->color[i] * coefficients[2] * phongValue;                               // specular
                    }
                }
            }

            // Process each spot light
            for (SpotLight *sl : spotLights)
            {
                // Check if intersection point is within spotlight cone
                // Cast ray from spotlight to intersection point
                Vector3D lightToPoint = intersectionPoint - sl->position;
                double lightDistance = lightToPoint.length();
                lightToPoint.normalize();

                // Check if point is within spotlight cone
                double angle = acos(lightToPoint.dot(sl->direction));
                if (angle <= sl->cutoffAngle)
                {
                    // the angle is less or equal to the cutoff angle, so we will proceed with the spotlight
                    // Create ray from spotlight to intersection point
                    Ray lightRay(sl->position, lightToPoint);

                    // Check shadow
                    bool inShadow = false;
                    for (Object *obj : objects)
                    {
                        // we now check all objects in the scene to see if any of them is between the spotlight and the intersection point
                        // if so, then the intersection point is shadowed
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
                        // Calculate diffuse component Lambert
                        // represents how much light hits the surface based on its angle to the light source
                        Vector3D toLight = sl->position - intersectionPoint;
                        toLight.normalize();
                        double lambertValue = max(0.0, normal.dot(toLight));

                        // Calculate specular component using Phong reflection model
                        // represents the shiny highlight on the surface based on the angle to the light source and viewer
                        Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                        reflected.normalize();
                        Vector3D toCamera = r.start - intersectionPoint;
                        toCamera.normalize();
                        double phongValue = max(0.0, reflected.dot(toCamera));
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
                Vector3D incident = r.dir; // Incident ray direction
                Vector3D reflected = incident - normal * (2.0 * incident.dot(normal));
                reflected.normalize();

                // Create reflected ray starting slightly forward from intersection point to avoid self-intersection
                Vector3D reflectedStart = intersectionPoint + reflected * 0.001;
                Ray reflectedRay(reflectedStart, reflected);

                // Find nearest intersecting object for the reflected ray
                double minT = numeric_limits<double>::max();
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
                    // nearestObject is the closest object in the reflected ray
                    // Initialize reflected color array
                    double reflectedColor[3] = {0, 0, 0};
                    nearestObject->intersect(reflectedRay, reflectedColor, level + 1); // we will increment the recursion level by 1

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
        // meaning the ray does not intersect the triangle in the positive direction, so any reflection calculation will not be needed
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
    double all_coefficients[10]; // A B C D E F G H I J coefficients

    GeneralQuadric(double A, double B, double C, double D, double E,
                   double F, double G, double H, double I, double J)
    {
        all_coefficients[0] = A;
        all_coefficients[1] = B;
        all_coefficients[2] = C;
        all_coefficients[3] = D;
        all_coefficients[4] = E;
        all_coefficients[5] = F;
        all_coefficients[6] = G;
        all_coefficients[7] = H;
        all_coefficients[8] = I;
        all_coefficients[9] = J;
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

        Vector3D o = r.start; // ray origin
        Vector3D d = r.dir; // ray direction

        // substitutes the ray equation into the quadric equation to derive the quadratic coefficients
        // a, b, c are coefficients of the quadratic equation
        double a = all_coefficients[0] * d.x * d.x + all_coefficients[1] * d.y * d.y + all_coefficients[2] * d.z * d.z +
                   all_coefficients[3] * d.x * d.y + all_coefficients[4] * d.x * d.z + all_coefficients[5] * d.y * d.z;

        double b = 2 * all_coefficients[0] * o.x * d.x + 2 * all_coefficients[1] * o.y * d.y + 2 * all_coefficients[2] * o.z * d.z +
                   all_coefficients[3] * (o.x * d.y + o.y * d.x) + all_coefficients[4] * (o.x * d.z + o.z * d.x) +
                   all_coefficients[5] * (o.y * d.z + o.z * d.y) + all_coefficients[6] * d.x + all_coefficients[7] * d.y + all_coefficients[8] * d.z;

        double c = all_coefficients[0] * o.x * o.x + all_coefficients[1] * o.y * o.y + all_coefficients[2] * o.z * o.z +
                   all_coefficients[3] * o.x * o.y + all_coefficients[4] * o.x * o.z + all_coefficients[5] * o.y * o.z +
                   all_coefficients[6] * o.x + all_coefficients[7] * o.y + all_coefficients[8] * o.z + all_coefficients[9];

        double discriminant = b * b - 4 * a * c;

        // check the discriminant to determine if there are real solutions
        if (discriminant < 0)
            return -1.0;

        double sqrt_discriminant = sqrt(discriminant);
        double t1 = (-b - sqrt_discriminant) / (2.0 * a);
        double t2 = (-b + sqrt_discriminant) / (2.0 * a);

        // Choose the smallest positive t value, if both are negative return -1
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
            Vector3D lightDir = intersectionPoint - pl->position; // Light direction
            double lightDistance = lightDir.length();
            lightDir.normalize();
            Ray lightRay(pl->position, lightDir);

            // Check if intersection point is in shadow
            bool inShadow = false;
            for (Object *obj : objects)
            {
                // Check if any object is between the light and the intersection point
                // if so, then the intersection point is shadowed
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
                double lambertValue = max(0.0, normal.dot(toLight));

                // Calculate specular component (Phong)
                Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                reflected.normalize();
                Vector3D toCamera = r.start - intersectionPoint;
                toCamera.normalize();
                double phongValue = max(0.0, reflected.dot(toCamera));
                phongValue = pow(phongValue, shine);

                // Add diffuse and specular contributions
                for (int i = 0; i < 3; i++)
                {
                    color[i] += pl->color[i] * coefficients[1] * lambertValue * intersectionPointColor[i]; // diffuse
                    color[i] += pl->color[i] * coefficients[2] * phongValue;                               // specular
                }
            }
        }

        // Process each spot light
        for (SpotLight *sl : spotLights)
        {
            // Check if intersection point is within spotlight cone
            // Cast ray from spotlight to intersection point
            Vector3D lightToPoint = intersectionPoint - sl->position;
            double lightDistance = lightToPoint.length();
            lightToPoint.normalize();

            // Check if point is within spotlight cone
            double angle = acos(lightToPoint.dot(sl->direction));
            // If the angle is less than or equal to the cutoff angle, we proceed with the spotlight
            // Create ray from spotlight to intersection point
            if (angle <= sl->cutoffAngle)
            {
                Ray lightRay(sl->position, lightToPoint);

                // Check shadow
                bool inShadow = false;
                for (Object *obj : objects)
                {
                    // Check if any object is between the spotlight and the intersection point
                    // if so, then the intersection point is shadowed
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
                    // Calculate diffuse component (Lambert)
                    Vector3D toLight = sl->position - intersectionPoint;
                    toLight.normalize();
                    double lambertValue = max(0.0, normal.dot(toLight));

                    // Calculate specular component (Phong)
                    Vector3D reflected = normal * (2.0 * normal.dot(toLight)) - toLight;
                    reflected.normalize();
                    Vector3D toCamera = r.start - intersectionPoint;
                    toCamera.normalize();
                    double phongValue = max(0.0, reflected.dot(toCamera));
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
            double minT = numeric_limits<double>::max();
            Object *nearestObject = nullptr;

            // Find the nearest object intersecting the reflected ray
            for (Object *obj : objects)
            {
                // Check intersection with the reflected ray
                // if the object is not the current object, we will check for intersection
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
                // nearestObject is the closest object in the reflected ray
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
        normal.x = 2 * all_coefficients[0] * point.x + all_coefficients[3] * point.y + all_coefficients[4] * point.z + all_coefficients[6];
        normal.y = 2 * all_coefficients[1] * point.y + all_coefficients[3] * point.x + all_coefficients[5] * point.z + all_coefficients[7];
        normal.z = 2 * all_coefficients[2] * point.z + all_coefficients[4] * point.x + all_coefficients[5] * point.y + all_coefficients[8];
        normal.normalize();
    }

    bool isWithinReferenceCube(const Vector3D &point) const
    {
        // Check if point is within the reference cube bounds
        if (length > 0 && fabs(point.x - reference_point.x) > length / 2)
            return false;
        if (width > 0 && fabs(point.y - reference_point.y) > width / 2)
            return false;
        if (height > 0 && fabs(point.z - reference_point.z) > height / 2)
            return false;
        return true;
    }
};

#endif // SCENE_OBJECTS_HPP