/**
 * Ray Tracing Assignment - Complete Implementation
 * 
 * A unified ray tracing application with OpenGL visualization.
 * Features:
 * - Scene loading from file with support for spheres, triangles, quadrics
 * - Real-time camera movement and controls
 * - Ray tracing with Phong illumination, shadows, and reflections
 * - Texture mapping for floors with toggle functionality
 * - Integrated camera system with comprehensive movement controls
 * - OpenGL visualization of objects and lights
 * 
 * Controls:
 * - Arrow keys: Move forward/back/left/right
 * - Page Up/Down: Move up/down
 * - 1/2: Rotate look-at point left/right
 * - 3/4: Rotate camera up/down
 * - 5/6: Tilt camera
 * - w/s: Move camera up/down
 * - 0: Capture ray-traced image
 * - t/T: Toggle floor texture (checkerboard/texture)
 * - a: Toggle coordinate axes
 * - ESC: Exit
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <limits>
#include <cmath>
#include "bitmap_image.hpp"
#include "2005090_Header.hpp" // Object and ray tracing class definitions

// Define M_PI if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// OpenGL / GLUT Headers
#ifdef __APPLE__
#include <GLUT/glut.h> // Use GLUT framework on macOS
#else
#include <GL/glut.h> // Use standard GLUT location on Linux/Windows
#endif

#define M_PI 3.14159265358979323846

const char* textureFileName = "texture2.bmp"; // Default texture file name

/*
 * TEXTURE MAPPING FEATURE:
 * - Press 't' or 'T' to toggle between checkerboard and texture modes for the floor
 * - The system loads "texture_1.bmp" by default
 * - To use different textures, call switchTexture("filename.bmp") 
 * - Supports common image formats (BMP, PNG, JPG, TGA, etc.)
 * - Texture coordinates are automatically mapped from world coordinates
 * - Falls back to checkerboard pattern if texture loading fails
 * - Texture repeats 10 times across the floor for better detail
 */

// Global variables
std::vector<Object *> objects;
std::vector<PointLight *> pointLights;
std::vector<SpotLight *> spotLights;

// Recursion level for reflections
int maxRecursionLevel = 3;

// TextureManager static members
unsigned char* TextureManager::textureData = nullptr;
int TextureManager::textureWidth = 0;
int TextureManager::textureHeight = 0;
int TextureManager::textureChannels = 0;
bool TextureManager::textureLoaded = false;

// Camera position and orientation variables
GLfloat eyex = 100, eyey = -150, eyez = 60;    // Camera position coordinates - positioned for good scene overview
GLfloat centerx = 20, centery = 20, centerz = 20; // Look-at point coordinates - looking towards scene center
GLfloat upx = 0, upy = 0, upz = 1;             // Up vector coordinates - Z is up

float movementSpeed = 5.0f; // Speed of camera movement - increased for larger scene
float rotationSpeed = 0.05f; // Speed of camera rotation - slightly decreased for smoother control

// Object visibility flags
bool isAxes = true; // Toggle for coordinate axes

// Camera parameters
double fovY = 45.0;  // Better field of view for scene overview
double aspectRatio = 1.0;
double nearPlane = 1.0;
double farPlane = 1000.0;

// Window dimensions
int windowWidth = 1280; // Increased for better resolution
int windowHeight = 720;
double worldWidth = 2.0;
double worldHeight = 2.0;

/**
 * Initialize OpenGL settings
 * Sets up background color and enables depth testing
 */
void initGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);

    // Set up projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, static_cast<double>(windowWidth) / windowHeight, 0.1f, 1000.0f);
}

// Scene loading function
void loadData()
{
    std::ifstream file("scene.txt");
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open scene.txt file!" << std::endl;
        return;
    }

    int recursionLevel, imageSize;
    file >> recursionLevel >> imageSize; // the first pair of inputs are the recursion level and image size
    
    // Set global recursion level
    maxRecursionLevel = recursionLevel;
    
    std::cout << "Loading scene..." << std::endl;
    std::cout << "Recursion Level: " << recursionLevel << std::endl;
    std::cout << "Image Size: " << imageSize << std::endl;

    // Update global parameters based on scene file
    windowHeight = imageSize;
    windowWidth = imageSize;
    worldWidth = 200.0;   // Increase world dimensions to match scene scale
    worldHeight = 200.0;

    int numObjects;
    file >> numObjects;
    std::cout << "Number of objects to load: " << numObjects << std::endl;

    // Clear existing objects
    for (Object *obj : objects)
        delete obj;
    objects.clear();

    // Load objects
    for (int i = 0; i < numObjects; i++)
    {
        std::string objectType;
        file >> objectType;

        if (objectType == "sphere")
        {
            double centerX, centerY, centerZ, radius;
            double colorR, colorG, colorB;
            double ambient, diffuse, specular, reflection;
            int shininess;

            file >> centerX >> centerY >> centerZ;
            file >> radius;
            file >> colorR >> colorG >> colorB;
            file >> ambient >> diffuse >> specular >> reflection;
            file >> shininess;

            // now we create the sphere with all the variables
            Sphere *sphere = new Sphere(Vector3D(centerX, centerY, centerZ), radius);
            sphere->setColor(colorR, colorG, colorB);
            sphere->setCoefficients(ambient, diffuse, specular, reflection);
            sphere->setShine(shininess);
            objects.push_back(sphere);
            
            std::cout << "Loaded sphere at (" << centerX << ", " << centerY << ", " << centerZ 
                      << ") with radius " << radius << std::endl;
        }
        else if (objectType == "triangle")
        {
            double x1, y1, z1, x2, y2, z2, x3, y3, z3;
            double colorR, colorG, colorB;
            double ambient, diffuse, specular, reflection;
            int shininess;

            file >> x1 >> y1 >> z1;
            file >> x2 >> y2 >> z2;
            file >> x3 >> y3 >> z3;
            file >> colorR >> colorG >> colorB;
            file >> ambient >> diffuse >> specular >> reflection;
            file >> shininess;

            // now we create the triangle with all the variables
            Triangle *triangle = new Triangle(Vector3D(x1, y1, z1), 
                                            Vector3D(x2, y2, z2), 
                                            Vector3D(x3, y3, z3));
            triangle->setColor(colorR, colorG, colorB);
            triangle->setCoefficients(ambient, diffuse, specular, reflection);
            triangle->setShine(shininess);
            objects.push_back(triangle);
            
            std::cout << "Loaded triangle with vertices (" << x1 << "," << y1 << "," << z1 
                      << "), (" << x2 << "," << y2 << "," << z2 
                      << "), (" << x3 << "," << y3 << "," << z3 << ")" << std::endl;
        }
        else if (objectType == "general")
        {
            double A, B, C, D, E, F, G, H, I, J;
            double refX, refY, refZ, length, width, height;
            double colorR, colorG, colorB;
            double ambient, diffuse, specular, reflection;
            int shininess;

            file >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
            file >> refX >> refY >> refZ >> length >> width >> height;
            file >> colorR >> colorG >> colorB;
            file >> ambient >> diffuse >> specular >> reflection;
            file >> shininess;

            // now we create the general quadric with all the variables
            // GeneralQuadric is a class that represents a quadric surface defined by the coefficients A, B, C, D, E, F, G, H, I, J
            // and the reference point (refX, refY, refZ) with dimensions (
            GeneralQuadric *quadric = new GeneralQuadric(A, B, C, D, E, F, G, H, I, J);
            quadric->reference_point = Vector3D(refX, refY, refZ);
            quadric->length = length;
            quadric->width = width;
            quadric->height = height;
            quadric->setColor(colorR, colorG, colorB);
            quadric->setCoefficients(ambient, diffuse, specular, reflection);
            quadric->setShine(shininess);
            objects.push_back(quadric);
            
            std::cout << "Loaded general quadric at (" << refX << "," << refY << "," << refZ 
                      << ") with dimensions " << length << "x" << width << "x" << height << std::endl;
        }
    }

    // Read floor parameters from scene file
    // The scene.txt mentions floor should be 1000 width with 20 tile width
    double floorWidth = 1000.0;
    double tileWidth = 20.0;
    
    std::cout << "Creating floor with width " << floorWidth << " and tile size " << tileWidth << std::endl;
    
    // Add floor (as mentioned in the scene description)
    Floor *floor = new Floor(floorWidth, tileWidth, false); // Start with checkerboard mode
    objects.push_back(floor);
    
    std::cout << "Floor created with checkerboard pattern. Use 't' key to toggle texture mode." << std::endl;

    // Load point lights
    int numPointLights;
    file >> numPointLights;

    // clear any existing light pointers
    for (PointLight *light : pointLights)
        delete light;
    pointLights.clear();

    for (int i = 0; i < numPointLights; i++)
    {
        double x, y, z;
        double r, g, b;
        file >> x >> y >> z;
        file >> r >> g >> b;

        PointLight *light = new PointLight(Vector3D(x, y, z), r, g, b);
        pointLights.push_back(light);
        
        std::cout << "Loaded point light at (" << x << ", " << y << ", " << z 
                  << ") with color (" << r << ", " << g << ", " << b << ")" << std::endl;
    }

    // Load spot lights
    int numSpotLights;
    file >> numSpotLights;

    // Clear existing spot lights pointers
    for (SpotLight *light : spotLights)
        delete light;
    spotLights.clear();

    for (int i = 0; i < numSpotLights; i++)
    {
        double posX, posY, posZ;
        double colorR, colorG, colorB;
        double dirX, dirY, dirZ;
        double cutoffAngle;

        file >> posX >> posY >> posZ;
        file >> colorR >> colorG >> colorB;
        file >> dirX >> dirY >> dirZ;
        file >> cutoffAngle;

        // Convert cutoff angle from degrees to radians
        cutoffAngle = cutoffAngle * M_PI / 180.0;

        SpotLight *spotlight = new SpotLight(Vector3D(posX, posY, posZ),
                                           Vector3D(dirX, dirY, dirZ),
                                           cutoffAngle,
                                           colorR, colorG, colorB);
        spotLights.push_back(spotlight);
        
        std::cout << "Loaded spotlight at (" << posX << ", " << posY << ", " << posZ 
                  << ") pointing (" << dirX << ", " << dirY << ", " << dirZ 
                  << ") with cutoff " << cutoffAngle << std::endl;
    }

    file.close();
}

/**
 * Switches the floor texture to a new image file
 * @param filename The name of the texture file to load
 */
void switchTexture(const char* filename)
{
    if (TextureManager::loadTexture(filename)) {
        std::cout << "Switched to texture: " << filename << std::endl;
        // Enable texture mode on all floors
        for (Object *obj : objects) {
            Floor *floor = dynamic_cast<Floor*>(obj);
            if (floor) {
                floor->setTextureMode(true);
            }
        }
    } else {
        std::cout << "Failed to load texture: " << filename << std::endl;
    }
}

/// @brief Captures the current window content to an image file
void capture()
{
    static int captureCount = 11; // Start from Output_11.bmp as per the requirement
    const int imageWidth = windowWidth;  // Use actual window dimensions
    const int imageHeight = windowHeight;
    
    // Initialize bitmap image and set background color
    bitmap_image image(imageWidth, imageHeight);

    image.clear(); // Clear the image to default color
    
    // Set a light background color instead of pure black
    for (int i = 0; i < imageWidth; i++)
    {
        for (int j = 0; j < imageHeight; j++)
        {
            image.set_pixel(i, j, 50, 50, 80); // Dark blue background
        }
    }
    
    // Use camera variables from the header file
    Vector3D eye(eyex, eyey, eyez);
    Vector3D lookAt(centerx, centery, centerz);
    Vector3D up(upx, upy, upz);
    
    // Calculate view direction and orthonormal basis
    Vector3D l = lookAt - eye;
    l.normalize();
    Vector3D r = l.cross(up);
    r.normalize();
    Vector3D u = r.cross(l);
    
    // Calculate plane distance based on field of view
    double viewAngle = fovY * (M_PI / 180.0); // Convert to radians
    double planeDistance = (windowHeight / 2.0) / tan(viewAngle / 2.0);
    
    // Calculate top-left corner of the image plane
    Vector3D topLeft = eye + l * planeDistance - r * (windowWidth / 2.0) + u * (windowHeight / 2.0);
    
    // Calculate pixel step sizes
    // pixel height and width in world coordinates
    double du = (double)windowWidth / imageWidth; 
    double dv = (double)windowHeight / imageHeight;
    
    // Choose middle of the grid cell
    topLeft = topLeft + r * (0.5 * du) - u * (0.5 * dv); // adjust to center pixel
    
    std::cout << "Capturing image " << captureCount << "..." << std::endl;
    
    // Ray tracing loop
    for (int i = 0; i < imageWidth; i++)
    {
        for (int j = 0; j < imageHeight; j++)
        {
            // Calculate current pixel position
            Vector3D curPixel = topLeft + r * (i * du) - u * (j * dv);
            
            // Cast ray from eye to current pixel
            Vector3D rayDirection = curPixel - eye;
            Ray ray(eye, rayDirection);
            
            double tMin = std::numeric_limits<double>::max();
            Object *nearestObject = nullptr;
            double dummyColor[3] = {0, 0, 0};
            
            // Find nearest intersection
            for (Object *obj : objects)
            {
                // we will check t value for each objects in the scene
                double t = obj->intersect(ray, dummyColor, 0);
                if (t > 0 && t < tMin)
                {
                    tMin = t;
                    nearestObject = obj;
                }
            }
            
            // Calculate final color
            double *color = new double[3]{0.05, 0.05, 0.1}; // Darker background color
            if (nearestObject)
            {
                // nearestObject is not null, so there is an intersection
                // Reset color for object intersection
                color[0] = color[1] = color[2] = 0.0;
                // Get actual color from the nearest object
                double t_min = nearestObject->intersect(ray, color, 1);
            }
            
            // Clamp color values to [0, 1] range
            for (int k = 0; k < 3; k++)
            {
                color[k] = std::min(1.0, std::max(0.0, color[k]));
            }
            
            // Update image pixel (i,j)
            image.set_pixel(i, j,
                            static_cast<unsigned char>(color[0] * 255),
                            static_cast<unsigned char>(color[1] * 255),
                            static_cast<unsigned char>(color[2] * 255));
            
            delete[] color; // Clean up dynamic memory
        }
    }
    
    // Save image with naming convention: Output_11.bmp, Output_12.bmp, etc.
    std::string filename = "Output_" + std::to_string(captureCount++) + ".bmp";
    image.save_image(filename);
    std::cout << "Image saved as " << filename << std::endl;
}

/**
 * Draw coordinate axes
 * X axis: red, Y axis: green, Z axis: blue
 */
void drawAxes()
{
    glLineWidth(3); // Set line thickness

    glBegin(GL_LINES);

    // X axis (red)
    glColor3f(1, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(1, 0, 0);

    // Y axis (green)
    glColor3f(0, 1, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 1, 0);

    // Z axis (blue)
    glColor3f(0, 0, 1);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 1);

    glEnd();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Use camera variables from the header file
    gluLookAt(eyex, eyey, eyez,
              centerx, centery, centerz,
              upx, upy, upz);

    // Draw all objects
    for (Object *obj : objects)
    {
        obj->draw();
    }
    // Draw coordinate axes if enabled
    if (isAxes)
    {
        drawAxes();
    }


    // Draw all lights
    for (PointLight *light : pointLights)
    {
        light->draw();
    }

    for (SpotLight *light : spotLights)
    {
        light->draw();
    }

    glutSwapBuffers();
}

void reshape(int width, int height)
{
    windowWidth = width;
    windowHeight = height;
    glViewport(0, 0, width, height);
    aspectRatio = static_cast<double>(width) / height;
    initGL();
}

/**
 * Window reshape callback from camera settings
 * Handles window resizing and maintains aspect ratio
 */
void reshapeListener(GLsizei width, GLsizei height)
{
    // Prevent division by zero
    if (height == 0)
        height = 1;

    // Calculate aspect ratio
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set viewport to cover entire window
    glViewport(0, 0, width, height);

    // Set up perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // 45-degree field of view, aspect ratio, near and far clipping planes
    gluPerspective(45.0f, aspect, 0.1f, 1000.0f);
}

// Draw coordinate axes


void keyboard(unsigned char key, int x, int y)
{
    // Calculate view direction vector for camera operations
    double lx = centerx - eyex;
    double lz = centerz - eyez;
    double s;

    switch (key)
    {
    case '0':
        capture();
        break;
    case 't':
    case 'T':
        // Toggle texture mode for floor
        for (Object *obj : objects) {
            Floor *floor = dynamic_cast<Floor*>(obj);
            if (floor) {
                bool currentMode = floor->useTexture;
                floor->setTextureMode(!currentMode);
                std::cout << "Floor texture mode: " << (floor->useTexture ? "ON (texture)" : "OFF (checkerboard)") << std::endl;
                if (floor->useTexture && !TextureManager::textureLoaded) {
                    std::cout << "Warning: Texture mode enabled but no texture loaded!" << std::endl;
                }
                break; // Only toggle the first floor found
            }
        }
        break;
    case '1':
    {
        // Rotate look-at point left around camera (yaw left) - simplified algorithm like keys 3-6
        float ROTATION_ANGLE = rotationSpeed;
        
        // Calculate view direction vector
        float viewDirx = centerx - eyex;
        float viewDiry = centery - eyey;
        float viewDirz = centerz - eyez;
        
        // Get the view distance to maintain it
        float viewLength = sqrt(viewDirx * viewDirx + viewDiry * viewDiry + viewDirz * viewDirz);
        
        // Normalize the view direction
        viewDirx /= viewLength;
        viewDiry /= viewLength;
        viewDirz /= viewLength;
        
        // Normalize the up vector
        float upLength = sqrt(upx * upx + upy * upy + upz * upz);
        float normUpx = upx / upLength;
        float normUpy = upy / upLength;
        float normUpz = upz / upLength;
        
        // Rotation around up vector (yaw left)
        float cosAngle = cos(ROTATION_ANGLE);
        float sinAngle = sin(ROTATION_ANGLE);
        
        // Rotate view direction around up vector
        float newViewDirx = viewDirx * cosAngle + (normUpy * viewDirz - normUpz * viewDiry) * sinAngle;
        float newViewDiry = viewDiry * cosAngle + (normUpz * viewDirx - normUpx * viewDirz) * sinAngle;
        float newViewDirz = viewDirz * cosAngle + (normUpx * viewDiry - normUpy * viewDirx) * sinAngle;
        
        // Update center position
        centerx = eyex + newViewDirx * viewLength;
        centery = eyey + newViewDiry * viewLength;
        centerz = eyez + newViewDirz * viewLength;
        break;
    }
    case '2':
    {
        // Rotate look-at point right around camera (yaw right) - simplified algorithm like keys 3-6
        float ROTATION_ANGLE = -rotationSpeed; // Negative for opposite direction
        
        // Calculate view direction vector
        float viewDirx = centerx - eyex;
        float viewDiry = centery - eyey;
        float viewDirz = centerz - eyez;
        
        // Get the view distance to maintain it
        float viewLength = sqrt(viewDirx * viewDirx + viewDiry * viewDiry + viewDirz * viewDirz);
        
        // Normalize the view direction
        viewDirx /= viewLength;
        viewDiry /= viewLength;
        viewDirz /= viewLength;
        
        // Normalize the up vector
        float upLength = sqrt(upx * upx + upy * upy + upz * upz);
        float normUpx = upx / upLength;
        float normUpy = upy / upLength;
        float normUpz = upz / upLength;
        
        // Rotation around up vector (yaw right)
        float cosAngle = cos(ROTATION_ANGLE);
        float sinAngle = sin(ROTATION_ANGLE);
        
        // Rotate view direction around up vector
        float newViewDirx = viewDirx * cosAngle + (normUpy * viewDirz - normUpz * viewDiry) * sinAngle;
        float newViewDiry = viewDiry * cosAngle + (normUpz * viewDirx - normUpx * viewDirz) * sinAngle;
        float newViewDirz = viewDirz * cosAngle + (normUpx * viewDiry - normUpy * viewDirx) * sinAngle;
        
        // Update center position
        centerx = eyex + newViewDirx * viewLength;
        centery = eyey + newViewDiry * viewLength;
        centerz = eyez + newViewDirz * viewLength;
        break;
    }
    case '3': // Rotate camera upwards
    {
        float ROTATION_ANGLE = rotationSpeed;
        //  the new vector vireDir is the vector from the camera to the reference point
        float viewDirx = centerx - eyex;
        float viewDiry = centery - eyey;
        float viewDirz = centerz - eyez;

        // Normalize the view direction
        // getting the magnitude of the viewDir vector
        float viewLength = sqrt(viewDirx * viewDirx + viewDiry * viewDiry + viewDirz * viewDirz);
        // now we are normalizing the viewDir vector as a unit vector
        viewDirx /= viewLength;
        viewDiry /= viewLength;
        viewDirz /= viewLength;

        // Normalize the up vector to get  the direction unit vector
        float upLength = sqrt(upx * upx + upy * upy + upz * upz);
        upx /= upLength;
        upy /= upLength;
        upz /= upLength;

        // we will create a rotation matrix around the up vector
        float cosAngle = cos(ROTATION_ANGLE);
        float sinAngle = sin(ROTATION_ANGLE);

        // Rotate the view direction, this is the new view direction from the camera to the reference point, here there is ROTATION ANGLE difference from now viewDir vector, and 90 - ROTATION ANGLE reference from the up vector, so we just add cos angle with viewDir anmd sin angle with up vector
        float newViewDirx = viewDirx * cosAngle + upx * sinAngle;
        float newViewDiry = viewDiry * cosAngle + upy * sinAngle;
        float newViewDirz = viewDirz * cosAngle + upz * sinAngle;

        // Rotate the up vector, this is the new up vector here there is ROTATION ANGLE difference from now up vector, and 90 + ROTATION ANGLE difference from the viewDir vector, so we just subtract cos angle with up vector and sin angle with viewDir vector
        float newUpx = upx * cosAngle - viewDirx * sinAngle;
        float newUpy = upy * cosAngle - viewDiry * sinAngle;
        float newUpz = upz * cosAngle - viewDirz * sinAngle;

        // Update the center position by adding the new view direction with the camera position to get the new center
        centerx = eyex + newViewDirx * viewLength;
        centery = eyey + newViewDiry * viewLength;
        centerz = eyez + newViewDirz * viewLength;

        // Update the up vector
        upx = newUpx;
        upy = newUpy;
        upz = newUpz;
        break;
    }
    case '4': // Rotate camera downwards
    {
        // the new vector vireDir is the vector from the camera to the reference point
        float ROTATION_ANGLE = -rotationSpeed;
        float viewDirx = centerx - eyex;
        float viewDiry = centery - eyey;
        float viewDirz = centerz - eyez;

        // Normalize the view direction vector as a unit vector
        float viewLength = sqrt(viewDirx * viewDirx + viewDiry * viewDiry + viewDirz * viewDirz);
        viewDirx /= viewLength;
        viewDiry /= viewLength;
        viewDirz /= viewLength;

        // Normalize the up vector as a unit vector
        float upLength = sqrt(upx * upx + upy * upy + upz * upz);
        upx /= upLength;
        upy /= upLength;
        upz /= upLength;

        // Create rotation matrix around the right vector
        float cosAngle = cos(ROTATION_ANGLE);
        float sinAngle = sin(ROTATION_ANGLE);

        // Rotate the view direction vector
        float newViewDirx = viewDirx * cosAngle + upx * sinAngle;
        float newViewDiry = viewDiry * cosAngle + upy * sinAngle;
        float newViewDirz = viewDirz * cosAngle + upz * sinAngle;

        // Rotate the up vector
        float newUpx = upx * cosAngle - viewDirx * sinAngle;
        float newUpy = upy * cosAngle - viewDiry * sinAngle;
        float newUpz = upz * cosAngle - viewDirz * sinAngle;

        // Update the center position
        centerx = eyex + newViewDirx * viewLength;
        centery = eyey + newViewDiry * viewLength;
        centerz = eyez + newViewDirz * viewLength;

        // Update the up vector
        upx = newUpx;
        upy = newUpy;
        upz = newUpz;
        break;
    }
    case '5':
    {
        float TILT_ANGLE = rotationSpeed; // Positive for clockwise, negative for counter-clockwise

        // Calculate the view direction vector from the camera to the reference point
        float viewDirx = centerx - eyex;
        float viewDiry = centery - eyey;
        float viewDirz = centerz - eyez;

        // Normalize the view direction vector as a unit vector
        float viewLength = sqrt(viewDirx * viewDirx + viewDiry * viewDiry + viewDirz * viewDirz);
        viewDirx /= viewLength;
        viewDiry /= viewLength;
        viewDirz /= viewLength;

        // Normalize the up vector as a unit vector
        float upLength = sqrt(upx * upx + upy * upy + upz * upz);
        upx /= upLength;
        upy /= upLength;
        upz /= upLength;

        // Calculate the right vector
        // cross product of view direction and up vector
        float rightx = viewDiry * upz - viewDirz * upy;
        float righty = viewDirz * upx - viewDirx * upz;
        float rightz = viewDirx * upy - viewDiry * upx;

        // Create rotation matrix using Rodrigues' rotation formula
        float cosAngle = cos(TILT_ANGLE);
        float sinAngle = sin(TILT_ANGLE);

        // Rotate the up vector around the view direction axis, as the TILT_ANGLE differencec from the current up vector, and 90 - TILT_ANGLE difference from the right vector, so we just add cos angle with up vector and sin angle with right vector
        float newUpx = upx * cosAngle + rightx * sinAngle;
        float newUpy = upy * cosAngle + righty * sinAngle;
        float newUpz = upz * cosAngle + rightz * sinAngle;

        // Update the up vector
        upx = newUpx;
        upy = newUpy;
        upz = newUpz;

        // We make the up vector again the unit vector for the next rotation
        upLength = sqrt(upx * upx + upy * upy + upz * upz);
        upx /= upLength;
        upy /= upLength;
        upz /= upLength;
        break;
    }
    case '6':
    {
        float TILT_ANGLE = -rotationSpeed; // Positive for clockwise, negative for counter-clockwise

        // Calculate the view direction vector from the camera to the reference point
        float viewDirx = centerx - eyex;
        float viewDiry = centery - eyey;
        float viewDirz = centerz - eyez;

        // Normalize the view direction vector as a unit vector
        float viewLength = sqrt(viewDirx * viewDirx + viewDiry * viewDiry + viewDirz * viewDirz);
        viewDirx /= viewLength;
        viewDiry /= viewLength;
        viewDirz /= viewLength;

        // Normalize the up vector as a unit vector
        float upLength = sqrt(upx * upx + upy * upy + upz * upz);
        upx /= upLength;
        upy /= upLength;
        upz /= upLength;

        // Calculate the right vector
        // cross product of view direction and up vector
        float rightx = viewDiry * upz - viewDirz * upy;
        float righty = viewDirz * upx - viewDirx * upz;
        float rightz = viewDirx * upy - viewDiry * upx;

        // Create rotation matrix using Rodrigues' rotation formula
        float cosAngle = cos(TILT_ANGLE);
        float sinAngle = sin(TILT_ANGLE);

        // Rotate the up vector around the view direction axis , as the TILT_ANGLE differencec from the current up vector, and 90 - TILT_ANGLE difference from the right vector, so we just add cos angle with up vector and sin angle with right vector
        float newUpx = upx * cosAngle + rightx * sinAngle;
        float newUpy = upy * cosAngle + righty * sinAngle;
        float newUpz = upz * cosAngle + rightz * sinAngle;

        // Update the up vector
        upx = newUpx;
        upy = newUpy;
        upz = newUpz;

        // We make the up vector again the unit vector for the next rotation
        upLength = sqrt(upx * upx + upy * upy + upz * upz);
        upx /= upLength;
        upy /= upLength;
        upz /= upLength;
        break;
    }
    case 'a':
        isAxes = !isAxes;
        break; // Toggle axes
    // --- Program Control ---
    case 27:
        exit(0);
        break; // ESC key: exit program
    default:
        break;
    }
    glutPostRedisplay();
}

/**
 * Special key input handler (arrow keys, function keys)
 * Provides camera orbit functionality
 */
void specialKeyListener(int key, int x, int y)
{
    // Calculate view direction vector
    double lx = centerx - eyex;
    double lz = centerz - eyez;
    double s;

    switch (key)
    {
    case GLUT_KEY_LEFT:
    {
        // Calculate the vector from the camera to the reference point
        double cx = centerx - eyex;
        double cy = centery - eyey;
        double cz = centerz - eyez;

        // Calculate the cross product of the head vector and the camera-reference point vector
        double leftx = upy * cz - upz * cy;
        double lefty = upz * cx - upx * cz;
        double leftz = upx * cy - upy * cx;
        // Now the cross product of the head vector and the camera-reference point vector is the left vector
        // as cross vector is at the perpendicular of two vectors

        // Normalize the left vector as we only need the unit vector (direction)
        double length = sqrt(leftx * leftx + lefty * lefty + leftz * leftz); // calculating the length
        leftx /= length;
        lefty /= length;
        leftz /= length;

        // Move the camera position in the left direction by the movement speed
        eyex += leftx * movementSpeed;
        eyey += lefty * movementSpeed;
        eyez += leftz * movementSpeed;

        // Move the look-at point in the left direction by the movement speed
        centerx += leftx * movementSpeed;
        centery += lefty * movementSpeed;
        centerz += leftz * movementSpeed;

        break;
    }

    case GLUT_KEY_RIGHT:
    {
        // Calculate the vector from the camera to the reference point
        double cx = centerx - eyex;
        double cy = centery - eyey;
        double cz = centerz - eyez;

        // Calculate the cross product of the head vector and the camera-reference point vector
        double rightx = upy * cz - upz * cy;
        double righty = upz * cx - upx * cz;
        double rightz = upx * cy - upy * cx;
        // Now the cross product of the head vector and the camera-reference point vector is the right vector
        // as cross vector is at the perpendicular of two vectors
        // Normalize the right vector as we only need the unit vector (direction)
        double length = sqrt(rightx * rightx + righty * righty + rightz * rightz); // calculating the length
        rightx /= length;
        righty /= length;
        rightz /= length;

        // Move the camera position in the right direction by the movement speed
        eyex -= rightx * movementSpeed;
        eyey -= righty * movementSpeed;
        eyez -= rightz * movementSpeed;

        // Move the look-at point in the right direction by the movement speed
        centerx -= rightx * movementSpeed;
        centery -= righty * movementSpeed;
        centerz -= rightz * movementSpeed;

        break;
    }

    case GLUT_KEY_UP:
    {
        // Calculate the vector from the camera to the reference point
        double cx = centerx - eyex;
        double cy = centery - eyey;
        double cz = centerz - eyez;

        // Normalize the vector as we only need the unit vector (direction)
        double length = sqrt(cx * cx + cy * cy + cz * cz);
        cx /= length;
        cy /= length;
        cz /= length;

        // Move the camera position towards the c vector by a single movementSpeed unit
        eyex += cx * movementSpeed;
        eyey += cy * movementSpeed;
        eyez += cz * movementSpeed;

        // Move the look-at point towards the c vector by a single movementSpeed unit
        centerx += cx * movementSpeed;
        centery += cy * movementSpeed;
        centerz += cz * movementSpeed;

        break;
    }

    case GLUT_KEY_DOWN:
    {
        // Calculate the vector from the camera to the reference point
        double cx = centerx - eyex;
        double cy = centery - eyey;
        double cz = centerz - eyez;

        // Normalize the vector as we only need the unit vector (direction)
        double length = sqrt(cx * cx + cy * cy + cz * cz);
        cx /= length;
        cy /= length;
        cz /= length;

        // Move the camera position away from the c vector by a single movementSpeed unit
        eyex -= cx * movementSpeed;
        eyey -= cy * movementSpeed;
        eyez -= cz * movementSpeed;

        // Move the look-at point away from the c vector by a single movementSpeed unit
        centerx -= cx * movementSpeed;
        centery -= cy * movementSpeed;
        centerz -= cz * movementSpeed;

        break;
    }
    case GLUT_KEY_PAGE_UP:
    {
        // here we know the up side is the side of the head vector, which in this case is the up vector
        // Move the camera position along the head vector by a single unit
        eyex += upx * movementSpeed;
        eyey += upy * movementSpeed;
        eyez += upz * movementSpeed;

        // Move the reference point towards the up vector by a single unit of movement speed
        centerx += upx * movementSpeed;
        centery += upy * movementSpeed;
        centerz += upz * movementSpeed;

        break;
    }

    case GLUT_KEY_PAGE_DOWN:
    {
        // here we know the down side is the other side of the head vector, which in this case is the up vector
        // Move the camera position along the opposite direction of the head vector by a single unit
        eyex -= upx * movementSpeed;
        eyey -= upy * movementSpeed;
        eyez -= upz * movementSpeed;

        // Move the reference point towards the opposite direction of the up vector by a single unit of movement speed
        centerx -= upx * movementSpeed;
        centery -= upy * movementSpeed;
        centerz -= upz * movementSpeed;

        break;
    }
    }

    glutPostRedisplay(); // Request a screen refresh
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(windowWidth, windowHeight);
    glutCreateWindow("3D Scene with Ray Tracing");

    initGL();
    
    // Print keyboard instructions
    std::cout << "\n=== RAY TRACER CONTROLS ===" << std::endl;
    std::cout << "Camera Movement:" << std::endl;
    std::cout << "  Arrow keys: Move forward/back/left/right" << std::endl;
    std::cout << "  Page Up/Down: Move up/down" << std::endl;
    std::cout << "  1/2: Rotate look-at point left/right" << std::endl;
    std::cout << "  3/4: Rotate camera up/down" << std::endl;
    std::cout << "  5/6: Tilt camera" << std::endl;
    std::cout << "  w/s: Move camera up/down" << std::endl;
    std::cout << "Rendering:" << std::endl;
    std::cout << "  0: Capture ray-traced image" << std::endl;
    std::cout << "  t/T: Toggle floor texture (checkerboard/texture)" << std::endl;
    std::cout << "  a: Toggle coordinate axes" << std::endl;
    std::cout << "  ESC: Exit" << std::endl;
    std::cout << "========================\n" << std::endl;

    // Load the scene from the file
    loadData();
    
    // Try to load texture (optional)
    if (TextureManager::loadTexture(textureFileName)) {
        std::cout << "Texture loaded successfully. Press 't' to toggle between checkerboard and texture." << std::endl;
    } else {
        std::cout << "Texture not found or failed to load. Using checkerboard pattern only." << std::endl;
    }

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(specialKeyListener); // Add special key handler for arrow keys

    glutMainLoop();

    // Clean up
    TextureManager::cleanup();
    for (Object *obj : objects)
        delete obj;
    for (PointLight *light : pointLights)
        delete light;
    for (SpotLight *light : spotLights)
        delete light;

    return 0;
}