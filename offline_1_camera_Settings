/**
 * OpenGL 3D Drawing Demo
 *
 * This program demonstrates basic 3D rendering with OpenGL and GLUT including:
 * - Camera positioning with gluLookAt
 * - Drawing 3D shapes (cube and pyramid)
 * - Keyboard navigation for camera control
 * - Perspective projection
 * - Object toggling
 */

// --- Includes ---
// Standard Headers
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <bits/stdc++.h>
#include <GL/glew.h>
using namespace std;

// OpenGL / GLUT Headers
#ifdef __APPLE__
#include <GLUT/glut.h> // Use GLUT framework on macOS
#else
#include <GL/glut.h> // Use standard GLUT location on Linux/Windows
#endif

/// @brief the force which is applied to the sphere towards land
float gravity = -9.8f;
/// @brief the coefficient of friction when it rolls on the ground
float friction = 0.98f;
/// @brief  the percentage of force that helps the ball to bounce back from the floor after a collision
float restitution = 0.8f;
/// @brief defines the size of the cube with checkered floor
float cubeSize = 20.0f;
/// @brief the value of pi used in calculations
float pi = 3.14159f;
/// @brief increase of ball velocity per plus key press
float increasePerPlus = 10.0f;
/// @brief How quickly tthe sphere's angular velocity adjusts towards its ideal value when it hits a surface
const float rollMatchFactor = 0.5f;
/// @brief Limits how fast the sphere is allowed to rotate
const float maxAngularSpeed = 50.0f;
/// @brief  quadric to draw the sphere
GLUquadricObj *quadric;

/// @brief Now below are the original values of the sphere

float originalSphereRadius = 0.5f;
float originalSphereMass = 1.0f;
float originalSphereVelocity[] = {1.0f, 2.0f, 2.0f};
float originalSpherePosition[] = {5.0f, 5.0f, 5.0f};
float originalSphereAngularVelocity[] = {0.0f, 0.0f, 0.0f};
float originalSphereRotationAngle[] = {0.0f, 0.0f, 0.0f};
float originalSphereColor[] = {0.8f, 0.2f, 0.2f};

// --- Global Variables ---
/// @brief animation speed or how much will be the period between two timer function calls
int animationSpeed = 10;
/// @brief Camera position and orientation
GLfloat eyex = 4, eyey = 4, eyez = 4;          // Camera position coordinates
GLfloat centerx = 0, centery = 0, centerz = 0; // Look-at point coordinates
GLfloat upx = 0, upy = 1, upz = 0;             // Up vector coordinates

float movementSpeed = 0.3f; // Speed of camera movement
float rotationSpeed = 0.1f; // Speed of camera rotation

// Object visibility flags
bool isAxes = true;     // Toggle for coordinate axes
bool isCube = false;    // Toggle for cube
bool isPyramid = false; // Toggle for pyramid
bool paused = true;     // Toggle for pausing the animation
bool showArrow = true;  // Toggle for showing the arrow

// --- Function Declarations ---
void initGL();
void display();
void reshapeListener(GLsizei width, GLsizei height);
void keyboardListener(unsigned char key, int x, int y);
void specialKeyListener(int key, int x, int y);
void drawAxes();
void drawCube();
void drawPyramid();
void drawCubeWithCheckeredFloor();

/**
 * Initialize OpenGL settings
 * Sets up background color and enables depth testing
 */
void initGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black background
    glEnable(GL_DEPTH_TEST);              // Enable depth testing for z-culling
}

/// @brief struct which represents a sphere
typedef struct
{
    float radius;             // Radius of the sphere
    float mass;               // Mass of the sphere
    float position[3];        // Position vector of the sphere
    float velocity[3];        // Velocity vector of the sphere
    float angularVelocity[3]; // Angular velocity vector of the sphere
    float rotationAngle[3];   // Rotation angle vector of the sphere

} Sphere;

/// @brief sphere object we are going to use
Sphere sphere;

/// @brief set all the values of the sphere to their original values defined above
void initSphere()
{
    sphere.radius = originalSphereRadius;
    sphere.mass = originalSphereMass;
    sphere.position[0] = originalSpherePosition[0];
    sphere.position[1] = originalSpherePosition[1];
    sphere.position[2] = originalSpherePosition[2];
    sphere.velocity[0] = originalSphereVelocity[0];
    sphere.velocity[1] = originalSphereVelocity[1];
    sphere.velocity[2] = originalSphereVelocity[2];
    sphere.angularVelocity[0] = originalSphereAngularVelocity[0];
    sphere.angularVelocity[1] = originalSphereAngularVelocity[1];
    sphere.angularVelocity[2] = originalSphereAngularVelocity[2];
    sphere.rotationAngle[0] = originalSphereRotationAngle[0];
    sphere.rotationAngle[1] = originalSphereRotationAngle[1];
    sphere.rotationAngle[2] = originalSphereRotationAngle[2];

    quadric = gluNewQuadric();
    gluQuadricNormals(quadric, GLU_SMOOTH);
}

/// @brief  to add the stripes to the sphere
/// @param color the two colors we want add to the sphere
/// @param angle the angle of the sphere colors
void stripeColor(GLfloat *color, float angle)
{
    // makes 18 stripes at the longitude of the sphere
    int stripe = (int)(angle / (2 * pi / 18)) % 2;
    if (stripe == 0)
    {
        color[0] = 1.0f; // Red
        color[1] = 0.0f;
        color[2] = 0.0f;
    }
    else
    { // Green
        color[0] = 0.0f;
        color[1] = 1.0f;
        color[2] = 0.0f;
    }
}

/// @brief draw the sphere with the given color stripes
void drawSphere()
{
    glPushMatrix(); // save the current GL state

    glTranslatef(sphere.position[0], sphere.position[1], sphere.position[2]); // position the sphere using the position vector
    glRotatef(sphere.rotationAngle[0], 1.0f, 0.0f, 0.0f);                     // apply rotation to the x axis
    glRotatef(sphere.rotationAngle[1], 0.0f, 1.0f, 0.0f);                     // apply rotation to the y axis
    glRotatef(sphere.rotationAngle[2], 0.0f, 0.0f, 1.0f);                     // apply rotation to the z axis

    glEnable(GL_COLOR_MATERIAL); // enable color material to track GLcolor
    gluQuadricCallback(quadric, GLU_ERROR, NULL);
    gluQuadricDrawStyle(quadric, GLU_FILL); // fills surface style
    gluQuadricNormals(quadric, GLU_SMOOTH); // configure a quadric object to draw smooth surfaces
    gluQuadricTexture(quadric, GL_TRUE);

    // Set up the color callback
    GLfloat color[3];
    glBegin(GL_TRIANGLE_STRIP);
    for (float phi = 0; phi < pi; phi += pi / 20)
    {
        // it goes vertically from the top to the bottom, in 20 steps
        for (float theta = 0; theta <= 2 * pi; theta += pi / 20)
        {
            // it goes around the sphere in 40 steps horizontally
            stripeColor(color, theta);
            glColor3fv(color);
            glVertex3f(
                sphere.radius * sin(phi) * cos(theta),
                sphere.radius * cos(phi),
                sphere.radius * sin(phi) * sin(theta));

            stripeColor(color, theta);
            glColor3fv(color);
            glVertex3f(
                sphere.radius * sin(phi + pi / 20) * cos(theta),
                sphere.radius * cos(phi + pi / 20),
                sphere.radius * sin(phi + pi / 20) * sin(theta));
        }
    }
    glEnd();

    glDisable(GL_COLOR_MATERIAL);
    glPopMatrix();
}

/// @brief this function draws the velocity arrow of the sphere using the velocity vector of the sphere
void drawVelocityArrow()
{
    float speed = sqrt(sphere.velocity[0] * sphere.velocity[0] +
                       sphere.velocity[1] * sphere.velocity[1] +
                       sphere.velocity[2] * sphere.velocity[2]); // calculating the speed magnitude

    if (speed == 0.0f)
        return; /// we wont show the vector if there is no speed

    glPushMatrix();
    glTranslatef(sphere.position[0], sphere.position[1], sphere.position[2]); // we are finding the position vector of the sphere to start the speed arrow

    // Calculate arrow direction (normalized velocity)
    float dirX = sphere.velocity[0] / speed;
    float dirY = sphere.velocity[1] / speed;
    float dirZ = sphere.velocity[2] / speed;

    // Scale arrow length proportional to velocity magnitude
    float arrowLength = speed * 0.5f; // Adjust scale factor as needed, and multiply with speed to show the arrow length proportional to speed
    float headLength = 0.2f;          // we want the head length to be a bit longer, but constant size
    float headWidth = 0.1f;           // we want the head width to be a bit smaller, but constant size

    // Set arrow color (red)
    glColor3f(1.0f, 0.0f, 0.0f);
    glLineWidth(2.0f); // making the line width bigger to make it more visible

    // Draw arrow shaft
    glBegin(GL_LINES);
    glVertex3f(0, 0, 0);                                                    // start at the ccenter of the sphere
    glVertex3f(dirX * arrowLength, dirY * arrowLength, dirZ * arrowLength); // end the arrow line at the direction of the velocity
    glEnd();

    // Draw arrow head (cone)
    glPushMatrix();
    glTranslatef(dirX * arrowLength, dirY * arrowLength, dirZ * arrowLength); // Now we are at the end of the arrow

    // Rotate to point in velocity direction
    float angle = acos(dirY) * 180.0f / pi;
    float axisX = -dirZ;
    float axisZ = dirX;

    if (fabs(angle) > 0.001f)
    {
        glRotatef(angle, axisX, 0.0f, axisZ);
    }

    glutSolidCone(headWidth, headLength, 8, 1); // make the arrow head a solid cone
    glPopMatrix();

    glPopMatrix();
}

/// @brief this function checks for collisions for the sphere in all directions. that could be on the roof, on the 4 walls or bouncing off the floor. This method also check for rotation and angular rotational increase when the sphere collides with a surface like in the real world.
void checkCollisions()
{

    // Floor collision (Y-axis)
    if (sphere.position[1] - sphere.radius <= 0.0f)
    {
        sphere.position[1] = sphere.radius;                     // Prevent sinking through the floor by positioning the center of the sphere minimum redius above the floor
        sphere.velocity[1] = -sphere.velocity[1] * restitution; // as the sphere hits the floor, its velocity changes in the opposite direction with fraction of restitution

        // Apply friction to horizontal velocities
        sphere.velocity[0] *= friction; // as the sphere hits the floor, its horizontal x velocity remains at the same direction but with a fraction of friction
        sphere.velocity[2] *= friction; // as the sphere hits the floor, its horizontal z velocity remains at the same direction but with a fraction of friction

        // sphere.angularVelocity[0] += sphere.velocity[2] / sphere.radius;
        // sphere.angularVelocity[2] += -sphere.velocity[0] / sphere.radius;
    }

    // Wall collisions (X-axis)
    if (sphere.position[0] - sphere.radius <= 0.0f)
    {
        sphere.position[0] = sphere.radius; // set the sphere position to radius away from the wall
        sphere.velocity[0] = -sphere.velocity[0] * restitution; // the velocity will be in the opposite direction

        // Smoothly adjust Y-axis rotation

        // sphere.angularVelocity[1] += -sphere.velocity[2] / sphere.radius;
    }
    if (sphere.position[0] + sphere.radius >= cubeSize)
    {
        sphere.position[0] = cubeSize - sphere.radius; // set the sphere position to radius away from the wall
        sphere.velocity[0] = -sphere.velocity[0] * restitution; // the velocity will be in the opposite direction

        // sphere.angularVelocity[1] += -sphere.velocity[2] / sphere.radius;
    }

    // Wall collisions (Z-axis)
    if (sphere.position[2] - sphere.radius <= 0.0f)
    {
        sphere.position[2] = sphere.radius;
        sphere.velocity[2] = -sphere.velocity[2] * restitution;

        // sphere.angularVelocity[1] += sphere.velocity[0] / sphere.radius;
    }
    if (sphere.position[2] + sphere.radius >= cubeSize)
    {
        sphere.position[2] = cubeSize - sphere.radius;
        sphere.velocity[2] = -sphere.velocity[2] * restitution;

        // sphere.angularVelocity[1] += sphere.velocity[0] / sphere.radius;
    }

    // Ceiling collision (Y-axis top)
    if (sphere.position[1] + sphere.radius >= cubeSize)
    {
        sphere.position[1] = cubeSize - sphere.radius;
        sphere.velocity[1] = -sphere.velocity[1] * restitution;
    }
}

/// @brief updatePhysics of the ball
/// @param deltaTime it is the millisecond time after when the function is called
void updatePhysics(int deltaTime)
{
    float dt = deltaTime / 1000.0f; // Convert to seconds

    // Apply gravity
    sphere.velocity[1] += gravity * dt; // in every moment, the sphere will be affected by gravity

    // Update position
    sphere.position[0] += sphere.velocity[0] * dt; // now the sphere is moving according to its velocity
    sphere.position[1] += sphere.velocity[1] * dt;
    sphere.position[2] += sphere.velocity[2] * dt;

    checkCollisions(); // check for collisions

    // Calculate angular velocity according to velocity and radius
    sphere.angularVelocity[0] = sphere.velocity[2] / sphere.radius; 
    sphere.angularVelocity[1] = 0;                     
    sphere.angularVelocity[2] = -sphere.velocity[0] / sphere.radius;

    // Dampen angular velocity
    // sphere.angularVelocity[0] *= friction;
    // sphere.angularVelocity[1] *= friction;
    // sphere.angularVelocity[2] *= friction;

    // Update rotation according to angular velocity
    sphere.rotationAngle[0] += sphere.angularVelocity[0] * dt * 180.0f / pi;
    sphere.rotationAngle[1] += sphere.angularVelocity[1] * dt * 180.0f / pi;
    sphere.rotationAngle[2] += sphere.angularVelocity[2] * dt * 180.0f / pi;

    // Just to find out the current position for the sphere
    // cout << "Sphere position: (" << sphere.position[0] << ", " << sphere.position[1] << ", " << sphere.position[2] << ")" << endl;
}
/**
 * Main display function
 * Sets up the camera and renders visible objects
 */
void display()
{
    // Clear color and depth buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Set up the model-view matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Position camera using the eye, center and up vectors
    gluLookAt(eyex, eyey, eyez,          // Camera position
              centerx, centery, centerz, // Look-at point
              upx, upy, upz);            // Up vector

    // Draw objects based on visibility flags
    drawCubeWithCheckeredFloor();
    drawSphere();
    if (showArrow)
    {
        drawVelocityArrow();
    }
    if (isAxes)
        drawAxes();

    // Swap buffers (double buffering)
    glutSwapBuffers();
}

/**
 * Window reshape callback
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
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}

/**
 * Keyboard input handler for standard keys
 * Manages camera position, object visibility, and program exit
 */
void keyboardListener(unsigned char key, int x, int y)
{

    // Calculate view direction vector
    double lx = centerx - eyex;
    double lz = centerz - eyez;
    double s;

    switch (key)
    {
    case '1':
    {
        // we are getting the vector which originates from the camera position to the reference point
        double caX = centerx - eyex;
        double caY = centery - eyey;
        double caZ = centerz - eyez;

        // we are calculating the cross product of the new vector and the head vector
        double crossX = upy * caZ - upz * caY;
        double crossY = upz * caX - upx * caZ;
        double crossZ = upx * caY - upy * caX;
        // as a result the new cross vector will be the perpendicular of the head and the camera to refernce point vector
        // so it might be on the left or on the right side of the current camera position

        // getting the magnitude value of the cross product vector
        double crossMag = sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);
        // Now, we are making the the new t vector as the unit vector for the direction of cross vector
        double tX = crossX / crossMag;
        double tY = crossY / crossMag;
        double tZ = crossZ / crossMag;
        // Now, we are incrementing the camera position in the direction of the t vector
        double daX = tX * rotationSpeed;
        double daY = tY * rotationSpeed;
        double daZ = tZ * rotationSpeed;

        centerx = centerx + daX;
        centery = centery + daY;
        centerz = centerz + daZ;
        break;
    }
    case '2':
    {
        // we are getting the vector which originates from the camera position to the reference point
        double caX = centerx - eyex;
        double caY = centery - eyey;
        double caZ = centerz - eyez;

        // we are calculating the cross product of the new vector and the head vector
        double crossX = upy * caZ - upz * caY;
        double crossY = upz * caX - upx * caZ;
        double crossZ = upx * caY - upy * caX;

        // getting the magnitude value of the cross product vector
        double crossMag = sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);
        double tX = crossX / crossMag;
        double tY = crossY / crossMag;
        double tZ = crossZ / crossMag;

        // Now, we are calculating the value to change 
        double daX = tX * rotationSpeed;
        double daY = tY * rotationSpeed;
        double daZ = tZ * rotationSpeed;
        // Now, we are decrementing the camera position in the direction of the t vector
        centerx = centerx - daX;
        centery = centery - daY;
        centerz = centerz - daZ;
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
        // here the rotation angle is the angle between the right vector and the view direction vector, so we just add cos angle with viewDir anmd sin angle with right vector
        // so as it is minus, it will rotate at the different direction
        // Rotate the view direction vector this is the new view direction from the camera to the reference point here there is ROTATION ANGLE difference from now viewDir vector, and 90 - ROTATION ANGLE reference from the up vector, so we just add cos angle with viewDir anmd sin angle with up vector
        float newViewDirx = viewDirx * cosAngle + upx * sinAngle;
        float newViewDiry = viewDiry * cosAngle + upy * sinAngle;
        float newViewDirz = viewDirz * cosAngle + upz * sinAngle;

        // Rotate the up vector this is the new up vector, here there is ROTATION ANGLE difference from now up vector, and 90 + ROTATION ANGLE difference from the viewDir vector, so we just subtract cos angle with up vector and sin angle with viewDir vector
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

   

    case 'w':
    {
        // Move camera up without changing the reference point
        eyey += movementSpeed;
        break;
    }
    case 's':
        eyey -= movementSpeed;
        break; // Move camera down withiout changing the refernce point

    case 't':
        centerz += movementSpeed;
        break; // Move look-at point forward
    case 'y':
        centerz -= movementSpeed;
        break; // Move look-at point backward
    case '*':
        movementSpeed += 0.1f;
        rotationSpeed += 0.1f;
        break;
    case '/':
        movementSpeed -= 0.1f;
        rotationSpeed -= 0.1f;
        break;
    case '+':
        sphere.velocity[0] += increasePerPlus;
        sphere.velocity[1] += increasePerPlus;
        sphere.velocity[2] += increasePerPlus;
        break;
    case '-':
        sphere.velocity[0] -= increasePerPlus;
        sphere.velocity[1] -= increasePerPlus;
        sphere.velocity[2] -= increasePerPlus;
        break;
    case ' ':
        paused = !paused;
        break;
    case 'r':
        if (paused)
        {
            initSphere(); // only reset key is appliable if paused currently
        }
        break;

    // --- Object Visibility Toggles ---
    case 'a':
        isAxes = !isAxes;
        break; // Toggle axes
    case 'c':
        isCube = !isCube;
        break; // Toggle cube
    case 'p':
        isPyramid = !isPyramid;
        break; // Toggle pyramid
    case 'v':
        showArrow = !showArrow; // toggle the arrow of the sphere velocity
        break;

    // --- Program Control ---
    case 27:
        exit(0);
        break; // ESC key: exit program
    }

    glutPostRedisplay(); // Request a screen refresh
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

// Simulate physics

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


/// @brief draws a checkered floor made of alternating black and white square tiles
/// @param size what will be the size of the floor
/// @param tiles how many tiles in each direction
void drawCheckeredFloor(float size, int tiles)
{
    // Set up the colors
    GLfloat white[] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat black[] = {0.0f, 0.0f, 0.0f, 1.0f};

    float tileSize = size / tiles;

    // Draw the floor on XZ plane starting from (0,0,0)
    glBegin(GL_QUADS);
    for (int x = 0; x < tiles; x++)
    {
        for (int z = 0; z < tiles; z++)
        {
            if ((x + z) % 2 == 0)
            {
                glColor4fv(white); // checking the color of the tile
            }
            else
            {
                glColor4fv(black);
            }

            // draws a tile of the floor which will be tilesize length long, and xz plane
            glVertex3f(x * tileSize, 0.0f, z * tileSize);
            glVertex3f((x + 1) * tileSize, 0.0f, z * tileSize);
            glVertex3f((x + 1) * tileSize, 0.0f, (z + 1) * tileSize);
            glVertex3f(x * tileSize, 0.0f, (z + 1) * tileSize);
        }
    }
    glEnd();
}

void drawCubeWithCheckeredFloor()
{
    // Set up the cube's dimensions

    int tiles = cubeSize;
    // Number of tiles in each direction

    // Set up the colors
    GLfloat wallColor[] = {0.5f, 1.0f, 0.5f, 1.0f}; // the color of the walls
    GLfloat ceilingColor[] = {0.5f, 0.5f, 0.5f, 1.0f}; // the color of the ceiling

    // Draw the checkered floor starting from origin
    drawCheckeredFloor(cubeSize, tiles);

    // Draw the walls
    glBegin(GL_QUADS);
    glColor4fv(wallColor);
    // Front wall with the wall color
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(cubeSize, 0.0f, 0.0f);
    glVertex3f(cubeSize, cubeSize, 0.0f);
    glVertex3f(0.0f, cubeSize, 0.0f);

    // drawing the back wall with the wall color
    glVertex3f(0.0f, 0.0f, cubeSize);
    glVertex3f(cubeSize, 0.0f, cubeSize);
    glVertex3f(cubeSize, cubeSize, cubeSize);
    glVertex3f(0.0f, cubeSize, cubeSize);

    // drawing the left wall with the wall color
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, cubeSize);
    glVertex3f(0.0f, cubeSize, cubeSize);
    glVertex3f(0.0f, cubeSize, 0.0f);

    // drawing the right wall with the wall color
    glVertex3f(cubeSize, 0.0f, 0.0f);
    glVertex3f(cubeSize, 0.0f, cubeSize);
    glVertex3f(cubeSize, cubeSize, cubeSize);
    glVertex3f(cubeSize, cubeSize, 0.0f);
    glEnd();

    // Drawing the ceiling with the color of the ceiling
    glBegin(GL_QUADS);
    glColor4fv(ceilingColor); // cahnge the color as the color of the ceiling is different from the wall color
    glVertex3f(0.0f, cubeSize, 0.0f);
    glVertex3f(cubeSize, cubeSize, 0.0f);
    glVertex3f(cubeSize, cubeSize, cubeSize);
    glVertex3f(0.0f, cubeSize, cubeSize);
    glEnd();
}

/**
 * Timer function for animation.
 * This demonstrates the use of a timer instead of idle function.
 * Timer functions provide better control over animation speed.
 *
 * @param value Value passed to the timer function (not used here)
 */
void timerFunction(int value)
{
    // Update animation values
    // cout << "Timer function called" << endl;
    time_t currentTime = time(NULL);
    struct tm *timeInfo = localtime(&currentTime);
    if (paused == false)
    {
        // we wont apply physics if the scene is paused
        updatePhysics((float)(animationSpeed));
    }
    // Request a redisplay
    glutPostRedisplay();

    // Register the timer again
    glutTimerFunc(animationSpeed, timerFunction, 0);
}

/**
 * Main function: Program entry point
 */
int main(int argc, char **argv)
{
    // Initialize GLUT
    glutInit(&argc, argv);

    // Configure display mode and window
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(640, 640);
    glutInitWindowPosition(50, 50);
    glutCreateWindow("OpenGL 3D Drawing");

    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutTimerFunc(animationSpeed, timerFunction, 0);
    initSphere(); // initialize the sphere object
    // Initialize OpenGL settings
    initGL();

    // Enter the GLUT event loop
    glutMainLoop();

    return 0;
}