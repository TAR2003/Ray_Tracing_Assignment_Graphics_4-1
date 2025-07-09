/**
 * OpenGL 3D Camera Movement Demo
 *
 * This program demonstrates basic camera movement with OpenGL and GLUT including:
 * - Camera positioning with gluLookAt
 * - Keyboard navigation for camera control
 * - Perspective projection
 * - Coordinate axes display
 */

// --- Includes ---
// Standard Headers
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

// OpenGL / GLUT Headers
#ifdef __APPLE__
#include <GLUT/glut.h> // Use GLUT framework on macOS
#else
#include <GL/glut.h> // Use standard GLUT location on Linux/Windows
#endif

// --- Global Variables ---
/// @brief Camera position and orientation
GLfloat eyex = 4, eyey = 4, eyez = 4;          // Camera position coordinates
GLfloat centerx = 0, centery = 0, centerz = 0; // Look-at point coordinates
GLfloat upx = 0, upy = 1, upz = 0;             // Up vector coordinates

float movementSpeed = 0.3f; // Speed of camera movement
float rotationSpeed = 0.1f; // Speed of camera rotation

// Object visibility flags
bool isAxes = true;     // Toggle for coordinate axes

// --- Function Declarations ---
void initGL();
void display();
void reshapeListener(GLsizei width, GLsizei height);
void keyboardListener(unsigned char key, int x, int y);
void specialKeyListener(int key, int x, int y);
void drawAxes();

/**
 * Initialize OpenGL settings
 * Sets up background color and enables depth testing
 */
void initGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black background
    glEnable(GL_DEPTH_TEST);              // Enable depth testing for z-culling
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

    // --- Object Visibility Toggles ---
    case 'a':
        isAxes = !isAxes;
        break; // Toggle axes

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
    glutCreateWindow("OpenGL 3D Camera Movement");

    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);

    glutDisplayFunc(display);
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    
    // Initialize OpenGL settings
    initGL();

    // Enter the GLUT event loop
    glutMainLoop();

    return 0;
}