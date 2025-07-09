#include "main.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#define M_PI 3.14159265358979323846

// Camera variables
GLdouble eyeX = 4, eyeY = 4, eyeZ = 4;
GLdouble lookX = 0, lookY = 0, lookZ = 0;
GLdouble upX = 0, upY = 1, upZ = 0;
float movementSpeed = 0.3f;
float rotationSpeed = 0.1f;
bool isAxes = true;

// Ray tracing global variables
std::vector<Object*> objects;
std::vector<PointLight> pointLights;
std::vector<SpotLight> spotLights;
int recursionLevel = 1;
int imageWidth = 768, imageHeight = 768;
int captureCount = 1;

// Floor object
Floor* floor_obj = nullptr;

// === SPHERE IMPLEMENTATION ===
void Sphere::draw() {
    glPushMatrix();
    glTranslatef(reference_point.x, reference_point.y, reference_point.z);
    glColor3f(color.r, color.g, color.b);
    glutSolidSphere(length, 50, 50);
    glPopMatrix();
}

double Sphere::intersect(Ray* r, Color* color, int level) {
    Vector3D oc = r->start - reference_point;
    double a = r->dir.dot(r->dir);
    double b = 2.0 * oc.dot(r->dir);
    double c = oc.dot(oc) - length * length;
    double discriminant = b * b - 4 * a * c;
    
    if (discriminant < 0) return -1;
    
    double t = (-b - sqrt(discriminant)) / (2.0 * a);
    if (t < 0) {
        t = (-b + sqrt(discriminant)) / (2.0 * a);
        if (t < 0) return -1;
    }
    
    if (level == 0) return t;
    
    Vector3D intersectionPoint = r->start + r->dir * t;
    Vector3D normal = getNormal(intersectionPoint);
    *color = calculateLighting(r, intersectionPoint, normal, this->color, this, level);
    return t;
}

Vector3D Sphere::getNormal(const Vector3D& point) const {
    return (point - reference_point).normalize();
}

Color Sphere::getColorAt(const Vector3D& point) const {
    return color;
}

// === TRIANGLE IMPLEMENTATION ===
void Triangle::draw() {
    glColor3f(color.r, color.g, color.b);
    glBegin(GL_TRIANGLES);
    glVertex3f(a.x, a.y, a.z);
    glVertex3f(b.x, b.y, b.z);
    glVertex3f(c.x, c.y, c.z);
    glEnd();
}

double Triangle::intersect(Ray* r, Color* color, int level) {
    Vector3D edge1 = b - a;
    Vector3D edge2 = c - a;
    Vector3D h = r->dir.cross(edge2);
    double det = edge1.dot(h);
    
    if (fabs(det) < 1e-6) return -1;
    
    double invDet = 1.0 / det;
    Vector3D s = r->start - a;
    double u = invDet * s.dot(h);
    
    if (u < 0.0 || u > 1.0) return -1;
    
    Vector3D q = s.cross(edge1);
    double v = invDet * r->dir.dot(q);
    
    if (v < 0.0 || u + v > 1.0) return -1;
    
    double t = invDet * edge2.dot(q);
    if (t < 0) return -1;
    
    if (level == 0) return t;
    
    Vector3D intersectionPoint = r->start + r->dir * t;
    Vector3D normal = getNormal(intersectionPoint);
    *color = calculateLighting(r, intersectionPoint, normal, this->color, this, level);
    return t;
}

Vector3D Triangle::getNormal(const Vector3D& point) const {
    return (b - a).cross(c - a).normalize();
}

Color Triangle::getColorAt(const Vector3D& point) const {
    return color;
}

// === FLOOR IMPLEMENTATION ===
Floor::Floor(double floorWidth, double tileWidth) : floorWidth(floorWidth), tileWidth(tileWidth) {
    reference_point = Vector3D(-floorWidth/2, -floorWidth/2, 0);
    length = tileWidth;
}

void Floor::draw() {
    bool colorToggle = false;
    for (double x = -floorWidth/2; x < floorWidth/2; x += tileWidth) {
        for (double y = -floorWidth/2; y < floorWidth/2; y += tileWidth) {
            if (colorToggle) 
                glColor3f(1, 1, 1);
            else 
                glColor3f(0, 0, 0);
                
            glBegin(GL_QUADS);
            glVertex3f(x, y, 0);
            glVertex3f(x + tileWidth, y, 0);
            glVertex3f(x + tileWidth, y + tileWidth, 0);
            glVertex3f(x, y + tileWidth, 0);
            glEnd();
            
            colorToggle = !colorToggle;
        }
        if ((int)(floorWidth/tileWidth) % 2 == 0) 
            colorToggle = !colorToggle;
    }
}

double Floor::intersect(Ray* r, Color* color, int level) {
    if (fabs(r->dir.z) < 1e-6) return -1;
    
    double t = -(r->start.z) / r->dir.z;
    if (t < 0) return -1;
    
    Vector3D intersectionPoint = r->start + r->dir * t;
    if (fabs(intersectionPoint.x) > floorWidth/2 || fabs(intersectionPoint.y) > floorWidth/2) 
        return -1;
    
    if (level == 0) return t;
    
    Vector3D normal = getNormal(intersectionPoint);
    Color floorColor = getColorAt(intersectionPoint);
    *color = calculateLighting(r, intersectionPoint, normal, floorColor, this, level);
    return t;
}

Vector3D Floor::getNormal(const Vector3D& point) const {
    return Vector3D(0, 0, 1);
}

Color Floor::getColorAt(const Vector3D& point) const {
    int x = (int)((point.x - reference_point.x) / tileWidth);
    int y = (int)((point.y - reference_point.y) / tileWidth);
    
    if ((x + y) % 2 == 0)
        return Color(1, 1, 1); // White
    else
        return Color(0, 0, 0); // Black
}

// === GENERAL QUADRIC IMPLEMENTATION ===
GeneralQuadric::GeneralQuadric(double A, double B, double C, double D, double E,
                               double F, double G, double H, double I, double J) {
    coefficients[0] = A; coefficients[1] = B; coefficients[2] = C;
    coefficients[3] = D; coefficients[4] = E; coefficients[5] = F;
    coefficients[6] = G; coefficients[7] = H; coefficients[8] = I; coefficients[9] = J;
}

void GeneralQuadric::draw() {
    // Not required for quadric surfaces
}

double GeneralQuadric::intersect(Ray* r, Color* color, int level) {
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
    if (discriminant < 0) return -1;
    
    double t1 = (-Bq - sqrt(discriminant)) / (2 * Aq);
    double t2 = (-Bq + sqrt(discriminant)) / (2 * Aq);
    
    double t = -1;
    if (t1 > 0) t = t1;
    if (t2 > 0 && (t2 < t1 || t < 0)) t = t2;
    if (t < 0) return -1;
    
    Vector3D intersectionPoint = r->start + r->dir * t;
    if (!isWithinReferenceCube(intersectionPoint)) return -1;
    
    if (level == 0) return t;
    
    Vector3D normal = getNormal(intersectionPoint);
    *color = calculateLighting(r, intersectionPoint, normal, this->color, this, level);
    return t;
}

bool GeneralQuadric::isWithinReferenceCube(const Vector3D& point) const {
    if (width > 0 && (point.x < reference_point.x || point.x > reference_point.x + width))
        return false;
    if (height > 0 && (point.y < reference_point.y || point.y > reference_point.y + height))
        return false;
    if (length > 0 && (point.z < reference_point.z || point.z > reference_point.z + length))
        return false;
    return true;
}

Vector3D GeneralQuadric::getNormal(const Vector3D& point) const {
    double nx = 2 * coefficients[0] * point.x + coefficients[3] * point.y + coefficients[4] * point.z + coefficients[6];
    double ny = 2 * coefficients[1] * point.y + coefficients[3] * point.x + coefficients[5] * point.z + coefficients[7];
    double nz = 2 * coefficients[2] * point.z + coefficients[4] * point.x + coefficients[5] * point.y + coefficients[8];
    return Vector3D(nx, ny, nz).normalize();
}

Color GeneralQuadric::getColorAt(const Vector3D& point) const {
    return color;
}

// === LIGHT IMPLEMENTATIONS ===
void PointLight::draw() {
    glPushMatrix();
    glTranslatef(light_pos.x, light_pos.y, light_pos.z);
    glColor3f(color.r, color.g, color.b);
    glutSolidSphere(0.5, 10, 10);
    glPopMatrix();
}

void SpotLight::draw() {
    pointLight.draw();
}

/**
 * Initialize OpenGL settings
 * Sets up background color and enables depth testing
 */

/**
 * Main display function
 * Sets up the camera and renders visible objects
 */

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
    gluPerspective(45.0f, aspect, 0.1f, 1000.0f);
}

/**
 * Keyboard input handler for standard keys
 * Manages camera position, object visibility, and program exit
 */
void keyboardListener(unsigned char key, int x, int y)
{

    // Calculate view direction vector
    double lx = lookX - eyeX;
    double lz = lookZ - eyeZ;
    double s;

    switch (key)
    {
    case '0': // Capture image
        capture();
        break;
    case '1': // Rotate left
        {
            Vector3D newLook = Vector3D(lookX - eyeX, lookY - eyeY, lookZ - eyeZ);
            double cosAngle = cos(rotationSpeed);
            double sinAngle = sin(rotationSpeed);
            lookX = eyeX + newLook.x * cosAngle + s * sinAngle;
            lookY = eyeY + newLook.y * cosAngle + s * sinAngle;
            lookZ = eyeZ + newLook.z * cosAngle + s * sinAngle;
        }
        break;
    case '2': // Rotate right
        {
            Vector3D newLook = Vector3D(lookX - eyeX, lookY - eyeY, lookZ - eyeZ);
            double cosAngle = cos(-rotationSpeed);
            double sinAngle = sin(-rotationSpeed);
            lookX = eyeX + newLook.x * cosAngle + s * sinAngle;
            lookY = eyeY + newLook.y * cosAngle + s * sinAngle;
            lookZ = eyeZ + newLook.z * cosAngle + s * sinAngle;
        }
        break;
    case '3': // Look up
        {
            Vector3D newLook = Vector3D(lookX - eyeX, lookY - eyeY, lookZ - eyeZ);
            double cosAngle = cos(rotationSpeed);
            double sinAngle = sin(rotationSpeed);
            lookX = eyeX + newLook.x * cosAngle + s * sinAngle;
            lookY = eyeY + newLook.y * cosAngle + s * sinAngle;
            lookZ = eyeZ + newLook.z * cosAngle + s * sinAngle;
        }
        break;
    case '4': // Look down
        {
            Vector3D newLook = Vector3D(lookX - eyeX, lookY - eyeY, lookZ - eyeZ);
            double cosAngle = cos(-rotationSpeed);
            double sinAngle = sin(-rotationSpeed);
            lookX = eyeX + newLook.x * cosAngle + s * sinAngle;
            lookY = eyeY + newLook.y * cosAngle + s * sinAngle;
            lookZ = eyeZ + newLook.z * cosAngle + s * sinAngle;
        }
        break;
    case '5': // Tilt clockwise
        {
            double cosAngle = cos(-rotationSpeed);
            double sinAngle = sin(-rotationSpeed);
            double newUpX = upX * cosAngle + s * sinAngle;
            double newUpY = upY * cosAngle + s * sinAngle;
            double newUpZ = upZ * cosAngle + s * sinAngle;
            upX = newUpX; upY = newUpY; upZ = newUpZ;
        }
        break;
    case '6': // Tilt counter-clockwise
        {
            double cosAngle = cos(rotationSpeed);
            double sinAngle = sin(rotationSpeed);
            double newUpX = upX * cosAngle + s * sinAngle;
            double newUpY = upY * cosAngle + s * sinAngle;
            double newUpZ = upZ * cosAngle + s * sinAngle;
            upX = newUpX; upY = newUpY; upZ = newUpZ;
        }
        break;
    case 'a': // Toggle axes
        isAxes = !isAxes;
        break;
    case 27: // ESC key
        exit(0);
        break;
    }
    glutPostRedisplay();
}

void specialKeyListener(int key, int x, int y) {
    Vector3D l = Vector3D(lookX - eyeX, lookY - eyeY, lookZ - eyeZ).normalize();
    Vector3D r = l.cross(Vector3D(upX, upY, upZ)).normalize();
    Vector3D u = r.cross(l).normalize();
    
    switch (key) {
    case GLUT_KEY_UP: // Move forward
        eyeX += l.x * movementSpeed;
        eyeY += l.y * movementSpeed;
        eyeZ += l.z * movementSpeed;
        lookX += l.x * movementSpeed;
        lookY += l.y * movementSpeed;
        lookZ += l.z * movementSpeed;
        break;
    case GLUT_KEY_DOWN: // Move backward
        eyeX -= l.x * movementSpeed;
        eyeY -= l.y * movementSpeed;
        eyeZ -= l.z * movementSpeed;
        lookX -= l.x * movementSpeed;
        lookY -= l.y * movementSpeed;
        lookZ -= l.z * movementSpeed;
        break;
    case GLUT_KEY_LEFT: // Move left
        eyeX -= r.x * movementSpeed;
        eyeY -= r.y * movementSpeed;
        eyeZ -= r.z * movementSpeed;
        lookX -= r.x * movementSpeed;
        lookY -= r.y * movementSpeed;
        lookZ -= r.z * movementSpeed;
        break;
    case GLUT_KEY_RIGHT: // Move right
        eyeX += r.x * movementSpeed;
        eyeY += r.y * movementSpeed;
        eyeZ += r.z * movementSpeed;
        lookX += r.x * movementSpeed;
        lookY += r.y * movementSpeed;
        lookZ += r.z * movementSpeed;
        break;
    case GLUT_KEY_PAGE_UP: // Move up
        eyeX += u.x * movementSpeed;
        eyeY += u.y * movementSpeed;
        eyeZ += u.z * movementSpeed;
        lookX += u.x * movementSpeed;
        lookY += u.y * movementSpeed;
        lookZ += u.z * movementSpeed;
        break;
    case GLUT_KEY_PAGE_DOWN: // Move down
        eyeX -= u.x * movementSpeed;
        eyeY -= u.y * movementSpeed;
        eyeZ -= u.z * movementSpeed;
        lookX -= u.x * movementSpeed;
        lookY -= u.y * movementSpeed;
        lookZ -= u.z * movementSpeed;
        break;
    }
    glutPostRedisplay();
}

void initGL() {
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, 1, 0.1, 1000.0);
}

void freeMemory() {
    for (Object* obj : objects) {
        delete obj;
    }
    objects.clear();
    pointLights.clear();
    spotLights.clear();
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 800);
    glutInitWindowPosition(50, 50);
    glutCreateWindow("Ray Tracing Assignment - CSE410");
    
    // Initialize OpenGL settings
    initGL();
    
    // Load scene data
    loadData();
    
    // Register callbacks
    glutDisplayFunc(display);
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    
    std::cout << "=== Ray Tracing Assignment ===" << std::endl;
    std::cout << "Controls:" << std::endl;
    std::cout << "Arrow Keys: Move camera forward/backward/left/right" << std::endl;
    std::cout << "Page Up/Down: Move camera up/down" << std::endl;
    std::cout << "1/2: Rotate left/right" << std::endl;
    std::cout << "3/4: Look up/down" << std::endl;
    std::cout << "5/6: Tilt clockwise/counter-clockwise" << std::endl;
    std::cout << "0: Capture image" << std::endl;
    std::cout << "a: Toggle axes" << std::endl;
    std::cout << "ESC: Exit" << std::endl;
    
    // Enter the GLUT event loop
    glutMainLoop();
    
    // Cleanup
    freeMemory();
    
    return 0;
}

// === LIGHTING CALCULATION ===
Color calculateLighting(Ray* ray, const Vector3D& point, const Vector3D& normal, 
                       const Color& objectColor, const Object* obj, int level) {
    Color finalColor = objectColor * obj->coEfficients[0]; // Ambient component
    
    // Point lights
    for (const auto& light : pointLights) {
        Vector3D lightDir = (light.light_pos - point).normalize();
        
        // Check for shadows
        Ray shadowRay(point + normal * 0.001, lightDir); // Slightly offset to avoid self-intersection
        bool inShadow = false;
        for (Object* object : objects) {
            if (object != obj && object->intersect(&shadowRay, nullptr, 0) > 0) {
                inShadow = true;
                break;
            }
        }
        
        if (!inShadow) {
            // Diffuse component
            double lambertValue = std::max(0.0, normal.dot(lightDir));
            finalColor = finalColor + light.color * objectColor * obj->coEfficients[1] * lambertValue;
            
            // Specular component
            Vector3D reflectedRay = (normal * 2.0 * normal.dot(lightDir) - lightDir).normalize();
            Vector3D viewDir = (ray->start - point).normalize();
            double phongValue = std::max(0.0, reflectedRay.dot(viewDir));
            phongValue = pow(phongValue, obj->shine);
            finalColor = finalColor + light.color * obj->coEfficients[2] * phongValue;
        }
    }
    
    // Spot lights
    for (const auto& spotlight : spotLights) {
        Vector3D lightDir = (spotlight.pointLight.light_pos - point).normalize();
        Vector3D spotDir = spotlight.light_dir.normalize();
        
        // Check if point is within spotlight cone
        double spotAngle = acos((-lightDir).dot(spotDir)) * 180.0 / M_PI;
        if (spotAngle <= spotlight.cutoff_angle) {
            // Check for shadows
            Ray shadowRay(point + normal * 0.001, lightDir);
            bool inShadow = false;
            for (Object* object : objects) {
                if (object != obj && object->intersect(&shadowRay, nullptr, 0) > 0) {
                    inShadow = true;
                    break;
                }
            }
            
            if (!inShadow) {
                // Diffuse component
                double lambertValue = std::max(0.0, normal.dot(lightDir));
                finalColor = finalColor + spotlight.pointLight.color * objectColor * obj->coEfficients[1] * lambertValue;
                
                // Specular component
                Vector3D reflectedRay = (normal * 2.0 * normal.dot(lightDir) - lightDir).normalize();
                Vector3D viewDir = (ray->start - point).normalize();
                double phongValue = std::max(0.0, reflectedRay.dot(viewDir));
                phongValue = pow(phongValue, obj->shine);
                finalColor = finalColor + spotlight.pointLight.color * obj->coEfficients[2] * phongValue;
            }
        }
    }
    
    // Reflection
    if (level < recursionLevel && obj->coEfficients[3] > 0) {
        Vector3D reflectedDir = (ray->dir - normal * 2.0 * ray->dir.dot(normal)).normalize();
        Ray reflectedRay(point + normal * 0.001, reflectedDir);
        
        // Find nearest intersecting object
        double tMin = -1;
        Object* nearestObject = nullptr;
        for (Object* object : objects) {
            double t = object->intersect(&reflectedRay, nullptr, 0);
            if (t > 0 && (tMin < 0 || t < tMin)) {
                tMin = t;
                nearestObject = object;
            }
        }
        
        if (nearestObject) {
            Color reflectedColor;
            nearestObject->intersect(&reflectedRay, &reflectedColor, level + 1);
            finalColor = finalColor + reflectedColor * obj->coEfficients[3];
        }
    }
    
    finalColor.clamp();
    return finalColor;
}

// === LOAD DATA FROM scene.txt ===
void loadData() {
    std::ifstream file("scene.txt");
    if (!file.is_open()) {
        std::cout << "Error: Could not open scene.txt" << std::endl;
        return;
    }
    
    file >> recursionLevel >> imageWidth;
    imageHeight = imageWidth; // Square image
    
    int numObjects;
    file >> numObjects;
    
    for (int i = 0; i < numObjects; i++) {
        std::string objectType;
        file >> objectType;
        
        if (objectType == "sphere") {
            double x, y, z, radius;
            file >> x >> y >> z >> radius;
            
            double r, g, b;
            file >> r >> g >> b;
            
            double amb, diff, spec, refl;
            file >> amb >> diff >> spec >> refl;
            
            int shine;
            file >> shine;
            
            Sphere* sphere = new Sphere(Vector3D(x, y, z), radius);
            sphere->setColor(r, g, b);
            sphere->setCoEfficients(amb, diff, spec, refl);
            sphere->setShine(shine);
            objects.push_back(sphere);
        }
        else if (objectType == "triangle") {
            double x1, y1, z1, x2, y2, z2, x3, y3, z3;
            file >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
            
            double r, g, b;
            file >> r >> g >> b;
            
            double amb, diff, spec, refl;
            file >> amb >> diff >> spec >> refl;
            
            int shine;
            file >> shine;
            
            Triangle* triangle = new Triangle(Vector3D(x1, y1, z1), Vector3D(x2, y2, z2), Vector3D(x3, y3, z3));
            triangle->setColor(r, g, b);
            triangle->setCoEfficients(amb, diff, spec, refl);
            triangle->setShine(shine);
            objects.push_back(triangle);
        }
        else if (objectType == "general") {
            double A, B, C, D, E, F, G, H, I, J;
            file >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
            
            double refX, refY, refZ, length, width, height;
            file >> refX >> refY >> refZ >> length >> width >> height;
            
            double r, g, b;
            file >> r >> g >> b;
            
            double amb, diff, spec, refl;
            file >> amb >> diff >> spec >> refl;
            
            int shine;
            file >> shine;
            
            GeneralQuadric* quadric = new GeneralQuadric(A, B, C, D, E, F, G, H, I, J);
            quadric->reference_point = Vector3D(refX, refY, refZ);
            quadric->length = length;
            quadric->width = width;
            quadric->height = height;
            quadric->setColor(r, g, b);
            quadric->setCoEfficients(amb, diff, spec, refl);
            quadric->setShine(shine);
            objects.push_back(quadric);
        }
    }
    
    // Add floor
    floor_obj = new Floor(1000, 20);
    floor_obj->setColor(0.5, 0.5, 0.5);
    floor_obj->setCoEfficients(0.4, 0.2, 0.1, 0.2);
    floor_obj->setShine(5);
    objects.push_back(floor_obj);
    
    // Load point lights
    int numPointLights;
    file >> numPointLights;
    
    for (int i = 0; i < numPointLights; i++) {
        double x, y, z, r, g, b;
        file >> x >> y >> z >> r >> g >> b;
        pointLights.push_back(PointLight(Vector3D(x, y, z), Color(r, g, b)));
    }
    
    // Load spot lights
    int numSpotLights;
    file >> numSpotLights;
    
    for (int i = 0; i < numSpotLights; i++) {
        double x, y, z, r, g, b, dx, dy, dz, cutoff;
        file >> x >> y >> z >> r >> g >> b >> dx >> dy >> dz >> cutoff;
        spotLights.push_back(SpotLight(Vector3D(x, y, z), Color(r, g, b), Vector3D(dx, dy, dz), cutoff));
    }
    
    file.close();
    std::cout << "Scene loaded successfully!" << std::endl;
}

// === CAPTURE FUNCTION ===
void capture() {
    std::cout << "Capturing image..." << std::endl;
    
    bitmap_image image(imageWidth, imageHeight);
    image.set_all_channels(0, 0, 0); // Black background
    
    double planeDistance = (imageHeight / 2.0) / tan((M_PI / 180.0) * 45.0 / 2.0); // 45 degree FOV
    
    Vector3D eye(eyeX, eyeY, eyeZ);
    Vector3D look(lookX, lookY, lookZ);
    Vector3D up(upX, upY, upZ);
    
    Vector3D l = (look - eye).normalize();
    Vector3D r = l.cross(up).normalize();
    Vector3D u = r.cross(l).normalize();
    
    Vector3D topleft = eye + l * planeDistance - r * (imageWidth / 2.0) + u * (imageHeight / 2.0);
    
    double du = 1.0;
    double dv = 1.0;
    
    topleft = topleft + r * (0.5 * du) - u * (0.5 * dv);
    
    for (int i = 0; i < imageWidth; i++) {
        for (int j = 0; j < imageHeight; j++) {
            Vector3D curPixel = topleft + r * (i * du) - u * (j * dv);
            Ray ray(eye, (curPixel - eye).normalize());
            
            // Find nearest intersecting object
            double tMin = -1;
            Object* nearestObject = nullptr;
            for (Object* object : objects) {
                double t = object->intersect(&ray, nullptr, 0);
                if (t > 0 && (tMin < 0 || t < tMin)) {
                    tMin = t;
                    nearestObject = object;
                }
            }
            
            if (nearestObject) {
                Color pixelColor;
                nearestObject->intersect(&ray, &pixelColor, 1);
                image.set_pixel(i, j, (unsigned char)(pixelColor.r * 255), 
                                     (unsigned char)(pixelColor.g * 255), 
                                     (unsigned char)(pixelColor.b * 255));
            }
        }
    }
    
    std::string filename = "Output_1" + std::to_string(captureCount) + ".bmp";
    image.save_image(filename);
    captureCount++;
    std::cout << "Image saved as " << filename << std::endl;
}

// === OPENGL FUNCTIONS ===
void drawAxes() {
    glLineWidth(3);
    glBegin(GL_LINES);
    
    // X-axis (Red)
    glColor3f(1, 0, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(100, 0, 0);
    
    // Y-axis (Green)
    glColor3f(0, 1, 0);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 100, 0);
    
    // Z-axis (Blue)
    glColor3f(0, 0, 1);
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 100);
    
    glEnd();
}

// Simple display function for camera demo
void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    // Set camera position
    gluLookAt(eyeX, eyeY, eyeZ,
              lookX, lookY, lookZ,
              upX, upY, upZ);
    
    // Draw axes if enabled
    if (isAxes)
        drawAxes();

    // Draw all objects
    for (Object* obj : objects) {
        obj->draw();
    }
    
    // Draw lights
    for (const auto& light : pointLights) {
        const_cast<PointLight&>(light).draw();
    }
    
    for (const auto& spotlight : spotLights) {
        const_cast<SpotLight&>(spotlight).draw();
    }
    
    glutSwapBuffers();
}

void reshapeListener(GLsizei width, GLsizei height) {
    if (height == 0) height = 1;
    GLfloat aspect = (GLfloat)width / (GLfloat)height;
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, aspect, 0.1f, 1000.0f);
}

void keyboardListener(unsigned char key, int x, int y) {
    Vector3D l = Vector3D(lookX - eyeX, lookY - eyeY, lookZ - eyeZ).normalize();
    Vector3D r = l.cross(Vector3D(upX, upY, upZ)).normalize();
    Vector3D u = r.cross(l).normalize();
    
    switch (key) {
    case '0': // Capture image
        capture();
        break;
    case '1': // Rotate left
        {
            Vector3D newLook = Vector3D(lookX - eyeX, lookY - eyeY, lookZ - eyeZ);
            double cosAngle = cos(rotationSpeed);
            double sinAngle = sin(rotationSpeed);
            lookX = eyeX + newLook.x * cosAngle + r.x * sinAngle;
            lookY = eyeY + newLook.y * cosAngle + r.y * sinAngle;
            lookZ = eyeZ + newLook.z * cosAngle + r.z * sinAngle;
        }
        break;
    case '2': // Rotate right
        {
            Vector3D newLook = Vector3D(lookX - eyeX, lookY - eyeY, lookZ - eyeZ);
            double cosAngle = cos(-rotationSpeed);
            double sinAngle = sin(-rotationSpeed);
            lookX = eyeX + newLook.x * cosAngle + r.x * sinAngle;
            lookY = eyeY + newLook.y * cosAngle + r.y * sinAngle;
            lookZ = eyeZ + newLook.z * cosAngle + r.z * sinAngle;
        }
        break;
    case '3': // Look up
        {
            Vector3D newLook = Vector3D(lookX - eyeX, lookY - eyeY, lookZ - eyeZ);
            double cosAngle = cos(rotationSpeed);
            double sinAngle = sin(rotationSpeed);
            lookX = eyeX + newLook.x * cosAngle + u.x * sinAngle;
            lookY = eyeY + newLook.y * cosAngle + u.y * sinAngle;
            lookZ = eyeZ + newLook.z * cosAngle + u.z * sinAngle;
        }
        break;
    case '4': // Look down
        {
            Vector3D newLook = Vector3D(lookX - eyeX, lookY - eyeY, lookZ - eyeZ);
            double cosAngle = cos(-rotationSpeed);
            double sinAngle = sin(-rotationSpeed);
            lookX = eyeX + newLook.x * cosAngle + u.x * sinAngle;
            lookY = eyeY + newLook.y * cosAngle + u.y * sinAngle;
            lookZ = eyeZ + newLook.z * cosAngle + u.z * sinAngle;
        }
        break;
    case '5': // Tilt clockwise
        {
            double cosAngle = cos(-rotationSpeed);
            double sinAngle = sin(-rotationSpeed);
            double newUpX = upX * cosAngle + r.x * sinAngle;
            double newUpY = upY * cosAngle + r.y * sinAngle;
            double newUpZ = upZ * cosAngle + r.z * sinAngle;
            upX = newUpX; upY = newUpY; upZ = newUpZ;
        }
        break;
    case '6': // Tilt counter-clockwise
        {
            double cosAngle = cos(rotationSpeed);
            double sinAngle = sin(rotationSpeed);
            double newUpX = upX * cosAngle + r.x * sinAngle;
            double newUpY = upY * cosAngle + r.y * sinAngle;
            double newUpZ = upZ * cosAngle + r.z * sinAngle;
            upX = newUpX; upY = newUpY; upZ = newUpZ;
        }
        break;
    case 'a': // Toggle axes
        isAxes = !isAxes;
        break;
    case 27: // ESC key
        exit(0);
        break;
    }
    glutPostRedisplay();
}

void specialKeyListener(int key, int x, int y) {
    Vector3D l = Vector3D(lookX - eyeX, lookY - eyeY, lookZ - eyeZ).normalize();
    Vector3D r = l.cross(Vector3D(upX, upY, upZ)).normalize();
    Vector3D u = r.cross(l).normalize();
    
    switch (key) {
    case GLUT_KEY_UP: // Move forward
        eyeX += l.x * movementSpeed;
        eyeY += l.y * movementSpeed;
        eyeZ += l.z * movementSpeed;
        lookX += l.x * movementSpeed;
        lookY += l.y * movementSpeed;
        lookZ += l.z * movementSpeed;
        break;
    case GLUT_KEY_DOWN: // Move backward
        eyeX -= l.x * movementSpeed;
        eyeY -= l.y * movementSpeed;
        eyeZ -= l.z * movementSpeed;
        lookX -= l.x * movementSpeed;
        lookY -= l.y * movementSpeed;
        lookZ -= l.z * movementSpeed;
        break;
    case GLUT_KEY_LEFT: // Move left
        eyeX -= r.x * movementSpeed;
        eyeY -= r.y * movementSpeed;
        eyeZ -= r.z * movementSpeed;
        lookX -= r.x * movementSpeed;
        lookY -= r.y * movementSpeed;
        lookZ -= r.z * movementSpeed;
        break;
    case GLUT_KEY_RIGHT: // Move right
        eyeX += r.x * movementSpeed;
        eyeY += r.y * movementSpeed;
        eyeZ += r.z * movementSpeed;
        lookX += r.x * movementSpeed;
        lookY += r.y * movementSpeed;
        lookZ += r.z * movementSpeed;
        break;
    case GLUT_KEY_PAGE_UP: // Move up
        eyeX += u.x * movementSpeed;
        eyeY += u.y * movementSpeed;
        eyeZ += u.z * movementSpeed;
        lookX += u.x * movementSpeed;
        lookY += u.y * movementSpeed;
        lookZ += u.z * movementSpeed;
        break;
    case GLUT_KEY_PAGE_DOWN: // Move down
        eyeX -= u.x * movementSpeed;
        eyeY -= u.y * movementSpeed;
        eyeZ -= u.z * movementSpeed;
        lookX -= u.x * movementSpeed;
        lookY -= u.y * movementSpeed;
        lookZ -= u.z * movementSpeed;
        break;
    }
    glutPostRedisplay();
}

void initGL() {
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, 1, 0.1, 1000.0);
}

void freeMemory() {
    for (Object* obj : objects) {
        delete obj;
    }
    objects.clear();
    pointLights.clear();
    spotLights.clear();
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 800);
    glutInitWindowPosition(50, 50);
    glutCreateWindow("Ray Tracing Assignment - CSE410");
    
    // Initialize OpenGL settings
    initGL();
    
    // Load scene data
    loadData();
    
    // Register callbacks
    glutDisplayFunc(display);
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    
    std::cout << "=== Ray Tracing Assignment ===" << std::endl;
    std::cout << "Controls:" << std::endl;
    std::cout << "Arrow Keys: Move camera forward/backward/left/right" << std::endl;
    std::cout << "Page Up/Down: Move camera up/down" << std::endl;
    std::cout << "1/2: Rotate left/right" << std::endl;
    std::cout << "3/4: Look up/down" << std::endl;
    std::cout << "5/6: Tilt clockwise/counter-clockwise" << std::endl;
    std::cout << "0: Capture image" << std::endl;
    std::cout << "a: Toggle axes" << std::endl;
    std::cout << "ESC: Exit" << std::endl;
    
    // Enter the GLUT event loop
    glutMainLoop();
    
    // Cleanup
    freeMemory();
    
    return 0;
}