# Ray Tracing Assignment - Computer Graphics Course

## Project Overview

This repository contains a comprehensive ray tracing application developed for CSE410 Offline Assignment 3. The project implements a full-featured ray tracer with real-time OpenGL visualization, supporting multiple geometric primitives, advanced lighting models, texture mapping, and recursive reflections. The application is designed to load scenes from configuration files and render high-quality images using physically-based ray tracing techniques.

## Technical Architecture

### Core Components

The ray tracing engine is built around several key architectural components:

- **Ray-Object Intersection Engine**: Implements precise mathematical intersection algorithms for spheres, triangles, and general quadric surfaces
- **Phong Illumination Model**: Supports ambient, diffuse, and specular lighting with configurable material properties
- **Recursive Reflection System**: Handles multiple levels of reflections with configurable recursion depth
- **Texture Mapping System**: Supports both procedural checkerboard patterns and bitmap textures with UV coordinate mapping
- **Camera System**: Implements a flexible camera model with position, look-at, and up vectors for scene navigation

### Supported Geometric Primitives

#### Spheres
- **Definition**: Center point, radius, color, and material coefficients
- **Intersection**: Analytical ray-sphere intersection using quadratic formula
- **Normal Calculation**: Surface normal computed as normalized vector from center to intersection point

#### Triangles
- **Definition**: Three vertices defining the triangle plane
- **Intersection**: Möller-Trumbore ray-triangle intersection algorithm with barycentric coordinates
- **Normal Calculation**: Cross product of edge vectors for consistent surface orientation

#### General Quadric Surfaces
- **Definition**: Arbitrary quadric surfaces defined by 10 coefficients (A, B, C, D, E, F, G, H, I, J)
- **Support**: Ellipsoids, paraboloids, hyperboloids, and other second-degree surfaces
- **Clipping**: Axis-aligned bounding box clipping with configurable dimensions
- **Intersection**: Solves quadratic equation in ray parameter space

#### Floor Plane
- **Definition**: Infinite XY-plane at Z=0 with configurable tile size
- **Rendering Modes**: Procedural checkerboard pattern or texture mapping
- **Texture Mapping**: World coordinate to texture coordinate transformation with repeating patterns

### Lighting System

#### Point Light Sources
- **Attributes**: 3D position and RGB color intensity
- **Illumination**: Omnidirectional light emission with distance-based attenuation
- **Shadow Casting**: Ray casting from intersection points to light sources for shadow computation

#### Spotlight Sources
- **Attributes**: Position, direction vector, cutoff angle, and color
- **Cone Geometry**: Angular attenuation based on angle from central axis
- **Shadow Casting**: Directional shadow computation within cone boundaries

### Material Properties

Each object supports the following material coefficients:
- **Ambient Coefficient**: Controls ambient light contribution
- **Diffuse Coefficient**: Controls Lambertian diffuse reflection
- **Specular Coefficient**: Controls specular highlight intensity
- **Reflection Coefficient**: Controls mirror-like reflection strength
- **Shininess**: Controls specular highlight size (Phong exponent)

## File Structure and Organization

### Core Implementation Files

#### `2005090_Main.cpp` (1021 lines)
Primary application file containing:
- **OpenGL/GLUT Setup**: Window creation, viewport configuration, and rendering context initialization
- **Scene Loading**: Parser for `scene.txt` file format with error handling
- **Camera System**: Interactive camera controls with movement and rotation
- **Ray Tracing Integration**: Interface between OpenGL display and ray tracing engine
- **User Interface**: Keyboard input handling and real-time parameter adjustment
- **Image Capture**: BMP image output with configurable resolution

Key Functions:
- `loadData()`: Scene file parser and object instantiation
- `capture()`: Ray tracing image generation and file output
- `display()`: OpenGL rendering loop
- `keyboard()`: Input event handling
- `specialKeyListener()`: Arrow key navigation

#### `2005090_Header.hpp` (988 lines)
Comprehensive header file containing:
- **Vector3D Class**: 3D vector mathematics with dot product, cross product, and normalization
- **Color Class**: RGB color representation with utility functions
- **TextureManager Class**: Static texture loading and sampling using stb_image library
- **Ray Class**: Ray representation with origin and direction vectors
- **Object Hierarchy**: Abstract base class and concrete implementations
- **Light Classes**: Point light and spotlight implementations

Core Classes:
- `Vector3D`: Complete 3D vector arithmetic implementation
- `Object`: Abstract base class defining intersection and rendering interface
- `Sphere`: Analytical sphere intersection and normal calculation
- `Triangle`: Barycentric coordinate-based triangle intersection
- `GeneralQuadric`: Quadric surface intersection with bounding box clipping
- `Floor`: Infinite plane with checkerboard and texture support
- `PointLight`/`SpotLight`: Light source implementations with shadow casting

### Configuration and Assets

#### `scene.txt` (129 lines)
Scene description file format:
```
recursion_level image_size
num_objects
[object_definitions]
num_point_lights
[point_light_definitions]
num_spot_lights
[spot_light_definitions]
```

Object Definition Formats:
- **Sphere**: `sphere center_x center_y center_z radius r g b ambient diffuse specular reflection shininess`
- **Triangle**: `triangle x1 y1 z1 x2 y2 z2 x3 y3 z3 r g b ambient diffuse specular reflection shininess`
- **General Quadric**: `general A B C D E F G H I J ref_x ref_y ref_z length width height r g b ambient diffuse specular reflection shininess`

#### Texture Assets
- `texture1.bmp`, `texture2.bmp`: Primary texture files for floor mapping
- `sample_texture.bmp`, `sample_texture2.bmp`: Additional texture options
- `sample.bmp`: Reference texture file

#### External Libraries

#### `bitmap_image.hpp` (4925 lines)
- **Purpose**: Platform-independent BMP image reading and writing
- **Features**: 24-bit RGB bitmap support with pixel manipulation
- **License**: MIT License
- **Usage**: Image output for ray-traced results

#### `stb_image.h`
- **Purpose**: Universal image loading library
- **Supported Formats**: BMP, PNG, JPEG, TGA, HDR, PIC, PNM
- **Integration**: Texture loading for floor mapping
- **Implementation**: Single-header library with STB_IMAGE_IMPLEMENTATION

### Build Configuration

#### `a.sh`
Build script for Windows with MinGW:
```bash
g++ main.cpp -o demo.exe -lfreeglut -lglew32 -lopengl32 -lglu32; start demo.exe
```

#### Dependencies
- **FreeGLUT**: Cross-platform windowing and input handling
- **GLEW**: OpenGL extension wrangler library
- **OpenGL32**: Core OpenGL library
- **GLU32**: OpenGL utility library

### Output Files

#### Generated Images
- `Output_11.bmp`, `Output_12.bmp`, `Output_13.bmp`, `Output_14.bmp`: Sample ray-traced output images
- **Naming Convention**: Incremental numbering starting from Output_11.bmp
- **Resolution**: Configurable via scene file (default 2000x2000 pixels)

#### Compiled Executable
- `demo.exe`: Windows executable with embedded dependencies
- **Alternative**: Executable in `2005090/` subdirectory

## Installation and Setup

### Prerequisites

#### System Requirements
- **Operating System**: Windows 10/11 (primary), Linux, or macOS
- **Compiler**: GCC (MinGW on Windows) or Visual Studio 2019+
- **Memory**: Minimum 4GB RAM (8GB recommended for high-resolution rendering)
- **Graphics**: OpenGL 3.0+ compatible graphics card

#### Required Libraries
- **FreeGLUT**: Cross-platform windowing toolkit
- **GLEW**: OpenGL Extension Wrangler Library
- **OpenGL**: Graphics rendering pipeline
- **GLU**: OpenGL Utility Library

### Windows Installation

#### Method 1: Using MinGW (Recommended)
1. Install MinGW-w64 from [mingw-w64.org](https://www.mingw-w64.org/)
2. Install FreeGLUT development libraries
3. Install GLEW development libraries
4. Clone the repository:
   ```bash
   git clone https://github.com/TAR2003/Ray_Tracing_Assignment_Graphics_4-1.git
   cd Ray_Tracing_Assignment_Graphics_4-1
   ```
5. Build using the provided script:
   ```bash
   ./a.sh
   ```

#### Method 2: Manual Compilation
```bash
g++ 2005090_Main.cpp -o demo.exe -lfreeglut -lglew32 -lopengl32 -lglu32
```

### Linux Installation

#### Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install freeglut3-dev libglew-dev libgl1-mesa-dev libglu1-mesa-dev
git clone https://github.com/TAR2003/Ray_Tracing_Assignment_Graphics_4-1.git
cd Ray_Tracing_Assignment_Graphics_4-1
g++ 2005090_Main.cpp -o demo -lglut -lGLEW -lGL -lGLU
./demo
```

#### CentOS/RHEL/Fedora:
```bash
sudo yum install freeglut-devel glew-devel mesa-libGL-devel mesa-libGLU-devel
# or for newer versions:
sudo dnf install freeglut-devel glew-devel mesa-libGL-devel mesa-libGLU-devel
g++ 2005090_Main.cpp -o demo -lglut -lGLEW -lGL -lGLU
./demo
```

### macOS Installation

```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install freeglut glew

# Build
g++ 2005090_Main.cpp -o demo -framework GLUT -framework OpenGL -lGLEW
./demo
```

## Usage Guide

### Basic Operation

1. **Launch Application**: Run the compiled executable (`demo.exe` on Windows, `./demo` on Unix systems)
2. **Scene Loading**: The application automatically loads `scene.txt` from the current directory
3. **Navigation**: Use keyboard controls to explore the 3D scene
4. **Rendering**: Press '0' to capture a ray-traced image
5. **Exit**: Press ESC to close the application

### Keyboard Controls

#### Camera Movement
| Key | Function |
|-----|----------|
| ↑ | Move camera forward |
| ↓ | Move camera backward |
| ← | Move camera left |
| → | Move camera right |
| Page Up | Move camera up |
| Page Down | Move camera down |
| W | Move camera up (alternative) |
| S | Move camera down (alternative) |

#### Camera Rotation
| Key | Function |
|-----|----------|
| 1 | Rotate look-at point left |
| 2 | Rotate look-at point right |
| 3 | Rotate camera up |
| 4 | Rotate camera down |
| 5 | Tilt camera left |
| 6 | Tilt camera right |

#### Rendering and Display
| Key | Function |
|-----|----------|
| 0 | Capture ray-traced image |
| T | Toggle floor texture (checkerboard ↔ texture) |
| A | Toggle coordinate axes display |
| ESC | Exit application |

### Scene File Configuration

#### Scene File Format (`scene.txt`)

The scene configuration follows a structured format:

```
recursion_level image_size
number_of_objects
[object_definitions...]
number_of_point_lights
[point_light_definitions...]
number_of_spot_lights
[spot_light_definitions...]
```

#### Object Definitions

**Sphere Format:**
```
sphere
center_x center_y center_z
radius
color_r color_g color_b
ambient_coeff diffuse_coeff specular_coeff reflection_coeff
shininess
```

**Triangle Format:**
```
triangle
vertex1_x vertex1_y vertex1_z
vertex2_x vertex2_y vertex2_z
vertex3_x vertex3_y vertex3_z
color_r color_g color_b
ambient_coeff diffuse_coeff specular_coeff reflection_coeff
shininess
```

**General Quadric Format:**
```
general
A B C D E F G H I J
reference_x reference_y reference_z length width height
color_r color_g color_b
ambient_coeff diffuse_coeff specular_coeff reflection_coeff
shininess
```

#### Light Definitions

**Point Light Format:**
```
position_x position_y position_z
color_r color_g color_b
```

**Spotlight Format:**
```
position_x position_y position_z
color_r color_g color_b
direction_x direction_y direction_z
cutoff_angle_degrees
```

### Advanced Configuration

#### Material Properties

- **Ambient Coefficient** (0.0-1.0): Controls ambient light contribution
- **Diffuse Coefficient** (0.0-1.0): Controls Lambertian diffuse reflection
- **Specular Coefficient** (0.0-1.0): Controls specular highlight intensity
- **Reflection Coefficient** (0.0-1.0): Controls mirror reflection strength
- **Shininess** (1-100+): Controls specular highlight size (higher = sharper highlights)

#### Rendering Parameters

- **Recursion Level**: Maximum reflection depth (1-10 recommended)
- **Image Size**: Output image resolution in pixels (square images only)
- **Floor Width**: Physical extent of the checkerboard floor
- **Tile Width**: Size of individual checkerboard tiles

### Texture Management

#### Supported Texture Formats
- BMP (Windows Bitmap)
- PNG (Portable Network Graphics)
- JPEG/JPG (Joint Photographic Experts Group)
- TGA (Truevision Graphics Adapter)
- HDR (High Dynamic Range)
- PIC (Softimage PIC)
- PNM (Portable Anymap)

#### Texture Loading
The application automatically attempts to load the default texture file specified in the source code. Alternative textures can be loaded by modifying the `textureFileName` variable and recompiling.

#### Texture Coordinate Mapping
- Floor textures use world coordinate to UV coordinate transformation
- Texture coordinates are automatically scaled and repeated for optimal visual quality
- The texture scale factor can be adjusted in the `Floor::getColorAt()` method

## Ray Tracing Algorithm Implementation

### Core Ray Tracing Pipeline

#### 1. Camera Ray Generation
For each pixel in the output image:
- Calculate the pixel's position in world space using camera parameters
- Generate a ray from the camera position through the pixel
- Normalize the ray direction vector

#### 2. Scene Intersection Testing
For each ray:
- Test intersection with all objects in the scene
- Find the closest valid intersection point
- Store intersection distance, object reference, and surface normal

#### 3. Lighting Calculation
At each intersection point:
- Initialize color with ambient lighting component
- For each light source in the scene:
  - Cast shadow ray from intersection point to light source
  - If path is unobstructed, calculate diffuse and specular contributions
  - Apply Phong illumination model with material properties

#### 4. Recursive Reflection
If reflection coefficient > 0 and recursion depth < maximum:
- Calculate reflection ray direction using surface normal
- Recursively trace reflection ray
- Blend reflected color with surface color using reflection coefficient

#### 5. Color Accumulation and Output
- Combine ambient, diffuse, specular, and reflection components
- Clamp color values to valid range [0, 1]
- Convert to 8-bit RGB values and store in output image

### Mathematical Foundations

#### Ray-Sphere Intersection
Given ray R(t) = O + tD and sphere with center C and radius r:
- Solve quadratic equation: ||O + tD - C||² = r²
- Discriminant determines intersection existence
- Choose closest positive t value

#### Ray-Triangle Intersection
Using Möller-Trumbore algorithm:
- Express intersection point using barycentric coordinates
- Solve linear system for ray parameter t and barycentric coordinates (u, v)
- Valid intersection if t > 0, u ≥ 0, v ≥ 0, and u + v ≤ 1

#### Ray-Quadric Intersection
For general quadric Ax² + By² + Cz² + Dxy + Exz + Fyz + Gx + Hy + Iz + J = 0:
- Substitute ray equation R(t) = O + tD
- Solve resulting quadratic equation in t
- Apply bounding box clipping for finite quadric surfaces

#### Phong Illumination Model
Total color = Ambient + Diffuse + Specular:
- **Ambient**: Ka × Ia
- **Diffuse**: Kd × Id × max(0, N · L)
- **Specular**: Ks × Is × max(0, R · V)^n

Where:
- Ka, Kd, Ks: Material coefficients
- Ia, Id, Is: Light intensities
- N: Surface normal
- L: Light direction
- R: Reflection direction
- V: View direction
- n: Shininess exponent

## Performance Optimization

### Computational Complexity
- **Time Complexity**: O(pixels × objects × lights × recursion_depth)
- **Space Complexity**: O(objects + lights + image_size²)
- **Bottlenecks**: Intersection testing and recursive reflection calculations

### Optimization Strategies

#### 1. Early Ray Termination
- Discard rays with minimal contribution (low reflection coefficients)
- Limit recursion depth to prevent excessive computation
- Use epsilon values to handle floating-point precision issues

#### 2. Efficient Intersection Testing
- Bounding box checks before detailed intersection calculations
- Analytical solutions for primitive geometries
- Spatial data structures (future enhancement)

#### 3. Memory Management
- Dynamic allocation for scene objects with proper cleanup
- Static texture management to avoid redundant loading
- Efficient color array handling in intersection calculations

### Scalability Considerations
- Linear scaling with image resolution
- Quadratic scaling with scene complexity
- Exponential scaling with recursion depth

## Troubleshooting

### Common Issues

#### Compilation Errors

**Missing Libraries:**
```
error: GL/glut.h: No such file or directory
```
**Solution**: Install FreeGLUT development packages

**Linker Errors:**
```
undefined reference to `glutInit'
```
**Solution**: Ensure all required libraries are linked (-lglut -lGLEW -lGL -lGLU)

#### Runtime Issues

**Black Output Images:**
- Verify scene file format and syntax
- Check camera position and orientation
- Ensure objects are within the camera's view frustum
- Validate light source positions and intensities

**Application Crashes:**
- Check scene file path and permissions
- Verify texture file accessibility
- Monitor memory usage for large scenes
- Validate OpenGL context creation

**Poor Performance:**
- Reduce recursion level for faster rendering
- Decrease image resolution for testing
- Optimize scene complexity
- Consider simplified lighting models

#### Debug Mode

Enable debug output by uncommenting debug statements in the source code:
- Object loading verification
- Intersection calculation details
- Lighting computation breakdown
- Memory allocation tracking

### Testing Procedures

#### Unit Testing
- Individual object intersection algorithms
- Light calculation accuracy
- Texture coordinate mapping
- Camera transformation matrices

#### Integration Testing
- Complete scene rendering pipeline
- Multi-object scene handling
- Complex lighting scenarios
- Recursive reflection accuracy

#### Performance Testing
- Benchmark different scene configurations
- Memory usage profiling
- Rendering time analysis
- Quality vs. performance trade-offs

## Development Guidelines

### Code Organization

#### Design Patterns
- **Abstract Factory**: Object creation from scene descriptions
- **Strategy Pattern**: Different intersection algorithms per object type
- **Template Method**: Common intersection calculation framework
- **Singleton**: Texture manager for global texture access

#### Coding Standards
- Consistent naming conventions (camelCase for variables, PascalCase for classes)
- Comprehensive inline documentation
- Error handling with meaningful messages
- Memory management with RAII principles

### Extension Points

#### Adding New Geometric Primitives
1. Inherit from `Object` base class
2. Implement `intersect()` method with ray-object mathematics
3. Implement `getNormal()` for surface normal calculation
4. Implement `draw()` for OpenGL visualization
5. Add scene file parser support in `loadData()`

#### Implementing New Lighting Models
1. Extend lighting calculation in `Object::intersectPoint()`
2. Add new light source classes inheriting appropriate base class
3. Implement attenuation and distribution functions
4. Update scene file format for new light parameters

#### Advanced Rendering Features
- **Anti-aliasing**: Multiple rays per pixel with averaging
- **Global Illumination**: Monte Carlo ray tracing for indirect lighting
- **Volumetric Rendering**: Participating media and atmospheric effects
- **Motion Blur**: Temporal sampling for moving objects
- **Depth of Field**: Camera lens simulation with focal distance

### Testing Framework

#### Automated Testing
- Unit tests for mathematical functions
- Integration tests for rendering pipeline
- Performance benchmarks for optimization
- Visual regression tests for rendering quality

#### Manual Testing
- Interactive scene navigation
- Parameter adjustment verification
- Cross-platform compatibility
- User interface responsiveness

## Academic Context

### Course Integration
This project serves as a comprehensive assignment for computer graphics courses, covering:
- **Ray Tracing Fundamentals**: Intersection algorithms and lighting models
- **3D Mathematics**: Vector operations and coordinate transformations
- **Computer Graphics Pipeline**: Rasterization vs. ray tracing comparison
- **Software Engineering**: Large-scale C++ project organization

### Learning Objectives
Students gain practical experience in:
- Implementing mathematical algorithms in code
- Managing complex software projects
- Understanding graphics hardware and software interaction
- Optimizing performance-critical applications

### Assessment Criteria
- **Correctness**: Accurate implementation of ray tracing algorithms
- **Completeness**: Support for all required geometric primitives and lighting
- **Quality**: Code organization, documentation, and error handling
- **Performance**: Reasonable rendering times for typical scenes
- **Innovation**: Creative scene designs and optimization techniques

## Future Enhancements

### Planned Features
- **Acceleration Structures**: Bounding Volume Hierarchy (BVH) for faster intersection testing
- **Advanced Materials**: Bump mapping, normal mapping, and procedural textures
- **Post-Processing Effects**: Tone mapping, color correction, and filtering
- **Animation Support**: Keyframe interpolation and motion blur
- **GPU Acceleration**: CUDA or OpenCL implementation for parallel processing

### Research Directions
- **Path Tracing**: Unbiased global illumination algorithms
- **Bidirectional Path Tracing**: Improved convergence for complex lighting
- **Photon Mapping**: Caustics and subsurface scattering simulation
- **Real-Time Ray Tracing**: Hardware-accelerated ray tracing with RTX

## Contributing

### Development Workflow
1. Fork the repository
2. Create a feature branch
3. Implement changes with appropriate testing
4. Submit pull request with detailed description
5. Code review and integration

### Contribution Guidelines
- Follow existing code style and conventions
- Include comprehensive documentation for new features
- Provide test cases for algorithmic changes
- Maintain backwards compatibility when possible

## License and Credits

### License
This project is developed as an academic assignment and is intended for educational purposes. The code may be used and modified for non-commercial educational activities.

### Third-Party Libraries
- **bitmap_image.hpp**: MIT License - Platform Independent Bitmap Image Reader Writer Library by Arash Partow
- **stb_image.h**: Public Domain - Image loading library by Sean Barrett
- **FreeGLUT**: X-Consortium License - Cross-platform windowing toolkit
- **GLEW**: Modified BSD License - OpenGL Extension Wrangler Library

### Acknowledgments
- Course instructors and teaching assistants for project guidance
- Computer graphics community for algorithm references and implementation techniques
- Open source library developers for providing essential tools

### Contact Information
For questions, issues, or contributions related to this project, please refer to the repository's issue tracker or contact the course instructors through official academic channels.