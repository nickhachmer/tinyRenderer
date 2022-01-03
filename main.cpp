#include <vector>
#include <cmath>
#include <algorithm>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include <iostream>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0, 255,   0,   255);
Model* model = NULL;
const int width  = 800;
const int height = 800;
const int depth  = 255;

// I dont fully understand what is going on here...
mat<4,4> myviewport(int x, int y, int w, int h) {
    mat<4,4> m = mat<4,4>::identity();
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}

void line(int x0, int y0, int x1, int y1, TGAImage &image, const TGAColor &color) {
    
    // ensure width is greater than the height
    bool steep = false;
    if ( std::abs(y1 - y0) > std::abs(x1 - x0) ) {
        // if height > width then transpose the points
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    
    // ensure we go left to righ
    if ( x0 > x1 ) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    int dx = x1 - x0; 
    int dy = y1 - y0; 

    int yIncrementor = y1 > y0 ? 1 : -1; 
    float derror =  std::abs(dy) * 2;
    float error = 0;

    int y = y0;
    for ( float x = x0; x <= x1; x++ ) { 
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        } 

        error += derror;
        if (error > dx) {
            y += yIncrementor;
            error -= dx * 2;
        }
    }
}

void line(vec2 t0, vec2 t1, TGAImage &image, const TGAColor &color) {
    line(t0.x, t0.y, t1.x, t1.y, image, color);
}

vec3 barycentric(vec2 p0, vec2 p1, vec2 p2, vec2 P) {
    // create vectors from the triangle points
    vec3 xcomp, ycomp;

    // x compoent of the vectors
    xcomp[0] = p1.x - p0.x;
    xcomp[1] = p2.x - p0.x;
    xcomp[2] = p0.x - P.x;

    // y component of the vectors
    ycomp[0] = p1.y - p0.y;
    ycomp[1] = p2.y - p0.y;
    ycomp[2] = p0.y - P.y;

    // with the x and y compoonents in a vector, he have a system of two linear equations

    // use the cross product to compute the intersection of the two lines
    vec3 result = cross(xcomp, ycomp);

    // triangle is degenerate (all points lay on a line) if ..
    if (std::abs(result[2]) < 0.01) return {-1, -1, -1};

    // simply the vector so that the z component is 1.
    double u = result.x / result.z;
    double v = result.y / result.z;

    return {u, v, 1 - u - v};
}

vec3 barycentric(vec3 p0, vec3 p1, vec3 p2, vec3 P) {
    return barycentric({p0.x, p0.y}, {p1.x, p1.y}, {p2.x, p2.y}, {P.x, P.y});
}

void triangle(vec3* pts, vec2* uvs, float *zbuffer, TGAImage &image, const double intensity) { 
    
    // create bounding box
    //     the minimum x and y in the first vector
    //     and the maximum x and y in the second vector
    double maxWidth = image.get_width() - 1;
    double maxHeight = image.get_height() - 1;
    vec2 bbx[2];

    bbx[0].x = std::min({pts[0].x, pts[1].x, pts[2].x});
    bbx[0].y = std::min({pts[0].y, pts[1].y, pts[2].y});
    bbx[1].x = std::min(maxWidth, std::max({pts[0].x, pts[1].x, pts[2].x}));
    bbx[1].y = std::min(maxHeight, std::max({pts[0].y, pts[1].y, pts[2].y}));

// this directive tells to compiler to parallelize this block of code
#pragma omp parallel for
    // iterate overall points in the bounding box to check if it is inside the triangle
    for (int x = (int) bbx[0].x; x <= (int) bbx[1].x; x++) {
        for (int y = (int) bbx[0].y; y <= (int) bbx[1].y; y++) {
            
            // dont draw points that are off screen
            if (x < 0 || y < 0) continue;
            
            vec3 P = { (double) x, (double) y, 0};
            vec3 bary = barycentric(pts[0], pts[1], pts[2], P);

            // check if point P is outside the triangle - don't draw it if it is
            if (bary[0] < 0 || bary[1] < 0 || bary[2] < 0) continue;

            // calculate the z coordinate of point P from the triangle
            float z = pts[0].z * bary[0] + pts[1].z * bary[1] + pts[2].z * bary[2];
            
            // use the z buffer as a way to determine if a point lies on top of another
            // if a point has a large z value than previous point that has been drawn then 
            // we draw on top of the point and store the new maximum z value.
            int idx = x + width * y;
            if (z > zbuffer[idx]) {

                zbuffer[idx] = z;

                // use the barycentric coordinates find the corresponding point in the uv triangle
                vec2 uv = uvs[0] + (uvs[1] - uvs[0]) * bary[0] + (uvs[2] - uvs[0]) * bary[1];
                
                image.set(x, y, model->diffuse(uv) * intensity);
            }
        }
    }
}

int main(int argc, char** argv) {

    if (argc == 2) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head/african_head.obj");
    }

    // initialize the z buffer
    float* zbuffer = new float[width * height];
    for (int i = width * height; i > 0; i--) {
        zbuffer[i] = -std::numeric_limits<float>::max();
    }

    // light going into the screen (I think....)
    const vec3 light(0, 0, 1);
    
    // create our transform matrices
    mat<4,4> viewPort   = myviewport(width/8, height/8, width*3/4, height*3/4);
    mat<4,4> projection = mat<4,4>::identity();
    // set specific element to give us perspective projection
    projection[3][2] = - 1.f / 1;

	TGAImage image(width, height, TGAImage::RGB);
    for ( int face = 0; face < model->nfaces(); face++ ) {
        
        // calculate the normal of the triangle
        vec3 normal = cross(model->vert(face, 1) - model->vert(face, 0), model->vert(face, 2) - model->vert(face, 0));
        normal.normalize();

        // use the normal to calculate how intense the light should be on the face
        float intensity = light * normal;

        // backface culling - removing triangles that are facing the opposite direction
        if (intensity < 0) continue;

        // obtain the vertices of the triangle in screen space
        vec3 screen_coords[3];
        vec2 uvs[3];
        for ( int vertex = 0; vertex < 3; vertex++ ) {
            vec3 world_coords = model->vert(face, vertex);

            vec4 result = viewPort * projection * embed<4>(world_coords);
            result[0] = result[0] / result[3];
            result[1] = result[1] / result[3];
            result[2] = result[2] / result[3];

            screen_coords[vertex] = proj<3>(result);
            uvs[vertex] = model->uv(face, vertex);
        }

        // rasterize the triangle (draw the pixels on the image)
        triangle(screen_coords, uvs, zbuffer, image, intensity);        
    }

    //image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    delete zbuffer;
	return 0;
}