#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0, 255,   0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;

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

// returns the minimum x and y in the first vector
//     and the maximum x and y in the second vector
vec2* boundingBox(vec2& t0, vec2& t1, vec2& t2, double maxWidth, double maxHeight) {
    vec2* bbx = new vec2[2];

    bbx[0].x = std::min(t0.x, std::min(t1.x, t2.x));
    bbx[0].y = std::min(t0.y, std::min(t1.y, t2.y));

    bbx[1].x = std::min(maxWidth, std::max(t0.x, std::max(t1.x, t2.x)));
    bbx[1].y = std::min(maxHeight, std::max(t0.y, std::max(t1.y, t2.y)));

    return bbx;
}

void triangle(vec2 t0, vec2 t1, vec2 t2, TGAImage &image, const TGAColor& color) { 
    
    vec2* bbx = boundingBox(t0, t1, t2, image.get_width() - 1, image.get_height() - 1);

    // iterate overall points in the bounding box to check if it is inside the triangle
    for (int x = bbx[0].x; x <= bbx[1].x; x++) {
        for (int y = bbx[0].y; y <= bbx[1].y; y++) {
            
            // create vectors from the triangle points
            vec3 xcomp, ycomp;

            // x compoent of the vectors
            xcomp[0] = t1.x - t0.x;
            xcomp[1] = t2.x - t0.x;
            xcomp[2] = t0.x - x;

            // y component of the vectors
            ycomp[0] = t1.y - t0.y;
            ycomp[1] = t2.y - t0.y;
            ycomp[2] = t0.y - y;

            // with the x and y compoonents in a vector, he have a system of two linear equations

            // use the cross product to compute the intersection of the two lines
            vec3 result = cross(xcomp, ycomp);

            if (std::abs(result[2]) < 1) continue;

            // simply the vector so that the z component is 1.
            double u = result.x / result.z;
            double v = result.y / result.z;

            // check if inside triangle
            if (u <= 1 && u >= 0
                && v <= 1 && v >= 0
                && u + v <= 1) {
                image.set(x, y, color);
            }
        }
    }

    delete bbx;
}

int main(int argc, char** argv) {

    if (argc == 2) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head/african_head.obj");
    }

    // light going into the screen (I think....)
    const vec3 light(0, 0, 1);

	TGAImage image(width, height, TGAImage::RGB);
    for ( int face = 0; face < model->nfaces(); face++ ) {
        
        // calculate the normal of the triangle
        vec3 normal = cross(model->vert(face, 1) - model->vert(face, 0), model->vert(face, 2) - model->vert(face, 0));
        normal.normalize();

        // use the normal to calculate how intense the light should be on the face
        float intensity = light * normal;;

        // backface culling - removing triangles that are facing the opposite direction
        if (intensity < 0) continue;

        vec2 screen_coords[3];
        for ( int vertex = 0; vertex < 3; vertex++ ) {            
            vec3 world_coords = model->vert(face, vertex);
            world_coords.x = (world_coords.x + 1.) / 2. * width;
            world_coords.y = (world_coords.y + 1.) / 2. * width;

            screen_coords[vertex] = vec2(world_coords.x, world_coords.y);
        }

        triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
    }

    //image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
	return 0;
}