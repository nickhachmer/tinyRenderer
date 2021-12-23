#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 1000;

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

int main(int argc, char** argv) {
    if (argc == 2) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head/african_head.obj");
    }

	TGAImage image(width, height, TGAImage::RGB);
    for ( int face = 0; face < model->nfaces(); face++ ) {
        for ( int vertex = 0; vertex < 3; vertex++ ) {
            vec3 v0 = model->vert(face, vertex);
            vec3 v1 = model->vert(face, (vertex + 1) % 3);

            int x0 = (v0.x + 1.) / 2. * width;
            int y0 = (v0.y + 1.) / 2. * height;
            int x1 = (v1.x + 1.) / 2. * width;
            int y1 = (v1.y + 1.) / 2. * height;
            
            line(x0, y0, x1, y1, image, white);
        }
    }

    //image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
	return 0;
}