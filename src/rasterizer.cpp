#include "rasterizer.h"

using namespace std;

namespace CGL {

    RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
        size_t width, size_t height,
        unsigned int sample_rate) {
        this->psm = psm;
        this->lsm = lsm;
        this->width = width;
        this->height = height;
        this->sample_rate = sample_rate;

        sample_buffer.resize(width * height * sample_rate, Color::White);
    }

    // Used by rasterize_point and rasterize_line
    void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
        // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
        // NOTE: You are not required to implement proper supersampling for points and lines
        // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
        sample_buffer[y * width + x] = c;
        // We create another function super_sample_fill_pixel() to replace this function!
    }

    // fill the pixel in supersample cases
    void RasterizerImp::super_sample_fill_pixel(float x, float y, int bias, Color c) {
        int sx = (int)floor(x);
        int sy = (int)floor(y);

        if (sx < 0 || sx >= width) return;
        if (sy < 0 || sy >= height) return;
        // The original pixel will occupy "sample_rate" number of positions instead of one.
        // bias is in [0, sample_rate - 1] 
        sample_buffer[sample_rate * (sy * width + sx) + bias] = c;
    }

    // Rasterize a point: simple example to help you start familiarizing
    // yourself with the starter code.
    void RasterizerImp::rasterize_point(float x, float y, Color color) {
        // fill in the nearest pixel
        int sx = (int)floor(x);
        int sy = (int)floor(y);

        // check bounds
        if (sx < 0 || sx >= width) return;
        if (sy < 0 || sy >= height) return;

        for (int bias = 0; bias < sample_rate; bias++) {
            super_sample_fill_pixel(sx, sy, bias, color);
        }
        return;
    }

    // Rasterize a line.
    void RasterizerImp::rasterize_line(float x0, float y0,
        float x1, float y1,
        Color color) {
        if (x0 > x1) {
            swap(x0, x1); swap(y0, y1);
        }

        float pt[] = { x0,y0 };
        float m = (y1 - y0) / (x1 - x0);
        float dpt[] = { 1,m };
        int steep = abs(m) > 1;
        if (steep) {
            dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
            dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
        }

        while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
            rasterize_point(pt[0], pt[1], color);
            pt[0] += dpt[0]; pt[1] += dpt[1];
        }
    }


    // Rasterize a triangle.
    void RasterizerImp::rasterize_triangle(float x0, float y0,
        float x1, float y1,
        float x2, float y2,
        Color color) {
        // TODO: Task 1: Implement basic triangle rasterization here, no supersampling

        // Ensure that these points are arranged in counterclockwise order. If not, exchange them
        if ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) < 0) {
            swap(x1, x2);
            swap(y1, y2);
        }

        // find the boundary of the triangle
        int xmax = (int)std::max({ x0, x1, x2 });
        int xmin = (int)std::min({ x0, x1, x2 });
        int ymax = (int)std::max({ y0, y1, y2 });
        int ymin = (int)std::min({ y0, y1, y2 });


        // TODO: Task 2: Update to implement super-sampled rasterization

        int scale = (int)sqrt(sample_rate);
        // sample the point 
        for (int y = ymin; y < ymax; y++) {
            for (int x = xmin; x < xmax; x++) {
                // break one pixel into scale * scale smaller pixels and detect whether they are in the triangle
                for (int j = 0; j < scale; j++) {
                    for (int i = 0; i < scale; i++) {
                        float px = x + (i + 0.5) / scale;
                        float py = y + (j + 0.5) / scale;
                        // Use L = V * N to detect the inside property
                        float L0 = -1 * (px - x0) * (y1 - y0) + (py - y0) * (x1 - x0);
                        float L1 = -1 * (px - x1) * (y2 - y1) + (py - y1) * (x2 - x1);
                        float L2 = -1 * (px - x2) * (y0 - y2) + (py - y2) * (x0 - x2);
                        if ((L0 >= 0) && (L1 >= 0) && (L2 >= 0)) {
                            super_sample_fill_pixel(x, y, i + j * 2, color);
                        }
                    }
                }
            }
        }


    }


    void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
        float x1, float y1, Color c1,
        float x2, float y2, Color c2)
    {
        // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
        // Hint: You can reuse code from rasterize_triangle
        if ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) < 0) {
            swap(x1, x2);
            swap(y1, y2);
        }

        int xmax = (int)std::max({ x0, x1, x2 });
        int xmin = (int)std::min({ x0, x1, x2 });
        int ymax = (int)std::max({ y0, y1, y2 });
        int ymin = (int)std::min({ y0, y1, y2 });


        int scale = (int)sqrt(sample_rate);
        for (int x = xmin; x < xmax; x++) {
            for (int y = ymin; y < ymax; y++) {
                // initialize the target color
                Color new_color = Color(0, 0, 0);
                for (int i = 0; i < scale; i++) {
                    for (int j = 0; j < scale; j++) {
                        float px = x + (i + 0.5) / scale;
                        float py = y + (j + 0.5) / scale;
                        float L0 = -1 * (px - x0) * (y1 - y0) + (py - y0) * (x1 - x0);
                        float L1 = -1 * (px - x1) * (y2 - y1) + (py - y1) * (x2 - x1);
                        float L2 = -1 * (px - x2) * (y0 - y2) + (py - y2) * (x0 - x2);
                        if ((L0 >= 0) && (L1 >= 0) && (L2 >= 0)) {
                            // calculate the barycentric coordinate
                            float alpha = (-(px - x1) * (y2 - y1) + (py - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                            float beta = (-(px - x2) * (y0 - y2) + (py - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                            float gamma = 1 - alpha - beta;
                            // use the weight to get the new color
                            new_color.r = alpha * c0.r + beta * c1.r + gamma * c2.r;
                            new_color.g = alpha * c0.g + beta * c1.g + gamma * c2.g;
                            new_color.b = alpha * c0.b + beta * c1.b + gamma * c2.b;
                            super_sample_fill_pixel(x, y, i + j * 2, new_color);
                        }
                    }
                }
            }
        }


    }


    void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
        float x1, float y1, float u1, float v1,
        float x2, float y2, float u2, float v2,
        Texture& tex)
    {
        // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
        // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
        // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
        if ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) < 0) {
            swap(x1, x2);
            swap(y1, y2);
        }

        // get the boundary and use ceil() / floor() to improve the accuracy
        int xmax = ceil(std::max({ x0, x1, x2 }));
        int xmin = floor(std::min({ x0, x1, x2 }));
        int ymax = ceil(std::max({ y0, y1, y2 }));
        int ymin = floor(std::min({ y0, y1, y2 }));


        int scale = (int)sqrt(sample_rate);
        // initialize the SampleParams
        SampleParams sp;
        sp.lsm = lsm;
        sp.psm = psm;
        Color c;

        // inintialize the transition matrix
        Matrix3x3 M(x0, x1, x2, y0, y1, y2, 1, 1, 1);
        Matrix3x3 inv_M = M.inv();

        // convert the point into Vector form
        Vector2D X0 = Vector2D(x0, y0);
        Vector2D X1 = Vector2D(x1, y1);
        Vector2D X2 = Vector2D(x2, y2);
        // triangle line
        Vector2D l0 = X1 - X0;
        Vector2D l1 = X2 - X1;
        Vector2D l2 = X0 - X2;
        // normal vector
        Vector2D N0 = Vector2D(-l0.y, l0.x);
        Vector2D N1 = Vector2D(-l1.y, l1.x);
        Vector2D N2 = Vector2D(-l2.y, l2.x);

        for (int y = ymin; y < ymax; y++) {
            for (int x = xmin; x < xmax; x++) {
                for (int j = 0; j < scale; j++) {
                    for (int i = 0; i < scale; i++) {
                        float px = (float)x + ((float)i + 0.5) / (float)scale;
                        float py = (float)y + ((float)j + 0.5) / (float)scale;
                        Vector2D p0 = Vector2D(px, py);
                        Vector2D V1 = p0 - X0;
                        Vector2D V2 = p0 - X1;
                        Vector2D V3 = p0 - X2;
                        // Use L = V * N to detect whether the point is inside the triangle
                        if ((dot(V1, N0) >= 0) && (dot(V2, N1) >= 0) && (dot(V3, N2) >= 0)) {
                            Vector3D p = Vector3D(px, py, 1);
                            // Use M^-1 * (x, y, 1) to get the barycentric coordinate
                            Vector3D weight_p = inv_M * p;
                            weight_p.z = 1 - weight_p.x - weight_p.y;
                            // Use the weight to calculate the uv coordinate of the sampled point
                            sp.p_uv = Vector2D(dot(weight_p, Vector3D(u0, u1, u2)), dot(weight_p, Vector3D(v0, v1, v2)));

                            Vector3D weight_p_dx = inv_M * Vector3D(px + 1., py, 1);
                            weight_p_dx.z = 1 - weight_p_dx.x - weight_p_dx.y;
                            // Since dx=1, then we just need to use subtraction to get d(u,v)/dx
                            sp.p_dx_uv = Vector2D(dot(weight_p_dx, Vector3D(u0, u1, u2)), dot(weight_p_dx, Vector3D(v0, v1, v2))) - sp.p_uv;

                            Vector3D weight_p_dy = inv_M * Vector3D(px, py + 1., 1);
                            weight_p_dy.z = 1 - weight_p_dy.x - weight_p_dy.y;
                            // Since dy=1, then we just need to use subtraction to get d(u,v)/dy
                            sp.p_dy_uv = Vector2D(dot(weight_p_dy, Vector3D(u0, u1, u2)), dot(weight_p_dy, Vector3D(v0, v1, v2))) - sp.p_uv;

                            // check the texture image and get the color
                            c = tex.sample(sp);
                            super_sample_fill_pixel(x, y, i + 2 * j, c);

                        }

                    }
                }

            }
        }

    }

    void RasterizerImp::set_sample_rate(unsigned int rate) {
        // TODO: Task 2: You may want to update this function for supersampling support

        this->sample_rate = rate;


        this->sample_buffer.resize(width * height * this->sample_rate, Color::White);
    }


    void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
        size_t width, size_t height)
    {
        // TODO: Task 2: You may want to update this function for supersampling support

        this->width = width;
        this->height = height;
        this->rgb_framebuffer_target = rgb_framebuffer;


        this->sample_buffer.resize(width * height * this->sample_rate, Color::White);
    }


    void RasterizerImp::clear_buffers() {
        std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
        std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
    }


    // This function is called at the end of rasterizing all elements of the
    // SVG file.  If you use a supersample buffer to rasterize SVG elements
    // for antialising, you could use this call to fill the target framebuffer
    // pixels from the supersample buffer data.
    //
    void RasterizerImp::resolve_to_framebuffer() {
        // TODO: Task 2: You will likely want to update this function for supersampling support

        int scale = (int)sqrt(sample_rate);
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                // Get all the small pixels related to the position, sum and average their colors to get a new color.
                Color colorsum = Color(0, 0, 0);
                for (int j = 0; j < scale; j++) {
                    for (int i = 0; i < scale; i++) {
                        colorsum += sample_buffer[sample_rate * (y * width + x) + i + j * 2];
                    }
                }
                colorsum = colorsum * (1.0 / sample_rate);
                for (int k = 0; k < 3; ++k) {
                    this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&colorsum.r)[k] * 255;
                }
            }
        }


    }

    Rasterizer::~Rasterizer() { }


    
}// CGL
