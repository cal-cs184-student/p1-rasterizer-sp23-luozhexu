#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

    Color Texture::sample(const SampleParams& sp) {
        // TODO: Task 6: Fill this in.
        float level = get_level(sp);
        if (sp.lsm == L_ZERO) {
            if (sp.psm == P_NEAREST) {
                return sample_nearest(sp.p_uv, 0);
            }
            else if (sp.psm == P_LINEAR) {
                return sample_bilinear(sp.p_uv, 0);
            }

        }
        else if (sp.lsm == L_NEAREST) {
            int level2 = round(level);
            // return magenta for invalid level
            if (level2 > mipmap.size() || level < 0) {
                return Color(1, 0, 1);
            }
            if (sp.psm == P_NEAREST) {
                return sample_nearest(sp.p_uv, level2);
            }
            else if (sp.psm == P_LINEAR) {
                return sample_bilinear(sp.p_uv, level2);
            }
        }
        else if (sp.lsm == L_LINEAR) {
            int level_up = ceil(level);
            int level_down = floor(level);
            // return magenta for invalid level
            if (level_up > mipmap.size() || level_down < 0) {
                return Color(1, 0, 1);
            }
            if (sp.psm == P_NEAREST) {
                // round up/down to get the adjecent neighbors
                Color c1 = sample_nearest(sp.p_uv, level_up);
                Color c2 = sample_nearest(sp.p_uv, level_down);
                // use the differences of current level from these two levels as weights
                Color c = (level_up - level) * c1 + (level - level_down) * c2;
                return c;
            }
            else if (sp.psm == P_LINEAR) {
                Color c1 = sample_bilinear(sp.p_uv, level_up);
                Color c2 = sample_bilinear(sp.p_uv, level_down);
                Color c = (level_up - level) * c1 + (level - level_down) * c2;
                return c;
            }
        }

    }

    float Texture::get_level(const SampleParams& sp) {
        // TODO: Task 6: Fill this in.
        // scale the derivatives according to the texture size
        Vector2D p_tx_uv = Vector2D(sp.p_dx_uv.x * (width), sp.p_dx_uv.y * (height));
        Vector2D p_ty_uv = Vector2D(sp.p_dy_uv.x * (width), sp.p_dy_uv.y * (height));
        // get the bigger norm
        float L = std::max(p_tx_uv.norm(), p_ty_uv.norm());
        float D = log2(L);
        return D;
    }


    Color MipLevel::get_texel(int tx, int ty) {
        return Color(&texels[tx * 3 + ty * width * 3]);
    }

    Color Texture::sample_nearest(Vector2D uv, int level) {
        // TODO: Task 5: Fill this in.
        auto& mip = mipmap[level];
        // return magenta for invalid level
        if (level > mipmap.size() || level < 0) {
            return Color(1, 0, 1);
        }
        // round up/down to get the nearest texel
        int x = round(uv.x * (mipmap[level].width));
        int y = round(uv.y * (mipmap[level].height));
        return mipmap[level].get_texel(x, y);

    }

    Color Texture::sample_bilinear(Vector2D uv, int level) {
        // TODO: Task 5: Fill this in.
        auto& mip = mipmap[level];
        //// return magenta for invalid level
        if (level > mipmap.size() || level < 0) {
            return Color(1, 0, 1);
        }

        int x_l = floor(uv.x * (mipmap[level].width));
        int x_r = ceil(uv.x * (mipmap[level].width));
        int y_d = floor(uv.y * (mipmap[level].height));
        int y_u = ceil(uv.y * (mipmap[level].height));

        // get the four nearest neighbors
        Color u01 = mipmap[level].get_texel(x_l, y_u);
        Color u11 = mipmap[level].get_texel(x_r, y_u);
        Color u00 = mipmap[level].get_texel(x_l, y_d);
        Color u10 = mipmap[level].get_texel(x_r, y_d);

        // implement lerp(v, t0, t1) = t0 + v(t1 - t0) function
        float s = uv.x * (mipmap[level].width) - x_l;

        Color u0 = u00 + s * (u10 + (-1) * u00);
        Color u1 = u01 + s * (u11 + (-1) * u01);

        float t = uv.y * (mipmap[level].height) - y_d;
        Color f = u0 + t * (u1 + (-1) * u0);
        return f;

    }



    /****************************************************************************/

    // Helpers

    inline void uint8_to_float(float dst[3], unsigned char* src) {
        uint8_t* src_uint8 = (uint8_t*)src;
        dst[0] = src_uint8[0] / 255.f;
        dst[1] = src_uint8[1] / 255.f;
        dst[2] = src_uint8[2] / 255.f;
    }

    inline void float_to_uint8(unsigned char* dst, float src[3]) {
        uint8_t* dst_uint8 = (uint8_t*)dst;
        dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
        dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
        dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
    }

    void Texture::generate_mips(int startLevel) {

        // make sure there's a valid texture
        if (startLevel >= mipmap.size()) {
            std::cerr << "Invalid start level";
        }

        // allocate sublevels
        int baseWidth = mipmap[startLevel].width;
        int baseHeight = mipmap[startLevel].height;
        int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

        numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
        mipmap.resize(startLevel + numSubLevels + 1);

        int width = baseWidth;
        int height = baseHeight;
        for (int i = 1; i <= numSubLevels; i++) {

            MipLevel& level = mipmap[startLevel + i];

            // handle odd size texture by rounding down
            width = max(1, width / 2);
            //assert (width > 0);
            height = max(1, height / 2);
            //assert (height > 0);

            level.width = width;
            level.height = height;
            level.texels = vector<unsigned char>(3 * width * height);
        }

        // create mips
        int subLevels = numSubLevels - (startLevel + 1);
        for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
            mipLevel++) {

            MipLevel& prevLevel = mipmap[mipLevel - 1];
            MipLevel& currLevel = mipmap[mipLevel];

            int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
            int currLevelPitch = currLevel.width * 3; // 32 bit RGB

            unsigned char* prevLevelMem;
            unsigned char* currLevelMem;

            currLevelMem = (unsigned char*)&currLevel.texels[0];
            prevLevelMem = (unsigned char*)&prevLevel.texels[0];

            float wDecimal, wNorm, wWeight[3];
            int wSupport;
            float hDecimal, hNorm, hWeight[3];
            int hSupport;

            float result[3];
            float input[3];

            // conditional differentiates no rounding case from round down case
            if (prevLevel.width & 1) {
                wSupport = 3;
                wDecimal = 1.0f / (float)currLevel.width;
            }
            else {
                wSupport = 2;
                wDecimal = 0.0f;
            }

            // conditional differentiates no rounding case from round down case
            if (prevLevel.height & 1) {
                hSupport = 3;
                hDecimal = 1.0f / (float)currLevel.height;
            }
            else {
                hSupport = 2;
                hDecimal = 0.0f;
            }

            wNorm = 1.0f / (2.0f + wDecimal);
            hNorm = 1.0f / (2.0f + hDecimal);

            // case 1: reduction only in horizontal size (vertical size is 1)
            if (currLevel.height == prevLevel.height) {
                //assert (currLevel.height == 1);

                for (int i = 0; i < currLevel.width; i++) {
                    wWeight[0] = wNorm * (1.0f - wDecimal * i);
                    wWeight[1] = wNorm * 1.0f;
                    wWeight[2] = wNorm * wDecimal * (i + 1);

                    result[0] = result[1] = result[2] = 0.0f;

                    for (int ii = 0; ii < wSupport; ii++) {
                        uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
                        result[0] += wWeight[ii] * input[0];
                        result[1] += wWeight[ii] * input[1];
                        result[2] += wWeight[ii] * input[2];
                    }

                    // convert back to format of the texture
                    float_to_uint8(currLevelMem + (3 * i), result);
                }

                // case 2: reduction only in vertical size (horizontal size is 1)
            }
            else if (currLevel.width == prevLevel.width) {
                //assert (currLevel.width == 1);

                for (int j = 0; j < currLevel.height; j++) {
                    hWeight[0] = hNorm * (1.0f - hDecimal * j);
                    hWeight[1] = hNorm;
                    hWeight[2] = hNorm * hDecimal * (j + 1);

                    result[0] = result[1] = result[2] = 0.0f;
                    for (int jj = 0; jj < hSupport; jj++) {
                        uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
                        result[0] += hWeight[jj] * input[0];
                        result[1] += hWeight[jj] * input[1];
                        result[2] += hWeight[jj] * input[2];
                    }

                    // convert back to format of the texture
                    float_to_uint8(currLevelMem + (currLevelPitch * j), result);
                }

                // case 3: reduction in both horizontal and vertical size
            }
            else {

                for (int j = 0; j < currLevel.height; j++) {
                    hWeight[0] = hNorm * (1.0f - hDecimal * j);
                    hWeight[1] = hNorm;
                    hWeight[2] = hNorm * hDecimal * (j + 1);

                    for (int i = 0; i < currLevel.width; i++) {
                        wWeight[0] = wNorm * (1.0f - wDecimal * i);
                        wWeight[1] = wNorm * 1.0f;
                        wWeight[2] = wNorm * wDecimal * (i + 1);

                        result[0] = result[1] = result[2] = 0.0f;

                        // convolve source image with a trapezoidal filter.
                        // in the case of no rounding this is just a box filter of width 2.
                        // in the general case, the support region is 3x3.
                        for (int jj = 0; jj < hSupport; jj++)
                            for (int ii = 0; ii < wSupport; ii++) {
                                float weight = hWeight[jj] * wWeight[ii];
                                uint8_to_float(input, prevLevelMem +
                                    prevLevelPitch * (2 * j + jj) +
                                    3 * (2 * i + ii));
                                result[0] += weight * input[0];
                                result[1] += weight * input[1];
                                result[2] += weight * input[2];
                            }

                        // convert back to format of the texture
                        float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
                    }
                }
            }
        }
    }

}
