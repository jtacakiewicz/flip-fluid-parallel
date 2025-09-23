#ifndef EMP_FLUID_HPP
#define EMP_FLUID_HPP
#include "SFML/Graphics/Image.hpp"
#include "SFML/Graphics/Sprite.hpp"
#include "SFML/Graphics/Texture.hpp"
#include "SFML/Graphics/RenderTarget.hpp"
#include "particle.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <unordered_map>
enum class eCellTypes {
    Fluid,
    Air,
    Solid
};
class Fluid {
    int m_width;
    int m_height;
    int m_num_cells;
    float m_cell_size;
    std::vector<Color> cell_color;
    std::vector<eCellTypes> cell_type;
    std::vector<float> solid;
    std::vector<float> pressure;
    std::vector<vec2f> prev_velocities;
    std::vector<vec2f> velocities;

    std::vector<vec2f> prev_air_velocities;
    std::vector<vec2f> air_velocities;

    std::vector<vec2f> velocities_diff;
    std::vector<float> particle_density;
    float particleRestDensity = 0;
    //x and y are not in grid coordinates but in global
    template<class T, class extracter>
    float sampleField(float x, float y, T* field, extracter getF, float dx_offset = NAN, float dy_offset = NAN) const;

    template<class T, class extrT>
    void advectAny(float dt, std::vector<T>& vec, std::vector<vec2f>& vels, extrT func, float dx, float dy) const;

public:
    inline int width() const {
        return m_width;
    }
    inline int height() const {
        return m_height;
    }
    inline float cell_size() const {
        return m_cell_size;
    }
    std::vector<float> smoke;
    float fluid_density = 1;
    float air_density = 0.001;
    float flipRatio = 0.9f;

    void updateParticleDensity(Particles& particles);
    void transferVelocitiesToGrid(float flipRatio, Particles& particles);
    void transferVelocitiesFromGrid(float flipRatio, Particles& particles);
    void solveIncompressibility(float dt, eCellTypes expected_type, std::vector<vec2f>& vels, std::function<bool(int)> solid, float density, int numIters, float overRelaxation, bool compensateDrift);

    std::map<std::string, float> simulate(Particles& particles, AABB sim_area, float dt, vec2f gravity, int numPressureIters, int numParticleIters, float overRelaxation, bool compensateDrift);
    void draw(AABB area, Particles& particles, sf::RenderTarget &window,
              std::unordered_map<eCellTypes, Color> color_table = {
              {eCellTypes::Air, Color(0, 0, 50)},
              {eCellTypes::Solid, Color(100, 100, 100)},
              {eCellTypes::Fluid, Color(70, 100, 220)}});
    Fluid(float cell_size, int width, int height);
};
template<class T, class extracter>
float Fluid::sampleField(float x, float y, T* field, extracter getF, float dx_offset, float dy_offset) const {
    float inv_csize = 1.0 / m_cell_size;
    float half_csize = 0.5 * m_cell_size;

    x = std::clamp(x, m_cell_size, (m_width - 1) * m_cell_size);
    y = std::clamp(y, m_cell_size, (m_height - 1) * m_cell_size);

    float dx = (isnan(dx_offset) ? half_csize : dx_offset);
    float dy = (isnan(dy_offset) ? half_csize : dy_offset);

    float x0 = std::fminf(std::floor((x-dx)*inv_csize), (float)m_width-1.f);
    float tx = (x - dx - x0*m_cell_size) * inv_csize;
    float x1 = std::fminf(x0 + 1, m_width-1.f);

    float y0 = std::fminf(std::floor((y-dy)*inv_csize), (float)m_height-1.f);
    float ty = (y - dy - y0*m_cell_size) * inv_csize;
    float y1 = std::fminf(y0 + 1, m_height-1.f);

    float sx = 1.0 - tx;
    float sy = 1.0 - ty;

    auto n = m_width;
    float val = sx*sy * getF(field[int(x0 + y0*n)]) +
        tx*sy * getF(field[int(x1 + y0*n)]) +
        tx*ty * getF(field[int(x1 + y1*n)]) +
        sx*ty * getF(field[int(x0 + y1*n)]);

    return val;
}

template<class T, class extrT>
void Fluid::advectAny(float dt, std::vector<T>& vec, std::vector<vec2f>& vels, extrT func, float dx, float dy) const {
    auto temporary = vec;
    auto half_size = 0.5 * m_cell_size;
    auto n = m_width;

    for (auto j = 0; j < m_height-1; j++) {
        for (auto i = 0; i < m_width-1; i++) {
            if (solid[i + j * m_width] != 0.0) {
                auto uu = (vels[i + j * n].x + vels[i+1 + j * n].x)*0.5;
                auto vv = (vels[i + j * n].y + vels[i + (j+1) * n].y)*0.5;
                auto x = i*m_cell_size + half_size - dt*uu;
                auto y = j*m_cell_size + half_size - dt*vv;

                func(temporary[i + j * n]) = sampleField<T, extrT>(x,y, vec.data(), func, dx, dy);
            }
        }	 
    }
    std::swap(vec, temporary);
}
#endif
