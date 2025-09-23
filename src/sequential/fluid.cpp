#include "fluid.hpp"
#include "cuda/particle.hpp"
#include "time.hpp"
#include <algorithm>
#include <functional>
#include <vector>
Fluid::Fluid(float cell_size, int width, int height) : m_cell_size(cell_size), m_width(width), m_height(height), m_num_cells(width * height){
    cell_color = std::vector<Color>(m_num_cells);
    cell_type = std::vector<eCellTypes>(m_num_cells, eCellTypes::Air);
    solid = std::vector<float>(m_num_cells, 1.f);
    particle_density = std::vector<float>(m_num_cells, 0.f);
    pressure = std::vector<float>(m_num_cells, 0.f);
    smoke = std::vector<float>(m_num_cells, 0.f);

    velocities = std::vector<vec2f>(m_num_cells, vec2f(0, 0));
    air_velocities = std::vector<vec2f>(m_num_cells, vec2f(0, 0));
    prev_velocities = std::vector<vec2f>(m_num_cells, vec2f(0, 0));
    velocities_diff = std::vector<vec2f>(m_num_cells, vec2f(0, 0));

    for(int i = 0; i < m_height; i++) {
        for(int j = 0; j < m_width; j++) {
            if(j == 0 || i == 0 || i == m_height - 1 || j == m_width - 1)
                solid[i * width + j] = 0.f;
            if(i > 10 && i < 20 && j > 10 && j < 20) {
                smoke[i * width + j] = 1.f;
            }
        }
    }
}

std::tuple<float, float, float> getCoords(float coord, float inv_cell, float limit) {
    auto v0 = fmin(floorf(coord*inv_cell), limit);
    auto share_v = (coord - v0/inv_cell) * inv_cell;
    auto v1 = fmin(v0 + 1, limit);
    return {share_v, v0, v1};
}
std::tuple<float, float, float, float> combineShares(float share0_x, float share0_y) {
    auto share1_x = 1.0 - share0_x;
    auto share1_y = 1.0 - share0_y;

    auto d0 = share1_x*share1_y;
    auto d1 = share0_x*share1_y;
    auto d2 = share0_x*share0_y;
    auto d3 = share1_x*share0_y;
    return {d0, d1, d2, d3};
}

void Fluid::updateParticleDensity(Particles& particles)
{
    int n = m_width;
    float h = m_cell_size;
    float inv_csize = 1.f / m_cell_size;
    float half_csize = 0.5 * h;

    std::fill(particle_density.begin(), particle_density.end(), 0.f);

    for (int i = 0; i < Particles::max_particle_count; i++) {
        auto x = particles.position[i].x;
        auto y = particles.position[i].y;

        auto [share0_x, x0, x1] = getCoords(x - half_csize, inv_csize, m_width - 1);
        auto [share0_y, y0, y1] = getCoords(y - half_csize, inv_csize, m_height - 1);
        auto [d0, d1, d2, d3] = combineShares(share0_x, share0_y);

        if(y1 >= m_height) y1 = y0;
        if(x1 >= m_width) x1 = x0;
        if (x0 < m_width && y0 < m_height) particle_density[x0 + y0 * n] += d0;
        if (x1 < m_width && y0 < m_height) particle_density[x1 + y0 * n] += d1;
        if (x1 < m_width && y1 < m_height) particle_density[x1 + y1 * n] += d2;
        if (x0 < m_width && y1 < m_height) particle_density[x0 + y1 * n] += d3;
    }

    if (particleRestDensity == 0.0) {
        auto sum = 0.0;
        auto numFluidCells = 0;

        for (int i = 0; i < m_num_cells; i++) {
            if (cell_type[i] == eCellTypes::Fluid) {
                sum += particle_density[i];
                numFluidCells++;
            }
        }

        if (numFluidCells > 0)
            particleRestDensity = sum / numFluidCells;
    }
}
void Fluid::transferVelocitiesToGrid(float flipRatio, Particles& particles) {
    auto inv_cell_size = 1.f / m_cell_size;
    auto half_cell_size = 0.5 * m_cell_size;

    prev_velocities = velocities;

    for (auto i = 0; i < m_num_cells; i++) 
        cell_type[i] = (solid[i] == 0.0 ? eCellTypes::Solid : eCellTypes::Air);

    for (auto i = 0; i < Particles::max_particle_count; i++) {
        auto x = particles.position[i].x;
        auto y = particles.position[i].y;
        auto xi = /* std::clamp( */floorf(x * inv_cell_size)/* , 0.f, m_width - 1.f) */;
        auto yi = /* std::clamp( */floorf(y * inv_cell_size)/* , 0.f, m_height - 1.f) */;
        auto cellNr = xi + yi * m_width;
        if (cell_type[cellNr] == eCellTypes::Air) {
            cell_type[cellNr] = eCellTypes::Fluid;
        }
    }

    std::fill(velocities.begin(), velocities.end(), vec2f(0, 0));
    std::fill(velocities_diff.begin(), velocities_diff.end(), vec2f(0, 0));
    auto vel = [&](int idx, int component) -> float&{
        if(component == 0) {
            return velocities[idx].x;
        }
        return velocities[idx].y;
    };
    auto prev_vel = [&](int idx, int component) -> float&{
        if(component == 0) {
            return prev_velocities[idx].x;
        }
        return prev_velocities[idx].y;
    };
    auto diff= [&](int idx, int component) -> float&{
        if(component == 0) {
            return velocities_diff[idx].x;
        }
        return velocities_diff[idx].y;
    };

    for (auto component = 0; component < 2; component++) {
        auto offset_x = component == 0 ? 0.0 : half_cell_size;
        auto offset_y = component == 0 ? half_cell_size : 0.0;

        for (auto i = 0; i < Particles::max_particle_count; i++) {
            auto x = particles.position[i].x;
            auto y = particles.position[i].y;

            auto [share0_x, x0, x1] = getCoords(x - offset_x, inv_cell_size, m_width - 1);
            auto [share0_y, y0, y1] = getCoords(y - offset_y, inv_cell_size, m_height - 1);

            auto [d0, d1, d2, d3] = combineShares(share0_x, share0_y);

            auto bl = x0 + y0 * m_width;
            auto br = x1 + y0 * m_width;
            auto tr = x1 + y1 * m_width;
            auto tl = x0 + y1 * m_width;

            //transfering weighted velocity
            auto pv = particles.velocity[i].x;
            if(component == 1)
                pv = particles.velocity[i].y;

            diff(bl, component) += d0;
            diff(br, component) += d1;
            diff(tr, component) += d2;
            diff(tl, component) += d3;

            vel(bl, component) += pv * d0;
            vel(br, component) += pv * d1;
            vel(tr, component) += pv * d2;
            vel(tl, component) += pv * d3;
        }
        for (auto i = 0; i < velocities.size(); i++) {
            if (diff(i, component) > 0.0) {
                vel(i, component) /= diff(i, component);
            }
        }
        for (auto j = 0; j < m_height; j++) {
            for (auto i = 0; i < m_width; i++) {
                auto solid = (cell_type[i + j * m_width] == eCellTypes::Solid);
                if(component == 0) {
                    if (solid || (i > 0 && cell_type[(i - 1) + j * m_width] == eCellTypes::Solid))
                        velocities[i + j * m_width].x = 0;
                }else {
                    if (solid || (j > 0 && cell_type[i + (j - 1) * m_width] == eCellTypes::Solid))
                        velocities[i + j * m_width].y = 0;
                }
            }
        }
    }
}
void Fluid::transferVelocitiesFromGrid(float flipRatio, Particles& particles) {
    auto inv_cell_size = 1.f / m_cell_size;
    auto half_cell_size = 0.5 * m_cell_size;

    for (auto component = 0; component < 2; component++) {
        auto offset_x = component == 0 ? 0.0 : half_cell_size;
        auto offset_y = component == 0 ? half_cell_size : 0.0;

        auto vel = [&](int idx) -> float&{
            if(component == 0) {
                return velocities[idx].x;
            }
            return velocities[idx].y;
        };
        auto prev_vel = [&](int idx) -> float&{
            if(component == 0) {
                return prev_velocities[idx].x;
            }
            return prev_velocities[idx].y;
        };
        auto diff= [&](int idx) -> float&{
            if(component == 0) {
                return velocities_diff[idx].x;
            }
            return velocities_diff[idx].y;
        };

        for (auto i = 0; i < Particles::max_particle_count; i++) {
            auto x = particles.position[i].x;
            auto y = particles.position[i].y;

            auto [share0_x, x0, x1] = getCoords(x - offset_x, inv_cell_size, m_width - 1);
            auto [share0_y, y0, y1] = getCoords(y - offset_y, inv_cell_size, m_height - 1);

            auto [d0, d1, d2, d3] = combineShares(share0_x, share0_y);

            auto bl = x0 + y0 * m_width;
            auto br = x1 + y0 * m_width;
            auto tr = x1 + y1 * m_width;
            auto tl = x0 + y1 * m_width;

            auto offset = (component == 0 ? 1 : m_width);
            auto valid0 = cell_type[bl] == eCellTypes::Fluid || cell_type[bl - offset] == eCellTypes::Fluid ? 1.0 : 0.0;
            auto valid1 = cell_type[br] == eCellTypes::Fluid || cell_type[br - offset] == eCellTypes::Fluid ? 1.0 : 0.0;
            auto valid2 = cell_type[tr] == eCellTypes::Fluid || cell_type[tr - offset] == eCellTypes::Fluid ? 1.0 : 0.0;
            auto valid3 = cell_type[tl] == eCellTypes::Fluid || cell_type[tl - offset] == eCellTypes::Fluid ? 1.0 : 0.0;

            auto v = particles.velocity[i].x;
            if(component == 1) 
                v = particles.velocity[i].y;
            auto d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

            if (d > 0.0) {

                auto picV = (valid0 * d0 * vel(bl) + valid1 * d1 * vel(br) + valid2 * d2 * vel(tr) + valid3 * d3 * vel(tl)) / d;
                auto corr = (valid0 * d0 * (vel(bl) - prev_vel(bl)) + valid1 * d1 * (vel(br) - prev_vel(br))
                    + valid2 * d2 * (vel(tr) - prev_vel(tr)) + valid3 * d3 * (vel(tl) - prev_vel(tl))) / d;
                auto flipV = v + corr;
                auto vel_comp = (1.0 - flipRatio) * picV + flipRatio * flipV;
                if(component == 0) {
                    particles.velocity[i].x = vel_comp;
                }else {
                    particles.velocity[i].y = vel_comp;
                }
            }
        }
    }
}
void Fluid::solveIncompressibility(float dt, eCellTypes expected_type, std::vector<vec2f>& vels, std::function<bool(int)> solid, float density, int numIters, float overRelaxation, bool compensateDrift) {
    auto n = m_width;
    auto cp = density * m_cell_size / dt;

    for (auto iter = 0; iter < numIters; iter++) {

        for (auto j = 0; j < m_height; j++) {
            for (auto i = 0; i < m_width; i++) {

                if (cell_type[i + j*n] != expected_type)
                    continue;

                auto center = i + j*n;
                auto left = (i - 1) + j*n;
                auto right = (i + 1) + j*n;
                auto bottom = i + (j - 1)*n;
                auto top = i + (j + 1)*n;

                auto sx0 = !solid(left);
                auto sx1 = !solid(right);
                auto sy0 = !solid(bottom);
                auto sy1 = !solid(top);
                auto s = sx0 + sx1 + sy0 + sy1;
                if (s == 0.0)
                    continue;

                auto div = vels[right].x - vels[center].x + 
                    vels[top].y - vels[center].y;

                if (particleRestDensity > 0.0 && compensateDrift) {
                    auto k = 0.6;
                    auto compression = particle_density[i + j*n] - particleRestDensity;
                    if (compression < 0.0)
                        div = div - k * compression;
                }

                auto dp = -div / s;
                dp *= overRelaxation;

                vels[center].x -= sx0 * dp;
                vels[right].x += sx1 * dp;
                vels[center].y -= sy0 * dp;
                vels[top].y += sy1 * dp;
            }
        }
    }
}
std::map<std::string, float> Fluid::simulate(Particles& particles, AABB sim_area, float dt, vec2f gravity, int numPressureIters, int numParticleIters, float overRelaxation, bool compensateDrift) {
    auto numSubSteps = 1;
    auto sdt = dt / numSubSteps;
    std::map<std::string, float> bench;

    sim_area.setSize(sim_area.size() - vec2f(m_cell_size, m_cell_size)*2.f);

    float col_time = 0.f;
    float fluid_time = 0.f;
    for (int step = 0; step < numSubSteps; step++) {
        Stopwatch local_stop;
        {
            ParticleSolveBlock solv(particles);
            for(int i = 0; i < numParticleIters; i++) {
                accelerate(particles, gravity);
                bench["particles::accelerate"] += local_stop.restart();
                integrate(particles, sdt / (float)numParticleIters);
                bench["particles::integrate"] += local_stop.restart();
                constraint(particles, sim_area);
                bench["particles::constraint"] += local_stop.restart();
                auto results = collide(particles, sim_area);
                for(auto [key, v] : results) bench[key] += v;
                bench["particles::collide"] += local_stop.restart();
            }
        }
        bench["fluid::advectSmoke"] += local_stop.restart();
        transferVelocitiesToGrid(1.f, particles);
        bench["fluid::transfer1"] += local_stop.restart();
        updateParticleDensity(particles);
        bench["fluid::density"] += local_stop.restart();
        prev_velocities = velocities;
        solveIncompressibility(sdt, eCellTypes::Fluid, velocities, [&](int idx){ return this->cell_type[idx] == eCellTypes::Solid; }, fluid_density, numPressureIters, overRelaxation, compensateDrift);
        advectAny(sdt, air_velocities, air_velocities, [&](vec2f& v)->float&{return v.x;}, 0.f, m_cell_size/2.f);
        advectAny(sdt, air_velocities, air_velocities, [&](vec2f& v)->float&{return v.y;}, m_cell_size/2.f, 0.f);
        float a = 0.03;
        for(int i = 0; i < m_num_cells; i++) {
            air_velocities[i] = (velocities[i] * a) + (air_velocities[i] * (1.f - a));
        }
        solveIncompressibility(sdt, eCellTypes::Air, air_velocities, [&](int idx){ return this->cell_type[idx] == eCellTypes::Solid || this->cell_type[idx] == eCellTypes::Fluid;}, air_density, numPressureIters, overRelaxation, false);
        advectAny(sdt, smoke, air_velocities, [&](float& f)->float&{return f;}, m_cell_size/2.f, m_cell_size/2.f);
        bench["fluid::incompressibility"] += local_stop.restart();
        bench["fluid::advectSmoke"] += local_stop.restart();
        transferVelocitiesFromGrid(flipRatio, particles);
        bench["fluid::transfer2"] += local_stop.restart();
    }
    return bench;
}
void Fluid::draw(AABB area, Particles& particles, sf::RenderTarget &window,
          std::unordered_map<eCellTypes, Color> color_table) {
    std::vector<float> vals(m_num_cells, 0);
    updateParticleDensity(particles);
    for (auto i = 0; i < Particles::max_particle_count; i++) {
        auto x = particles.position[i].x;
        auto y = particles.position[i].y;

        auto [share0_x, x0, x1] = getCoords(x, 1.f / m_cell_size, m_width - 1);
        auto [share0_y, y0, y1] = getCoords(y, 1.f / m_cell_size, m_height - 1);

        auto [d0, d1, d2, d3] = combineShares(share0_x, share0_y);

        auto bl = x0 + y0 * m_width;
        auto br = x1 + y0 * m_width;
        auto tr = x1 + y1 * m_width;
        auto tl = x0 + y1 * m_width;

        vals[bl] += d0;
        vals[br] += d1;
        vals[tr] += d2;
        vals[tl] += d3;
    }
    float max_v = *std::max_element(vals.begin(), vals.end());

    sf::Image img(sf::Vector2u(m_width, m_height));
    for (int i = 0; i < m_height; i++) {
        for (int j = 0; j < m_width; j++) {
            auto type = cell_type[i * m_width + j];
            auto coord = sf::Vector2u(j, m_height - i - 1);
            if (!color_table.contains(type)) {
                img.setPixel(coord, sf::Color(255, 0, 255));
            } else if (type == eCellTypes::Fluid) {
                float val = (particle_density[j + i*m_width] / particleRestDensity);
                val = std::clamp(val, 0.1f, 1.f);
                uint8_t r = 105 * (val * 0.7f + 0.3f);
                uint8_t g = 120 * (val * 0.7f + 0.3f);
                uint8_t b = 220 * (val * 0.7f + 0.3f);
                img.setPixel(coord, sf::Color(r, g, b));
            }else if (smoke[i * m_width + j] != 0.f) {
                uint8_t val = 255 *smoke[i * m_width + j]; 
                img.setPixel(coord, sf::Color(val, val, val));
            }else {
                img.setPixel(coord, color_table.at(type));
            }
        }
    }
    sf::Texture tex;
    auto success = tex.loadFromImage(img);
    sf::Sprite spr(tex);
    vec2f scale = vec2f(area.size().x / m_width, area.size().y / m_height);
    spr.setScale(scale);
    spr.setPosition(area.bl());
    window.draw(spr);
}
