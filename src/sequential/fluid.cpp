#include "fluid.hpp"
#include "cuda/particle.hpp"
#include "time.hpp"
#include <algorithm>
#include <vector>
Fluid::Fluid(float cell_size, int width, int height) : m_cell_size(cell_size), m_width(width), m_height(height), m_num_cells(width * height){
    cell_color = std::vector<Color>(m_num_cells);
    cell_type = std::vector<eCellTypes>(m_num_cells, eCellTypes::Air);
    solid = std::vector<float>(m_num_cells, 1.f);
    particle_density = std::vector<float>(m_num_cells, 0.f);
    pressure = std::vector<float>(m_num_cells, 0.f);
    smoke = std::vector<float>(m_num_cells, 0.f);

    velocities = std::vector<vec2f>(m_num_cells, vec2f(0, 0));
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

        auto x0 = floorf((x - half_csize) * inv_csize);
        auto tx = ((x - half_csize) - x0 * h) * inv_csize;
        auto x1 = fmin(x0 + 1, m_width-2);
        
        auto y0 = floorf((y-half_csize)*inv_csize);
        auto ty = ((y - half_csize) - y0*h) * inv_csize;
        auto y1 = fmin(y0 + 1, m_height-2);

        auto sx = 1.0 - tx;
        auto sy = 1.0 - ty;

        if (x0 < m_width && y0 < m_height) particle_density[x0 + y0 * n] += sx * sy;
        if (x1 < m_width && y0 < m_height) particle_density[x1 + y0 * n] += tx * sy;
        if (x1 < m_width && y1 < m_height) particle_density[x1 + y1 * n] += tx * ty;
        if (x0 < m_width && y1 < m_height) particle_density[x0 + y1 * n] += sx * ty;
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
std::tuple<float, float, float, float> combineShares(float share0_x, float share0_y) {
    auto share1_x = 1.0 - share0_x;
    auto share1_y = 1.0 - share0_y;

    auto d0 = share1_x*share1_y;
    auto d1 = share0_x*share1_y;
    auto d2 = share0_x*share0_y;
    auto d3 = share1_x*share0_y;
    return {d0, d1, d2, d3};
}
std::tuple<float, float, float> getCoords(float coord, float inv_cell, float limit) {
    auto v0 = fmin(floorf(coord*inv_cell), limit);
    auto share_v = (coord - v0/inv_cell) * inv_cell;
    auto v1 = fmin(v0 + 1, limit);
    return {share_v, v0, v1};
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

    std::fill(velocities_diff.begin(), velocities_diff.end(), vec2f(0, 0));
    std::fill(velocities.begin(), velocities.end(), vec2f(0, 0));
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

            // x = std::clamp(x, h, (m_width - 1) * h);
            // y = std::clamp(y, h, (m_height - 1) * h);

            auto [share0_x, x0, x1] = getCoords(x - offset_x, inv_cell_size, m_width - 2);
            auto [share0_y, y0, y1] = getCoords(y - offset_y, inv_cell_size, m_height - 2);

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
                        velocities[i + j * m_width].x = prev_velocities[i + j * m_width].x;
                }else {
                    if (solid || (j > 0 && cell_type[i + (j - 1) * m_width] == eCellTypes::Solid))
                        velocities[i + j * m_width].y = prev_velocities[i + j * m_width].y;
                }
            }
        }
    }
    for(int i = 0; i < m_num_cells; i++) {
        for(int c = 0; c < 2; c++) {
            if(diff(i, c) <= 0.f)
                vel(i, c) = prev_vel(i, c);
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

            auto [share0_x, x0, x1] = getCoords(x - offset_x, inv_cell_size, m_width - 2);
            auto [share0_y, y0, y1] = getCoords(y - offset_y, inv_cell_size, m_height - 2);

            auto [d0, d1, d2, d3] = combineShares(share0_x, share0_y);

            auto bl = x0 + y0 * m_width;
            auto br = x1 + y0 * m_width;
            auto tr = x1 + y1 * m_width;
            auto tl = x0 + y1 * m_width;

            auto offset = (component == 0 ? 1 : m_width);
            auto valid0 = cell_type[bl] != eCellTypes::Air || cell_type[bl - offset] != eCellTypes::Air ? 1.0 : 0.0;
            auto valid1 = cell_type[br] != eCellTypes::Air || cell_type[br - offset] != eCellTypes::Air ? 1.0 : 0.0;
            auto valid2 = cell_type[tr] != eCellTypes::Air || cell_type[tr - offset] != eCellTypes::Air ? 1.0 : 0.0;
            auto valid3 = cell_type[tl] != eCellTypes::Air || cell_type[tl - offset] != eCellTypes::Air ? 1.0 : 0.0;

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
                float change = corr;
                if(component == 0) {
                    particles.velocity[i].x = vel_comp;
                }else {
                    particles.velocity[i].y = vel_comp;
                }
            }
        }
    }
}
void Fluid::solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift) {
    std::fill(pressure.begin(), pressure.end(), 0.f);
    prev_velocities = velocities;

    auto n = m_width;

    for (auto iter = 0; iter < numIters; iter++) {

        for (auto j = 1; j < m_height-1; j++) {
            for (auto i = 1; i < m_width-1; i++) {

                if (cell_type[i + j*n] == eCellTypes::Solid)
                    continue;

                auto center = i + j*n;
                auto left = (i - 1) + j*n;
                auto right = (i + 1) + j*n;
                auto bottom = i + (j - 1)*n;
                auto top = i + (j + 1)*n;
                float rho_center = (cell_type[center] == eCellTypes::Fluid) ? fluid_density : air_density;
                float rho_left   = (cell_type[left]   == eCellTypes::Fluid) ? fluid_density : air_density;
                float rho_right  = (cell_type[right]  == eCellTypes::Fluid) ? fluid_density : air_density;
                float rho_bottom = (cell_type[bottom] == eCellTypes::Fluid) ? fluid_density : air_density;
                float rho_top    = (cell_type[top]    == eCellTypes::Fluid) ? fluid_density : air_density;

                // Inverse densities on faces
                float inv_rho_x0 = (solid[left] == 1.0f)
                    ? 0.5f * (1.0f / rho_center + 1.0f / rho_left)
                    : 0.0f;
                float inv_rho_x1 = (solid[right] == 1.0f)
                    ? 0.5f * (1.0f / rho_center + 1.0f / rho_right)
                    : 0.0f;

                float inv_rho_y0 = (solid[bottom] == 1.0f)
                    ? 0.5f * (1.0f / rho_center + 1.0f / rho_bottom)
                    : 0.0f;
                float inv_rho_y1 = (solid[top] == 1.0f)
                    ? 0.5f * (1.0f / rho_center + 1.0f / rho_top)
                    : 0.0f;

                auto div = velocities[right].x - velocities[center].x + 
                    velocities[top].y - velocities[center].y;

                div /= m_cell_size;

                if (particleRestDensity > 0.0f && compensateDrift) {
                    float compression = particle_density[center] - particleRestDensity;
                    if (compression > 0.0f) {
                        float k = 1.0f;
                        div -= k * compression;
                    }
                }

                // Build diagonal coefficient
                float A_center = (inv_rho_x0 + inv_rho_x1 +
                    inv_rho_y0 + inv_rho_y1);

                if (A_center == 0.0f)
                    continue;

                float dp = -div / A_center;

                // Over-relaxation
                dp *= overRelaxation;

                pressure[center] += dp;

                // Correct velocities
                velocities[center].x -= dp * inv_rho_x0;
                velocities[right].x += dp * inv_rho_x1;
                velocities[center].y -= dp * inv_rho_y0;
                velocities[top].y += dp * inv_rho_y1;
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
        advectAny(sdt, velocities, [&](vec2f& v)->float&{return v.x;}, 0.f, m_cell_size/2.f);
        advectAny(sdt, velocities, [&](vec2f& v)->float&{return v.y;}, m_cell_size/2.f, 0.f);
        bench["fluid::advectSmoke"] += local_stop.restart();
        transferVelocitiesToGrid(1.f, particles);
        bench["fluid::transfer1"] += local_stop.restart();
        updateParticleDensity(particles);
        bench["fluid::density"] += local_stop.restart();
        solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
        bench["fluid::incompressibility"] += local_stop.restart();
        advectAny(sdt, smoke, [&](float& f)->float&{return f;});
        bench["fluid::advectSmoke"] += local_stop.restart();
        transferVelocitiesFromGrid(flipRatio, particles);
        bench["fluid::transfer2"] += local_stop.restart();
    }
    return bench;
}
void Fluid::draw(AABB area, Particles& particles, sf::RenderTarget &window,
          std::unordered_map<eCellTypes, Color> color_table) {
    std::vector<float> vals(m_num_cells, 0);
    for (auto i = 0; i < Particles::max_particle_count; i++) {
        auto x = particles.position[i].x;
        auto y = particles.position[i].y;

        auto [share0_x, x0, x1] = getCoords(x, 1.f / m_cell_size, m_width - 2);
        auto [share0_y, y0, y1] = getCoords(y, 1.f / m_cell_size, m_height - 2);

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
                float val = (vals[i * m_width + j] / max_v);
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
