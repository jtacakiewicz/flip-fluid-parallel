#include "fluid.hpp"
#include "benchmark/benchmark.hpp"
#include "cuda/particle.hpp"
#include "time.hpp"
#include <SFML/System/Vector2.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
Fluid::Fluid(float cell_size, int width, int height)
    : m_cell_size(cell_size)
    , m_width(width)
    , m_height(height)
    , m_num_cells(width * height)
{
    cell_color = std::vector<Color>(m_num_cells);
    cell_type = std::vector<eCellTypes>(m_num_cells, eCellTypes::Air);
    solid = std::vector<float>(m_num_cells, 0.f);
    particle_density = std::vector<float>(m_num_cells, 0.f);
    pressure = std::vector<float>(m_num_cells, 0.f);
    smoke = std::vector<float>(m_num_cells, 0.f);

    velocities = std::vector<vec2f>(m_num_cells, vec2f(0, 0));
    air_velocities = std::vector<vec2f>(m_num_cells, vec2f(0, 0));
    prev_velocities = std::vector<vec2f>(m_num_cells, vec2f(0, 0));
    velocities_diff = std::vector<vec2f>(m_num_cells, vec2f(0, 0));

    for(int i = 0; i < m_height; i++) {
        for(int j = 0; j < m_width; j++) {
            if(j == 0 || i == 0 || i == m_height - 1 || j == m_width - 1) {
                solid[i * width + j] = 1.f;
            }
            if(i > 10 && i < 20 && j > 10 && j < 20) {
                smoke[i * width + j] = 1.f;
            }
        }
    }
}

std::tuple<float, float, float> getCoords(float coord, float inv_cell, float limit)
{
    auto v0 = fmin(floorf(coord * inv_cell), limit);
    auto share_v = (coord - v0 / inv_cell) * inv_cell;
    auto v1 = fmin(v0 + 1, limit);
    return { share_v, v0, v1 };
}
std::tuple<float, float, float, float> combineShares(float share0_x, float share0_y)
{
    auto share1_x = 1.0 - share0_x;
    auto share1_y = 1.0 - share0_y;

    auto d0 = share1_x * share1_y;
    auto d1 = share0_x * share1_y;
    auto d2 = share0_x * share0_y;
    auto d3 = share1_x * share0_y;
    return { d0, d1, d2, d3 };
}

void Fluid::collideWithGrid(Particles &particles, float dt)
{
    EMP_BENCHMARK_FUNC()
    const std::vector<eCellTypes> collision_blocks = { eCellTypes::Solid };
    for(int i = 0; i < particles.max_particle_count; i++) {
        auto pos = particles.position[i];
        auto vel = particles.velocity[i];
        auto collision = findCollision(pos, vel * dt, collision_blocks);
        int ii = pos.x / m_cell_size;
        int jj = pos.y / m_cell_size;

        if(collision.has_value()) {
            auto n = collision->normal;
            particles.position[i] = collision->collision_point + n;
            particles.velocity[i] -= n * dot(n, vel);
        }
    }
}
std::optional<Fluid::CollisionDetection> Fluid::findCollision(vec2f origin, vec2f dir, const std::vector<eCellTypes> &col_types)
{
    auto origin_type = cell_type[floorf(origin.x / m_cell_size) + floorf(origin.y / m_cell_size) * m_width];
    auto backingUp = std::find(col_types.begin(), col_types.end(), origin_type) != col_types.end();
    if(backingUp) {
        dir *= -100.f;
    }
    if(qlen(dir) < 0.0001f) {
        return {};
    }
    float inv_c_size = 1.f / m_cell_size;
    float c_size = m_cell_size;
    int n = m_width;

    vec2f end = origin + dir;
    float slope = dir.y / dir.x;
    if(dir.x == 0.f) {
        slope = INFINITY;
    }
    float constant = origin.y - slope * origin.x;
    if(dir.x == 0.f) {
        constant = 0.f;
    }

    float x_adv_dir = std::copysign(1.f, dir.x);
    float y_adv_dir = std::copysign(1.f, dir.y);

    float cur_x = floorf(origin.x * inv_c_size) * c_size + c_size * (dir.x > 0.f);
    float cur_y = floorf(origin.y * inv_c_size) * c_size + c_size * (dir.y > 0.f);
    vec2f last_meas = origin;
    float summaric_len = 0.f;
    while(cur_x * x_adv_dir < end.x * x_adv_dir + c_size && cur_y * y_adv_dir < end.y * y_adv_dir + c_size) {
        float x_step = (cur_y - constant) / slope;
        if(slope == INFINITY) {
            x_step = origin.x;
        }
        float y_step = slope * cur_x + constant;
        vec2f y_axis(x_step, cur_y);
        vec2f x_axis(cur_x, y_step);
        vec2f norm;
        float len1 = qlen(last_meas - y_axis);
        float len2 = qlen(last_meas - x_axis);
        int i, j;
        if(len1 < len2) {
            last_meas = y_axis;
            norm = { 0.f, -y_adv_dir };
            i = x_step * inv_c_size;
            j = (cur_y + y_adv_dir * c_size * 0.5f) * inv_c_size;

            summaric_len += len1;
            cur_y += c_size * y_adv_dir;
        } else {
            last_meas = x_axis;
            norm = { -x_adv_dir, 0.f };
            i = (cur_x + x_adv_dir * c_size * 0.5f) * inv_c_size;
            j = y_step * inv_c_size;

            summaric_len += len2;
            cur_x += c_size * x_adv_dir;
        }
        if(qlen(last_meas - origin) > qlen(dir)) {
            break;
        }
        if(i >= m_width || j >= m_height || i < 0 || j < 0) {
            continue;
        }
        auto type = cell_type[i + j * n];
        auto hasCollided = std::find(col_types.begin(), col_types.end(), type) != col_types.end();
        if(hasCollided && !backingUp) {
            CollisionDetection result {
                .collision_point = last_meas,
                .normal = norm,
                .grid_idx = { i, j }
            };
            return result;
        } else if(!hasCollided && backingUp) {
            CollisionDetection result {
                .collision_point = last_meas,
                .normal = -norm,
                .grid_idx = { i, j }
            };
            return result;
        }
    }
    return {};
}

void Fluid::updateParticleDensity(Particles &particles)
{
    EMP_BENCHMARK_FUNC()
    int n = m_width;
    float h = m_cell_size;
    float inv_csize = 1.f / m_cell_size;
    float half_csize = 0.5 * h;

    std::fill(particle_density.begin(), particle_density.end(), 0.f);

    for(int i = 0; i < Particles::max_particle_count; i++) {
        auto x = particles.position[i].x;
        auto y = particles.position[i].y;

        auto [share0_x, x0, x1] = getCoords(x - half_csize, inv_csize, m_width - 1);
        auto [share0_y, y0, y1] = getCoords(y - half_csize, inv_csize, m_height - 1);
        auto [d0, d1, d2, d3] = combineShares(share0_x, share0_y);

        if(y1 >= m_height) {
            y1 = y0;
        }
        if(x1 >= m_width) {
            x1 = x0;
        }
        if(x0 < m_width && y0 < m_height) {
            particle_density[x0 + y0 * n] += d0;
        }
        if(x1 < m_width && y0 < m_height) {
            particle_density[x1 + y0 * n] += d1;
        }
        if(x1 < m_width && y1 < m_height) {
            particle_density[x1 + y1 * n] += d2;
        }
        if(x0 < m_width && y1 < m_height) {
            particle_density[x0 + y1 * n] += d3;
        }
    }

    if(particleRestDensity == 0.0) {
        auto sum = 0.0;
        auto numFluidCells = 0;

        for(int i = 0; i < m_num_cells; i++) {
            if(cell_type[i] == eCellTypes::Fluid) {
                sum += particle_density[i];
                numFluidCells++;
            }
        }

        if(numFluidCells > 0) {
            particleRestDensity = sum / numFluidCells;
        }
    }
}
void Fluid::transferVelocitiesToGrid(float flipRatio, Particles &particles)
{
    for(auto i = 0; i < m_num_cells; i++) {                                     
        cell_type[i] = (solid[i] == 1.0 ? eCellTypes::Solid : eCellTypes::Air); 
    }
    EMP_BENCHMARK_FUNC();
    auto inv_cell_size = 1.f / m_cell_size;
    auto half_cell_size = 0.5 * m_cell_size;

    prev_velocities = velocities;

    for(auto i = 0; i < Particles::max_particle_count; i++) {
        auto x = particles.position[i].x;
        auto y = particles.position[i].y;
        auto xi = /* std::clamp( */ floorf(x * inv_cell_size) /* , 0.f, m_width - 1.f) */;
        auto yi = /* std::clamp( */ floorf(y * inv_cell_size) /* , 0.f, m_height - 1.f) */;
        auto cellNr = xi + yi * m_width;
        if(cell_type[cellNr] == eCellTypes::Air) {
            cell_type[cellNr] = eCellTypes::Fluid;
        }
    }

    std::fill(velocities.begin(), velocities.end(), vec2f(0, 0));
    std::fill(velocities_diff.begin(), velocities_diff.end(), vec2f(0, 0));
    auto vel = [&](int idx, int component) -> float & {
        if(component == 0) {
            return velocities[idx].x;
        }
        return velocities[idx].y;
    };
    auto prev_vel = [&](int idx, int component) -> float & {
        if(component == 0) {
            return prev_velocities[idx].x;
        }
        return prev_velocities[idx].y;
    };
    auto diff = [&](int idx, int component) -> float & {
        if(component == 0) {
            return velocities_diff[idx].x;
        }
        return velocities_diff[idx].y;
    };

    for(auto component = 0; component < 2; component++) {
        auto offset_x = component == 0 ? 0.0 : half_cell_size;
        auto offset_y = component == 0 ? half_cell_size : 0.0;

        for(auto i = 0; i < Particles::max_particle_count; i++) {
            auto x = particles.position[i].x;
            auto y = particles.position[i].y;

            auto [share0_x, x0, x1] = getCoords(x - offset_x, inv_cell_size, m_width - 1);
            auto [share0_y, y0, y1] = getCoords(y - offset_y, inv_cell_size, m_height - 1);

            auto [d0, d1, d2, d3] = combineShares(share0_x, share0_y);

            auto bl = x0 + y0 * m_width;
            auto br = x1 + y0 * m_width;
            auto tr = x1 + y1 * m_width;
            auto tl = x0 + y1 * m_width;

            //  transfering weighted velocity
            auto pv = particles.velocity[i].x;
            if(component == 1) {
                pv = particles.velocity[i].y;
            }

            diff(bl, component) += d0;
            diff(br, component) += d1;
            diff(tr, component) += d2;
            diff(tl, component) += d3;

            vel(bl, component) += pv * d0;
            vel(br, component) += pv * d1;
            vel(tr, component) += pv * d2;
            vel(tl, component) += pv * d3;
        }
        for(auto i = 0; i < velocities.size(); i++) {
            if(diff(i, component) > 0.0) {
                vel(i, component) /= diff(i, component);
            }
        }
        for(auto j = 0; j < m_height; j++) {
            for(auto i = 0; i < m_width; i++) {
                auto solid = (cell_type[i + j * m_width] == eCellTypes::Solid);
                if(component == 0) {
                    if(solid || (i > 0 && cell_type[(i - 1) + j * m_width] == eCellTypes::Solid)) {
                        velocities[i + j * m_width].x = 0;
                    }
                } else {
                    if(solid || (j > 0 && cell_type[i + (j - 1) * m_width] == eCellTypes::Solid)) {
                        velocities[i + j * m_width].y = 0;
                    }
                }
            }
        }
    }
}
void Fluid::transferVelocitiesFromGrid(float flipRatio, Particles &particles)
{
    EMP_BENCHMARK_FUNC()
    auto inv_cell_size = 1.f / m_cell_size;
    auto half_cell_size = 0.5 * m_cell_size;

    for(auto component = 0; component < 2; component++) {
        auto offset_x = component == 0 ? 0.0 : half_cell_size;
        auto offset_y = component == 0 ? half_cell_size : 0.0;

        auto vel = [&](int idx) -> float & {
            if(component == 0) {
                return velocities[idx].x;
            }
            return velocities[idx].y;
        };
        auto prev_vel = [&](int idx) -> float & {
            if(component == 0) {
                return prev_velocities[idx].x;
            }
            return prev_velocities[idx].y;
        };
        auto diff = [&](int idx) -> float & {
            if(component == 0) {
                return velocities_diff[idx].x;
            }
            return velocities_diff[idx].y;
        };

        for(auto i = 0; i < Particles::max_particle_count; i++) {
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
            if(component == 1) {
                v = particles.velocity[i].y;
            }
            auto d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

            if(d > 0.0) {

                auto picV = (valid0 * d0 * vel(bl) + valid1 * d1 * vel(br) + valid2 * d2 * vel(tr) + valid3 * d3 * vel(tl)) / d;
                auto corr = (valid0 * d0 * (vel(bl) - prev_vel(bl)) + valid1 * d1 * (vel(br) - prev_vel(br)) +
                             valid2 * d2 * (vel(tr) - prev_vel(tr)) + valid3 * d3 * (vel(tl) - prev_vel(tl))) /
                            d;
                auto flipV = v + corr;
                auto vel_comp = (1.0 - flipRatio) * picV + flipRatio * flipV;
                if(component == 0) {
                    particles.velocity[i].x = vel_comp;
                } else {
                    particles.velocity[i].y = vel_comp;
                }
            }
        }
    }
}
void Fluid::transferBetweenGrids(std::vector<vec2f> &vel1, eCellTypes type1, std::vector<vec2f> &vel2, eCellTypes type2,
                                 float ratio)
{
    EMP_BENCHMARK_FUNC()
    auto n = m_width;
    auto oneOfTypes = [&](eCellTypes type) -> bool { return type == type1 || type == type2; };
    for(auto j = 0; j < m_height - 1; j++) {
        for(auto i = 0; i < m_width - 1; i++) {
            int cur = i + j * n;
            if(!oneOfTypes(cell_type[cur])) {
                continue;
            }

            for(int component = 0; component < 2; component++) {
                int other = (i + (1 - component)) + (j + component) * n;
                if(!oneOfTypes(cell_type[other])) {
                    continue;
                }
                if(cell_type[cur] == cell_type[other]) {
                    continue;
                }
                if(component == 0) {
                    auto combined = (vel2[cur].x - vel1[cur].x) * ratio;
                    vel1[cur].x += combined;
                } else {
                    auto combined = (vel2[cur].y - vel1[cur].y) * ratio;
                    vel1[cur].y += combined;
                }
            }
        }
    }
}
void Fluid::solveIncompressibility(float dt, eCellTypes expected_type, std::vector<vec2f> &vels, std::function<float(int)> solid,
                                   float density, int numIters, float overRelaxation, bool compensateDrift)
{
    EMP_BENCHMARK_FUNC()
    auto n = m_width;
    auto cp = density * m_cell_size / dt;

    for(auto iter = 0; iter < numIters; iter++) {

        for(auto j = 0; j < m_height; j++) {
            for(auto i = 0; i < m_width; i++) {

                if(cell_type[i + j * n] != expected_type) {
                    continue;
                }

                auto center = i + j * n;
                auto left = (i - 1) + j * n;
                auto right = (i + 1) + j * n;
                auto bottom = i + (j - 1) * n;
                auto top = i + (j + 1) * n;

                auto sx0 = 1.f - solid(left);
                auto sx1 = 1.f - solid(right);
                auto sy0 = 1.f - solid(bottom);
                auto sy1 = 1.f - solid(top);
                auto s = sx0 + sx1 + sy0 + sy1;
                if(s == 0.0) {
                    continue;
                }

                auto div = vels[right].x - vels[center].x + vels[top].y - vels[center].y;

                if(particleRestDensity > 0.0 && compensateDrift) {
                    auto k = 0.6;
                    auto compression = particle_density[i + j * n] - particleRestDensity;
                    if(compression < 0.0) {
                        div = div - k * compression;
                    }
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
void Fluid::simulate(Particles &particles, AABB sim_area, float dt, vec2f gravity, int numPressureIters, int numParticleIters,
                     float overRelaxation, bool compensateDrift)
{
    auto numSubSteps = 1;
    auto sdt = dt / numSubSteps;

    sim_area.setSize(sim_area.size() - vec2f(m_cell_size, m_cell_size) * 2.f);

    float col_time = 0.f;
    float fluid_time = 0.f;
    for(int step = 0; step < numSubSteps; step++) {
        Stopwatch local_stop;
        for(auto i = 0; i < m_num_cells; i++) {
            cell_type[i] = (solid[i] == 1.0 ? eCellTypes::Solid : eCellTypes::Air);
        }
        {
            for(int i = 0; i < numParticleIters; i++) {
                accelerate(particles, gravity);
                integrate(particles, sdt / (float)numParticleIters);
                collideWithGrid(particles, sdt / (float)numParticleIters);
                if(i == 0) {
                    collide(particles, sim_area);
                }
                constraint(particles, sim_area);
            }
        }
        transferVelocitiesToGrid(flipRatio, particles);
        updateParticleDensity(particles);
        prev_velocities = velocities;
        solveIncompressibility(
            sdt, eCellTypes::Fluid, velocities, [&](int idx) { return this->cell_type[idx] == eCellTypes::Solid; }, fluid_density,
            numPressureIters, overRelaxation, compensateDrift);
        advectAny(sdt, air_velocities, air_velocities, [&](vec2f &v) -> float & { return v.x; }, 0.f, m_cell_size / 2.f);
        advectAny(sdt, air_velocities, air_velocities, [&](vec2f &v) -> float & { return v.y; }, m_cell_size / 2.f, 0.f);

        float transfer_coef = (fluid_density + air_density) / fluid_density;
        transferBetweenGrids(air_velocities, eCellTypes::Air, velocities, eCellTypes::Fluid, transfer_coef);
        solveIncompressibility(
            sdt, eCellTypes::Air, air_velocities,
            [&](int idx) {
                if(this->cell_type[idx] == eCellTypes::Solid) {
                    return 1.f;
                }
                if(this->cell_type[idx] == eCellTypes::Fluid) {
                    float perc = particle_density[idx] / particleRestDensity;
                    return std::clamp<float>(perc, 0.f, 1.f);
                }
                return 0.f;
            },
            air_density, numPressureIters, overRelaxation, false);
        advectAny(sdt, smoke, air_velocities, [&](float &f) -> float & { return f; }, m_cell_size / 2.f, m_cell_size / 2.f);
        transferBetweenGrids(velocities, eCellTypes::Fluid, air_velocities, eCellTypes::Air, 1.f - transfer_coef);
        transferVelocitiesFromGrid(flipRatio, particles);
    }
}
void Fluid::draw(AABB area, Particles &particles, sf::RenderTarget &window, std::unordered_map<eCellTypes, Color> color_table)
{
    EMP_BENCHMARK_FUNC()
    std::vector<float> vals(m_num_cells, 0);
    updateParticleDensity(particles);
    for(auto i = 0; i < Particles::max_particle_count; i++) {
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
    for(int i = 0; i < m_height; i++) {
        for(int j = 0; j < m_width; j++) {
            auto type = cell_type[i * m_width + j];
            auto coord = sf::Vector2u(j, m_height - i - 1);
            if(!color_table.contains(type)) {
                img.setPixel(coord, sf::Color(255, 0, 255));
            } else if(type == eCellTypes::Fluid) {
                float val = (particle_density[j + i * m_width] / particleRestDensity);
                val = std::clamp(val, 0.1f, 1.f);
                uint8_t r = 105 * (val * 0.7f + 0.3f);
                uint8_t g = 120 * (val * 0.7f + 0.3f);
                uint8_t b = 220 * (val * 0.7f + 0.3f);
                img.setPixel(coord, sf::Color(r, g, b));
            } else if(type == eCellTypes::Solid) {
                img.setPixel(coord, color_table.at(type));
            } else if(smoke[i * m_width + j] > 0.15f) {
                uint8_t val = 255 * smoke[i * m_width + j];
                img.setPixel(coord, sf::Color(val, val, val));
            } else {
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
