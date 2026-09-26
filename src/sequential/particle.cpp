#include "particle.hpp"
#include "benchmark/benchmark.hpp"
#include "geometry_func.hpp"
#include "time.hpp"
#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/PrimitiveType.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <execution>
#include <iostream>
#include <stdexcept>
#include <unordered_set>
#include <vector>

float Particles::radius = 6.f;
float Particles::diameter = Particles::radius * 2.f;
uint32_t Particles::max_particle_count = 2500;

typedef std::vector<uint32_t> CompactVec;

struct Grid {
    std::vector<CompactVec> cells;
    uint32_t cols = 0, rows = 0;
    float cell_size = 1.f;
    vec2f origin{ 0, 0 };

    void reset(uint32_t c, uint32_t r, float cs, vec2f o)
    {
        cols = c;
        rows = r;
        cell_size = cs;
        origin = o;
        if(cells.size() != (size_t)cols * rows) {
            cells.assign((size_t)cols * rows, {});
        } else {
            for(auto &cell : cells) cell.clear();
        }
    }
    CompactVec       &at(int row, int col)       { return cells[(size_t)row * cols + col]; }
    const CompactVec &at(int row, int col) const { return cells[(size_t)row * cols + col]; }
};

static Grid s_grid;

static void buildGrid(Particles &p, AABB area, float cell_size)
{
    uint32_t inner_cols = (uint32_t)std::floor(area.size().x / cell_size) + 1;
    uint32_t inner_rows = (uint32_t)std::floor(area.size().y / cell_size) + 1;

    s_grid.reset(inner_cols + 2, inner_rows + 2, cell_size, area.min);

    for(uint32_t i = 0; i < Particles::max_particle_count; ++i) {
        int col = (int)((p.position[i].x - area.min.x) / cell_size) + 1;
        int row = (int)((p.position[i].y - area.min.y) / cell_size) + 1;
        col = std::clamp(col, 1, (int)s_grid.cols - 2);
        row = std::clamp(row, 1, (int)s_grid.rows - 2);
        s_grid.at(row, col).push_back(i);
    }
}

void cleanup(Particles &particles)
{
    delete[] particles.acceleration;
    delete[] particles.velocity;
    delete[] particles.position;
    delete[] particles.previous_position;
}

void accelerate(Particles &particles, vec2f gravity)
{
    EMP_BENCHMARK_FUNC();
    for(uint32_t i = 0; i < Particles::max_particle_count; i++) {
        particles.acceleration[i] += gravity;
    }
}

void integrate(Particles &particles, float dt)
{
    EMP_BENCHMARK_FUNC();
    for(uint32_t i = 0; i < Particles::max_particle_count; i++) {
        particles.previous_position[i] = particles.position[i];
        particles.position[i] += particles.velocity[i] * dt;
        particles.previous_position[i] -= particles.acceleration[i] * dt * dt;
        particles.velocity[i] += particles.acceleration[i] * dt;
        particles.acceleration[i] = { 0, 0 };
    }
}

void deriveVelocities(Particles &particles, float dt)
{
    EMP_BENCHMARK_FUNC();
    if(dt == 0.f) return;
    for(uint32_t i = 0; i < Particles::max_particle_count; i++) {
        particles.velocity[i] = (particles.position[i] - particles.previous_position[i]) / dt;
    }
}

void processCollision(Particles &particles, uint32_t i, uint32_t ii)
{
    vec2f diff = particles.position[ii] - particles.position[i];
    const float min_dist = particles.diameter;
    float r2 = qlen(diff);
    if(r2 < min_dist * min_dist && r2 > 1e-10f) {
        float l = std::sqrt(r2);
        vec2f n = diff / l;
        float c = (min_dist - l) * 0.5f;
        particles.position[i] -= n * c;
        particles.position[ii] += n * c;
    }
}

void collide(Particles &particles, AABB sim_area)
{
    EMP_BENCHMARK_FUNC();
    buildGrid(particles, sim_area, particles.diameter);

    for(uint32_t row = 1; row + 1 < s_grid.rows; ++row) {
        for(uint32_t col = 1; col + 1 < s_grid.cols; ++col) {
            const auto &cell = s_grid.at(row, col);
            if(cell.empty()) continue;

            static const int dr[4] = { 0, 1, 1, 1 };
            static const int dc[4] = { 1, 0, 1, -1 };

            for(size_t i = 0; i < cell.size(); ++i) {
                for(size_t j = i + 1; j < cell.size(); ++j) {
                    processCollision(particles, cell[i], cell[j]);
                }
                for(int k = 0; k < 4; ++k) {
                    const auto &other = s_grid.at(row + dr[k], col + dc[k]);
                    for(uint32_t j : other) {
                        processCollision(particles, cell[i], j);
                    }
                }
            }
        }
    }
}

void collide(Particles &particles)
{
    EMP_BENCHMARK_FUNC();
    for(uint32_t i = 0; i < Particles::max_particle_count; i++) {
        for(uint32_t ii = i + 1; ii < Particles::max_particle_count; ii++) {
            processCollision(particles, i, ii);
        }
    }
}
void constraint(Particles &particles, AABB area)
{
    EMP_BENCHMARK_FUNC();
    for(uint32_t i = 0; i < Particles::max_particle_count; i++) {
        particles.position[i].x = std::clamp(particles.position[i].x, area.min.x, area.max.x);
        particles.position[i].y = std::clamp(particles.position[i].y, area.min.y, area.max.y);
    }
}
void draw(Particles &particles, sf::RenderTarget &window, sf::Color color)
{
    EMP_BENCHMARK_FUNC();
    sf::VertexArray quads(sf::PrimitiveType::Triangles, 6 * Particles::max_particle_count);
    float winH = (float)window.getSize().y;

    for(uint32_t i = 0; i < Particles::max_particle_count; i++) {
        vec2f pos = particles.position[i];
        pos.y = winH - pos.y;

        quads[i * 6 + 0].position = { pos.x,                   pos.y + particles.radius };
        quads[i * 6 + 1].position = { pos.x + particles.radius, pos.y };
        quads[i * 6 + 2].position = { pos.x,                   pos.y - particles.radius };
        quads[i * 6 + 3].position = { pos.x,                   pos.y - particles.radius };
        quads[i * 6 + 4].position = { pos.x - particles.radius, pos.y };
        quads[i * 6 + 5].position = { pos.x,                   pos.y + particles.radius };

        for(int j = 0; j < 6; j++) {
            quads[i * 6 + j].color = color;
        }
    }
    window.draw(quads);
}
void init(Particles &particles, AABB screen_area, float spacing, int seed)
{
    particles.position          = new vec2f[Particles::max_particle_count];
    particles.previous_position = new vec2f[Particles::max_particle_count];
    particles.velocity          = new vec2f[Particles::max_particle_count];
    particles.acceleration      = new vec2f[Particles::max_particle_count];

    srand(seed);
    float w = screen_area.size().x - particles.radius * 2.f;
    int width = std::max(1, (int)(w / (particles.radius * 2.f * spacing)));
    for(uint32_t i = 0; i < Particles::max_particle_count; i++) {
        particles.position[i].x = (i % width) * particles.radius * 2.f * spacing + screen_area.min.x;
        particles.position[i].y = (i / width) * particles.radius * 2.f * spacing + screen_area.min.y;
        particles.previous_position[i] = particles.position[i];
        particles.velocity[i]     = { 0, 0 };
        particles.acceleration[i] = { 0, 0 };
    }
}
