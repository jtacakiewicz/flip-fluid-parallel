#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/PrimitiveType.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/Audio.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <cmath>
#include <iomanip>
#include <numeric>
#include "benchmark/benchmark.hpp"
#include "fluid.hpp"
#include "particle.hpp"
#include "geometry_func.hpp"
#include "time.hpp"

#include <SFML/Window/Mouse.hpp>
#include <iostream>
#include <stdexcept>
using namespace sf;

int main(int argc, char **argv)
{
    float w = 1280;
    float h = 720;
    int numParticleIters = 16;
    int numFluidIters = 32;
    float overrelaxation = 1.9f;
    float cell_size_scale = 3.f;
    float raporting_interval = 1.f;
    float spacing = 1.f;
    bool pushOut = true;
    Particles particles;
    std::string help_msg = R"""(
Usage: fluid-sim [OPTIONS]
Options:
        -w --width           [width] 
        -q --height          [height] 
        -c --cell-size       [cell scale] 
        -r --radius          [particle radius] 
        -n --num             [amount of particles]
        -i --particle-iters  [num of iterations for particle solver]
        -f --fluid-iters     [num of iterations for fluid solver]
        -x --overrelaxation  [overrelaxation coef]
        -d --correct-drift   [1/0 should drift be corrected]
        -t --raport-interval [raporting time in seconds]
        -s --spacing         [spacing between particles (scale of radius)]
When using application:
        Click and hold mouse to interact,
        Press G to toggle grid view,
        Press P to toggle particle view
)""";
    try {
        for(int i = 1; i < argc; i += 2) {
            std::string flag = argv[i];
            std::string arg = i + 1 < argc ? argv[i + 1] : "";
            if(flag == "-w" || flag == "--width") {
                w = stoi(arg);
            } else if(flag == "-q" || flag == "--height") {
                h = stoi(arg);
            } else if(flag == "-c" || flag == "--cell-size") {
                cell_size_scale = stof(arg);
            } else if(flag == "-r" || flag == "--radius") {
                particles.radius = stof(arg);
            } else if(flag == "-n" || flag == "--num") {
                particles.max_particle_count = stoi(arg);
            } else if(flag == "-i" || flag == "--particle-iters") {
                numParticleIters = stoi(arg);
            } else if(flag == "-f" || flag == "--fluid-iters") {
                numFluidIters = stoi(arg);
            } else if(flag == "-x" || flag == "--overrelaxation") {
                overrelaxation = stof(arg);
            } else if(flag == "-d" || flag == "--correct-drift") {
                pushOut = stoi(arg);
            } else if(flag == "-t" || flag == "--raport-interval") {
                raporting_interval = stof(arg);
            } else if(flag == "-s" || flag == "--spacing") {
                spacing = stof(arg);
            } else if(flag == "-h" || flag == "--help") {
                std::cout << help_msg;
                return 0;
            } else {
                throw std::invalid_argument("unrecognized flag");
            }
        }
    } catch(...) {
        printf("incorrect arguments were given!");
        printf("%s", help_msg.c_str());
        return 1;
    }
    auto fluid_cell_size = particles.diameter * cell_size_scale;
    AABB screen_area = AABB::CreateMinSize({ 0, 0 }, { w, h });
    RenderWindow window(VideoMode(Vector2u(w, h)), "demo");
    auto area = screen_area;
    area.setSize(area.size() * 0.8f);
    init(particles, area, spacing);

    auto fluid_size = screen_area.size() / fluid_cell_size;
    Fluid fluid(fluid_cell_size, fluid_size.x, fluid_size.y);

    auto dispNameValue = [&](std::string name, auto value, bool isLast = false) { std::cout << name << ": " << value << "\n"; };
    dispNameValue("particle radius", particles.radius);
    dispNameValue("particle count", Particles::max_particle_count);
    dispNameValue("num of particle iters", numParticleIters);
    dispNameValue("num of fluid iters", numFluidIters);
    dispNameValue("overrelaxation", overrelaxation);

    float total_time = 0;
    Clock deltaClock;
    bool shouldReport = true;
    Stopwatch report_clock;
    uint32_t sample_count;
    report_clock.restart();
    while(window.isOpen()) {
        while(const std::optional event = window.pollEvent()) {
            if(event->is<sf::Event::Closed>()) {
                window.close();
            }
        }
        float deltaTime = deltaClock.restart().asSeconds();

        static vec2f last_mouse_pos;
        auto posi = sf::Mouse::getPosition(window);
        vec2f mouse_pos = { (float)posi.x, (float)posi.y };
        mouse_pos.y = window.getSize().y - mouse_pos.y;
        vec2f mouse_dir = mouse_pos - last_mouse_pos;
        const float brush_size = 50.f;
        if(sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)) {
            for(int i = 0; i < Particles::max_particle_count; i++) {
                auto scalar = length(mouse_dir) / deltaTime;
                vec2f norm;
                if(qlen(mouse_dir) == 0.f) {
                    norm = { 0, 0 };
                } else {
                    norm = normal(mouse_dir);
                }
                if(length(mouse_pos - particles.position[i]) < brush_size) {
                    particles.velocity[i] = norm * std::clamp(scalar, 0.f, 1000.f);
                }
            }
        }

        total_time += deltaTime;

        fluid.simulate(particles, screen_area, deltaTime, vec2f(0, -1000.f), numFluidIters, numParticleIters, overrelaxation,
                       pushOut);
        sample_count += 1;
        if(report_clock.getElapsedTime() > raporting_interval && shouldReport) {
            auto total_time = BenchmarkRegistry().get().getMeasurement("ROOT");
            std::cout << "FPS:\t" << std::setprecision(2) << 1.0 / (total_time / sample_count) << '\n';
            BenchmarkRegistry().get().print(sample_count);
            BenchmarkRegistry().get().reset();
            report_clock.restart();
            sample_count = 0;
        }

        window.clear();
        static bool drawParticles = false;
        static bool drawGrid = true;
        static bool pressed = 0;
        if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::P)) {
            if(!pressed) {
                drawParticles = !drawParticles;
            }
            pressed = true;
        } else if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::G)) {
            if(!pressed) {
                drawGrid = !drawGrid;
            }
            pressed = true;
        } else {
            pressed = false;
        }
        if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::S)) {
            int x = mouse_pos.x / fluid.cell_size();
            int y = mouse_pos.y / fluid.cell_size();
            fluid.solid[y * fluid.width() + x] = 1.f;
        }

        if(drawGrid) {
            fluid.draw(screen_area, particles, window);
        }
        if(drawParticles) {
            draw(particles, window, Color(70, 70, 250));
        }
        sf::CircleShape cs(brush_size);
        cs.setOrigin({ brush_size, brush_size });
        cs.setPosition({ mouse_pos.x, screen_area.size().y - mouse_pos.y });
        cs.setFillColor(Color(0, 0, 0, 0));
        cs.setOutlineColor(Color(255, 255, 255));
        cs.setOutlineThickness(2.f);
        window.draw(cs);
        vec2f origin = { 500, 500 };
        auto collision = fluid.findCollision(origin, mouse_pos - origin, { eCellTypes::Solid });
        sf::CircleShape circ(5.f);
        circ.setOrigin({ 5.f, 5.f });
        circ.setPosition(vec2f(origin.x, h - origin.y));
        window.draw(circ);
        if(collision) {
            auto cp = collision->collision_point;

            circ.setPosition(vec2f(cp.x, h - cp.y));
            sf::Vertex verts[2];
            for(int i = 0; i < 2; i++) {
                verts[i].color = sf::Color(0xffffff);
                cp = cp + collision->normal * 10.f;
                verts[i].position = vec2f(cp.x, h - cp.y);
            }
            window.draw(circ);
            window.draw(verts, 2U, sf::PrimitiveType::Lines);
        }
        window.display();
        last_mouse_pos = mouse_pos;
    }
    cleanup(particles);

    return 0;
}
