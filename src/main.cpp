#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/PrimitiveType.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/Audio.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include <SFML/Window/Keyboard.hpp>
#include <any>
#include <cmath>
#include <iomanip>
#include <numeric>
#include "benchmark/benchmark.hpp"
#include "fluid.hpp"
#include "particle.hpp"
#include "geometry_func.hpp"
#include "time.hpp"
#include "vec2.hpp"

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
        Press 1 for grabbing, 2 for smoke generation and 3 for solid block drawing.
        Press G to toggle grid view.
        Press P to toggle particle view.
        Press Q to toggle pressure view.
)""";
    enum MOUSE_MODE { HOLD, SMOKE_GEN, SOLID_GEN };
    MOUSE_MODE cur_mouse_mode = HOLD;
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

    std::cout << "particle radius: " << particles.radius << "\n";
    std::cout << "particle count: " << Particles::max_particle_count << "\n";
    std::cout << "num of particle iters: " << numParticleIters << "\n";
    std::cout << "num of fluid iters: " << numFluidIters << "\n";
    std::cout << "overrelaxation: " << overrelaxation << "\n";
    bool drawParticles = false;
    bool drawGrid = true;
    bool drawPressure = false;
    bool shortcut_pressed = 0;
    std::map<sf::Keyboard::Key, bool *> display_shortcuts = {
        { sf::Keyboard::Key::P, &drawParticles },
        { sf::Keyboard::Key::G, &drawGrid      },
        { sf::Keyboard::Key::X, &drawPressure  },
    };

    float total_time = 0;
    Clock deltaClock;
    bool shouldReport = true;
    Stopwatch report_clock;
    uint32_t sample_count;
    report_clock.restart();
    float brush_size = fluid.cell_size();

    while(window.isOpen()) {
        float deltaTime = deltaClock.restart().asSeconds();
        while(const std::optional event = window.pollEvent()) {
            if(event->is<sf::Event::Closed>()) {
                window.close();
            }
            if(event->is<sf::Event::MouseWheelScrolled>()) {
                auto wheel = event->getIf<sf::Event::MouseWheelScrolled>();
                const float scroll_speed = 10.f;  //  X cells per second
                if(wheel->delta == 0.f) {
                } else if(wheel->delta > 0.f) {
                    brush_size += fluid.cell_size() * deltaTime * scroll_speed;
                } else {
                    brush_size -= fluid.cell_size() * deltaTime * scroll_speed;
                }
                brush_size = std::max(brush_size, 1.f);
            }
        }

        //  mouse controls
        static vec2f last_mouse_pos;
        auto posi = sf::Mouse::getPosition(window);
        vec2f mouse_pos = { (float)posi.x, (float)posi.y };
        mouse_pos.y = window.getSize().y - mouse_pos.y;
        vec2f mouse_dir = mouse_pos - last_mouse_pos;
        if(sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)) {
            switch(cur_mouse_mode) {
                case MOUSE_MODE::HOLD:
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
                    break;
                case MOUSE_MODE::SOLID_GEN:
                case MOUSE_MODE::SMOKE_GEN: {
                    int x = mouse_pos.x / fluid.cell_size();
                    int y = mouse_pos.y / fluid.cell_size();
                    float hsize = brush_size / 2.f / fluid.cell_size() + 1;
                    for(int j = y - hsize; j < y + hsize; j++) {
                        for(int i = x - hsize; i < x + hsize; i++) {
                            if(length(mouse_pos - vec2f(i + 0.5f, j + 0.5f) * fluid.cell_size()) > brush_size) {
                                continue;
                            }
                            if(cur_mouse_mode == MOUSE_MODE::SOLID_GEN) {
                                fluid.solid[j * fluid.width() + i] = 1.f;
                            } else if(cur_mouse_mode == MOUSE_MODE::SMOKE_GEN) {
                                fluid.smoke[j * fluid.width() + i] = 1.f;
                            }
                        }
                    }
                } break;
            }
        }
        //  other controls
        bool any_shortcut_pressed = false;
        for(auto m : display_shortcuts) {
            if(sf::Keyboard::isKeyPressed(m.first)) {
                if(!shortcut_pressed) {
                    *m.second = !*m.second;
                }
                any_shortcut_pressed = true;
                shortcut_pressed = true;
            }
        }
        if(!any_shortcut_pressed) {
            shortcut_pressed = false;
        }
        if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num1)) {
            cur_mouse_mode = MOUSE_MODE::HOLD;
        }
        if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num2)) {
            cur_mouse_mode = MOUSE_MODE::SMOKE_GEN;
        }
        if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Num3)) {
            cur_mouse_mode = MOUSE_MODE::SOLID_GEN;
        }

        //  simulation
        total_time += deltaTime;

        fluid.simulate(particles, screen_area, deltaTime, vec2f(0, -1000.f), numFluidIters, numParticleIters, overrelaxation,
                       pushOut);
        sample_count += 1;
        if(report_clock.getElapsedTime() > raporting_interval && shouldReport && EMP_BENCHMARK) {
            auto total_time = BenchmarkRegistry().get().getMeasurement("ROOT");
            std::cout << "FPS:\t" << std::setprecision(2) << 1.0 / (total_time / sample_count) << '\n';
            BenchmarkRegistry().get().print(sample_count);
            BenchmarkRegistry().get().reset();
            report_clock.restart();
            sample_count = 0;
        }

        window.clear();

        //  drawing
        if(drawGrid) {
            std::unordered_map<eCellTypes, Color> color_table = {
                { eCellTypes::Air,   Color(0,   0,   50)  },
                { eCellTypes::Solid, Color(100, 100, 100) },
                { eCellTypes::Fluid, Color(70,  100, 220) },
            };
            fluid.draw(screen_area, particles, window, color_table, drawPressure);
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
        window.display();
        last_mouse_pos = mouse_pos;
    }
    cleanup(particles);

    return 0;
}
