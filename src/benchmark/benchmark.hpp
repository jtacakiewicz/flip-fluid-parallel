#ifndef EMP_BENCHMARK_HPP
#define EMP_BENCHMARK_HPP
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
struct MeasurementGroup {
    double sum = 0.0;
    std::unordered_map<std::string, MeasurementGroup> subgroups;
    void print_children(double root_t, uint32_t sample_count = 1, std::ostream &os = std::cout, std::string prefix = "") const
    {
        for(const auto &g : subgroups) {
            double lpercentage = g.second.sum / this->sum * 100.0;
            double gpercentage = g.second.sum / root_t * 100.0;
            os << prefix << g.first << ":\t" << g.second.sum / sample_count << "s (l: " << std::setprecision(2) << lpercentage
               << "%; g: " << gpercentage << "%)\n";
            g.second.print_children(root_t, sample_count, os, prefix + '\t');
        }
    }
    void updateSum()
    {
        for(auto &g : subgroups) {
            g.second.updateSum();
            sum += g.second.sum;
        }
    }
};
class BenchmarkRegistry : private MeasurementGroup {
public:
    static BenchmarkRegistry &get()
    {
        static BenchmarkRegistry s_instance;
        return s_instance;
    }
    void record(std::string name, double value)
    {
        auto group = this->getSmallestGroup(name);
        if(group) {
            group->sum += value;
        }
        m_changed = true;
    }
    void print(uint32_t sample_count = 1, std::ostream &os = std::cout)
    {
        calcChange();
        os << "ROOT:\t" << this->sum / sample_count << "s\n";
        std::string prefix = "\t";
        this->print_children(this->sum, sample_count, os, prefix);
    }
    double getMeasurement(std::string name)
    {
        calcChange();
        if(name == "ROOT") {
            return this->sum;
        }
        auto group = this->cgetSmallestGroup(name);
        if(!group) {
            return -1;
        }
        return group->sum;
    }
    void reset()
    {
        this->sum = 0.0;
        this->subgroups.clear();
    }

private:
    bool m_changed = false;
    void calcChange()
    {
        if(!m_changed) {
            return;
        }
        m_changed = false;
        this->updateSum();
    }
    const MeasurementGroup *cgetSmallestGroup(std::string name) const
    {
        size_t left = 0;
        size_t right = left;
        const MeasurementGroup *group = this;
        do {
            right = name.find("::", left);
            auto subgroup = name.substr(left, right - left);
            if(!group->subgroups.contains(subgroup)) {
                return nullptr;
            }
            group = &group->subgroups.at(subgroup);
            left = right + 2;
        } while(right != std::string::npos);
        return group;
    }
    MeasurementGroup *getSmallestGroup(std::string name)
    {
        size_t left = 0;
        size_t right = left;
        MeasurementGroup *group = this;
        do {
            right = name.find("::", left);
            auto subgroup = name.substr(left, right - left);
            group = &group->subgroups[subgroup];
            left = right + 2;
        } while(right != std::string::npos);
        return group;
    }
};

class BenchmarkStopwatch {
    typedef std::chrono::duration<double> duration_t;
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> time_point_t;
    time_point_t m_start;
    time_point_t m_stop;
    std::string m_name;

public:
    BenchmarkStopwatch(std::string name)
        : m_name(name)
    {
    }
    void start() { m_start = std::chrono::high_resolution_clock::now(); }
    void stop()
    {
        m_stop = std::chrono::high_resolution_clock::now();
        auto dur = m_stop - m_start;
        double elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(dur).count()) / 1e9;
        BenchmarkRegistry::get().record(m_name, elapsed);
    }
};
class BenchmarkScoped : private BenchmarkStopwatch {
public:
    BenchmarkScoped(std::string name)
        : BenchmarkStopwatch(name)
    {
        this->start();
    }
    ~BenchmarkScoped() { this->stop(); }
};

#ifndef EMP_BENCHMARK 
#define EMP_BENCHMARK 1
#endif

#if EMP_BENCHMARK
#define EMP_CONCAT(a, b)                       EMP_CONCAT_INNER(a, b)
#define EMP_CONCAT_INNER(a, b)                 a##b
#define EMP_UNIQUE_NAME(base)                  EMP_CONCAT(EMP_CONCAT(base, __COUNTER__), EMP)
#define EMP_BENCHMARK_FUNC_INNER(uname)        BenchmarkScoped uname(__FILE__ + std::string("::") + __func__);
#define EMP_BENCHMARK_FUNCn_INNER(name, uname) BenchmarkScoped uname(__FILE__ + std::string("::") + __func__ + "::" + name);

#define EMP_BENCHMARK_FUNCn(name) EMP_BENCHMARK_FUNCn_INNER(name, EMP_UNIQUE_NAME(__func__))
#define EMP_BENCHMARK_FUNC()      EMP_BENCHMARK_FUNC_INNER(EMP_UNIQUE_NAME(__func__))
#else
#define EMP_CONCAT(a, b)
#define EMP_CONCAT_INNER(a, b)
#define EMP_UNIQUE_NAME(base)
#define EMP_BENCHMARK_FUNC_INNER(uname)
#define EMP_BENCHMARK_FUNCn_INNER(name, uname)

#define EMP_BENCHMARK_FUNCn(name)
#define EMP_BENCHMARK_FUNC()
#endif
#endif
