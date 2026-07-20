#include "io/io_utilities.h"

#include <fmt/core.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <string>
#include <ctime>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

bool stdin_is_interactive()
{
#ifdef _WIN32
    return _isatty(_fileno(stdin)) != 0;
#else
    return isatty(fileno(stdin)) != 0;
#endif
}

bool stdout_is_interactive()
{
#ifdef _WIN32
    return _isatty(_fileno(stdout)) != 0;
#else
    return isatty(fileno(stdout)) != 0;
#endif
}

namespace {
using SteadyClock = std::chrono::steady_clock;
SteadyClock::time_point progress_started;
SteadyClock::time_point progress_last_rendered;
index_t last_progress_step = static_cast<index_t>(-1);
bool progress_active = false;
}

std::string format_duration(double seconds)
{
    if (!std::isfinite(seconds) || seconds < 0.0)
        return "--:--";
    const auto total = static_cast<long long>(std::llround(seconds));
    const auto hours = total / 3600;
    const auto mins = (total % 3600) / 60;
    const auto secs = total % 60;
    return hours > 0 ? fmt::format("{}:{:02}:{:02}", hours, mins, secs)
                     : fmt::format("{:02}:{:02}", mins, secs);
}

void print_console_status(const std::string& label, const std::string& detail)
{
    fmt::print("  {:<12}{}{}\n", label, detail.empty() ? "" : " ", detail);
}

void print_console_error(const std::string& message)
{
    finish_simulation_progress();
    fmt::print(stderr, "\n  ERROR       {}\n", message);
}

void print_simulation_progress(index_t step, index_t total_steps)
{
    if (!should_print_step_marker() || total_steps <= 0)
        return;
    if (!progress_active) {
        progress_started = SteadyClock::now();
        progress_last_rendered = progress_started;
        progress_active = true;
        last_progress_step = static_cast<index_t>(-1);
        fmt::print("\n  Simulation\n");
    }

    // Keep rendering off the hot path. Interactive output is capped at 4 Hz;
    // redirected output emits at most about 20 records for the whole run.
    // All expensive formatting happens only after this gate.
    static const bool terminal = stdout_is_interactive();
    const auto now = SteadyClock::now();
    if (step != total_steps && step != 1) {
        if (terminal) {
            constexpr auto refresh_period = std::chrono::milliseconds(250);
            if (now - progress_last_rendered < refresh_period)
                return;
        } else {
            const index_t interval = std::max<index_t>(1, total_steps / 20);
            if (last_progress_step != static_cast<index_t>(-1) && step - last_progress_step < interval)
                return;
        }
    }

    const double elapsed = std::chrono::duration<double>(now - progress_started).count();
    const double fraction = std::clamp(static_cast<double>(step) / total_steps, 0.0, 1.0);
    const bool estimate_ready = elapsed >= 0.5;
    const double eta = estimate_ready && step > 0 ? elapsed * (total_steps - step) / step : -1.0;
    const std::string rate = estimate_ready ? fmt::format("{:5.1f} step/s", step / elapsed)
                                            : "   -- step/s";
    constexpr int width = 28;
    const int filled = static_cast<int>(std::round(fraction * width));
    const std::string bar = std::string(filled, '=') +
                            (filled < width ? ">" + std::string(width - filled - 1, ' ') : "");
    const std::string line = fmt::format(
        "  [{:<{}}] {:6.2f}%  {}/{}  elapsed {}  ETA {}  {}",
        bar, width, fraction * 100.0, step, total_steps,
        format_duration(elapsed), format_duration(eta), rate);
    if (terminal)
        fmt::print("\r{:<140}", line);
    else
        fmt::print("{}\n", line);
    if (terminal)
        std::fflush(stdout);
    last_progress_step = step;
    progress_last_rendered = now;
}

void finish_simulation_progress()
{
    if (progress_active && stdout_is_interactive())
        fmt::print("\n");
    progress_active = false;
}

void print_openpic_banner()
{
    fmt::print("{}", R"(
  OPEN-PIC  3.0
  3D Hybrid Particle-in-Cell Plasma Simulator
  ------------------------------------------------------------
  Institute of Radiophysics and Electronics, NAS RA
  irphe.am
)");
}

time_t elapsed_seconds()
{
    static time_t start = time(0);
    return (time(0) - start);
}

std::string create_out_file_name(const std::string& prefix, const std::string& type_name, index_t step)
{
    static size_t time_steps_width = fmt::formatted_size("{}", PIC::Config::time_steps());
    return fmt::format("{}_{}_{:0{}}.dat", prefix, type_name, step, time_steps_width);
}
