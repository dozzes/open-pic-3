#include "config/config.h"
#include "lua/use_lua.h"
#include "lua/bind_to_lua.h"
#include "io/io_utilities.h"
#include "particles/particles.h"
#include "grid/grid.h"
#include "core/simulate.h"
#include "mpi/mpi_runtime.h"

#include <fmt/core.h>

#include <stdexcept>
#include <string>

namespace {

struct CommandLine
{
    bool        mpi_enabled = false;
    const char* lua_script  = nullptr;
};

CommandLine parse_command_line(int argc, char* argv[])
{
    CommandLine command_line;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-mpi") {
            command_line.mpi_enabled = true;
        } else if (!command_line.lua_script) {
            command_line.lua_script = argv[i];
        } else {
            throw std::runtime_error("Too many command-line arguments.");
        }
    }

    return command_line;
}

void print_usage(bool wait_for_key)
{
    std::cerr << "Usage: open-pic [-mpi] <setup-script.lua>\n"
                 "\nOptions:\n"
                 "  -mpi    Run under MPI (launch with mpiexec/mpirun)\n";
    if (wait_for_key && stdin_is_interactive()) {
        // Only meaningful for a double-clicked serial run, where the console
        // window would otherwise close before the message can be read; an
        // MPI launch (mpiexec) has no single interactive console to wait on,
        // and a batch job with redirected stdin must never block on a prompt.
        std::cerr << "Press any key to exit ...";
        std::cin.get();
    }
}

} // namespace

int main(int argc, char* argv[])
{
    bool mpi_requested = false;
    int  exit_code     = 0;

    try {
        CommandLine command_line = parse_command_line(argc, argv);
        mpi_requested            = command_line.mpi_enabled;

        if (!command_line.lua_script) {
            print_usage(!mpi_requested);
            return -1;
        }

        if (mpi_requested) {
            // MPI_Init must precede all other MPI calls, including any that print output.
            PIC::Mpi::init(&argc, &argv);
            PIC::Config::parameters().process_rank = PIC::Mpi::rank();
            PIC::Config::parameters().process_num  = PIC::Mpi::size();
        }

        // In replicated-grid MPI mode every rank runs this whole function, so
        // all console/log output not specific to a single rank is gated on
        // is_root() to avoid N-way duplicated banners, prompts, and logs.
        if (PIC::Mpi::is_root())
            print_openpic_banner();

        auto& lua = get_lua_state();

        Particles particles;
        Grid      grid;
        if (mpi_requested)
            // Partition must be set before bind_to_lua so that Lua injection callbacks
            // see the correct process_idx and skip particle creation on non-owner ranks.
            particles.use_rank_partition(PIC::Mpi::rank(), PIC::Mpi::size());

        if (PIC::Mpi::is_root())
            print_console_status("Config", command_line.lua_script);

        if (!bind_to_lua(lua, command_line.lua_script, grid, particles)) {
            if (PIC::Mpi::is_root()) {
                print_console_error("Lua configuration failed.");
                // Interactive pause only makes sense for a serial, console-attached
                // run; under MPI, Mpi::abort below tears down all ranks instead,
                // and with redirected stdin the prompt would hang a batch job.
                if (!mpi_requested && stdin_is_interactive()) {
                    std::cerr << "\nPress any key to exit ...";
                    std::cin.get();
                }
            }
            if (mpi_requested)
                PIC::Mpi::abort(-1);
            return -1;
        }

        if (PIC::Mpi::is_root()) {
            print_console_status("Mode", mpi_requested
                ? fmt::format("MPI ({} ranks)", PIC::Config::parameters().process_num)
                : "serial");
            PIC::Config::to_stream(std::cout);
            std::ofstream ofs_params("opic_config.txt");
            if (ofs_params) {
                ofs_params << std::boolalpha;
                PIC::Config::to_stream(ofs_params);
            }
        }

        time_t start_time = elapsed_seconds();

        // simulate() handles both modes itself via Mpi::enabled().
        PIC::simulate(grid, particles);

        time_t run_time = (elapsed_seconds() - start_time);

        if (PIC::Mpi::is_root()) {
            finish_simulation_progress();
            fmt::print("\n  SUCCESS     Simulation completed in {}\n",
                       format_duration(static_cast<double>(run_time)));
            PIC::Config::logger() << "\nRun time: " << run_time << " sec.\nFinished!";
        }
    } catch (const sol::error& e) {
        // sol::error derives from std::exception, so it must be caught first;
        // e.what() already holds the full Lua error message (sol2 copies it
        // into the exception at throw time), so no separate Lua-stack read
        // is needed or safe this far from the original throw site.
        if (PIC::Mpi::is_root()) {
            log(e.what(), false);
            print_console_error(e.what());
        }
        // Mpi::abort broadcasts the error across all ranks; a bare return would
        // leave other ranks spinning indefinitely in MPI_Allreduce.
        if (mpi_requested)
            PIC::Mpi::abort(-1);
        exit_code = -1;
    } catch (const std::exception& e) {
        if (PIC::Mpi::is_root()) {
            log(e.what(), false);
            print_console_error(e.what());
        }
        if (mpi_requested)
            PIC::Mpi::abort(-1);
        exit_code = -1;
    } catch (...) {
        if (PIC::Mpi::is_root())
            print_console_error("An unknown error occurred.");
        if (mpi_requested)
            PIC::Mpi::abort(-1);
        exit_code = -1;
    }

    // All catch blocks fall through to here rather than returning directly,
    // so MPI_Finalize always runs for a rank that reached MPI_Init — skipping
    // it leaves the process hanging or corrupts the next MPI job on shared
    // hardware. The early "no script given" path returns before MPI_Init, so
    // it does not need this and returns directly above.
    if (mpi_requested)
        PIC::Mpi::finalize();

    return exit_code;
}
