#pragma once

#include "config/config.h"

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <omp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>

/************************************************************************/
/* print OpenPIC banner                                                 */
/************************************************************************/
void print_openpic_banner();

/************************************************************************/
/* timing functions                                                     */
/************************************************************************/
time_t elapsed_seconds();

// True when stdin is a live terminal. Gates every "press any key" style
// prompt: under a batch scheduler or redirected stdin such a prompt would
// block the job forever instead of informing anyone.
bool stdin_is_interactive();
bool stdout_is_interactive();

void print_console_status(const std::string& label, const std::string& detail = {});
void print_console_error(const std::string& message);
void print_simulation_progress(index_t step, index_t total_steps);
void finish_simulation_progress();
std::string format_duration(double seconds);

/************************************************************************/
/* messaging functions                                                  */
/************************************************************************/
// Profiling output levels:
// 0 = silent (no output)
// 1 = throttled progress with ETA and throughput
// 2 = progress + intermediate profiling timers
inline int verbosity_level()
{
    return PIC::Config::parameters().verbosity_level;
}

inline bool should_print_tm()
{
    // For intermediate timers (inside each step)
    return verbosity_level() >= 2;
}

inline bool should_print_step_marker()
{
    // For step markers (Step X of Y) at the start of each step
    return verbosity_level() >= 1;
}

template <class MsgT> void print_tm(const MsgT& msg, std::ostream& os = std::cout)
{
    if (PIC::Config::process_rank() != 0 || !should_print_tm())
        return;

    // Determine the width once for efficiency
    static size_t time_steps_width = fmt::formatted_size("{}", PIC::Config::time_steps());

    // {0:>{1}} -> Argument 0 (seconds), right-aligned, width from Argument 1 (width)
    fmt::print(os, "{:>{}} sec: {}\n", elapsed_seconds(), time_steps_width, msg);
}

template <class MsgPrefixT, class MsgSuffixT>
void print_tm(const MsgPrefixT& msg_prefix, const MsgSuffixT& msg_suffix, std::ostream& os = std::cout)
{
    if (PIC::Config::process_rank() != 0 || !should_print_tm())
        return;

    static size_t time_steps_width = fmt::formatted_size("{}", PIC::Config::time_steps());
    fmt::print(os, "{:>{}} sec: {}{}\n", elapsed_seconds(), time_steps_width, msg_prefix, msg_suffix);
}

template <class MsgT> void print(const MsgT& msg, std::ostream& os = std::cout)
{
    os << "\n" << msg << std::endl;
}

/************************************************************************/
/* Create data file name for saving                                     */
/************************************************************************/
std::string create_out_file_name(const std::string& prefix, const std::string& type_name, index_t step);

/************************************************************************/
/* Logging function                                                     */
/************************************************************************/
template <class MsgT> void log(const MsgT& msg, bool duplicate_on_stdout = true)
{
    const int thread_num = omp_get_thread_num();

    char log_file_name[50];
    sprintf(log_file_name, "opic_thread_%d.log", thread_num);
    std::ofstream ofs_log(log_file_name, std::ios_base::app);

    if (!ofs_log) {
        char error_msg[200];
        sprintf(error_msg, "Can't open \"%s\" file.", log_file_name);
        std::cout << error_msg;
        return;
    }

    ofs_log << std::setiosflags(std::ios_base::scientific);

    print(msg, ofs_log);

    if (duplicate_on_stdout) {
        print(msg);
    }
}
