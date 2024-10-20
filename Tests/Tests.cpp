/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */
// Doctest header
#ifdef AMREX_USE_SYCL
// Intel's GPU runtime uses SIGSEGV to trigger migration of managed memory
#define DOCTEST_CONFIG_NO_POSIX_SIGNALS
#endif
#define DOCTEST_CONFIG_IMPLEMENT
#define DOCTEST_CONFIG_NO_UNPREFIXED_OPTIONS
#include "doctest.h"

#include "TestCases.hpp" // Test cases are defined here
#include "doctestCLIArgs.hpp"
#include "doctestOutput.hpp"

// system headers
#include <iomanip>
#include <iostream>

#include "AMReX.H"
#include "AMReX_REAL.H"
#include "AMReX_ccse-mpi.H"

namespace doctest
{
// Unfortunately the following has to be global and non-const
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
CLIArgs cli_args;
} // namespace doctest

// NOLINTBEGIN(bugprone-exception-escape)
int main(int argc, char *argv[])
{
#ifdef BL_USE_MPI
    // We can only initialize and finalize MPI once so do it here
    MPI_Init(&argc, &argv);
#endif
    doctest::cli_args.set(argv);

    doctest::Context doctest_context(argc, argv);
#ifdef BL_USE_MPI
    doctest_context.setCout(
        &doctest::hide_output_from_non_zero_ranks(std::cout));
#endif

    // Default AMReX verbosity to 0 to avoid the "Initialized", "Finalized" and
    // memory usage messages
    amrex::system::verbose = 0;

    // increase output precision
    constexpr int output_precision = 17;
    std::cout << std::setprecision(output_precision);

    // also increase precision of doctest's stringmaker
    std::ostream *doctest_stringstream = doctest::detail::tlssPush();
    doctest_stringstream->precision(output_precision);
    doctest::detail::tlssPop();

    int result = doctest_context.run();

#ifdef BL_USE_MPI
    MPI_Finalize();
#endif

    return result;
}
// NOLINTEND(bugprone-exception-escape)