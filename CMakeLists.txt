cmake_minimum_required(VERSION 3.14)
project(GenomeAssembler)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. Create a library from your source code
# This library will contain compare.cpp and any other logic files.
add_library(assembler_lib
        src/compare.cpp
        # Add other .cpp files from src/ here, but NOT main.cpp
)

# 2. Tell the library where to find its own headers
target_include_directories(assembler_lib PUBLIC src)

# 3. Create your main application executable
add_executable(assembler_app
        src/main.cpp
)

# 4. Link your main app against your library
target_link_libraries(assembler_app PRIVATE assembler_lib)


# ----------------------------------------------------
# Testing Configuration
# ----------------------------------------------------

# 5. Enable testing for the project
enable_testing()

# 6. Find the GoogleTest package
# You may need to install this (e.g., sudo apt-get install libgtest-dev)
find_package(GTest REQUIRED)

# 7. Create a new directory for test files
# This command is not strictly needed, but assumes you created tests/
file(GLOB TEST_FILES "tests/*.cpp")

# 8. Add the test executable
add_executable(run_tests
        ${TEST_FILES}
)

# 9. Link the test executable against GTest and your assembler library
target_link_libraries(run_tests PRIVATE
        assembler_lib     # <-- Link to your code
        GTest::gtest      # <-- Link to GTest
        GTest::gtest_main # <-- Link to GTest main
)

# 10. Add your test to CTest
include(GoogleTest)
gtest_discover_tests(run_tests)