# Add the minimum required CMake version at the very top
cmake_minimum_required(VERSION 3.10)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(GenomeAssembler)

# 1. Create a library from your source code
add_library(assembler_lib
        src/compare.cpp
        src/convert.cpp
        src/bitvector.cpp
        src/graphdbj.cpp
)

# Move definitions to after target creation
target_compile_definitions(assembler_lib PUBLIC -D_GLIBCXX_USE_CXX11_ABI=1)

# 2. Tell the library where to find its own headers
target_include_directories(assembler_lib PUBLIC src)

# 3. Create your main application executable
add_executable(assembler_app
        src/main.cpp
)

# Move definitions to after target creation
target_compile_definitions(assembler_app PRIVATE -D_GLIBCXX_USE_CXX11_ABI=1)

# 4. Link your main app against your library
target_link_libraries(assembler_app PRIVATE assembler_lib)


# ----------------------------------------------------
# Testing Configuration
# ----------------------------------------------------

# 5. Enable testing for the project
enable_testing()

# 6. Use FetchContent to obtain GoogleTest if not available
include(FetchContent)

FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/heads/main.zip
        DOWNLOAD_EXTRACT_TIMESTAMP true
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

# 7. Create a new directory for test files
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

# Define macro with path to test data directory so tests can find files reliably
target_compile_definitions(run_tests PRIVATE TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/Tests_Et_Ref")

# 10. Add your test to CTest
include(GoogleTest)
gtest_discover_tests(run_tests)