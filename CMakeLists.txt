cmake_minimum_required(VERSION 3.16)
project(NanoKMC VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(NANOKMC_PAPER_BUILD
       "Use the compiler settings frozen for the NanoKMC paper benchmark"
       OFF)
option(NANOKMC_WARNINGS
       "Enable useful compiler warnings for NanoKMC sources"
       ON)

set(NANOKMC_SOURCES
    src/main.cpp
    src/CalcClass.cpp
    src/SystemClass.cpp
    src/SystemClassEval.cpp
    src/SystemKMCHomogenous.cpp
    src/SystemKMCClassical.cpp
    src/SystemKMCActiveFilteredBase.cpp
    src/SystemKMCActiveFilteredGeneric.cpp
    src/SystemKMCActiveFilteredBinaryNN.cpp
    src/SystemKMCPartialFilterOptimized.cpp
    src/SystemKMCExactClassOptimized.cpp
    src/SystemKMCRateCategoryOptimized.cpp
    src/internal/SystemKMCBinaryRateBase.cpp
    src/internal/SystemKMCPartialFilterCore.cpp
    src/internal/SystemKMCExactClassCore.cpp
)

add_executable(nanokmc ${NANOKMC_SOURCES})
target_include_directories(nanokmc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_compile_definitions(nanokmc PRIVATE
    NANOKMC_VERSION="${PROJECT_VERSION}"
    $<$<BOOL:${NANOKMC_PAPER_BUILD}>:NANOKMC_PAPER_BUILD_ENABLED=1>
)

if (NANOKMC_WARNINGS)
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(nanokmc PRIVATE
            -Wall -Wextra -Wpedantic
            -Wno-unused-parameter
            -Wno-unused-variable
            -Wno-unused-but-set-variable
            -Wno-sign-compare
            -Wno-parentheses)
    elseif (MSVC)
        target_compile_options(nanokmc PRIVATE /W3)
    endif()
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" AND NANOKMC_PAPER_BUILD)
    # Frozen single-core benchmark setting. Keep this explicit so a CMake
    # Release default such as -O3 does not silently change the paper build.
    target_compile_options(nanokmc PRIVATE -O1)
    target_compile_definitions(nanokmc PRIVATE NDEBUG)
elseif (MSVC AND NANOKMC_PAPER_BUILD)
    target_compile_options(nanokmc PRIVATE /O1)
    target_compile_definitions(nanokmc PRIVATE NDEBUG)
endif()

include(CTest)
if (BUILD_TESTING)
    add_test(NAME nanokmc_help COMMAND nanokmc --help)
    add_test(NAME nanokmc_version COMMAND nanokmc --version)
    add_test(NAME nanokmc_list_solvers COMMAND nanokmc --list-solvers)

    find_package(Python3 COMPONENTS Interpreter QUIET)
    if (Python3_Interpreter_FOUND)
        add_test(
            NAME nanokmc_repository_contract
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/validation/check_repository.py)
        add_test(
            NAME nanokmc_solver_smoke
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/validation/run_smoke.py
                    --exe $<TARGET_FILE:nanokmc>)
        add_test(
            NAME nanokmc_output_smoke
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/validation/run_output_smoke.py
                    --exe $<TARGET_FILE:nanokmc>)
        add_test(
            NAME nanokmc_profile_axis_smoke
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/validation/run_profile_axis_smoke.py
                    --exe $<TARGET_FILE:nanokmc>)
        add_test(
            NAME nanokmc_rate_category_limits
            COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/validation/run_rate_category_limits.py
                    --exe $<TARGET_FILE:nanokmc>)
    endif()
endif()

include(GNUInstallDirs)
install(TARGETS nanokmc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
