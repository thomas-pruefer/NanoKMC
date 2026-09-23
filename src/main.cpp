#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

#include "CalcClass.h"
#include "SystemClass.h"
#include "SystemKMCActiveFilteredBinaryNN.h"
#include "SystemKMCActiveFilteredGeneric.h"
#include "SystemKMCClassical.h"
#include "SystemKMCExactClassOptimized.h"
#include "SystemKMCPartialFilterOptimized.h"
#include "SystemKMCRateCategoryOptimized.h"

namespace {

constexpr const char* kSolverIds[] = {
    "KMCClassical",
    "KMCActiveFilteredGeneric",
    "KMCActiveFilteredBinaryNN",
    "KMCPartialFilterOptimized",
    "KMCRateCategoryOptimized",
    "KMCExactClassOptimized",
};

void PrintSupportedSolvers(std::ostream& out) {
    for (const char* id : kSolverIds) {
        out << "  " << id << '\n';
    }
}

std::unique_ptr<SystemClass> CreateSystem(const std::string& systemId) {
    if (systemId == "KMCClassical")
        return std::make_unique<SystemKMCClassical>();
    if (systemId == "KMCActiveFilteredGeneric")
        return std::make_unique<SystemKMCActiveFilteredGeneric>();
    if (systemId == "KMCActiveFilteredBinaryNN")
        return std::make_unique<SystemKMCActiveFilteredBinaryNN>();
    if (systemId == "KMCPartialFilterOptimized")
        return std::make_unique<SystemKMCPartialFilterOptimized>();
    if (systemId == "KMCRateCategoryOptimized")
        return std::make_unique<SystemKMCRateCategoryOptimized>();
    if (systemId == "KMCExactClassOptimized")
        return std::make_unique<SystemKMCExactClassOptimized>();

    std::cerr << "Unsupported SystemID: " << systemId << "\n"
              << "Supported SystemID values:\n";
    PrintSupportedSolvers(std::cerr);
    return nullptr;
}

int ChangeDirectory(const char* path) {
#ifdef _WIN32
    return _chdir(path);
#else
    return chdir(path);
#endif
}

void PrintVersion() {
    std::cout << "NanoKMC " << NANOKMC_VERSION;
#ifdef NANOKMC_PAPER_BUILD_ENABLED
    std::cout << " (manuscript-benchmark compiler policy)";
#else
    std::cout << " (standard compiler policy)";
#endif
    std::cout << '\n';
}

void PrintHelp() {
    PrintVersion();
    std::cout << "Usage: nanokmc [run-directory]\n\n"
              << "NanoKMC reads nanokmc.in and output/CalcData.csv from the "
                 "selected run directory (or the current directory when no "
                 "directory is supplied).\n\n"
              << "Options:\n"
              << "  -h, --help          Show this help text\n"
              << "  --version           Show the software version\n"
              << "  --list-solvers      List supported SystemID values\n\n"
              << "Supported solvers:\n";
    PrintSupportedSolvers(std::cout);
    std::cout << "\nKMCRateCategoryOptimized: RateCategoryCount=1|2|4|8 "
                 "(default validated setting: 4).\n\n"
              << "Simulation / morphology EvalParam tokens:\n"
              << "  Benchmark\n"
              << "  ClusterDistribution\n"
              << "  AxialCompositionProfile\n"
              << "  CylindricalCompositionProfile\n"
              << "  SphericalCompositionProfile\n"
              << "  ProjectedCompositionXZ\n\n"
              << "State / visualization EvalParam tokens:\n"
              << "  Rasmol\n"
              << "  OnefileRasmol\n"
              << "  BlenderSimple\n"
              << "  WriteCSV\n\n"
              << "Profile axes: AxialCompositionProfileAxis=X|Y|Z and "
                 "CylindricalCompositionProfileAxis=X|Y|Z.\n";
}

bool CheckRunDirectory() {
    namespace fs = std::filesystem;
    bool ok = true;
    if (!fs::is_regular_file("nanokmc.in")) {
        std::cerr << "Missing required input file: nanokmc.in\n";
        ok = false;
    }
    if (!fs::is_regular_file(fs::path("output") / "CalcData.csv")) {
        std::cerr << "Missing required checkpoint file: output/CalcData.csv\n";
        ok = false;
    }
    return ok;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc > 2) {
        std::cerr << "Too many arguments.\n";
        PrintHelp();
        return 2;
    }

    if (argc == 2) {
        const std::string arg = argv[1];
        if (arg == "--help" || arg == "-h") {
            PrintHelp();
            return 0;
        }
        if (arg == "--version") {
            PrintVersion();
            return 0;
        }
        if (arg == "--list-solvers") {
            PrintSupportedSolvers(std::cout);
            return 0;
        }
        if (ChangeDirectory(argv[1]) != 0) {
            std::cerr << "Could not change into run directory: " << argv[1]
                      << '\n';
            return 1;
        }
    }

    if (!CheckRunDirectory()) {
        std::cerr << "See docs/quickstart.md and docs/input-output.md for the required "
                     "run-directory layout.\n";
        return 2;
    }

    CalcClass calc;
    calc.Init();

    auto system = CreateSystem(calc.SystemID);
    if (!system) return 2;

    system->Init(calc.GetStepData(0), calc.ColumnParam);

    for (int i = 0; i < calc.RecordNumber; ++i) {
        if (calc.Recorded[i] == 0) {
            calc.CalcData[i][2] = i;
            if (i > 0 && calc.Recorded[i - 1] == 1) {
                system = CreateSystem(calc.SystemID);
                if (!system) return 2;
                system->Init(calc.GetStepData(i - 1), calc.ColumnParam);
            }
            system->RunMC(calc.GetStepData(i));
        } else if (i > 0) {
            system = CreateSystem(calc.SystemID);
            if (!system) return 2;
            system->Init(calc.GetStepData(i), calc.ColumnParam);
        }
        calc.AddCalcData(system->GetStepData(), i);
        system->RunEval(calc.GetStepData(i));
    }

    calc.WriteCalcData();
    system->Zip();
    return 0;
}
