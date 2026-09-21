#include "SystemKMCActiveFilteredGeneric.h"
SystemKMCActiveFilteredGeneric::SystemKMCActiveFilteredGeneric() {
}

SystemKMCActiveFilteredGeneric::~SystemKMCActiveFilteredGeneric() {
}

double SystemKMCActiveFilteredGeneric::AcceptanceProbabilityFromEnvironment(
    const LocalPairEnvironment& env) {
    const int a = env.firstSpecies;
    const int b = env.secondSpecies;
    if (a == b) {
        std::cerr << "ActiveFilteredGeneric: inactive pair passed to evaluator"
                  << std::endl;
        std::exit(2);
    }

    // Symmetric homogeneous NN benchmark Hamiltonian:
    //   E_bond = -Ea/2 for equal species, 0 for unlike species.
    // The mutual selected bond is unlike both before and after the exchange,
    // so it cancels. The explicit ordered environment is nevertheless retained
    // in full, allowing future derived systems to use directional/motif terms.
    const int firstMutualDirection = env.selectedDirection;
    const int secondMutualDirection = OppositeDirection(firstMutualDirection);

    int sameBefore = 0;
    int sameAfter = 0;
    for (int d = 0; d < CoordinationNumber; ++d) {
        if (d != firstMutualDirection) {
            const int s = env.firstNeighborSpecies[d];
            sameBefore += (s == a) ? 1 : 0;
            sameAfter  += (s == b) ? 1 : 0;
        }
        if (d != secondMutualDirection) {
            const int s = env.secondNeighborSpecies[d];
            sameBefore += (s == b) ? 1 : 0;
            sameAfter  += (s == a) ? 1 : 0;
        }
    }

    const double before = -0.5 * Ea * static_cast<double>(sameBefore);
    const double after  = -0.5 * Ea * static_cast<double>(sameAfter);
    const double deltaE = after - before;
    if (deltaE <= 0.0) return 1.0;
    return std::exp(-deltaE / kT);
}

double SystemKMCActiveFilteredGeneric::PairAcceptanceProbability(
    std::uint32_t firstSite, std::uint32_t secondSite,
    int canonicalDirection) {
    LocalPairEnvironment env;
    GatherLocalPairEnvironment(firstSite, secondSite, canonicalDirection, env);
    return AcceptanceProbabilityFromEnvironment(env);
}
