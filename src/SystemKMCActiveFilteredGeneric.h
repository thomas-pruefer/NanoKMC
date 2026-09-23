#pragma once

#include "SystemKMCActiveFilteredBase.h"

// General Active-Filtered local-environment backend.
//
// Candidate filtering, selection, timing and active-list updates are inherited
// unchanged from SystemKMCActiveFilteredBase. This backend reconstructs the
// ordered nearest-neighbour environment only after a candidate has been
// selected. The bundled default Hamiltonian is a symmetric homogeneous model:
// an equal-species NN bond contributes -Ea/2 and an unlike bond contributes 0.
// The local-environment representation preserves arbitrary species labels; the
// bundled initializer/examples are binary, while derived models can provide
// true multi-species initialization and acceptance physics.
//
// More complex local Hamiltonians can derive from this class and override
// AcceptanceProbabilityFromEnvironment() without changing the Active-Filtered
// structural algorithm. The LocalPairEnvironment retains directional species
// information rather than compressing it to coordination counts.
class SystemKMCActiveFilteredGeneric : public SystemKMCActiveFilteredBase {

    public:
        SystemKMCActiveFilteredGeneric();
        ~SystemKMCActiveFilteredGeneric() override;

    protected:
        const char* AcceptanceBackendName() const override {
            return "KMCActiveFilteredGeneric";
        }

        double PairAcceptanceProbability(std::uint32_t firstSite,
                                         std::uint32_t secondSite,
                                         int canonicalDirection) override;

        virtual double AcceptanceProbabilityFromEnvironment(
            const LocalPairEnvironment& env);
};
