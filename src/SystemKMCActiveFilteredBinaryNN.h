#pragma once

#include "SystemKMCActiveFilteredBase.h"

#include <array>
#include <vector>
#include <cstdint>
#include <cmath>

// Binary nearest-neighbour optimized backend for NanoKMC Active-Filtered.
//
// This class uses exactly the same structurally filtered candidate population,
// uniform selection and common-MCS clock as ActiveFilteredGeneric. It exploits
// only the symmetric two-species nearest-neighbour benchmark Hamiltonian:
// one cached integer per lattice site (the number of species-1 neighbours) is
// sufficient to determine the selected pair's Metropolis factor in O(1).
//
// No energetic information is used for event preselection. The optimization is
// confined to post-selection acceptance evaluation and local cache maintenance.
class SystemKMCActiveFilteredBinaryNN : public SystemKMCActiveFilteredBase {

    public:
        SystemKMCActiveFilteredBinaryNN();
        ~SystemKMCActiveFilteredBinaryNN() override;

    protected:
        const char* AcceptanceBackendName() const override {
            return "KMCActiveFilteredBinaryNN";
        }

        double PairAcceptanceProbability(std::uint32_t firstSite,
                                         std::uint32_t secondSite,
                                         int canonicalDirection) override;
        void InitializeAcceptanceState() override;
        void UpdateAcceptanceStateAfterExchange(std::uint32_t firstSite,
                                                std::uint32_t secondSite,
                                                int oldFirstSpecies,
                                                int oldSecondSpecies) override;
        void ValidateAcceptanceStateOrAbort(const char* context) override;

    private:
        static constexpr int RateClassCount = 8;
        std::vector<unsigned char> speciesOneNeighborCount;
        std::array<double, RateClassCount> classProbability{};

        void InitializeNeighborCounts();
        void UpdateNeighborCountsAfterExchange(std::uint32_t firstSite,
                                               std::uint32_t secondSite,
                                               int oldFirstSpecies,
                                               int oldSecondSpecies);
        double ReferenceProbabilityFromExplicitCounts(
            std::uint32_t firstSite,
            std::uint32_t secondSite) const;
};
