#pragma once

#include "SystemKMCHomogenous.h"

#include <array>
#include <vector>
#include <cstdint>
#include <limits>
#include <cmath>
#include <algorithm>

// Shared structural engine for the NanoKMC Active-Filtered method.
//
// The method implemented here is deliberately independent of the local
// Hamiltonian:
//   * maintain the complete population of structurally active unlike
//     nearest-neighbour bonds;
//   * select uniformly from that filtered population;
//   * evaluate a post-selection acceptance probability supplied by a derived
//     local-physics backend;
//   * advance the FCC common-MCS clock by 6/B for every selected active bond;
//   * refresh only candidate bonds incident on the exchanged sites.
//
// The present repository instantiates FCC topology (z=12, six canonical
// undirected bonds/site). The separation between structural engine and
// acceptance backend is intentional: future lattice/topology layers can reuse
// the same method without changing the local-physics interface.
//
// Scientific boundary: active/event lists, reverse addressing and local list
// maintenance are established prior art. The NanoKMC-specific design choice is
// to maintain the complete *structural* unlike-bond population while leaving
// energetics as a post-selection Metropolis test. See docs/literature.md,
// especially BKL (1975), Sadiq (1984), Shida & Henriques (1997), and Schulze
// (2002, 2008).
class SystemKMCActiveFilteredBase : public SystemKMCHomogenous {

    public:
        SystemKMCActiveFilteredBase();
        ~SystemKMCActiveFilteredBase() override;

    protected:
        static constexpr int CoordinationNumber = 12;
        static constexpr int CanonicalDirectionCount = 6;
        static constexpr std::uint32_t InvalidRow32 =
            std::numeric_limits<std::uint32_t>::max();
        static constexpr std::uint32_t InvalidSite32 =
            std::numeric_limits<std::uint32_t>::max();

        struct LocalPairEnvironment {
            std::uint32_t firstSite = 0;
            std::uint32_t secondSite = 0;
            int selectedDirection = 0;
            int firstSpecies = 0;
            int secondSpecies = 0;
            std::array<std::uint32_t, CoordinationNumber> firstNeighborSite{};
            std::array<std::uint32_t, CoordinationNumber> secondNeighborSite{};
            std::array<int, CoordinationNumber> firstNeighborSpecies{};
            std::array<int, CoordinationNumber> secondNeighborSpecies{};
        };

        void InitSubSys() override;
        void DetermineBonds() override;
        void JumpAttempt() override;
        void EvalSysStep(unsigned long long n) override;
        double MCSIncrementPerAttempt() const override;

        // Acceptance-backend interface. Candidate selection remains purely
        // structural in all Active-Filtered implementations.
        virtual const char* AcceptanceBackendName() const = 0;
        virtual double PairAcceptanceProbability(std::uint32_t firstSite,
                                                 std::uint32_t secondSite,
                                                 int canonicalDirection) = 0;

        // Optional derived-cache lifecycle hooks.
        virtual void InitializeAcceptanceState();
        virtual void UpdateAcceptanceStateAfterExchange(std::uint32_t firstSite,
                                                        std::uint32_t secondSite,
                                                        int oldFirstSpecies,
                                                        int oldSecondSpecies);
        virtual void ValidateAcceptanceStateOrAbort(const char* context);

        // Fixed FCC topology and site state. Site indices do not move when
        // atoms exchange positions.
        std::vector<int> fixedX;
        std::vector<int> fixedY;
        std::vector<int> fixedZ;
        std::vector<std::uint32_t> coordToSite;
        std::vector<std::uint32_t> neighborSite; // [site*12 + direction]
        std::vector<unsigned char> siteSpecies;

        // Dense canonical active-bond list and reverse addressing.
        std::vector<std::uint32_t> activeKeys;
        std::vector<std::uint32_t> activeRow;

        // Stamp-based de-duplication for the incident-bond union.
        std::vector<std::uint32_t> bondTouchStamp;
        std::uint32_t bondStampCounter;
        std::vector<std::uint32_t> affectedBondScratch;

        unsigned long long CoordLinearIndex(int x, int y, int z) const;
        void NeighborCoords(int x, int y, int z, int direction,
                            int& xn, int& yn, int& zn) const;
        void BuildFixedTopology();
        void InitializeSiteSpeciesCache();
        void InitializeFastActiveList();

        std::uint32_t CanonicalBondKeyFromIncident(std::uint32_t site,
                                                   int direction) const;
        void CollectAffectedBondKeys(std::uint32_t firstSite,
                                     std::uint32_t secondSite);
        void AddActiveKey(std::uint32_t key);
        void RemoveActiveKey(std::uint32_t key);
        void RefreshActiveKey(std::uint32_t key);

        void GatherLocalPairEnvironment(std::uint32_t firstSite,
                                        std::uint32_t secondSite,
                                        int canonicalDirection,
                                        LocalPairEnvironment& env) const;

        int OppositeDirection(int direction) const;
        void ValidateBaseStateOrAbort(const char* context);
};
