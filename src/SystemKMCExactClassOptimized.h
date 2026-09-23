#pragma once

#include "internal/SystemKMCExactClassCore.h"
#include <array>
#include <vector>
#include <cstdint>
#include <limits>
#include <cmath>

// Optimized exact finite-rate-class reference solver for the bundled
// binary FCC exchange model. The design belongs to the established rejection-free/rate-class
// lineage: BKL/n-fold (1975), active-bond exact classes in Sadiq (1984), and
// finite-class list implementations such as Schulze (2002).
// Sadiq DOI: 10.1016/0021-9991(84)90028-7
// Schulze DOI: 10.1103/PhysRevE.65.036704
// See docs/literature.md for the claim/provenance boundary.
//
// The internal Core class supplies the validated class-list mechanics; this
// public solver adds cached local counts and differential class refresh. This is
// a controlled literature/reference architecture, not a NanoKMC method.
class SystemKMCExactClassOptimized : public SystemKMCExactClassCore {

    public:
        SystemKMCExactClassOptimized();
        ~SystemKMCExactClassOptimized();

    protected:
        void InitSubSys() override;
        void DetermineBonds() override;
        void JumpAttempt() override;
        double MCSIncrementPerAttempt() const override;
        bool RecomputeMCSIncrementEveryAttempt() const override;
        void EvalSysStep(unsigned long long n) override;

    private:
        // Number of species-0 nearest neighbours at each fixed FCC site.
        // This reduces exact rate classification from two 12-neighbour scans
        // per bond to O(1) arithmetic.
        std::vector<unsigned char> zeroNeighborCount;

        // Sparse-delta scratch for the <=20 sites touched by an A<->B exchange.
        std::vector<signed char> zeroNeighborDelta;
        std::vector<std::uint32_t> siteTouchStamp;
        std::uint32_t siteStampCounter;
        std::vector<unsigned long long> touchedSiteScratch;
        std::vector<unsigned long long> changedSiteScratch;

        // Stamp-based bond de-duplication avoids allocating/sorting ~O(10^2)
        // descriptors on every exact-class event.
        std::vector<std::uint32_t> bondTouchStamp;
        std::uint32_t bondStampCounter;
        std::vector<unsigned long long> affectedBondScratch;

        // Current 6/Q clock value.  The base RunMC still requests it before
        // each event, but the hot call becomes a load rather than a division.
        double optimizedClockIncrement;
        unsigned long long optimizedEventCount;

        void InitializeZeroNeighborCounts();
        signed char FastRateClassForKey(unsigned long long key);
        void BeginSiteDeltaTransaction();
        void AccumulateZeroNeighborDelta(unsigned long long site, int delta);
        void BuildAffectedSitesAndBonds(int x, int y, int z,
                                        int xn, int yn, int zn,
                                        int s, int sn);
        void ApplyZeroNeighborDeltas();
        void RefreshBondKeyDifferential(unsigned long long key);
        void UpdateOptimizedClock();
        void ValidateOptimizedStateOrAbort(const char* context);
};
