#pragma once

#include "internal/SystemKMCPartialFilterCore.h"
#include <vector>
#include <cstdint>
#include <limits>
#include <cmath>

// Optimized partial structural-filter reference solver used by the paper
// benchmark. Its candidate population follows the intermediate binary-alloy
// construction described in Appendix A of Bortz, Kalos & Lebowitz (1975):
// maintain A sites with at least one B neighbour, select one eligible A, then
// select a neighbour and retain residual structural/Metropolis rejection.
// DOI: 10.1016/0021-9991(75)90060-1. See docs/literature.md.
//
// The internal Core class supplies the validated eligible-site mechanics; this
// public solver adds cached local counts and avoids unused active-table work.
// This is a controlled literature/reference architecture, not a NanoKMC method.
class SystemKMCPartialFilterOptimized : public SystemKMCPartialFilterCore {

    public:
        SystemKMCPartialFilterOptimized();
        ~SystemKMCPartialFilterOptimized();

    protected:
        void InitSubSys() override;
        void DetermineBonds() override;
        void JumpAttempt() override;
        void EvalSysStep(unsigned long long n) override;

    private:
        // Species-1 nearest-neighbour count on the fixed lattice coordinate
        // grid.  Only FCC sites are queried.  This simultaneously provides an
        // O(1) eligibility test and O(1) Metropolis rate-class evaluation.
        std::vector<unsigned char> speciesOneNeighborCount;

        // Stamp-based de-duplication for the <=20 atom identities whose
        // eligibility may change after an accepted exchange.
        std::vector<std::uint32_t> atomTouchStamp;
        std::uint32_t atomStampCounter;
        std::vector<unsigned long long> affectedAtomFastScratch;

        unsigned long long CoordLinearIndexFast(int x, int y, int z) const;
        unsigned char SpeciesOneCountAt(int x, int y, int z) const;
        void SetSpeciesOneCountAt(int x, int y, int z, int value);
        void InitializeSpeciesOneNeighborCounts();
        double FastPairMetropolisProbability(int x, int y, int z,
                                             int xn, int yn, int zn);
        void UpdateSpeciesOneNeighborCountsAfterExchange(int x, int y, int z,
                                                         int xn, int yn, int zn,
                                                         int s, int sn);
        void CollectAffectedAtomsFast(int x, int y, int z,
                                      int xn, int yn, int zn);
        void RefreshEligibilityFast(unsigned long long atomIndex);
        void UpdateBondNumberAnalytically(unsigned char bAtA,
                                          unsigned char bAtB);
        void ValidateOptimizedStateOrAbort(const char* context);
        void ValidateFastProbabilitySampleOrAbort(const char* context);
};
