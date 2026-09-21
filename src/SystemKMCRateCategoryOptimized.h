#pragma once

#include "internal/SystemKMCBinaryRateBase.h"
#include <array>
#include <vector>
#include <cstdint>
#include <cmath>
#include <limits>

// Coarse rate-category rejection reference solver.
//
// Scientific role:
//   active unlike bonds -> coarse rate category -> category/member selection
//   -> residual acceptance using exact_rate/category_upper_bound.
//
// This follows the inverted-list / rate-envelope design lineage discussed by
// Schulze (2008) and demonstrated empirically by Saum, Schulze & Ratsch (2009).
// It is a controlled in-code literature/reference architecture, not a NanoKMC
// Active-Filtered method.  See docs/literature.md.
//
// For the bundled binary FCC NN benchmark there are eight exact Metropolis
// classes c=0..7.  RateCategoryCount groups those classes contiguously:
//   M=1: {0..7}                    (single majorant population)
//   M=2: {0..3}, {4..7}
//   M=4: {0,1}, {2,3}, {4,5}, {6,7}   [release default]
//   M=8: one exact class per category (rejection-free limiting case)
class SystemKMCRateCategoryOptimized : public SystemKMCBinaryRateBase {

    public:
        SystemKMCRateCategoryOptimized();
        ~SystemKMCRateCategoryOptimized();

    protected:
        void InitSubSys() override;
        void DetermineBonds() override;
        void JumpAttempt() override;
        double MCSIncrementPerAttempt() const override;
        bool RecomputeMCSIncrementEveryAttempt() const override;
        void EvalSysStep(unsigned long long n) override;

    private:
        static constexpr int MaxCategoryCount = RateClassCount;
        static constexpr signed char InactiveCategory = -1;

        int rateCategoryCount;
        int rateCategoryWidth;
        std::array<std::vector<unsigned long long>, MaxCategoryCount> categoryBonds;
        std::array<double, MaxCategoryCount> categoryUpperProbability;
        std::vector<signed char> bondCategory;
        std::vector<unsigned long long> categoryRow;
        double totalMajorantWeight;
        double categoryClockIncrement;
        unsigned long long categoryEventCount;

        // Same O(1) binary-NN classifier/cache strategy as the exact-class
        // reference solver.  Only the coarse category is persisted per bond;
        // the exact class is evaluated for the selected candidate to perform
        // the residual rejection test.
        std::vector<unsigned char> zeroNeighborCount;
        std::vector<signed char> zeroNeighborDelta;
        std::vector<std::uint32_t> siteTouchStamp;
        std::uint32_t siteStampCounter;
        std::vector<unsigned long long> touchedSiteScratch;
        std::vector<unsigned long long> changedSiteScratch;

        std::vector<std::uint32_t> bondTouchStamp;
        std::uint32_t bondStampCounter;
        std::vector<unsigned long long> affectedBondScratch;

        int CategoryForRateClass(int rateClass) const;
        void InitializeCategoryBounds();
        void InitializeZeroNeighborCounts();
        signed char FastRateClassForKey(unsigned long long key);
        void AddCategoryBond(unsigned long long key, signed char category);
        void RemoveCategoryBond(unsigned long long key);
        void RecalculateMajorantWeight();
        void UpdateClockIncrement();

        void BeginSiteDeltaTransaction();
        void AccumulateZeroNeighborDelta(unsigned long long site, int delta);
        void BuildAffectedSitesAndBonds(int x, int y, int z,
                                        int xn, int yn, int zn,
                                        int s, int sn);
        void ApplyZeroNeighborDeltas();
        void RefreshBondKeyDifferential(unsigned long long key);
        void ValidateOptimizedStateOrAbort(const char* context);
};
