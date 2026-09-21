#pragma once

#include "SystemKMCBinaryRateBase.h"
#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>

// Internal mechanics for the exact finite-rate-class reference solver. Exact
// energetic classes and rejection-free class/member selection are established
// prior art; see docs/literature.md for BKL, Sadiq and Schulze references.
class SystemKMCExactClassCore : public SystemKMCBinaryRateBase {

    public:
        SystemKMCExactClassCore();
        ~SystemKMCExactClassCore();

    protected:
        void InitSubSys() override;
        void DetermineBonds() override;
        void JumpAttempt() override;
        double MCSIncrementPerAttempt() const override;
        bool RecomputeMCSIncrementEveryAttempt() const override;
        void EvalSysStep(unsigned long long n) override;

        std::array<std::vector<unsigned long long>, RateClassCount> classBonds;
        std::vector<signed char> bondClass;
        std::vector<unsigned long long> bondRow;
        double totalRateWeight;

        void AddBondKey(unsigned long long key, signed char rateClass);
        void RemoveBondKey(unsigned long long key);
        void CollectAffectedBondKeys(int x, int y, int z,
                                     int xn, int yn, int zn,
                                     std::vector<unsigned long long>& keys);
        void RecalculateTotalRateWeight();
        void ValidateClassPopulationOrAbort(const char* context) const;
};
