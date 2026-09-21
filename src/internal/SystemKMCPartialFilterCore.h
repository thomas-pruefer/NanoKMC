#pragma once

#include "../SystemKMCHomogenous.h"

#include <array>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdint>

class SystemKMCPartialFilterCore : public SystemKMCHomogenous {

    public:
        SystemKMCPartialFilterCore();
        ~SystemKMCPartialFilterCore();

    protected:
        void InitSubSys() override;
        void CodeCalculation() override;
        void DetermineBonds() override;
        void JumpAttempt() override = 0;
        double MCSIncrementPerAttempt() const override;
        bool RecomputeMCSIncrementEveryAttempt() const override;
        void EvalSysStep(unsigned long long n) override;

    protected:
        static constexpr int RateClassCount = 8;

        // BKL (1975), Appendix-A style population: species-0 (A) atoms that currently
        // have at least one species-1 (B) nearest neighbour.  Atom indices are
        // stable identities in NanoKMC because ExchangeSites moves xpr/ypr/zpr with
        // the exchanged atom.
        std::vector<unsigned long long> eligibleA;
        std::vector<unsigned long long> eligibleRow;
        unsigned long long invalidRow;

        // Same finite FCC Metropolis factors as the active and exact-class
        // homogeneous benchmark solvers.  This avoids reconstructing the
        // 18-site bit pattern for a proposal that has already survived the
        // partial structural filter.
        std::array<double, RateClassCount> classProbability;

        // Precomputed mean common-MCS increment for the current eligible-site
        // population.  RunMC asks every attempt because the denominator is not
        // BondNumber, but this value itself is updated only when M_A changes.
        double partialClockIncrement;

        void NeighborCoords(int x, int y, int z, int direction,
                            int& xn, int& yn, int& zn) const;
        bool HasSpeciesOneNeighbor(unsigned long long atomIndex);
        bool IsEligibleA(unsigned long long atomIndex);
        void AddEligible(unsigned long long atomIndex);
        void RemoveEligible(unsigned long long atomIndex);
        void RefreshEligibility(unsigned long long atomIndex);
        void UpdateClockIncrement();
        void ValidateEligiblePopulationOrAbort(const char* context);

        int CountSpeciesZeroNeighborsExcept(int x, int y, int z,
                                            int ex, int ey, int ez);
        double PairMetropolisProbability(int x, int y, int z,
                                         int xn, int yn, int zn);
};
