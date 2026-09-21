#pragma once

#include "SystemKMCHomogenous.h"
#include <array>
#include <cmath>

// Classical blind-proposal Kawasaki solver for the binary FCC benchmark.
// Each attempt selects one lattice site uniformly and then one of its 12
// nearest neighbours uniformly. Same-species proposals are structural nulls;
// unlike pairs retain the same Metropolis acceptance rule used by the other
// public benchmark solvers. This class is the deliberately unfiltered reference
// against which the proposal-frequency mapping of the filtered solvers is
// defined; it is not presented as a new Monte Carlo algorithm.
class SystemKMCClassical : public SystemKMCHomogenous {

    public:
        SystemKMCClassical();
        ~SystemKMCClassical();

    protected:
        void InitSubSys() override;
        void CodeCalculation() override;
        void DetermineBonds() override;
        void JumpAttempt() override;
        double MCSIncrementPerAttempt() const override;

    private:
        static constexpr int CoordinationNumber = 12;
        static constexpr int RateClassCount = 8;
        std::array<double, RateClassCount> classProbability;

        void NeighborCoords(int x, int y, int z, int direction,
                            int& xn, int& yn, int& zn) const;
        int CountSpeciesOneNeighbors(int x, int y, int z);
        double PairMetropolisProbability(int x, int y, int z,
                                         int xn, int yn, int zn,
                                         int s, int sn,
                                         int& bAtA, int& bAtB);
        void UpdateBondNumberAfterAcceptedExchange(int bAtA, int bAtB);
};
