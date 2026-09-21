#pragma once

#include "../SystemKMCHomogenous.h"

#include <array>
#include <vector>
#include <limits>
#include <cmath>

// Internal utility base shared by the binary rate-informed reference solvers.
// It defines the fixed FCC topology and the eight exact Metropolis factors for
// the bundled symmetric binary nearest-neighbour benchmark, but it does not
// prescribe how those rates are stored or selected.
class SystemKMCBinaryRateBase : public SystemKMCHomogenous {

    public:
        SystemKMCBinaryRateBase();
        ~SystemKMCBinaryRateBase();

    protected:
        void InitSubSys() override;
        void CodeCalculation() override;

        static constexpr int RateClassCount = 8;
        static constexpr signed char InactiveClass = -1;
        std::array<double, RateClassCount> classProbability;

        std::vector<int> fixedX;
        std::vector<int> fixedY;
        std::vector<int> fixedZ;
        std::vector<unsigned long long> fixedSiteAtCoord;
        unsigned long long invalidFixedSite;

        void BuildFixedSiteMap();
        unsigned long long CoordLinearIndex(int x, int y, int z) const;
        unsigned long long FixedSiteIndex(int x, int y, int z) const;
        void NeighborCoords(int x, int y, int z, int direction,
                            int& xn, int& yn, int& zn) const;
        unsigned long long CanonicalBondKeyFromIncident(
            unsigned long long fixedSite, int direction);
        void ResolveCanonicalBond(unsigned long long key,
                                  int& x, int& y, int& z,
                                  int& xn, int& yn, int& zn) const;
        int CountSpeciesZeroNeighborsExcept(int x, int y, int z,
                                            int ex, int ey, int ez);
        signed char RateClassForKey(unsigned long long key);
};
