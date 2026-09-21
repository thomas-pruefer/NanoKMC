#include "SystemKMCBinaryRateBase.h"
SystemKMCBinaryRateBase::SystemKMCBinaryRateBase()
    : invalidFixedSite(std::numeric_limits<unsigned long long>::max()) {
    classProbability.fill(0.0);
}

SystemKMCBinaryRateBase::~SystemKMCBinaryRateBase() {
}

void SystemKMCBinaryRateBase::InitSubSys() {
    SystemKMCHomogenous::InitSubSys();
    if (NSpecies != 2) {
        std::cerr << "Binary rate-informed reference solvers require exactly two species."
                  << std::endl;
        std::exit(2);
    }
}

void SystemKMCBinaryRateBase::CodeCalculation() {
    SystemClass::CodeCalculation();
    classProbability[0] = 1.0;
    for (int c = 1; c < RateClassCount; ++c) {
        classProbability[c] = std::exp(-Ea / kT * static_cast<double>(c));
    }
}

unsigned long long SystemKMCBinaryRateBase::CoordLinearIndex(
    int x, int y, int z) const {
    return (static_cast<unsigned long long>(x) * static_cast<unsigned long long>(ly1)
          + static_cast<unsigned long long>(y)) * static_cast<unsigned long long>(lz1)
          + static_cast<unsigned long long>(z);
}

unsigned long long SystemKMCBinaryRateBase::FixedSiteIndex(
    int x, int y, int z) const {
    return fixedSiteAtCoord[CoordLinearIndex(x, y, z)];
}

void SystemKMCBinaryRateBase::BuildFixedSiteMap() {
    fixedX.clear(); fixedY.clear(); fixedZ.clear();
    fixedX.reserve(TotalAtoms); fixedY.reserve(TotalAtoms); fixedZ.reserve(TotalAtoms);

    const unsigned long long fullGridSize =
        static_cast<unsigned long long>(lx1) *
        static_cast<unsigned long long>(ly1) *
        static_cast<unsigned long long>(lz1);
    fixedSiteAtCoord.assign(fullGridSize, invalidFixedSite);

    unsigned long long site = 0;
    for (int x = 0; x < lx1; ++x) {
        for (int y = 0; y < ly1; ++y) {
            for (int z = 0; z < lz1; ++z) {
                if (((x ^ y ^ z) & 1) == 0) {
                    fixedX.push_back(x); fixedY.push_back(y); fixedZ.push_back(z);
                    fixedSiteAtCoord[CoordLinearIndex(x, y, z)] = site;
                    ++site;
                }
            }
        }
    }
    if (site != static_cast<unsigned long long>(TotalAtoms)) {
        std::cerr << "SystemKMCBinaryRateBase: fixed-site map size mismatch: "
                  << site << " vs TotalAtoms=" << TotalAtoms << std::endl;
        std::exit(2);
    }
}

void SystemKMCBinaryRateBase::NeighborCoords(
    int x, int y, int z, int direction,
    int& xn, int& yn, int& zn) const {

    const int xp1 = ((x + 1) & lx);
    const int yp1 = ((y + 1) & ly);
    const int zp1 = ((z + 1) & lz);
    const int xm1 = ((x + lx) & lx);
    const int ym1 = ((y + ly) & ly);
    const int zm1 = ((z + lz) & lz);

    switch (direction) {
        case  0: xn=xp1; yn=  y; zn=zp1; break;
        case  1: xn=  x; yn=ym1; zn=zp1; break;
        case  2: xn=xm1; yn=  y; zn=zp1; break;
        case  3: xn=  x; yn=yp1; zn=zp1; break;
        case  4: xn=xp1; yn=yp1; zn=  z; break;
        case  5: xn=xm1; yn=yp1; zn=  z; break;
        case  6: xn=  x; yn=yp1; zn=zm1; break;
        case  7: xn=xp1; yn=  y; zn=zm1; break;
        case  8: xn=xm1; yn=  y; zn=zm1; break;
        case  9: xn=  x; yn=ym1; zn=zm1; break;
        case 10: xn=xp1; yn=ym1; zn=  z; break;
        case 11: xn=xm1; yn=ym1; zn=  z; break;
        default: xn=x; yn=y; zn=z; break;
    }
}

unsigned long long SystemKMCBinaryRateBase::CanonicalBondKeyFromIncident(
    unsigned long long fixedSite, int direction) {

    if (direction < 6) {
        return fixedSite * 6ULL + static_cast<unsigned long long>(direction);
    }

    int xn, yn, zn;
    NeighborCoords(fixedX[fixedSite], fixedY[fixedSite], fixedZ[fixedSite],
                   direction, xn, yn, zn);
    const unsigned long long other = FixedSiteIndex(xn, yn, zn);
    if (other == invalidFixedSite) {
        std::cerr << "SystemKMCBinaryRateBase: invalid FCC neighbor in canonicalization"
                  << std::endl;
        std::exit(2);
    }
    const int canonicalDirection = xyz2nn(direction);
    return other * 6ULL + static_cast<unsigned long long>(canonicalDirection);
}

void SystemKMCBinaryRateBase::ResolveCanonicalBond(
    unsigned long long key,
    int& x, int& y, int& z,
    int& xn, int& yn, int& zn) const {

    const unsigned long long site = key / 6ULL;
    const int direction = static_cast<int>(key % 6ULL);
    x = fixedX[site]; y = fixedY[site]; z = fixedZ[site];
    NeighborCoords(x, y, z, direction, xn, yn, zn);
}

int SystemKMCBinaryRateBase::CountSpeciesZeroNeighborsExcept(
    int x, int y, int z, int ex, int ey, int ez) {

    int count = 0;
    for (int d = 0; d < 12; ++d) {
        int xn, yn, zn;
        NeighborCoords(x, y, z, d, xn, yn, zn);
        if (xn == ex && yn == ey && zn == ez) continue;
        if ((this->*GetSpecies)(xn, yn, zn) == 0) ++count;
    }
    return count;
}

signed char SystemKMCBinaryRateBase::RateClassForKey(unsigned long long key) {
    int x, y, z, xn, yn, zn;
    ResolveCanonicalBond(key, x, y, z, xn, yn, zn);

    const int s = (this->*GetSpecies)(x, y, z);
    const int sn = (this->*GetSpecies)(xn, yn, zn);
    if (s == sn) return InactiveClass;

    int x0, y0, z0, x1, y1, z1;
    if (s == 0 && sn == 1) {
        x0=x; y0=y; z0=z; x1=xn; y1=yn; z1=zn;
    } else if (s == 1 && sn == 0) {
        x0=xn; y0=yn; z0=zn; x1=x; y1=y; z1=z;
    } else {
        std::cerr << "SystemKMCBinaryRateBase requires binary species 0/1."
                  << std::endl;
        std::exit(2);
    }

    const int nini = CountSpeciesZeroNeighborsExcept(x0, y0, z0, x1, y1, z1);
    const int nfin = CountSpeciesZeroNeighborsExcept(x1, y1, z1, x0, y0, z0);
    const int delta = nini - nfin;
    const int rateClass = (delta > 0) ? delta : 0;
    if (rateClass < 0 || rateClass >= RateClassCount) {
        std::cerr << "SystemKMCBinaryRateBase: invalid FCC rate class "
                  << rateClass << " from nini=" << nini << " nfin=" << nfin
                  << std::endl;
        std::exit(2);
    }
    return static_cast<signed char>(rateClass);
}
