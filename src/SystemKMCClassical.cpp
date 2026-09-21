#include "SystemKMCClassical.h"
SystemKMCClassical::SystemKMCClassical() {
    classProbability.fill(0.0);
}

SystemKMCClassical::~SystemKMCClassical() {
}

void SystemKMCClassical::InitSubSys() {
    SystemKMCHomogenous::InitSubSys();
    if (NSpecies != 2) {
        std::cerr << "KMCClassical requires exactly two species." << std::endl;
        std::exit(2);
    }
    txtstream << "System: KMCClassical";
    logging(txtstream.str());
}

void SystemKMCClassical::CodeCalculation() {
    classProbability[0] = 1.0;
    for (int c = 1; c < RateClassCount; ++c) {
        classProbability[c] = std::exp(-Ea / kT * static_cast<double>(c));
    }
}

void SystemKMCClassical::NeighborCoords(
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

int SystemKMCClassical::CountSpeciesOneNeighbors(int x, int y, int z) {
    int count = 0;
    for (int d = 0; d < CoordinationNumber; ++d) {
        int xn, yn, zn;
        NeighborCoords(x, y, z, d, xn, yn, zn);
        if ((this->*GetSpecies)(xn, yn, zn) == 1) ++count;
    }
    return count;
}

void SystemKMCClassical::DetermineBonds() {
    BondNumber = 0;
    for (unsigned long long atom = 0;
         atom < static_cast<unsigned long long>(TotalAtoms); ++atom) {
        const int x = xpr[atom];
        const int y = ypr[atom];
        const int z = zpr[atom];
        const int s = (this->*GetSpecies)(x, y, z);
        for (int d = 0; d < 6; ++d) {
            int xn, yn, zn;
            NeighborCoords(x, y, z, d, xn, yn, zn);
            if (s != (this->*GetSpecies)(xn, yn, zn)) ++BondNumber;
        }
    }

    txtstream << "Classical initialization: B=" << BondNumber;
    logging(txtstream.str());
}

double SystemKMCClassical::MCSIncrementPerAttempt() const {
    if (TotalAtoms <= 0) return 0.0;
    // One classical MCS is N random site-neighbour proposals.
    return 1.0 / static_cast<double>(TotalAtoms);
}

double SystemKMCClassical::PairMetropolisProbability(
    int x, int y, int z, int xn, int yn, int zn,
    int s, int sn, int& bAtA, int& bAtB) {

    ++BenchmarkProbabilityEvaluationCount;
    if (s == sn || s < 0 || s > 1 || sn < 0 || sn > 1) {
        std::cerr << "KMCClassical: invalid unlike binary proposal" << std::endl;
        std::exit(2);
    }

    if (s == 0) {
        bAtA = CountSpeciesOneNeighbors(x, y, z);
        bAtB = CountSpeciesOneNeighbors(xn, yn, zn);
    } else {
        bAtA = CountSpeciesOneNeighbors(xn, yn, zn);
        bAtB = CountSpeciesOneNeighbors(x, y, z);
    }

    // For A=0, B=1, the selected A-B bond remains unlike after exchange and
    // cancels. The exact FCC NN Metropolis class is b(B)-b(A)+1.
    const int delta = bAtB - bAtA + 1;
    const int rateClass = (delta > 0) ? delta : 0;
    if (rateClass < 0 || rateClass >= RateClassCount) {
        std::cerr << "KMCClassical: invalid Metropolis class " << rateClass
                  << " delta=" << delta << std::endl;
        std::exit(2);
    }
    return classProbability[rateClass];
}

void SystemKMCClassical::UpdateBondNumberAfterAcceptedExchange(
    int bAtA, int bAtB) {

    // Exact change in the number of unlike undirected bonds for a binary A-B
    // exchange. The selected A-B bond itself remains unlike.
    const long long delta = 2LL * (1LL - static_cast<long long>(bAtA)
                                        + static_cast<long long>(bAtB));
    const long long updated = static_cast<long long>(BondNumber) + delta;
    if (updated < 0) {
        std::cerr << "KMCClassical: BondNumber underflow" << std::endl;
        std::exit(2);
    }
    BondNumber = static_cast<unsigned long long>(updated);
}

void SystemKMCClassical::JumpAttempt() {
    const unsigned long long selectedAtom =
        RandomIndex(static_cast<unsigned long long>(TotalAtoms));
    const int x = xpr[selectedAtom];
    const int y = ypr[selectedAtom];
    const int z = zpr[selectedAtom];
    const int s = (this->*GetSpecies)(x, y, z);

    const int direction = static_cast<int>(RandomIndex(CoordinationNumber));
    int xn, yn, zn;
    NeighborCoords(x, y, z, direction, xn, yn, zn);
    const int sn = (this->*GetSpecies)(xn, yn, zn);

    if (s == sn) {
        ++BenchmarkInactiveRejectedCount;
        return;
    }

    int bAtA = 0;
    int bAtB = 0;
    const double p = PairMetropolisProbability(
        x, y, z, xn, yn, zn, s, sn, bAtA, bAtB);
    if (!RandomAccept(p)) return;

    UpdateBondNumberAfterAcceptedExchange(bAtA, bAtB);

    // Classical blind proposal does not require atom-identity bookkeeping: xpr
    // remains the fixed site enumeration and only the site species are swapped.
    // Fixed-site enumeration keeps the proposal stream deterministic for a fixed seed.
    (this->*SetSpecies)(x, y, z, sn);
    (this->*SetSpecies)(xn, yn, zn, s);
    ++NAccepted;
}
