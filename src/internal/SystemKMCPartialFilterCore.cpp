#include "SystemKMCPartialFilterCore.h"
SystemKMCPartialFilterCore::SystemKMCPartialFilterCore()
    : invalidRow(std::numeric_limits<unsigned long long>::max()),
      partialClockIncrement(0.0) {
    classProbability.fill(0.0);
}

SystemKMCPartialFilterCore::~SystemKMCPartialFilterCore() {
}

void SystemKMCPartialFilterCore::InitSubSys() {
    SystemKMCHomogenous::InitSubSys();
    if (NSpecies != 2) {
        std::cerr << "KMCPartialFilterCore requires exactly two species." << std::endl;
        std::exit(2);
    }
}

void SystemKMCPartialFilterCore::CodeCalculation() {
    // Same finite FCC Metropolis factors as the other binary benchmark solvers.
    classProbability[0] = 1.0;
    for (int c = 1; c < RateClassCount; ++c) {
        classProbability[c] = std::exp(-Ea / kT * static_cast<double>(c));
    }
}

bool SystemKMCPartialFilterCore::RecomputeMCSIncrementEveryAttempt() const {
    // The clock denominator is M_A (eligible A sites), not BondNumber.  M_A can
    // change after an accepted exchange while BondNumber remains unchanged.
    // The returned value is precomputed, so this virtual call does not perform
    // a division on rejected proposals.
    return true;
}

double SystemKMCPartialFilterCore::MCSIncrementPerAttempt() const {
    return partialClockIncrement;
}

void SystemKMCPartialFilterCore::EvalSysStep(unsigned long long n) {
    SystemKMCHomogenous::EvalSysStep(n);
}

void SystemKMCPartialFilterCore::NeighborCoords(
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

bool SystemKMCPartialFilterCore::HasSpeciesOneNeighbor(
    unsigned long long atomIndex) {

    const int x = xpr[atomIndex];
    const int y = ypr[atomIndex];
    const int z = zpr[atomIndex];
    for (int d = 0; d < 12; ++d) {
        int xn, yn, zn;
        NeighborCoords(x, y, z, d, xn, yn, zn);
        if ((this->*GetSpecies)(xn, yn, zn) == 1) return true;
    }
    return false;
}

bool SystemKMCPartialFilterCore::IsEligibleA(
    unsigned long long atomIndex) {

    const int x = xpr[atomIndex];
    const int y = ypr[atomIndex];
    const int z = zpr[atomIndex];
    if ((this->*GetSpecies)(x, y, z) != 0) return false;
    return HasSpeciesOneNeighbor(atomIndex);
}

void SystemKMCPartialFilterCore::AddEligible(
    unsigned long long atomIndex) {

    if (eligibleRow[atomIndex] != invalidRow) return;
    eligibleRow[atomIndex] = static_cast<unsigned long long>(eligibleA.size());
    eligibleA.push_back(atomIndex);
}

void SystemKMCPartialFilterCore::RemoveEligible(
    unsigned long long atomIndex) {

    const unsigned long long row = eligibleRow[atomIndex];
    if (row == invalidRow) return;

    const unsigned long long lastAtom = eligibleA.back();
    eligibleA[row] = lastAtom;
    eligibleRow[lastAtom] = row;
    eligibleA.pop_back();
    eligibleRow[atomIndex] = invalidRow;
}

void SystemKMCPartialFilterCore::RefreshEligibility(
    unsigned long long atomIndex) {

    const bool shouldBeEligible = IsEligibleA(atomIndex);
    const bool isListed = (eligibleRow[atomIndex] != invalidRow);
    if (shouldBeEligible && !isListed) AddEligible(atomIndex);
    else if (!shouldBeEligible && isListed) RemoveEligible(atomIndex);
}

void SystemKMCPartialFilterCore::UpdateClockIncrement() {
    if (eligibleA.empty()) {
        partialClockIncrement = 0.0;
    } else {
        // Classical reference: N site-neighbour attempts per MCS.  A specific
        // undirected A-B bond is proposed with intensity 1/6 per MCS because
        // it can be selected from either orientation.  The partial filter picks
        // only the A orientation with probability 1/(12 M_A), hence
        // Delta MCS = 1/(2 M_A) restores the same per-bond proposal intensity.
        partialClockIncrement =
            1.0 / (2.0 * static_cast<double>(eligibleA.size()));
    }
}

void SystemKMCPartialFilterCore::ValidateEligiblePopulationOrAbort(
    const char* context) {

    unsigned long long bruteCount = 0;
    for (unsigned long long atom = 0;
         atom < static_cast<unsigned long long>(TotalAtoms); ++atom) {
        const bool expected = IsEligibleA(atom);
        const bool listed = (eligibleRow[atom] != invalidRow);
        if (expected != listed) {
            std::cerr << "KMCPartialFilterCore: eligibility mismatch at "
                      << context << " for atom " << atom
                      << " expected=" << expected << " listed=" << listed
                      << std::endl;
            std::exit(2);
        }
        if (expected) ++bruteCount;
    }
    if (bruteCount != eligibleA.size()) {
        std::cerr << "KMCPartialFilterCore: eligible population mismatch at "
                  << context << ": brute=" << bruteCount
                  << " list=" << eligibleA.size() << std::endl;
        std::exit(2);
    }
}

void SystemKMCPartialFilterCore::DetermineBonds() {
    // Build the eligible-A population and structural unlike-bond count directly.
    eligibleA.clear();
    eligibleA.reserve(static_cast<std::size_t>(AtomNumber[0]));
    eligibleRow.assign(static_cast<std::size_t>(TotalAtoms), invalidRow);

    BondNumber = 0;
    for (unsigned long long atom = 0;
         atom < static_cast<unsigned long long>(TotalAtoms); ++atom) {
        const int x = xpr[atom];
        const int y = ypr[atom];
        const int z = zpr[atom];
        if (IsEligibleA(atom)) AddEligible(atom);
        for (int d = 0; d < 6; ++d) {
            int xn, yn, zn;
            NeighborCoords(x, y, z, d, xn, yn, zn);
            if ((this->*GetSpecies)(x, y, z) !=
                (this->*GetSpecies)(xn, yn, zn)) {
                ++BondNumber;
            }
        }
    }

    UpdateClockIncrement();
    ValidateEligiblePopulationOrAbort("initialization");

    txtstream << "Partial-filter initialization: B=" << BondNumber
              << " eligible_A=" << eligibleA.size()
              << " dMCS=" << std::setprecision(16) << partialClockIncrement;
    logging(txtstream.str());
}

int SystemKMCPartialFilterCore::CountSpeciesZeroNeighborsExcept(
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

double SystemKMCPartialFilterCore::PairMetropolisProbability(
    int x, int y, int z, int xn, int yn, int zn) {

    ++BenchmarkProbabilityEvaluationCount;

    // JumpAttempt always orients the proposal from selected A (species 0) to
    // its randomly chosen neighbour.  Only A-B proposals reach this function.
    const int s = (this->*GetSpecies)(x, y, z);
    const int sn = (this->*GetSpecies)(xn, yn, zn);
    if (s != 0 || sn != 1) {
        std::cerr << "KMCPartialFilterCore: probability requested for non A-B pair"
                  << std::endl;
        std::exit(2);
    }

    const int nini = CountSpeciesZeroNeighborsExcept(x, y, z, xn, yn, zn);
    const int nfin = CountSpeciesZeroNeighborsExcept(xn, yn, zn, x, y, z);
    const int delta = nini - nfin;
    const int rateClass = (delta > 0) ? delta : 0;
    if (rateClass < 0 || rateClass >= RateClassCount) {
        std::cerr << "KMCPartialFilterCore: invalid rate class " << rateClass
                  << " from nini=" << nini << " nfin=" << nfin << std::endl;
        std::exit(2);
    }
    return classProbability[rateClass];
}
