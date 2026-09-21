#include "SystemKMCPartialFilterOptimized.h"
SystemKMCPartialFilterOptimized::SystemKMCPartialFilterOptimized()
    : atomStampCounter(0) {
}

SystemKMCPartialFilterOptimized::~SystemKMCPartialFilterOptimized() {
}

void SystemKMCPartialFilterOptimized::InitSubSys() {
    SystemKMCPartialFilterCore::InitSubSys();
    txtstream << "System partial-filter optimized: KMCPartialFilterOptimized";
    logging(txtstream.str());
}

unsigned long long SystemKMCPartialFilterOptimized::CoordLinearIndexFast(
    int x, int y, int z) const {
    return (static_cast<unsigned long long>(x) * static_cast<unsigned long long>(ly1)
          + static_cast<unsigned long long>(y)) * static_cast<unsigned long long>(lz1)
          + static_cast<unsigned long long>(z);
}

unsigned char SystemKMCPartialFilterOptimized::SpeciesOneCountAt(
    int x, int y, int z) const {
    return speciesOneNeighborCount[CoordLinearIndexFast(x, y, z)];
}

void SystemKMCPartialFilterOptimized::SetSpeciesOneCountAt(
    int x, int y, int z, int value) {
    if (value < 0 || value > 12) {
        std::cerr << "KMCPartialFilterOptimized: neighbour-count range error"
                  << std::endl;
        std::exit(2);
    }
    speciesOneNeighborCount[CoordLinearIndexFast(x, y, z)] =
        static_cast<unsigned char>(value);
}

void SystemKMCPartialFilterOptimized::InitializeSpeciesOneNeighborCounts() {
    const unsigned long long fullGridSize =
        static_cast<unsigned long long>(lx1) *
        static_cast<unsigned long long>(ly1) *
        static_cast<unsigned long long>(lz1);
    speciesOneNeighborCount.assign(fullGridSize, 0);

    for (unsigned long long atom = 0;
         atom < static_cast<unsigned long long>(TotalAtoms); ++atom) {
        const int x = xpr[atom];
        const int y = ypr[atom];
        const int z = zpr[atom];
        int count = 0;
        for (int d = 0; d < 12; ++d) {
            int xn, yn, zn;
            NeighborCoords(x, y, z, d, xn, yn, zn);
            if ((this->*GetSpecies)(xn, yn, zn) == 1) ++count;
        }
        SetSpeciesOneCountAt(x, y, z, count);
    }
}

void SystemKMCPartialFilterOptimized::DetermineBonds() {
    // Let the internal core construct the validated initial eligible-A population
    // and structural unlike-bond count.
    SystemKMCPartialFilterCore::DetermineBonds();
    InitializeSpeciesOneNeighborCounts();

    atomTouchStamp.assign(static_cast<std::size_t>(TotalAtoms), 0U);
    affectedAtomFastScratch.reserve(24);
    atomStampCounter = 0;

    // Verify the cached neighbour counts against the core eligibility definition
    // once at initialization. This is outside the evolution timer.
    for (unsigned long long atom = 0;
         atom < static_cast<unsigned long long>(TotalAtoms); ++atom) {
        const int x = xpr[atom];
        const int y = ypr[atom];
        const int z = zpr[atom];
        const bool fastEligible = ((this->*GetSpecies)(x, y, z) == 0)
                               && (SpeciesOneCountAt(x, y, z) > 0);
        const bool listed = (eligibleRow[atom] != invalidRow);
        if (fastEligible != listed) {
            std::cerr << "KMCPartialFilterOptimized: initial eligibility mismatch"
                      << " atom=" << atom << std::endl;
            std::exit(2);
        }
    }

    ValidateFastProbabilitySampleOrAbort("initialization");

    txtstream << "Partial-filter optimized initialization: B=" << BondNumber
              << " eligible_A=" << eligibleA.size()
              << " dMCS=" << std::setprecision(16) << partialClockIncrement;
    logging(txtstream.str());
}

double SystemKMCPartialFilterOptimized::FastPairMetropolisProbability(
    int x, int y, int z, int xn, int yn, int zn) {

    ++BenchmarkProbabilityEvaluationCount;
    const int s = (this->*GetSpecies)(x, y, z);
    const int sn = (this->*GetSpecies)(xn, yn, zn);
    if (s != 0 || sn != 1) {
        std::cerr << "KMCPartialFilterOptimized: rate requested for non A-B pair"
                  << std::endl;
        std::exit(2);
    }

    // For binary A=0/B=1: n0 = 12-n1.  Excluding the A exchange partner at
    // the B endpoint yields delta = n1(B)-n1(A)+1.
    const int delta = static_cast<int>(SpeciesOneCountAt(xn, yn, zn))
                    - static_cast<int>(SpeciesOneCountAt(x, y, z)) + 1;
    const int rateClass = (delta > 0) ? delta : 0;
    if (rateClass < 0 || rateClass >= RateClassCount) {
        std::cerr << "KMCPartialFilterOptimized: invalid rate class "
                  << rateClass << " delta=" << delta << std::endl;
        std::exit(2);
    }
    return classProbability[rateClass];
}

void SystemKMCPartialFilterOptimized::UpdateSpeciesOneNeighborCountsAfterExchange(
    int x, int y, int z, int xn, int yn, int zn, int s, int sn) {

    const int deltaFirst = (sn == 1 ? 1 : 0) - (s == 1 ? 1 : 0);
    const int deltaSecond = (s == 1 ? 1 : 0) - (sn == 1 ? 1 : 0);
    for (int d = 0; d < 12; ++d) {
        int ax, ay, az;
        NeighborCoords(x, y, z, d, ax, ay, az);
        SetSpeciesOneCountAt(ax, ay, az,
            static_cast<int>(SpeciesOneCountAt(ax, ay, az)) + deltaFirst);
        NeighborCoords(xn, yn, zn, d, ax, ay, az);
        SetSpeciesOneCountAt(ax, ay, az,
            static_cast<int>(SpeciesOneCountAt(ax, ay, az)) + deltaSecond);
    }
}

void SystemKMCPartialFilterOptimized::CollectAffectedAtomsFast(
    int x, int y, int z, int xn, int yn, int zn) {

    if (++atomStampCounter == 0) {
        std::fill(atomTouchStamp.begin(), atomTouchStamp.end(), 0U);
        atomStampCounter = 1;
    }
    affectedAtomFastScratch.clear();

    const auto addAtomAt = [&](int sx, int sy, int sz) {
        const unsigned long long atom = xyzpointer[sx][sy][sz][0];
        if (atomTouchStamp[atom] != atomStampCounter) {
            atomTouchStamp[atom] = atomStampCounter;
            affectedAtomFastScratch.push_back(atom);
        }
    };

    addAtomAt(x, y, z);
    addAtomAt(xn, yn, zn);
    for (int d = 0; d < 12; ++d) {
        int ax, ay, az;
        NeighborCoords(x, y, z, d, ax, ay, az);
        addAtomAt(ax, ay, az);
        NeighborCoords(xn, yn, zn, d, ax, ay, az);
        addAtomAt(ax, ay, az);
    }
}

void SystemKMCPartialFilterOptimized::RefreshEligibilityFast(
    unsigned long long atomIndex) {

    const int x = xpr[atomIndex];
    const int y = ypr[atomIndex];
    const int z = zpr[atomIndex];
    const bool shouldBeEligible = ((this->*GetSpecies)(x, y, z) == 0)
                               && (SpeciesOneCountAt(x, y, z) > 0);
    const bool isListed = (eligibleRow[atomIndex] != invalidRow);
    if (shouldBeEligible && !isListed) AddEligible(atomIndex);
    else if (!shouldBeEligible && isListed) RemoveEligible(atomIndex);
}

void SystemKMCPartialFilterOptimized::UpdateBondNumberAnalytically(
    unsigned char bAtA, unsigned char bAtB) {

    // Before the exchange x is A and y is B.  The selected A-B bond stays
    // unlike.  For the remaining 22 incident bonds the exact interface-count
    // change is 2*(1 - nB(x) + nB(y)), using pre-exchange neighbour counts.
    const long long delta = 2LL * (1LL - static_cast<long long>(bAtA)
                                        + static_cast<long long>(bAtB));
    const long long updated = static_cast<long long>(BondNumber) + delta;
    if (updated < 0) {
        std::cerr << "KMCPartialFilterOptimized: BondNumber underflow"
                  << std::endl;
        std::exit(2);
    }
    BondNumber = static_cast<unsigned long long>(updated);
}

void SystemKMCPartialFilterOptimized::JumpAttempt() {
    if (eligibleA.empty()) return;

    const unsigned long long selectedAtom =
        eligibleA[RandomIndex(static_cast<unsigned long long>(eligibleA.size()))];
    const int x = xpr[selectedAtom];
    const int y = ypr[selectedAtom];
    const int z = zpr[selectedAtom];

#ifndef NDEBUG
    if ((this->*GetSpecies)(x, y, z) != 0 || SpeciesOneCountAt(x, y, z) == 0) {
        std::cerr << "KMCPartialFilterOptimized: eligible-list invariant failure"
                  << std::endl;
        std::exit(2);
    }
#endif

    const int direction = static_cast<int>(RandomIndex(12ULL));
    int xn, yn, zn;
    NeighborCoords(x, y, z, direction, xn, yn, zn);
    const int sn = (this->*GetSpecies)(xn, yn, zn);
    if (sn != 1) {
        ++BenchmarkInactiveRejectedCount;
        return;
    }

    const double jumpProbability = FastPairMetropolisProbability(x, y, z, xn, yn, zn);
    if (!RandomAccept(jumpProbability)) return;

    const unsigned char bAtA = SpeciesOneCountAt(x, y, z);
    const unsigned char bAtB = SpeciesOneCountAt(xn, yn, zn);
    UpdateBondNumberAnalytically(bAtA, bAtB);

    const int s = 0;
    ExchangeSites(x, y, z, xn, yn, zn, s, sn);
    UpdateSpeciesOneNeighborCountsAfterExchange(x, y, z, xn, yn, zn, s, sn);
    ++NAccepted;
    ++BenchmarkActiveTableUpdateCount;

    CollectAffectedAtomsFast(x, y, z, xn, yn, zn);
    for (unsigned long long atom : affectedAtomFastScratch) {
        RefreshEligibilityFast(atom);
    }
    UpdateClockIncrement();

    // BondNumber is maintained exactly by the local binary-exchange identity
    // above and independently recounted by the Benchmark evaluator at checkpoints.
}

void SystemKMCPartialFilterOptimized::ValidateOptimizedStateOrAbort(
    const char* context) {

    unsigned long long eligibleCount = 0;
    for (unsigned long long atom = 0;
         atom < static_cast<unsigned long long>(TotalAtoms); ++atom) {
        const int x = xpr[atom];
        const int y = ypr[atom];
        const int z = zpr[atom];
        int bruteB = 0;
        for (int d = 0; d < 12; ++d) {
            int ax, ay, az;
            NeighborCoords(x, y, z, d, ax, ay, az);
            if ((this->*GetSpecies)(ax, ay, az) == 1) ++bruteB;
        }
        const int cachedB = static_cast<int>(SpeciesOneCountAt(x, y, z));
        if (bruteB != cachedB) {
            std::cerr << "KMCPartialFilterOptimized: cached neighbour-count mismatch at "
                      << context << " atom=" << atom
                      << " cached=" << cachedB << " brute=" << bruteB << std::endl;
            std::exit(2);
        }
        const bool expectedEligible = ((this->*GetSpecies)(x, y, z) == 0) && bruteB > 0;
        const bool listed = (eligibleRow[atom] != invalidRow);
        if (expectedEligible != listed) {
            std::cerr << "KMCPartialFilterOptimized: eligibility mismatch at "
                      << context << " atom=" << atom << std::endl;
            std::exit(2);
        }
        if (expectedEligible) ++eligibleCount;
    }
    if (eligibleCount != eligibleA.size()) {
        std::cerr << "KMCPartialFilterOptimized: eligible-count mismatch at "
                  << context << " cached=" << eligibleA.size()
                  << " brute=" << eligibleCount << std::endl;
        std::exit(2);
    }
}

void SystemKMCPartialFilterOptimized::EvalSysStep(unsigned long long n) {
    if (SysEvalParam.find(" OptimizationAudit ") != pos) {
        ValidateOptimizedStateOrAbort("evaluation checkpoint");
        ValidateFastProbabilitySampleOrAbort("evaluation checkpoint");
    }
    SystemKMCPartialFilterCore::EvalSysStep(n);
}

void SystemKMCPartialFilterOptimized::ValidateFastProbabilitySampleOrAbort(
    const char* context) {

    unsigned long long checked = 0;
    constexpr unsigned long long maxChecks = 4096ULL;
    for (unsigned long long atom = 0;
         atom < static_cast<unsigned long long>(TotalAtoms) && checked < maxChecks; ++atom) {
        const int x = xpr[atom];
        const int y = ypr[atom];
        const int z = zpr[atom];
        if ((this->*GetSpecies)(x, y, z) != 0) continue;
        for (int d = 0; d < 12 && checked < maxChecks; ++d) {
            int xn, yn, zn;
            NeighborCoords(x, y, z, d, xn, yn, zn);
            if ((this->*GetSpecies)(xn, yn, zn) != 1) continue;
            const double fast = FastPairMetropolisProbability(x, y, z, xn, yn, zn);
            const double reference = PairMetropolisProbability(x, y, z, xn, yn, zn);
            const double tol = 1e-14 * std::max(1.0, std::abs(reference));
            if (std::abs(fast - reference) > tol) {
                std::cerr << "KMCPartialFilterOptimized: fast probability mismatch at "
                          << context << " fast=" << std::setprecision(17) << fast
                          << " reference=" << reference << std::endl;
                std::exit(2);
            }
            ++checked;
        }
    }
}
