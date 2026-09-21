#include "SystemKMCExactClassOptimized.h"
SystemKMCExactClassOptimized::SystemKMCExactClassOptimized()
    : siteStampCounter(0), bondStampCounter(0), optimizedClockIncrement(0.0),
      optimizedEventCount(0) {
}

SystemKMCExactClassOptimized::~SystemKMCExactClassOptimized() {
}

void SystemKMCExactClassOptimized::InitSubSys() {
    SystemKMCExactClassCore::InitSubSys();
    txtstream << "System exact-class optimized: KMCExactClassOptimized";
    logging(txtstream.str());
}

bool SystemKMCExactClassOptimized::RecomputeMCSIncrementEveryAttempt() const {
    // Q changes after essentially every executed event.  RunMC must request the
    // current value before each event, but it is precomputed after the previous
    // local update so this call itself is only a load.
    return true;
}

double SystemKMCExactClassOptimized::MCSIncrementPerAttempt() const {
    return optimizedClockIncrement;
}

void SystemKMCExactClassOptimized::InitializeZeroNeighborCounts() {
    zeroNeighborCount.assign(static_cast<std::size_t>(TotalAtoms), 0);
    for (unsigned long long site = 0;
         site < static_cast<unsigned long long>(TotalAtoms); ++site) {
        int count = 0;
        const int x = fixedX[site];
        const int y = fixedY[site];
        const int z = fixedZ[site];
        for (int d = 0; d < 12; ++d) {
            int xn, yn, zn;
            NeighborCoords(x, y, z, d, xn, yn, zn);
            if ((this->*GetSpecies)(xn, yn, zn) == 0) ++count;
        }
        zeroNeighborCount[site] = static_cast<unsigned char>(count);
    }
}

signed char SystemKMCExactClassOptimized::FastRateClassForKey(
    unsigned long long key) {

    const unsigned long long siteA = key / 6ULL;
    const int direction = static_cast<int>(key % 6ULL);
    const int x = fixedX[siteA];
    const int y = fixedY[siteA];
    const int z = fixedZ[siteA];
    int xn, yn, zn;
    NeighborCoords(x, y, z, direction, xn, yn, zn);
    const unsigned long long siteB = FixedSiteIndex(xn, yn, zn);

    const int sA = (this->*GetSpecies)(x, y, z);
    const int sB = (this->*GetSpecies)(xn, yn, zn);
    if (sA == sB) return InactiveClass;

    unsigned long long site0, site1;
    if (sA == 0 && sB == 1) {
        site0 = siteA;
        site1 = siteB;
    } else if (sA == 1 && sB == 0) {
        site0 = siteB;
        site1 = siteA;
    } else {
        std::cerr << "KMCExactClassOptimized requires binary species 0/1."
                  << std::endl;
        std::exit(2);
    }

    // For an A(0)-B(1) bond, the B endpoint counts the A partner as one of
    // its species-0 neighbours while the A endpoint does not count the B
    // partner.  Excluding the exchange partner therefore gives
    // delta = n0(A) - (n0(B)-1).
    const int delta = static_cast<int>(zeroNeighborCount[site0])
                    - static_cast<int>(zeroNeighborCount[site1]) + 1;
    const int rateClass = (delta > 0) ? delta : 0;
    if (rateClass < 0 || rateClass >= RateClassCount) {
        std::cerr << "KMCExactClassOptimized: invalid rate class "
                  << rateClass << " delta=" << delta << std::endl;
        std::exit(2);
    }
    return static_cast<signed char>(rateClass);
}

void SystemKMCExactClassOptimized::BeginSiteDeltaTransaction() {
    if (++siteStampCounter == 0) {
        std::fill(siteTouchStamp.begin(), siteTouchStamp.end(), 0U);
        siteStampCounter = 1;
    }
    touchedSiteScratch.clear();
    changedSiteScratch.clear();
}

void SystemKMCExactClassOptimized::AccumulateZeroNeighborDelta(
    unsigned long long site, int delta) {

    if (siteTouchStamp[site] != siteStampCounter) {
        siteTouchStamp[site] = siteStampCounter;
        zeroNeighborDelta[site] = 0;
        touchedSiteScratch.push_back(site);
    }
    zeroNeighborDelta[site] = static_cast<signed char>(
        static_cast<int>(zeroNeighborDelta[site]) + delta);
}

void SystemKMCExactClassOptimized::BuildAffectedSitesAndBonds(
    int x, int y, int z, int xn, int yn, int zn, int s, int sn) {

    BeginSiteDeltaTransaction();

    const int deltaAtFirstSite = (sn == 0 ? 1 : 0) - (s == 0 ? 1 : 0);
    const int deltaAtSecondSite = (s == 0 ? 1 : 0) - (sn == 0 ? 1 : 0);

    for (int d = 0; d < 12; ++d) {
        int ax, ay, az;
        NeighborCoords(x, y, z, d, ax, ay, az);
        AccumulateZeroNeighborDelta(FixedSiteIndex(ax, ay, az), deltaAtFirstSite);
        NeighborCoords(xn, yn, zn, d, ax, ay, az);
        AccumulateZeroNeighborDelta(FixedSiteIndex(ax, ay, az), deltaAtSecondSite);
    }

    // A bond class depends only on the endpoint species and the species-0
    // neighbour counts at its endpoints.  Common neighbours of the exchanged
    // sites receive +1 and -1 and therefore need no reclassification unless
    // they are themselves an exchanged endpoint (for which the species changes).
    const unsigned long long firstSite = FixedSiteIndex(x, y, z);
    const unsigned long long secondSite = FixedSiteIndex(xn, yn, zn);
    for (unsigned long long site : touchedSiteScratch) {
        if (zeroNeighborDelta[site] != 0 || site == firstSite || site == secondSite) {
            changedSiteScratch.push_back(site);
        }
    }

    if (++bondStampCounter == 0) {
        std::fill(bondTouchStamp.begin(), bondTouchStamp.end(), 0U);
        bondStampCounter = 1;
    }
    affectedBondScratch.clear();
    for (unsigned long long site : changedSiteScratch) {
        for (int d = 0; d < 12; ++d) {
            const unsigned long long key = CanonicalBondKeyFromIncident(site, d);
            if (bondTouchStamp[key] != bondStampCounter) {
                bondTouchStamp[key] = bondStampCounter;
                affectedBondScratch.push_back(key);
            }
        }
    }
}

void SystemKMCExactClassOptimized::ApplyZeroNeighborDeltas() {
    for (unsigned long long site : touchedSiteScratch) {
        const int updated = static_cast<int>(zeroNeighborCount[site])
                          + static_cast<int>(zeroNeighborDelta[site]);
        if (updated < 0 || updated > 12) {
            std::cerr << "KMCExactClassOptimized: neighbour-count range error"
                      << std::endl;
            std::exit(2);
        }
        zeroNeighborCount[site] = static_cast<unsigned char>(updated);
    }
}

void SystemKMCExactClassOptimized::RefreshBondKeyDifferential(
    unsigned long long key) {

    const signed char oldClass = bondClass[key];
    const signed char newClass = FastRateClassForKey(key);
    if (oldClass == newClass) return;

    if (oldClass != InactiveClass) RemoveBondKey(key);
    if (newClass != InactiveClass) AddBondKey(key, newClass);
}

void SystemKMCExactClassOptimized::UpdateOptimizedClock() {
    optimizedClockIncrement = (totalRateWeight > 0.0)
        ? 6.0 / totalRateWeight : 0.0;
}

void SystemKMCExactClassOptimized::DetermineBonds() {
    SystemKMCExactClassCore::DetermineBonds();
    InitializeZeroNeighborCounts();

    zeroNeighborDelta.assign(static_cast<std::size_t>(TotalAtoms), 0);
    siteTouchStamp.assign(static_cast<std::size_t>(TotalAtoms), 0U);
    bondTouchStamp.assign(static_cast<std::size_t>(TotalAtoms) * 6ULL, 0U);
    touchedSiteScratch.reserve(24);
    changedSiteScratch.reserve(20);
    affectedBondScratch.reserve(192);
    siteStampCounter = 0;
    bondStampCounter = 0;
    optimizedEventCount = 0;
    UpdateOptimizedClock();

    // One-time consistency check: the O(1) classifier must reproduce the
    // reference classifier for every canonical bond at initialization.
    const unsigned long long canonicalBondCount =
        static_cast<unsigned long long>(TotalAtoms) * 6ULL;
    for (unsigned long long key = 0; key < canonicalBondCount; ++key) {
        if (FastRateClassForKey(key) != bondClass[key]) {
            std::cerr << "KMCExactClassOptimized: initial fast-class mismatch at key "
                      << key << std::endl;
            std::exit(2);
        }
    }

    txtstream << "Exact-class optimized initialization: B=" << BondNumber
              << " Q=" << std::setprecision(16) << totalRateWeight
              << " dMCS=" << optimizedClockIncrement;
    logging(txtstream.str());
}

void SystemKMCExactClassOptimized::JumpAttempt() {
    if (BondNumber == 0 || totalRateWeight <= 0.0) return;

    const double target = RandomUnit() * totalRateWeight;
    double cumulative = 0.0;
    int selectedClass = -1;
    for (int c = 0; c < RateClassCount; ++c) {
        cumulative += static_cast<double>(classBonds[c].size()) * classProbability[c];
        if (!classBonds[c].empty() && target < cumulative) {
            selectedClass = c;
            break;
        }
    }
    if (selectedClass < 0) {
        for (int c = RateClassCount - 1; c >= 0; --c) {
            if (!classBonds[c].empty()) { selectedClass = c; break; }
        }
    }
    if (selectedClass < 0) return;

    const std::vector<unsigned long long>& selectedList = classBonds[selectedClass];
    const unsigned long long key = selectedList[RandomIndex(selectedList.size())];

    int x, y, z, xn, yn, zn;
    ResolveCanonicalBond(key, x, y, z, xn, yn, zn);
    const int s = (this->*GetSpecies)(x, y, z);
    const int sn = (this->*GetSpecies)(xn, yn, zn);
    if (s == sn) {
        std::cerr << "KMCExactClassOptimized: selected class bond became inactive"
                  << std::endl;
        std::exit(2);
    }

    BuildAffectedSitesAndBonds(x, y, z, xn, yn, zn, s, sn);

    ExchangeSites(x, y, z, xn, yn, zn, s, sn);
    ApplyZeroNeighborDeltas();
    ++NAccepted;
    ++BenchmarkActiveTableUpdateCount;

    for (unsigned long long affectedKey : affectedBondScratch) {
        RefreshBondKeyDifferential(affectedKey);
    }

    ++optimizedEventCount;
    // Incremental +/- updates are exact enough for event selection, but a rare
    // eight-term re-sum bounds floating-point drift over very long runs.
    if ((optimizedEventCount & 0xFFFFULL) == 0ULL) {
        RecalculateTotalRateWeight();
    }
    UpdateOptimizedClock();
}

void SystemKMCExactClassOptimized::ValidateOptimizedStateOrAbort(
    const char* context) {

    // Full O(N) audit used only at evaluation/checkpoint boundaries, outside
    // EvolutionWallSecondsCumulative.  It validates both the cached neighbour
    // counts and every exact class against the original reference classifier.
    for (unsigned long long site = 0;
         site < static_cast<unsigned long long>(TotalAtoms); ++site) {
        int brute = 0;
        const int x = fixedX[site];
        const int y = fixedY[site];
        const int z = fixedZ[site];
        for (int d = 0; d < 12; ++d) {
            int ax, ay, az;
            NeighborCoords(x, y, z, d, ax, ay, az);
            if ((this->*GetSpecies)(ax, ay, az) == 0) ++brute;
        }
        if (brute != static_cast<int>(zeroNeighborCount[site])) {
            std::cerr << "KMCExactClassOptimized: cached neighbour-count mismatch at "
                      << context << " site=" << site
                      << " cached=" << static_cast<int>(zeroNeighborCount[site])
                      << " brute=" << brute << std::endl;
            std::exit(2);
        }
    }

    unsigned long long activeCount = 0;
    double q = 0.0;
    const unsigned long long canonicalBondCount =
        static_cast<unsigned long long>(TotalAtoms) * 6ULL;
    for (unsigned long long key = 0; key < canonicalBondCount; ++key) {
        const signed char expected = RateClassForKey(key);
        if (expected != bondClass[key]) {
            std::cerr << "KMCExactClassOptimized: class mismatch at "
                      << context << " key=" << key
                      << " cached=" << static_cast<int>(bondClass[key])
                      << " expected=" << static_cast<int>(expected) << std::endl;
            std::exit(2);
        }
        if (expected != InactiveClass) {
            ++activeCount;
            q += classProbability[static_cast<int>(expected)];
        }
    }
    if (activeCount != BondNumber) {
        std::cerr << "KMCExactClassOptimized: active-count mismatch at "
                  << context << " cached=" << BondNumber
                  << " brute=" << activeCount << std::endl;
        std::exit(2);
    }
    const double tol = 1e-10 * std::max(1.0, std::abs(q));
    if (std::abs(q - totalRateWeight) > tol) {
        std::cerr << "KMCExactClassOptimized: Q mismatch at " << context
                  << " cached=" << std::setprecision(17) << totalRateWeight
                  << " brute=" << q << std::endl;
        std::exit(2);
    }
}

void SystemKMCExactClassOptimized::EvalSysStep(unsigned long long n) {
    if (SysEvalParam.find(" OptimizationAudit ") != pos) {
        ValidateOptimizedStateOrAbort("evaluation checkpoint");
    }
    SystemKMCExactClassCore::EvalSysStep(n);
}
