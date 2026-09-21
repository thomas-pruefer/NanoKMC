#include "SystemKMCRateCategoryOptimized.h"
SystemKMCRateCategoryOptimized::SystemKMCRateCategoryOptimized()
    : rateCategoryCount(4), rateCategoryWidth(2), totalMajorantWeight(0.0),
      categoryClockIncrement(0.0), categoryEventCount(0),
      siteStampCounter(0), bondStampCounter(0) {
    categoryUpperProbability.fill(0.0);
}

SystemKMCRateCategoryOptimized::~SystemKMCRateCategoryOptimized() {
}

void SystemKMCRateCategoryOptimized::InitSubSys() {
    SystemKMCBinaryRateBase::InitSubSys();
    if (NSpecies != 2) {
        std::cerr << "KMCRateCategoryOptimized requires exactly two species."
                  << std::endl;
        std::exit(2);
    }

    const int requested = Str2Int(InputParam("RateCategoryCount", "nanokmc.in"));
    rateCategoryCount = (requested < 0) ? 4 : requested;
    if (rateCategoryCount != 1 && rateCategoryCount != 2 &&
        rateCategoryCount != 4 && rateCategoryCount != 8) {
        std::cerr << "KMCRateCategoryOptimized: RateCategoryCount must be "
                     "1, 2, 4, or 8 (default: 4)." << std::endl;
        std::exit(2);
    }
    rateCategoryWidth = RateClassCount / rateCategoryCount;

    txtstream << "System rate-category optimized: KMCRateCategoryOptimized"
              << " - RateCategoryCount=" << rateCategoryCount
              << " - exact classes/category=" << rateCategoryWidth;
    logging(txtstream.str());
}

bool SystemKMCRateCategoryOptimized::RecomputeMCSIncrementEveryAttempt() const {
    // Qhat changes after accepted exchanges.  The value itself is cached, so
    // requesting it per proposal is only a load.
    return true;
}

double SystemKMCRateCategoryOptimized::MCSIncrementPerAttempt() const {
    return categoryClockIncrement;
}

int SystemKMCRateCategoryOptimized::CategoryForRateClass(int rateClass) const {
    if (rateClass < 0 || rateClass >= RateClassCount) return -1;
    return rateClass / rateCategoryWidth;
}

void SystemKMCRateCategoryOptimized::InitializeCategoryBounds() {
    categoryUpperProbability.fill(0.0);
    for (int g = 0; g < rateCategoryCount; ++g) {
        const int firstClass = g * rateCategoryWidth;
        // p_c decreases monotonically with c, hence the first exact class in
        // each contiguous category is its valid majorant probability.
        categoryUpperProbability[g] = classProbability[firstClass];
    }
}

void SystemKMCRateCategoryOptimized::InitializeZeroNeighborCounts() {
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

signed char SystemKMCRateCategoryOptimized::FastRateClassForKey(
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
        std::cerr << "KMCRateCategoryOptimized requires binary species 0/1."
                  << std::endl;
        std::exit(2);
    }

    const int delta = static_cast<int>(zeroNeighborCount[site0])
                    - static_cast<int>(zeroNeighborCount[site1]) + 1;
    const int rateClass = (delta > 0) ? delta : 0;
    if (rateClass < 0 || rateClass >= RateClassCount) {
        std::cerr << "KMCRateCategoryOptimized: invalid exact class "
                  << rateClass << " delta=" << delta << std::endl;
        std::exit(2);
    }
    return static_cast<signed char>(rateClass);
}

void SystemKMCRateCategoryOptimized::AddCategoryBond(
    unsigned long long key, signed char category) {

    if (category == InactiveCategory) return;
    if (bondCategory[key] != InactiveCategory) {
        std::cerr << "KMCRateCategoryOptimized: duplicate category insertion"
                  << std::endl;
        std::exit(2);
    }

    const int g = static_cast<int>(category);
    categoryRow[key] = categoryBonds[g].size();
    categoryBonds[g].push_back(key);
    bondCategory[key] = category;
    ++BondNumber;
    totalMajorantWeight += categoryUpperProbability[g];
}

void SystemKMCRateCategoryOptimized::RemoveCategoryBond(
    unsigned long long key) {

    const signed char category = bondCategory[key];
    if (category == InactiveCategory) return;

    const int g = static_cast<int>(category);
    std::vector<unsigned long long>& list = categoryBonds[g];
    const unsigned long long row = categoryRow[key];
    const unsigned long long lastKey = list.back();

    list[row] = lastKey;
    categoryRow[lastKey] = row;
    list.pop_back();

    bondCategory[key] = InactiveCategory;
    categoryRow[key] = 0;
    --BondNumber;
    totalMajorantWeight -= categoryUpperProbability[g];
}

void SystemKMCRateCategoryOptimized::RecalculateMajorantWeight() {
    totalMajorantWeight = 0.0;
    for (int g = 0; g < rateCategoryCount; ++g) {
        totalMajorantWeight += static_cast<double>(categoryBonds[g].size())
                             * categoryUpperProbability[g];
    }
}

void SystemKMCRateCategoryOptimized::UpdateClockIncrement() {
    categoryClockIncrement = (totalMajorantWeight > 0.0)
        ? 6.0 / totalMajorantWeight : 0.0;
}

void SystemKMCRateCategoryOptimized::BeginSiteDeltaTransaction() {
    if (++siteStampCounter == 0) {
        std::fill(siteTouchStamp.begin(), siteTouchStamp.end(), 0U);
        siteStampCounter = 1;
    }
    touchedSiteScratch.clear();
    changedSiteScratch.clear();
}

void SystemKMCRateCategoryOptimized::AccumulateZeroNeighborDelta(
    unsigned long long site, int delta) {

    if (siteTouchStamp[site] != siteStampCounter) {
        siteTouchStamp[site] = siteStampCounter;
        zeroNeighborDelta[site] = 0;
        touchedSiteScratch.push_back(site);
    }
    zeroNeighborDelta[site] = static_cast<signed char>(
        static_cast<int>(zeroNeighborDelta[site]) + delta);
}

void SystemKMCRateCategoryOptimized::BuildAffectedSitesAndBonds(
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

void SystemKMCRateCategoryOptimized::ApplyZeroNeighborDeltas() {
    for (unsigned long long site : touchedSiteScratch) {
        const int updated = static_cast<int>(zeroNeighborCount[site])
                          + static_cast<int>(zeroNeighborDelta[site]);
        if (updated < 0 || updated > 12) {
            std::cerr << "KMCRateCategoryOptimized: neighbour-count range error"
                      << std::endl;
            std::exit(2);
        }
        zeroNeighborCount[site] = static_cast<unsigned char>(updated);
    }
}

void SystemKMCRateCategoryOptimized::RefreshBondKeyDifferential(
    unsigned long long key) {

    const signed char oldCategory = bondCategory[key];
    const signed char exactClass = FastRateClassForKey(key);
    const signed char newCategory = (exactClass == InactiveClass)
        ? InactiveCategory
        : static_cast<signed char>(CategoryForRateClass(static_cast<int>(exactClass)));

    if (oldCategory == newCategory) return;
    if (oldCategory != InactiveCategory) RemoveCategoryBond(key);
    if (newCategory != InactiveCategory) AddCategoryBond(key, newCategory);
}

void SystemKMCRateCategoryOptimized::DetermineBonds() {
    BuildFixedSiteMap();
    InitializeCategoryBounds();
    for (int g = 0; g < MaxCategoryCount; ++g) categoryBonds[g].clear();

    const unsigned long long canonicalBondCount =
        static_cast<unsigned long long>(TotalAtoms) * 6ULL;
    bondCategory.assign(canonicalBondCount, InactiveCategory);
    categoryRow.assign(canonicalBondCount, 0ULL);

    InitializeZeroNeighborCounts();
    zeroNeighborDelta.assign(static_cast<std::size_t>(TotalAtoms), 0);
    siteTouchStamp.assign(static_cast<std::size_t>(TotalAtoms), 0U);
    bondTouchStamp.assign(canonicalBondCount, 0U);
    touchedSiteScratch.reserve(24);
    changedSiteScratch.reserve(20);
    affectedBondScratch.reserve(192);
    siteStampCounter = 0;
    bondStampCounter = 0;

    BondNumber = 0;
    totalMajorantWeight = 0.0;
    categoryEventCount = 0;

    for (unsigned long long key = 0; key < canonicalBondCount; ++key) {
        const signed char exactClass = FastRateClassForKey(key);
        if (exactClass == InactiveClass) continue;
        const signed char category = static_cast<signed char>(
            CategoryForRateClass(static_cast<int>(exactClass)));
        AddCategoryBond(key, category);
    }

    RecalculateMajorantWeight();
    UpdateClockIncrement();

    txtstream << "Rate-category initialization: B=" << BondNumber
              << " M=" << rateCategoryCount
              << " Qhat=" << std::setprecision(16) << totalMajorantWeight
              << " dMCS=" << categoryClockIncrement;
    logging(txtstream.str());
    for (int g = 0; g < rateCategoryCount; ++g) {
        const int firstClass = g * rateCategoryWidth;
        const int lastClass = firstClass + rateCategoryWidth - 1;
        txtstream << "Rate category " << g
                  << " exact_classes=" << firstClass << ".." << lastClass
                  << " N=" << categoryBonds[g].size()
                  << " p_hat=" << categoryUpperProbability[g];
        logging(txtstream.str());
    }

    ValidateOptimizedStateOrAbort("initialization");
}

void SystemKMCRateCategoryOptimized::JumpAttempt() {
    if (BondNumber == 0 || totalMajorantWeight <= 0.0) return;

    const double target = RandomUnit() * totalMajorantWeight;
    double cumulative = 0.0;
    int selectedCategory = -1;
    for (int g = 0; g < rateCategoryCount; ++g) {
        cumulative += static_cast<double>(categoryBonds[g].size())
                    * categoryUpperProbability[g];
        if (!categoryBonds[g].empty() && target < cumulative) {
            selectedCategory = g;
            break;
        }
    }
    if (selectedCategory < 0) {
        for (int g = rateCategoryCount - 1; g >= 0; --g) {
            if (!categoryBonds[g].empty()) { selectedCategory = g; break; }
        }
    }
    if (selectedCategory < 0) return;

    const std::vector<unsigned long long>& selectedList =
        categoryBonds[selectedCategory];
    const unsigned long long key =
        selectedList[RandomIndex(static_cast<unsigned long long>(selectedList.size()))];

    const signed char exactClass = FastRateClassForKey(key);
    if (exactClass == InactiveClass) {
        std::cerr << "KMCRateCategoryOptimized: selected category bond became inactive"
                  << std::endl;
        std::exit(2);
    }
    if (CategoryForRateClass(static_cast<int>(exactClass)) != selectedCategory) {
        std::cerr << "KMCRateCategoryOptimized: selected bond/category mismatch"
                  << std::endl;
        std::exit(2);
    }

    ++BenchmarkProbabilityEvaluationCount;
    const double exactProbability = classProbability[static_cast<int>(exactClass)];
    const double upperProbability = categoryUpperProbability[selectedCategory];
    const double residualAcceptance = exactProbability / upperProbability;
    if (!RandomAccept(residualAcceptance)) return;

    int x, y, z, xn, yn, zn;
    ResolveCanonicalBond(key, x, y, z, xn, yn, zn);
    const int s = (this->*GetSpecies)(x, y, z);
    const int sn = (this->*GetSpecies)(xn, yn, zn);
    if (s == sn) {
        std::cerr << "KMCRateCategoryOptimized: selected unlike bond became equal"
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

    ++categoryEventCount;
    if ((categoryEventCount & 0xFFFFULL) == 0ULL) {
        RecalculateMajorantWeight();
    }
    UpdateClockIncrement();
}

void SystemKMCRateCategoryOptimized::ValidateOptimizedStateOrAbort(
    const char* context) {

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
            std::cerr << "KMCRateCategoryOptimized: cached neighbour-count mismatch at "
                      << context << " site=" << site
                      << " cached=" << static_cast<int>(zeroNeighborCount[site])
                      << " brute=" << brute << std::endl;
            std::exit(2);
        }
    }

    unsigned long long activeCount = 0;
    double qhat = 0.0;
    const unsigned long long canonicalBondCount =
        static_cast<unsigned long long>(TotalAtoms) * 6ULL;
    for (unsigned long long key = 0; key < canonicalBondCount; ++key) {
        const signed char exactClass = RateClassForKey(key);
        const signed char expectedCategory = (exactClass == InactiveClass)
            ? InactiveCategory
            : static_cast<signed char>(
                CategoryForRateClass(static_cast<int>(exactClass)));
        if (expectedCategory != bondCategory[key]) {
            std::cerr << "KMCRateCategoryOptimized: category mismatch at "
                      << context << " key=" << key
                      << " cached=" << static_cast<int>(bondCategory[key])
                      << " expected=" << static_cast<int>(expectedCategory)
                      << std::endl;
            std::exit(2);
        }
        if (expectedCategory != InactiveCategory) {
            ++activeCount;
            qhat += categoryUpperProbability[static_cast<int>(expectedCategory)];
        }
    }

    if (activeCount != BondNumber) {
        std::cerr << "KMCRateCategoryOptimized: active-count mismatch at "
                  << context << " cached=" << BondNumber
                  << " brute=" << activeCount << std::endl;
        std::exit(2);
    }
    const double tol = 1e-10 * std::max(1.0, std::abs(qhat));
    if (std::abs(qhat - totalMajorantWeight) > tol) {
        std::cerr << "KMCRateCategoryOptimized: Qhat mismatch at " << context
                  << " cached=" << std::setprecision(17) << totalMajorantWeight
                  << " brute=" << qhat << std::endl;
        std::exit(2);
    }

    unsigned long long listPopulation = 0;
    for (int g = 0; g < rateCategoryCount; ++g) {
        listPopulation += categoryBonds[g].size();
    }
    if (listPopulation != BondNumber) {
        std::cerr << "KMCRateCategoryOptimized: list population mismatch at "
                  << context << std::endl;
        std::exit(2);
    }
}

void SystemKMCRateCategoryOptimized::EvalSysStep(unsigned long long n) {
    if (SysEvalParam.find(" OptimizationAudit ") != pos) {
        ValidateOptimizedStateOrAbort("evaluation checkpoint");
    }
    SystemKMCHomogenous::EvalSysStep(n);
}
