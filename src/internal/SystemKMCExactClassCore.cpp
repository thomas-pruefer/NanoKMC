#include "SystemKMCExactClassCore.h"
SystemKMCExactClassCore::SystemKMCExactClassCore()
    : totalRateWeight(0.0) {
}

SystemKMCExactClassCore::~SystemKMCExactClassCore() {
}

void SystemKMCExactClassCore::InitSubSys() {
    SystemKMCBinaryRateBase::InitSubSys();
    totalRateWeight = 0.0;
}

bool SystemKMCExactClassCore::RecomputeMCSIncrementEveryAttempt() const {
    return true;
}

double SystemKMCExactClassCore::MCSIncrementPerAttempt() const {
    if (totalRateWeight <= 0.0) return 0.0;
    return 6.0 / totalRateWeight;
}

void SystemKMCExactClassCore::EvalSysStep(unsigned long long n) {
    SystemKMCHomogenous::EvalSysStep(n);
}

void SystemKMCExactClassCore::AddBondKey(
    unsigned long long key, signed char rateClass) {

    if (rateClass == InactiveClass) return;
    if (bondClass[key] != InactiveClass) {
        std::cerr << "KMCExactClassCore: attempted duplicate bond insertion" << std::endl;
        std::exit(2);
    }
    const int c = static_cast<int>(rateClass);
    bondRow[key] = classBonds[c].size();
    classBonds[c].push_back(key);
    bondClass[key] = rateClass;
    ++BondNumber;
    totalRateWeight += classProbability[c];
}

void SystemKMCExactClassCore::RemoveBondKey(unsigned long long key) {
    const signed char rateClass = bondClass[key];
    if (rateClass == InactiveClass) return;

    const int c = static_cast<int>(rateClass);
    std::vector<unsigned long long>& list = classBonds[c];
    const unsigned long long row = bondRow[key];
    const unsigned long long lastKey = list.back();

    list[row] = lastKey;
    bondRow[lastKey] = row;
    list.pop_back();

    bondClass[key] = InactiveClass;
    bondRow[key] = 0;
    --BondNumber;
    totalRateWeight -= classProbability[c];
}

void SystemKMCExactClassCore::CollectAffectedBondKeys(
    int x, int y, int z, int xn, int yn, int zn,
    std::vector<unsigned long long>& keys) {

    std::vector<unsigned long long> sites;
    sites.reserve(26);
    const auto addSite = [&](int sx, int sy, int sz) {
        const unsigned long long site = FixedSiteIndex(sx, sy, sz);
        if (site != invalidFixedSite) sites.push_back(site);
    };

    addSite(x, y, z); addSite(xn, yn, zn);
    for (int d = 0; d < 12; ++d) {
        int ax, ay, az;
        NeighborCoords(x, y, z, d, ax, ay, az); addSite(ax, ay, az);
        NeighborCoords(xn, yn, zn, d, ax, ay, az); addSite(ax, ay, az);
    }

    std::sort(sites.begin(), sites.end());
    sites.erase(std::unique(sites.begin(), sites.end()), sites.end());

    keys.clear(); keys.reserve(sites.size() * 12ULL);
    for (unsigned long long site : sites) {
        for (int d = 0; d < 12; ++d) {
            keys.push_back(CanonicalBondKeyFromIncident(site, d));
        }
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
}

void SystemKMCExactClassCore::RecalculateTotalRateWeight() {
    totalRateWeight = 0.0;
    for (int c = 0; c < RateClassCount; ++c) {
        totalRateWeight += static_cast<double>(classBonds[c].size())
                         * classProbability[c];
    }
}

void SystemKMCExactClassCore::ValidateClassPopulationOrAbort(
    const char* context) const {
    unsigned long long sum = 0;
    for (int c = 0; c < RateClassCount; ++c) sum += classBonds[c].size();
    if (sum != BondNumber) {
        std::cerr << "KMCExactClassCore: class population mismatch at " << context
                  << ": sum=" << sum << " BondNumber=" << BondNumber << std::endl;
        std::exit(2);
    }
}

void SystemKMCExactClassCore::DetermineBonds() {
    BuildFixedSiteMap();
    for (int c = 0; c < RateClassCount; ++c) classBonds[c].clear();

    const unsigned long long canonicalBondCount =
        static_cast<unsigned long long>(TotalAtoms) * 6ULL;
    bondClass.assign(canonicalBondCount, InactiveClass);
    bondRow.assign(canonicalBondCount, 0ULL);
    BondNumber = 0;
    totalRateWeight = 0.0;

    for (unsigned long long key = 0; key < canonicalBondCount; ++key) {
        const signed char c = RateClassForKey(key);
        if (c != InactiveClass) AddBondKey(key, c);
    }
    ValidateClassPopulationOrAbort("initialization");
    RecalculateTotalRateWeight();

    txtstream << "Exact-class initialization: B=" << BondNumber
              << " Q=" << std::setprecision(16) << totalRateWeight;
    logging(txtstream.str());
    for (int c = 0; c < RateClassCount; ++c) {
        txtstream << "Exact class " << c << ": N=" << classBonds[c].size()
                  << " p=" << classProbability[c];
        logging(txtstream.str());
    }
}

void SystemKMCExactClassCore::JumpAttempt() {
    if (BondNumber == 0 || totalRateWeight <= 0.0) return;

    const double target = RandomUnit() * totalRateWeight;
    double cumulative = 0.0;
    int selectedClass = -1;
    for (int c = 0; c < RateClassCount; ++c) {
        cumulative += static_cast<double>(classBonds[c].size()) * classProbability[c];
        if (!classBonds[c].empty() && target < cumulative) {
            selectedClass = c; break;
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
        std::cerr << "KMCExactClassCore: selected class bond became inactive" << std::endl;
        std::exit(2);
    }

    std::vector<unsigned long long> affected;
    CollectAffectedBondKeys(x, y, z, xn, yn, zn, affected);
    for (unsigned long long affectedKey : affected) RemoveBondKey(affectedKey);

    ExchangeSites(x, y, z, xn, yn, zn, s, sn);
    ++NAccepted;
    ++BenchmarkActiveTableUpdateCount;

    for (unsigned long long affectedKey : affected) {
        const signed char c = RateClassForKey(affectedKey);
        if (c != InactiveClass) AddBondKey(affectedKey, c);
    }
    RecalculateTotalRateWeight();
    ValidateClassPopulationOrAbort("local update");
}
