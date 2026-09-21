#include "SystemKMCActiveFilteredBase.h"
SystemKMCActiveFilteredBase::SystemKMCActiveFilteredBase()
    : bondStampCounter(0) {
}

SystemKMCActiveFilteredBase::~SystemKMCActiveFilteredBase() {
}

void SystemKMCActiveFilteredBase::InitSubSys() {
    SystemKMCHomogenous::InitSubSys();
    if (NSpecies < 2) {
        std::cerr << "Active-Filtered KMC requires at least two species." << std::endl;
        std::exit(2);
    }
    txtstream << "System: " << AcceptanceBackendName();
    logging(txtstream.str());
}

unsigned long long SystemKMCActiveFilteredBase::CoordLinearIndex(
    int x, int y, int z) const {
    return (static_cast<unsigned long long>(x) *
            static_cast<unsigned long long>(ly1)
          + static_cast<unsigned long long>(y)) *
            static_cast<unsigned long long>(lz1)
          + static_cast<unsigned long long>(z);
}

void SystemKMCActiveFilteredBase::NeighborCoords(
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

int SystemKMCActiveFilteredBase::OppositeDirection(int direction) const {
    static constexpr int opposite[CoordinationNumber] =
        {8,6,7,9,11,10,1,2,0,3,5,4};
    if (direction < 0 || direction >= CoordinationNumber) return -1;
    return opposite[direction];
}

void SystemKMCActiveFilteredBase::BuildFixedTopology() {
    const auto nSites = static_cast<unsigned long long>(TotalAtoms);
    if (nSites * static_cast<unsigned long long>(CanonicalDirectionCount) >=
        static_cast<unsigned long long>(InvalidRow32)) {
        std::cerr << "Active-Filtered: lattice too large for 32-bit canonical keys"
                  << std::endl;
        std::exit(2);
    }

    fixedX.resize(nSites);
    fixedY.resize(nSites);
    fixedZ.resize(nSites);

    const unsigned long long fullGridSize =
        static_cast<unsigned long long>(lx1) *
        static_cast<unsigned long long>(ly1) *
        static_cast<unsigned long long>(lz1);
    coordToSite.assign(static_cast<std::size_t>(fullGridSize), InvalidSite32);

    for (std::uint32_t site = 0;
         site < static_cast<std::uint32_t>(nSites); ++site) {
        const int x = xpr[site];
        const int y = ypr[site];
        const int z = zpr[site];
        fixedX[site] = x;
        fixedY[site] = y;
        fixedZ[site] = z;
        coordToSite[CoordLinearIndex(x,y,z)] = site;
    }

    neighborSite.resize(static_cast<std::size_t>(nSites) * CoordinationNumber);
    for (std::uint32_t site = 0;
         site < static_cast<std::uint32_t>(nSites); ++site) {
        for (int d = 0; d < CoordinationNumber; ++d) {
            int xn, yn, zn;
            NeighborCoords(fixedX[site], fixedY[site], fixedZ[site], d,
                           xn, yn, zn);
            const std::uint32_t nb = coordToSite[CoordLinearIndex(xn,yn,zn)];
            if (nb == InvalidSite32) {
                std::cerr << "Active-Filtered: invalid FCC neighbour mapping"
                          << std::endl;
                std::exit(2);
            }
            neighborSite[static_cast<std::size_t>(site) * CoordinationNumber + d] = nb;
        }
    }
}

void SystemKMCActiveFilteredBase::InitializeSiteSpeciesCache() {
    if (NSpecies > 256) {
        std::cerr << "Active-Filtered currently supports at most 256 species in its fixed-site cache." << std::endl;
        std::exit(2);
    }
    siteSpecies.resize(static_cast<std::size_t>(TotalAtoms));
    for (std::uint32_t site = 0;
         site < static_cast<std::uint32_t>(TotalAtoms); ++site) {
        const int species =
            (this->*GetSpecies)(fixedX[site], fixedY[site], fixedZ[site]);
        if (species < 0 || species >= NSpecies) {
            std::cerr << "Active-Filtered: invalid species in fixed-site cache"
                      << std::endl;
            std::exit(2);
        }
        siteSpecies[site] = static_cast<unsigned char>(species);
    }
}

void SystemKMCActiveFilteredBase::AddActiveKey(std::uint32_t key) {
    if (activeRow[key] != InvalidRow32) return;
    const std::uint32_t row = static_cast<std::uint32_t>(activeKeys.size());
    activeRow[key] = row;
    activeKeys.push_back(key);
}

// Dense-list deletion is O(1): move the last key into the removed row and
// repair its reverse row address. This swap-delete pattern is established
// event-list machinery; its role here is to keep the structural population
// compact without imposing ordering semantics.
void SystemKMCActiveFilteredBase::RemoveActiveKey(std::uint32_t key) {
    const std::uint32_t row = activeRow[key];
    if (row == InvalidRow32) return;
    const std::uint32_t lastKey = activeKeys.back();
    activeKeys[row] = lastKey;
    activeRow[lastKey] = row;
    activeKeys.pop_back();
    activeRow[key] = InvalidRow32;
}

void SystemKMCActiveFilteredBase::InitializeFastActiveList() {
    const std::size_t canonicalCount =
        static_cast<std::size_t>(TotalAtoms) * CanonicalDirectionCount;
    activeRow.assign(canonicalCount, InvalidRow32);
    activeKeys.clear();
    activeKeys.reserve(static_cast<std::size_t>(BondNumber));

    for (std::uint32_t site = 0;
         site < static_cast<std::uint32_t>(TotalAtoms); ++site) {
        const std::size_t base =
            static_cast<std::size_t>(site) * CoordinationNumber;
        for (int d = 0; d < CanonicalDirectionCount; ++d) {
            const std::uint32_t nb = neighborSite[base+d];
            if (siteSpecies[site] != siteSpecies[nb]) {
                AddActiveKey(site * CanonicalDirectionCount +
                             static_cast<std::uint32_t>(d));
            }
        }
    }
    BondNumber = static_cast<unsigned long long>(activeKeys.size());
    bondTouchStamp.assign(canonicalCount, 0U);
    affectedBondScratch.clear();
    affectedBondScratch.reserve(2 * CoordinationNumber);
    bondStampCounter = 0;
}

std::uint32_t SystemKMCActiveFilteredBase::CanonicalBondKeyFromIncident(
    std::uint32_t site, int direction) const {
    if (direction < CanonicalDirectionCount) {
        return site * CanonicalDirectionCount +
               static_cast<std::uint32_t>(direction);
    }
    const std::uint32_t nb =
        neighborSite[static_cast<std::size_t>(site) * CoordinationNumber + direction];
    const int opposite = OppositeDirection(direction);
    if (opposite < 0 || opposite >= CanonicalDirectionCount) {
        std::cerr << "Active-Filtered: invalid reverse canonical direction"
                  << std::endl;
        std::exit(2);
    }
    return nb * CanonicalDirectionCount +
           static_cast<std::uint32_t>(opposite);
}

// Only bonds incident on either exchanged FCC endpoint can change from like to
// unlike or vice versa. With z=12, the two incident sets share the selected
// bond, giving 2*z-1 = 23 unique undirected canonical bonds. The stamp array
// de-duplicates these keys without heap allocation or sorting in the hot path.
void SystemKMCActiveFilteredBase::CollectAffectedBondKeys(
    std::uint32_t firstSite, std::uint32_t secondSite) {
    if (++bondStampCounter == 0) {
        std::fill(bondTouchStamp.begin(), bondTouchStamp.end(), 0U);
        bondStampCounter = 1;
    }
    affectedBondScratch.clear();
    const auto addIncident = [&](std::uint32_t site) {
        for (int d = 0; d < CoordinationNumber; ++d) {
            const std::uint32_t key = CanonicalBondKeyFromIncident(site, d);
            if (bondTouchStamp[key] != bondStampCounter) {
                bondTouchStamp[key] = bondStampCounter;
                affectedBondScratch.push_back(key);
            }
        }
    };
    addIncident(firstSite);
    addIncident(secondSite);
}

void SystemKMCActiveFilteredBase::RefreshActiveKey(std::uint32_t key) {
    const std::uint32_t site = key / CanonicalDirectionCount;
    const int d = static_cast<int>(key % CanonicalDirectionCount);
    const std::uint32_t nb =
        neighborSite[static_cast<std::size_t>(site) * CoordinationNumber + d];
    const bool shouldBeActive = siteSpecies[site] != siteSpecies[nb];
    const bool isActive = activeRow[key] != InvalidRow32;
    if (shouldBeActive && !isActive) AddActiveKey(key);
    else if (!shouldBeActive && isActive) RemoveActiveKey(key);
}

void SystemKMCActiveFilteredBase::GatherLocalPairEnvironment(
    std::uint32_t firstSite, std::uint32_t secondSite,
    int canonicalDirection, LocalPairEnvironment& env) const {
    env.firstSite = firstSite;
    env.secondSite = secondSite;
    env.selectedDirection = canonicalDirection;
    env.firstSpecies = static_cast<int>(siteSpecies[firstSite]);
    env.secondSpecies = static_cast<int>(siteSpecies[secondSite]);

    const std::size_t firstBase =
        static_cast<std::size_t>(firstSite) * CoordinationNumber;
    const std::size_t secondBase =
        static_cast<std::size_t>(secondSite) * CoordinationNumber;
    for (int d = 0; d < CoordinationNumber; ++d) {
        const std::uint32_t a = neighborSite[firstBase+d];
        const std::uint32_t b = neighborSite[secondBase+d];
        env.firstNeighborSite[d] = a;
        env.secondNeighborSite[d] = b;
        env.firstNeighborSpecies[d] = static_cast<int>(siteSpecies[a]);
        env.secondNeighborSpecies[d] = static_cast<int>(siteSpecies[b]);
    }
}

void SystemKMCActiveFilteredBase::InitializeAcceptanceState() {
}

void SystemKMCActiveFilteredBase::UpdateAcceptanceStateAfterExchange(
    std::uint32_t, std::uint32_t, int, int) {
}

void SystemKMCActiveFilteredBase::ValidateAcceptanceStateOrAbort(const char*) {
}

void SystemKMCActiveFilteredBase::DetermineBonds() {
    // Build only the structures used by the Active-Filtered engine.
    BuildFixedTopology();
    InitializeSiteSpeciesCache();
    InitializeFastActiveList();
    InitializeAcceptanceState();

    txtstream << "Active-Filtered initialization: B=" << BondNumber
              << " flat_keys=" << activeKeys.size();
    logging(txtstream.str());
}

// Classical FCC Kawasaki sampling performs N directed site-neighbour proposals
// per MCS, so each undirected NN bond is proposed 1/6 times per MCS in
// expectation. Selecting uniformly among B unlike undirected bonds therefore
// requires 6/B MCS per filtered selection to preserve that proposal frequency.
// See docs/literature.md for the filtered-pair normalization boundary.
double SystemKMCActiveFilteredBase::MCSIncrementPerAttempt() const {
    if (BondNumber == 0) return 0.0;
    return static_cast<double>(CanonicalDirectionCount) /
           static_cast<double>(BondNumber);
}

// One Active-Filtered attempt: uniform structural candidate selection,
// post-selection energetics, then bounded local maintenance only after an
// accepted exchange. A rejected Metropolis test leaves all maintained state
// unchanged.
void SystemKMCActiveFilteredBase::JumpAttempt() {
    if (activeKeys.empty()) return;

    const std::uint32_t key = activeKeys[
        RandomIndex(static_cast<unsigned long long>(activeKeys.size()))];
    const std::uint32_t firstSite = key / CanonicalDirectionCount;
    const int direction = static_cast<int>(key % CanonicalDirectionCount);
    const std::uint32_t secondSite =
        neighborSite[static_cast<std::size_t>(firstSite) * CoordinationNumber + direction];

#ifndef NDEBUG
    if (siteSpecies[firstSite] == siteSpecies[secondSite]) {
        std::cerr << "Active-Filtered: inactive key selected" << std::endl;
        std::exit(2);
    }
#endif

    ++BenchmarkProbabilityEvaluationCount;
    const double jumpProbability =
        PairAcceptanceProbability(firstSite, secondSite, direction);
    if (!RandomAccept(jumpProbability)) return;

    const int firstSpecies = static_cast<int>(siteSpecies[firstSite]);
    const int secondSpecies = static_cast<int>(siteSpecies[secondSite]);
    CollectAffectedBondKeys(firstSite, secondSite);

    ExchangeSites(fixedX[firstSite], fixedY[firstSite], fixedZ[firstSite],
                  fixedX[secondSite], fixedY[secondSite], fixedZ[secondSite],
                  firstSpecies, secondSpecies);

    // Derived cached descriptors see the pre-exchange species values explicitly.
    UpdateAcceptanceStateAfterExchange(firstSite, secondSite,
                                       firstSpecies, secondSpecies);
    siteSpecies[firstSite] = static_cast<unsigned char>(secondSpecies);
    siteSpecies[secondSite] = static_cast<unsigned char>(firstSpecies);

    ++NAccepted;
    ++BenchmarkActiveTableUpdateCount;

    for (std::uint32_t affectedKey : affectedBondScratch) {
        RefreshActiveKey(affectedKey);
    }
    BondNumber = static_cast<unsigned long long>(activeKeys.size());
}

void SystemKMCActiveFilteredBase::ValidateBaseStateOrAbort(const char* context) {
    unsigned long long bruteActive = 0;
    for (std::uint32_t site = 0;
         site < static_cast<std::uint32_t>(TotalAtoms); ++site) {
        const int realSpecies =
            (this->*GetSpecies)(fixedX[site], fixedY[site], fixedZ[site]);
        if (realSpecies != static_cast<int>(siteSpecies[site])) {
            std::cerr << "Active-Filtered species-cache mismatch at " << context
                      << " site=" << site << std::endl;
            std::exit(2);
        }
        const std::size_t base =
            static_cast<std::size_t>(site) * CoordinationNumber;
        for (int d = 0; d < CanonicalDirectionCount; ++d) {
            const std::uint32_t key =
                site * CanonicalDirectionCount + static_cast<std::uint32_t>(d);
            const std::uint32_t nb = neighborSite[base+d];
            const bool expected = siteSpecies[site] != siteSpecies[nb];
            const bool listed = activeRow[key] != InvalidRow32;
            if (expected != listed) {
                std::cerr << "Active-Filtered active-list mismatch at " << context
                          << " key=" << key << std::endl;
                std::exit(2);
            }
            if (expected) ++bruteActive;
        }
    }

    if (bruteActive != activeKeys.size() || bruteActive != BondNumber) {
        std::cerr << "Active-Filtered active-count mismatch at " << context
                  << " list=" << activeKeys.size()
                  << " BondNumber=" << BondNumber
                  << " brute=" << bruteActive << std::endl;
        std::exit(2);
    }
    for (std::uint32_t row = 0; row < activeKeys.size(); ++row) {
        const std::uint32_t key = activeKeys[row];
        if (activeRow[key] != row) {
            std::cerr << "Active-Filtered reverse-row mismatch at " << context
                      << " row=" << row << std::endl;
            std::exit(2);
        }
    }
    ValidateAcceptanceStateOrAbort(context);
}

void SystemKMCActiveFilteredBase::EvalSysStep(unsigned long long n) {
    if (SysEvalParam.find(" OptimizationAudit ") != pos ||
        SysEvalParam.find(" ActiveFilteredAudit ") != pos) {
        ValidateBaseStateOrAbort("evaluation checkpoint");
    }
    SystemKMCHomogenous::EvalSysStep(n);
}
