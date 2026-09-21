#include "SystemKMCActiveFilteredBinaryNN.h"
SystemKMCActiveFilteredBinaryNN::SystemKMCActiveFilteredBinaryNN() {
}

SystemKMCActiveFilteredBinaryNN::~SystemKMCActiveFilteredBinaryNN() {
}

void SystemKMCActiveFilteredBinaryNN::InitializeNeighborCounts() {
    speciesOneNeighborCount.assign(static_cast<std::size_t>(TotalAtoms), 0U);
    for (std::uint32_t site = 0;
         site < static_cast<std::uint32_t>(TotalAtoms); ++site) {
        int count = 0;
        const std::size_t base =
            static_cast<std::size_t>(site) * CoordinationNumber;
        for (int d = 0; d < CoordinationNumber; ++d) {
            const std::uint32_t nb = neighborSite[base+d];
            count += (siteSpecies[nb] == 1) ? 1 : 0;
        }
        speciesOneNeighborCount[site] = static_cast<unsigned char>(count);
    }
}

void SystemKMCActiveFilteredBinaryNN::InitializeAcceptanceState() {
    if (NSpecies != 2) {
        std::cerr << "KMCActiveFilteredBinaryNN requires NSpecies=2. "
                  << "Use KMCActiveFilteredGeneric for other species counts."
                  << std::endl;
        std::exit(2);
    }

    classProbability[0] = 1.0;
    for (int c = 1; c < RateClassCount; ++c) {
        classProbability[c] = std::exp(-Ea / kT * static_cast<double>(c));
    }
    InitializeNeighborCounts();

    txtstream << "Active-Filtered BinaryNN descriptor: one cached species-1 "
                 "nearest-neighbour count per site; 8 precomputed Metropolis factors";
    logging(txtstream.str());
}

double SystemKMCActiveFilteredBinaryNN::PairAcceptanceProbability(
    std::uint32_t firstSite, std::uint32_t secondSite, int) {
    const int firstSpecies = siteSpecies[firstSite];
    const int secondSpecies = siteSpecies[secondSite];
    if (firstSpecies == secondSpecies ||
        firstSpecies < 0 || firstSpecies > 1 ||
        secondSpecies < 0 || secondSpecies > 1) {
        std::cerr << "KMCActiveFilteredBinaryNN: invalid active pair"
                  << std::endl;
        std::exit(2);
    }

    const int firstCount =
        static_cast<int>(speciesOneNeighborCount[firstSite]);
    const int secondCount =
        static_cast<int>(speciesOneNeighborCount[secondSite]);

    // Symmetric homogeneous benchmark.  The selected A-B bond and
    // the four shared FCC neighbours cancel.  The remaining exchange-energy
    // class is determined exactly by the two cached coordination counts.
    int delta;
    if (firstSpecies == 0) {
        delta = secondCount - firstCount + 1;
    } else {
        delta = firstCount - secondCount + 1;
    }
    const int rateClass = (delta > 0) ? delta : 0;
    if (rateClass < 0 || rateClass >= RateClassCount) {
        std::cerr << "KMCActiveFilteredBinaryNN: invalid Metropolis class "
                  << rateClass << " (delta=" << delta << ")" << std::endl;
        std::exit(2);
    }
    return classProbability[rateClass];
}

void SystemKMCActiveFilteredBinaryNN::UpdateNeighborCountsAfterExchange(
    std::uint32_t firstSite, std::uint32_t secondSite,
    int oldFirstSpecies, int oldSecondSpecies) {
    const int deltaFirst =
        (oldSecondSpecies == 1 ? 1 : 0) - (oldFirstSpecies == 1 ? 1 : 0);
    const int deltaSecond =
        (oldFirstSpecies == 1 ? 1 : 0) - (oldSecondSpecies == 1 ? 1 : 0);

    const std::size_t firstBase =
        static_cast<std::size_t>(firstSite) * CoordinationNumber;
    const std::size_t secondBase =
        static_cast<std::size_t>(secondSite) * CoordinationNumber;

    for (int d = 0; d < CoordinationNumber; ++d) {
        const std::uint32_t a = neighborSite[firstBase+d];
        const int va = static_cast<int>(speciesOneNeighborCount[a]) + deltaFirst;
        if (va < 0 || va > CoordinationNumber) {
            std::cerr << "KMCActiveFilteredBinaryNN: neighbour-count range error"
                      << std::endl;
            std::exit(2);
        }
        speciesOneNeighborCount[a] = static_cast<unsigned char>(va);

        const std::uint32_t b = neighborSite[secondBase+d];
        const int vb = static_cast<int>(speciesOneNeighborCount[b]) + deltaSecond;
        if (vb < 0 || vb > CoordinationNumber) {
            std::cerr << "KMCActiveFilteredBinaryNN: neighbour-count range error"
                      << std::endl;
            std::exit(2);
        }
        speciesOneNeighborCount[b] = static_cast<unsigned char>(vb);
    }
}

void SystemKMCActiveFilteredBinaryNN::UpdateAcceptanceStateAfterExchange(
    std::uint32_t firstSite, std::uint32_t secondSite,
    int oldFirstSpecies, int oldSecondSpecies) {
    UpdateNeighborCountsAfterExchange(firstSite, secondSite,
                                      oldFirstSpecies, oldSecondSpecies);
}

double SystemKMCActiveFilteredBinaryNN::ReferenceProbabilityFromExplicitCounts(
    std::uint32_t firstSite, std::uint32_t secondSite) const {
    int firstCount = 0;
    int secondCount = 0;
    const std::size_t firstBase =
        static_cast<std::size_t>(firstSite) * CoordinationNumber;
    const std::size_t secondBase =
        static_cast<std::size_t>(secondSite) * CoordinationNumber;
    for (int d = 0; d < CoordinationNumber; ++d) {
        firstCount += siteSpecies[neighborSite[firstBase+d]] == 1 ? 1 : 0;
        secondCount += siteSpecies[neighborSite[secondBase+d]] == 1 ? 1 : 0;
    }

    int delta;
    if (siteSpecies[firstSite] == 0) {
        delta = secondCount - firstCount + 1;
    } else {
        delta = firstCount - secondCount + 1;
    }
    const int rateClass = (delta > 0) ? delta : 0;
    return classProbability[rateClass];
}

void SystemKMCActiveFilteredBinaryNN::ValidateAcceptanceStateOrAbort(
    const char* context) {
    for (std::uint32_t site = 0;
         site < static_cast<std::uint32_t>(TotalAtoms); ++site) {
        int brute = 0;
        const std::size_t base =
            static_cast<std::size_t>(site) * CoordinationNumber;
        for (int d = 0; d < CoordinationNumber; ++d) {
            brute += siteSpecies[neighborSite[base+d]] == 1 ? 1 : 0;
        }
        if (brute != static_cast<int>(speciesOneNeighborCount[site])) {
            std::cerr << "KMCActiveFilteredBinaryNN: cached-neighbour mismatch at "
                      << context << " site=" << site << std::endl;
            std::exit(2);
        }
    }

    unsigned long long checked = 0;
    constexpr unsigned long long maxChecks = 4096ULL;
    for (std::uint32_t key : activeKeys) {
        if (checked >= maxChecks) break;
        const std::uint32_t firstSite = key / CanonicalDirectionCount;
        const int direction = static_cast<int>(key % CanonicalDirectionCount);
        const std::uint32_t secondSite =
            neighborSite[static_cast<std::size_t>(firstSite) *
                         CoordinationNumber + direction];
        const double cached = PairAcceptanceProbability(
            firstSite, secondSite, direction);
        const double reference = ReferenceProbabilityFromExplicitCounts(
            firstSite, secondSite);
        const double tol = 1e-14 * std::max(1.0, std::abs(reference));
        if (std::abs(cached-reference) > tol) {
            std::cerr << "KMCActiveFilteredBinaryNN: probability mismatch at "
                      << context << " key=" << key << std::endl;
            std::exit(2);
        }
        ++checked;
    }
}
