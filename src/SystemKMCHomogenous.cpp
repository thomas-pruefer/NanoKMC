#include "SystemKMCHomogenous.h"
SystemKMCHomogenous::SystemKMCHomogenous() {

}

SystemKMCHomogenous::~SystemKMCHomogenous() {

}

void SystemKMCHomogenous::InitSubSys() {
    Ea=Str2Double(InputParam("Ea", "nanokmc.in"));
    dev=Str2Double(InputParam("dev", "nanokmc.in"));
    clvl=Str2Double(InputParam("clvl", "nanokmc.in"));

    if (kT <= 0.0) {
        std::cerr << "kT must be > 0." << std::endl;
        std::exit(2);
    }
    if (clvl < 0.0 || clvl > 1.0) {
        std::cerr << "clvl must satisfy 0 <= clvl <= 1." << std::endl;
        std::exit(2);
    }

    txtstream << "Ea: " << Ea; logging(txtstream.str());
    txtstream << "kT: " << kT; logging(txtstream.str());
    txtstream << "clvl: " << clvl; logging(txtstream.str());
    txtstream << "dev: " << dev; logging(txtstream.str());
    txtstream << "Base model: homogeneous FCC nearest-neighbour exchange"; logging(txtstream.str());
}

double SystemKMCHomogenous::DistributionFunction(double, double, double, int l) {
    double dpart;

    if (l==0) {
        dpart=clvl;
    } else {
        dpart=1;
    }

    return dpart;
}


// Homogeneous FCC Active-Filtered clock:
// one selected undirected active A-B pair represents z/(2B)=6/B
// classical MCS for an FCC lattice with coordination z=12.

double SystemKMCHomogenous::MCSIncrementPerAttempt() const {
    if (BondNumber == 0) return 0.0;
    return 6.0 / static_cast<double>(BondNumber);
}

void SystemKMCHomogenous::EvalSysStep(unsigned long long) {

}
