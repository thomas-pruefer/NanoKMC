#include "SystemKMCHomogenous.h"
SystemKMCHomogenous::SystemKMCHomogenous() {

}

SystemKMCHomogenous::~SystemKMCHomogenous() {

}

void SystemKMCHomogenous::InitSubSys() {
    Ea=Str2Double(InputParam("Ea", "nanokmc.in"));
    dev=Str2Double(InputParam("dev", "nanokmc.in"));
    clvl=Str2Double(InputParam("clvl", "nanokmc.in"));

    txtstream << "Ea: " << Ea; logging(txtstream.str());
    txtstream << "kT: " << kT; logging(txtstream.str());
    txtstream << "clvl: " << clvl; logging(txtstream.str());
    txtstream << "dev: " << dev; logging(txtstream.str());
    txtstream << "Base model: homogeneous FCC Kawasaki benchmark"; logging(txtstream.str());
}

double SystemKMCHomogenous::DistributionFunction(double x, double y, double z, int l) {
    double dpart;

    if (l==0) {
        dpart=clvl;
    } else {
        dpart=1;
    }

    return dpart;
}


// Homogeneous Active-Filtered benchmark clock:
// one selected undirected active A-B pair represents z/(2B)=6/B
// classical MCS for an FCC lattice with coordination z=12.

double SystemKMCHomogenous::MCSIncrementPerAttempt() const {
    if (BondNumber == 0) return 0.0;
    return 6.0 / static_cast<double>(BondNumber);
}

void SystemKMCHomogenous::EvalSysStep(unsigned long long n) {

}
