#include "SystemClass.h"
// NanoKMC public evaluation/output layer.
//
// The public surface is intentionally split into:
//   (1) simulation / morphology observables, and
//   (2) state / visualization exports.
//
// The public evaluator provides generic, solver-independent spatial profiles
// and state/visualization exports as self-describing outputs.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

std::string NormalizeAxis(std::string axis, const char* fallback = "X") {
    if (axis.empty() || axis == "-1") axis = fallback;
    const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(axis[0])));
    if (c == 'X' || c == 'Y' || c == 'Z') return std::string(1, c);
    return fallback;
}

void EnsureEvaluationDirectory(const std::string& name) {
    std::filesystem::create_directories(std::filesystem::path("evaluation") / name);
}

bool IsEmptyOrMissing(const std::filesystem::path& path) {
    std::error_code ec;
    return !std::filesystem::exists(path, ec) || std::filesystem::file_size(path, ec) == 0;
}

int AxisCoordinate(const std::string& axis, int x, int y, int z) {
    if (axis == "Y") return y;
    if (axis == "Z") return z;
    return x;
}

int AxisLength(const std::string& axis, int lx1, int ly1, int lz1) {
    if (axis == "Y") return ly1;
    if (axis == "Z") return lz1;
    return lx1;
}

void PerpendicularCoordinates(const std::string& axis, int x, int y, int z,
                              int& u, int& v) {
    if (axis == "Y") { u = x; v = z; return; }
    if (axis == "Z") { u = x; v = y; return; }
    u = y; v = z;
}

void PerpendicularLengths(const std::string& axis, int lx1, int ly1, int lz1,
                          int& lu, int& lv) {
    if (axis == "Y") { lu = lx1; lv = lz1; return; }
    if (axis == "Z") { lu = lx1; lv = ly1; return; }
    lu = ly1; lv = lz1;
}

} // namespace

void SystemClass::CalcClustDist(int s, int sn) {
    long ncid;
    int clusteratoms, xx, yy, zz, at, i, k, nclusterbonds;
    long *clusterbonds, *surfaceatoms;
    int *clusters;
    clusters = new int [TotalAtoms];
    surfaceatoms = new long [TotalAtoms];
    clusterbonds = new long [TotalAtoms];

    unsigned long long x_minus, y_minus, z_minus, x__plus, y__plus, z__plus, x__zero, y__zero, z__zero;
    atx = new int [TotalAtoms];
    atz = new int [TotalAtoms];
    aty = new int [TotalAtoms];

    // calculating clusters
    ncid=0;
    for(zz=0;zz<lz1;zz++) {
        for(yy=0;yy<ly1;yy++) {
            for(xx=0;xx<lx1;xx++) {
                if(((xx^yy^zz)&LongOne)==0) {
                    if ((this->*GetSpecies)(xx, yy, zz)==s) {
                        (this->*SetSpecies)(xx, yy, zz, sn);
                        kk      =  0;
                        at		=  0;
                        atx[0]	= xx;
                        aty[0]	= yy;
                        atz[0]	= zz;
                        do {
                            top(atx[at], aty[at], atz[at], s, sn);
                            at++;
                        } while(at<=kk);
                        clusters[ncid]=at;

                        for (i=0;i<at;i++) {
                            (this->*SetSpecies)(atx[i], aty[i], atz[i], s);
                        }

                        nclusterbonds=0;
                        surfaceatoms[ncid]=0;
                        for (i=0;i<at;i++) {
                            x_minus = (lx+atx[i]) & lx, y_minus = (ly+aty[i]) & ly, z_minus = (lz+atz[i]) & lz,
                            x__plus = ( atx[i]+1) & lx, y__plus = ( aty[i]+1) & ly, z__plus = ( atz[i]+1) & lz,
                            x__zero =   atx[i]        , y__zero =   aty[i]        , z__zero =   atz[i]        ;
                            bool hasSameSpeciesNeighbor = false;
                            if ((this->*GetSpecies)(x_minus,y__zero,z__plus)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if ((this->*GetSpecies)(x__zero,y__plus,z__plus)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if ((this->*GetSpecies)(x__plus,y__zero,z__plus)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if ((this->*GetSpecies)(x__zero,y_minus,z__plus)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if ((this->*GetSpecies)(x_minus,y__zero,z_minus)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if ((this->*GetSpecies)(x__zero,y_minus,z_minus)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if ((this->*GetSpecies)(x__zero,y__plus,z_minus)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if ((this->*GetSpecies)(x__plus,y__zero,z_minus)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if ((this->*GetSpecies)(x_minus,y_minus,z__zero)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if ((this->*GetSpecies)(x_minus,y__plus,z__zero)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if ((this->*GetSpecies)(x__plus,y_minus,z__zero)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if ((this->*GetSpecies)(x__plus,y__plus,z__zero)==s) {nclusterbonds++; hasSameSpeciesNeighbor=true;}
                            if (hasSameSpeciesNeighbor) {
                                surfaceatoms[ncid]++;
                            }
                        }

                        for (i=0;i<at;i++) {
                            (this->*SetSpecies)(atx[i], aty[i], atz[i], sn);
                        }

                        clusterbonds[ncid]=nclusterbonds;
                        ncid++;
                    }
                }
            }
        }
    }
    delete[] atx;
    delete[] aty;
    delete[] atz;

    // Calculate ClusterDistr
    for (i=0; i<1000; i++) {
        clusterDistr[i]=0;
        clusterbondsDistr[i]=0;
        surfaceatomsDistr[i]=0;
    }

    clusteratoms=0;
    for (i=0; i<ncid; i++) {
        clusteratoms=clusteratoms+clusters[i];
        for (k=0; k<999; k++) {
            if (clusters[i]<=10*(k+1)) {
                clusterDistr[k]++;

                clusterbondsDistr[k]=clusterbondsDistr[k]+clusterbonds[k];
                surfaceatomsDistr[k]=surfaceatomsDistr[k]+surfaceatoms[k];
                break;
            }
        }
        if (clusters[i]>9990) {
            clusterDistr[999]++;
            clusterbondsDistr[999]=clusterbondsDistr[999]+clusterbonds[999];
            surfaceatomsDistr[999]=surfaceatomsDistr[999]+surfaceatoms[999];
        }
    }

    delete[] surfaceatoms;
    delete[] clusters;
    delete[] clusterbonds;

}


void SystemClass::BenchmarkEvaluation(unsigned long long n) {

    std::fstream BenchmarkFile;
    std::ifstream ExistingBenchmarkFile("evaluation/Benchmark.csv");
    bool writeHeader = (ExistingBenchmarkFile.peek() == std::ifstream::traits_type::eof());
    ExistingBenchmarkFile.close();

    BenchmarkFile.open("evaluation/Benchmark.csv", std::ios::out | std::ios::app);
    BenchmarkFile << std::setprecision(15);

    const double coordinationNumber = 12.0;
    const double totalUndirectedBonds = 0.5 * coordinationNumber * static_cast<double>(TotalAtoms);

    // Independently recount unequal nearest-neighbor bonds at checkpoints.
    // All public solvers maintain BondNumber as this same structural quantity,
    // so the scan is a solver-independent consistency audit.
    unsigned long long unequalBondCount = 0;
    for (long atom=0; atom<TotalAtoms; atom++) {
        int x = xpr[atom];
        int y = ypr[atom];
        int z = zpr[atom];
        int xn = x;
        int yn = y;
        int zn = z;
        int s0 = (this->*GetSpecies)(x, y, z);

        for (int j=1; j<7; j++) {
            switch (j) {
                case  1: xn=((x+ 1) & lx); yn=( y         ); zn=((z+ 1) & lz); break;
                case  2: xn=( x         ); yn=((y+ly) & ly); zn=((z+ 1) & lz); break;
                case  3: xn=((x+lx) & lx); yn=( y         ); zn=((z+ 1) & lz); break;
                case  4: xn=( x         ); yn=((y+ 1) & ly); zn=((z+ 1) & lz); break;
                case  5: xn=((x+ 1) & lx); yn=((y+ 1) & ly); zn=( z         ); break;
                case  6: xn=((x+lx) & lx); yn=((y+ 1) & ly); zn=( z         ); break;
            }
            if (s0 != (this->*GetSpecies)(xn, yn, zn)) {
                unequalBondCount++;
            }
        }
    }

    double interfaceDensity = 0.0;
    if (totalUndirectedBonds > 0.0) {
        interfaceDensity = static_cast<double>(unequalBondCount) / totalUndirectedBonds;
    }

    double attemptsPerAccepted = 0.0;
    double probabilityEvaluationsPerAccepted = 0.0;
    double acceptedFraction = 0.0;
    double inactiveRejectedFraction = 0.0;
    double probabilityEvaluationFraction = 0.0;

    if (NAccepted > 0) {
        attemptsPerAccepted = static_cast<double>(NJumps) / static_cast<double>(NAccepted);
        probabilityEvaluationsPerAccepted = static_cast<double>(BenchmarkProbabilityEvaluationCount) / static_cast<double>(NAccepted);
    }
    if (NJumps > 0) {
        acceptedFraction = static_cast<double>(NAccepted) / static_cast<double>(NJumps);
        inactiveRejectedFraction = static_cast<double>(BenchmarkInactiveRejectedCount) / static_cast<double>(NJumps);
        probabilityEvaluationFraction = static_cast<double>(BenchmarkProbabilityEvaluationCount) / static_cast<double>(NJumps);
    }

    if (writeHeader) {
        BenchmarkFile << "StepCount"
                      << ";CalcDataStep"
                      << ";RecordCount"
                      << ";BondNumberInternal"
                      << ";UnequalBondCount"
                      << ";InterfaceDensityFCC"
                      << ";ActiveFractionFCC"
                      << ";NJumpsRecord"
                      << ";NAcceptedRecord"
                      << ";InactiveRejectedCountRecord"
                      << ";ProbabilityEvaluationCountRecord"
                      << ";ActiveTableUpdatesRecord"
                      << ";AttemptsPerAccepted"
                      << ";ProbabilityEvaluationsPerAccepted"
                      << ";AcceptedFraction"
                      << ";InactiveRejectedFraction"
                      << ";ProbabilityEvaluationFraction"
                      << ";TimeMsCumulative"
                      << ";CommonMCSExact"
                      << ";EvolutionWallSecondsCumulative"
                      << ";TotalAtoms"
                      << ";CoordinationNumber"
                      << ";TotalUndirectedBonds";
        for (int s=0; s<NSpecies; s++) {
            BenchmarkFile << ";AtomNumber" << s;
        }
        BenchmarkFile << std::endl;
    }

    BenchmarkFile << StepCount
                  << ";" << n
                  << ";" << RecordCount
                  << ";" << BondNumber
                  << ";" << unequalBondCount
                  << ";" << interfaceDensity
                  << ";" << interfaceDensity
                  << ";" << NJumps
                  << ";" << NAccepted
                  << ";" << BenchmarkInactiveRejectedCount
                  << ";" << BenchmarkProbabilityEvaluationCount
                  << ";" << BenchmarkActiveTableUpdateCount
                  << ";" << attemptsPerAccepted
                  << ";" << probabilityEvaluationsPerAccepted
                  << ";" << acceptedFraction
                  << ";" << inactiveRejectedFraction
                  << ";" << probabilityEvaluationFraction
                  << ";" << StepData[3]
                  << ";" << (static_cast<double>(StepCount) + MCStepAccount)
                  << ";" << EvolutionWallSecondsCumulative
                  << ";" << TotalAtoms
                  << ";" << coordinationNumber
                  << ";" << totalUndirectedBonds;

    for (int s=0; s<NSpecies; s++) {
        BenchmarkFile << ";" << AtomNumber[s];
    }
    BenchmarkFile << std::endl;
    BenchmarkFile.close();
}


void SystemClass::EvalClustDist(unsigned long long n) {

    int k, s, sn;
    std::string loadfile;
    std::vector<std::fstream> surfaceatomsFile(NSpecies);
    std::vector<std::fstream> clusterbondsFile(NSpecies);

    // Cluster analysis is intentionally non-destructive. CalcClustDist()
    // temporarily rewrites species labels while traversing connected
    // components. The original lattice is restored in memory after each
    // species so evaluation does not depend on external archive round-trips and cannot leave
    // the live KMC lattice inconsistent if the reload fails. Keep an exact
    // in-memory copy of the packed lattice instead and restore from it.
    std::vector<unsigned long long> gitBackup;
    gitBackup.reserve(static_cast<size_t>(cubx) * static_cast<size_t>(cuby) * static_cast<size_t>(cubz));
    for (int gx=0; gx<cubx; ++gx) {
        for (int gy=0; gy<cuby; ++gy) {
            for (int gz=0; gz<cubz; ++gz) {
                gitBackup.push_back(git[gx][gy][gz]);
            }
        }
    }

    for (s=0;s<NSpecies;s++) {
        strcpy(st2, "evaluation/surfaceatoms_S");
        strcat(st2, Int2Str(s).c_str());
        strcat(st2, ".csv");
        surfaceatomsFile[s].open(st2, std::ios::out | std::ios::app);
    }

    for (s=0;s<NSpecies;s++) {
        strcpy(st2, "evaluation/clusterbonds_S");
        strcat(st2, Int2Str(s).c_str());
        strcat(st2, ".csv");
        clusterbondsFile[s].open(st2, std::ios::out | std::ios::app);
    }

    for (s=0; s<NSpecies; s++) {
        if (s==0) {
            sn=NSpecies-1;
        } else {
            sn=s-1;
        }

        // Calculate Cluster Distribution
        CalcClustDist(s, sn);

        // Plot ClusterDistribution for Species s
        strcpy(st2, GetBitfileName(n).c_str());
        ClusterDistributionFile[s] << st2;
        surfaceatomsFile[s] << st2;
        clusterbondsFile[s] << st2;
        for (k=0; k<1000; k++) {
            ClusterDistributionFile[s] << ";" << clusterDistr[k];
            surfaceatomsFile[s] << ";" << surfaceatomsDistr[k];
            clusterbondsFile[s] << ";" << clusterbondsDistr[k];
        }
        ClusterDistributionFile[s] << ";" << std::endl;
        surfaceatomsFile[s] << ";" << std::endl;
        clusterbondsFile[s] << ";" << std::endl;

        // Restore the exact pre-evaluation lattice in memory. Solver-specific
        // cached structures remain valid because the restored lattice is identical
        // to the pre-evaluation state.
        size_t gi = 0;
        for (int gx=0; gx<cubx; ++gx) {
            for (int gy=0; gy<cuby; ++gy) {
                for (int gz=0; gz<cubz; ++gz) {
                    git[gx][gy][gz] = gitBackup[gi++];
                }
            }
        }

    }

    for (s=0;s<NSpecies;s++) {
        surfaceatomsFile[s].close();
    }

    for (s=0;s<NSpecies;s++) {
        clusterbondsFile[s].close();
    }
}


void SystemClass::count_erase(unsigned long long aa, unsigned long long bb, unsigned long long cc, int s, int sn){
    if ((this->*GetSpecies)(aa, bb, cc)==s) {
        (this->*SetSpecies)(aa, bb, cc, sn);
        kk++;
        atx[kk]=aa;
        aty[kk]=bb;
        atz[kk]=cc;
    }
}


void SystemClass::top (unsigned long long x, unsigned long long y, unsigned long long z, int s, int sn) {

unsigned long long	x_minus = (lx+x) & lx, y_minus = (ly+y) & ly, z_minus = (lz+z) & lz,
		x__plus = ( x+1) & lx, y__plus = ( y+1) & ly, z__plus = ( z+1) & lz,
		x__zero =   x        , y__zero =   y        , z__zero =   z        ;

count_erase(x_minus,y__zero,z__plus,s,sn);count_erase(x__zero,y__plus,z__plus,s,sn);
count_erase(x__plus,y__zero,z__plus,s,sn);count_erase(x__zero,y_minus,z__plus,s,sn);
count_erase(x_minus,y__zero,z_minus,s,sn);count_erase(x__zero,y_minus,z_minus,s,sn);
count_erase(x__zero,y__plus,z_minus,s,sn);count_erase(x__plus,y__zero,z_minus,s,sn);
count_erase(x_minus,y_minus,z__zero,s,sn);count_erase(x_minus,y__plus,z__zero,s,sn);
count_erase(x__plus,y_minus,z__zero,s,sn);count_erase(x__plus,y__plus,z__zero,s,sn);

}


void SystemClass::WriteRasmol (std::string rfile) {

    logging("WriteRasmol: begin");

    std::vector<std::fstream> outputfile(NSpecies);
    std::fstream rasmolscript;
    std::vector<int> ReducedAtomNumber(NSpecies, 0);
    int i, s, f;
    const char* localdir = "Rasmol";

    // Public RasMol export writes all selected species at the checkpoint.
    const int r = 0;
    EnsureEvaluationDirectory(localdir);

    logging("WriteRasmol: opening XYZ files");
    for (s=0; s<NSpecies; s++) {
        strcpy(st2, "evaluation/");
        strcat(st2, localdir);
        strcat(st2, "/");
        strcat(st2, rfile.c_str());
        strcat(st2, "_S");
        strcat(st2, Int2Str(s).c_str());
        strcat(st2, ".xyz");
        outputfile[s].open(st2, std::ios::out);
        if (!outputfile[s].is_open()) {
            txtstream << "WriteRasmol: ERROR opening " << st2;
            logging(txtstream.str());
        }
    }

    logging("WriteRasmol: counting atoms");
    for (i=0; i<TotalAtoms; i++) {
        f=CheckFreedom (xpr[i], ypr[i], zpr[i]);
        if (f<(13-r) && f>=r) {
            const int species=(this->*GetSpecies)(xpr[i], ypr[i], zpr[i]);
            if (species < 0 || species >= NSpecies) {
                txtstream << "WriteRasmol: invalid species " << species
                          << " at atom " << i << " / " << TotalAtoms;
                logging(txtstream.str());
                continue;
            }
            ReducedAtomNumber[species]++;
        }
    }

    for (s=0; s<NSpecies; s++) {
        outputfile[s] << ReducedAtomNumber[s] << std::endl;
        outputfile[s] << std::endl;
    }

    logging("WriteRasmol: writing atoms");
    for (i=0; i<TotalAtoms; i++) {
        f=CheckFreedom (xpr[i], ypr[i], zpr[i]);
        if (f<(13-r) && f>=r) {
            s=(this->*GetSpecies)(xpr[i], ypr[i], zpr[i]);
            if (s < 0 || s >= NSpecies) {
                continue;
            }
            if (RasmolSpeciesPlot[s]==1) {
                outputfile[s] << SpeciesName[s] << "\t ";
                outputfile[s] << xpr[i] << "\t ";
                outputfile[s] << ypr[i] << "\t ";
                outputfile[s] << zpr[i] << "\t ";
                outputfile[s] << "0" << std::endl;
            }
        }
    }
    for (s=0; s<NSpecies; s++) {
        outputfile[s].close();
    }

    logging("WriteRasmol: writing script");
    strcpy(st2, "evaluation/");
    strcat(st2, localdir);
    strcat(st2, "/");
    strcat(st2, rfile.c_str());
    strcat(st2, ".rsm");
    rasmolscript.open(st2, std::ios::out);

    if (rasmolscript.is_open()) {
        rasmolscript << "background [255,255,255]" << std::endl;
        for (s=0; s<NSpecies; s++) {
            if (RasmolSpeciesPlot[s]==1) {
                rasmolscript << "load xyz " << rfile << "_S" << s << ".xyz" << std::endl;
                rasmolscript << "colour atoms " << SpeciesColor[s] << std::endl;
                rasmolscript << "spacefill 200" << std::endl;
                rasmolscript << "rotate z -90.00" << std::endl;
            }
        }
        rasmolscript << "rotate all" << std::endl;
        rasmolscript.close();
    } else {
        txtstream << "WriteRasmol: ERROR opening script " << st2;
        logging(txtstream.str());
    }

    logging("WriteRasmol: end");
}


void SystemClass::WriteOnefileRasmol (std::string rfile) {

    const int r = 1;

    std::fstream outputfile;
    std::fstream rasmolscript;
    int i, s, f, ReducedAtomNumber, NPlots;
    const char* st1 = "OnefileRasmol";
    EnsureEvaluationDirectory(st1);

    strcpy(st2, "evaluation/");
    strcat(st2, st1);
    strcat(st2, "/");
    strcat(st2, rfile.c_str());
    strcat(st2, ".xyz");
    outputfile.open(st2, std::ios::out);
    ReducedAtomNumber=0;
    NPlots=0;
    for (s=0; s<NSpecies; s++) {
        NPlots=NPlots+RasmolSpeciesPlot[s];
    }

    for (i=0;i<TotalAtoms;i++) {
        f=CheckFreedom (xpr[i], ypr[i], zpr[i]);
        s=(this->*GetSpecies)(xpr[i], ypr[i], zpr[i]);
        if (RasmolSpeciesPlot[s]==1) {
            if (f < (13-r) && f >= r) {
                        ReducedAtomNumber++;
            }
        }
    }
    outputfile << ReducedAtomNumber << std::endl;
    outputfile << std::endl;


    for (i=0;i<TotalAtoms;i++) {
        f=CheckFreedom (xpr[i], ypr[i], zpr[i]);
        if (f < (13-r) && f >= r) {
            s=(this->*GetSpecies)(xpr[i], ypr[i], zpr[i]);
            if (RasmolSpeciesPlot[s]==1) {
                    outputfile << SpeciesName[s] << "	 ";
                    outputfile << xpr[i] << "	 ";
                    outputfile << ypr[i] << "	 ";
                    outputfile << zpr[i] << "	 ";
                    if (NPlots==1) {
                        outputfile << 13-f;
                    } else {
                        outputfile << s;
                    }
                    outputfile << std::endl;
            }
        }
    }

    outputfile.close();

    strcpy(st2, "evaluation/");
    strcat(st2, st1);
    strcat(st2, "/");
    strcat(st2, rfile.c_str());
    strcat(st2, ".rsm");
    rasmolscript.open(st2, std::ios::out);

    rasmolscript << "background [255,255,255]" << std::endl;
    strcpy(st2, "load xyz ");
    strcat(st2, rfile.c_str());
    strcat(st2, ".xyz");
    rasmolscript << st2 << std::endl;
    rasmolscript << "spacefill 200" << std::endl;
    rasmolscript << "rotate z -90.00" << std::endl;
    rasmolscript << "rotate all" << std::endl;
    rasmolscript.close();
}


void SystemClass::WriteBlenderSimple (std::string rfile) {

    const int r = 1;

    std::fstream outputfile;
    std::fstream rasmolscript;
    int i, s, f, ReducedAtomNumber, NPlots;
    const char* st1 = "BlenderSimple";
    EnsureEvaluationDirectory(st1);

    strcpy(st2, "evaluation/");
    strcat(st2, st1);
    strcat(st2, "/");
    strcat(st2, rfile.c_str());
    strcat(st2, ".xyz");
    outputfile.open(st2, std::ios::out);
    ReducedAtomNumber=0;
    NPlots=0;
    for (s=0; s<NSpecies; s++) {
        NPlots=NPlots+RasmolSpeciesPlot[s];
    }

    for (i=0;i<TotalAtoms;i++) {
        f=CheckFreedom (xpr[i], ypr[i], zpr[i]);
        s=(this->*GetSpecies)(xpr[i], ypr[i], zpr[i]);
        if (RasmolSpeciesPlot[s]==1) {
            if ((f<=(12-r))) {
                        ReducedAtomNumber++;
            }
        }
    }
    outputfile << ReducedAtomNumber << std::endl;
    outputfile << std::endl;


    for (i=0;i<TotalAtoms;i++) {
        f=CheckFreedom (xpr[i], ypr[i], zpr[i]);
        if ((f<=(12-r)) ) {
            s=(this->*GetSpecies)(xpr[i], ypr[i], zpr[i]);
            if (RasmolSpeciesPlot[s]==1) {

                    outputfile << "Si" << "	 ";
                    outputfile << xpr[i] << "	 ";
                    outputfile << ypr[i] << "	 ";
                    outputfile << zpr[i] << "	 ";
                    if (NPlots==1) {
                        outputfile << 12-f;
                    } else {
                        outputfile << s;
                    }
                    outputfile << std::endl;
            }
        }
    }

    outputfile.close();

    strcpy(st2, "evaluation/");
    strcat(st2, st1);
    strcat(st2, "/");
    strcat(st2, rfile.c_str());
    strcat(st2, ".rsm");
    rasmolscript.open(st2, std::ios::out);

    rasmolscript << "background [255,255,255]" << std::endl;
    strcpy(st2, "load xyz ");
    strcat(st2, rfile.c_str());
    strcat(st2, ".xyz");
    rasmolscript << st2 << std::endl;
    rasmolscript << "spacefill 200" << std::endl;
    rasmolscript << "rotate z -90.00" << std::endl;
    rasmolscript << "rotate all" << std::endl;
    rasmolscript.close();
}


void SystemClass::WriteCSV (std::string rfile) {

    std::fstream outputfile;
    int i, s;
    const char* st1 = "CSV";
    EnsureEvaluationDirectory(st1);

    strcpy(st2, "evaluation/");
    strcat(st2, st1);
    strcat(st2, "/");
    strcat(st2, rfile.c_str());
    strcat(st2, ".csv");
    outputfile.open(st2, std::ios::out);

    for (i=0;i<TotalAtoms;i++) {
        s=(this->*GetSpecies)(xpr[i], ypr[i], zpr[i]);
        if (RasmolSpeciesPlot[s]==1) {
            outputfile << xpr[i] << ";";
            outputfile << ypr[i] << ";";
            outputfile << zpr[i] << ";";
            outputfile << s << ";";
            outputfile << std::endl;
        }
    }

    outputfile.close();
}


int SystemClass::CheckFreedom(int x, int y, int z) {
    int i, xn=0, yn=0, zn=0, j, s;

    j=0;
    s=(this->*GetSpecies)(x, y, z);

    for (i=0;i<12;i++) {

        switch(i){
        case  0: xn=((x+ 1) & lx); yn=  y; zn=((z+ 1) & lz); break;
        case  1: xn=  x; yn=((y+ly) & ly); zn=((z+ 1) & lz); break;
        case  2: xn=((x+lx) & lx); yn=  y; zn=((z+ 1) & lz); break;
        case  3: xn=  x; yn=((y+ 1) & ly); zn=((z+ 1) & lz); break;
        case  4: xn=((x+ 1) & lx); yn=((y+ 1) & ly); zn=  z; break;
        case  5: xn=((x+lx) & lx); yn=((y+ 1) & ly); zn=  z; break;
        case  6: xn=  x; yn=((y+ 1) & ly); zn=((z+lz) & lz); break;
        case  7: xn=((x+ 1) & lx); yn=  y; zn=((z+lz) & lz); break;
        case  8: xn=((x+lx) & lx); yn=  y; zn=((z+lz) & lz); break;
        case  9: xn=  x; yn=((y+ly) & ly); zn=((z+lz) & lz); break;
        case 10: xn=((x+ 1) & lx); yn=((y+ly) & ly); zn=  z; break;
        case 11: xn=((x+lx) & lx); yn=((y+ly) & ly); zn=  z; break;}
        if ((this->*GetSpecies)(xn, yn, zn)==s) {
            j++;
        }
    }
    return j;
}

// Species composition as a function of one configurable Cartesian coordinate.
// Bins use stored FCC half-grid indices; physical coordinates use lc/2.
void SystemClass::AxialCompositionProfile() {
    const std::string axis = NormalizeAxis(InputParam("AxialCompositionProfileAxis", "nanokmc.in"));
    const int length = AxisLength(axis, lx1, ly1, lz1);
    std::vector<long long> total(length, 0);
    std::vector<std::vector<long long>> counts(NSpecies, std::vector<long long>(length, 0));

    for (long i = 0; i < TotalAtoms; ++i) {
        const int q = AxisCoordinate(axis, xpr[i], ypr[i], zpr[i]);
        const int species = (this->*GetSpecies)(xpr[i], ypr[i], zpr[i]);
        if (q >= 0 && q < length && species >= 0 && species < NSpecies) {
            ++total[q];
            ++counts[species][q];
        }
    }

    EnsureEvaluationDirectory("AxialCompositionProfile");
    const std::filesystem::path path = std::filesystem::path("evaluation") /
        "AxialCompositionProfile" / ("axis_" + axis + ".csv");
    const bool writeHeader = IsEmptyOrMissing(path);
    std::ofstream out(path, std::ios::app);
    out << std::setprecision(17);
    if (writeHeader) {
        out << "checkpoint_mcs;axis;coordinate_index;coordinate;species;species_name;count;site_count;fraction\n";
    }
    for (int q = 0; q < length; ++q) {
        if (total[q] == 0) continue;
        for (int s = 0; s < NSpecies; ++s) {
            out << StepCount << ';' << axis << ';' << q << ';' << (0.5 * lc * q) << ';'
                << s << ';' << SpeciesName[s] << ';' << counts[s][q] << ';' << total[q] << ';'
                << (static_cast<double>(counts[s][q]) / static_cast<double>(total[q])) << '\n';
        }
    }
}

// Box-centred cylindrical composition map c_s(q,r).  This is a spatial
// composition profile, not a pair-correlation radial distribution function.
void SystemClass::CylindricalCompositionProfile() {
    const std::string axis = NormalizeAxis(InputParam("CylindricalCompositionProfileAxis", "nanokmc.in"));
    const int axialLength = AxisLength(axis, lx1, ly1, lz1);
    int lu = 0, lv = 0;
    PerpendicularLengths(axis, lx1, ly1, lz1, lu, lv);
    const double centerU = 0.5 * static_cast<double>(lu - 1);
    const double centerV = 0.5 * static_cast<double>(lv - 1);
    const int radialBins = std::max(1, static_cast<int>(std::floor(
        std::sqrt(centerU * centerU + centerV * centerV))) + 1);

    const std::size_t cells = static_cast<std::size_t>(axialLength) * radialBins;
    std::vector<long long> total(cells, 0);
    std::vector<std::vector<long long>> counts(NSpecies, std::vector<long long>(cells, 0));

    for (long i = 0; i < TotalAtoms; ++i) {
        const int axial = AxisCoordinate(axis, xpr[i], ypr[i], zpr[i]);
        int u = 0, v = 0;
        PerpendicularCoordinates(axis, xpr[i], ypr[i], zpr[i], u, v);
        const int radial = static_cast<int>(std::floor(std::sqrt(
            (u - centerU) * (u - centerU) + (v - centerV) * (v - centerV))));
        const int species = (this->*GetSpecies)(xpr[i], ypr[i], zpr[i]);
        if (axial < 0 || axial >= axialLength || radial < 0 || radial >= radialBins ||
            species < 0 || species >= NSpecies) continue;
        const std::size_t cell = static_cast<std::size_t>(axial) * radialBins + radial;
        ++total[cell];
        ++counts[species][cell];
    }

    EnsureEvaluationDirectory("CylindricalCompositionProfile");
    const std::filesystem::path path = std::filesystem::path("evaluation") /
        "CylindricalCompositionProfile" / ("axis_" + axis + ".csv");
    const bool writeHeader = IsEmptyOrMissing(path);
    std::ofstream out(path, std::ios::app);
    out << std::setprecision(17);
    if (writeHeader) {
        out << "checkpoint_mcs;axis;axial_index;axial_coordinate;radial_bin;radial_coordinate;species;species_name;count;site_count;fraction\n";
    }
    for (int axial = 0; axial < axialLength; ++axial) {
        for (int radial = 0; radial < radialBins; ++radial) {
            const std::size_t cell = static_cast<std::size_t>(axial) * radialBins + radial;
            if (total[cell] == 0) continue;
            for (int s = 0; s < NSpecies; ++s) {
                out << StepCount << ';' << axis << ';' << axial << ';' << (0.5 * lc * axial) << ';'
                    << radial << ';' << (0.5 * lc * (radial + 0.5)) << ';' << s << ';'
                    << SpeciesName[s] << ';' << counts[s][cell] << ';' << total[cell] << ';'
                    << (static_cast<double>(counts[s][cell]) / static_cast<double>(total[cell])) << '\n';
            }
        }
    }
}

// Box-centred spherical-shell composition profile.  Shells are defined by
// Euclidean distance in stored half-grid coordinates and reported in lc/2 units.
void SystemClass::SphericalCompositionProfile() {
    const double cx = 0.5 * static_cast<double>(lx1 - 1);
    const double cy = 0.5 * static_cast<double>(ly1 - 1);
    const double cz = 0.5 * static_cast<double>(lz1 - 1);
    const int radialBins = std::max(1, static_cast<int>(std::floor(
        std::sqrt(cx * cx + cy * cy + cz * cz))) + 1);

    std::vector<long long> total(radialBins, 0);
    std::vector<std::vector<long long>> counts(NSpecies, std::vector<long long>(radialBins, 0));

    for (long i = 0; i < TotalAtoms; ++i) {
        const double dx = xpr[i] - cx;
        const double dy = ypr[i] - cy;
        const double dz = zpr[i] - cz;
        const int radial = static_cast<int>(std::floor(std::sqrt(dx * dx + dy * dy + dz * dz)));
        const int species = (this->*GetSpecies)(xpr[i], ypr[i], zpr[i]);
        if (radial < 0 || radial >= radialBins || species < 0 || species >= NSpecies) continue;
        ++total[radial];
        ++counts[species][radial];
    }

    EnsureEvaluationDirectory("SphericalCompositionProfile");
    const std::filesystem::path path = std::filesystem::path("evaluation") /
        "SphericalCompositionProfile" / "profile.csv";
    const bool writeHeader = IsEmptyOrMissing(path);
    std::ofstream out(path, std::ios::app);
    out << std::setprecision(17);
    if (writeHeader) {
        out << "checkpoint_mcs;radial_bin;radial_coordinate;species;species_name;count;site_count;fraction\n";
    }
    for (int radial = 0; radial < radialBins; ++radial) {
        if (total[radial] == 0) continue;
        for (int s = 0; s < NSpecies; ++s) {
            out << StepCount << ';' << radial << ';' << (0.5 * lc * (radial + 0.5)) << ';'
                << s << ';' << SpeciesName[s] << ';' << counts[s][radial] << ';' << total[radial] << ';'
                << (static_cast<double>(counts[s][radial]) / static_cast<double>(total[radial])) << '\n';
        }
    }
}

// XZ composition projection obtained by integrating lattice occupancy over Y.
// Each XZ bin reports species counts and fractions over all valid Y sites.
void SystemClass::ProjectedCompositionXZ() {
    const std::size_t cells = static_cast<std::size_t>(lx1) * lz1;
    std::vector<long long> total(cells, 0);
    std::vector<std::vector<long long>> counts(NSpecies, std::vector<long long>(cells, 0));

    for (long i = 0; i < TotalAtoms; ++i) {
        const int x = xpr[i];
        const int z = zpr[i];
        const int species = (this->*GetSpecies)(xpr[i], ypr[i], zpr[i]);
        if (x < 0 || x >= lx1 || z < 0 || z >= lz1 || species < 0 || species >= NSpecies) continue;
        const std::size_t cell = static_cast<std::size_t>(x) * lz1 + z;
        ++total[cell];
        ++counts[species][cell];
    }

    EnsureEvaluationDirectory("ProjectedCompositionXZ");
    const std::filesystem::path path = std::filesystem::path("evaluation") /
        "ProjectedCompositionXZ" / "projection.csv";
    const bool writeHeader = IsEmptyOrMissing(path);
    std::ofstream out(path, std::ios::app);
    out << std::setprecision(17);
    if (writeHeader) {
        out << "checkpoint_mcs;x_index;x_coordinate;z_index;z_coordinate;species;species_name;count;site_count;fraction\n";
    }
    for (int x = 0; x < lx1; ++x) {
        for (int z = 0; z < lz1; ++z) {
            const std::size_t cell = static_cast<std::size_t>(x) * lz1 + z;
            if (total[cell] == 0) continue;
            for (int s = 0; s < NSpecies; ++s) {
                out << StepCount << ';' << x << ';' << (0.5 * lc * x) << ';'
                    << z << ';' << (0.5 * lc * z) << ';' << s << ';' << SpeciesName[s] << ';'
                    << counts[s][cell] << ';' << total[cell] << ';'
                    << (static_cast<double>(counts[s][cell]) / static_cast<double>(total[cell])) << '\n';
            }
        }
    }
}
