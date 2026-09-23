#include "SystemClass.h"

#include <charconv>
#include <cstdlib>
#include <system_error>

std::mt19937_64 SystemClass::RandomEngine;
bool SystemClass::RandomEngineSeeded = false;
std::uint64_t SystemClass::RandomEngineSeed = 0;

void SystemClass::InitializeRandom(std::uint64_t seed) {
    // A NanoKMC process normally uses one Seed value.  Keep the engine state
    // across System object re-creation so segmented calculations continue the
    // same random stream.
    if (!RandomEngineSeeded || RandomEngineSeed != seed) {
        RandomEngine.seed(seed);
        RandomEngineSeed = seed;
        RandomEngineSeeded = true;
    }
}

std::uint64_t SystemClass::RandomIndex(std::uint64_t upperExclusive) {
    if (upperExclusive <= 1) return 0;

    // Rejection removes modulo bias while retaining one raw RNG call in the
    // overwhelmingly common case.  Unsigned wrap-around makes threshold equal
    // to 2^64 mod upperExclusive.
    const unsigned long long threshold = (0ULL - upperExclusive) % upperExclusive;
    unsigned long long value;
    do {
        value = RandomEngine();
    } while (value < threshold);
    return value % upperExclusive;
}

double SystemClass::RandomUnit() {
    // 53 random bits fill the precision of an IEEE-754 double and guarantee
    // 0 <= u < 1.  This is used for distribution functions that need one
    // shared uniform variate across several cumulative species thresholds.
    return static_cast<double>(RandomEngine() >> 11) * (1.0 / 9007199254740992.0);
}

bool SystemClass::RandomAccept(double probability) {
    if (probability >= 1.0) return true;
    if (probability <= 0.0) return false;

    // The probability is already a double from the scientific model. Scaling
    // by exactly 2^64 maps that binary probability to the corresponding number
    // of accepted 64-bit RNG outcomes.  This is an integer-threshold Metropolis
    // test: one raw RNG call and one integer comparison in the hot path.
    const double scaled = probability * 18446744073709551616.0; // 2^64
    const unsigned long long threshold = static_cast<unsigned long long>(scaled);
    if (threshold == 0) return false;
    return RandomEngine() < threshold;
}

SystemClass::SystemClass() {
    BenchmarkProbabilityEvaluationCount = 0;
    BenchmarkInactiveRejectedCount = 0;
    BenchmarkActiveTableUpdateCount = 0;
}

SystemClass::~SystemClass() {

    int i, j, s;

    delete[] xpr;
    delete[] ypr;
    delete[] zpr;
    delete[] npr;
    delete[] StepData;
    delete[] AtomNumber;

    SCIFile.close();
    if (EvalParam.find(" ClusterDistribution ")!=pos) {
        for (s=0;s<NSpecies;s++) {
            ClusterDistributionFile[s].close();
        }
    }
    Logfile.close();
    for (i=0; i<cubx; i++) {
		for (j=0; j<cuby; j++) {
            delete[] git[i][j];
		}
		delete[] git[i];
	}
    delete[] git;

    for (i=0; i<lx1; i++) {
        for (j=0; j<ly1; j++) {
            delete[] xyzpointer[i][j];
        }
        delete[] xyzpointer[i];
    }
    delete[] xyzpointer;
    if (EvalParam.find(" ClusterDistribution ")!=pos) {
        delete[] ClusterDistributionFile;
    }
    delete[] SpeciesName;
    delete[] SpeciesColor;
    delete[] RasmolSpeciesPlot;

}

void SystemClass::Init(long long* SystemInit, int ColumnParam) {

    int i, j, k, s;
    std::string txt;

    time1=0;
    tstart = clock();
    MCStepAccount = 0.0;
    EvolutionWallSecondsCumulative = 0.0;
    // Initialize system
    StepCount=SystemInit[0];
    BondNumber=SystemInit[1];
    RecordCount=SystemInit[2];
    StepData=new long long[ColumnParam];
    for (i=0; i<ColumnParam; i++) {
        StepData[i]=SystemInit[i];
    }

    // Open logfile
    Logfile.open("log.dat", std::ios::out | std::ios::app);
    logging("===========================================");
    logging("=========  System Initialization  =========");
    txtstream << "=========  BondNumber: " << BondNumber ; logging(txtstream.str());
    txtstream << "=========  StepCount: " << StepCount ; logging(txtstream.str());
    txtstream << "=========  RecordCount: " << RecordCount ; logging(txtstream.str());
    logging("===========================================");

    // load system settings
    logging("load system settings...");
    InputSystem();

    // Set species number
    AtomNumber=new long [NSpecies];
    if (NSpecies <= 2) {
        SetSpecies=&SystemClass::SetSpecies2S;
        GetSpecies=&SystemClass::GetSpecies2S;
    } else if (NSpecies <= 4) {
        SetSpecies=&SystemClass::SetSpecies4S;
        GetSpecies=&SystemClass::GetSpecies4S;
    } else {
        // InputSystem() validates NSpecies <= 16.
        SetSpecies=&SystemClass::SetSpecies16S;
        GetSpecies=&SystemClass::GetSpecies16S;
    }

    // calculate dimensions
    lx = (LongOne << nx) - LongOne;
    ly = (LongOne << ny) - LongOne;
    lz = (LongOne << nz) - LongOne;
    lx1 = lx + LongOne;
    ly1 = ly + LongOne;
    lz1 = lz + LongOne;
    cubx = LongOne << (nx-LongTwo);
    cuby = LongOne << (ny-LongTwo);
    cubz = LongOne << (nz-LongTwo);
    TotalAtoms = static_cast<long>(1ULL << (nx + ny + nz - 1));
    txtstream << "nx: \t " << nx; logging(txtstream.str());
    txtstream << "ny: \t " << ny; logging(txtstream.str());
    txtstream << "nz: \t " << nz; logging(txtstream.str());
    txtstream << "lx: \t " << lx; logging(txtstream.str());
    txtstream << "ly: \t " << ly; logging(txtstream.str());
    txtstream << "lz: \t " << lz; logging(txtstream.str());

    // Initialize Subsystem
    logging("Initialize Subsystem...");
    InitSubSys();

    if (NSpecies>4) {
        cubx=cubx*2;
    }

    // Initialize arrays and pointers
	logging("Initialize arrays and pointers...");
    git = new unsigned long long **[cubx];
	for (i=0; i<cubx; i++) {
		git[i] = new unsigned long long *[cuby];
		for (j=0; j<cuby; j++) {
		    git[i][j] = new unsigned long long [cubz];
		}
	}

    xyzpointer = new unsigned long ***[lx1];
    for (i=0; i<lx1; i++) {
        xyzpointer[i] = new unsigned long **[ly1];
        for (j=0; j<ly1; j++) {
            xyzpointer[i][j] = new unsigned long *[lz1];
        }
    }

	for(i=0;i<cubx;i++) {
		for(j=0;j<cuby;j++) {
			for(k=0;k<cubz;k++) {
				git[i][j][k]=0;
			}
		}
	}

    logging("Initialize bitfile...");
    if (BondNumber>0) {

        // Loading existing bitfile
        LoadBitfile(GetBitfileName(StepCount));

        CountAtoms();

    } else {

        // Generating bitfile
        DistributeAtoms ();

    }

    txtstream << "StepCount: " << StepCount; logging(txtstream.str());
    txtstream << "TotalAtoms: " << TotalAtoms; logging(txtstream.str());

    long long SumAtoms=0;
    for (i=0;i<NSpecies;i++) {
        SumAtoms=SumAtoms+AtomNumber[i];
    }
    txtstream << "SumAtoms: " << SumAtoms; logging(txtstream.str());

    // Solver-specific precomputation.
    CodeCalculation();

    // Determine number of active bonds and active bonds array
    logging("Determine bonds...");
    DetermineBonds();
    StepData[1]=BondNumber;
    txtstream << "BondNumber: " << BondNumber; logging(txtstream.str());

    // other
    SCIFile.open("evaluation/eval.sce", std::ios::out | std::ios::app);

    // Initialize evaluation
    if (EvalParam.find(" ClusterDistribution ")!=pos) {
        ClusterDistributionFile = new std::fstream[NSpecies];
        for (s=0;s<NSpecies;s++) {
            strcpy(st2, "evaluation/clusters_S");
            strcat(st2, Int2Str(s).c_str());
            strcat(st2, ".csv");
            ClusterDistributionFile[s].open(st2, std::ios::out | std::ios::app);
        }
    }

}


void SystemClass::CodeCalculation() {
    // Optional solver-specific precomputation hook.
}

void SystemClass::SetSpecies16S (int x, int y, int z, unsigned long long s) {
    git[xyz2cubx16S(x)][xyz2cubyz16S(y)][xyz2cubyz16S(z)]=(~(Long15<<(xyz2bi16S(x, y, z))) & git[xyz2cubx16S(x)][xyz2cubyz16S(y)][xyz2cubyz16S(z)]);
    if (s>0) {
        git[xyz2cubx16S(x)][xyz2cubyz16S(y)][xyz2cubyz16S(z)]=((s<<(xyz2bi16S(x, y, z))) | git[xyz2cubx16S(x)][xyz2cubyz16S(y)][xyz2cubyz16S(z)]);
    }
}

void SystemClass::SetSpecies4S (int x, int y, int z, unsigned long long s) {
    switch (s) {
        case   0: git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]=(~(LongThree<<(xyz2bi4S(x, y, z))) & git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]); break;
        case   1: git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]=((LongOne<<(xyz2bi4S(x, y, z))) | git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]);
                  git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]=(~(LongOne<<(1+xyz2bi4S(x, y, z))) & git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]);break;
        case   2: git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]=(~(LongOne<<(xyz2bi4S(x, y, z))) & git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]);
                  git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]=((LongOne<<(1+xyz2bi4S(x, y, z))) | git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]);break;
        case   3: git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]=((LongThree<<(xyz2bi4S(x, y, z))) | git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)]); break;}
}

void SystemClass::SetSpecies2S (int x, int y, int z, unsigned long long s) {
    switch (s) {
        case   0: git[xyz2cubxyz2S(x)][xyz2cubxyz2S(y)][xyz2cubxyz2S(z)]=(git[xyz2cubxyz2S(x )][xyz2cubxyz2S(y )][xyz2cubxyz2S(z )]&(~(LongOne<<xyz2bi2S( x, y, z)))); break;
        case   1: git[xyz2cubxyz2S(x)][xyz2cubxyz2S(y)][xyz2cubxyz2S(z)]=(git[xyz2cubxyz2S(x )][xyz2cubxyz2S(y )][xyz2cubxyz2S(z )]|((LongOne<<xyz2bi2S( x, y, z)))); break;}
}

void SystemClass::InputSystem() {
    int s;

    const std::uint64_t seedInput = Str2UInt64(InputParam("Seed", "nanokmc.in"));
    InitializeRandom(seedInput);
    txtstream << "Random engine: mt19937_64, Seed: " << seedInput; logging(txtstream.str());
    kT=Str2Double(InputParam("kT", "nanokmc.in"));
    lc=Str2Double(InputParam("lc", "nanokmc.in"));
    EvalParam=InputParam("EvalParam", "nanokmc.in");
    SysEvalParam=InputParam("SysEvalParam", "nanokmc.in");
    NSpecies=Str2Int(InputParam("NSpecies", "nanokmc.in"));
    nx=Str2Int(InputParam("knx", "nanokmc.in"));
    ny=Str2Int(InputParam("kny", "nanokmc.in"));
    nz=Str2Int(InputParam("knz", "nanokmc.in"));
    // Backward compatibility with older example files that used nx/ny/nz.
    if (nx < 0) nx=Str2Int(InputParam("nx", "nanokmc.in"));
    if (ny < 0) ny=Str2Int(InputParam("ny", "nanokmc.in"));
    if (nz < 0) nz=Str2Int(InputParam("nz", "nanokmc.in"));

    if (NSpecies < 2 || NSpecies > 16) {
        std::cerr << "NSpecies must be between 2 and 16 for the current lattice "
                     "storage implementation." << std::endl;
        std::exit(2);
    }
    if (nx < 2 || ny < 2 || nz < 2) {
        std::cerr << "knx, kny and knz must each be >= 2 for the packed FCC "
                     "lattice representation." << std::endl;
        std::exit(2);
    }
    if (nx >= std::numeric_limits<int>::digits ||
        ny >= std::numeric_limits<int>::digits ||
        nz >= std::numeric_limits<int>::digits) {
        std::cerr << "knx, kny and knz are too large for the integer coordinate "
                     "representation on this platform." << std::endl;
        std::exit(2);
    }
    const int totalAtomExponent = nx + ny + nz - 1;
    if (totalAtomExponent >= std::numeric_limits<long>::digits) {
        std::cerr << "Requested FCC lattice is too large for TotalAtoms on this "
                     "platform." << std::endl;
        std::exit(2);
    }

    txtstream << "NSpecies: \t" << NSpecies; logging(txtstream.str());
    SpeciesName= new std::string [NSpecies];
    for (s=0; s<NSpecies; s++) {
        strcpy(st2, "SpeciesName(");
        strcat(st2, Int2Str(s+1).c_str());
        strcat(st2, ")");
        SpeciesName[s]=InputParam(st2, "nanokmc.in");
        txtstream << "Species " << s << ": \t" << SpeciesName[s]; logging(txtstream.str());
    }

    txtstream << "RasmolSpeciesPlot "; logging(txtstream.str());
	RasmolSpeciesPlot = new int [NSpecies];
    for (s=0; s<NSpecies; s++) {
        strcpy(st2, "RasmolSpeciesPlot(");
        strcat(st2, Int2Str(s+1).c_str());
        strcat(st2, ")");
        RasmolSpeciesPlot[s]=Str2Int(InputParam(st2, "nanokmc.in"));
        txtstream << "RasmolSpeciesPlot Species " << s << ": \t" << RasmolSpeciesPlot[s]; logging(txtstream.str());
    }

    SpeciesColor = new std::string [NSpecies];
    for (s=0; s<NSpecies; s++) {
        strcpy(st2, "SpeciesColor(");
        strcat(st2, Int2Str(s+1).c_str());
        strcat(st2, ")");
        SpeciesColor[s]=InputParam(st2, "nanokmc.in");
    }

    txtstream << "nx: \t" << nx; logging(txtstream.str());
    txtstream << "ny: \t" << ny; logging(txtstream.str());
    txtstream << "nz: \t" << nz; logging(txtstream.str());
    txtstream << "kT: \t" << kT; logging(txtstream.str());
    txtstream << "lc: \t" << lc; logging(txtstream.str());
    txtstream << "EvalParam: \t" << EvalParam; logging(txtstream.str());
    txtstream << "SysEvalParam: \t" << SysEvalParam; logging(txtstream.str());


}

std::string SystemClass::GetBitfileName(unsigned long n) {
    char s[255];
    char jst[50];
    sprintf(jst, "%lu", n);
    strcpy(s, "");

    if(n==0) strcat(s, "00000000");
    else {
        double lo = log10((double)(n));
        int nullen = 1+(int)(lo);
        int maxnullen = 8 - nullen;
        for(int zahl=1; zahl<=maxnullen; zahl++) strcat(s, "0");
        strcat(s, jst);
    }
    return s;
}

void SystemClass::CountAtoms() {
    int i, j, k, l, n;

    // make xyz arrays
    xpr = new int [TotalAtoms];
    ypr = new int [TotalAtoms];
    zpr = new int [TotalAtoms];
    npr = new unsigned long [TotalAtoms];
    for (i=0;i<NSpecies;i++) {
        AtomNumber[i]=0;
    }

    n=0;
    for (i=0; i<lx1; i++){
        for (j=0; j<ly1; j++){
            for (k=0; k<lz1; k++){
                if (((i ^ j ^ k) & LongOne)==0){
                    npr[n]=n;
                    xpr[n]=i;
                    ypr[n]=j;
                    zpr[n]=k;
                    l=(this->*GetSpecies)(i, j, k);
                    AtomNumber[l]++;
                    xyzpointer[i][j][k]=&npr[n];
                    n++;
                }
            }
        }
    }

    for (i=0; i<NSpecies; i++) {
        StepData[5+i]=AtomNumber[i];
    }

    txtstream << "Counted atoms: "; logging(txtstream.str());
    for (i=0; i<NSpecies; i++) {
        txtstream << "Species " << i << ": \t " << AtomNumber[i]; logging(txtstream.str());
    }
}

void SystemClass::LoadBitfile(std::string Filename) {
    int i,j,k;
    char cdln1[255];
    std::fstream bitfile;

    // unpack file
    strcpy(cdln1, "unzip -o -q output/bit/");
    strcat(cdln1, Filename.c_str());
    strcat(cdln1, ".zip ");
    system(cdln1);

    // load bitfile
    strcpy(cdln1, Filename.c_str());
    strcat(cdln1, ".dat");
    bitfile.open(cdln1, std::ios::in | std::ios::binary);
    for (i=0; i<cubx; i++){
        for (j=0; j<cuby; j++){
            for (k=0; k<cubz; k++){
                bitfile.read((char *) &git[i][j][k], 8);
            }
        }
    }
    bitfile.close();

    // delete .dat file
    #ifdef _WIN32
        strcpy(cdln1, "del " );
    #else
        strcpy(cdln1, "rm " );
    #endif
    strcat(cdln1, Filename.c_str());
    strcat(cdln1, ".dat");
    system(cdln1);
    txtstream << "Bitfile loaded: \t " << cdln1; logging(txtstream.str());

}

void SystemClass::DistributeAtoms () {
    int i, j, k, l, n;
    double rnd, x, y, z;

    // make xyz arrays
    xpr = new int [TotalAtoms];
    ypr = new int [TotalAtoms];
    zpr = new int [TotalAtoms];
    npr = new unsigned long [TotalAtoms];
    for (i=0;i<NSpecies;i++) {
        AtomNumber[i]=0;
    }

    n=0;
    for (i=0; i<lx1; i++){
        for (j=0; j<ly1; j++){
            for (k=0; k<lz1; k++){
                if (((i ^ j ^ k) & LongOne)==0){
                    x=lc*i/2;
                    y=lc*j/2;
                    z=lc*k/2;
                    rnd=RandomUnit();
                    for (l=0; l<NSpecies; l++) {
                        if (rnd<DistributionFunction(x, y, z, l)) {
                            npr[n]=n;
                            xpr[n]=i;
                            ypr[n]=j;
                            zpr[n]=k;
                            AtomNumber[l]++;
                            (this->*SetSpecies)(i, j, k, l);
                            xyzpointer[i][j][k]=&npr[n];
                            n++;
                            break;
                        }

                    }
                }
            }
        }
    }

    txtstream << "Placed atoms: "; logging(txtstream.str());
    for (i=0; i<NSpecies; i++) {
        txtstream << "Species " << i << ": \t " << AtomNumber[i]; logging(txtstream.str());
    }

    SaveBitFile(0);
    //RunEval(StepData);
}

int SystemClass::xyz2nn(int i) {
    switch(i){
        case 0: return 8; break;
        case 1: return 6; break;
        case 2: return 7; break;
        case 3: return 9; break;
        case 4: return 11; break;
        case 5: return 10; break;
        case 6: return 1; break;
        case 7: return 2; break;
        case 8: return 0; break;
        case 9: return 3; break;
        case 10: return 5; break;
        case 11: return 4; break;
        case 12: return 8; break;
        case 13: return 6; break;
        case 14: return 7; break;
        case 15: return 9; break;
        case 16: return 11; break;
        case 17: return 10; break;
        case 18: return 1; break;
        case 19: return 2; break;
        case 20: return 0; break;
        case 21: return 3; break;
        case 22: return 5; break;
        case 23: return 4;
        default: return -1;
    }
}

/////////////////////////////////////////////////////////////////////////
void SystemClass::ExchangeSites(int x, int y, int z, int xn, int yn, int zn, int s, int sn) {
    long long iiin, iii;

    iiin=xyzpointer[xn][yn][zn][0];
    iii=xyzpointer[x][y][z][0];
    xpr[iiin]=x;
    ypr[iiin]=y;
    zpr[iiin]=z;
    xyzpointer[x][y][z]=&npr[iiin];

    (this->*SetSpecies)(x, y, z, sn);
    (this->*SetSpecies)(xn, yn, zn, s);

    xpr[iii]=xn;
    ypr[iii]=yn;
    zpr[iii]=zn;
    xyzpointer[xn][yn][zn]=&npr[iii];
}

void SystemClass::RunMC(long long* in) {
    int i;
    const long long n = in[0];
    RecordCount=in[2];
    NJumps=0;
    NAccepted=0;
    BenchmarkProbabilityEvaluationCount=0;
    BenchmarkInactiveRejectedCount=0;
    BenchmarkActiveTableUpdateCount=0;

    // Clean benchmark timer: only time spent evolving the KMC state is counted.
    // Initialization, SaveBitFile, RunEval, morphology calculations and file output
    // outside this loop are deliberately excluded.
    const auto evolutionStart = std::chrono::steady_clock::now();

    // All public solvers advance the shared common-MCS coordinate through an
    // explicit per-attempt/event increment.
    unsigned long long cachedBondNumber = LongN;
    double cachedMCSIncrement = 0.0;

    while (StepCount < static_cast<unsigned long long>(n)) {
        // Recompute only when the state-dependent denominator changes.
        if (RecomputeMCSIncrementEveryAttempt() ||
            BondNumber != cachedBondNumber || cachedMCSIncrement <= 0.0) {
            cachedMCSIncrement = MCSIncrementPerAttempt();
            cachedBondNumber = BondNumber;
        }

        // A zero increment denotes an absorbing state for catalogue-based
        // solvers. The classical solver always returns 1/N.
        if (cachedMCSIncrement <= 0.0) {
            StepCount = static_cast<unsigned long long>(n);
            MCStepAccount = 0.0;
            break;
        }

        // The increment belongs to the pre-attempt state. JumpAttempt may
        // change BondNumber on acceptance.
        JumpAttempt();
        NJumps++;
        MCStepAccount += cachedMCSIncrement;

        while (MCStepAccount >= 1.0) {
            MCStepAccount -= 1.0;
            StepCount++;
            EvalStep();
        }
    }

    const auto evolutionEnd = std::chrono::steady_clock::now();
    EvolutionWallSecondsCumulative +=
        std::chrono::duration<double>(evolutionEnd - evolutionStart).count();

    RecordCount++;
    txtstream << "RunMC - StepCount: " << StepCount
              << " - CommonMCSExact: " << (static_cast<double>(StepCount) + MCStepAccount)
              << " - BondNumber: " << BondNumber
              << " - RecordCount: " << RecordCount; logging(txtstream.str());
    txtstream << "EvolutionWallSecondsCumulative: " << EvolutionWallSecondsCumulative;
    logging(txtstream.str());

    // Retain the process-CPU field in CalcData for compatibility.
    // It is NOT the clean benchmark timer because it starts in Init().
    time1 = clock() - tstart;
    time1 = time1/CLOCKS_PER_SEC;
    txtstream << "Legacy process Time : " << time1; logging(txtstream.str());

    // save results
    // 0: requested MCStep record
    // 1: BondNumber
    // 2: RecordCount
    // 3: legacy process time in ms (compatibility only)
    // 4: jump attempts during record
    // 5+: atoms of each species

    StepData[0]=n;
    StepData[1]=BondNumber;
    StepData[2]=RecordCount;
    StepData[3]=abs(time1*1000);
    StepData[4]=NJumps;
    for (i=0; i<NSpecies; i++) {
        StepData[5+i]=AtomNumber[i];
    }
    SaveBitFile(n);
}

void SystemClass::RunEval(long long* in) {
    const unsigned long long checkpoint = static_cast<unsigned long long>(in[0]);
    RecordCount=in[2];

    txtstream << "RunEval - StepCount: " << StepCount << " - BondNumber: " << BondNumber
              << " - RecordCount: " << RecordCount;
    logging(txtstream.str());

    // ---------------------------------------------------------------------
    // Simulation / morphology observables
    // ---------------------------------------------------------------------
    if (EvalParam.find(" Benchmark ")!=pos) {
        logging("Running Benchmark");
        BenchmarkEvaluation(checkpoint);
    }
    if (EvalParam.find(" ClusterDistribution ")!=pos) {
        logging("Running ClusterDistribution");
        EvalClustDist(checkpoint);
    }
    if (EvalParam.find(" AxialCompositionProfile ")!=pos) {
        logging("Running AxialCompositionProfile");
        AxialCompositionProfile();
    }
    if (EvalParam.find(" CylindricalCompositionProfile ")!=pos) {
        logging("Running CylindricalCompositionProfile");
        CylindricalCompositionProfile();
    }
    if (EvalParam.find(" SphericalCompositionProfile ")!=pos) {
        logging("Running SphericalCompositionProfile");
        SphericalCompositionProfile();
    }
    if (EvalParam.find(" ProjectedCompositionXZ ")!=pos) {
        logging("Running ProjectedCompositionXZ");
        ProjectedCompositionXZ();
    }

    // ---------------------------------------------------------------------
    // State / visualization exports
    // ---------------------------------------------------------------------
    const std::string checkpointName = GetBitfileName(checkpoint);
    if (EvalParam.find(" Rasmol ")!=pos) {
        logging("Running Rasmol");
        WriteRasmol(checkpointName);
    }
    if (EvalParam.find(" OnefileRasmol ")!=pos) {
        logging("Running OnefileRasmol");
        WriteOnefileRasmol(checkpointName);
    }
    if (EvalParam.find(" BlenderSimple ")!=pos) {
        logging("Running BlenderSimple");
        WriteBlenderSimple(checkpointName);
    }
    if (EvalParam.find(" WriteCSV ")!=pos) {
        logging("Running WriteCSV");
        WriteCSV(checkpointName);
    }

    EvalSysStep(checkpoint);
}

void SystemClass::SaveBitFile(unsigned long n){

    int i,j,k;
    char jst[50];
    char st2[255], cdln3[255];
    std::fstream bitfile;

    sprintf(jst, "%lu", n);
    strcpy(st2, "");

    if(n==0) strcat(st2, "00000000");

    else {
        double lo = log10((double)(n));
        int nullen = 1+(int)(lo);
        int maxnullen = 8 - nullen;
        for(int zahl=1; zahl<=maxnullen; zahl++) strcat(st2, "0");
        strcat(st2, jst);
    }

    strcpy(cdln3, "zip -o -q -m output/bit/");
    strcat(cdln3, st2);
    strcat(cdln3, ".zip ");

    strcat(cdln3, st2);
    strcat(cdln3, ".dat");

    strcat(st2, ".dat");

    bitfile.open(st2, std::ios::out | std::ios::binary);
    for (i=0; i<cubx; i++){
        for (j=0; j<cuby; j++){
            for (k=0; k<cubz; k++){
                bitfile.write((char *) &git[i][j][k], 8);
            }
        }
    }
    bitfile.close();
    txtstream << "Bitfile saved: " << st2; logging(txtstream.str());
    system(cdln3);
}

void SystemClass::InitSubSys() {

}

void SystemClass::EvalSysStep(unsigned long long) {

}

bool SystemClass::RecomputeMCSIncrementEveryAttempt() const {
    return false;
}

void SystemClass::EvalStep() {

}

long long* SystemClass::GetStepData() {
    return StepData;
}

double SystemClass::Str2Double (const std::string &str) {
    std::stringstream ss(str);
    double n;
    double m;
    ss >> n;
    if (str.find("D") != std::string::npos) {
        m=Str2Double(str.substr(str.find("D")+1,str.find("D")+2));
        n=n*pow(10,m);
    }
    return n;
}

int SystemClass::Str2Int (const std::string &str) {
    std::stringstream ss(str);
    int n;
    double m;
    ss >> n;
    if (str.find("D") != std::string::npos) {
        m=Str2Int(str.substr(str.find("D")+1,str.find("D")+2));
        n=n*pow(10,m);
    }
    return n;
}

std::uint64_t SystemClass::Str2UInt64(const std::string &str) {
    if (str.empty() || str.front() == '-') {
        std::cerr << "Expected a non-negative uint64_t value, got: " << str
                  << std::endl;
        std::exit(2);
    }

    std::uint64_t value = 0;
    const char* begin = str.data();
    const char* end = begin + str.size();
    const auto result = std::from_chars(begin, end, value, 10);
    if (result.ec != std::errc() || result.ptr != end) {
        std::cerr << "Invalid uint64_t value: " << str << std::endl;
        std::exit(2);
    }
    return value;
}

std::string SystemClass::Int2Str (int n) {

    std::ostringstream convert;

    convert << n;
    return convert.str();
}

std::string SystemClass::InputParam(std::string paramname, std::string filename) {
    std::string line;
    const std::string key = paramname + "=";

    std::fstream inputfile(filename.c_str(), std::ios::in);
    while (std::getline(inputfile, line)) {
        const std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        if (line.compare(first, 2, "//") == 0) continue;
        if (line.compare(first, key.size(), key) != 0) continue;

        const std::size_t valueBegin = first + key.size();
        std::size_t semi = line.find(";", valueBegin);
        if (semi == std::string::npos) semi = line.size();
        std::string param = line.substr(valueBegin, semi - valueBegin);

        const std::size_t valueFirst = param.find_first_not_of(" \t\r\n");
        const std::size_t valueLast = param.find_last_not_of(" \t\r\n");
        if (valueFirst == std::string::npos) {
            param.clear();
        } else {
            param = param.substr(valueFirst, valueLast - valueFirst + 1);
        }

        if (param.size() >= 2 && param.front() == '"' && param.back() == '"') {
            param = param.substr(1, param.size() - 2);
        }
        return param;
    }
    return "-1";
}

void SystemClass::Zip() {
    // Consolidate per-checkpoint packed lattice states into one archive.

    #ifdef _WIN32
        system("if exist output\\bit\\0*.zip zip -0 -o -q -m output/bit/data.zip output/bit/0*.zip");
    #else
        system("ls output/bit/0*.zip >/dev/null 2>&1 && zip -0 -o -q -m output/bit/data.zip output/bit/0*.zip");
    #endif
    logging("Bitfiles zipped");
    SysZip();
}

void SystemClass::SysZip() {

}

void SystemClass::logging(std::string txt) {
    Logfile << txt << std::endl;
    std::cout << txt << std::endl;
    txtstream.str(std::string());;
}
