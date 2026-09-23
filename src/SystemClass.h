#pragma once

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <chrono>
#include <iomanip>
#include <random>
#include <limits>
#include <cstdint>

class SystemClass {

    public:
        SystemClass();
        virtual ~SystemClass();
        void RunMC(long long* in);

        void Init(long long* SystemInit, int ColumnParam);
        void Zip();
        long long* GetStepData();
        void RunEval(long long* in);
        unsigned long long BondNumber;
        int lx1, ly1, lz1, lx, ly, lz, nx, ny, nz;
        unsigned long long NJumps;
        long long NAccepted;
        long long BenchmarkProbabilityEvaluationCount;
        long long BenchmarkInactiveRejectedCount;
        long long BenchmarkActiveTableUpdateCount;
        long TotalAtoms, *AtomNumber;
        // Compact packed lattice storage. Species are stored per fixed FCC site;
        // the packing density depends on NSpecies (1, 2 or 4 bits per site).
        unsigned long long *** git;
        unsigned long long StepCount;
        unsigned long long LongOne=1, LongTwo=2, LongThree=3, LongZero=0,
                           Long15=15, LongN=std::numeric_limits<unsigned long long>::max();
        long long *StepData;
        void logging(std::string txt);

    protected:
        virtual void JumpAttempt() = 0;
        virtual double MCSIncrementPerAttempt() const = 0;
        virtual bool RecomputeMCSIncrementEveryAttempt() const;
        // Shared, deterministic RNG for all NanoKMC systems.  The engine is
        // process-wide so re-creating a System object during a segmented run does
        // not restart the random stream.
        static std::mt19937_64 RandomEngine;
        static bool RandomEngineSeeded;
        static std::uint64_t RandomEngineSeed;
        void InitializeRandom(std::uint64_t seed);
        std::uint64_t RandomIndex(std::uint64_t upperExclusive);
        double RandomUnit();
        bool RandomAccept(double probability);
        virtual void SysZip();
        std::string GetBitfileName(unsigned long n);
        void LoadBitfile(std::string Filename);
        virtual void InitSubSys();
        virtual void DetermineBonds() = 0;
        virtual void DistributeAtoms();
        int xyz2nn(int i);
        void SaveBitFile(unsigned long n);
        void InputSystem();
        virtual void EvalStep();
        virtual void EvalSysStep(unsigned long long n);
        void EvalClustDist(unsigned long long n);
        void BenchmarkEvaluation(unsigned long long n);
        void count_erase(unsigned long long aa, unsigned long long bb, unsigned long long cc, int s, int sn);
        void top(unsigned long long x, unsigned long long y, unsigned long long z, int s, int sn);
        double Str2Double (const std::string &str);
        int Str2Int (const std::string &str);
        std::uint64_t Str2UInt64(const std::string &str);
        virtual double DistributionFunction(double x, double y, double z, int l) = 0;
        void WriteRasmol (std::string rfile);
        void WriteBlenderSimple (std::string rfile);
        void WriteOnefileRasmol (std::string rfile);
        void WriteCSV (std::string rfile);
        // Generic spatial-composition observables.  These routines operate on
        // the current lattice checkpoint and emit self-describing tidy CSVs.
        void AxialCompositionProfile();
        void CylindricalCompositionProfile();
        void SphericalCompositionProfile();
        void ProjectedCompositionXZ();
        void CountAtoms();
        void SetSpecies16S(int x, int y, int z, unsigned long long  s);
        void SetSpecies4S(int x, int y, int z, unsigned long long s);
        void SetSpecies2S(int x, int y, int z, unsigned long long s);
        void CalcClustDist(int s, int sn);
        std::string InputParam(std::string paramname, std::string filename);
        std::string Int2Str (int n);
        int CheckFreedom (int x, int y, int z);
        virtual void CodeCalculation();
        void ExchangeSites(int x, int y, int z, int xn, int yn, int zn, int s, int sn);

        int cubx, cuby, cubz;
        // Fixed-site coordinate tables and reverse site lookup. These preserve
        // the compact legacy lattice layout while solver code works with stable
        // site indices.
        int *xpr, *ypr, *zpr;
        unsigned long ****xyzpointer, *npr;
        std::string name, EvalParam, SysEvalParam, *SpeciesName;
        double kT, lc, Ea;
        std::fstream *ClusterDistributionFile, SCIFile, Logfile;
        int kk, *atx, *aty, *atz;
        char st2[255];
        int clusterDistr[1000];
        int clusterbondsDistr[1000];
        int surfaceatomsDistr[1000];
        int RecordCount;
        // Fractional remainder of the common MCS clock. StepCount stores whole MCS.
        double MCStepAccount;
        // Clean elapsed wall time spent only inside KMC evolution loops.
        double EvolutionWallSecondsCumulative;
        std::stringstream txtstream;
        double time1, tstart, dev;
        int NSpecies;
        void (SystemClass::*SetSpecies)(int, int, int, unsigned long long);
        int (SystemClass::*GetSpecies)(int, int, int);

        std::size_t pos=-1;
        std::string *SpeciesColor;
        int *RasmolSpeciesPlot;

        inline unsigned long long xyz2cubxyz2S(unsigned long long xyz){return (xyz>>LongTwo);}
        inline unsigned long long xyz2bi2S(unsigned long long x,unsigned long long y,unsigned long long z){return((x & 3)+4*(y & 3)+16*(z & 3));}
        inline unsigned long long xyz2cubxyz4S(unsigned long long xyz){return (xyz>>LongTwo);}
        inline unsigned long long xyz2bi4S(unsigned long long x,unsigned long long y,unsigned long long z){return((x & 3)-(((y ^ z) & LongOne))+4*(y & 3)+16*(z & 3));}

        inline unsigned long long xyz2cubyz16S(unsigned long long xyz){return (xyz>>LongTwo);}
        inline unsigned long long xyz2cubx16S(unsigned long long xyz){return (xyz>>LongOne);}
        inline unsigned long long xyz2bi16S(unsigned long long, unsigned long long y, unsigned long long z){return(4*(y & 3)+16*(z & 3));}
        inline int GetSpecies16S(int x, int y, int z) {
            return (((Long15<<(xyz2bi16S(x, y, z))) & git[xyz2cubx16S(x)][xyz2cubyz16S(y)][xyz2cubyz16S(z)])>>(xyz2bi16S(x, y, z)));
        }

        inline int GetSpecies4S(int x, int y, int z) {return ((LongThree<<(xyz2bi4S(x, y, z))) & git[xyz2cubxyz4S(x)][xyz2cubxyz4S(y)][xyz2cubxyz4S(z)])>>(xyz2bi4S(x, y, z));}
        inline int GetSpecies2S(int x, int y, int z) {return ((LongOne<<(xyz2bi2S(x, y, z))) & git[xyz2cubxyz2S(x)][xyz2cubxyz2S(y)][xyz2cubxyz2S(z)])>>(xyz2bi2S(x, y, z));}
};
