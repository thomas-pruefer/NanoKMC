#pragma once

#include <string>
#include <fstream>

class CalcClass {

    public:
        CalcClass();
        ~CalcClass();

        void Init();
        long long* GetStepData(int i);
        int Str2Int (const std::string &str);
        double Str2Double (const std::string &str);
        void AddCalcData(long long* SystemData, int i);
        std::string InputParam(std::string paramname, std::string filename);
        void WriteCalcData();

        long long* Recorded = nullptr;
        std::string SystemID;
        long long **CalcData = nullptr;
        int RecordCount = 0, RecordNumber = 0, ColumnParam = 0;

    protected:
        int NSpecies;
        std::fstream CalcDatafile;
        const double lc=0.4338;
};
