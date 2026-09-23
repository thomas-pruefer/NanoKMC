#include "CalcClass.h"
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>


CalcClass::CalcClass() {

}

CalcClass::~CalcClass() {
    if (CalcData != nullptr) {
        for (int i = 0; i < RecordNumber; ++i) {
            delete[] CalcData[i];
        }
        delete[] CalcData;
    }
    delete[] Recorded;
}

double CalcClass::Str2Double (const std::string &str) {
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

int CalcClass::Str2Int (const std::string &str) {
    std::stringstream ss(str);
    int n;
    ss >> n;
    return n;
}

void CalcClass::Init() {
    // Ensure the standard run folders exist and clear existing evaluations.
    #ifdef _WIN32
        system("if not exist evaluation mkdir evaluation");
        system("if not exist output mkdir output");
        system("if not exist output\\bit mkdir output\\bit");
        system("del /Q evaluation\\*.* 2>nul");
    #else
        system("mkdir -p evaluation output/bit");
        system("rm -f evaluation/*");
    #endif

    // Extract existing lattice checkpoints when this run is resumed from
    // output/bit/data.zip. The archive stores lattice state only; RNG state is
    // intentionally not part of the checkpoint format.
    std::ifstream datazip("output/bit/data.zip");
    if (datazip.good()) {
        system("unzip -o -q output/bit/data.zip");
    }
    datazip.close();

    NSpecies = Str2Int(InputParam("NSpecies", "nanokmc.in"));
    if (NSpecies < 2 || NSpecies > 16) {
        std::cerr << "NSpecies must be between 2 and 16 for the packed lattice."
                  << std::endl;
        std::exit(2);
    }
    ColumnParam = 5 + NSpecies;

    // Blank lines are ignored deliberately. Each nonblank row is one checkpoint
    // record; accepting an empty row as an all-zero record can silently create an
    // unintended extra checkpoint.
    std::ifstream input("output/CalcData.csv");
    if (!input) {
        std::cerr << "Could not open output/CalcData.csv" << std::endl;
        std::exit(2);
    }
    std::vector<std::string> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
        rows.push_back(line);
    }
    input.close();

    RecordNumber = static_cast<int>(rows.size());
    if (RecordNumber == 0) {
        std::cerr << "output/CalcData.csv contains no checkpoint records" << std::endl;
        std::exit(2);
    }

    CalcData = new long long *[RecordNumber];
    for (int j = 0; j < RecordNumber; ++j) {
        CalcData[j] = new long long[ColumnParam]();
    }

    Recorded = new long long[RecordNumber]();
    RecordCount = RecordNumber;

    for (int j = 0; j < RecordNumber; ++j) {
        std::stringstream row(rows[j]);
        std::string field;
        int column = 0;
        while (std::getline(row, field, ';') && column < ColumnParam) {
            const std::size_t first = field.find_first_not_of(" \t\r\n");
            const std::size_t last = field.find_last_not_of(" \t\r\n");
            if (first == std::string::npos) {
                field = "0";
            } else {
                field = field.substr(first, last - first + 1);
            }
            try {
                std::size_t parsed = 0;
                const long long value = std::stoll(field, &parsed, 10);
                if (parsed != field.size()) throw std::invalid_argument("trailing characters");
                CalcData[j][column++] = value;
            } catch (const std::exception&) {
                std::cerr << "Invalid integer in output/CalcData.csv row " << (j + 1)
                          << ", column " << (column + 1) << ": " << field
                          << std::endl;
                std::exit(2);
            }
        }
        if (column < ColumnParam) {
            std::cerr << "output/CalcData.csv row " << (j + 1) << " has "
                      << column << " fields; expected at least " << ColumnParam
                      << std::endl;
            std::exit(2);
        }
        if (CalcData[j][0] < 0) {
            std::cerr << "output/CalcData.csv row " << (j + 1)
                      << " has a negative requested MCS checkpoint" << std::endl;
            std::exit(2);
        }
        if (j > 0 && CalcData[j][0] < CalcData[j - 1][0]) {
            std::cerr << "output/CalcData.csv checkpoints must be nondecreasing; "
                      << "row " << (j + 1) << " is smaller than the previous row"
                      << std::endl;
            std::exit(2);
        }

        if (CalcData[j][2] == 0 && RecordNumber == RecordCount) {
            RecordCount = j;
        }
        if (CalcData[j][1] != 0) {
            Recorded[j] = 1;
        }
    }

    SystemID = InputParam("SystemID", "nanokmc.in");
}

std::string CalcClass::InputParam(std::string paramname, std::string filename) {
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

void CalcClass::AddCalcData(long long* SystemData, int i) {
    for (int j=0; j<ColumnParam; j++) {
        CalcData[i][j]=SystemData[j];
    }
}

void CalcClass::WriteCalcData() {
    int j, i;
    std::string str, str2, index;
    CalcDatafile.open("output/CalcData.csv", std::fstream::out);
    for (i=0; i<RecordNumber; i++) {
        for (j=0; j<ColumnParam; j++) {
            CalcDatafile << CalcData[i][j] << ";";
        }
        CalcDatafile << std::endl;
    }
    CalcDatafile.close();

    // write index parameter settings
    std::fstream evalfile("evaluation/eval.sce", std::ios::out);
    index=InputParam("Index", "nanokmc.in");
    while (index.find("_") != std::string::npos) {
        j=index.find("_");
        str=index.substr(0,j);
        index=index.substr(j+1,index.size());
        if (str.substr(0,3)=="Int") {
            str=str.substr(3,str.size());
            if (str.substr(0,3)=="Sys") {
                str2=str.substr(3,str.size());
            } else {
                str2=str.substr(4,str.size());
            }
            str2=InputParam(str2, "nanokmc.in");
            evalfile << "setIntMeta(\"" << str << "\"," << str2 << ");" << std::endl;
        }
    }
    if (index.substr(0,3)=="Int") {
        index=index.substr(3,index.size());
        if (index.substr(0,3)=="Sys") {
            str2=index.substr(3,index.size());
        } else {
            str2=index.substr(4,index.size());
        }
        str2=InputParam(str2, "nanokmc.in");
        evalfile << "setIntMeta(\"" << index << "\"," << str2 << ");" << std::endl;
    }
    str2=InputParam("Fluence", "nanokmc.in");
    evalfile << "setIntMeta(\"" << "SysFluence" << "\"," << str2 << ");" << std::endl;
    str2=InputParam("Seed", "nanokmc.in");
    evalfile << "setIntMeta(\"" << "CalcSeed" << "\"," << str2 << ");" << std::endl;
    evalfile.close();
}

long long* CalcClass::GetStepData(int i) {

    return CalcData[i];

}
