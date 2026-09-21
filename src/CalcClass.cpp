#include "CalcClass.h"
#include <sstream>
#include <fstream>
#include <cmath>


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
    if (str.find("D")!=-1) {
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
    int ColumnParam0, l1, l2;
    std::fstream CalcInputfile;
    char s[255];
    std::string str,inputvalue[15];
    std::string inputs[15];
    int i, j, k, l;
//std::cout << " sd" << std::endl;
//exit(0);
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

    // Extract existing bit files when this run is resumed from output/bit/data.zip.
    std::ifstream datazip("output/bit/data.zip");
    if (datazip.good()) {
        system("unzip -o -q output/bit/data.zip");
    }
    datazip.close();

    // Count rows and colums
    NSpecies=Str2Int(InputParam("NSpecies", "nanokmc.in"));
    ColumnParam=5+NSpecies;
    CalcDatafile.open("output/CalcData.csv", std::fstream::in);

    RecordNumber=0;
    while (CalcDatafile.getline(s,255)) {
        str=s;
        RecordNumber++;
    }
    CalcDatafile.close();

    // initialize CalcData array
    CalcData = new long long *[RecordNumber];
    for (j=0; j<RecordNumber; j++) {
        CalcData[j] = new long long[ColumnParam];
        for (i=0; i<ColumnParam; i++) {
            CalcData[j][i]=0;
        }
    }

    // Initialize Recorded
    Recorded = new long long [RecordNumber];
    for (i=0; i<RecordNumber; i++) {
        Recorded[i]=0;
    }
    //Recorded[0]=1;

    // Read CalcData.dat
    CalcDatafile.open("output/CalcData.csv", std::fstream::in);
    RecordCount=RecordNumber;

    for (j=0; j<RecordNumber; j++) {
        CalcDatafile.getline(s,255);
        str=s;
        l1=0;
        l2=str.find(";",l1);
        k=0;
        while (l2>-1) {
            CalcData[j][k]=Str2Int(str.substr(l1,l2-l1));
            k=k+1;
            l1=l2+1;
            l2=str.find(";",l1);
        }
        //std::cout << RecordCount << std::endl;
        if (CalcData[j][2]==0 && RecordNumber==RecordCount) {
            RecordCount=j;
        }
        if (CalcData[j][1]!=0) {
            Recorded[j]=1;
        }
    }
    CalcDatafile.close();

    SystemID=InputParam("SystemID", "nanokmc.in");
}

std::string CalcClass::InputParam(std::string paramname, std::string filename) {
    std::string line;
    std::string param = "-1";
    std::string key = paramname + "=";

    std::fstream inputfile(filename.c_str(), std::ios::in);
    while (std::getline(inputfile, line)) {
        std::size_t j = line.find(key);
        if (j == std::string::npos) continue;

        std::size_t eq = line.find("=", j);
        std::size_t semi = line.find(";", eq);
        if (eq == std::string::npos) continue;
        if (semi == std::string::npos) semi = line.size();

        param = line.substr(eq + 1, semi - eq - 1);

        // trim whitespace
        std::size_t first = param.find_first_not_of(" \t\r\n");
        std::size_t last = param.find_last_not_of(" \t\r\n");
        if (first == std::string::npos) {
            param = "";
        } else {
            param = param.substr(first, last - first + 1);
        }

        // strip optional quotes
        if (param.size() >= 2 && param.front() == '"' && param.back() == '"') {
            param = param.substr(1, param.size() - 2);
        }
        break;
    }
    inputfile.close();
    return param;
}

void CalcClass::AddCalcData(long long* SystemData, int i) {
    int j, k;

    for (j=0; j<ColumnParam; j++) {
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
    while (index.find("_")!=-1) {
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
