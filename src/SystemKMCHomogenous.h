#pragma once

#include "SystemClass.h"

#include <string.h>

class SystemKMCHomogenous : public SystemClass {

    public:
        SystemKMCHomogenous();
        ~SystemKMCHomogenous();

    protected:
        void InitSubSys() override;
        void EvalSysStep(unsigned long long n) override;
        double DistributionFunction(double x, double y, double z, int l) override;
        double MCSIncrementPerAttempt() const override;
        double clvl;

};
