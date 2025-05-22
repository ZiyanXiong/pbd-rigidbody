#pragma once

#include "Common.h"
#include "Utils.h"

namespace _2psp {
    class Model;
    class SimEnvGenerator {
    public:
        static Model* createGroundTest(std::string solver = "TGS");
    };

}