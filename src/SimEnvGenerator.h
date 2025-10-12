#pragma once

#include "Common.h"
#include "Utils.h"

namespace _2psp {
    class Model;
    class SimEnvGenerator {
    public:
        static Model* createScene(int SceneId, std::string solver = "TGS", int substeps = 150, dtype h = 1.0/60);
    };

}