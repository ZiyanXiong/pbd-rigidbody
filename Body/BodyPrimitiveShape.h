#pragma once
#include "Body/Body.h"

namespace _2psp {
    class Shape;

    class BodyPrimitiveShape : public Body {
    public:
        BodyPrimitiveShape(Model* sim, Joint* joint, dtype density)
            :Body(sim, joint, density) { }
    };

}