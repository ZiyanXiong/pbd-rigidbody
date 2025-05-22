#pragma once
#include "Body/BodyPrimitiveShape.h"

namespace _2psp {

    class BodyCuboid : public BodyPrimitiveShape {
    public:
        Vector3 _length;

        BodyCuboid(Model* sim, Joint* joint, Vector3 length, dtype density);

        void compute_mass_inertial();
    };

}