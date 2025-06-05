#pragma once
#include "Body/BodyPrimitiveShape.h"

namespace _2psp {

    class BodyCuboid : public BodyPrimitiveShape {
    public:
        Vector3 _length;
        bool _is_infinite_mass;

        BodyCuboid(Model* sim, Joint* joint, Vector3 length, dtype density, bool is_infinite_mass = false);

        void compute_mass_inertial();
    };

}