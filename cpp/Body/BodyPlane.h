#pragma once
#include "Body/BodyPrimitiveShape.h"

namespace _2psp {

    class BodyPlane : public BodyPrimitiveShape {
    public:
        Vector3 _length;

        BodyPlane(Model* sim, const Matrix4 E_g);

        void compute_mass_inertial();
    };

}