#include "Body/BodyPlane.h"
#include "Utils.h"

namespace _2psp {


    BodyPlane::BodyPlane(Model* sim, const Matrix4 E_g) :
        BodyPrimitiveShape(sim, nullptr, 0)
    {
        _E_0i = E_g;
        _E_i0 = math::Einv(_E_0i);
    }

    void BodyPlane::compute_mass_inertial()
    {
        // Compute the mass
        _mass_inv = 0;

        // Compute the inertia
        _Inertia_inv.segment<3>(0) = Vector3::Zero();
        _Inertia_inv.segment<3>(3) = Vector3::Ones() * _mass_inv;
    }

}
