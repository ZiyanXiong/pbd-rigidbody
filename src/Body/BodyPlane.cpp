#include "Body/BodyPlane.h"
#include "Utils.h"

namespace _2psp {


    BodyPlane::BodyPlane(Model* sim, const Matrix4 E_g) :
        BodyPrimitiveShape(sim, nullptr, 0)
    {
        _E_i0 = E_g;
        _E_0i = math::Einv(_E_i0);
    }

    void BodyPlane::init() 
    {
        compute_mass_inertial();
        _phi.setZero();
        _phi_0.setZero();
        _phi_dt.setZero();
        _delta_phi.setZero();
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
