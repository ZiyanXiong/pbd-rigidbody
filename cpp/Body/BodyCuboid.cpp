#include "Body/BodyCuboid.h"
#include "Utils.h"

namespace _2psp {


    BodyCuboid::BodyCuboid(Model* sim, Joint* joint, Vector3 length, dtype density) :
        BodyPrimitiveShape(sim, joint, density)
    {
        _length = length;
    }

    void BodyCuboid::compute_mass_inertial()
	{
		// Compute the mass
		dtype mass = _density * _length.x() * _length.y() * _length.z();
        _mass_inv = 1.0 / mass;

        // Compute the inertia
        Vector6 Inertia;
		Inertia(0) = (mass / 12.0) * (_length.y() * _length.y() + _length.z() * _length.z());
		Inertia(1) = (mass / 12.0) * (_length.x() * _length.x() + _length.z() * _length.z());
		Inertia(2) = (mass / 12.0) * (_length.x() * _length.x() + _length.y() * _length.y());
        Inertia.segment<3>(3) = Vector3::Ones() * mass;
        _Inertia_inv = Inertia.cwiseInverse();
	}

}
