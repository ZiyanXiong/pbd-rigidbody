#include "Body/Body.h"
#include "Utils.h"
#include "Model.h"

namespace _2psp {


    Body::Body(Model* sim, Joint* joint, dtype density):
        _mass_inv(0.0),
        _index(-1),
        _layer(-1),
        _mu(0.5),
        _q(nullptr, 0),
        _dq(nullptr, 0)
    {
        _sim = sim;

        _joint = joint;

        _density = density;

        _color << 0.25f, 0.148f, 0.06475f;

        _use_texture = false;

        _contact_bodies_next.clear();

    }


    Body::~Body() {

    }


    void Body::init() {
        compute_mass_inertial();
        _E_0i.setIdentity(4,4);
        _E_i0.setIdentity(4,4);
        _phi.setZero();
        _phi_0.setZero();
        _phi_dt.setZero();
        _delta_phi.setZero();
    }

    void Body::set_color(Vector3 color) {
        _color = color.cast<float>();
        _use_texture = false;
    }

    void Body::set_texture(std::string texture_path) {
        _texture_path = texture_path;
        _use_texture = true;
    }
    /*
    void Body::set_q(const VectorX& q) {
        _E_0i.topLeftCorner(3, 3) = Quat(q.segment<4>(0)).toRotationMatrix();
        _E_0i.topRightCorner(3, 1) = q.segment<3>(4);
        _E_i0.topLeftCorner(3, 3) = _E_0i.topLeftCorner(3, 3).transpose();
        _E_i0.topRightCorner(3, 1) = -q.segment<3>(4);
    }

    void Body::set_dq(const VectorX& dq) {
        AngleAxis dAa(Quat(dq.segment<4>(0)));
        _phi.segment<3>(0) = dAa.angle() * dAa.axis();
        _phi.segment<3>(3) = dq.segment<3>(4);
    }
    */
    void Body::update() {
        _E_i0.topLeftCorner(3, 3) = Quat(_q.segment<4>(0)).toRotationMatrix();
        _E_i0.topRightCorner(3, 1) = _q.segment<3>(4);
        _E_0i.topLeftCorner(3, 3) = _E_i0.topLeftCorner(3, 3).transpose();
        _E_0i.topRightCorner(3, 1) = -_q.segment<3>(4);
        AngleAxis dAa(Quat(_dq.segment<4>(0)));
        _phi.segment<3>(0) = dAa.angle() * dAa.axis();
        _phi.segment<3>(3) = _dq.segment<3>(4);
    }

    void Body::step_unconstrained(MapVectorX& fm) {
        Matrix3 I_inv = _E_i0.topLeftCorner(3, 3) * _Inertia_inv.head<3>().asDiagonal() * _E_i0.topLeftCorner(3, 3).transpose();
        Matrix3 I;
        if (_Inertia_inv.head<3>().norm() > math::eps_big) {
			I = _E_i0.topLeftCorner(3, 3) * _Inertia_inv.head<3>().cwiseInverse().asDiagonal() * _E_i0.topLeftCorner(3, 3).transpose();
		}
        else {
			I.setIdentity(3, 3);
		}
        dtype mass;
        if (_mass_inv > math::eps_big) {
			mass = 1.0 / _mass_inv;
		}
        else {
			mass = 1.0;
		}
        Vector3 t = fm.tail<3>() - _phi.segment<3>(0).cross(I * _phi.segment<3>(0));
        Vector3 f = fm.head<3>() + _sim->_options->_gravity * mass;

        _phi.segment<3>(0) += _sim->_options->_h * I_inv * t;
        _phi.segment<3>(3) += _sim->_options->_h * _mass_inv * f;
        _phi_0 = _phi; // save the unconstrained step result
    }

    Vector3 Body::transform_point(const Vector3& p) {
        // Transform point p from body frame to world frame
        return _E_i0.topLeftCorner(3,3)* p + _E_i0.topRightCorner(3,1);
    }

    void Body::update_substep_states(dtype dt) {
        _phi_dt += _phi * dt;
        _delta_phi += _phi * dt; // accumulate the delta phi for velocity solve
    }

    void Body::interagate_state() {
        //std::cout << "Body States Phi:\n" << _phi.transpose() << std::endl;
        //std::cout << "Body States Phi_dt:\n" << _phi_dt.transpose() << std::endl;

        Vector3 dtheta_axis = _delta_phi.segment<3>(0).normalized();
        AngleAxis dAa(_delta_phi.segment<3>(0).norm(), dtheta_axis);
        //std::cout << "dAa:\n" << Quat(dAa) << std::endl;
        //std::cout << "q0:\n" << _q.segment<4>(0) << std::endl;
        _q.segment<4>(0) = (dAa * Quat(_q.segment<4>(0))).normalized().coeffs();
        _q.segment<3>(4) += _delta_phi.segment<3>(3);
        //std::cout << "q1:\n" << _q.segment<4>(0) << std::endl;
        //std::cout << "q1_new:\n" << (Quat(dAa) * Quat(_q.segment<4>(0))).normalized().coeffs() << std::endl;

        _dq.segment<4>(0) = Quat(AngleAxis(_phi.segment<3>(0).norm(), _phi.segment<3>(0).normalized())).coeffs();
        _dq.segment<3>(4) = _phi.segment<3>(3);
        _phi.setZero();
        _phi_dt.setZero();
        _delta_phi.setZero();
        _layer = -1;
        _contact_bodies_next.clear();
    }

}