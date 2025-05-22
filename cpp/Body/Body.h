#pragma once
#include "Common.h"
#include "Utils.h"

namespace _2psp {

    class Joint;
    class Model;
    class Body;

    class Body {
    public:
        // simulation
        Model* _sim;

        // name
        std::string _name = "";

        // structure
        Joint* _joint;
        //Body* _parent;

        std::vector<Body*> _contact_bodies_next;      // Contating bodies in next layer
        std::vector<Body*> _contact_bodies_current;    // Contating bodies in current layer

        // index
       int _index;                          // index for body in body list
       int _layer;                          // layer number for body in contact graph

        // data
        dtype _density;
        dtype _mass_inv;
        Vector6 _Inertia_inv;                   // diagonal of mass matrix M_i in body frame

        // constants
        SE3 _E_ij, _E_ji;                   // transformation between body i and its joint j, constant

        // variables
        SE3 _E_0i, _E_i0;                   // transformation between body i and world frame

        // body states
        MapVectorX _q;
        MapVectorX _dq;

        // phi
        se3 _phi;                           // spatial velocity
        se3 _phi_dt;                       // temporary variable for gauss-seidel solver

        dtype _mu;

        // rendering
        Vector3f _color;
        std::string _texture_path;
        bool _use_texture;

        Body(Model* sim, Joint* joint, dtype density);

        ~Body();

        // init body
        void init();

        // set rendering color
        void set_color(Vector3 color);
        void set_texture(std::string texture_path);

        // update
        void update();
        // void set_q(const VectorX& q);
        // void set_dq(const VectorX& dq);
        void interagate_state();

        Vector3 transform_point(const Vector3& p);

        // unconstrained step
        void step_unconstrained(se3& fm);
        void update_substep_states(dtype dt);

        void virtual compute_mass_inertial() = 0;
    };

}