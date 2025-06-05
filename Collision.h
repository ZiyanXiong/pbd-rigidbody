#pragma once
#pragma once
#include "Common.h"
#include "Utils.h"
#include "CollisionDetection/Contact.h"

namespace _2psp {
    class Body;
    class Contact;

    class Collision {
    public:
        Body* _body1;
        Body* _body2;
        vector<Contact> _contacts;
        //JacobianMatrixVector _contact_frames;
        //JacobianMatrixVector _rxn_w1;
        //JacobianMatrixVector _rxn_w2;
        JacobianMatrixVector _J1;
        JacobianMatrixVector _J2;
        JacobianMatrixVector _J_div_m1;
        JacobianMatrixVector _J_div_m2;
        Matrix3X _w1;
        Matrix3X _w2;
        Matrix3X _d;
        Matrix3X _c;
        Matrix3X _c0;
        Matrix3X _lambdas_prev;
        Matrix3X _lambdas;
        VectorX _dlambdas_nor;

        dtype _mu;
        bool _shock_porpagate;
        bool _is_converged;
        bool _is_sp_valid;

        Collision(Body* body1, Body* body2, vector<Contact> contacts);

        void init();
        void compute_c(dtype h);
        void solve_nor(dtype h);
        void solve_tan(dtype h);
        void solve_nor_2psp(dtype h);
        void solve_tan_2psp(dtype h);
        void solve_vel_nor(dtype h);
        void solve_vel_tan(dtype h);

        void apply_accumulated_impulse();
        void apply_mass_averaged_impulse_nor();

        bool is_converged(dtype h, dtype tol);

        //void reset_time_report();
        //void print_time_report();
    };

}