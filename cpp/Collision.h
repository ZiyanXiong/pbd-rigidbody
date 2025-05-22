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
        Matrix3X _lambdas;
        dtype _mu;

        Collision(Body* body1, Body* body2, vector<Contact> contacts);

        void init();
        void solve_nor(dtype h);
        void solve_tan(dtype h);

        //void reset_time_report();
        //void print_time_report();
    };

}