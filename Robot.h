#pragma once
#include "Common.h"
#include "Utils.h"
#include "Collision.h"

namespace _2psp {

    //class Joint;
    class Body;
    //class Force;
    //class Actuator;
    class Collision;
    //class Constraint;

    class Robot {
    public:
        //vector<Joint*> _joints;
        vector<Body*> _bodies;
        //vector<Force*> _forces;
       // vector<Actuator*> _actuators;
       // vector<Constraint*> _constraints;
        vector<Collision> _collisions;
        vector<vector<size_t>> _collision_layer; 

        int _ndof_m, _ndof_u;
        int _ind_m, _ind_u;
        int _n_c; // number of constraints


        Robot(): _ndof_m(-1), _ndof_u(-1), _ind_m(-1), _ind_u(-1) {
            //_joints.clear();
            _bodies.clear();
            //_forces.clear();
            //_actuators.clear();
            _collisions.clear();
            _collision_layer.clear();
        }

        void add_body(Body* body);

       // void add_actuator(Actuator* actuator);

        //void add_constraint(Constraint* constraint);

        //void contact_culling(int k);

        // init robot
        void init(int ind_m);
        void init_collisions();
        void reset();
        //void construct_dfs_order(Joint* now);

        void construct_collision_order();
        void solve_collisions(dtype h);
        bool solve_collisions_2psp(dtype h, int& solve_count, int sp_iter_max, dtype tol);
        void solve_velocity(dtype h);

        void step_unconstrained(VectorX& f_t);

        // set state variables
        // void set_q(const VectorX q);
        // void set_dq(const VectorX qdot);
        void update();
        void set_state_memory(dtype* const qp_start, dtype* const dqp_start);
        void interagate_state();

        // get state variables
        //VectorX get_phi();

        //void reset_time_report();
        //void print_time_report();
    };

}