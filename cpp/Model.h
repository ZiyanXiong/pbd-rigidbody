#pragma once
#include "Common.h"
#include "Utils.h"

namespace _2psp {

    class Robot;
    class BodyPlane;

    class Model {
    private:
        // state history
        std::vector<VectorX> _q_his, _dq_his;

        // compute M and f matrices
        //void computeMatrices(MatrixX& M, VectorX& f); // evaluate force and mass matrix
        //void computeMatrices(MatrixX& M, VectorX& f, JacobianMatrixVector& dM_dq, MatrixX& K, MatrixX& D); // evaluate derivatives to q
        //void computeMatrices(MatrixX& M, VectorX& f, JacobianMatrixVector& dM_dq, MatrixX& K, MatrixX& D, MatrixX& df_du); // evaluate derivatives to q and to optimization parameters
        //void computeMatrices(
        //    MatrixX& M, VectorX& f,
        //    JacobianMatrixVector& dM_dq, MatrixX& K, MatrixX& D,
        //    MatrixX& df_du,
        //    JacobianMatrixVector& dM_dp, MatrixX& df_dp); // evaluate derivatives to q and to optimization parameters

        //// compute variables
        //void computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq);
        //void computeVariablesWithDerivative(VectorX& variables, MatrixX& dvar_dq, MatrixX& dvar_dp);

        // temporary variables for evaluation functions
        // VectorX _q1, _qdot1, _q0, _qdot0, _q_alpha, _qdot_alpha;

        // xml file related
        std::string _asset_folder;

    public:
        // -------------------- Constants -----------------------
        class Options {
        public:
            Vector3 _gravity;
            dtype _h;
            int _substep; // substep for temporal gauss seidel
            string _solver; // [2PSP, TGS]
            string _unit; // ["cm-g", "m-kg"]

            Options(Vector3 gravity = -980. * Vector3::UnitZ(), dtype h = 0.02, int substep = 100, string solver = "TGS", string unit = "cm-g") :
                _gravity(gravity), _h(h), _substep(substep), _solver(solver), _unit(unit) {}
        };

        Options* _options;


        class ViewerOptions {
        public:
            int _fps;
            dtype _speed;
            Vector3 _camera_pos;
            Vector3 _camera_up;
            Vector3 _camera_lookat;
            bool _ground;
            Matrix4 _E_g;
            bool _record;
            std::string _record_folder;
            bool _loop;     // whether loop the replay
            bool _infinite; // whether to replay until close the window

            ViewerOptions() {
                _fps = 30;
                _speed = 1.;
                _camera_pos = Vector3(4., -5., 3.);
                _camera_up = Vector3::UnitZ();
                _camera_lookat = Vector3::Zero();
                _ground = true;
                _E_g.topLeftCorner(3, 3).setIdentity();
                _E_g.topRightCorner(3, 1) = Vector3::UnitZ() * -2.;
                _record = false;
                _loop = true;
                _infinite = true;
            }
        };

        ViewerOptions* _viewer_options;

        class TimeReport {
        public:
            long long _time_solver, _time_save_backward, _time_backward;
            long long _time_compute_matrices, _time_compose_matrices;
            long long _time_compute_dJ, _time_compute_df;
            long long _time_dM_dp, _time_df_dp;
            long long _time_dM_dp1, _time_dM_dp2, _time_dM_dp4;

            void reset() {
                _time_solver = _time_save_backward = _time_backward = 0;
                _time_compute_matrices = _time_compose_matrices = 0;
                _time_compute_dJ = _time_compute_df = 0;
                _time_dM_dp = _time_df_dp = 0;
                _time_dM_dp1 = _time_dM_dp2 = _time_dM_dp4 = 0;
            }
        };

        TimeReport _time_report;

        // -------------------- forward dynamics related -------------------
        std::string _name;

        // robot related
        vector<Robot*> _robots;

        bool _ground;
        BodyPlane* _ground_plane;

        int  _ndof_m, _ndof_u;

        // controller parameters
        VectorX _phi;

        // states
        VectorX _q_init, _dq_init;
        VectorX _q, _dq;

        // verbose output
        bool _verbose;

        // constructors
        Model(Options* options, std::string name = "");
        //Simulation(Options* options, ViewerOptions* viewer_options, std::string name = "");
        //Simulation(std::string xml_file_path, bool verbose = false);

        // destructor
        ~Model();

        //Joint* parse_from_xml_file(pugi::xml_node root, pugi::xml_node node, \
        //    Joint* parent_joint, int& joint_cnt, bool verbose = false);

        void add_robot(Robot* robot) {
            _robots.push_back(robot);
        }

        // init simulation
        void init();

        // init states
        void set_state_init(const VectorX q_init, const VectorX dq_init);
        void set_q_init(const VectorX q_init);
        void set_dq_init(const VectorX dq_init);
        //const VectorX get_q_init();
        //const VectorX get_dq_init();

        // states
        void set_state(const VectorX q, const VectorX dq);
        void set_q(const VectorX q);
        void set_dq(const VectorX dq);
        //const VectorX get_q();
        //const VectorX get_qdot();

        // control variables
        //void set_u(const VectorX& u);
        //void get_ctrl_range(VectorX& ctrl_min, VectorX& ctrl_max);
        //void print_ctrl_info();

        void set_ground_plane(BodyPlane* ground_plane);

        // updata robot to propagate the state
        void update_robot_states_memory();

        // reset the simulation
        //void reset();

        // verbose defines if log state history
        void forward(int num_steps, bool save_history = false);

        void collision_detection();
        void step_unconstrained();
        void temporal_gauss_seidel();
        void update_robot();

        // export simulation replay to a folder
        void export_replay(std::string folder);

        //void print_time_report();
    };

}