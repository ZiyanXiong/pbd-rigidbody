#include "SimEnvGenerator.h"
#include "Model.h"
#include "Robot.h"
#include "Body/BodyCuboid.h"
#include "Body/BodyPlane.h"

#include "Utils.h"
#include "Common.h"

namespace _2psp{
	Model* SimEnvGenerator::createGroundTest(std::string solver) {
        // simulation options
        Model::Options* options = new Model::Options();
        options->_gravity = -980. * Vector3::UnitZ();
        options->_h = 0.01;
        options->_solver = solver;

        // construct simulation
        Model* sim = new Model(options, "Scene 0");

        // define robot
        Robot* robot = new Robot();
        Vector3 length(4.0, 4.0, 4.0);
        BodyCuboid* body1 = new BodyCuboid(sim, nullptr, length, (dtype)1.0);
        BodyCuboid* body2 = new BodyCuboid(sim, nullptr, length, (dtype)1.0);

        // add ground contact
        //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
        Matrix4 E_g = Matrix4::Identity();
        //E_g.topLeftCorner(3, 3) = R;
        VectorX q(14);
        q << 0, 0, 0, 1, 0, 0, 2,
            0, 0, 0, 1, 0, 0, 6;
        VectorX dq(14);
        dq << 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0;

        robot->add_body(body1);
        robot->add_body(body2);
        sim->add_robot(robot);
        sim->set_ground_plane(new BodyPlane(sim, E_g));
        sim->set_state_init(q, dq);
        sim->init();

        return sim;
	}
}