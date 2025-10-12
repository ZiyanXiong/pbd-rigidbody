#include "SimEnvGenerator.h"
#include "Model.h"
#include "Robot.h"
#include "Body/BodyCuboid.h"
#include "Body/BodyPlane.h"
#include "Body/BodyMesh.h"

#include "Utils.h"
#include "Common.h"

namespace _2psp{
	Model* SimEnvGenerator::createScene(int SceneId, std::string solver, int substeps, dtype h) {
        // simulation options
        Model::Options* options = new Model::Options();
        options->_gravity = -980. * Vector3::UnitZ();
        options->_h = h;
        options->_substep = substeps; // substep for temporal gauss seidel
        options->_solver = solver;

        switch (SceneId) {
        case 0:
        {
            // Stacking 10 boxes
            // construct simulation
            Model* sim = new Model(options, "Scene 0");

            // define robot
            Robot* robot = new Robot();
            Vector3 length(4.0, 4.0, 4.0);

            int n_bodies = 2; // number of bodies
            VectorX q(7 * n_bodies);
            VectorX dq(7 * n_bodies);

            for (int i = 0; i < n_bodies; ++i) {
                robot->add_body(new BodyCuboid(sim, nullptr, length, (dtype)1.0));
                q.segment<7>(7 * i) << 0., 0., 0., 1., 0., 0., 4. * i + 2.; // orientation and position
                dq.segment<7>(7 * i) << 0., 0., 0., 1., 0., 0., 0.; // velocity
            }


            // add ground contact
            //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
            Matrix4 E_g = Matrix4::Identity();

            sim->add_robot(robot);
            sim->set_ground_plane(new BodyPlane(sim, E_g));
            sim->set_state_init(q, dq);
            sim->init();
            return sim;
            break;
        }
        case 4: 
        {
            // Stacking 10 boxes on a slope
            // construct simulation
            options->_2psp_iter_max = 35;
            Model* sim = new Model(options, "Scene 4");

            // define robot
            Robot* robot = new Robot();
            Vector3 length(4.0, 4.0, 4.0);
            dtype angle = constants::pi / 12.;
            Quat rotation(AngleAxis(angle, Vector3::UnitY()));

            int n_bodies = 10; // number of bodies
            VectorX q(7 * n_bodies);
            VectorX dq(7 * n_bodies);

            for (int i = 0; i < n_bodies; ++i) {
                robot->add_body(new BodyCuboid(sim, nullptr, length, (dtype)1.0));
                dtype x = -i;
                dtype y = 0;
                dtype z = 4. * i + 2.;
                Vector3 pos = rotation * Vector3(x, y, z);
                q.segment<7>(7 * i) << rotation.coeffs(), pos; // orientation and position
                dq.segment<7>(7 * i) << 0., 0., 0., 1., 0., 0., 0.; // velocity
            }

            // add ground contact
            //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
            Matrix4 E_g = Matrix4::Identity();
            E_g.topLeftCorner(3, 3) = rotation.toRotationMatrix();

            sim->add_robot(robot);
            sim->set_ground_plane(new BodyPlane(sim, E_g));
            sim->set_state_init(q, dq);
            sim->init();
            return sim;
            break;
        }
        case 5:
        {
            // Stacking: Inverted Tower
            Model* sim = new Model(options, "Scene 5");

            // define robot
            Robot* robot = new Robot();

            int n_bodies = 5; // number of bodies
            VectorX q(7 * n_bodies);
            VectorX dq(7 * n_bodies);

            for (int i = 0; i < n_bodies; ++i) {
                dtype scale = 4.0 * pow(2.0, i);
                robot->add_body(new BodyCuboid(sim, nullptr, Vector3(scale, scale, scale), (dtype)1.0));
                q.segment<7>(7 * i) << 0., 0., 0., 1., 0., 0., 1.5 * scale - 4; // orientation and position
                dq.segment<7>(7 * i) << 0., 0., 0., 1., 0., 0., 0.; // velocity
            }

            // add ground contact
            //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
            Matrix4 E_g = Matrix4::Identity();

            sim->add_robot(robot);
            sim->set_ground_plane(new BodyPlane(sim, E_g));
            sim->set_state_init(q, dq);
            sim->init();
            return sim;
            break;
        }
        case 8:
        {
            // Stacking: Unbalanced goal post
            // construct simulation
            options->_2psp_iter_max = 175;

            Model* sim = new Model(options, "Scene 8");

            // define robot
            Robot* robot = new Robot();
            Vector3 length(4.0, 4.0, 4.0);

            int n_bodies = 20; // number of bodies
            VectorX q(7 * n_bodies);
            VectorX dq(7 * n_bodies);

            robot->add_body(new BodyCuboid(sim, nullptr, length, (dtype)1.0));
            q.segment<7>(0) << 0., 0., 0., 1., 0., 0., 2.; // orientation and position
            dq.segment<7>(0) << 0., 0., 0., 1., 0., 0., 0.; // velocity

            for (int i = 1; i < 9; ++i) {
                robot->add_body(new BodyCuboid(sim, nullptr, length, (dtype)1.0));
                q.segment<7>(7 * i) << 0., 0., 0., 1., -20., 0., 4. * i + 6.; // orientation and position
                dq.segment<7>(7 * i) << 0., 0., 0., 1., 0., 0., 0.; // velocity
            }

            for (int i = 9; i < 18; ++i) {
                robot->add_body(new BodyCuboid(sim, nullptr, length, (dtype)1.0));
                q.segment<7>(7 * i) << 0., 0., 0., 1., 20., 0., 4. * (i - 8) + 6.; // orientation and position
                dq.segment<7>(7 * i) << 0., 0., 0., 1., 0., 0., 0.; // velocity
            }

            robot->add_body(new BodyCuboid(sim, nullptr, Vector3(11. * 4., 4., 4.), (dtype)1.0));
            q.segment<7>(7 * 18) << 0., 0., 0., 1., 0., 0., 6.; // orientation and position
            dq.segment<7>(7 * 18) << 0., 0., 0., 1., 0., 0., 0.; // velocity

            robot->add_body(new BodyCuboid(sim, nullptr, length, (dtype)2.7));
            q.segment<7>(7 * 19) << 0., 0., 0., 1., -20., 0., 46.; // orientation and position
            dq.segment<7>(7 * 19) << 0., 0., 0., 1., 0., 0., 0.; // velocity

            // add ground contact
            //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
            Matrix4 E_g = Matrix4::Identity();

            sim->add_robot(robot);
            sim->set_ground_plane(new BodyPlane(sim, E_g));
            sim->set_state_init(q, dq);
            sim->init();
            return sim;
            break;
        }
        case 9:
        {
            // Stacking: Earthquake
            // construct simulation
            options->_2psp_iter_max = 75;
            options->_2psp_tol = math::eps_big;
            Model* sim = new Model(options, "Scene 9");

            // define robot
            Robot* robot = new Robot();
            Vector3 length(4.0, 4.0, 4.0);

            int n_bodies = 16; // number of bodies
            VectorX q(7 * n_bodies);
            VectorX dq(7 * n_bodies);

            robot->add_body(new BodyCuboid(sim, nullptr, Vector3(40,40,0.2), (dtype)1.0, true));
            q.segment<7>(0) << 0., 0., 0., 1., 0., 0., 1.; // orientation and position
            dq.segment<7>(0) << 0., 0., 0., 1., 0., 0., 0.; // velocity

            for (int i = 1; i < n_bodies; ++i) {
                robot->add_body(new BodyCuboid(sim, nullptr, length, (dtype)1.0));
                q.segment<7>(7 * i) << 0., 0., 0., 1., 0., 0., 4. * (i - 1) + 2. + 1.1; // orientation and position
                dq.segment<7>(7 * i) << 0., 0., 0., 1., 0., 0., 0.; // velocity
            }


            // add ground contact
            //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
            Matrix4 E_g = Matrix4::Identity();

            sim->add_robot(robot);
            sim->set_ground_plane(new BodyPlane(sim, E_g));
            sim->set_state_init(q, dq);
            sim->init();
            return sim;
            break;
        }
        case 10:
        {
            // Stacing: Jenga Add
            // construct simulation
            options->_2psp_iter_max = 75;
            options->_2psp_tol = math::eps_big;
            Model* sim = new Model(options, "Scene 10");
            // define robot
            Robot* robot = new Robot();
            Vector3 length(8.0, 24.0, 4.0);
            dtype density = 0.6; // density of the cuboid bodies
            dtype mu = 0.2; // friction coefficient

            int layers = 25; // number of layers
            int n_bodies = 53; // number of bodies
            VectorX q(7 * n_bodies);
            VectorX dq(7 * n_bodies);
            Eigen::VectorXi layer_pattern(25);
            layer_pattern << 4, 10, 8, 6, 2, 2, 1, 9, 7, 8, 1, 10,
                9, 3, 2, 2, 4, 6, 5, 3, 7, 2, 3, 4, 4;
            int body_ind = 0;
            for (int l = 0; l < layers; ++l) {
                for (int i = 0; i < 3; ++i) {
                    if(layer_pattern(l) < 3 && (i == 0 || i == 2))
						continue; 
                    if (layer_pattern(l) == 3 && (i == 0))
                        continue;
                    if(layer_pattern(l) == 4 && (i == 1))
						continue;
                    if (layer_pattern(l) == 5 && (i == 2))
                        continue;

                    robot->add_body(new BodyCuboid(sim, nullptr, length, density));
                    robot->_bodies.back()->_mu = mu; // set friction coefficient
                    dtype angle = constants::pi / 2. * (l % 2);
                    dtype x = -4.03 * 4 + 2.05 * 4 * (i + 1);
                    dtype y = 0;
                    dtype z = 4. * l + 2.;
                    if (l == layers - 1 && i == 0)
                        z += 3.0 * 4;
                    Quat rotation(AngleAxis(angle, Vector3::UnitZ()));
                    Vector3 pos = rotation * Vector3(x, y, z);
                    q.segment<7>(7 * body_ind) << rotation.coeffs(), pos; // orientation and position
                    dq.segment<7>(7 * body_ind) << 0., 0., 0., 1., 0., 0., 0.; // velocity
                    body_ind++;
                }
            }

            // add ground contact
            //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
            Matrix4 E_g = Matrix4::Identity();

            sim->add_robot(robot);
            sim->set_ground_plane(new BodyPlane(sim, E_g));
            sim->set_state_init(q, dq);
            sim->init();
            return sim;
            break;
        }
        case 11:
        {
            // Stacing: Jenga Add
            // construct simulation
            options->_2psp_iter_max = 75;
            options->_2psp_tol = math::eps_big;
            Model* sim = new Model(options, "Scene 11");
            // define robot
            Robot* robot = new Robot();
            Vector3 length(8.0, 24.0, 4.0);
            dtype density = 0.6; // density of the cuboid bodies
            dtype mu = 0.2; // friction coefficient

            int layers = 25; // number of layers
            int n_bodies = 53; // number of bodies
            VectorX q(7 * n_bodies);
            VectorX dq(7 * n_bodies);
            Eigen::VectorXi layer_pattern(25);
            layer_pattern << 4, 10, 8, 6, 2, 2, 1, 9, 7, 8, 1, 10,
                9, 3, 2, 2, 4, 6, 5, 3, 7, 2, 3, 4, 5;
            int body_ind = 0;
            for (int l = 0; l < layers; ++l) {
                for (int i = 0; i < 3; ++i) {
                    if (layer_pattern(l) < 3 && (i == 0 || i == 2))
                        continue;
                    if (layer_pattern(l) == 3 && (i == 0))
                        continue;
                    if (layer_pattern(l) == 4 && (i == 1))
                        continue;
                    if (layer_pattern(l) == 5 && (i == 2))
                        continue;

                    robot->add_body(new BodyCuboid(sim, nullptr, length, density));
                    robot->_bodies.back()->_mu = mu; // set friction coefficient
                    dtype angle = constants::pi / 2. * (l % 2);
                    dtype x = -4.03 * 4 + 2.05 * 4 * (i + 1);
                    dtype y = 0;
                    dtype z = 4. * l + 2.;

                    Quat rotation(AngleAxis(angle, Vector3::UnitZ()));
                    Vector3 pos = rotation * Vector3(x, y, z);
                    q.segment<7>(7 * body_ind) << rotation.coeffs(), pos; // orientation and position
                    dq.segment<7>(7 * body_ind) << 0., 0., 0., 1., 0., 0., 0.; // velocity
                    body_ind++;
                }
            }

            // add ground contact
            //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
            Matrix4 E_g = Matrix4::Identity();

            sim->add_robot(robot);
            sim->set_ground_plane(new BodyPlane(sim, E_g));
            sim->set_state_init(q, dq);
            sim->init();
            return sim;
            break;
        }
        case 12:{
            // Stacking: Bowls
            // construct simulation
            options->_2psp_iter_max = 45;
            //options->_2psp_tol = math::eps;
            Model* sim = new Model(options, "Scene 12");

            // define robot
            Robot* robot = new Robot();
            dtype mu = 0.6; // friction coefficient

            dtype angle;
            int n_bodies = 10; // number of bodies
            VectorX q(7 * n_bodies);
            VectorX dq(7 * n_bodies);
            std::vector<std::string> file_names = { "../ShapeFiles/bowl/bowl.obj",
                "../ShapeFiles/bowl/bowl_part2.obj",
                "../ShapeFiles/bowl/bowl_part3.obj",
                "../ShapeFiles/bowl/bowl_part4.obj",
                "../ShapeFiles/bowl/bowl_part5.obj",
                "../ShapeFiles/bowl/bowl_part6.obj" };

            for (int i = 0; i < n_bodies; ++i) {
                robot->add_body(new BodyMesh(sim, nullptr, file_names, (dtype)1.0));
                robot->_bodies.back()->_mu = mu; // set friction coefficient
                dtype x = 0;
                dtype y = 0;
                dtype z = (4. * i + 2.) * 0.49 + 0.3;
                if (i % 2 == 1) { 
                    angle = constants::pi / 60.;
                    x = -0.03 * 4.;
                }
                else {
                    angle = -constants::pi / 60.;
                    x = 0.03 * 4.;
                }
                if (i == 0) {
                    angle = 0;
                    x = 0;
                }
                Quat rotation(AngleAxis(angle, Vector3::UnitY()));
                Vector3 pos = Vector3(x, y, z);
                q.segment<7>(7 * i) << rotation.coeffs(), pos; // orientation and position
                dq.segment<7>(7 * i) << 0., 0., 0., 1., 0., 0., 0.; // velocity
            }

            // add ground contact
            //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
            Matrix4 E_g = Matrix4::Identity();

            sim->add_robot(robot);
            sim->set_ground_plane(new BodyPlane(sim, E_g));
            sim->set_state_init(q, dq);
            sim->init();
            return sim;
            break;
        }
        case 13:
        {
            // Circular Tower
            // construct simulation
            options->_2psp_iter_max = 45;
            options->_2psp_tol = math::eps;
            options->_2psp_stable_tol = 1;

            Model* sim = new Model(options, "Scene 13");

            // define robot
            Robot* robot = new Robot();
            Vector3 length(8.0, 14.0, 4.0);

            int layers = 25; // number of layers
            int n_bodies = 5 * layers + 1; // number of bodies
            VectorX q(7 * n_bodies);
            VectorX dq(7 * n_bodies);

            for (int l = 0; l < layers; ++l) {
                for (int i = 0; i < 5; ++i) {
                    robot->add_body(new BodyCuboid(sim, nullptr, length, (dtype)1.0));
                    robot->_bodies.back()->_mu = 0.4; // set friction coefficient
					dtype angle = 2 * constants::pi / 5. * (i+1) + constants::pi / 5. * (l+1);
					dtype x = -4.5 * 4;
					dtype y = 0;
                    dtype z = 4. * l + 2.;
                    Quat rotation(AngleAxis(angle, Vector3::UnitZ()));
                    Vector3 pos = rotation * Vector3(x, y, z);
					q.segment<7>(7 * (5 * l + i)) << rotation.coeffs(), pos; // orientation and position
					dq.segment<7>(7 * (5 * l + i)) << 0., 0., 0., 1., 0., 0., 0.; // velocity
				}
			}

            // add the moving body
            robot->add_body(new BodyCuboid(sim, nullptr, length, (dtype)1.0));
            robot->_bodies.back()->_mu = 0.4; // set friction coefficient
            dtype angle = 2 * constants::pi / 5. + constants::pi / 5. * 24 + constants::pi / 2.;
            dtype x = 0;
            dtype y = 59;
            dtype z = 4. * 24 - 2.;
            Quat rotation(AngleAxis(angle, Vector3::UnitZ()));
            Vector3 pos = rotation * Vector3(x, y, z);
            q.segment<7>(7 * (n_bodies-1)) << rotation.coeffs(), pos; // orientation and position
            Vector3 velocity = rotation * Vector3(0, -450, 20);
            //Vector3 velocity =  Vector3(100, 0, 0);
            dq.segment<7>(7 * (n_bodies-1)) << 0., 0., 0., 1., velocity; // velocity
            

            // add ground contact
            //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
            Matrix4 E_g = Matrix4::Identity();

            sim->add_robot(robot);
            sim->set_ground_plane(new BodyPlane(sim, E_g));
            sim->set_state_init(q, dq);
            sim->init();
            return sim;
            break;
        }
        case 14: 
        {
            // Different height Tower
            // construct simulation
            options->_2psp_iter_max = 75;
            options->_2psp_tol = math::eps_big;
            if(options->_substep > 495)
                options->_substep = 495; // limit the substep due to numerical issue with TGS for this scene
            Model* sim = new Model(options, "Scene 14");

            // define robot
            Robot* robot = new Robot();
            Vector3 length1(4.0, 4.0, 4.0);
            Vector3 length2(4.0, 4.0, 8.0);

            int layers = 2; // number of layers
            int n_bodies = 16 * layers; // number of bodies
            VectorX q(7 * n_bodies);
            VectorX dq(7 * n_bodies);
            dtype mu = 0.5;

            for (int l = 0; l < layers; ++l) {
                for (int i = 0; i < 10; ++i) {
                    robot->add_body(new BodyCuboid(sim, nullptr, length1, (dtype)1.0));
                    robot->_bodies.back()->_mu = mu; // set friction coefficient
                    dtype angle = 0.0;
                    dtype x;
                    dtype y = 0;
                    dtype z;
                    if (l == 0) {
                        x = -7.5 * 4.;
                        z = 4. * i + 2.;
                    }
                    else {
                        x = 7.5 * 4.;
                        z = 4. * i + 2. + 11. * 4.;
                    }

                    Quat rotation(AngleAxis(angle, Vector3::UnitZ()));
                    Vector3 pos = rotation * Vector3(x, y, z);
                    q.segment<7>(7 * (16 * l + i)) << rotation.coeffs(), pos; // orientation and position
                    dq.segment<7>(7 * (16 * l + i)) << 0., 0., 0., 1., 0., 0., 0.; // velocity
                }

                for (int i = 0; i < 5; i++) {
                    robot->add_body(new BodyCuboid(sim, nullptr, length2, (dtype)1.0));
                    robot->_bodies.back()->_mu = mu; // set friction coefficient
                    dtype angle = 0.0;
                    dtype x;
                    dtype y = 0;
                    dtype z;
                    if (l == 0) {
                        x = 7.5 * 4.;
                        z = (4. * i + 2.) * 2.;
                    }
                    else {
                        x = -7.5 * 4.;
                        z = (4. * i + 2.) * 2. + 11. * 4.;
                    }

                    Quat rotation(AngleAxis(angle, Vector3::UnitZ()));
                    Vector3 pos = rotation * Vector3(x, y, z);
                    q.segment<7>(7 * (16 * l + 10 + i)) << rotation.coeffs(), pos; // orientation and position
                    dq.segment<7>(7 * (16 * l + 10 + i)) << 0., 0., 0., 1., 0., 0., 0.; // velocity
                }

                // add the horiziontal body
                robot->add_body(new BodyCuboid(sim, nullptr, Vector3(16 * 4., 4., 4.), (dtype)1.0));
                robot->_bodies.back()->_mu = mu; // set friction coefficient
                dtype angle = 0;
                dtype x = 0;
                dtype y = 0;
                dtype z;
                if (l == 0) {
                    z = 10.5 * 4.;
                }
                else {
                    z = 21.5 * 4.;
                }
                Quat rotation(AngleAxis(angle, Vector3::UnitZ()));
                Vector3 pos = rotation * Vector3(x, y, z);
                q.segment<7>(7 *(16 * l + 15)) << rotation.coeffs(), pos; // orientation and position
                dq.segment<7>(7 * (16 * l + 15)) << 0., 0., 0., 1., 0., 0., 0.; // velocity
            }

            // add ground contact
            //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
            Matrix4 E_g = Matrix4::Identity();

            sim->add_robot(robot);
            sim->set_ground_plane(new BodyPlane(sim, E_g));
            sim->set_state_init(q, dq);
            sim->init();
            return sim;
            break;
        }
        case 15: 
        {
            // Stacking: Drop
            // construct simulation
            options->_2psp_iter_max = 75;
            options->_2psp_tol = math::eps_big;
            options->_h = 0.01;
            Model* sim = new Model(options, "Scene 15");
            // define robot
            Robot* robot = new Robot();
            Vector3 length1(12.0, 18.0, 12.0);
            Vector3 length2(12.0, 18.0, 6.0);
            dtype density = 0.6; // density of the cuboid bodies
            dtype mu = 0.5; // friction coefficient

            int layers = 7; // number of layers
            int n_bodies = 21; // number of bodies
            VectorX q(7 * n_bodies);
            VectorX dq(7 * n_bodies);
            Eigen::VectorXi random_num(28);
            random_num << -45, 15, 11, 9, -124, 2, 0, 13, 37, 11, 0, 15,
                120, 3, 2, 2, -71, 8, 6, 4, 40, 2, 4, 5, -16, 12, 3, 8;
            int body_ind = 0;
            for (int l = 0; l < layers; ++l) {
                Quat rz(AngleAxis(constants::pi/360.*random_num(4*l), Vector3::UnitZ()));
                for (int i = 0; i < 2; ++i) {
                    if (i == 0) {
                        robot->add_body(new BodyCuboid(sim, nullptr, length1, density));
                        robot->_bodies.back()->_mu = mu; // set friction coefficient
                        dtype x = 1.15 * 6;
                        dtype y = 0;
                        dtype z = 6. * (l + 0.65) * 2.65;
                        Quat rotation = AngleAxis(constants::pi / 180. * random_num(4 * l + 1), Vector3::UnitY()) * rz;
                        Vector3 pos = rz * Vector3(x, y, z);
                        q.segment<7>(7 * body_ind) << rotation.coeffs(), pos; // orientation and position
                        dq.segment<7>(7 * body_ind) << 0., 0., 0., 1., 0., 0., 0.; // velocity
                        body_ind++;
                    }
                    else {
                        {
                            robot->add_body(new BodyCuboid(sim, nullptr, length2, density));
                            robot->_bodies.back()->_mu = mu; // set friction coefficient
                            dtype x = -1.15 * 6;
                            dtype y = 0;
                            dtype z = 6. * (l + 0.35) * 2.65;
                            Quat rotation = AngleAxis(constants::pi / 180. * random_num(4 * l + 2), Vector3::UnitY()) * rz;
                            Vector3 pos = rz * Vector3(x, y, z);
                            q.segment<7>(7 * body_ind) << rotation.coeffs(), pos; // orientation and position
                            dq.segment<7>(7 * body_ind) << 0., 0., 0., 1., 0., 0., 0.; // velocity
                            body_ind++;
                        }

                        {
                            robot->add_body(new BodyCuboid(sim, nullptr, length2, density));
                            robot->_bodies.back()->_mu = mu; // set friction coefficient
                            dtype x = -1.25 * 6;
                            dtype y = 0;
                            dtype z = 6. * (l + 0.85) * 2.65;
                            Quat rotation = AngleAxis(constants::pi / 180. * random_num(4 * l + 3), Vector3::UnitY()) * rz;
                            Vector3 pos = rz * Vector3(x, y, z);
                            q.segment<7>(7 * body_ind) << rotation.coeffs(), pos; // orientation and position
                            dq.segment<7>(7 * body_ind) << 0., 0., 0., 1., 0., 0., 0.; // velocity
                            body_ind++;
                        }
                    }
                }
            }

            // add ground contact
            //Matrix3 R = Eigen::AngleAxis<dtype>(-constants::pi / 2., Vector3::UnitX()).matrix();
            Matrix4 E_g = Matrix4::Identity();

            sim->add_robot(robot);
            sim->set_ground_plane(new BodyPlane(sim, E_g));
            sim->set_state_init(q, dq);
            sim->init();
            return sim;
            break;
        }
			default:
				std::cerr << "Unknown SceneId: " << SceneId << std::endl;
				return nullptr;
		}

	}
}