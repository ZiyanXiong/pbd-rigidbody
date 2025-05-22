#include "Model.h"
#include "Robot.h"
#include "Body/BodyCuboid.h"
#include "Body/BodyPlane.h"
#include "CollisionDetection/CollisionDetection.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>

namespace _2psp {
	Model::Model(Options* options, std::string name) : 
		_options(options),
		_name(name),
		_ndof_m(0)
	{
		_q_his.clear();
		_dq_his.clear();
		_robots.clear();
	}

	Model::~Model()
	{

	}

	void Model::init() {
		assert(_robots.size() >= 1 && "The number of robots can't be smaller than 1 before initlizating the model");
		for (size_t i = 0; i < _robots.size(); i++) {
			_robots[i]->init(_ndof_m);
			_ndof_m += _robots[i]->_ndof_m;
		}

		assert(_ndof_m == _q_init.size() && "The size of q_init is not equal to ndof_m" );
		assert(_ndof_m == _dq_init.size() && "The size of dq_init is not equal to ndof_m");
		set_state(_q_init, _dq_init);
		update_robot_states_memory();

		if (_ground) {
			_ground_plane->init();
		}
	}

	void Model::set_state_init(const VectorX q_init, const VectorX dq_init) {
		_q_init = q_init;
		_dq_init = dq_init;
	}

	void Model::set_q_init(const VectorX q_init) {
		_q_init = q_init;
	}

	void Model::set_dq_init(const VectorX dq_init) {
		_dq_init = dq_init;
	}

	void Model::set_state(const VectorX q, const VectorX dq) {
		_q = q;
		_dq = dq;
	}

	void Model::set_q(const VectorX q) {
		_q = q;
	}

	void Model::set_dq(const VectorX dq) {
		_dq = dq;
	}

	void Model::update_robot_states_memory() {
		for (size_t i = 0; i < _robots.size(); i++) {
			_robots[i]->set_state_memory(&_q(_robots[i]->_ind_m), &_dq(_robots[i]->_ind_m));
		}
	}

	void Model::update_robot() {
		for (size_t i = 0; i < _robots.size(); i++) {
			_robots[i]->update();
		}
	}

	void Model::set_ground_plane(BodyPlane* ground_plane) {
		_ground_plane = ground_plane;
		_ground = true;
	}

	void Model::collision_detection() {
		// Ground collision
		std::vector<size_t> ground_collision;
		for (size_t i = 0; i < _robots.size(); i++) {
			for (size_t j = 0; j < _robots[i]->_bodies.size(); j++) {
				if(collision_detection_ground_cuboid(_ground_plane, static_cast<BodyCuboid*>(_robots[i]->_bodies[j]), _robots[i]->_collisions))
					ground_collision.push_back(_robots[i]->_collisions.size() - 1);
			}
		}
		_robots[0]->_collision_layer.push_back(ground_collision);

		
		// Free objects collision
		for (size_t i = 0; i < _robots[0]->_bodies.size(); i++) {
			for (size_t j = i + 1; j < _robots[0]->_bodies.size(); j++) {
				collision_detection_cuboid_cuboid(static_cast<BodyCuboid*>(_robots[0]->_bodies[i]), static_cast<BodyCuboid*>(_robots[0]->_bodies[j]), _robots[0]->_collisions);
			}
		}
		_robots[0]->construct_collision_order();
	
	}

	void Model::forward(int num_steps, bool save_history) {
		for(int i = 0; i < num_steps; i++) {
			update_robot();
			//std::cout << "Body States:\n" << _q.transpose() << std::endl;
			collision_detection();
			step_unconstrained();
			temporal_gauss_seidel();
			if(save_history) {
				_q_his.push_back(_q);
				_dq_his.push_back(_dq);
			}
			//std::cout << "Body States:\n" << _q.transpose() << std::endl;
		}
	}

	void Model::step_unconstrained() {
		for (size_t i = 0; i < _robots.size(); i++) {
			_robots[i]->step_unconstrained();
		}
	}

	void Model::temporal_gauss_seidel() {
		dtype hs = _options->_h / _options->_substep;
		for (size_t i = 0; i < _robots.size(); i++) {
			_robots[i]->init_collisions();
		}
		for (int i = 0; i < _options->_substep; i++) {
			for (size_t j = 0; j < _robots.size(); j++) {
				_robots[j]->solve_collisions(hs);
			}
		}
		for (size_t i = 0; i < _robots.size(); i++) {
			_robots[i]->interagate_state();
		}
	}

	void Model::export_replay(std::string folder) {
		if (folder.back() != '/') {
			folder += '/';
		}
		//create folder if it does not exist
		if (!std::filesystem::exists(folder)) {
			std::filesystem::create_directories(folder);
		}
		std::ofstream result_file;
		result_file.open(folder + "Body_States_" + _options->_solver +"_" + std::to_string(_options->_substep) + ".txt");
		size_t body_num = 0;
		for (size_t i = 0; i < _robots.size(); i++) {
			body_num += _robots[i]->_bodies.size();
		}

		if (result_file.is_open()) {
			result_file << "#Body number: " << body_num << std::endl;
			result_file << "#Step number: " << _q_his.size()  << std::endl;
			for (size_t i = 0; i < _q_his.size(); i++) {
			result_file << _q_his[i].transpose() << std::endl;
			}
			result_file.close();
			std::cout << "Results written to file successfully." << std::endl;
		}
		else {
			std::cerr << "Unable to open file." << std::endl;
		}
	}

}