#include "Robot.h"
#include "Body/Body.h"
#include <queue>

namespace _2psp {
	void Robot::init(int ind_m) {
		for (size_t i = 0; i < _bodies.size(); i++) {
			_bodies[i]->init();
		}
		_ind_m = ind_m;
		_ndof_m = 7 * _bodies.size();
	}

	void Robot::add_body(Body* body) {
		_bodies.push_back(body);
	}

	/*
	void Robot::set_q(const VectorX q) {
		for (size_t i = 0; i < _bodies.size() * 6; i += 6) {
			_bodies[i]->set_q(q.segment<6>(i));
		}
	}

	void Robot::set_dq(const VectorX dq) {
		for (size_t i = 0; i < _bodies.size() * 6; i += 6) {
			_bodies[i]->set_q(dq.segment<6>(i));
		}
	}
	*/

	void Robot::update() {
		for (size_t i = 0; i < _bodies.size(); i++) {
			_bodies[i]->update();
		}
	}

	void Robot::set_state_memory(dtype* const qp_start, dtype* const dqp_start) {
		dtype* current_qp = qp_start;
		dtype* current_dqp = dqp_start;
		for (size_t i = 0; i < _bodies.size(); i ++) {
			new (&_bodies[i]->_q) MapVectorX(current_qp, 7);
			new (&_bodies[i]->_dq) MapVectorX(current_dqp, 7);
			current_qp += 7;
			current_dqp += 7;
		}
	}

	void Robot::step_unconstrained() {
		se3 f_t = se3::Zero();
		for (size_t i = 0; i < _bodies.size(); i++) {
			_bodies[i]->step_unconstrained(f_t);
		}
	}

	void Robot::init_collisions() {
		for (size_t i = 0; i < _collisions.size(); i++) {
			_collisions[i].init();
		}
	}

	void Robot::solve_collisions(dtype h) {
		for (size_t layer = 0; layer < _collision_layer.size(); layer++) {
			for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
				_collisions[_collision_layer[layer][i]].solve_nor(h);
				_collisions[_collision_layer[layer][i]].solve_tan(h);
			}
		}
		for (size_t i = 0; i < _bodies.size(); i++) {
			_bodies[i]->update_substep_states(h);
		}
	}

	void Robot::interagate_state() {
		for (size_t i = 0; i < _bodies.size(); i++) {
			_bodies[i]->interagate_state();
		}

		//if (_collisions.size() > 0) {
		//	std::cout << "lambdas:\n" << _collisions[0]._lambdas << std::endl;
		//}

		_collisions.clear();
		_collision_layer.clear();
	}

	void Robot::construct_collision_order() {
		for(size_t i = 0; i < _collision_layer[0].size(); i++) {
			size_t ind = _collision_layer[0][i];
			_collisions[ind]._body1->_layer = 0;
			std::queue<Body*> body_queue;
			body_queue.push(_collisions[ind]._body1);
			while(!body_queue.empty()) {
				Body* current_body = body_queue.front();
				body_queue.pop();
				for (size_t j = 0; j < current_body->_contact_bodies_current.size(); j++) {
					Body* collision_body = current_body->_contact_bodies_current[j];
					if (collision_body->_layer < current_body->_layer) {
						collision_body->_layer = current_body->_layer;
						body_queue.push(collision_body);
					}
				}
				for (size_t j = 0; j < current_body->_contact_bodies_next.size(); j++) {
					Body* collision_body = current_body->_contact_bodies_next[j];
					if (collision_body->_layer < current_body->_layer + 1) {
						collision_body->_layer = current_body->_layer + 1;
						body_queue.push(collision_body);
					}
				}
			}
		}

		for (size_t i = _collision_layer[0].size(); i < _collisions.size(); i++) {
			if( _collisions[i]._body2->_layer > _collision_layer.size())
				for(size_t j = 0; j <_collisions[i]._body2->_layer - _collision_layer.size(); j++)
					_collision_layer.push_back(std::vector<size_t>());
			_collision_layer[_collisions[i]._body2->_layer].push_back(i);
		}
	}
}