#include "Robot.h"
#include "Body/Body.h"
#include <queue>

namespace _2psp {
	void Robot::init(int ind_m) {
		for (size_t i = 0; i < _bodies.size(); i++) {
			_bodies[i]->init();
			_bodies[i]->_index = i;
		}
		_ind_m = ind_m;
		_ndof_m = 7 * _bodies.size();
		_n_c = 0; // Reset the number of constraints
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

	void Robot::reset() {
		for (size_t i = 0; i < _bodies.size(); i++) {
			_bodies[i]->_phi = _bodies[i]->_phi_0;
			_bodies[i]->_phi_dt.setZero();
			_bodies[i]->_delta_phi.setZero();
		}
		for (size_t i = 0; i < _collisions.size(); i++) {
			_collisions[i]._lambdas.setZero();
		}
	}

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

	void Robot::step_unconstrained(VectorX& f_t) {
		for (size_t i = 0; i < _bodies.size(); i++) {
			MapVectorX f_t_i(f_t.data() + i * 6, 6);
			_bodies[i]->step_unconstrained(f_t_i);

			//std::cout << "Body " << i << " Force: " << f_t_i.transpose() << std::endl;
		}
	}

	void Robot::init_collisions() {
		_n_c = 0; // Reset the number of constraints
		for (size_t i = 0; i < _collisions.size(); i++) {
			_collisions[i].init();
			_n_c += _collisions[i]._contacts.size();
		}
	}

	void Robot::solve_collisions(dtype h) {
		for (size_t layer = 0; layer < _collision_layer.size(); layer++) {
			for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
				_collisions[_collision_layer[layer][i]].solve_nor(h);
			}
			for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
				_collisions[_collision_layer[layer][i]].solve_tan(h);
			}
		}

		for (size_t i = 0; i < _bodies.size(); i++) {
			_bodies[i]->update_substep_states(h);
		}

		//for (size_t i = 0; i < _bodies.size(); i++) {
		//	std::cout << "Body " << i << " States Phi:\n" << _bodies[i]->_phi.transpose() << std::endl;
		//	std::cout << "Body " << i << " States Phi_dt:\n" << _bodies[i]->_phi_dt.transpose() << std::endl;
		//}
		//std::cout << "Lambdas: ";
		//for (size_t i = 0; i < _collisions.size(); i++) {
		//	std::cout << _collisions[i]._lambdas.transpose();
		//}
		//std::cout << std::endl;
	}

	bool Robot::solve_collisions_2psp(dtype h, int& solve_count, int sp_iter_max, dtype tol, dtype stable_tol) {
		int up_iter_max = sp_iter_max;
		int down_iter_max = sp_iter_max;
		bool upward_success = true;
		bool downward_success = true;

		solve_count = 0;

		for (size_t layer = 0; layer < _collision_layer.size(); layer++) {
			for (int iter = 0; iter < up_iter_max; iter++) {
				for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
					if (_collisions[_collision_layer[layer][i]]._shock_porpagate) {
						_collisions[_collision_layer[layer][i]].solve_nor_2psp(h);
						_collisions[_collision_layer[layer][i]].apply_mass_averaged_impulse_nor();
					}
					else {
						_collisions[_collision_layer[layer][i]].solve_nor(h);
					}
					solve_count += _collisions[_collision_layer[layer][i]]._contacts.size();
				}

				for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
					if(_collisions[_collision_layer[layer][i]]._shock_porpagate)
						_collisions[_collision_layer[layer][i]].solve_tan_2psp(h);
					else
						_collisions[_collision_layer[layer][i]].solve_tan(h);
				}

				bool is_converged = true;
				//for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
				//	is_converged &= _collisions[_collision_layer[layer][i]]._is_converged;
				//	_collisions[_collision_layer[layer][i]]._is_converged = true;
				//}

				for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
					is_converged = is_converged && _collisions[_collision_layer[layer][i]].is_converged(h, tol);
					//is_converged &= _collisions[_collision_layer[layer][i]].is_converged(h);
				}

				if (layer == 1 && iter > 100) {
					//std::cout << "Layer " << layer << " Iteration " << iter << " is_converged: " << is_converged << std::endl;
					//_collisions[_collision_layer[layer][0]].comute_c(h);
					//std::cout << "rs_normal 1:\n" << _collisions[_collision_layer[layer][0]]._c << std::endl;
					//_collisions[_collision_layer[layer][1]].comute_c(h);
					//std::cout << "rs_normal 2:\n" << _collisions[_collision_layer[layer][1]]._c << std::endl;

					//std::cout << "lambda 1:\n" << _collisions[_collision_layer[layer][0]]._lambdas << std::endl;
					//std::cout << "lambda 2:\n" << _collisions[_collision_layer[layer][1]]._lambdas << std::endl;
					//std::cout << "dlambda_nor 1:\n" << _collisions[_collision_layer[layer][0]]._dlambdas_nor << std::endl;
					//std::cout << "dlambda_nor 2:\n" << _collisions[_collision_layer[layer][1]]._dlambdas_nor << std::endl;

					//std::cout << "dlambda_nor 1:\n" << _collisions[_collision_layer[layer][0]]._dlambdas_nor << std::endl;
					//std::cout << "dlambda_nor 2:\n" << _collisions[_collision_layer[layer][1]]._dlambdas_nor << std::endl;

					//std::cout << "J1:\n" << _collisions[_collision_layer[layer][0]]._J1 << _collisions[_collision_layer[layer][1]]._J1 << std::endl;
					//std::cout << "_J_div_m1:\n" << _collisions[_collision_layer[layer][0]]._J_div_m1 << _collisions[_collision_layer[layer][1]]._J_div_m1 << std::endl;


					//std::cout << "body phi1:\n" << _collisions[_collision_layer[layer][1]]._body1->_phi << std::endl;
					//std::cout << "body phi2:\n" << _collisions[_collision_layer[layer][1]]._body1->_phi << std::endl;
				}

				if (is_converged)
				{
					//std::cout << "iter: " << iter << std::endl;
					break;
				}
			}

			for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
				if (_collisions[_collision_layer[layer][i]]._shock_porpagate) {
					_collisions[_collision_layer[layer][i]].compute_c(h);
					if ((_collisions[_collision_layer[layer][i]]._c.row(0).array() > stable_tol).any()) {
						//std::cout << "rs_normal :" << _collisions[_collision_layer[layer][i]]._c.row(0) << std::endl;
						//std::cout << "body1 layer:" << _collisions[_collision_layer[layer][i]]._body1->_layer << std::endl;
						//std::cout << "body2 layer:" << _collisions[_collision_layer[layer][i]]._body2->_layer << std::endl;
						upward_success = false;
					}
				}
			}


			//for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
			//	upward_success &= _collisions[_collision_layer[layer][i]]._is_sp_valid;
			//	_collisions[_collision_layer[layer][i]]._is_sp_valid = true;
			//}

			if (!upward_success)
				break;
		}


		//std::cout << "Iteration number :" << solve_count / constraint_count << std::endl;
		//for (size_t i = 0; i < _bodies.size(); i++) {
		//	std::cout << "Body " << i << " States Phi:\n" << _bodies[i]->_phi.transpose() << std::endl;
		//	std::cout << "Body " << i << " States Phi_dt:\n" << _bodies[i]->_phi_dt.transpose() << std::endl;
		//}
		
		if (upward_success) {
			for (int layer = _collision_layer.size() - 1; layer >= 0; layer--) {
				for (int iter = 0; iter < up_iter_max; iter++) {
					for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
						if (_collisions[_collision_layer[layer][i]]._shock_porpagate) {
							_collisions[_collision_layer[layer][i]].solve_nor_2psp(h);
							_collisions[_collision_layer[layer][i]].apply_mass_averaged_impulse_nor();
						}
						else {
							_collisions[_collision_layer[layer][i]].solve_nor(h);
						}
						solve_count += _collisions[_collision_layer[layer][i]]._contacts.size();
					}

					for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
						if (_collisions[_collision_layer[layer][i]]._shock_porpagate)
							_collisions[_collision_layer[layer][i]].solve_tan_2psp(h);
						else
							_collisions[_collision_layer[layer][i]].solve_tan(h);
					}

					bool is_converged = true;
					for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
						is_converged &= _collisions[_collision_layer[layer][i]]._is_converged;
						_collisions[_collision_layer[layer][i]]._is_converged = true;
					}

					for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
						is_converged = is_converged && _collisions[_collision_layer[layer][i]].is_converged(h, tol);
						//is_converged &= _collisions[_collision_layer[layer][i]].is_converged(h);
					}

					if (is_converged)
					{
						//std::cout << "iter: " << iter << std::endl;
						break;
					}
				}


				for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
					if (_collisions[_collision_layer[layer][i]]._shock_porpagate) {
						_collisions[_collision_layer[layer][i]].compute_c(h);
						if ((_collisions[_collision_layer[layer][i]]._c.row(0).array() > stable_tol).any()) {
							//std::cout << "rs_normal :" << _collisions[_collision_layer[layer][i]]._c.row(0) << std::endl;
							//std::cout << "body1 layer:" << _collisions[_collision_layer[layer][i]]._body1->_layer << std::endl;
							//std::cout << "body1 phi:" << _collisions[_collision_layer[layer][i]]._body1->_phi.transpose() << std::endl;
							//std::cout << "body2 layer:" << _collisions[_collision_layer[layer][i]]._body2->_layer << std::endl;
							//std::cout << "body2 phi:" << _collisions[_collision_layer[layer][i]]._body2->_phi.transpose() << std::endl;
							downward_success = false;
						}
					}
				}


				if (!downward_success)
					break;

				for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
					if(_collisions[_collision_layer[layer][i]]._shock_porpagate)
						_collisions[_collision_layer[layer][i]].apply_accumulated_impulse();
				}
			}
		}
		

		//int iter_total = solve_count / constraint_count;
		//std::cout << "Iteration number :" << iter_total << std::endl;

		//for (size_t i = 0; i < _bodies.size(); i++) {
		//	std::cout << "Body " << i << " States Phi:\n" << _bodies[i]->_phi.transpose() << std::endl;
		//	std::cout << "Body " << i << " States Phi_dt:\n" << _bodies[i]->_phi_dt.transpose() << std::endl;
		//}

		for (size_t i = 0; i < _bodies.size(); i++) {
			_bodies[i]->update_substep_states(h);
		}
		//return true;
		return upward_success && downward_success;
	}

	void Robot::solve_velocity(dtype h) {
		for (size_t layer = 0; layer < _collision_layer.size(); layer++) {
			for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
				_collisions[_collision_layer[layer][i]].solve_vel_nor(h);
			}
			for (size_t i = 0; i < _collision_layer[layer].size(); i++) {
				_collisions[_collision_layer[layer][i]].solve_vel_tan(h);
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
			_collisions[ind]._body1->_layer = 1;
			_collisions[ind]._body2->_layer = 0;
			//std::cout << "body1:" << _collisions[ind]._body1->_index << ", body2:" << _collisions[ind]._body2->_index << std::endl;
			//std::cout << "body1 layer:" << _collisions[ind]._body1->_layer << ", body2 layer:" << _collisions[ind]._body2->_layer << std::endl;
			std::queue<Body*> body_queue;
			body_queue.push(_collisions[ind]._body1);
			while(!body_queue.empty()) {
				Body* current_body = body_queue.front();
				body_queue.pop();

				for (size_t j = 0; j < current_body->_contact_bodies_next.size(); j++) {
					Body* collision_body = current_body->_contact_bodies_next[j];
					if (collision_body->_layer < current_body->_layer + 1) {
						collision_body->_layer = current_body->_layer + 1;
						body_queue.push(collision_body);
					}
				}
			}
		}

		int max_layer = 0;
		for (size_t i = 0; i < _bodies.size(); i++) {
			if (_bodies[i]->_layer > max_layer) {
				max_layer = _bodies[i]->_layer;
			}
		}
		if (max_layer < 1) {
			max_layer = 1; // Ensure at least one layer exists
		}

		for (size_t i = 0; i < max_layer - 1; i++) {
			_collision_layer.push_back(std::vector<size_t>()); // Initialize each layer with an empty vector
		}	

		for (size_t i = 0; i < _collisions.size(); i++) {
			if (_collisions[i]._body2->_layer == 0) {
				continue;
			}

			if (_collisions[i]._body1->_layer < 0) 
				_collisions[i]._body1->_layer = max_layer;
			if (_collisions[i]._body2->_layer < 0)
				_collisions[i]._body2->_layer = max_layer;

			if (_collisions[i]._body1->_layer == _collisions[i]._body2->_layer)
				_collisions[i]._shock_porpagate = false;

			if(_collisions[i]._shock_porpagate)
				_collision_layer[_collisions[i]._body1->_layer-1].push_back(i);
			else {
				//std::cout << "Collision index:" << i << " Body1 layer :" << _collisions[i]._body1->_layer << " Body2 layer :" << _collisions[i]._body2->_layer << std::endl;
				_collision_layer[min(_collisions[i]._body1->_layer, _collisions[i]._body2->_layer) - 1].push_back(i);
			}
		}
	}
}