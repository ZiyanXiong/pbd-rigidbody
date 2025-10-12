#include "Collision.h"
#include "Body/Body.h"

namespace _2psp {

	Collision::Collision(Body* body1, Body* body2, vector<Contact> contacts) : _body1(body1), _body2(body2) {
		_contacts = contacts;
		size_t n = _contacts.size();

		_w1 = Matrix3X::Zero(3,n);
		_w2 = Matrix3X::Zero(3,n);
		_d = Matrix3X::Zero(3,n);
		_c = Matrix3X::Zero(3,n);
		_c0 = Matrix3X::Zero(3, n);
		_J1 = JacobianMatrixVector(3, 6, n);
		_J2 = JacobianMatrixVector(3, 6, n);
		_J_div_m1 = JacobianMatrixVector(3, 6, n);
		_J_div_m2 = JacobianMatrixVector(3, 6, n);
		_lambdas_prev = Matrix3X::Zero(3, n);
		_lambdas = Matrix3X::Zero(3,n);
		_dlambdas_nor = VectorX::Zero(n);

		_mu = 0.5 * (_body1->_mu + _body2->_mu);
		_shock_porpagate = true;
		_is_converged = true;
		_is_sp_valid = true;
	}

	void Collision::init() {
		Matrix3 rxn_l = Matrix3::Zero();
		Vector3 Iinv1 = _body1->_Inertia_inv.head(3);
		Vector3 Iinv2 =_body2->_Inertia_inv.head(3);

		//std::cout << "E1:\n" << _body1->_E_0i << std::endl;
		//std::cout << "E2:\n" << _body2->_E_0i << std::endl;

		for (size_t i = 0; i < _contacts.size(); ++i) {
			Vector3 xw1 = _body1->transform_point(_contacts[i]._xi1);
			Vector3 xw2 = _body2->transform_point(_contacts[i]._xi2);

			//std::cout << "xw1:" << xw1.transpose() << " xl1:" << _contacts[i]._xi1.transpose() << std::endl;
			//std::cout << "xw2:" << xw2.transpose() << " xl2:" << _contacts[i]._xi2.transpose() << std::endl;
			//std::cout << "raxnw:" << math::get_tangent(_contacts[i]._normal).transpose() * math::gamma(_body1->_E_i0.topLeftCorner(3, 3) * _contacts[i]._xi1) << std::endl;
			//std::cout << "raxnl:" << _body1->_E_i0.topLeftCorner(3,3)*((_body1->_E_0i.topLeftCorner(3,3)*math::get_tangent(_contacts[i]._normal)).transpose() * math::gamma(_contacts[i]._xi1)) << std::endl;	
			
			_d.col(i) = xw1 - xw2;
			//std::cout << "d_i:\n" << _d.col(i).transpose() << std::endl;
			_J1(i) = math::get_tangent(_contacts[i]._normal).transpose() * math::gamma(_body1->_E_i0.topLeftCorner(3, 3) * _contacts[i]._xi1);
			_J2(i) =_J1(i).topRightCorner(3, 3) * math::gamma(_body2->_E_i0.topLeftCorner(3, 3) * _contacts[i]._xi2);


			_J_div_m1(i).topRightCorner(3, 3) = _J1(i).topRightCorner(3, 3).array() * _body1->_mass_inv;
			_J_div_m2(i).topRightCorner(3, 3) = _J2(i).topRightCorner(3, 3).array() * _body2->_mass_inv;

			rxn_l = math::skew(_contacts[i]._xi1) * (_body1->_E_0i.topLeftCorner(3, 3) * _J1(i).topRightCorner(3,3).transpose());

			//std::cout << "rxn_l:\n" << rxn_l << std::endl;
			//std::cout << "nl:" << (_body1->_E_0i.topLeftCorner(3, 3) * _J1(i).topRightCorner(3, 3).transpose()) << std::endl;
			//std::cout << "rxnl_0:\n" << math::skew(_contacts[i]._xi1) * (_body1->_E_0i.topLeftCorner(3, 3) * _J1(i).row(0).tail(3).transpose()) << std::endl;
			//std::cout << "E_0i*Ei0:\n" << _body1->_E_0i.topLeftCorner(3, 3) * _body1->_E_i0.topLeftCorner(3,3) << std::endl;

			_w1.col(i) = (rxn_l.transpose() * Iinv1.asDiagonal() * rxn_l).diagonal() + Vector3::Ones() * _body1->_mass_inv;
			_J_div_m1(i).topLeftCorner(3, 3) = (_body1->_E_i0.topLeftCorner(3, 3) * Iinv1.asDiagonal() * rxn_l).transpose();

			rxn_l = math::skew(_contacts[i]._xi2) * (_body2->_E_0i.topLeftCorner(3, 3) * _J2(i).topRightCorner(3,3).transpose());
			_w2.col(i) = (rxn_l.transpose() * Iinv2.asDiagonal() * rxn_l).diagonal() + Vector3::Ones() * _body2->_mass_inv;
			_J_div_m2(i).topLeftCorner(3, 3) = (_body2->_E_i0.topLeftCorner(3, 3) * Iinv2.asDiagonal() * rxn_l).transpose();
		}

		//std::cout << "w1:\n" << _w1 << std::endl;
		//std::cout << "w2:\n" << _w2 << std::endl;
		//std::cout << "J1:\n" << _J1 << std::endl;
		//std::cout << "J1/M1:\n" << _J_div_m1 << std::endl;
		//std::cout << "J2/M2:\n" << _J_div_m2 << std::endl;
		//std::cout << "I_inv1:\n" << Iinv1 << std::endl;
		//std::cout << "I_inv2:\n" << Iinv2 << std::endl;
		//std::cout << "Contact Frame:\n" << math::get_tangent(_contacts[0]._normal) << std::endl;
	}

	void Collision::compute_c(dtype h) {
		for (size_t i = 0; i < _contacts.size(); ++i) {
			_c.col(i) = _J1(i) * _body1->_phi - _J2(i) * _body2->_phi + _J1(i).topRightCorner(3,3) * _d.col(i) / h;
		}
	}

	bool Collision::is_converged(dtype h, dtype tol) {
		//std::cout << "delta lambdas:\n" << _lambdas - _lambdas_prev << std::endl;
		//return (_lambdas-_lambdas_prev).norm() < math::eps;
		//_c0 = _c;
		//compute_c(h);
		//std::cout << "_c:\n" << _c.transpose()  << std::endl;
		//std::cout << "dc_norm:\n" << (_c - _c0).norm() << std::endl;
		return (_c-_c0).norm() < tol;
	}

	void Collision::solve_nor(dtype h) {
		//std::cout << "Body1 Phi:\n" << _body1->_phi.transpose() << std::endl;
		//std::cout << "Body2 Phi:\n" << _body2->_phi.transpose() << std::endl;
		//std::cout << "Body1 Phi_dt:\n" << _body1->_phi_dt.transpose() << std::endl;
		//std::cout << "Body2 Phi_dt:\n" << _body2->_phi_dt.transpose() << std::endl;
		//std::cout << "d:\n" << _d << std::endl;
		for (size_t i = 0; i < _contacts.size(); ++i) {
			dtype dlambda_nor;
			dtype c;
			c = _J1(i).row(0).dot(_body1->_phi + _body1->_phi_dt / h) - _J2(i).row(0).dot(_body2->_phi + _body2->_phi_dt / h) + _J1(i).row(0).tail(3).dot(_d.col(i)) / h;
			//std::cout << "_J2(i).row(0): " << _J2(i).row(0) << std::endl;
			//std::cout << "Body2 Phi: " << _body2->_phi.transpose() << std::endl;
			//std::cout << "Body2 Phi_dt: " << _body2->_phi_dt.transpose() << std::endl;
			//std::cout << "c: "<< c << " = " << _J1(i).row(0).dot(_body1->_phi + _body1->_phi_dt / h) << " - " << _J2(i).row(0).dot(_body2->_phi + _body2->_phi_dt / h) << " + " << _J1(i).row(0).tail(3).dot(_d.col(i)) / h << std::endl;
			//c += _J1(i).row(0).dot(_body1->_phi_dt / h) - _J2(i).row(0).dot(_body2->_phi_dt / h) + _J1(i).row(0).tail(3).dot(_d.col(i)) / h;

			dlambda_nor = -c / (_w1.col(i)(0) + _w2.col(i)(0));
			if(_lambdas.col(i)(0) + dlambda_nor < 0) {
				dlambda_nor = -_lambdas.col(i)(0);
			}
			_lambdas.col(i)(0) += dlambda_nor;
			_body1->_phi += _J_div_m1(i).row(0).transpose() * dlambda_nor;
			_body2->_phi -= _J_div_m2(i).row(0).transpose() * dlambda_nor;

			//dtype dlambda_tan1, dlambda_tan2, c1, c2, lambda_tan_norm;
			//c1 = _J1(i).row(1).dot(_body1->_phi + _body1->_phi_dt / h) - _J2(i).row(1).dot(_body2->_phi + _body2->_phi_dt / h) + _J1(i).row(1).tail(3).dot(_d.col(i)) / h;
			////c1 += _J1(i).row(1).dot(_body1->_phi_dt / h) - _J2(i).row(1).dot(_body2->_phi_dt / h) + _J1(i).row(1).tail(3).dot(_d.col(i)) / h;
			//c2 = _J1(i).row(2).dot(_body1->_phi + _body1->_phi_dt / h) - _J2(i).row(2).dot(_body2->_phi + _body2->_phi_dt / h) + _J1(i).row(2).tail(3).dot(_d.col(i)) / h;
			////c2 += _J1(i).row(2).dot(_body1->_phi_dt / h) - _J2(i).row(2).dot(_body2->_phi_dt / h) + _J1(i).row(2).tail(3).dot(_d.col(i)) / h;

			//dlambda_tan1 = -c1 / (_w1.col(i)(1) + _w2.col(i)(1));
			//dlambda_tan2 = -c2 / (_w1.col(i)(2) + _w2.col(i)(2));
			//Vector2 lambda_tan(dlambda_tan1 + _lambdas.col(i)(1), dlambda_tan2 + _lambdas.col(i)(2));
			//lambda_tan_norm = lambda_tan.norm();
			//if (lambda_tan_norm > _mu * _lambdas.col(i)(0) && lambda_tan_norm > math::eps_big) {
			//	dlambda_tan1 = _mu * _lambdas.col(i)(0) * lambda_tan(0) / lambda_tan_norm - _lambdas.col(i)(1);
			//	dlambda_tan2 = _mu * _lambdas.col(i)(0) * lambda_tan(1) / lambda_tan_norm - _lambdas.col(i)(2);
			//}
			//_lambdas.col(i)(1) += dlambda_tan1;
			//_lambdas.col(i)(2) += dlambda_tan2;
			////std::cout << "lambdas :\n" << _lambdas.col(i) << std::endl;
			//_body1->_phi += _J_div_m1(i).row(1).transpose() * dlambda_tan1;
			//_body1->_phi += _J_div_m1(i).row(2).transpose() * dlambda_tan2;
			//_body2->_phi -= _J_div_m2(i).row(1).transpose() * dlambda_tan1;
			//_body2->_phi -= _J_div_m2(i).row(2).transpose() * dlambda_tan2;
		}
	}

	void Collision::solve_tan(dtype h) {
		for (size_t i = 0; i < _contacts.size(); ++i) {
			dtype dlambda_tan1, dlambda_tan2, c1, c2, lambda_tan_norm;
			c1 = _J1(i).row(1).dot(_body1->_phi + _body1->_phi_dt / h) - _J2(i).row(1).dot(_body2->_phi + _body2->_phi_dt / h) + _J1(i).row(1).tail(3).dot(_d.col(i)) / h;
			//c1 += _J1(i).row(1).dot(_body1->_phi_dt / h) - _J2(i).row(1).dot(_body2->_phi_dt / h) + _J1(i).row(1).tail(3).dot(_d.col(i)) / h;
			c2 = _J1(i).row(2).dot(_body1->_phi + _body1->_phi_dt / h) - _J2(i).row(2).dot(_body2->_phi + _body2->_phi_dt / h) + _J1(i).row(2).tail(3).dot(_d.col(i)) / h;
			//c2 += _J1(i).row(2).dot(_body1->_phi_dt / h) - _J2(i).row(2).dot(_body2->_phi_dt / h) + _J1(i).row(2).tail(3).dot(_d.col(i)) / h;

			dlambda_tan1 = -c1 / (_w1.col(i)(1) + _w2.col(i)(1));
			dlambda_tan2 = -c2 / (_w1.col(i)(2) + _w2.col(i)(2));
			Vector2 lambda_tan(dlambda_tan1 + _lambdas.col(i)(1), dlambda_tan2 + _lambdas.col(i)(2));
			lambda_tan_norm = lambda_tan.norm();
			if (lambda_tan_norm > _mu * _lambdas.col(i)(0) && lambda_tan_norm > math::eps_big) {
				dlambda_tan1 = _mu * _lambdas.col(i)(0) * lambda_tan(0) / lambda_tan_norm - _lambdas.col(i)(1);
				dlambda_tan2 = _mu * _lambdas.col(i)(0) * lambda_tan(1) / lambda_tan_norm - _lambdas.col(i)(2);
			}
			_lambdas.col(i)(1) += dlambda_tan1;
			_lambdas.col(i)(2) += dlambda_tan2;
			//std::cout << "lambdas :\n" << _lambdas.col(i) << std::endl;
			_body1->_phi += _J_div_m1(i).row(1).transpose() * dlambda_tan1;
			_body1->_phi += _J_div_m1(i).row(2).transpose() * dlambda_tan2;
			_body2->_phi -= _J_div_m2(i).row(1).transpose() * dlambda_tan1;
			_body2->_phi -= _J_div_m2(i).row(2).transpose() * dlambda_tan2;
		}
	}

	void Collision::solve_nor_2psp(dtype h) {
		for (size_t i = 0; i < _contacts.size(); ++i) {
			dtype dlambda_nor;
			dtype c;
			c = _J1(i).row(0).dot(_body1->_phi) - _J2(i).row(0).dot(_body2->_phi);
			c += _J1(i).row(0).tail(3).dot(_d.col(i)) / h;
			_c0.col(i)(0) = _c.col(i)(0);
			_c.col(i)(0) = c;
			dlambda_nor = -c / _w1(0,i);
			//if (_lambdas.col(i)(0) + dlambda_nor < 0) {
			//	dlambda_nor = -_lambdas.col(i)(0);
			//}
			//_lambdas.col(i)(0) += dlambda_nor;
			//_body1->_phi += _J_div_m1(i).row(0).transpose() * dlambda_nor;
			//std::cout << "Body1 Phi:\n" << _body1->_phi.transpose() << std::endl;
			_dlambdas_nor(i) = dlambda_nor;

			//if (abs(c)>math::eps_big)
			//	_is_converged = false;

			//if (_shock_porpagate && c > 1)
			//	_is_sp_valid = false;

			//_body2->_phi -= _J_div_m2(i).row(0).transpose() * dlambda_nor;
		}
	}

	void Collision::apply_mass_averaged_impulse_nor() {
		//dtype w_sum = 0;
		//for (int i = 0; i < _contacts.size(); ++i) {
		//	w_sum += _w1(0,i);
		//}

		for (size_t i = 0; i < _contacts.size(); ++i) {
			// Apply the impulse to body1
			dtype dlambda_nor = _dlambdas_nor(i) * 0.5;
			if (_lambdas.col(i)(0) + dlambda_nor < 0) {
				dlambda_nor = -_lambdas.col(i)(0);
			}
			_lambdas.col(i)(0) += dlambda_nor;
			_body1->_phi += _J_div_m1(i).row(0).transpose() * dlambda_nor;
		}
	}


	void Collision::solve_tan_2psp(dtype h) {
		for (size_t i = 0; i < _contacts.size(); ++i) {
			dtype dlambda_tan1, dlambda_tan2, c1, c2, lambda_tan_norm;
			c1 = _J1(i).row(1).dot(_body1->_phi) - _J2(i).row(1).dot(_body2->_phi);
			c1 +=  _J1(i).row(1).tail(3).dot(_d.col(i)) / h;
			c2 = _J1(i).row(2).dot(_body1->_phi) - _J2(i).row(2).dot(_body2->_phi);
			c2 +=  _J1(i).row(2).tail(3).dot(_d.col(i)) / h;
			_c0.col(i)(1) = _c.col(i)(1);
			_c0.col(i)(2) = _c.col(i)(2);
			_c.col(i)(1) = c1;
			_c.col(i)(2) = c2;

			dlambda_tan1 = -c1 / _w1(1,i);
			dlambda_tan2 = -c2 / _w1(2,i);
			Vector2 lambda_tan(dlambda_tan1 + _lambdas.col(i)(1), dlambda_tan2 + _lambdas.col(i)(2));
			lambda_tan_norm = lambda_tan.norm();
			if (lambda_tan_norm > _mu * _lambdas.col(i)(0) && lambda_tan_norm > math::eps) {
				dlambda_tan1 = _mu * _lambdas.col(i)(0) * lambda_tan(0) / lambda_tan_norm - _lambdas.col(i)(1);
				dlambda_tan2 = _mu * _lambdas.col(i)(0) * lambda_tan(1) / lambda_tan_norm - _lambdas.col(i)(2);
			}
			//_lambdas_prev.col(i) = _lambdas.col(i);
			_lambdas.col(i)(1) += dlambda_tan1;
			_lambdas.col(i)(2) += dlambda_tan2;
			_body1->_phi += _J_div_m1(i).row(1).transpose() * dlambda_tan1;
			_body1->_phi += _J_div_m1(i).row(2).transpose() * dlambda_tan2;
			//std::cout << "Body1 Phi:\n" << _body1->_phi.transpose() << std::endl;

			//if (abs(c1) > math::eps_big || abs(c2) > math::eps_big)
			//	_is_converged = false;

			//_body2->_phi -= _J_div_m2(i).row(1).transpose() * dlambda_tan1;
			//_body2->_phi -= _J_div_m2(i).row(2).transpose() * dlambda_tan2;
		}
	}

	void Collision::apply_accumulated_impulse() {
		for (size_t i = 0; i < _contacts.size(); ++i) {
			_body2->_phi -= _J_div_m2(i).transpose() * _lambdas.col(i);
		}
		//std::cout << "lambdas:\n" << _lambdas << std::endl;
	}

	void Collision::solve_vel_nor(dtype h) {
		for (size_t i = 0; i < _contacts.size(); ++i) {
			dtype dlambda_nor;
			dtype c;
			c = _J1(i).row(0).dot(_body1->_phi) - _J2(i).row(0).dot(_body2->_phi);
			//c = _J1(i).row(0).dot(_body1->_phi + _body1->_phi_dt / h) - _J2(i).row(0).dot(_body2->_phi + _body2->_phi_dt / h);
			dlambda_nor = -c / (_w1.col(i)(0) + _w2.col(i)(0));
			if (_lambdas.col(i)(0) + dlambda_nor < 0) {
				dlambda_nor = -_lambdas.col(i)(0);
			}
			_lambdas.col(i)(0) += dlambda_nor;
			_body1->_phi += _J_div_m1(i).row(0).transpose() * dlambda_nor;
			_body2->_phi -= _J_div_m2(i).row(0).transpose() * dlambda_nor;

			//dtype dlambda_tan1, dlambda_tan2, c1, c2, lambda_tan_norm;
			//c1 = _J1(i).row(1).dot(_body1->_phi) - _J2(i).row(1).dot(_body2->_phi);
			//c2 = _J1(i).row(2).dot(_body1->_phi) - _J2(i).row(2).dot(_body2->_phi);
			////c1 = _J1(i).row(1).dot(_body1->_phi + _body1->_phi_dt / h) - _J2(i).row(1).dot(_body2->_phi + _body2->_phi_dt / h);
			////c2 = _J1(i).row(2).dot(_body1->_phi + _body1->_phi_dt / h) - _J2(i).row(2).dot(_body2->_phi + _body2->_phi_dt / h);

			//dlambda_tan1 = -c1 / (_w1.col(i)(1) + _w2.col(i)(1));
			//dlambda_tan2 = -c2 / (_w1.col(i)(2) + _w2.col(i)(2));
			//Vector2 lambda_tan(dlambda_tan1 + _lambdas.col(i)(1), dlambda_tan2 + _lambdas.col(i)(2));
			//lambda_tan_norm = lambda_tan.norm();
			//if (lambda_tan_norm > _mu * _lambdas.col(i)(0)) {
			//	dlambda_tan1 = _mu * _lambdas.col(i)(0) * lambda_tan(0) / lambda_tan_norm - _lambdas.col(i)(1);
			//	dlambda_tan2 = _mu * _lambdas.col(i)(0) * lambda_tan(1) / lambda_tan_norm - _lambdas.col(i)(2);
			//}
			//_lambdas.col(i)(1) += dlambda_tan1;
			//_lambdas.col(i)(2) += dlambda_tan2;
			//_body1->_phi += _J_div_m1(i).row(1).transpose() * dlambda_tan1;
			//_body1->_phi += _J_div_m1(i).row(2).transpose() * dlambda_tan2;
			//_body2->_phi -= _J_div_m2(i).row(1).transpose() * dlambda_tan1;
			//_body2->_phi -= _J_div_m2(i).row(2).transpose() * dlambda_tan2;
		}
	}

	void Collision::solve_vel_tan(dtype h) {
		for (size_t i = 0; i < _contacts.size(); ++i) {
			dtype dlambda_tan1, dlambda_tan2, c1, c2, lambda_tan_norm;
			c1 = _J1(i).row(1).dot(_body1->_phi) - _J2(i).row(1).dot(_body2->_phi);
			c2 = _J1(i).row(2).dot(_body1->_phi) - _J2(i).row(2).dot(_body2->_phi);
			//c1 = _J1(i).row(1).dot(_body1->_phi + _body1->_phi_dt / h) - _J2(i).row(1).dot(_body2->_phi + _body2->_phi_dt / h);
			//c2 = _J1(i).row(2).dot(_body1->_phi + _body1->_phi_dt / h) - _J2(i).row(2).dot(_body2->_phi + _body2->_phi_dt / h);

			dlambda_tan1 = -c1 / (_w1.col(i)(1) + _w2.col(i)(1));
			dlambda_tan2 = -c2 / (_w1.col(i)(2) + _w2.col(i)(2));
			Vector2 lambda_tan(dlambda_tan1 + _lambdas.col(i)(1), dlambda_tan2 + _lambdas.col(i)(2));
			lambda_tan_norm = lambda_tan.norm();
			if (lambda_tan_norm > _mu * _lambdas.col(i)(0)) {
				dlambda_tan1 = _mu * _lambdas.col(i)(0) * lambda_tan(0) / lambda_tan_norm - _lambdas.col(i)(1);
				dlambda_tan2 = _mu * _lambdas.col(i)(0) * lambda_tan(1) / lambda_tan_norm - _lambdas.col(i)(2);
			}
			_lambdas.col(i)(1) += dlambda_tan1;
			_lambdas.col(i)(2) += dlambda_tan2;
			_body1->_phi += _J_div_m1(i).row(1).transpose() * dlambda_tan1;
			_body1->_phi += _J_div_m1(i).row(2).transpose() * dlambda_tan2;
			_body2->_phi -= _J_div_m2(i).row(1).transpose() * dlambda_tan1;
			_body2->_phi -= _J_div_m2(i).row(2).transpose() * dlambda_tan2;
		}
	}
}