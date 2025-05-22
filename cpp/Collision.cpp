#include "Collision.h"
#include "Body/Body.h"

namespace _2psp {

	Collision::Collision(Body* body1, Body* body2, vector<Contact> contacts) : _body1(body1), _body2(body2) {
		_contacts = contacts;
		size_t n = _contacts.size();

		_w1 = Matrix3X::Zero(3,n);
		_w2 = Matrix3X::Zero(3,n);
		_d = Matrix3X::Zero(3,n);
		_J1 = JacobianMatrixVector(3, 6, n);
		_J2 = JacobianMatrixVector(3, 6, n);
		_J_div_m1 = JacobianMatrixVector(3, 6, n);
		_J_div_m2 = JacobianMatrixVector(3, 6, n);
		_lambdas = Matrix3X::Zero(3,n);
		_mu = 0.5 * (_body1->_mu + _body2->_mu);
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


			_d.col(i) = xw1 - xw2;
			_J1(i) = math::get_tangent(_contacts[i]._normal).transpose() * math::gamma(_body1->_E_i0.topLeftCorner(3, 3) * _contacts[i]._xi1);
			_J2(i) =_J1(i).topRightCorner(3, 3) * math::gamma(_body2->_E_i0.topLeftCorner(3, 3) * _contacts[i]._xi2);


			_J_div_m1(i).topRightCorner(3, 3) = _J1(i).topRightCorner(3, 3).array() * _body1->_mass_inv;
			_J_div_m2(i).topRightCorner(3, 3) = _J2(i).topRightCorner(3, 3).array() * _body2->_mass_inv;

			rxn_l = _body1->_E_0i.topLeftCorner(3, 3)* _J1(i).topLeftCorner(3,3).transpose();
			_w1.col(i) = (rxn_l.transpose() * Iinv1.asDiagonal() * rxn_l).diagonal() + Vector3::Ones() * _body1->_mass_inv;
			_J_div_m1(i).topLeftCorner(3, 3) = (_body1->_E_i0.topLeftCorner(3, 3) * Iinv1.asDiagonal() * rxn_l).transpose();

			rxn_l = _body2->_E_0i.topLeftCorner(3, 3) * _J2(i).topLeftCorner(3, 3).transpose();
			_w2.col(i) = (rxn_l.transpose() * Iinv2.asDiagonal() * rxn_l).diagonal() + Vector3::Ones() * _body2->_mass_inv;
			_J_div_m2(i).topLeftCorner(3, 3) = (_body2->_E_i0.topLeftCorner(3, 3) * Iinv2.asDiagonal() * rxn_l).transpose();
		}

		//std::cout << "w1:\n" << _w1 << std::endl;
		//std::cout << "w2:\n" << _w2 << std::endl;
		//std::cout << "J1:\n" << _J1 << std::endl;
		//std::cout << "J2:\n" << _J2 << std::endl;
		//std::cout << "I_inv1:\n" << Iinv1 << std::endl;
		//std::cout << "I_inv2:\n" << Iinv2 << std::endl;
		//std::cout << "Contact Frame:\n" << math::get_tangent(_contacts[0]._normal) << std::endl;
	}

	void Collision::solve_nor(dtype h) {
		//std::cout << "Body1 Phi:\n" << _body1->_phi.transpose() << std::endl;
		//std::cout << "Body2 Phi:\n" << _body2->_phi.transpose() << std::endl;
		//std::cout << "Body1 Phi_dt:\n" << _body1->_phi_dt.transpose() << std::endl;
		//std::cout << "Body2 Phi_dt:\n" << _body2->_phi_dt.transpose() << std::endl;
		for (size_t i = 0; i < _contacts.size(); ++i) {
			dtype dlambda_nor;
			dtype c;
			c = _J1(i).row(0).dot(_body1->_phi) - _J2(i).row(0).dot(_body2->_phi);
			c += _J1(i).row(0).dot(_body1->_phi_dt / h) - _J2(i).row(0).dot(_body2->_phi_dt / h) + _J1(i).row(0).tail(3).dot(_d.col(i));
			dlambda_nor = -c / (_w1.col(i)(0) + _w2.col(i)(0));
			if(_lambdas.col(i)(0) + dlambda_nor < 0) {
				dlambda_nor = -_lambdas.col(i)(0);
			}
			_lambdas.col(i)(0) += dlambda_nor;
			_body1->_phi += _J_div_m1(i).row(0).transpose() * dlambda_nor;
			_body2->_phi -= _J_div_m2(i).row(0).transpose() * dlambda_nor;
		}
	}
	/*
	void Collision::solve_tan() {
		for (size_t i = 0; i < _contacts.size(); ++i) {
			Vector2 dlambda_tan;
			Vector2 c;
			c = _J1(i)(seq(1, 2), all) * _body1->_phi - _J2(i)(seq(1, 2), all) * _body2->_phi;
			c -= _J1(i)(seq(1, 2), all) * _body1->_phi_dt - _J2(i)(seq(1, 2), all) * _body2->_phi_dt + _J1(i)(seq(1, 2), seq(3, 5)) * _d.col(i);
			dlambda_tan = -c.array() / (_w1.col(i).tail(2) + _w2.col(i).tail(2)).array();

			Vector2 dlambda_tan_max = (_lambdas.col(i).tail(2) + dlambda_tan);
			dtype lambda_tan_max_norm = dlambda_tan_max.norm();
			if (lambda_tan_max_norm > _mu * _lambdas.col(i)(0)) {
				dlambda_tan_max = _mu * _lambdas.col(i)(0) * dlambda_tan_max / lambda_tan_max_norm;
				dlambda_tan = dlambda_tan_max - _lambdas.col(i).tail(2);
			}
			_lambdas.col(i).tail(2) += dlambda_tan;
			_body1->_phi += _J_div_m1(i).bottomRows(2).transpose() * dlambda_tan;
			_body2->_phi -= _J_div_m2(i).bottomRows(2).transpose() * dlambda_tan;
		}
	}
	*/
	void Collision::solve_tan(dtype h) {
		for (size_t i = 0; i < _contacts.size(); ++i) {
			dtype dlambda_tan1, dlambda_tan2, c1, c2, lambda_tan_norm;
			c1 = _J1(i).row(1).dot(_body1->_phi) - _J2(i).row(1).dot(_body2->_phi);
			c1 += _J1(i).row(1).dot(_body1->_phi_dt / h) - _J2(i).row(1).dot(_body2->_phi_dt / h) + _J1(i).row(1).tail(3).dot(_d.col(i));
			c2 = _J1(i).row(2).dot(_body1->_phi) - _J2(i).row(2).dot(_body2->_phi);
			c2 += _J1(i).row(2).dot(_body1->_phi_dt / h) - _J2(i).row(2).dot(_body2->_phi_dt / h) + _J1(i).row(2).tail(3).dot(_d.col(i));

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