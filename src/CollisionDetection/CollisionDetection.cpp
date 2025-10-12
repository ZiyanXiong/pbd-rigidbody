#include "CollisionDetection/CollisionDetection.h"
#include "CollisionDetection/contact.h"
#include "Body/BodyCuboid.h"
#include "Body/BodyPlane.h"
#include "Body/BodyMesh.h"
#include "ode/odeBoxBox.h"
#include "coal/coalMeshMesh.h"

namespace _2psp {
	
	bool collision_detection_ground_cuboid(BodyPlane* ground, BodyCuboid* body, std::vector<Collision>& collisions) {
		bool collision = false;
		Vector3 half_length = 0.5 * body->_length;
		MatrixX xl(4, 8);
		xl << -half_length(0), -half_length(0), -half_length(0), -half_length(0), half_length(0), half_length(0), half_length(0), half_length(0),
			-half_length(1), -half_length(1), half_length(1), half_length(1), -half_length(1), -half_length(1), half_length(1), half_length(1),
			-half_length(2), half_length(2), -half_length(2), half_length(2), -half_length(2), half_length(2), -half_length(2), half_length(2),
			1, 1, 1, 1, 1, 1, 1, 1;
		MatrixX xg = ground->_E_0i * body->_E_i0 * xl;
		Vector3 nw = ground->_E_i0.col(2).head(3);

		//std::cout << "Ei0_g: \n" << ground->_E_i0 << std::endl;
		//std::cout << "E0i_g: \n" << ground->_E_0i << std::endl;
		//std::cout << "x_g: \n" << xg << std::endl;

		vector<Contact> contacts;
		contacts.clear();
		for (int i = 0; i < 8; i++) {
			dtype d = xg(2, i);
			if (d < 0.2) {
				collision = true;
				Vector4 xg_i = xg.col(i);
				xg_i(2) = 0;
				// We assume the order of bodies in the collison class is body1(Top) <- body2(Bottom)
				contacts.push_back(Contact(xl.col(i).head(3), xg_i.head(3), d, nw));
			}
		}
		if (contacts.size() > 0) {
			// We assume the order of bodies in the collision is body1 <- body2
			collisions.push_back(Collision(body, ground, contacts));
		}

		return collision;
	}

	bool collision_detection_cuboid_cuboid(BodyCuboid* cuboid1, BodyCuboid* cuboid2, std::vector<Collision>& collisions) {
		Eigen::Matrix4d E1 = cuboid1->_E_i0.cast<double>();
		Eigen::Matrix4d E2 = cuboid2->_E_i0.cast<double>();
		Eigen::Vector3d s1 = cuboid1->_length.cast<double>();
		Eigen::Vector3d s2 = cuboid2->_length.cast<double>();

		bool collision = false;
		ode::Contacts results = ode::odeBoxBox(E1, s1, E2, s2);
		//if (results.count == 1) {
		//	std::cout << "Collision Detection: " << results.count << std::endl;
		//	std::cout << "E1: " << E1 << std::endl;
		//	std::cout << "E2: " << E2 << std::endl;
		//	std::cout << "results: \n normal: " << results.normal.transpose() << std::endl;
		//	std::cout << "depth: " << results.depths[0] << std::endl;
		//	std::cout << "pos: " << results.positions[0] << std::endl;
		//	ode::Contacts results_new = ode::odeBoxBox(E1, s1, E2, s2);
		//}

		if (results.count > 0) {
			collision = true;
			bool swap_body = results.normal.dot(Eigen::Vector3d::UnitZ()) > 0;
			if (cuboid1->_is_infinite_mass && cuboid2->_is_infinite_mass) {
				// If both bodies are infinite mass, we can't to solve the collision
				return false;
			}
			if(cuboid1->_is_infinite_mass)
				swap_body = true; // If body1 is infinite mass, we always swap the body

			std::vector<Contact> contacts;
			for (int i = 0; i < results.count; i++) {
				Vector3d r1 = E1.topLeftCorner(3, 3).transpose() * (results.positions[i].cast<dtype>() - E1.topRightCorner(3, 1));
				Vector3d r2 = E2.topLeftCorner(3, 3).transpose() * (results.positions[i].cast<dtype>() - static_cast<dtype>(results.depths[i]) * results.normal.cast<dtype>() - E2.topRightCorner(3, 1));
				// The normal from conllider is pointng from body1 to body2
				// We assume the order of bodies in the collison class is body1(Top) <- body2(Bottom)
				if (swap_body)
					contacts.push_back(Contact(r2, r1, static_cast<dtype>(results.depths[i]), results.normal.cast<dtype>()));
				else
					contacts.push_back(Contact(r1, r2, static_cast<dtype>(results.depths[i]), -results.normal.cast<dtype>()));
			}

			if (swap_body)
				collisions.push_back(Collision(cuboid2, cuboid1, contacts));
			else
				collisions.push_back(Collision(cuboid1, cuboid2, contacts));

			if (abs(results.normal.dot(Eigen::Vector3d::UnitZ())) < 0.25)
			{
				collisions.back()._shock_porpagate = false; // If the normal is close to vertical, we don't use shock porpagation
			}
			else if (swap_body) {
				cuboid1->_contact_bodies_next.push_back(cuboid2);
				collisions.back()._shock_porpagate = true;
			}
			else {
				cuboid2->_contact_bodies_next.push_back(cuboid1);
				collisions.back()._shock_porpagate = true;
			}

			//for (size_t i = 0; i < contacts.size(); i++)
			//{
			//	std::cout << "contact : " << i;
			//	std::cout << " depth: " << contacts[i]._d;
			//	std::cout << " normal: " << contacts[i]._normal.transpose();
			//	std::cout << " x1: " << contacts[i]._xi1.transpose();
			//	std::cout << " x2: " << contacts[i]._xi2.transpose() << std::endl;
			//}

		}

		return collision;
	}

	bool collision_detection_ground_mesh(BodyPlane* ground, BodyMesh* body, std::vector<Collision>& collisions) {
		bool collision = false;
		Matrix3X& xl = body->_V;
		Matrix4 E_ig = ground->_E_0i * body->_E_i0;
		Matrix3 R_ig = E_ig.topLeftCorner(3, 3);
		Vector3 p_ig = E_ig.topRightCorner(3, 1);
		Matrix3X xg = (R_ig * xl).colwise() + p_ig;
		Vector3 nw = ground->_E_i0.col(2).head(3);

		//std::cout << "Ei0_g: \n" << ground->_E_i0 << std::endl;
		//std::cout << "E0i_g: \n" << ground->_E_0i << std::endl;
		//std::cout << "x_g: \n" << xg << std::endl;
		dtype depth_max = xg.row(2).minCoeff();
		vector<Contact> contacts;
		contacts.clear();
		std::vector<int> indices;
		if (depth_max < 0.2) {
			collision = true;
			//std::cout << "Collision Detection: " << depth_max << std::endl;
			for (int i = 0; i < xg.cols(); i++) {
				if (xg(2, i) < depth_max + 5e-2) {
					indices.push_back(i);
				}
			}
		};

		if (indices.size() > 8)
		{
			std::vector<int> indices_new(4, 0);
			std::vector<dtype> min_max_value(4, 0);
			min_max_value[0] = min_max_value[1] = xg(0, 0);
			min_max_value[2] = min_max_value[3] = xg(1, 0);
			for (int i = 0; i < indices.size(); i++) {
				if (xg(0, indices[i]) < min_max_value[0])
				{
					min_max_value[0] = xg(0, indices[i]);
					indices_new[0] = indices[i];
				}

				if (xg(0, indices[i]) > min_max_value[1])
				{
					min_max_value[1] = xg(0, indices[i]);
					indices_new[1] = indices[i];
				}

				if (xg(1, indices[i]) < min_max_value[2])
				{
					min_max_value[2] = xg(1, indices[i]);
					indices_new[2] = indices[i];
				}
				if (xg(1, indices[i]) > min_max_value[3])
				{
					min_max_value[3] = xg(1, indices[i]);
					indices_new[3] = indices[i];
				}
			}
			indices = indices_new;
		}
		if (indices.size() > 0) {
			for (int i = 0; i < indices.size(); i++) {
				Vector3 xg_i = xg.col(indices[i]);
				xg_i(2) = 0;
				// We assume the order of bodies in the collison class is body1(Top) <- body2(Bottom)
				contacts.push_back(Contact(xl.col(indices[i]), xg_i, xg(2, indices[i]), nw));
			}
			// We assume the order of bodies in the collision is body1 <- body2
			collisions.push_back(Collision(body, ground, contacts));
		}
		//for (size_t i = 0; i < contacts.size(); i++)
		//{
		//	std::cout << "contact : " << i;
		//	std::cout << " depth: " << contacts[i]._d;
		//	std::cout << " normal: " << contacts[i]._normal.transpose();
		//	std::cout << " x1: " << contacts[i]._xi1.transpose();
		//	std::cout << " x2: " << contacts[i]._xi2.transpose() << std::endl;
		//}
		return collision;
	}

	bool collision_detection_mesh_mesh(BodyMesh* body1, BodyMesh* body2, std::vector<Collision>& collisions) {
		Eigen::Matrix4d E1 = body1->_E_i0.cast<double>();
		Eigen::Matrix4d E2 = body2->_E_i0.cast<double>();

		bool collision = false;
		std::vector<Contact> contacts;
		if (body1->_convex_shapes.size() > 1 && body2->_convex_shapes.size() > 1) {
			for (size_t i = 1; i < body2->_convex_shapes.size(); i+=2)
			{
				coal::Contacts results = coal::coalMeshMesh(E1 * body1->_E_io, body1->_convex_shapes[i], E2 * body2->_E_io, body2->_convex_shapes[0]);
				//coal::Contacts results = coal::coalMeshMesh(E1 * body1->_E_io, body1->_file_names[i], E2 * body2->_E_io, body2->_file_names[0]);
				for (int j = 0; j < min(results.count, 2); j++) {
					Vector3d r1 = E1.topLeftCorner(3, 3).transpose() * (results.positions[j].cast<dtype>() - static_cast<dtype>(results.depths[j] * 0.5) * results.normal.cast<dtype>() - E1.topRightCorner(3, 1));
					Vector3d r2 = E2.topLeftCorner(3, 3).transpose() * (results.positions[j].cast<dtype>() + static_cast<dtype>(results.depths[j] * 0.5) * results.normal.cast<dtype>() - E2.topRightCorner(3, 1));
					// The normal from conllider is pointng from body1 to body2
					// We assume the order of bodies in the collison class is body1(Top) <- body2(Bottom)
					contacts.push_back(Contact(r2, r1, static_cast<dtype>(results.depths[j]), results.normal.cast<dtype>()));
					//std::cout << "Collision Detection: " << results.count << std::endl;
					//std::cout << "E1: " << E1 << std::endl;
					//std::cout << "E2: " << E2 << std::endl;
					//std::cout << "E1_io: " << body1->_E_io << std::endl;
					//std::cout << "E2_io: " << body2->_E_io << std::endl;
					//std::cout << "results: \n normal: " << results.normal.transpose() << std::endl;
					//std::cout << "depth: " << results.depths[0] << std::endl;
					//std::cout << "pos: " << results.positions[0].transpose() << std::endl;
				}
			}
		}
		if (contacts.size() > 0) {
			collision = true;
			collisions.push_back(Collision(body2, body1, contacts));
			body1->_contact_bodies_next.push_back(body2);
			collisions.back()._shock_porpagate = true;
			//std::cout << "Colliding body: " << body1->_index << ", " << body2->_index << ", contact num: " << contacts.size() << std::endl;
		}

		//for (size_t i = 0; i < contacts.size(); i++)
		//{
		//	std::cout << "contact : " << i;
		//	std::cout << " depth: " << contacts[i]._d;
		//	std::cout << " normal: " << contacts[i]._normal.transpose();
		//	std::cout << " x1: " << contacts[i]._xi1.transpose();
		//	std::cout << " x2: " << contacts[i]._xi2.transpose() << std::endl;
		//}
		return collision;
	}
}