#include "CollisionDetection/CollisionDetection.h"
#include "CollisionDetection/contact.h"
#include "Body/BodyCuboid.h"
#include "Body/BodyPlane.h"
#include "ode/odeBoxBox.h"

namespace _2psp {
	
	bool collision_detection_ground_cuboid(BodyPlane* ground, BodyCuboid* body, std::vector<Collision>& collisions) {
		bool collision = false;
		Vector3 half_length = 0.5 * body->_length;
		MatrixX xl(4, 8);
		xl << -half_length(0), -half_length(0), half_length(0), half_length(0), -half_length(0), -half_length(0), half_length(0), half_length(0),
			half_length(1), -half_length(1), -half_length(1), half_length(1), half_length(1), -half_length(1), -half_length(1), half_length(1),
			-half_length(2), -half_length(2), -half_length(2), -half_length(2), half_length(2), half_length(2), half_length(2), half_length(2),
			1, 1, 1, 1, 1, 1, 1, 1;
		MatrixX xg = ground->_E_i0 * body->_E_0i * xl;
		Vector3 nw = ground->_E_0i.col(2).head(3);

		//std::cout << "xl: \n" << xl << std::endl;

		vector<Contact> contacts;
		contacts.clear();
		for (int i = 0; i < 8; i++) {
			dtype d = xg(2, i);
			if (d < 0.2) {
				collision = true;
				Vector4 xg_i = xg.col(i);
				xg_i(2) = 0;
				// We assume the order of bodies in the collison class is body1(Top) <- body2(Bottom)
				contacts.push_back(Contact(xl.col(i).head(3), (ground->_E_0i * xg_i).head(3), d, nw));
			}
		}
		if (contacts.size() > 0) {
			// We assume the order of bodies in the collision is body1 <- body2
			collisions.push_back(Collision(body, ground, contacts));
		}

		return collision;
	}

	bool collision_detection_cuboid_cuboid(BodyCuboid* cuboid1, BodyCuboid* cuboid2, std::vector<Collision>& collisions) {
		Eigen::Matrix4d E1 = cuboid1->_E_0i.cast<double>();
		Eigen::Matrix4d E2 = cuboid2->_E_0i.cast<double>();
		Eigen::Vector3d s1 = cuboid1->_length.cast<double>();
		Eigen::Vector3d s2 = cuboid2->_length.cast<double>();

		bool collision = false;
		ode::Contacts results = ode::odeBoxBox(E1, s1, E2, s2);

		if (results.count > 0) {
			collision = true;
			bool swap_body = results.normal.dot(Eigen::Vector3d::UnitZ()) > 0;
			std::vector<Contact> contacts;
			for (int i = 0; i < results.count; i++) {
				Vector3d r1 = cuboid1->_E_0i.topLeftCorner(3, 3) * (results.positions[i].cast<dtype>() - cuboid1->_E_0i.topRightCorner(3, 1));
				Vector3d r2 = cuboid2->_E_0i.topLeftCorner(3, 3) * (results.positions[i].cast<dtype>() - static_cast<dtype>(results.depths[i]) * results.normal.cast<dtype>() - cuboid2->_E_0i.topRightCorner(3, 1));
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
				cuboid1->_contact_bodies_current.push_back(cuboid2);
				cuboid2->_contact_bodies_current.push_back(cuboid1);
			}
			else if(swap_body)
				cuboid1->_contact_bodies_next.push_back(cuboid2);
			else
				cuboid2->_contact_bodies_next.push_back(cuboid1);

		}

		return collision;
	}
}