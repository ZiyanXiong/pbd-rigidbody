#pragma once
#include "Common.h"
#include "Utils.h"
#include "Collision.h"

namespace _2psp {
	class Contact;
	class Body;
	class BodyCuboid;
	class BodyPlane;
	class Collision;
	class BodyMesh;

	// detect the collision between ground and a cuboid body
	// return the list of the contact points represented in body frame.
	bool collision_detection_ground_cuboid(BodyPlane* ground, BodyCuboid* body, std::vector<Collision>& collisions);
	bool collision_detection_ground_mesh(BodyPlane* ground, BodyMesh* body, std::vector<Collision>& collisions);

	bool collision_detection_cuboid_cuboid(BodyCuboid* cuboid1, BodyCuboid* cuboid2, std::vector<Collision>& collisions);
	bool collision_detection_mesh_mesh(BodyMesh* body1, BodyMesh* body2, std::vector<Collision>& collisions);

	// detect the collision between a general body and a primitive body
	// the general body should be able to give a list of contact points on the surface.
	// the primitive body is supposed to have an anlytical distance field
	//void collision_detection_general_primitive(const Body* contact_body, const Body* primitive_body, std::vector<Contact>& contacts);

	// detect the collision between a general body and a SDF body
	//void collision_detection_general_SDF(const Body* contact_body, const BodySDFObj* SDF_body, std::vector<Contact>& contacts);
}