#pragma once
#ifndef _COALMESHMESH_
#define _COALMESHMESH_

#include <Eigen/Dense>

namespace coal {
	struct Contacts
	{
		// Number of contacts
		int count;
		// Maximum penetration depth
		double depthMax;
		// Penetration depths
		double depths[8];
		// Contact points in world space
		Eigen::Vector3d positions[8];
		// Contact normal (same for all points)
		Eigen::Vector3d normal;
	};

	Contacts coalMeshMesh(const Eigen::Matrix4d& M1, const std::shared_ptr<coal::ConvexBase> shape1,
		const Eigen::Matrix4d& M2, const std::shared_ptr<coal::ConvexBase> shape2);

	Contacts coalMeshMesh(const Eigen::Matrix4d& M1, const std::string shape1,
		const Eigen::Matrix4d& M2, const std::string shape2);

}
#endif
