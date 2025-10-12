#include "coal/math/transform.h"
#include "coal/mesh_loader/loader.h"
#include "coal/BVH/BVH_model.h"
#include "coal/collision.h"
#include "coal/collision_data.h"
#include "coal/contact_patch.h"
#include "coalMeshMesh.h"

std::shared_ptr<coal::ConvexBase> loadConvexMesh(const std::string& file_name) {
    coal::NODE_TYPE bv_type = coal::BV_AABB;
    coal::MeshLoader loader(bv_type);
    coal::BVHModelPtr_t bvh = loader.load(file_name);
    bvh->buildConvexHull(true, "Qt");
    return bvh->convex;
}

coal::Contacts coal::coalMeshMesh(const Eigen::Matrix4d& M1,
    std::shared_ptr<coal::ConvexBase> shape1,
    const Eigen::Matrix4d& M2,
    std::shared_ptr<coal::ConvexBase> shape2) {

    coal::Transform3s T1;
    T1.setRotation(M1.topLeftCorner(3, 3));
    T1.setTranslation(M1.topRightCorner(3, 1));
    coal::Transform3s T2;
    T2.setRotation(M2.topLeftCorner(3, 3));
    T2.setTranslation(M2.topRightCorner(3, 1));

    coal::CollisionRequest col_req;
    col_req.security_margin = 1e-1;
    coal::CollisionResult col_res;

    // Collision call
    coal::collide(shape1.get(), T1, shape2.get(), T2, col_req, col_res);

    coal::ContactPatchRequest patch_req;
    coal::ContactPatchResult patch_res;
    patch_req.setPatchTolerance(5e-2);
    patch_req.setNumSamplesCurvedShapes(8);
    coal::computeContactPatch(shape1.get(), T1, shape2.get(), T2, col_res,
        patch_req, patch_res);

    coal::Contacts results;
    if (patch_res.numContactPatches() > 0 && col_res.isCollision()) {
        coal::ContactPatch contactpatch = patch_res.getContactPatch(0);

        results.depthMax = contactpatch.penetration_depth;
        results.count = contactpatch.size();
        if (results.count > 8)
            results.count = 8;

        for (size_t i = 0; i < contactpatch.size() && i < 8; ++i) {
            results.positions[i] = contactpatch.getPoint(i);
            results.depths[i] = contactpatch.penetration_depth;
        }
        results.normal << contactpatch.getNormal();
    }

    col_res.clear();

    return results;
}

coal::Contacts coal::coalMeshMesh(const Eigen::Matrix4d& M1,
    const std::string shape1_path,
    const Eigen::Matrix4d& M2,
    const std::string shape2_path) {

    std::shared_ptr<coal::ConvexBase> shape1 = loadConvexMesh(shape1_path);
    std::shared_ptr<coal::ConvexBase> shape2 = loadConvexMesh(shape2_path);

    coal::Transform3s T1;
    T1.setRotation(M1.topLeftCorner(3, 3));
    T1.setTranslation(M1.topRightCorner(3, 1));
    coal::Transform3s T2;
    T2.setRotation(M2.topLeftCorner(3, 3));
    T2.setTranslation(M2.topRightCorner(3, 1));

    coal::CollisionRequest col_req;
    col_req.security_margin = 1e-1;
    coal::CollisionResult col_res;

    // Collision call
    coal::collide(shape1.get(), T1, shape2.get(), T2, col_req, col_res);

    coal::ContactPatchRequest patch_req;
    coal::ContactPatchResult patch_res;
    patch_req.setPatchTolerance(5e-2);
    patch_req.setNumSamplesCurvedShapes(8);
    coal::computeContactPatch(shape1.get(), T1, shape2.get(), T2, col_res,
        patch_req, patch_res);

    coal::Contacts results;
    if (patch_res.numContactPatches() > 0 && col_res.isCollision()) {
        coal::ContactPatch contactpatch = patch_res.getContactPatch(0);

        results.depthMax = contactpatch.penetration_depth;
        results.count = contactpatch.size();
        if (results.count > 8)
            results.count = 8;

        for (size_t i = 0; i < contactpatch.size() && i < 8; ++i) {
            results.positions[i] = contactpatch.getPoint(i);
            results.depths[i] = contactpatch.penetration_depth;
        }
        results.normal << contactpatch.getNormal();
    }

    col_res.clear();
    
    return results;
}