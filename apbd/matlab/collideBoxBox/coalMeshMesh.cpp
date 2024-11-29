// Command for compiling
// mex COPTIMFLAGS='-O3 -DNDEBUG' -I"C:\Users\dsdsx\anaconda3\envs\pyPBD\Library\include" -I"C:\Users\dsdsx\anaconda3\envs\pyPBD\Library\include\eigen3" -L"C:\Users\dsdsx\anaconda3\envs\pyPBD\Library\lib" -lcoal coalMeshMesh.cpp
#pragma warning(disable : 4996)

#include "coal/math/transform.h"
#include "coal/mesh_loader/loader.h"
#include "coal/BVH/BVH_model.h"
#include "coal/collision.h"
#include "coal/collision_data.h"
#include "coal/contact_patch.h"

#include <iostream>
#include <memory>

//#define DEBUG_MAIN
#define MATLAB_MEX_BUILD

/** If instricuted, compile a mex function for Matlab.  */
#ifdef MATLAB_MEX_BUILD
#include "mex.h"
#define A(i, j) A[i + j * M]
#else
#define mexPrintf printf
#endif

// Function to load a convex mesh from a `.obj`, `.stl` or `.dae` file.
//
// This function imports the object inside the file as a BVHModel, i.e. a point
// cloud which is hierarchically transformed into a tree of bounding volumes.
// The leaves of this tree are the individual points of the point cloud
// stored in the `.obj` file.
// This BVH can then be used for collision detection.
//
// For better computational efficiency, we sometimes prefer to work with
// the convex hull of the point cloud. This insures that the underlying object
// is convex and thus very fast collision detection algorithms such as
// GJK or EPA can be called with this object.
// Consequently, after creating the BVH structure from the point cloud, this
// function also computes its convex hull.
std::shared_ptr<coal::ConvexBase> loadConvexMesh(const std::string& file_name) {
  coal::NODE_TYPE bv_type = coal::BV_AABB;
  coal::MeshLoader loader(bv_type);
  coal::BVHModelPtr_t bvh = loader.load(file_name);
  bvh->buildConvexHull(true, "Qt");
  return bvh->convex;
}

#ifdef DEBUG_MAIN
    int main() {
      // Create the coal shapes.
      // Coal supports many primitive shapes: boxes, spheres, capsules, cylinders,
      // ellipsoids, cones, planes, halfspace and convex meshes (i.e. convex hulls
      // of clouds of points). It also supports BVHs (bounding volumes hierarchies),
      // height-fields and octrees.
      std::shared_ptr<coal::ConvexBase> shape1 = loadConvexMesh("./cube1.obj");
      std::shared_ptr<coal::ConvexBase> shape2 = loadConvexMesh("./cube1.obj");

      // Define the shapes' placement in 3D space
      coal::Transform3s T1;
      T1.setQuatRotation(coal::Quaternion3f::Identity());
      T1.setTranslation(coal::Vec3s(0.0, 0.0, 0.0));
      coal::Transform3s T2 = coal::Transform3s::Identity();
      // T2.setQuatRotation(coal::Quaternion3f::UnitRandom());
      T2.setTranslation(coal::Vec3s(0.2, 0.0, 2.1));

      // Define collision requests and results.
      //
      // The collision request allows to set parameters for the collision pair.
      // For example, we can set a positive or negative security margin.
      // If the distance between the shapes is less than the security margin, the
      // shapes will be considered in collision. Setting a positive security margin
      // can be usefull in motion planning, i.e to prevent shapes from getting too
      // close to one another. In physics simulation, allowing a negative security
      // margin may be usefull to stabilize the simulation.
      coal::CollisionRequest col_req;
      col_req.security_margin = 1e-1;
      // A collision result stores the result of the collision test (signed distance
      // between the shapes, witness points location, normal etc.)
      coal::CollisionResult col_res;

      // Collision call
      coal::collide(shape1.get(), T1, shape2.get(), T2, col_req, col_res);

      coal::ContactPatchRequest patch_req;
      coal::ContactPatchResult patch_res;
      coal::computeContactPatch(shape1.get(), T1, shape2.get(), T2, col_res,
                                patch_req, patch_res);

      // We can access the collision result once it has been populated
      std::cout << "Collision? " << col_res.isCollision() << "\n";
      if (col_res.isCollision()) {
        coal::Contact contact = col_res.getContact(0);
        // The penetration depth does **not** take into account the security margin.
        // Consequently, the penetration depth is the true signed distance which
        // separates the shapes. To have the distance which takes into account the
        // security margin, we can simply add the two together.
        std::cout << "Penetration depth: " << contact.penetration_depth << "\n";
        std::cout << "Distance between the shapes including the security margin: "
                  << contact.penetration_depth + col_req.security_margin << "\n";
        std::cout << "Witness point on shape1: "
                  << contact.nearest_points[0].transpose() << "\n";
        std::cout << "Witness point on shape2: "
                  << contact.nearest_points[1].transpose() << "\n";
        std::cout << "Normal: " << contact.normal.transpose() << "\n";
      }

      // We can access the collision result once it has been populated
      std::cout << "Contact patch number: " << patch_res.numContactPatches()
                << "\n";
      if (patch_res.numContactPatches() > 0 && col_res.isCollision()) {
        coal::ContactPatch contactpatch = patch_res.getContactPatch(0);

        std::cout << "Penetration depth: " << contactpatch.penetration_depth
                  << "\n";
        std::cout << "Distance between the shapes including the security margin: "
                  << contactpatch.penetration_depth + col_req.security_margin
                  << "\n";
        for (size_t i = 0; i < contactpatch.size(); ++i) {
          std::cout << "Witness point on shape1: "
                    << (contactpatch.getPoint(i) +
                        0.5 * contactpatch.penetration_depth *
                            contactpatch.getNormal())
                           .transpose()
                    << "\n";
          std::cout << "Witness point on shape2: "
                    << (contactpatch.getPoint(i) -
                        0.5 * contactpatch.penetration_depth *
                            contactpatch.getNormal())
                           .transpose()
                    << "\n";
        }

        std::cout << "Normal: " << contactpatch.getNormal().transpose() << "\n";
      }

      // Before calling another collision test, it is important to clear the
      // previous results stored in the collision result.
      col_res.clear();

      return 0;
    }
#endif  // DEBUG


#ifdef MATLAB_MEX_BUILD
    struct Contacts {
      // Number of contacts
      int count = 0;
      // Maximum penetration depth
      double depthMax = 0;
      // Penetration depths
      double depths[8];
      // Contact points in world space
      Eigen::Vector3d positions[8];
      // Contact normal (same for all points)
      Eigen::Vector3d normal = Eigen::Vector3d(0,0,0);
    };

    Contacts coalMeshMesh(const Eigen::Matrix4d& M1,
                          const std::string &meshPath1,
                          const Eigen::Matrix4d& M2,
                          const std::string &meshPath2) {

      std::shared_ptr<coal::ConvexBase> shape1 = loadConvexMesh(meshPath1);
      std::shared_ptr<coal::ConvexBase> shape2 = loadConvexMesh(meshPath2);

      coal::Transform3s T1;
      T1.setRotation(M1.topLeftCorner(3,3));
      T1.setTranslation(M1.topRightCorner(3,1));
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

      Contacts results;
      if (patch_res.numContactPatches() > 0 && col_res.isCollision()) {
        coal::ContactPatch contactpatch = patch_res.getContactPatch(0);

        results.depthMax = contactpatch.penetration_depth;
        results.count = contactpatch.size();
        if(results.count > 8)
          results.count = 8;

        for (size_t i = 0; i < contactpatch.size() && i < 8; ++i) {
          results.positions[i] = contactpatch.getPoint(i);
          results.depths[i] = contactpatch.penetration_depth;
        }
        results.normal << contactpatch.getNormal();
      }

      return results;
    }

    void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                     const mxArray *prhs[]) {
      // check for proper number of arguments
      if (nrhs != 4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "4 inputs required.");
      }
      if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "1 output required.");
      }

      enum { RHS_E1 = 0, RHS_MPATH1, RHS_E2, RHS_MPATH2, RHS_COUNT };

      if (!mxIsDouble(prhs[RHS_E1]) ||
          mxGetNumberOfElements(prhs[RHS_E1]) != 16) {
        mexErrMsgTxt("E1 must be a mat4.");
      }
      if (!mxIsChar(prhs[RHS_MPATH1])) {
        mexErrMsgTxt("meshPath1 must be a valid file path.");
      }
      if (!mxIsDouble(prhs[RHS_E2]) ||
          mxGetNumberOfElements(prhs[RHS_E2]) != 16) {
        mexErrMsgTxt("E2 must be a mat4.");
      }
      if (!mxIsChar(prhs[RHS_MPATH2])) {
        mexErrMsgTxt("meshPath2 must be a valid file path.");
      }

      mwSize M;
      mwSize N;
      mxChar *pth;
      double *A;
      char file_path[120];

      // Convert E1
      M = mxGetM(prhs[RHS_E1]);
      A = mxGetPr(prhs[RHS_E1]);
      Eigen::Matrix4d E1;
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          E1(i, j) = A(i, j);
        }
      }

      // Convert mesh path 1
      N = mxGetN(prhs[RHS_MPATH1]);
      pth = mxGetChars(prhs[RHS_MPATH1]);
      for (int i = 0; i < N && i < 120; ++i) {
        file_path[i] = pth[i];
      }
      file_path[N] = '\0';
      std::string meshPath1(file_path);
      //mexPrintf(file_path);
      //mexPrintf(" File Name length %d.\n", N);

      // Convert E2
      M = mxGetM(prhs[RHS_E2]);
      A = mxGetPr(prhs[RHS_E2]);
      Eigen::Matrix4d E2;
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          E2(i, j) = A(i, j);
        }
      }

      // Convert mesh path 1
      N = mxGetN(prhs[RHS_MPATH2]);
      pth = mxGetChars(prhs[RHS_MPATH2]);
      for (int i = 0; i < N && i < 120; ++i) {
        file_path[i] = pth[i];
      }
      file_path[N] = '\0';
      std::string meshPath2(file_path);

      // Call collision detector
      Contacts contacts = coalMeshMesh(E1, meshPath1, E2, meshPath2);

      // Convert collisions
      int nfields = 5;
      const char **fnames;
      fnames = (const char **)mxCalloc(nfields, sizeof(*fnames));
      fnames[0] = "count";
      fnames[1] = "depthMax";
      fnames[2] = "nor";
      fnames[3] = "pos";
      fnames[4] = "depth";
      mxArray *c = mxCreateStructMatrix(1, 1, nfields, fnames);
      mxSetField(c, 0, "count", mxCreateDoubleScalar(contacts.count));
      mxSetField(c, 0, "depthMax", mxCreateDoubleScalar(contacts.depthMax));
      mxArray *mat;
      mat = mxCreateDoubleMatrix(3, 1, mxREAL);
      M = mxGetM(mat);
      A = mxGetPr(mat);
      for (int i = 0; i < 3; ++i) {
        A(i, 0) = contacts.normal(i);
      }
      mxSetField(c, 0, "nor", mat);
      //mexPrintf("contacts.count: %d\n", contacts.count);
      mat = mxCreateDoubleMatrix(3, contacts.count, mxREAL);
      M = mxGetM(mat);
      A = mxGetPr(mat);
      for (int k = 0; k < contacts.count; ++k) {
        for (int i = 0; i < 3; ++i) {
          A(i, k) = contacts.positions[k](i);
        }
      }
      mxSetField(c, 0, "pos", mat);
      //mexPrintf("contacts.count: %d\n", contacts.count);
      mat = mxCreateDoubleMatrix(1, contacts.count, mxREAL);
      M = mxGetM(mat);
      A = mxGetPr(mat);
      for (int k = 0; k < contacts.count; ++k) {
        A(0, k) = contacts.depths[k];
      }
      mxSetField(c, 0, "depth", mat);
      mxFree((void *)fnames);
      plhs[0] = c;
      //mexPrintf("Finish.\n");
    }
#endif