#pragma once
#include "Body/Body.h"
#include "coal/mesh_loader/loader.h"
#include "coal/BVH/BVH_model.h"

namespace _2psp {

    class BodyMesh : public Body {
    public:
        std::string _filename;  // .obj filename (path)
        Matrix3X _V;            // vertices
        Matrix3Xi _F;            // face elements
        bool _is_infinite_mass = false; // if true, the body has infinite mass and inertia
        SE3 _E_io;

        std::vector<std::shared_ptr<coal::ConvexBase>> _convex_shapes;
        std::vector<std::string> _file_names;
        BodyMesh(Model* sim, Joint* joint, std::vector<std::string> filenames, dtype density, bool is_infinite_mass = false);

        void compute_mass_inertial();

    private:
        void load_mesh(std::string filename);

        void process_mesh();

        void compute_mass_property(const Matrix3X& V, const Matrix3Xi& F, /*input*/
            dtype& mass, Vector3& COM,          /*output*/
            Matrix3& I);

        void VolumeIntegration(const Matrix3X& V, const Matrix3Xi& F, /*input*/
            dtype& T0, Vector3& T1, Vector3& T2, Vector3& TP /*output*/);
    };

}