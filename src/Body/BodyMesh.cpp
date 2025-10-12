#include "Body/BodyMesh.h"
#include "Utils.h"
#include "Body/tiny_obj_loader.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

namespace _2psp {

    std::shared_ptr<coal::ConvexBase> loadConvexMesh(const std::string& file_name) {
        coal::NODE_TYPE bv_type = coal::BV_AABB;
        coal::MeshLoader loader(bv_type);
        coal::BVHModelPtr_t bvh = loader.load(file_name);
        bvh->buildConvexHull(true, "Qt");
        return bvh->convex;
    }

    BodyMesh::BodyMesh(Model* sim, Joint* joint, std::vector<std::string> filenames, dtype density, bool is_infinite_mass) :
        Body(sim, joint, density)
    {
        _filename = filenames[0];
        for (size_t i = 0; i < filenames.size(); ++i) {
			_convex_shapes.push_back(loadConvexMesh(filenames[i]));
		}
        _file_names = filenames;
        _is_infinite_mass = is_infinite_mass;

    }

    void BodyMesh::compute_mass_inertial()
	{
        if (_is_infinite_mass) {
			_mass_inv = 0.0;
			_Inertia_inv.setZero();
			return;
		}

        load_mesh(_filename);
        process_mesh();
	}

    /*
    void BodyMesh::load_mesh(std::string filename) {
        std::vector<tinyobj::shape_t> obj_shape;
        std::vector<tinyobj::material_t> obj_material;
        tinyobj::attrib_t attrib;
        std::string err;
        tinyobj::LoadObj(&attrib, &obj_shape, &obj_material, &err, filename.c_str());

        int num_vertices = (int)attrib.vertices.size() / 3;
        _V.resize(3, num_vertices);
        for (int i = 0; i < num_vertices; i++) {
            _V.col(i) = Vector3(attrib.vertices[i * 3],
                attrib.vertices[i * 3 + 1],
                attrib.vertices[i * 3 + 2]);
        }

        int num_elements = (int)obj_shape[0].mesh.indices.size() / 3;
        _F.resize(3, num_elements);
        for (int i = 0; i < num_elements; i++) {
            _F.col(i) = Vector3i(obj_shape[0].mesh.indices[i * 3].vertex_index,
                obj_shape[0].mesh.indices[i * 3 + 1].vertex_index,
                obj_shape[0].mesh.indices[i * 3 + 2].vertex_index);
        }
    }
    */

    void BodyMesh::load_mesh(std::string filename) {
        std::vector<Eigen::Vector3d> vertices;
        std::vector<Eigen::Vector3i> faces;

        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
        }

        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string line_type;
            ss >> line_type;

            if (line_type == "v") {
                double x, y, z;
                ss >> x >> y >> z;
                vertices.push_back(Eigen::Vector3d(x, y, z));
            }
            else if (line_type == "f") {
                Eigen::Vector3i face;
                int v_idx;
                char slash;

                for (int i = 0; i < 3; ++i) {
                    ss >> v_idx;
                    // .obj files are 1-indexed, so we subtract 1 for 0-based indexing
                    face(i) = v_idx - 1;

                    // Handle different face formats (v/vt, v//vn, v/vt/vn)
                    if (ss.peek() == '/') {
                        ss.get(slash); // Consume the '/'
                        if (ss.peek() == '/') {
                            ss.get(slash); // Consume the second '/' for v//vn
                            ss >> v_idx; // Unused normal index
                        }
                        else {
                            ss >> v_idx; // Unused texture index
                            if (ss.peek() == '/') {
                                ss.get(slash); // Consume the '/' for v/vt/vn
                                ss >> v_idx; // Unused normal index
                            }
                        }
                    }
                }
                faces.push_back(face);
            }
        }

        // Convert std::vector to Eigen::Matrix
        _V.resize(3, vertices.size());
        for (size_t i = 0; i < vertices.size(); ++i) {
            _V.col(i) = vertices[i];
        }

        _F.resize(3, faces.size());
        for (size_t i = 0; i < faces.size(); ++i) {
            _F.col(i) = faces[i];
        }
    }


    void BodyMesh::process_mesh() {
        // compute mass properties
        dtype T0; Vector3 T1; Vector3 T2; Vector3 TP;
        Vector3 r;
        Matrix3 J;
        VolumeIntegration(_V, _F, T0, T1, T2, TP);

        //std::cout << "T0: " << T0 << std::endl;
        //std::cout << "T1: " << T1 << std::endl;
        //std::cout << "T2: " << T2 << std::endl;
        //std::cout << "TP: " << TP << std::endl;

        // compute mass
        dtype mass = T0 * _density;
        _mass_inv = 1.0 / mass;
        r = T1 / T0;

        J(0, 0) = _density * (T2(1) + T2(2)); // Ixx
        J(1, 1) = _density * (T2(2) + T2(0)); // Iyy
        J(2, 2) = _density * (T2(0) + T2(1)); // Izz

        J(0, 1) = -_density * TP(0); // Ixy
        J(1, 2) = -_density * TP(1); // Iyz
        J(2, 0) = -_density * TP(2); // Izx

        // Symmetrize the matrix
        J(1, 0) = J(0, 1);
        J(2, 1) = J(1, 2);
        J(0, 2) = J(2, 0);

        // --- Translate inertia tensor to center of mass using Parallel Axis Theorem ---
        // J_cm = J_origin - m * ([r]^T[r]I - r*r^T)
        // where [r] is the skew-symmetric matrix of r.
        J(0, 0) -= mass * (r(1) * r(1) + r(2) * r(2));
        J(1, 1) -= mass * (r(2) * r(2) + r(0) * r(0));
        J(2, 2) -= mass * (r(0) * r(0) + r(1) * r(1));

        J(0, 1) += mass * r(0) * r(1);
        J(1, 2) += mass * r(1) * r(2);
        J(2, 0) += mass * r(2) * r(0);

        // Symmetrize the translated matrix
        J(1, 0) = J(0, 1);
        J(2, 1) = J(1, 2);
        J(0, 2) = J(2, 0);

        //std::cout << "J: " << J << std::endl;
        // get the principal axes for inertia tensor by eigenvalue decomposition
        // https://en.wikipedia.org/wiki/Moment_of_inertia#Principal_axes
        _Inertia_inv.setZero();
        Eigen::SelfAdjointEigenSolver<Matrix3> eigensolver(J);
        Vector3 eig_values = eigensolver.eigenvalues();
        Matrix3 eig_vectors = eigensolver.eigenvectors();
        _Inertia_inv.head(3) = eig_values.cwiseInverse();
        _Inertia_inv(3) = _Inertia_inv(4) = _Inertia_inv(5) = 1 / mass;
        Matrix4 E = Matrix4::Identity();
        E.topLeftCorner(3, 3) = -eig_vectors;
        E.topRightCorner(3, 1) = r;

        // check for right-handedness
        Vector3 x = E.block(0, 0, 3, 1);
        Vector3 y = E.block(0, 1, 3, 1);
        Vector3 z = E.block(0, 2, 3, 1);
        if (x.cross(y).dot(z) < 0.0)
            E.block(0, 2, 3, 1) *= -1;

        // check if the rotation part is valid
        Matrix3 res = E.topLeftCorner(3, 3) * E.topLeftCorner(3, 3).transpose();
        if ((res - Matrix3::Identity()).norm() > 1e-6) {
            std::cerr << "invalid rotational part: " << std::endl << E.topLeftCorner(3, 3) << std::endl;
        }

        // transform the mesh into body frame
        _E_io = math::Einv(E);

        //std::cout << "E: " << E << std::endl;
        //std::cout << "VOLUME: " << T0 << std::endl;
        //std::cout << "COM: " << r.transpose() << std::endl;


        Matrix3 R_io = _E_io.topLeftCorner(3, 3);
        Vector3 p_io = _E_io.topRightCorner(3, 1);

        _V = (R_io * _V).colwise() + p_io;

        // std::cerr << "body " << _name << ", I = " << _Inertia.transpose() << std::endl;
    }
    
    void BodyMesh::compute_mass_property(
        const Matrix3X& V, const Matrix3Xi& F, /*input*/
        dtype& volume, Vector3& COM,          /*output*/
        Matrix3& I) {

        int num_vertices = V.cols();
        int num_faces = F.cols();

        // compute COM and volume
        volume = 0.;
        COM = Vector3::Zero();
        for (int i = 0; i < num_faces; i++) {
            Matrix3 A;
            A.col(0) = V.col(F(0, i));
            A.col(1) = V.col(F(1, i));
            A.col(2) = V.col(F(2, i));
            dtype vol = A.determinant();

            volume += vol;
            COM += vol * (A.col(0) + A.col(1) + A.col(2));
        }

        COM /= volume * 4.;
        volume /= 6.;

        // compute inertia tensor
        // assume mass = 1.
        Vector3 diag = Vector3::Zero();
        Vector3 offd = Vector3::Zero();
        for (int i = 0; i < num_faces; i++) {
            Matrix3 A;
            A.col(0) = V.col(F(0, i)) - COM;
            A.col(1) = V.col(F(1, i)) - COM;
            A.col(2) = V.col(F(2, i)) - COM;
            A.transposeInPlace();
            dtype d = A.determinant();

            for (int j = 0; j < 3; j++) {
                int j1 = (j + 1) % 3;
                int j2 = (j + 2) % 3;
                diag[j] += (A(0, j) * A(1, j) + A(1, j) * A(2, j) + A(2, j) * A(0, j) +
                    A(0, j) * A(0, j) + A(1, j) * A(1, j) + A(2, j) * A(2, j)) * d; // divide by 60.0f later;
                offd[j] += (A(0, j1) * A(1, j2) + A(1, j1) * A(2, j2) + A(2, j1) * A(0, j2) +
                    A(0, j1) * A(2, j2) + A(1, j1) * A(0, j2) + A(2, j1) * A(1, j2) +
                    A(0, j1) * A(0, j2) * 2 + A(1, j1) * A(1, j2) * 2 + A(2, j1) * A(2, j2) * 2) * d; // divide by 120.0f later
            }
        }

        diag /= volume * 60.;
        offd /= volume * 120.;
        I = (Matrix3() << diag(1) + diag(2), -offd(2), -offd(1),
            -offd(2), diag(0) + diag(2), -offd(0),
            -offd(1), -offd(0), diag(0) + diag(1)).finished();
    }

    /**
     * @brief Computes volume integrals of a triangulated mesh.
     *
     * This function translates a MATLAB script for computing volume integrals of a mesh
     * represented by vertices and faces into C++ using the Eigen library. These integrals
     * can be used to compute properties like volume, centroid, and moments of inertia.
     *
     * @param V An n-by-3 matrix of vertex coordinates, where n is the number of vertices.
     * @param F An m-by-3 matrix of face-vertex indices, where m is the number of faces.
     * Indices are assumed to be 0-based.
     * @param T0 A reference to a double to store the computed volume.
     * @param T1 A reference to an Eigen::Vector3d to store the first-order moments.
     * @param T2 A reference to an Eigen::Vector3d to store the second-order moments.
     * @param TP A reference to an Eigen::Vector3d to store the product of inertia terms.
     */
    void BodyMesh::VolumeIntegration(
        const Matrix3X& V,
        const Matrix3Xi& F,
        dtype& T0,
        Vector3& T1,
        Vector3& T2,
        Vector3& TP)
    {
        T0 = 0.0;
        T1.setZero();
        T2.setZero();
        TP.setZero();

        for (int i = 0; i < F.cols(); ++i) {
            // Get the vertices of the triangle
            Vector3 v0 = V.col(F(0, i));
            Vector3 v1 = V.col(F(1, i));
            Vector3 v2 = V.col(F(2, i));

            // Compute face normal
            Vector3 d10 = v1 - v0;
            Vector3 d20 = v2 - v0;
            Vector3 normal = d10.cross(d20);

            if (normal.norm() < 1e-9) {
                std::cerr << "Skipping bad triangle " << i << std::endl;
                continue;
            }
            Vector3 Normal = normal.normalized();

            int C;
            dtype nx = std::abs(Normal(0));
            dtype ny = std::abs(Normal(1));
            dtype nz = std::abs(Normal(2));

            if (nx > ny && nx > nz) {
                C = 0; // Project onto yz-plane
            }
            else if (ny > nz) {
                C = 1; // Project onto xz-plane
            }
            else {
                C = 2; // Project onto xy-plane
            }

            int A = (C + 1) % 3;
            int B = (A + 1) % 3;

            // Calculate the offset 'w'
            dtype w = -Normal.dot(v0);

            // Compute projection integrals
            dtype Pa = 0, Pb = 0, P1 = 0, Paa = 0, Pab = 0, Pbb = 0;
            dtype Paaa = 0, Paab = 0, Pabb = 0, Pbbb = 0;

            for (int j = 0; j < 3; ++j) {
                dtype a0 = V(A, F(j, i));
                dtype b0 = V(B, F(j, i));
                dtype a1 = V(A, F((j + 1) % 3, i));
                dtype b1 = V(B, F((j + 1) % 3, i));

                dtype da = a1 - a0;
                dtype db = b1 - b0;

                dtype a0_2 = a0 * a0, a0_3 = a0_2 * a0, a0_4 = a0_3 * a0;
                dtype b0_2 = b0 * b0, b0_3 = b0_2 * b0, b0_4 = b0_3 * b0;
                dtype a1_2 = a1 * a1, a1_3 = a1_2 * a1;
                dtype b1_2 = b1 * b1, b1_3 = b1_2 * b1;

                dtype C1 = a1 + a0;
                dtype Ca = a1 * C1 + a0_2;
                dtype Caa = a1 * Ca + a0_3;
                dtype Caaa = a1 * Caa + a0_4;
                dtype Cb = b1 * (b1 + b0) + b0_2;
                dtype Cbb = b1 * Cb + b0_3;
                dtype Cbbb = b1 * Cbb + b0_4;
                dtype Cab = 3 * a1_2 + 2 * a1 * a0 + a0_2;
                dtype Kab = a1_2 + 2 * a1 * a0 + 3 * a0_2;
                dtype Caab = a0 * Cab + 4 * a1_3;
                dtype Kaab = a1 * Kab + 4 * a0_3;
                dtype Cabb = 4 * b1_3 + 3 * b1_2 * b0 + 2 * b1 * b0_2 + b0_3;
                dtype Kabb = b1_3 + 2 * b1_2 * b0 + 3 * b1 * b0_2 + 4 * b0_3;

                P1 += db * C1;
                Pa += db * Ca;
                Paa += db * Caa;
                Paaa += db * Caaa;
                Pb += da * Cb;
                Pbb += da * Cbb;
                Pbbb += da * Cbbb;
                Pab += db * (b1 * Cab + b0 * Kab);
                Paab += db * (b1 * Caab + b0 * Kaab);
                Pabb += da * (a1 * Cabb + a0 * Kabb);
            }

            P1 /= 2.0;
            Pa /= 6.0;
            Paa /= 12.0;
            Paaa /= 20.0;
            Pb /= -6.0;
            Pbb /= -12.0;
            Pbbb /= -20.0;
            Pab /= 24.0;
            Paab /= 60.0;
            Pabb /= -60.0;

            // Compute face integrals
            dtype k1 = 1.0 / Normal(C);
            dtype k2 = k1 * k1;
            dtype k3 = k2 * k1;
            dtype k4 = k3 * k1;

            dtype Fa = k1 * Pa;
            dtype Fb = k1 * Pb;
            dtype Fc = -k2 * (Normal(A) * Pa + Normal(B) * Pb + w * P1);

            dtype Faa = k1 * Paa;
            dtype Fbb = k1 * Pbb;
            dtype Fcc = k3 * (pow(Normal(A), 2) * Paa + 2 * Normal(A) * Normal(B) * Pab + pow(Normal(B), 2) * Pbb +
                w * (2 * (Normal(A) * Pa + Normal(B) * Pb) + w * P1));

            dtype Faaa = k1 * Paaa;
            dtype Fbbb = k1 * Pbbb;
            dtype Fccc = -k4 * (pow(Normal(A), 3) * Paaa + 3 * pow(Normal(A), 2) * Normal(B) * Paab +
                3 * Normal(A) * pow(Normal(B), 2) * Pabb + pow(Normal(B), 3) * Pbbb +
                3 * w * (pow(Normal(A), 2) * Paa + 2 * Normal(A) * Normal(B) * Pab + pow(Normal(B), 2) * Pbb) +
                w * w * (3 * (Normal(A) * Pa + Normal(B) * Pb) + w * P1));

            dtype Faab = k1 * Paab;
            dtype Fbbc = -k2 * (Normal(A) * Pabb + Normal(B) * Pbbb + w * Pbb);
            dtype Fcca = k3 * (pow(Normal(A), 2) * Paaa + 2 * Normal(A) * Normal(B) * Paab + pow(Normal(B), 2) * Pabb +
                w * (2 * (Normal(A) * Paa + Normal(B) * Pab) + w * Pa));

            dtype Part;
            if (A == 0) Part = Fa; // A is x
            else if (B == 0) Part = Fb; // B is x
            else Part = Fc; // C is x

            T0 += Normal(0) * Part;

            Vector3 F1_vec, F2_vec, Fp_vec;
            F1_vec(A) = Faa; F1_vec(B) = Fbb; F1_vec(C) = Fcc;
            F2_vec(A) = Faaa; F2_vec(B) = Fbbb; F2_vec(C) = Fccc;
            Fp_vec(A) = Faab; Fp_vec(B) = Fbbc; Fp_vec(C) = Fcca;

            for (int j = 0; j < 3; ++j) {
                T1(j) += Normal(j) * F1_vec(j);
                T2(j) += Normal(j) * F2_vec(j);
                TP(j) += Normal(j) * Fp_vec(j);
            }
        }

        T1 /= 2.0;
        T2 /= 3.0;
        TP /= 2.0;
    }
}
