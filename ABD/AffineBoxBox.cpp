// To compile on my machine:

// mex COPTIMFLAGS='-O3 -DNDEBUG' -I"D:\ZiyanXiong\Cloth_Sim_Env\eigen-3.4.0" -I"C:\Program Files (x86)\libccd\include" -I"C:\Program Files (x86)\fcl\include" -L"C:\Program Files (x86)\libccd\lib" -L"C:\Program Files (x86)\fcl\lib" AffineBoxBox.cpp fcl.lib ccd.lib


#include <iostream>
#include "mex.h"
#include "fcl/fcl.h"

#define A(i,j) A[i+j*M]
using namespace std;
using namespace fcl;

//==============================================================================
template <typename S>
void test_Box_Box(Transform3<S>& tf0, Vector3<S>& whd0, Transform3<S>& tf1, Vector3<S>& whd1, vector<Contact<S>>& contacts)
{
    std::shared_ptr<Box<S>> box0(new Box<S>(whd0));
    std::shared_ptr<Box<S>> box1(new Box<S>(whd1));

    //  GJKSolver_indep solver;
    detail::GJKSolver_libccd<S> solver;

    static const int num_max_contacts = std::numeric_limits<int>::max();
    static const bool enable_contact = true;
    fcl::CollisionResult<S> result;
    fcl::CollisionRequest<S> request(num_max_contacts,
        enable_contact);

    CollisionObject<S> co0(box0, tf0);
    CollisionObject<S> co1(box1, tf1);

    fcl::collide(&co0, &co1, request, result);

    result.getContacts(contacts);
    for(int i = 0; i < contacts.size(); i++) {
        mexPrintf("Contact%d:\n", i);
        mexPrintf("pos:%lf %lf %lf\n", contacts[i].pos[0], contacts[i].pos[1], contacts[i].pos[2]);
        mexPrintf("normal:%lf %lf %lf\n", contacts[i].normal[0], contacts[i].normal[1], contacts[i].normal[2]);
    }
}

template <typename S>
void get_box(std::vector<Vector3<S>>& vertices, std::vector<Triangle>& triangles, Matrix3<S> A, Vector3<S> t)
{
    vertices.clear();
    triangles.clear();

    vertices.push_back(A * Vector3<S>(1.000000, 1.000000, -1.000000) + t);
    vertices.push_back(A * Vector3<S>(1.000000, - 1.000000, - 1.000000) + t);
    vertices.push_back(A * Vector3<S>(1.000000, 1.000000, 1.000000) + t);
    vertices.push_back(A * Vector3<S>(1.000000, - 1.000000, 1.000000) + t);
    vertices.push_back(A * Vector3<S>(-1.000000, 1.000000, - 1.000000) + t);
    vertices.push_back(A * Vector3<S>(-1.000000, - 1.000000, - 1.000000) + t);
    vertices.push_back(A * Vector3<S>(-1.000000, 1.000000, 1.000000) + t);
    vertices.push_back(A * Vector3<S>(-1.000000, - 1.000000, 1.000000) + t);

    triangles.push_back(Triangle(4, 2, 0));
    triangles.push_back(Triangle(2, 7, 3));
    triangles.push_back(Triangle(6, 5, 7));
    triangles.push_back(Triangle(1, 7, 5));
    triangles.push_back(Triangle(0, 3, 1));
    triangles.push_back(Triangle(4, 1, 5));
    triangles.push_back(Triangle(4, 6, 2));
    triangles.push_back(Triangle(2, 6, 7));
    triangles.push_back(Triangle(6, 4, 5));
    triangles.push_back(Triangle(1, 3, 7));
    triangles.push_back(Triangle(0, 2, 3));
    triangles.push_back(Triangle(4, 0, 1));
}

template<typename BV>
int mesh_mesh_collison(const Transform3<typename BV::S>& pose1, const Transform3<typename BV::S>& pose2,
    const std::vector<Vector3<typename BV::S>>& vertices1, const std::vector<Triangle>& triangles1,
    const std::vector<Vector3<typename BV::S>>& vertices2, const std::vector<Triangle>& triangles2, detail::SplitMethodType split_method, std::vector<Contact<typename BV::S>>& contacts)
{
    using S = typename BV::S;

    BVHModel<BV> m1;
    BVHModel<BV> m2;
    m1.bv_splitter.reset(new detail::BVSplitter<BV>(split_method));
    m2.bv_splitter.reset(new detail::BVSplitter<BV>(split_method));

    m1.beginModel();
    m1.addSubModel(vertices1, triangles1);
    m1.endModel();

    m2.beginModel();
    m2.addSubModel(vertices2, triangles2);
    m2.endModel();

    // Transform3<S> pose1(tf);
    // Transform3<S> pose2 = Transform3<S>::Identity();
    static const int num_max_contacts = std::numeric_limits<int>::max();
    static const bool enable_contact = true;

    //std::vector<Contact<S>> contacts;

    CollisionRequest<S> request(num_max_contacts, enable_contact);
    CollisionResult<S> result;
    int num_contacts = collide(&m1, pose1, &m2, pose2, request, result);

    result.getContacts(contacts);

    return contacts.size();
}

template <typename S>
void convex_convex_collison(const std::shared_ptr<const std::vector<Vector3<S>>>& vertices0, const std::shared_ptr<const std::vector<int>>& faces0,
    const std::shared_ptr<const std::vector<Vector3<S>>>& vertices1, const std::shared_ptr<const std::vector<int>>& faces1, std::vector<Contact<S>>& contacts)
{
    std::shared_ptr<Convex<S>> convex0 = std::make_shared<Convex<S>>(vertices0, 6, faces0);
    std::shared_ptr<Convex<S>> convex1 = std::make_shared<Convex<S>>(vertices1, 6, faces1);

    Transform3<S> tf0, tf1;
    tf0.setIdentity();
    //tf0.translation() = Vector3<S>(.9, 0, 0);
    //tf0.linear() = Quaternion<S>(.6, .8, 0, 0).toRotationMatrix();
    tf1.setIdentity();

    static const int num_max_contacts = std::numeric_limits<int>::max();
    static const bool enable_contact = true;
    fcl::CollisionResult<S> result;
    fcl::CollisionRequest<S> request(num_max_contacts,
        enable_contact);

    CollisionObject<S> co0(convex0, tf0);
    CollisionObject<S> co1(convex1, tf1);

    fcl::collide(&co0, &co1, request, result);
    result.getContacts(contacts);

}

template <typename S>
void get_convex_box(std::shared_ptr<std::vector<Vector3<S>>>& vertices, std::shared_ptr<std::vector<int>>& faces, Matrix3<S> A, Vector3<S> t)
{
    vertices->clear();
    faces->clear();

    vertices->push_back(A * Vector3<S>(1.000000, 1.000000, -1.000000) + t);
    vertices->push_back(A * Vector3<S>(1.000000, -1.000000, -1.000000) + t);
    vertices->push_back(A * Vector3<S>(1.000000, 1.000000, 1.000000) + t);
    vertices->push_back(A * Vector3<S>(1.000000, -1.000000, 1.000000) + t);
    vertices->push_back(A * Vector3<S>(-1.000000, 1.000000, -1.000000) + t);
    vertices->push_back(A * Vector3<S>(-1.000000, -1.000000, -1.000000) + t);
    vertices->push_back(A * Vector3<S>(-1.000000, 1.000000, 1.000000) + t);
    vertices->push_back(A * Vector3<S>(-1.000000, -1.000000, 1.000000) + t);

    faces->insert(faces->end(), { 4, 0, 4, 6, 2 });
    faces->insert(faces->end(), { 4, 3, 2, 6, 7 });
    faces->insert(faces->end(), { 4, 7, 6, 4, 5 });
    faces->insert(faces->end(), { 4, 5, 1, 3, 7 });
    faces->insert(faces->end(), { 4, 1, 0, 2, 3 });
    faces->insert(faces->end(), { 4, 5, 4, 0, 1 });
}

template <typename S>
void test_convex()
{
    // set mesh triangles and vertice indices
    std::shared_ptr<std::vector<Vector3<S>>> p1, p2;
    std::shared_ptr<std::vector<int>> t1, t2;

    p1 = std::make_shared<std::vector<Vector3<S>>>();
    p2 = std::make_shared<std::vector<Vector3<S>>>();
    t1 = std::make_shared<std::vector<int>>();
    t2 = std::make_shared<std::vector<int>>();

    Matrix3<S> A1, A2;
    Vector3<S> tl1, tl2;
    A1.setIdentity();
    A2.setIdentity();
    tl1 = Vector3<S>(0.0, 0.1, -1.5);
    tl2.setZero();

    get_convex_box<S>(p1, t1, A1, tl1);
    get_convex_box<S>(p2, t2, A2, tl2);

    vector<Contact<S>> contacts;

    convex_convex_collison<S>(p1, t1, p2, t2, contacts);

    cout << contacts.size() << " contacts found" << endl;
    for(int i = 0; i < contacts.size(); i++) {
        mexPrintf("Contact%d:\n", i);
        mexPrintf("pos:%lf %lf %lf\n", contacts[i].pos[0], contacts[i].pos[1], contacts[i].pos[2]);
        mexPrintf("normal:%lf %lf %lf\n", contacts[i].normal[0], contacts[i].normal[1], contacts[i].normal[2]);
    }
}

vector<Contact<double>> fclAffineBoxBox(Eigen::Matrix4d& E1, Eigen::Vector3d& whd1, Eigen::Matrix4d& E2, Eigen::Vector3d& whd2){
    // set mesh triangles and vertice indices
    std::shared_ptr<std::vector<Vector3<double>>> p1, p2;
    std::shared_ptr<std::vector<int>> t1, t2;

    p1 = std::make_shared<std::vector<Vector3<double>>>();
    p2 = std::make_shared<std::vector<Vector3<double>>>();
    t1 = std::make_shared<std::vector<int>>();
    t2 = std::make_shared<std::vector<int>>();

    get_convex_box<double>(p1, t1, 0.5 * E1.topLeftCorner(3, 3) * whd1.asDiagonal(), E1.block<3, 1>(0, 3));
    get_convex_box<double>(p2, t2, 0.5 * E2.topLeftCorner(3, 3) * whd2.asDiagonal(), E2.block<3, 1>(0, 3));

    vector<Contact<double>> contacts;
    convex_convex_collison<double>(p1, t1, p2, t2, contacts);

    return contacts;
}

vector<Contact<double>> fclAffineBoxBoxMesh(Eigen::Matrix4d& E1, Eigen::Vector3d& whd1, Eigen::Matrix4d& E2, Eigen::Vector3d& whd2){
    vector<Contact<double>> contacts;
    Transform3<double> tf1, tf2;

    tf1.setIdentity();
    tf2.setIdentity();

    // set mesh triangles and vertice indices
    std::vector<Vector3<double>> p1, p2;
    std::vector<Triangle> t1, t2;

    get_box<double>(p1, t1, 0.5 * E1.topLeftCorner(3, 3) * whd1.asDiagonal(), E1.block<3, 1>(0, 3));
    get_box<double>(p2, t2, 0.5 * E2.topLeftCorner(3, 3) * whd2.asDiagonal(), E2.block<3, 1>(0, 3));
    mesh_mesh_collison<OBBRSS<double>>(tf1, tf2, p1, t1, p2, t2, detail::SPLIT_METHOD_MEAN, contacts);

    return contacts;
}

vector<Contact<double>> fclBoxBox(Eigen::Matrix4d& E1, Eigen::Vector3d& whd1, Eigen::Matrix4d& E2, Eigen::Vector3d& whd2){
    vector<Contact<double>> contacts;
    Transform3<double> tf1, tf2;
    tf1.setIdentity();
    tf1.translation() = E1.block<3, 1>(0, 3);
    tf1.linear() = E1.topLeftCorner(3, 3);
    tf2.setIdentity();
    tf2.translation() = E2.block<3, 1>(0, 3);
    tf2.linear() = E2.topLeftCorner(3, 3);
    test_Box_Box<double>(tf1, whd1, tf2, whd2, contacts);
    return contacts;
}

//==============================================================================
void mexFunction( int nlhs, mxArray *plhs[],
				  int nrhs, const mxArray *prhs[])
{
	// check for proper number of arguments
	if(nrhs != 4) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","4 inputs required.");
	}
	if(nlhs != 1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","1 output required.");
	}

	enum
	{
		RHS_E1 = 0,
		RHS_WHD1,
		RHS_E2,
		RHS_WHD2,
		RHS_COUNT
	};

	if(!mxIsDouble(prhs[RHS_E1]) || mxGetNumberOfElements(prhs[RHS_E1]) != 16) {
		mexErrMsgTxt("E1 must be a mat4.");
	}
	if(!mxIsDouble(prhs[RHS_WHD1]) || mxGetNumberOfElements(prhs[RHS_WHD1]) != 3) {
		mexErrMsgTxt("whd1 must be a vec3.");
	}
	if(!mxIsDouble(prhs[RHS_E2]) || mxGetNumberOfElements(prhs[RHS_E2]) != 16) {
		mexErrMsgTxt("E1 must be a mat4.");
	}
	if(!mxIsDouble(prhs[RHS_WHD2]) || mxGetNumberOfElements(prhs[RHS_WHD2]) != 3) {
		mexErrMsgTxt("whd1 must be a vec3.");
	}

	mwSize M;
	double *A;

	// Convert E1
	M = mxGetM(prhs[RHS_E1]); A = mxGetPr(prhs[RHS_E1]);
	Eigen::Matrix4d E1;
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < 4; ++j) {
			E1(i,j) = A(i,j);
		}
	}

	// Convert whd1
	M = mxGetM(prhs[RHS_WHD1]); A = mxGetPr(prhs[RHS_WHD1]);
	Eigen::Vector3d whd1;
	for(int i = 0; i < 3; ++i) {
		whd1(i) = A(i,0);
	}

	// Convert E2
	M = mxGetM(prhs[RHS_E2]); A = mxGetPr(prhs[RHS_E2]);
	Eigen::Matrix4d E2;
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j < 4; ++j) {
			E2(i,j) = A(i,j);
		}
	}

	// Convert whd2
	M = mxGetM(prhs[RHS_WHD2]); A = mxGetPr(prhs[RHS_WHD2]);
	Eigen::Vector3d whd2;
	for(int i = 0; i < 3; ++i) {
		whd2(i) = A(i,0);
	}

	// Call collision detector
	//vector<Contact<double>> contacts = fclBoxBox(E1, whd1, E2, whd2);
    vector<Contact<double>> contacts = fclAffineBoxBox(E1, whd1, E2, whd2);

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
	mxSetField(c, 0, "count", mxCreateDoubleScalar(contacts.size()));
	mxSetField(c, 0, "depthMax", mxCreateDoubleScalar(0.0));
	mxArray *mat;
	mat = mxCreateDoubleMatrix(3, contacts.size(), mxREAL); M = mxGetM(mat); A = mxGetPr(mat);
	for(int k = 0; k < contacts.size(); ++k) {
	    for(int i = 0; i < 3; ++i) {
		    A(i,k) = contacts[k].normal(i);
	    }
	}
	mxSetField(c, 0, "nor", mat);
	mat = mxCreateDoubleMatrix(3, contacts.size(), mxREAL); M = mxGetM(mat); A = mxGetPr(mat);
	for(int k = 0; k < contacts.size(); ++k) {
		for(int i = 0; i < 3; ++i) {
			A(i,k) = contacts[k].pos(i);
		}
	}
	mxSetField(c, 0, "pos", mat);
	mat = mxCreateDoubleMatrix(1, contacts.size(), mxREAL); M = mxGetM(mat); A = mxGetPr(mat);
	for(int k = 0; k < contacts.size(); ++k) {
		A(0,k) = contacts[k].penetration_depth;
	}
	mxSetField(c, 0, "depth", mat);
	mxFree((void *)fnames);
	plhs[0] = c;
}

