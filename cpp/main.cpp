#include "SimEnvGenerator.h"
#include "Model.h"
#include "Robot.h"
#include "Body/BodyCuboid.h"
#include "Body/BodyPlane.h"

#include "Utils.h"
#include "Common.h"
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

int main() {
	_2psp::Model* sim =_2psp::SimEnvGenerator::createGroundTest("TGS");
	auto t0 = Clock::now();
	sim->forward(100, true);
	auto t1 = Clock::now();
	std::cerr << "time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms" << std::endl;
	sim->export_replay("../Results/");
}
