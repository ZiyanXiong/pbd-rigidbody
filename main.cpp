#include "SimEnvGenerator.h"
#include "Model.h"
#include "Robot.h"
#include "Body/BodyCuboid.h"
#include "Body/BodyPlane.h"

#include "Utils.h"
#include "Common.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
typedef std::chrono::high_resolution_clock Clock;

long long test_speed(int model_id, std::string solver, int substeps) {
	_2psp::Model* sim =_2psp::SimEnvGenerator::createScene(model_id, solver, substeps);
	std::chrono::milliseconds duration(0);
	for (int i = 0; i < 10; i++) {
		if (model_id == 11) {
			_2psp::VectorX f_tm = _2psp::VectorX::Zero(sim->_ndof_m / 7 * 6);
			f_tm.segment<6>(45 * 6) << 0., -800000, 0., 0., 0., 0.;
			auto t0 = Clock::now();
			for (int i = 0; i < 300; i++) {
				if (i < 20)
					sim->_f_tm = f_tm;
				else
					sim->_f_tm.setZero();
				sim->forward(1, true);
			}
			auto t1 = Clock::now();
			std::cerr << "Running" << sim->_name << " " << solver << " " << substeps << ", time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms" << std::endl;
			duration += std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
		}
		else {
			auto t0 = Clock::now();
			sim->forward(300, true);
			auto t1 = Clock::now();
			std::cerr <<"Running" << sim->_name << " " << solver << " " << substeps << ", time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms" << std::endl;
			duration += std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
		}
		sim->export_replay("../Results/" + std::to_string(model_id));
		sim->reset();
	}
	std::cerr << "Average time for " << sim->_name << " = " << duration.count() / 10 << "ms" << std::endl;
	return duration.count() / 10;
}

void save_time_report(std::string folder, Eigen::VectorXi& model_ids, std::vector<std::vector<long long>>& durations) {
	if (folder.back() != '/') {
		folder += '/';
	}
	//create folder if it does not exist
	if (!std::filesystem::exists(folder)) {
		std::filesystem::create_directories(folder);
	}
	std::ofstream result_file;
	result_file.open(folder + "time_report.txt");

	if (result_file.is_open()) {
		result_file << "# Each line is formatted as: Model id, TGS50, TGS150, TGS500, 2PSP" << std::endl;
		for (size_t i = 0; i < durations.size(); i++) {
			result_file << model_ids(i) << ", " << durations[i][0] << ", " << durations[i][1] << ", " << durations[i][2] << ", " << durations[i][3] << std::endl;
		}
		result_file.close();
		std::cout << "Time results written to file successfully." << std::endl;
	}
	else {
		std::cerr << "Unable to open file." << std::endl;
	}
}

int main() {
	Eigen::VectorXi model_ids(10);
	model_ids << 0, 4, 5, 8, 9, 10, 11, 13, 14, 15;
	//Eigen::VectorXi model_ids(5);
	//model_ids << 10, 11, 13, 14, 15;
	std::vector<std::vector<long long>> durations(model_ids.size(), std::vector<long long>(4, 0));

	for (int i = 0; i < model_ids.size(); i++) {
		std::cerr << "Testing model " << model_ids(i) << std::endl;
		durations[i][0] = test_speed(model_ids(i), "TGS", 50);
		durations[i][1] = test_speed(model_ids(i), "TGS", 150);
		durations[i][2] = test_speed(model_ids(i), "TGS", 500);
		durations[i][3] = test_speed(model_ids(i), "2PSP", 150);
	}

	save_time_report("../Results/", model_ids, durations);

	//int model_id = 9;
	//_2psp::Model* sim =_2psp::SimEnvGenerator::createScene(model_id, "TGS", 500);
	//if (model_id == 11) {
	//	_2psp::VectorX f_tm = _2psp::VectorX::Zero(sim->_ndof_m / 7 * 6);
	//	f_tm.segment<6>(45 * 6) << 0., -800000, 0., 0., 0., 0.;
	//	auto t0 = Clock::now();
	//	for (int i = 0; i < 300; i++) {
	//		if (i < 20)
	//			sim->_f_tm = f_tm;
	//		else
	//			sim->_f_tm.setZero();
	//		sim->forward(1, true);
	//	}
	//	auto t1 = Clock::now();
	//	std::cerr << "time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms" << std::endl;
	//}
	//else {
	//	auto t0 = Clock::now();
	//	sim->forward(300, true);
	//	auto t1 = Clock::now();
	//	std::cerr << "time = " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << "ms" << std::endl;
	//}
	//sim->export_replay("../Results/");
}
