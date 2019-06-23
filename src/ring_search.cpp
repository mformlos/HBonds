#include <sys/time.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <iostream>
#include <csignal>
#include <string>
#include <vector>
#include <math.h>
#include <list>
#include "../lib/Eigen/Eigen/Dense"
#include "../lib/Water.h"
#include "../lib/HelperFunctions.h"
using namespace Eigen;




int main(int argc, char* argv[]) {
	std::string ConfigFileStart{}, ConfigFileName{}, Directory{}, OutputDir{}, OutputFileName{}, OutputFileStart{};
	double CutOff{3.5}, MaxAngle{M_PI/6.};
	unsigned NumberOfMolecules{}, StartStep {}, SamplingStep{}, Step{};
	std::array<unsigned,3> Cells;
	std::array<double,3> BoxSize, CellSideLength;
	std::vector<Water> Molecules;
	std::vector<std::vector<std::vector<std::forward_list<Water*>>>> CellList;

	if (argc != 4) {
		std::cout << "usage: ./ring_search RUNPARENTDIR STARTSTEP SAMPLINGSTEP NMOLECULES" << std::endl;
	}

	Directory=argv[1];
	StartStep = std::stoi(argv[2]);
	SamplingStep = std::stoi(argv[3]);
	NumberOfMolecules = std::stoi(argv[4]);

	OutputDir = Directory+"/Hbonds";
	int dir_err {mkdir(OutputDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)};
	/*if (-1 == dir_err) {
		printf("Error creating directory!");
		return EXIT_FAILURE;
	}*/
	OutputFileStart=OutputDir+"/Hbonds";
	ConfigFileStart = Directory+"/conf";

	for (unsigned i = 0; i < NumberOfMolecules; i++ ) {
		Molecules.push_back(i);
	}
	timeval start {}, end {};
	gettimeofday(&start, NULL);

	Step = StartStep;
	while (true) {
		ConfigFileName = ConfigFileStart+std::to_string(Step)+".pdb";
		if (!initializePositions(Molecules, BoxSize, ConfigFileName)) {
			std::cout << "reached last step" << std::endl;
			break;
		}
		std::cout << "BoxSize: " << BoxSize[0] << " " << BoxSize[1] << " " << BoxSize[2] << std::endl;
		std::cout << Molecules[0].PositionO.transpose() << std::endl;
		std::cout << Molecules[0].PositionH1.transpose() << std::endl;
		std::cout << Molecules[0].PositionH2.transpose() << std::endl;

		CellList = updateCellList(Molecules, BoxSize, CellSideLength, Cells, CutOff);
		std::cout << "updated cell list" << std::endl;
		findNeighbours(Molecules, CellList, BoxSize, CellSideLength, Cells, CutOff, MaxAngle);
		/*for (unsigned i = 0; i < NumberOfMolecules; i++) {
			for (unsigned j = 0; j < NumberOfMolecules; j++) {
				if (HydrogenBond(Molecules[i], Molecules[j], BoxSize, CutOff, MaxAngle)) {
					Molecules[i].HBonds.push_front(&Molecules[j]);
					//other -> HBonds.push_front(&mol);
				}
			}
		}*/
		std::cout << "found neighbours" << std::endl;
		for (auto& mol : Molecules[0].HBonds ) {
			std::cout << mol -> Index << std::endl;

		}
		unsigned count {0};
		for (auto& mol : Molecules) {
			for (auto& bond : mol.HBonds ) {
				count++;
			}
		}
		std::cout << "total H-bonds: " << count << std::endl;

		/// RING SEARCH

		Step += SamplingStep;
	}
}
