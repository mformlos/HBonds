#include <sys/time.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <iostream>
#include <csignal>
#include <algorithm>
#include <string>
#include <vector>
#include <math.h>
#include <list>
#include <map>
#include <forward_list>
#include "../lib/Eigen/Eigen/Dense"
#include "../lib/Water.h"
#include "../lib/HelperFunctions.h"
using namespace Eigen;




int main(int argc, char* argv[]) {
	std::string ConfigFileStart{}, ConfigFileName{}, Directory{}, OutputDir{}, OutputFileName{}, OutputFileStart{};
	double CutOff{3.5}, MaxAngle{M_PI/6.}, AverageHbonds{0.0};
	unsigned NumberOfMolecules{}, StartStep {}, SamplingStep{}, Step{}, NSteps{0}, MaxRing{};
	std::array<unsigned,3> Cells;
	std::array<double,3> BoxSize, CellSideLength;
	std::vector<Water> Molecules;
	std::vector<std::vector<std::vector<std::forward_list<Water*>>>> CellList;
	std::map<unsigned, double> ring_dist;


	if (argc != 5) {
		std::cout << "usage: ./ring_search RUNPARENTDIR STARTSTEP SAMPLINGSTEP NMOLECULES MAXRING" << std::endl;
	}

	Directory=argv[1];
	StartStep = std::stoi(argv[2]);
	SamplingStep = std::stoi(argv[3]);
	NumberOfMolecules = std::stoi(argv[4]);
	MaxRing = std::stoi(argv[5]);

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

	for (unsigned i = 0; i <= MaxRing; i++) {
		ring_dist[i] = 0.0;
	}

	timeval start {}, end {};
	gettimeofday(&start, NULL);

	Step = StartStep;
	while (true) {
		ConfigFileName = ConfigFileStart+std::to_string(Step)+".pdb";
		if (!initializePositionsXYZ(Molecules, BoxSize, ConfigFileName)) {
			std::cout << "reached last step" << std::endl;
			break;
		}
		std::cout << "Step: " <<  Step << " , BoxSize: " << BoxSize[0] << " " << BoxSize[1] << " " << BoxSize[2] << std::endl;


		CellList = updateCellList(Molecules, BoxSize, CellSideLength, Cells, CutOff);
		//std::cout << "updated cell list" << std::endl;
		findNeighbours(Molecules, CellList, BoxSize, CellSideLength, Cells, CutOff, MaxAngle);
		/*for (unsigned i = 0; i < NumberOfMolecules; i++) {
			for (unsigned j = 0; j < NumberOfMolecules; j++) {
				if (HydrogenBond(Molecules[i], Molecules[j], BoxSize, CutOff, MaxAngle)) {
					Molecules[i].HBonds.push_front(&Molecules[j]);
					//other -> HBonds.push_front(&mol);
				}
			}
		}*/
		//std::cout << "found neighbours" << std::endl;

		unsigned count {0};
		for (auto& mol : Molecules) {
			for (auto& bond : mol.HBonds ) {
				count++;
			}
		}
		std::cout << "total H-bonds: " << count << std::endl;
		AverageHbonds += count;

		/// RING SEARCH
		for (unsigned i = 0; i < NumberOfMolecules; i++) {
			//std::cout << "start node: " << i << std::endl;
			searchRings(Molecules, i, NumberOfMolecules, MaxRing);
			//if (Molecules[i].MinRing < 4) std::cout << "start node: " << i << "Min Ring: " << Molecules[i].MinRing << std::endl;
		}
		for (unsigned i = 0; i < NumberOfMolecules; i++) {
			if (Molecules[i].MinRing == 100) Molecules[i].MinRing = 0;
			ring_dist[Molecules[i].MinRing]++;
			Molecules[i].MinRing = 100;
		}
		NSteps++;
		Step += SamplingStep;
	}
	std::ofstream output(Directory+"/ring_dist.dat");
	output << "# ring_size    P(ring_size)" << std::endl;
	for (auto& ring_size : ring_dist) {
		//output << ring_size.first << " " << ring_size.second << std::endl;
		output << ring_size.first << " " << ring_size.second/(NSteps*NumberOfMolecules) << std::endl;
	}
	output.close();
	std::cout << "Average Hbonds per molecule: " << AverageHbonds/(NumberOfMolecules*NSteps) << std::endl;
	gettimeofday(&end, NULL);
	double realTime { ((end.tv_sec - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6 };
	std::cout << "total time: " << realTime << " , time per step: " << realTime/NSteps<< std::endl;
}
