#include <fstream>
#include <istream>
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

using namespace Eigen;

inline Vector3d image(const Vector3d& pos, const std::array<double, 3>& BoxSize) {
    Vector3d newpos {pos};
    newpos(0) -= BoxSize[0]*floor(pos(0)/BoxSize[0]);
    newpos(1) -= BoxSize[1]*floor(pos(1)/BoxSize[1]);
    newpos(2) -= BoxSize[2]*floor(pos(2)/BoxSize[2]);
    return newpos;
}

inline Vector3d relPos(const Vector3d& posfirst, const Vector3d& possecond, const std::array<double, 3>& BoxSize) {
    Vector3d dist {possecond - posfirst};
    dist(0) -= BoxSize[0]*round(dist(0)/BoxSize[0]);
    dist(1) -= BoxSize[1]*round(dist(1)/BoxSize[1]);
    dist(2) -= BoxSize[2]*round(dist(2)/BoxSize[2]);
    return dist;
}

bool initializePositionsXYZ(std::vector<Water>& vec, std::array<double,3>& BoxSize, std::string filename) {
	std::ifstream file {filename};
	std::string line{};
	if (!file.is_open()) return false;
	std::string dump, Type;
	double x, y, z;
	unsigned count {0}, Natom {}, Ncurr{0}, NH1{0}, NH2{0};
	bool H1 {true};
	file >> std::skipws >> Natom;
	file >> std::skipws >> BoxSize[0] >> BoxSize[1] >> BoxSize[2];
	std::cout << Natom << " " << BoxSize[0] << " " << BoxSize[1] << " " << BoxSize[2] << std::endl;
	while (file >> std::skipws >> Type >> x >> y >> z) {
		if (Type == "O") {
    		vec[count].PositionO(0) = x;
    		vec[count].PositionO(1) = y;
    		vec[count].PositionO(2) = z;
		}
		else if (Type == "H" && H1) {
			vec[count].PositionH1(0) = x;
			vec[count].PositionH1(1) = y;
			vec[count].PositionH1(2) = z;
			H1=false;
			NH1++;
		}
		else if (Type == "H" && !H1) {
			vec[count].PositionH2(0) = x;
			vec[count].PositionH2(1) = y;
			vec[count].PositionH2(2) = z;
			H1=true;
			NH2++;
			count++;
		}
		Ncurr++;
	}
	if (count != vec.size()) {
		std::cout << "only " << count << " water molecules were initialized" << std::endl;
		return false;
	}
	if (Natom != Ncurr || NH2 != vec.size() || NH1 != vec.size()) {
		std::cout << "only " << Ncurr << " atoms, " << NH1 << "H1 atoms and " << NH2 << "H2 atoms were initialized" << std::endl;
		return false;
	}
	file.close();
	return true;
}

bool initializePositions(std::vector<Water>& vec, std::array<double,3>& BoxSize, std::string filename) {
    std::ifstream file {filename};
    std::string line{};
    if (!file.is_open()) return false;
    std::string dump, Type;
    double x, y, z;
    unsigned count {0}, Natom {}, Ncurr{};
    if (file.is_open()) {
    	getline(file, line);
    	getline(file, line);
    	getline(file, line);
    	file >> dump >> BoxSize[0] >> BoxSize[1] >> BoxSize[2] >> dump >> dump >> dump >> dump >> dump >> dump;
    	getline(file, line);
        while(getline(file, line)) {
            std::istringstream iss(line);
            if (iss  >> dump >> dump >> Type >> dump >> Ncurr >> x >> y >> z >> dump >> dump >> dump)
            {
            	if (Type == "OW") {
            		vec[count].PositionO(0) = x;
            		vec[count].PositionO(1) = y;
            		vec[count].PositionO(2) = z;
            	}
            	else if (Type == "HW1") {
            		vec[count].PositionH1(0) = x;
            		vec[count].PositionH1(1) = y;
            		vec[count].PositionH1(2) = z;
            	}
            	else if (Type == "HW2") {
					vec[count].PositionH2(0) = x;
					vec[count].PositionH2(1) = y;
					vec[count].PositionH2(2) = z;
					count++;
            	}
			}
            /*else if (iss >> x >> y >> z) {
            	BoxSize[0] = x;
            	BoxSize[1] = y;
            	BoxSize[2] = z;
            }*/
		}
        if (count != vec.size()) {
        	std::cout << "only " << count << " water molecules were initialized" << std::endl;
            return false;
        }
    }
    file.close();
    return true;
}

bool HydrogenBond (Water& First, Water& Second, std::array<double,3> BoxSize, double rCut, double angleMax) {
	bool bonded {false};
	Vector3d rOO {Vector3d::Zero()};
	Vector3d rOH {Vector3d::Zero()};
	double distanceOO {}, distanceOH{}, angle{};
	rOO = relPos(First.PositionO, Second.PositionO, BoxSize);
	distanceOO = rOO.norm();
	if (distanceOO > rCut) return false;
	/// first hydrogen
	rOH = relPos(First.PositionO, First.PositionH1, BoxSize);
	distanceOH = rOH.norm();
	angle = acos(rOH.dot(rOO)/(distanceOH*distanceOO));

	if (angle > 0.0 && angle < angleMax) return true;
	/// second hydrogen
	rOH = relPos(First.PositionO, First.PositionH2, BoxSize);
	distanceOH = rOH.norm();
	angle = acos(rOH.dot(rOO)/(distanceOH*distanceOO));
	if (angle > 0.0 && angle < angleMax) return true;
	return false;
}


std::vector<std::vector<std::vector<std::forward_list<Water*>>>> updateCellList(std::vector<Water>& Molecules, std::array<double,3>& BoxSize,  std::array<double,3>& CellSideLength, std::array<unsigned, 3>& Cells, double CutOff) {
	std::array<int, 3> CellNumber{};
	Cells[0] = unsigned(BoxSize[0]/CutOff);
	Cells[1] = unsigned(BoxSize[1]/CutOff);
	Cells[2] = unsigned(BoxSize[2]/CutOff);
	CellSideLength[0] = BoxSize[0]/(double)Cells[0];
	CellSideLength[1] = BoxSize[1]/(double)Cells[1];
	CellSideLength[2] = BoxSize[2]/(double)Cells[2];
	std::vector<std::vector<std::vector<std::forward_list<Water*>>>> CellList {std::vector<std::vector<std::vector<std::forward_list<Water*>>>>(Cells[0], std::vector<std::vector<std::forward_list<Water*>>>(Cells[1], std::vector<std::forward_list<Water*>>(Cells[2], std::forward_list<Water*>())))};
	for (auto& sheet : CellList) {
		for (auto& row : sheet) {
			for (auto& list : row) {
				list.clear();
			}
		}
	}
	for (auto& mol : Molecules) {
		Vector3d imPos;
		imPos = image(mol.PositionO, BoxSize);
		mol.HBonds.clear();
		for (unsigned i = 0; i < 3; i++) {
			CellNumber[i] = (int)(imPos(i)/CellSideLength[i]);
		}
		try {
			CellList[CellNumber[0]][CellNumber[1]][CellNumber[2]].push_front(&mol);
		}
		catch (std::exception& e) {
			std::cout << "Standard exception: " << e.what() << std::endl;
			std::cout << "Molecule: " << mol.Index << std::endl;
		}
	}
	return CellList;
}

bool findNeighbours(std::vector<Water>& Molecules, std::vector<std::vector<std::vector<std::forward_list<Water*>>>>& CellList, std::array<double,3>& BoxSize,  std::array<double,3>& CellSideLength, std::array<unsigned, 3>& Cells, double CutOff, double angleMax) {
	std::array<int, 3> CellNumber{};
	Vector3d imPos {};
    unsigned l, m, n;
	for (auto& mol : Molecules) {
		imPos = image(mol.PositionO, BoxSize);
		for (unsigned i = 0; i < 3; i++) {
			CellNumber[i] = (int)(imPos(i)/CellSideLength[i]);
		}
		for (int i = CellNumber[0]-2; i < CellNumber[0]+3; i++) {

			for (int j = CellNumber[1]-2; j < CellNumber[1]+3; j++) {

				for (int k = CellNumber[2]-2; k < CellNumber[2]+3; k++) {
					l = i - floor((double)i/Cells[0])*Cells[0];
					m = j - floor((double)j/Cells[1])*Cells[1];
					n = k - floor((double)k/Cells[2])*Cells[2];
					for (auto& other : CellList[l][m][n]) {
						if (other == &mol) continue;
						if (HydrogenBond(mol, *other, BoxSize, CutOff, angleMax)) {
							mol.HBonds.push_front(other -> Index);
							other -> HBonds.push_front(mol.Index);
						}
					}
				}
			}
		}
	}
}


/*void DFSUtil(unsigned u, unsigned v, std::vector<bool>& SptSet, std::vector<unsigned>& Distances)
{
    // Mark the current node as visited and
    // print it
    SptSet[v] = true;

    cout << v << " ";

    // Recur for all the vertices adjacent
    // to this vertex
    std::list<int>::iterator i;
    for (i = Molecules[v].HBonds.begin(); i !=  Molecules[v].HBonds.end(); ++i) {
        if (!SptSet[i->Index] ) DFSUtil(i->Index, SptSet);
    }
}*/

void loop (unsigned source, unsigned v, std::vector<Water>& Molecules, std::vector<unsigned>& Path, std::vector<std::vector<unsigned>>& AllPaths, unsigned MaxN) {
	//std::cout << Path.size() << " " << v << std::endl;

	for (auto& j : Molecules[v].HBonds) {
		std::vector<unsigned>::iterator PathLoc {std::find(Path.begin(), Path.end(), v)};
		if (PathLoc!=Path.end()) Path.erase(PathLoc+1, Path.end());
		bool stop {false};
		bool found_ring {false};
		//if (SptSet[j]) stop=true;
		if (j < source ) {
			//std::cout << "1" << std::endl;
			stop=true;
		}
		else if (Path.size()>MaxN) {
			//std::cout << "2" << std::endl;
			stop = true;
		}
		else if (std::find(Path.begin(), Path.end(), j) != Path.end()) {
			//std::cout << "3" << std::endl;
			stop=true;
		}
		else if (j == source) {
			if (Path.size()>1) {
				found_ring = true;
				Path.push_back(j);
			}
			stop=true;
			//std::cout << "4" << std::endl;
		}
		for (auto& bond : Molecules[j].HBonds) {
			if (std::find(Path.begin(), Path.end()-1, bond) != Path.end()-1) {
				//std::cout << "5" << std::endl;
				stop=true;
			}
		}
		if (found_ring) {
			AllPaths.push_back(Path);
			//std::cout << "ring found : ";
			for (auto& member : Path ) {
				//std::cout << member << " ";
				if (Path.size() < Molecules[member].MinRing) Molecules[member].MinRing = Path.size();
			}
			//std::cout << std::endl;
			//return true;
		}
		if (!stop) {
			Path.push_back(j);
			loop(source, j, Molecules, Path, AllPaths, MaxN);
		}
		//else return false;
	}
}

unsigned searchRings(std::vector<Water>& Molecules, unsigned Source, unsigned Nmol, unsigned MaxN) {
	//std::vector<unsigned> Distances(Nmol);
	//std::vector<bool> SptSet(Nmol);
	std::vector<unsigned> Path;
	std::vector<std::vector<unsigned>> AllPaths;
	//std::cout << "nearest neighbours: ";
	//for (auto& neighbour : Molecules[Source].HBonds) {
	//	std::cout << neighbour << " ";
	//}
	//std::cout << std::endl;
	for (auto& neighbour : Molecules[Source].HBonds) {
		Path.clear();
		//Path.push_back(Source);
		Path.push_back(neighbour);
		loop(Source, neighbour, Molecules, Path, AllPaths, MaxN);
	}

}

/*unsigned searchRing(std::vector<Water*>& Molecules, unsigned Source, unsigned Nmol, unsigned MaxN) {
	std::vector<unsigned> Distances(Nmol);
	std::vector<bool> SptSet(Nmol);
	Distances.fill(MaxN+1);
	Distances[Source]=MaxN;
	SptSet.fill(false);
	unsigned MinRing{MaxN}, CurrentMol{Source}, CurrentPath{0};
	while (true) {
		for (auto& Neighbour : Molecules[CurrentMol].HBonds) {
			unsigned v {Neighbour -> Index};
			unsigned newDist = CurrentPath+1;
			if (v == Source && newDist< MinRing && newDist > 2 ) {
				MinRing = newDist;
				if (MinRing == 3 ) return 3;
			}
			if (!SptSet(v) && newDist < Distance[v]) {
				Distance[v] = newDist;
				CurrentMol=Neighbour.Index;
				CurrentPath=Distance[CurrentMol];
			}
			if (CurrentMol == Source && )
				break;
			}
		}
		SptSet[CurrentMol]=true;
	}
	for (unsigned Target = 0; Target < Nmol; Target++) {
		// Pick the minimum distance vertex from the set of vertices not
		// yet processed. u is always equal to src in the first iteration.
		unsigned min = MaxN, min_index {};
		for (unsigned u = 0; u < Nmol; u++) {
		     if (sptSet[u] == false && Distance[u] <= min) {
		    	 min = Distance[u];
		    	 min_index = u;
		     }
		}
		//unsigned u = std::min_element(Distances.begin(), Distances.end())-Distances.begin();
		// Mark the picked vertex as processed

		// Update dist value of the adjacent vertices of the picked vertex.
	    for (auto& Neighbour : Molecules[u].HBonds) {
	    	unsigned v {Neighbour -> Index};
	    	unsigned newDist = Distance[u]+1;
	    	if (v == Source && newDist< MinRing && newDist > 2 ) MinRing = newDist;
	    	if (!SptSet[v] && newDist < Distance[v]) Distance[v] = newDist;
	    }
	    if (u != Source) SptSet[u] = true;
		else if (v == Source) {
			if (Distance[v] == 3) return true;

		}

		}
		for (int v = 0; v < V; v++) {

		 // Update dist[v] only if is not in sptSet, there is an edge from
		 // u to v, and total weight of path from src to  v through u is
		 // smaller than current value of dist[v]
		 if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX
									   && dist[u]+graph[u][v] < dist[v])
			dist[v] = dist[u] + graph[u][v];
	 }


}*/

