//#include <../lib/Eigen/Eigen/Dense>
using namespace Eigen;

struct Water {
	unsigned Index;
	Vector3d PositionO;
	Vector3d PositionH1;
	Vector3d PositionH2;
	std::list<unsigned> HBonds;
	unsigned MinRing;
	Water(unsigned I) :
		Index{I},
		PositionO{Vector3d::Zero()},
		PositionH1{Vector3d::Zero()},
		PositionH2{Vector3d::Zero()},
		MinRing{100} {}
};
