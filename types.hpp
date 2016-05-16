#pragma once
/** This is geometrical types for Moroz Anton Diploma */


#include <vector>
#include <cmath>

#include <Eigen/Eigen>

namespace Dimploma
{
	typedef Eigen::VectorXd Point;
	typedef Eigen::MatrixXd PointsVec;
	typedef PointsVec HyperCube;

	double dist(Point A, Point B)
	{
		return (A - B).norm();
	}

} //namespace Dimploma
