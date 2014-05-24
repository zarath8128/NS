#include <iostream>
#include <glsc.h>
#include "Staggered.h"
#include "Diffusion.h"
#include "Laplace.h"
#include "Coordinate.h"
#include <ZNAC/LA/LEQSolver.h>

using namespace NS;
using namespace ZNAC::LA;

int main()
{
	constexpr unsigned int N = 1024;
	constexpr unsigned int Nx = N;
	constexpr unsigned int Ny = N;
	constexpr unsigned int margin = 1;
	constexpr double Re = 10;
	constexpr double a = 0;
	constexpr double b = 1;
	constexpr double width = 1;
	constexpr double dx = 2*width/Nx;
	constexpr double Dx = 1/(Re*dx*dx);
	constexpr double height = 1;
	constexpr double dy = 2*height/Ny;
	constexpr double Dy = 1/(Re*dy*dy);
	constexpr double dt = 10;
	constexpr double ite_eps = 1e-15;

	Coordinate ux(-1, dx);
	auto top_u = [&](double x){return a*(width - x)*(width + x) + b;};
	auto bottom_u = [&](double x){return 0;};
	auto left_v = [&](double x){return 0;};
	auto right_v = [&](double x){return 0;};

	auto initial = [&](auto &v, auto &i){v[i] = 0;};

	Staggered nssys(Nx, Ny, initial, initial, initial);

	Grid utmp(Nx + 1, Ny, margin);
	Laplace laplace_u(Dx, Dy, Nx + 1, Ny, margin);
	Diffusion diffusion_u(Dx, Dy, dt, Nx + 1, Ny, margin);
	CG<double> cg(Nx, Nx, ite_eps);

	auto boundary_u = 
		[&](auto &U)
		{
			for(auto &u:Area(Nx + 1, Ny, margin, AreaIndex::POSITIVE, AreaIndex::ZERO))
				U[u] = top_u(ux(u.xi));
				//U[u(-1, 0)] = 0, U[u] = - U[u(-2, 0)];
				//U[u] = U[u(-1, 0)];// = 0;
				//U[u] = 0;
			for(auto &u:Area(Nx + 1, Ny, margin, AreaIndex::NEGATIVE, AreaIndex::ZERO))
				//U[u(1, 0)] = 0, U[u] = - U[u(2, 0)];
				U[u] = U[u(1, 0)];// = 0;
				//U[u] = 0;
			for(auto &u:Area(Nx + 1, Ny, margin, AreaIndex::ZERO, AreaIndex::POSITIVE))
				U[u] = 2*top_u(ux(u.xi)) - U[u(0, -1)];// + 4/Dy*(top(ux(u.xi))*(top(ux(u.xi + 1)) - top(ux(u.xi - 1)))/(2*dx) - (top(ux(u.xi + 1)) + top(ux(u.xi - 1)) - 2*top(ux(u.xi)))*Dx);
			for(auto &u:Area(Nx + 1, Ny, margin, AreaIndex::ZERO, AreaIndex::NEGATIVE))
				U[u] = U[u(0, 1)];
				//U[u] = 1;//U[u(0, 1)];
		};
	//auto boundary_v

	boundary_u(nssys.u);


	boundary_u(utmp);
	for(int i = 0; i < 1; ++i)
	{
		cg(diffusion_u, utmp, nssys.u, boundary_u, nssys.u.core);
//		for(auto &u:utmp.global)
//		{
//			if(u.xi == -1)
//				std::cout << std::endl;
//			std::cout << utmp[u] << "\t";
//		}
//		std::cout << "\n\n";
		for(auto &u:nssys.u.global)
			nssys.u[u] = utmp[u];
	}
//	boundary_condition(nssys.u);
//	for(auto &u:nssys.u.global)
//	{
//		if(u.xi == -1)
//			std::cout << std::endl;
//		std::cout << nssys.u[u] << "\t";
//	}
//	std::cout << "\n\n";
	return 0;
}
