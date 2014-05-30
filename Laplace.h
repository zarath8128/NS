#ifndef ZARATH_NS_LAPLACE_H
#define ZARATH_NS_LAPLACE_H

#include <ZNAC/LA/Matrix.h>
#include "Iterator.h"

namespace NS
{
	struct Laplace
	{
		const double Dx, Dy;
		const unsigned int Nx, Ny, margin;
		Laplace(double Dx, double Dy, unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Dx(Dx), Dy(Dy), Nx(Nx), Ny(Ny), margin(margin){}

		void operator()(ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod)const
		{
			for(auto &U:Core(Nx, Ny, margin))
				cod[U] = Dx*(dom[U(1, 0)] + dom[U(-1, 0)]) + Dy*(dom[U(0, 1)] + dom[U(0, -1)]) - 2*(Dx + Dy)*dom[U];
		}
	};

	struct Laplace_d
	{
		const double Dx, Dy;
		const unsigned int Nx, Ny, margin;
		Laplace_d(double Dx, double Dy, unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Dx(Dx), Dy(Dy), Nx(Nx), Ny(Ny), margin(margin){}

		void operator()(ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod)const
		{
		/*	for(auto &U:Core(Nx, Ny, margin))
				cod[U] = 
					Dx*(U.xi == 0 ? dom[U(1, 0)] - dom[U] : U.xi == (int)(Nx - 1) ? dom[U(-1, 0)] - dom[U] : dom[U(1, 0)] + dom[U(-1, 0)] - 2*dom[U]) + 
					Dy*(U.yi == 0 ? dom[U(0, 1)] - dom[U] : U.yi == (int)(Ny - 1) ? dom[U(0, -1)] - dom[U] : dom[U(0, 1)] + dom[U(0, -1)] - 2*dom[U]);
		*/	
			for(auto &U : Domain(Nx, Ny, margin, 1, Nx - 1, 1, Ny - 1))
				cod[U] = Dx*(dom[U(1, 0)] + dom[U(-1, 0)] - 2*dom[U]) + Dy*(dom[U(0, 1)] + dom[U(0, -1)] - 2*dom[U]);
			for(auto &U : Domain(Nx, Ny, margin, 0, 1, 0, 1))
				cod[U] = Dx*(dom[U(1, 0)] - dom[U]) + Dy*(dom[U(0, 1)] - dom[U]);
			for(auto &U : Domain(Nx, Ny, margin, 0, 1, 1, Ny - 1))
				cod[U] = Dx*(dom[U(1, 0)] - dom[U]) + Dy*(dom[U(0, 1)] + dom[U(0, -1)] - 2*dom[U]);
			for(auto &U : Domain(Nx, Ny, margin, 0, 1, Ny - 1, Ny))
				cod[U] = Dx*(dom[U(1, 0)] - dom[U]) + Dy*(dom[U(0, -1)] - dom[U]);
			for(auto &U : Domain(Nx, Ny, margin, Nx - 1, Nx, 0, 1))
				cod[U] = Dx*(dom[U(-1, 0)] - dom[U]) + Dy*(dom[U(0, 1)] - dom[U]);
			for(auto &U : Domain(Nx, Ny, margin, Nx - 1, Nx, 1, Ny - 1))
				cod[U] = Dx*(dom[U(-1, 0)] - dom[U]) + Dy*(dom[U(0, 1)] + dom[U(0, -1)] - 2*dom[U]);
			for(auto &U : Domain(Nx, Ny, margin, Nx - 1, Nx, Ny - 1, Ny))
				cod[U] = Dx*(dom[U(-1, 0)] - dom[U]) + Dy*(dom[U(0, -1)] - dom[U]);
			for(auto &U : Domain(Nx, Ny, margin, 1, Nx - 1, 0, 1))
				cod[U] = Dx*(dom[U(1, 0)] + dom[U(-1, 0)] - 2*dom[U]) + Dy*(dom[U(0, 1)] - dom[U]);
			for(auto &U : Domain(Nx, Ny, margin, 1, Nx - 1, Ny - 1, Ny))
				cod[U] = Dx*(dom[U(1, 0)] + dom[U(-1, 0)] - 2*dom[U]) + Dy*(dom[U(0, -1)] - dom[U]);
	
		}
	};
}

#endif
