#ifndef ZARATH_NS_DIFFUSION_H
#define ZARATH_NS_DIFFUSION_H

#include <ZNAC/LA/Matrix.h>
#include "Iterator.h"

namespace NS
{
	struct Diffusion
	{
		const double Dx, Dy, dt;
		const unsigned int Nx, Ny, margin;
		Diffusion(double Dx, double Dy, double dt, unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Dx(Dx), Dy(Dy), dt(dt), Nx(Nx), Ny(Ny), margin(margin){}

		template<class B>
		void operator()(ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod, const B& boundary_condition)const
		{
			//boundary_condition(dom);
			for(auto &U:Core(Nx, Ny, margin))
				cod[U] = dom[U] - dt*(Dx*(dom[U(1, 0)] + dom[U(-1, 0)] - 2*dom[U]) + Dy*(dom[U(0, 1)] + dom[U(0, -1)] - 2*dom[U]));
		}
	};
}

#endif
