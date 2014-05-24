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

		template<class B>
		void operator()(ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod, const B& boundary_condition)const
		{
			boundary_condition(dom);
			for(auto U:Core(Nx, Ny, margin))
				cod[U] = Dx*(cod[U(1, 0)] + cod[U(-1, 0)]) + Dy*(cod[U(0, 1)] + cod[U(0, -1)]) - 2*(Dx + Dy)*cod[U];
		}
	};
}

#endif
