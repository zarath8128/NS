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

		void operator()(ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod)const
		{
			for(auto &U:Core(Nx, Ny, margin))
				cod[U] = dom[U] - dt*(Dx*(dom[U(1, 0)] + dom[U(-1, 0)] - 2*dom[U]) + Dy*(dom[U(0, 1)] + dom[U(0, -1)] - 2*dom[U]));
		}
	};

	struct Bound_Flow
	{
		virtual ~Bound_Flow(){}
		virtual double operator()(int i) const {return i - i;}
	};

	struct Parabola
		:Bound_Flow
	{
		const double a, b; 
		constexpr Parabola(double a, double b, double left, double right, unsigned int Nx, unsigned int margin)
			:a(a), b(b), dx((left - right)/(Nx + 2*margin)), left(left), right(right){}
		virtual double operator()(int i) const
		{
			double x = left + (i + 1)*dx;
			return -a*(x - left)*(x - right) + b;
		}
	private:
		const double dx, left, right;
	};

	struct Diffusion_U
	{
		const double Dx, Dy, dt;
		const unsigned int Nx, Ny, margin;
		Diffusion_U(double Dx, double Dy, double dt, unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Dx(Dx), Dy(Dy), dt(dt), Nx(Nx), Ny(Ny), margin(margin){}

		void operator()(ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod)const
		{
		/*	for(auto &U:Core(Nx, Ny, margin))
				cod[U] = dom[U] - dt*(Dx*(dom[U(1, 0)] + dom[U(-1, 0)] - 2*dom[U]) + Dy*(dom[U(0, 1)] + dom[U(0, -1)] - 2*dom[U]));
		*/
			for(auto &U : Domain(Nx, Ny, margin, 0, Nx, 0, 1))
				cod[U] = dom[U] - dt*(Dx*(dom[U(1, 0)] + dom[U(-1, 0)] - 2*dom[U]) + Dy*(dom[U(0, 1)] + dom[U(0, -1)] - 3*dom[U]));
			for(auto &U : Domain(Nx, Ny, margin, 0, Nx, 1, Ny - 1))
				cod[U] = dom[U] - dt*(Dx*(dom[U(1, 0)] + dom[U(-1, 0)] - 2*dom[U]) + Dy*(dom[U(0, 1)] + dom[U(0, -1)] - 2*dom[U]));
			for(auto &U : Domain(Nx, Ny, margin, 0, Nx, Ny - 1, Ny))
				cod[U] = dom[U] - dt*(Dx*(dom[U(1, 0)] + dom[U(-1, 0)] - 2*dom[U]) + Dy*(dom[U(0, 1)] + dom[U(0, -1)] - 3*dom[U]));
		}
	};

	struct Diffusion_V
	{
		const double Dx, Dy, dt;
		const unsigned int Nx, Ny, margin;
		Diffusion_V(double Dx, double Dy, double dt, unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Dx(Dx), Dy(Dy), dt(dt), Nx(Nx), Ny(Ny), margin(margin){}

		void operator()(ZNAC::LA::IVector<double> &dom, ZNAC::LA::IVector<double> &cod)const
		{
			for(auto &U:Core(Nx, Ny, margin))
				cod[U] = dom[U] - dt*(Dx*(dom[U(1, 0)] + dom[U(-1, 0)] - 2*dom[U]) + Dy*(dom[U(0, 1)] + dom[U(0, -1)] - 2*dom[U]));
		}
	};
}

#endif
