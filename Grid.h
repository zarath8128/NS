#ifndef ZARATH_NS_GRID_H
#define ZARATH_NS_GRID_H

#include <ZNAC/LA/Vector.h>
#include "Iterator.h"

namespace NS
{
	struct Grid
		:public ZNAC::LA::Vector<double>
	{
	public:
		const unsigned int NX, NY, Nx, Ny, margin;
		const Global global;
		const Core core;
		const Boundary boundary;
		const AreaFactory area;
		const DomainFactory domain;

		template<class INIT>
		Grid(unsigned int Nx, unsigned int Ny, unsigned int margin, const INIT &init)
			:Vector((Nx + 2*margin)*(Ny + 2*margin)), NX(Nx + 2*margin), NY(Ny + 2*margin), Nx(Nx), Ny(Ny), margin(margin), 
			global(Nx, Ny, margin), core(Nx, Ny, margin), boundary(Nx, Ny, margin), area(Nx, Ny, margin), domain(Nx, Ny, margin), n(0, 0, Nx, Ny, margin)
		{
			for(auto &i:global)
				init(*this, i);
		}
		Grid(unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Grid(Nx, Ny, margin, [](auto &v, auto &i){v[i] = 0;}){}

		double &operator()(int xi, int yi){return (*this)[n.n(xi, yi)];}
		const double &operator()(int xi, int yi)const{return (*this)[n.n(xi, yi)];}
		bool operator==(const Grid &g) const {return Nx == g.Nx && Ny == g.Ny && margin == g.margin;}

	private:
		Index n;
	};
}

#endif
