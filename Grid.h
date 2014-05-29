#ifndef ZARATH_NS_GRID_H
#define ZARATH_NS_GRID_H

#include <ZNAC/LA/Vector.h>
#include "Iterator.h"

namespace NS
{
	struct Grid
		:public ZNAC::LA::Vector<double>
	{
		struct Global
		{
			constexpr Global(Grid &g):g(g){}
			constexpr GlobalIterator begin()const{return GlobalIterator::begin(g.Nx, g.Ny, g.margin);}
			constexpr GlobalIterator end()const{return GlobalIterator::end(g.Nx, g.Ny, g.margin);}
		private:
			Grid &g;
		};
		struct Boundary
		{
			constexpr Boundary(Grid &g):g(g){}
			constexpr BoundaryIterator begin() const {return BoundaryIterator::begin(g.Nx, g.Ny, g.margin);}
			constexpr BoundaryIterator end() const {return BoundaryIterator::end(g.Nx, g.Ny, g.margin);}
		private:
			Grid &g;
		};
		struct Area
		{
			constexpr Area(AreaIndex axi, AreaIndex ayi, Grid &g):axi(axi), ayi(ayi), g(g){}
			constexpr AreaIterator begin() const {return AreaIterator::begin(g.Nx, g.Ny, g.margin, axi, ayi);}
			constexpr AreaIterator end() const {return AreaIterator::end(g.Nx, g.Ny, g.margin, axi, ayi);}
		private:
			const AreaIndex axi, ayi;
			Grid &g;
		};
		struct AreaFactory
		{
			constexpr AreaFactory(Grid &g):g(g){}
			constexpr Area operator()(AreaIndex axi, AreaIndex ayi){return Area(axi, ayi, g);}
		private:
			Grid &g;
		};
	public:
		const unsigned int NX, NY, Nx, Ny, margin;
		const Global global;
		const Core core;
		const Boundary boundary;
		const AreaFactory area;

		template<class INIT>
		Grid(unsigned int Nx, unsigned int Ny, unsigned int margin, const INIT &init)
			:Vector((Nx + 2*margin)*(Ny + 2*margin)), NX(Nx + 2*margin), NY(Ny + 2*margin), Nx(Nx), Ny(Ny), margin(margin), 
			global(*this), core(Nx, Ny, margin), boundary(*this), area(*this), n(0, 0, Nx, Ny, margin)
		{
			for(auto &i:global)
				init(*this, i);
		}
		Grid(unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Grid(Nx, Ny, margin, [](auto &v, auto &i){v[i] = 0;}){}

		double &operator()(int xi, int yi){return (*this)[n.n(xi, yi)];}
		const double &operator()(int xi, int yi)const{return (*this)[n.n(xi, yi)];}

	private:
		Index n;
	};
}

#endif
