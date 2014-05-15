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
			constexpr GlobalIterator<double> begin()const{return GlobalIterator<double>::begin(g.Nx, g.Ny, g.margin, g);}
			constexpr GlobalIterator<double> end()const{return GlobalIterator<double>::end(g.Nx, g.Ny, g.margin, g);}
		private:
			Grid &g;
		};
		struct Core
		{
			constexpr Core(Grid &g):g(g){}
			constexpr CoreIterator<double> begin() const {return CoreIterator<double>::begin(g.Nx, g.Ny, g.margin, g);}
			constexpr CoreIterator<double> end() const {return CoreIterator<double>::end(g.Nx, g.Ny, g.margin, g);}
		private:
			Grid &g;
		};
		struct Boundary
		{
			constexpr Boundary(Grid &g):g(g){}
			constexpr BoundaryIterator<double> begin() const {return BoundaryIterator<double>::begin(g.Nx, g.Ny, g.margin, g);}
			constexpr BoundaryIterator<double> end() const {return BoundaryIterator<double>::end(g.Nx, g.Ny, g.margin, g);}
		private:
			Grid &g;
		};
		struct Area
		{
			constexpr Area(AreaIndex axi, AreaIndex ayi, Grid &g):axi(axi), ayi(ayi), g(g){}
			constexpr AreaIterator<double> begin() const {return AreaIterator<double>::begin(g.Nx, g.Ny, g.margin, g, axi, ayi);}
			constexpr AreaIterator<double> end() const {return AreaIterator<double>::end(g.Nx, g.Ny, g.margin, g, axi, ayi);}
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
			:Vector<double>((Nx + 2*margin)*(Ny + 2*margin)), NX(Nx + 2*margin), NY(Ny + 2*margin), Nx(Nx), Ny(Ny), margin(margin), 
			global(*this), core(*this), boundary(*this), area(*this)
		{
			for(auto i:global)
				init(i);
		}
	};
}

#endif
