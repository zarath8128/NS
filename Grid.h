#ifndef ZARATH_NS_GRID_H
#define ZARATH_NS_GRID_H

#include <ZNAC/LA/Vector.h>

namespace NS
{
	class Grid
		:public ZNAC::LA::Vector<double>
	{
		unsigned int n(int xi, int yi) const {return (yi + margin)*NX + xi + margin;}

		public:
		class Accessor
		{
		public:
			int xi, yi;

			Accessor(int xi, int yi, Grid &g):xi(xi), yi(yi), g(g){}

			double &operator=(const double &val){assert(g.n(xi, yi) < g.N()); return g[g.n(xi, yi)] = val;}
			double &operator()(int dxi, int dyi){assert(g.n(xi + dxi, yi + dyi) < g.N()); return g[g.n(xi + dxi, yi + dyi)];}
			operator double &(){assert(g.n(xi, yi) < g.N()); return g[g.n(xi, yi)];}

		private:
			Grid &g;
		};

		class Iterator
		{
		public:
			Iterator(int xi, int yi, Grid &g):xi(xi), yi(yi), g(g){}
			Accessor operator*(){return Accessor(xi, yi, g);}
			bool operator!=(const Iterator &i)const {return xi != i.xi || yi != i.yi;}
			virtual void operator++() const  = 0;

		protected:
			mutable int xi, yi;
			Grid &g;
		};

		class GlobalIterator
			:public Iterator
		{
		public:
			GlobalIterator(int xi, int yi, Grid &g):Iterator(xi, yi, g){}
			void operator++() const {if(++xi == (int)(g.NX - g.margin))xi = -g.margin, ++yi;}
		};

		const struct Global
		{
			constexpr Global(Grid &g):g(g){}
			GlobalIterator begin()const{return GlobalIterator(-g.margin, -g.margin, g);}
			GlobalIterator end()const{return GlobalIterator(-g.margin, g.NY - g.margin, g);}
		private:
			Grid &g;
		}global;

		class CoreIterator
			:public Iterator
		{
		public:
			CoreIterator(int xi, int yi, Grid &g):Iterator(xi, yi, g){}
			void operator++() const {if(++xi == (int)g.Nx) xi = 0, ++yi;}
		};

		const struct Core
		{
			constexpr Core(Grid &g):g(g){}
			CoreIterator begin() const {return CoreIterator(0, 0, g);}
			CoreIterator end() const {return CoreIterator(0, g.Ny, g);}
		private:
			Grid &g;
		}core;

	public:
		const unsigned int NX, NY, Nx, Ny, margin;

		template<class INIT>
		Grid(unsigned int Nx, unsigned int Ny, unsigned int margin, const INIT &init)
			:Vector<double>((Nx + 2*margin)*(Ny + 2*margin)), global(*this), core(*this), NX(Nx + 2*margin), NY(Ny + 2*margin), Nx(Nx), Ny(Ny), margin(margin)
		{
			for(auto i:global)
				init(i);
		}
	};
}

#endif
