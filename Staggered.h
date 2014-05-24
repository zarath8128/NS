#ifndef ZARATH_NS_STAGGERED_H
#define ZARATH_NS_STAGGERED_H

#include "Grid.h"

namespace NS
{
	class Staggered
	{
	public:
		Grid u, v, p;
		template<class INITU, class INITV, class INITP>
		Staggered(unsigned int Nx, unsigned int Ny, const INITU &init_u, const INITV &init_v, const INITP &init_p)
			:u(Nx + 1, Ny, 1, init_u), v(Nx, Ny + 1, 1, init_v), p(Nx, Ny, 1, init_p){}
		Staggered(unsigned int Nx, unsigned int Ny):Staggered(Nx, Ny, [](auto &v, auto &i){v[i] = 0;}, [](auto &v, auto &i){v[i] = 0;}, [](auto &v, auto &i){v[i] = 0;}){}
	};
}

#endif
