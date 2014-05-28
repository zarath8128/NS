#ifndef ZARATH_NS_ADVECTION_H
#define ZARATH_NS_ADVECTION_H

#include "Grid.h"

namespace NS
{
	void Advection(const Grid &u, const Grid &v, Grid &du, Grid &dv, double Dx, double Dy)
	{
		for(auto &i : u.core)
			du[i] = Dx*u[i]*(u[i(1, 0)] - u[i(-1, 0)]) + 0.25*Dy*(v(i.xi, i.yi - 1) + v(i.xi + 1, i.yi - 1) + v(i.xi, i.yi) + v(i.xi + 1, i.yi))*(u[i(0, 1)] - u[i(0, -1)]);
		for(auto &i : v.core)
			dv[i] = 0.25*Dx*(u(i.xi - 1, i.yi) + u(i.xi - 1, i.yi + 1) + u(i.xi, i.yi) + u(i.xi, i.yi + 1))*(v[i(1, 0)] - v[i(-1, 0)]) + Dy*v[i]*(v[i(0, 1)] - v[i(0, -1)]);
	}
}

#endif
