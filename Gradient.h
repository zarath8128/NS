#ifndef ZARATH_NS_GRADIENT_H
#define ZARATH_NS_GRADIENT_H

#include "Grid.h"

namespace NS
{
	void Gradient(const Grid &p, Grid &du, Grid &dv, double Dx, double Dy)
	{
		for(auto &i : du.core)
			du[i] = Dx*(p(i.xi + 1, i.yi) - p(i.xi, i.yi));
		for(auto &i : dv.core)
			dv[i] = Dy*(p(i.xi, i.yi + 1) - p(i.xi, i.yi));
	}
}

#endif
