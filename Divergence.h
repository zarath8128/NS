#ifndef ZARATH_NS_DIVERGENCE_H
#define ZARATH_NS_DIVERGENCE_H

#include "Grid.h"

namespace NS
{
	void Divergence(const Grid &u, const Grid &v, Grid &dp, double Dx, double Dy)
	{
		for(auto &i : dp.core)
			dp[i] = Dx*(u(i.xi, i.yi) - u(i.xi - 1, i.yi)) + Dy*(v(i.xi, i.yi) - v(i.xi, i.yi - 1));
	}
}

#endif
