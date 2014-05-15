#ifndef ZARATH_NS_INDEX_H
#define ZARATH_NS_INDEX_H

#include <ZNAC/LA/Vector.h>

namespace NS
{
	struct Index
	{
		int &xi, &yi;
		const unsigned int NX, NY, margin;

		constexpr Index(int &xi, int &yi, unsigned int NX, unsigned int NY, unsigned int margin)
			:xi(xi), yi(yi), NX(NX), NY(NY), margin(margin){}
		operator unsigned int()const{assert((yi + margin)*NX + xi + margin < NX*NY); return (yi + margin)*NX + xi + margin;}
		unsigned int operator()(int dxi, int dyi)const{assert((yi + dyi + margin)*NX + xi + dxi + margin < NX*NY); return (xi + dxi + margin)*NX + yi + dyi + margin;}
	};

	struct WrapperIndex
		:Index
	{
		constexpr WrapperIndex(int &xi, int &yi, unsigned int NX, unsigned int NY, unsigned int margin, double *buf)
			:Index(xi, yi, NX, NY, margin), buf(buf){}

		WrapperIndex &operator=(const double &val){assert((unsigned int)*this < NX*NY); buf[(unsigned int)*this] = val; return *this;}
		operator double()const{assert((unsigned int)*this < NX*NY); return buf[(unsigned int)*this];}
		double &operator()(int dxi, int dyi)const {assert((*(Index*)this)(dxi, dyi) < NX*NY); return buf[(*(Index*)this)(dxi, dyi)];}

	private:
		double *&buf;
	};
}

#endif
