#ifndef ZARATH_NS_COORDINATE_H
#define ZARATH_NS_COORDINATE_H

namespace NS
{
	class Coordinate
	{
	public:
		constexpr Coordinate(double origin, double delta)
			:origin(origin), delta(delta){}

		constexpr double operator()(int i) const {return origin + i*delta;}

	private:
		const double origin, delta;
	};
}

#endif
