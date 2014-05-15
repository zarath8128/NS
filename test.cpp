#include <iostream>
#include <glsc.h>
//#include "Index.h"
#include "Grid.h"
#include "Iterator.h"

using namespace NS;

int main()
{
	Grid g(8, 8, 1, [](auto &i){i = i.xi + 10*i.yi;});
	Accessor<double> a(0, 0, 8, 8, 1, g);
	for(auto x:g.global)
	{
		if(x.xi == -1)
			std::cout << std::endl;
		std::cout << "\t" << x;
	}
	std::cout << std::endl;
	for(auto x:g.core)
	{
		if(x.xi == 0)
			std::cout << std::endl;
		std::cout << "\t" << x;
	}
	std::cout << std::endl;
	for(auto x:g.boundary)
	{
		if(x.xi == -1)
			std::cout << std::endl;
		std::cout << "\t" << x;
	}
	for(int i = -1; i < 2; ++i)
		for(int j = -1; j < 2; ++j)
		{
			std::cout << std::endl << std::endl;
			for(auto x:g.area(AreaIndex(i), AreaIndex(j)))
				std::cout << x << std::endl;
		}
	return 0;
}
