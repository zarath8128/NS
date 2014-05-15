#include <iostream>
#include <glsc.h>
#include "Index.h"
#include "Grid.h"

using namespace NS;

int main()
{
	Grid g(16, 16, 1, [](auto &i){i = i.xi + i.yi;});
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
	return 0;
}
