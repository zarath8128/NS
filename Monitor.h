#ifndef ZARATH_NS_MONITOR_H
#define ZARATH_NS_MONITOR_H

#include <glsc.h>
#include "Grid.h"
#include <ZNAC/LA/LEQSolver.h>

namespace NS
{
	namespace Monitor
	{
		void Initialize(unsigned int size, double left, double right, double bottom, double top, unsigned int Nx, unsigned int Ny);
		void Finalize();

		void Bound();
		void VectorField(const Grid &u, const Grid &v);
		void PressureContln(const Grid &p);
		void ShowData(const Grid &u, const Grid &v, const Grid &p, double Dx, double Dy);
		void ShowCG(ZNAC::LA::CG<double> &cg);
		void PrintWord(const char *format, ...);
		void NewLine();
		void ResetLine();
		void ResetTab();
		void LineTab(int tab, unsigned int tab_stop = 20);
		void SetLineOffset(unsigned int x, unsigned int y);
	}
}

#endif
