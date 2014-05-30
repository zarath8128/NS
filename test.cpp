#include <iostream>
#include <glsc.h>
#include "Staggered.h"
#include "Diffusion.h"
#include "Laplace.h"
#include "Coordinate.h"
#include "Monitor.h"
#include "Advection.h"
#include "Gradient.h"
#include "Divergence.h"
#include <ZNAC/LA/LEQSolver.h>
#include <cfloat>
#include <cmath>

using namespace NS;
using namespace ZNAC::LA;
using namespace Monitor;

int main()
{
	//Definition of values;
	const unsigned int N = 128;
	const unsigned int Nx = N;
	const unsigned int Ny = N;
	const unsigned int margin = 1;
	const double Re = 3000;
	const double a = 0;
	const double b = 1;
	const double left = -1;
	const double right = 1;
	const double bottom = -1;
	const double top = 1;
	const double width = right - left;
	const double dx = width/Nx;
	const double Dx = 1/(Re*dx*dx);
	const double height = top - bottom;
	const double dy = height/Ny;
	const double Dy = 1/(Re*dy*dy);
	const double dt = 0.005;
	const double ite_eps = 1e-15;
	const unsigned int size = 150;//mm
	
	double t = 0;
	double T = 200;
	int Cut = 0.2/dt;
	int count = 0;

	bool All_in_One = false;
	bool EffectDiffusion = true && !All_in_One;
	bool ImplicitDiffusion = true;
	bool EffectAdvection = true && !All_in_One;
	bool EffectPressure = true && !All_in_One;
	bool Auto = true;
	bool Skip = false;
	bool Capture = false;

	Coordinate ux(-1 + dx, dx);
	auto top_u = [&](double x){return a*(width - x)*(width + x) + b;};
//	auto top_u = [&](double x){return a*cos(x*M_PI) + 1 + b;};

	Grid u(Nx - 1, Ny, margin);
	Grid v(Nx, Ny - 1, margin);
	Grid p(Nx, Ny, margin);
	Grid phi(Nx, Ny, margin);
	Grid utmp(Nx - 1, Ny, margin);
	Grid vtmp(Nx, Ny - 1, margin);
	Grid ptmp(Nx, Ny, margin);
//	Diffusion diffusion_u(Dx, Dy, dt, Nx - 1, Ny, margin);
	Diffusion_U diffusion_u(Dx, Dy, dt, Nx - 1, Ny, margin);
	Diffusion diffusion_v(Dx, Dy, dt, Nx, Ny - 1, margin);
	//Diffusion diffusion_p(Dx, Dy, 0.01, Nx, Ny, margin);
	Laplace laplace_u(Dx, Dy, Nx - 1, Ny, margin);
	Laplace laplace_v(Dx, Dy, Nx, Ny - 1, margin);
	Laplace_d laplace_phi(-1/(dx*dx), -1/(dy*dy), Nx, Ny, margin);
	CG<double> cg_phi(15, N, ite_eps);
	CG<double> cg_u(15, N, ite_eps);
	CG<double> cg_v(15, N, ite_eps);

	auto boundary_u = 
		[&](auto &U)
		{
			for(auto &i:Area(Nx - 1, Ny, margin, AreaIndex::POSITIVE, AreaIndex::ZERO))
				U[i] = 0;
			for(auto &i:Area(Nx - 1, Ny, margin, AreaIndex::NEGATIVE, AreaIndex::ZERO))
				U[i] = 0;
			for(auto &i:Area(Nx - 1, Ny, margin, AreaIndex::ZERO, AreaIndex::POSITIVE))
				U[i] = 2*top_u(ux(i.xi)) - U[i(0, -1)];// + 4/Dy*(top_u(ux(i.xi))*(top_u(ux(i.xi + 1)) - top_u(ux(i.xi - 1)))/(2*dx) - (top_u(ux(i.xi + 1)) + top_u(ux(i.xi - 1)) - 2*top_u(ux(i.xi)))*Dx);
			//	U[i] = 2 - U[i(0, -1)];// + 4/Dy*(top_u(ux(i.xi))*(top_u(ux(i.xi + 1)) - top_u(ux(i.xi - 1)))/(2*dx) - (top_u(ux(i.xi + 1)) + top_u(ux(i.xi - 1)) - 2*top_u(ux(i.xi)))*Dx);
			for(auto &i:Area(Nx - 1, Ny, margin, AreaIndex::ZERO, AreaIndex::NEGATIVE))
				U[i] = -U[i(0, 1)];
		};
	auto diffusion_boundary_u = 
		[&](auto &U)
		{
			for(auto &i : Domain(Nx - 1, Ny, margin, 0, Nx - 1, - 1, 0))
				U[i] = 0;
			for(auto &i : Domain(Nx - 1, Ny, margin, 0, Nx - 1, Ny, Ny+ 1))
				U[i] = 2*top_u(ux(i.xi));
		};

	auto boundary_v = 
		[&](auto &V)
		{
			for(auto &i:Area(Nx, Ny - 1, margin, AreaIndex::POSITIVE, AreaIndex::ZERO))
				V[i] = - V[i(-1, 0)];
				//V[i] = 0;
			for(auto &i:Area(Nx, Ny - 1, margin, AreaIndex::NEGATIVE, AreaIndex::ZERO))
				V[i] = - V[i(1, 0)];
				//V[i] = 0;
			for(auto &i:Area(Nx, Ny - 1, margin, AreaIndex::ZERO, AreaIndex::POSITIVE))
				V[i] = 0;//2 - V[i(0, -1)];
			for(auto &i:Area(Nx, Ny - 1, margin, AreaIndex::ZERO, AreaIndex::NEGATIVE))
				V[i] = 0;//-V[i(0, 1)];
		};

	auto boundary_p = 
		[&](auto &P)
		{
			double min = DBL_MAX;
			for(auto &i : P.core)
				min = min > P[i] ? P[i] : min;
			for(auto &i : P.core)
				P[i] -= min;
			for(auto &i : Area(Nx, Ny, margin, AreaIndex::POSITIVE, AreaIndex::ZERO))
				P[i] = P[i(-1, 0)] + Dx*dx*(2*u(i.xi - 1, i.yi) + 5*u(i.xi - 2, i.yi) + 4*u(i.xi - 3, i.yi) + u(i.xi - 4, i.yi));
			for(auto &i : Area(Nx, Ny, margin, AreaIndex::NEGATIVE, AreaIndex::ZERO))
				P[i] = P[i( 1, 0)] - Dx*dx*(2*u(i.xi, i.yi) - 5*u(i.xi + 1, i.yi) + 4*u(i.xi + 2, i.yi) - u(i.xi + 3, i.yi));
			for(auto &i : Area(Nx, Ny, margin, AreaIndex::ZERO, AreaIndex::POSITIVE))
				P[i] = P[i(0, -1)] + Dy*dy*(2*v(i.xi, i.yi - 1) + 5*v(i.xi, i.yi - 2) + 4*v(i.xi, i.yi - 3) + v(i.xi, i.yi - 4));
			for(auto &i : Area(Nx, Ny, margin, AreaIndex::ZERO, AreaIndex::NEGATIVE))
				P[i] = P[i(0,  1)] - Dy*dy*(2*v(i.xi, i.yi) - 5*v(i.xi, i.yi + 1) + 4*v(i.xi, i.yi + 2) - v(i.xi, i.yi + 3));
		};

	//initialize
	boundary_u(u);
	boundary_v(v);
	boundary_p(p);

	Initialize(size, -width, width, -width, width, Nx, Ny);

	if(Capture)
		g_capture_set("");

	while(t < T)
	{
		if(count++ % Cut == 0)
		{
			g_cls();

			g_line_color(G_BLACK);
			Bound();
			PressureContln(p);
			g_line_color(G_BLACK);
			VectorField(u, v);


			ResetLine();
			PrintWord("t      ");
			PrintWord(": %5.2f", t);
			PrintWord("dt     ");
			PrintWord(": %5.2e", dt);
			NewLine();
			PrintWord("Count");
			PrintWord(": %d", count);
			PrintWord("Cut");
			PrintWord(": %d", Cut);
			NewLine();
			PrintWord("Re     ");
			PrintWord(": %5.2f", Re);
			NewLine();
			PrintWord("a      ");
			PrintWord(": %5.2f", a);
			PrintWord("b      ");
			PrintWord(": %5.2f", b);
			NewLine();
			PrintWord("width  ");
			PrintWord(": %5.2f", width);
			PrintWord("height ");
			PrintWord(": %5.2f", height);
			NewLine();
			PrintWord("Nx     ");
			PrintWord(": %d", Nx);
			PrintWord("Ny     ");
			PrintWord(": %d", Ny);
			NewLine();
			PrintWord("dx     ");
			PrintWord(": %5.2f", dx);
			PrintWord("dy     ");
			PrintWord(": %5.2f", dy);
			NewLine();
			ShowData(u, v, p, 1/dx, 1/dy);
			NewLine();
			PrintWord("U Diffusion CG param");
			NewLine();
			ShowCG(cg_u);
			NewLine();
			PrintWord("V Diffusion CG param");
			NewLine();
			ShowCG(cg_v);
			NewLine();
			PrintWord("P Laplace CG param");
			NewLine();
			ShowCG(cg_phi);
			NewLine();

			if(Capture)
				g_capture();

			if(Auto)
				g_sleep(Skip ? 0.001 : 0.16);
			else
				g_sleep(-1);
		}

		if(EffectDiffusion)
		{
			if(ImplicitDiffusion)
			{
				diffusion_boundary_u(u);
				cg_u(diffusion_u, u, u, u.core);
				cg_v(diffusion_v, v, v, v.core);
			}
			else
			{
				laplace_u(u, utmp);
				laplace_v(v, vtmp);
				for(auto &i : u.core)
					u[i] += dt*utmp[i];
				for(auto &i : v.core)
					v[i] += dt*vtmp[i];
			}
			boundary_u(u);
			boundary_v(v);
		}

		if(EffectAdvection)
		{
			Advection(u, v, utmp, vtmp, 1/(2*dx), 1/(2*dy));
			for(auto &i : u.core)
				u[i] -= dt*utmp[i];
			for(auto &i : v.core)
				v[i] -= dt*vtmp[i];
			boundary_u(u);
			boundary_v(v);
		}

		if(EffectPressure)
		{
			Gradient(p, utmp, vtmp, 1/dx, 1/dy);
			for(auto &i : u.core)
				u[i] -= dt*utmp[i];
			for(auto &i : v.core)
				v[i] -= dt*vtmp[i];
			boundary_u(u);
			boundary_v(v);
			
			Divergence(u, v, ptmp, 1/dx, 1/dy);
			cg_phi(laplace_phi, phi, ptmp, phi.core);

			Gradient(phi, utmp, vtmp, 1/dx, 1/dy);
			for(auto &i : u.core)
				u[i] += utmp[i];
			for(auto &i : v.core)
				v[i] += vtmp[i];
			for(auto &i : p.core)
				p[i] -= 1/dt*phi[i];
			boundary_u(u);
			boundary_v(v);
			boundary_p(p);
		}

		if(All_in_One)
		{
			Advection(u, v, utmp, vtmp, 1/(2*dx), 1/(2*dy));
			cg_u(diffusion_u, u, u, u.core);
			cg_v(diffusion_v, v, v, v.core);
			for(auto &i : u.core)
				u[i] -= dt*utmp[i];
			for(auto &i : v.core)
				v[i] -= dt*vtmp[i];
			Gradient(p, utmp, vtmp, 1/dx, 1/dy);
			for(auto &i : u.core)
				u[i] -= dt*utmp[i];
			for(auto &i : v.core)
				v[i] -= dt*vtmp[i];
			boundary_u(u);
			boundary_v(v);
			Divergence(u, v, ptmp, 1/dx, 1/dy);
			cg_phi(laplace_phi, phi, ptmp, phi.core);
			//cg(diffusion_p, phi, ptmp, phi.core);

			Gradient(phi, utmp, vtmp, 1/dx, 1/dy);
			for(auto &i : u.core)
				u[i] += utmp[i];
			for(auto &i : v.core)
				v[i] += vtmp[i];
			for(auto &i : p.core)
				p[i] -= 1/dt*phi[i];
			boundary_u(u);
			boundary_v(v);
			boundary_p(p);
		}

		t += dt;
	}

	Finalize();

	return 0;
}
