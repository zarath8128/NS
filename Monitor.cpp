#include "Monitor.h"
#include <cmath>
#include <cstdarg>
#include <cfloat>

namespace NS
{
	namespace Monitor
	{
		static constexpr unsigned int contln_num = 200;

		static struct 
		{
			double left, right, bottom, top;
		} bound;

		static struct
		{
			unsigned int Nx, Ny;
		} n;

		static int line = 0;
		static int column = 0;
		static constexpr unsigned int line_height = 8;
		static constexpr unsigned int column_width = 24;
		static unsigned int line_offset_x;
		static unsigned int line_offset_y;

		static inline double x(int i)
		{
			return bound.left + (i + 0.5)*(bound.right - bound.left)/n.Nx;
		}

		static inline double y(int i)
		{
			return bound.bottom + (i + 0.5)*(bound.top - bound.bottom)/n.Ny;
		}

		static inline void arrow(double x, double y, double vx, double vy)
		{
			g_arrow(x, y, vx, vy, 0.9*(bound.right - bound.left)/n.Nx*sqrt(vx*vx + vy*vy), 0.02);
		}

		static constexpr double sq(double x){return x*x;}

		static inline void color(double &r, double &g, double &b, double val, double min, double max)
		{
			const double center = 0.5*(min + max);

			if(val < center)
			{
				/*b = 1 - (val - min)/(center - min);
				g = (val - min)/(center - min);
				r = 0;
				*/
				b = cos((val - min)*0.5*M_PI/(center - min));
				r = 0;
			}
			else if(val < max)
			{
				/*b = 0;
				g = 1 - (val - center)/(max - center);
				r = (val - center)/(max - center);
				*/
				b = 0;
				r = cos(0.5*M_PI*(val - max)/(max - center));
			}
				g = cos(M_PI*(val - center)/(max - min));

/*			b = exp(-a*sq((val - min)));
			g = exp(-a*sq((val - center)));
			r = exp(-a*sq((val - max)));
			double len = 1/sqrt(r*r + b*b + g*g);
			b *= len;
			g *= len;
			r *= len;
*/		}

		void Initialize(unsigned int size, double left, double right, double bottom, double top, unsigned int Nx, unsigned int Ny)
		{
			bound.left 		 = left;
			bound.right		 = right;
			bound.bottom	 = bottom;
			bound.top		 = top;

			n.Nx 			 = Nx;
			n.Ny			 = Ny;

			line_offset_x	 = size;
			line_offset_y	 = line_height;

			g_init("dammy", 1.6*size, size);
			g_device(G_DISP);
			g_def_scale(0, left, right, bottom, top, 0.05*size, 0.05*size, size*0.9, size*0.9);
			g_sel_scale(0);

			g_def_text(0, G_BLACK, G_FONT_TIMES_18);
			g_sel_text(0);
		}

		void Finalize()
		{
			g_term();
		}

		void Bound()
		{
			g_move(bound.right, bound.top);
			g_plot(bound.right, bound.bottom);
			g_plot(bound.left, bound.bottom);
			g_plot(bound.left, bound.top);
		}

		void VectorField(const Grid &u, const Grid &v)
		{
			for(unsigned int j = 0; j < n.Ny; ++j)
				for(unsigned int i = 0; i < n.Nx; ++i)
					arrow(x(i), y(j), 0.5*(u(i - 1, j) + u(i, j)), 0.5*(v(i, j - 1) + v(i, j)));
		}

		void PressureContln(const Grid &p)
		{
			double max = -DBL_MAX;
			for(auto &i : p.core)
				max = max < p[i] ? p[i] : max;

			double *P = (double *)malloc(sizeof(double)*n.Nx*n.Ny);

			for(auto &i : p.core)
				P[i.xi*n.Ny + i.yi] = p[i];

			double r, g, b;

			for(unsigned int i = 0; i < contln_num; ++i)
			{
				double lev = max*i/(contln_num - 1);
				color(r, g, b, lev, 0, max);
				g_line_color(g_rgb_color(r, g, b));
				g_contln(bound.left, bound.right, bound.bottom, bound.top, (double *)P, n.Nx, n.Ny, max*i/(contln_num - 1));
			}
			free(P);
		}

		void ShowData(const Grid &u, const Grid &v, const Grid &p, double Dx, double Dy)
		{
			double div_max = -DBL_MAX, div_min = DBL_MAX;
			double p_max = 0;
			double div_ave = 0, div_all = 0;
			for(auto &i : p.core)
			{
				p_max = p[i] > p_max ? p[i] : p_max;
				double div = Dx*(u(i.xi, i.yi) - u(i.xi - 1, i.yi)) + Dy*(v(i.xi, i.yi) - v(i.xi, i.yi - 1));
				div_max = div > div_max ? div : div_max;
				div_min = div < div_min ? div : div_min;
				div_all += ZNAC::ABS(div);
				div_ave += div;
			}

			PrintWord("Div Max");
			PrintWord(": %e", div_max);
			PrintWord("Div Min");
			PrintWord(": %e", div_min);
			NewLine();
			PrintWord("Div Ave");
			PrintWord(": %e", div_ave);
			PrintWord("Div Err");
			PrintWord(": %e", div_all);
			NewLine();
			PrintWord("Pre Max");
			PrintWord(": %5.2f", p_max);

		}

		void ShowCG(ZNAC::LA::CG<double> &cg)
		{
			PrintWord("repeat");
			PrintWord(": %d", cg.rep);
			PrintWord("total");
			PrintWord(": %d", cg.n);
			NewLine();
			PrintWord("state");
			PrintWord(": %s", cg ? "success" : "fail");
		}

		void G_Text(const char *text)
		{
			g_text(line_offset_x + column*column_width + column/2*column_width/2, line_offset_y + line*line_height, text);
		}

		void PrintLine(const char *format, ...)
		{
			va_list param;
			va_start(param, format);
			char buf[128];
			vsnprintf(buf, 127, format, param);
			va_end(param);
			G_Text(buf);
			line++;
			column = 0;
		}

		void PrintWord(const char *format, ...)
		{
			va_list param;
			va_start(param, format);
			char buf[128];
			vsnprintf(buf, 127, format, param);
			va_end(param);
			G_Text(buf);
			column++;
		}

		void NewLine()
		{
			line++;
			column = 0;
		}

		void ResetLine()
		{
			column = 0;
			line = 0;
		}

		void SetLineOffset(unsigned int x, unsigned int y)
		{
			line_offset_x = x;
			line_offset_y = y;
		}
	}
}
