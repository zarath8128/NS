#ifndef ZARATH_NS_ITERATOR_H
#define ZARATH_NS_ITERATOR_H

#include <ZNAC/LA/Vector.h>

namespace NS
{
	struct Index
	{
		mutable int xi, yi;
		const unsigned int Nx, Ny, margin;
		
		constexpr Index(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin)
			:xi(xi), yi(yi), Nx(Nx), Ny(Ny), margin(margin){}
		constexpr Index(const Index &i)
			:xi(i.xi), yi(i.yi), Nx(i.Nx), Ny(i.Ny), margin(i.margin){}

		constexpr operator unsigned int() const {return index();}
		constexpr unsigned int index(int dxi = 0, int dyi = 0)const{return (yi + dyi + margin)*(Nx + 2*margin) + xi + dxi + margin;}
		constexpr unsigned int operator()(int dxi, int dyi)const{return index(dxi, dyi);}
		constexpr unsigned int n(int xi, int yi)const{return (yi + margin)*(Nx + 2*margin) + xi + margin;}
	};

	static inline std::ostream &operator<<(std::ostream &dest, const Index& i)
	{
		dest << "xi : " << i.xi << "\tyi : " << i.yi << "\ti : " << (unsigned int)i;
		return dest;
	}

	struct Iterator
		:Index
	{
		constexpr Iterator(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Index(xi, yi, Nx, Ny, margin){}

		constexpr const Iterator& operator*()const{return *this;}
		constexpr bool operator!=(const Index &i)const{return xi != i.xi || yi != i.yi;}
		virtual void operator++() const = 0;
	};

	struct GlobalIterator
		:Iterator
	{
		constexpr GlobalIterator(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Iterator(xi, yi, Nx, Ny, margin){}
		void operator++() const {if(++xi == (int)(Nx + margin))xi = -margin, ++yi;}
		static constexpr GlobalIterator begin(unsigned int Nx, unsigned int Ny, unsigned int margin){return GlobalIterator(-margin, -margin, Nx, Ny, margin);}
		static constexpr GlobalIterator end(unsigned int Nx, unsigned int Ny, unsigned int margin){return GlobalIterator(-margin, Ny + margin, Nx, Ny, margin);}
	};

	struct Global
	{
		const unsigned int Nx, Ny, margin;
		constexpr Global(unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Nx(Nx), Ny(Ny), margin(margin){}
		constexpr GlobalIterator begin() const {return GlobalIterator::begin(Nx, Ny, margin);}
		constexpr GlobalIterator end() const {return GlobalIterator::end(Nx, Ny, margin);}
	};

	struct CoreIterator
		:Iterator
	{
		constexpr CoreIterator(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Iterator(xi, yi, Nx, Ny, margin){}
		void operator++() const {if(++xi == (int)Nx)xi = 0, ++yi;}
		static constexpr CoreIterator begin(unsigned int Nx, unsigned int Ny, unsigned int margin){return CoreIterator(0, 0, Nx, Ny, margin);}
		static constexpr CoreIterator end(unsigned int Nx, unsigned int Ny, unsigned int margin){return CoreIterator(0, Ny, Nx, Ny, margin);}
	};

	struct Core
	{
		const unsigned int Nx, Ny, margin;
		constexpr Core(unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Nx(Nx), Ny(Ny), margin(margin){}
		constexpr CoreIterator begin()const{return CoreIterator::begin(Nx, Ny, margin);}
		constexpr CoreIterator end()const{return CoreIterator::end(Nx, Ny, margin);}
	};

	struct BoundaryIterator
		:Iterator
	{
		constexpr BoundaryIterator(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Iterator(xi, yi, Nx, Ny, margin){}
		void operator++() const
		{
			if(yi < 0 || yi > (int)Ny - 1)
				if(++xi == int(Nx + margin))
					xi = -margin, ++yi;
				else{}
			else
				if(++xi == 0)
					xi = Nx;
				else if(xi == int(Nx + margin))
					xi = -margin, ++yi;
		}
		static constexpr BoundaryIterator begin(unsigned int Nx, unsigned int Ny, unsigned int margin){return BoundaryIterator(-margin, -margin, Nx, Ny, margin);}
		static constexpr BoundaryIterator end(unsigned int Nx, unsigned int Ny, unsigned int margin){return BoundaryIterator(-margin, Ny + margin, Nx, Ny, margin);}
	};

	struct Boundary
	{
		const unsigned int Nx, Ny, margin;
		constexpr Boundary(unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Nx(Nx), Ny(Ny), margin(margin){}
		constexpr BoundaryIterator begin()const{return BoundaryIterator::begin(Nx, Ny, margin);}
		constexpr BoundaryIterator end()const{return BoundaryIterator::end(Nx, Ny, margin);}
	};

	enum struct AreaIndex:int
	{
		NEGATIVE 	= -1,
		ZERO 		= 0,
		POSITIVE	= 1,
	};

	struct AreaIterator
		:Iterator
	{
		const AreaIndex axi, ayi;
		constexpr AreaIterator(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin, AreaIndex axi, AreaIndex ayi)
			:Iterator(xi, yi, Nx, Ny, margin), axi(axi), ayi(ayi){}
		void operator++() const
		{if(++xi == XEnd())xi = XBegin(), ++yi;}
		static constexpr AreaIterator begin(unsigned int Nx, unsigned int Ny, unsigned int margin, AreaIndex axi, AreaIndex ayi)
		{
			return AreaIterator(XBegin(axi, Nx, margin), YBegin(ayi, Ny, margin), Nx, Ny, margin, axi, ayi);
		}
		static constexpr AreaIterator end(unsigned int Nx, unsigned int Ny, unsigned int margin, AreaIndex axi, AreaIndex ayi)
		{
			return AreaIterator(XBegin(axi, Nx, margin), YEnd(ayi, Ny, margin), Nx, Ny, margin, axi, ayi);
		}

	private:
		constexpr int XBegin()const{return axi == AreaIndex::NEGATIVE ? -margin : axi == AreaIndex::ZERO ? 0 : Nx;}
		constexpr int XEnd()const{return axi == AreaIndex::NEGATIVE ? 0 : axi == AreaIndex::ZERO ? Nx : Nx + margin;}
		constexpr int YBegin()const{return ayi == AreaIndex::NEGATIVE ? -margin : ayi == AreaIndex::ZERO ? 0 : Ny;}
		constexpr int YEnd()const{return ayi == AreaIndex::NEGATIVE ? 0 : ayi == AreaIndex::ZERO ? Ny : Ny + margin;}
		constexpr static int XBegin(AreaIndex axi, unsigned int Nx, unsigned int margin){return axi == AreaIndex::NEGATIVE ? -margin : axi == AreaIndex::ZERO ? 0 : Nx;}
		constexpr static int XEnd(AreaIndex axi, unsigned int Nx, unsigned int margin){return axi == AreaIndex::NEGATIVE ? 0 : axi == AreaIndex::ZERO ? Nx : Nx + margin;}
		constexpr static int YBegin(AreaIndex ayi, unsigned int Ny, unsigned int margin){return ayi == AreaIndex::NEGATIVE ? -margin : ayi == AreaIndex::ZERO ? 0 : Ny;}
		constexpr static int YEnd(AreaIndex ayi, unsigned int Ny, unsigned int margin){return ayi == AreaIndex::NEGATIVE ? 0 : ayi == AreaIndex::ZERO ? Ny : Ny + margin;}
	};

	struct Area
	{
		const unsigned int Nx, Ny, margin;
		const AreaIndex axi, ayi;
		constexpr Area(unsigned int Nx, unsigned int Ny, unsigned int margin, AreaIndex axi, AreaIndex ayi)
			:Nx(Nx), Ny(Ny), margin(margin), axi(axi), ayi(ayi){}
		constexpr AreaIterator begin()const{return AreaIterator::begin(Nx, Ny, margin, axi, ayi);}
		constexpr AreaIterator end()const{return AreaIterator::end(Nx, Ny, margin, axi, ayi);}
	};

	struct AreaFactory
	{
		const unsigned int Nx, Ny, margin;
		constexpr AreaFactory(unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Nx(Nx), Ny(Ny), margin(margin){}
		constexpr Area operator()(AreaIndex axi, AreaIndex ayi){return Area(Nx, Ny, margin, axi, ayi);}
	};

	struct DomainIterator
		: Iterator
	{
		//x_begin <= xi < x_end
		//y_begin <= yi < y_end
		const int x_begin, x_end, y_begin, y_end;
		constexpr DomainIterator(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin, int x_begin, int x_end, int y_begin, int y_end)
			:Iterator(xi, yi, Nx, Ny, margin), x_begin(x_begin), x_end(x_end), y_begin(y_begin), y_end(y_end){}
		void operator++() const {if(++xi == x_end)xi = x_begin, ++yi;}
		static constexpr DomainIterator begin(unsigned int Nx, unsigned int Ny, unsigned int margin, int x_begin, int x_end, int y_begin, int y_end)
		{return DomainIterator(x_begin, y_begin, Nx, Ny, margin, x_begin, x_end, y_begin, y_end);}
		static constexpr DomainIterator end(unsigned int Nx, unsigned int Ny, unsigned int margin, int x_begin, int x_end, int y_begin, int y_end)
		{return DomainIterator(x_begin, y_end, Nx, Ny, margin, x_begin, x_end, y_begin, y_end);}
	};

	struct Domain
	{
		const unsigned int Nx, Ny, margin;
		const int x_begin, x_end, y_begin, y_end;

		constexpr Domain(unsigned int Nx, unsigned int Ny, unsigned int margin, int x_begin, int x_end, int y_begin, int y_end)
			:Nx(Nx), Ny(Ny), margin(margin), x_begin(x_begin), x_end(x_end), y_begin(y_begin), y_end(y_end){}
		constexpr DomainIterator begin() const {return DomainIterator::begin(Nx, Ny, margin, x_begin, x_end, y_begin, y_end);}
		constexpr DomainIterator end() const {return DomainIterator::end(Nx, Ny, margin, x_begin, x_end, y_begin, y_end);}
	};

	struct DomainFactory
	{
		const unsigned int Nx, Ny, margin;
		constexpr DomainFactory(unsigned int Nx, unsigned int Ny, unsigned int margin)
			:Nx(Nx), Ny(Ny), margin(margin){}
		Domain operator()(int x_begin, int x_end, int y_begin, int y_end)const{assert(x_begin < x_end); assert(y_begin < y_end); return Domain(Nx, Ny, margin, x_begin, x_end, y_begin, y_end);}
	};
}

#endif
