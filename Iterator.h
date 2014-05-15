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

		constexpr unsigned int index(int dxi = 0, int dyi = 0)const{return (yi + dyi + margin)*(Nx + 2*margin) + xi + dxi + margin;}
	};

	template<class T>
	struct Accessor
		:Index
	{

		constexpr Accessor(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<T> &v)
			:Index(xi, yi, Nx, Ny, margin), v(v){}
		constexpr Accessor(const Index &i, ZNAC::LA::IVector<T> &v)
			:Index(i), v(v){}

		constexpr operator double &(){return v[index()];}
		constexpr T &operator=(const T &val){return v[index()] = val;}
		constexpr T &operator()(int dxi, int dyi){return v[index(dxi, dyi)];}


	private:
		ZNAC::LA::IVector<T> &v;
	};

	template<class T>
	struct Iterator
		:Index
	{
		constexpr Iterator(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<T> &v)
			:Index(xi, yi, Nx, Ny, margin), v(v){}

		constexpr Accessor<T> operator*()const{return Accessor<T>(*this, v);}
		constexpr bool operator!=(const Index &i)const{return xi != i.xi || yi != i.yi;}
		virtual void operator++() const = 0;

	private:
		ZNAC::LA::IVector<T> &v;
	};

	template<class T>
	struct GlobalIterator
		:Iterator<T>
	{
		constexpr GlobalIterator(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<double> &v)
			:Iterator<T>(xi, yi, Nx, Ny, margin, v){}
		constexpr void operator++() const {if(++Iterator<T>::xi == (int)(Iterator<T>::Nx + Iterator<T>::margin))Iterator<T>::xi = -Iterator<T>::margin, ++Iterator<T>::yi;}
		static constexpr GlobalIterator begin(unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<T> &v){return GlobalIterator(-margin, -margin, Nx, Ny, margin, v);}
		static constexpr GlobalIterator end(unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<T> &v){return GlobalIterator(-margin, Ny + margin, Nx, Ny, margin, v);}
	};

	template<class T>
	struct CoreIterator
		:Iterator<T>
	{
		constexpr CoreIterator(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<double> &v)
			:Iterator<T>(xi, yi, Nx, Ny, margin, v){}
		constexpr void operator++() const {if(++Iterator<T>::xi == (int)Iterator<T>::Nx)Iterator<T>::xi = 0, ++Iterator<T>::yi;}
		static constexpr CoreIterator begin(unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<T> &v){return CoreIterator(0, 0, Nx, Ny, margin, v);}
		static constexpr CoreIterator end(unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<T> &v){return CoreIterator(0, Ny, Nx, Ny, margin, v);}
	};

	template<class T>
	struct BoundaryIterator
		:Iterator<T>
	{
		constexpr BoundaryIterator(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<double> &v)
			:Iterator<T>(xi, yi, Nx, Ny, margin, v){}
		constexpr void operator++() const
		{
			if(Iterator<T>::yi < 0 || Iterator<T>::yi > (int)Iterator<T>::Ny - 1)
				if(++Iterator<T>::xi == int(Iterator<T>::Nx + Iterator<T>::margin))
					Iterator<T>::xi = -Iterator<T>::margin, ++Iterator<T>::yi;
				else{}
			else
				if(++Iterator<T>::xi == 0)
					Iterator<T>::xi = Iterator<T>::Nx;
				else if(Iterator<T>::xi == int(Iterator<T>::Nx + Iterator<T>::margin))
					Iterator<T>::xi = -Iterator<T>::margin, ++Iterator<T>::yi;
		}
		static constexpr BoundaryIterator begin(unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<T> &v){return BoundaryIterator(-margin, -margin, Nx, Ny, margin, v);}
		static constexpr BoundaryIterator end(unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<T> &v){return BoundaryIterator(-margin, Ny + margin, Nx, Ny, margin, v);}
	};

	enum struct AreaIndex:int
	{
		NEGATIVE 	= -1,
		ZERO 		= 0,
		POSITIVE	= 1,
	};

	template<class T>
	struct AreaIterator
		:Iterator<T>
	{
		const AreaIndex axi, ayi;
		constexpr AreaIterator(int xi, int yi, unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<T> &v, AreaIndex axi, AreaIndex ayi)
			:Iterator<T>(xi, yi, Nx, Ny, margin, v), axi(axi), ayi(ayi){}
		constexpr void operator++() const
		{if(++Iterator<T>::xi == XEnd())Iterator<T>::xi = XBegin(), ++Iterator<T>::yi;}
		static constexpr AreaIterator begin(unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<T> &v, AreaIndex axi, AreaIndex ayi)
		{
			return AreaIterator(XBegin(axi, Nx, margin), YBegin(ayi, Ny, margin), Nx, Ny, margin, v, axi, ayi);
		}
		static constexpr AreaIterator end(unsigned int Nx, unsigned int Ny, unsigned int margin, ZNAC::LA::IVector<T> &v, AreaIndex axi, AreaIndex ayi)
		{
			return AreaIterator(XBegin(axi, Nx, margin), YEnd(ayi, Ny, margin), Nx, Ny, margin, v, axi, ayi);
		}

	private:
		constexpr int XBegin()const{return axi == AreaIndex::NEGATIVE ? -Iterator<T>::margin : axi == AreaIndex::ZERO ? 0 : Iterator<T>::Nx;}
		constexpr int XEnd()const{return axi == AreaIndex::NEGATIVE ? 0 : axi == AreaIndex::ZERO ? Iterator<T>::Nx : Iterator<T>::Nx + Iterator<T>::margin;}
		constexpr int YBegin()const{return ayi == AreaIndex::NEGATIVE ? -Iterator<T>::margin : ayi == AreaIndex::ZERO ? 0 : Iterator<T>::Ny;}
		constexpr int YEnd()const{return ayi == AreaIndex::NEGATIVE ? 0 : ayi == AreaIndex::ZERO ? Iterator<T>::Ny : Iterator<T>::Ny + Iterator<T>::margin;}
		constexpr static int XBegin(AreaIndex axi, unsigned int Nx, unsigned int margin){return axi == AreaIndex::NEGATIVE ? -margin : axi == AreaIndex::ZERO ? 0 : Nx;}
		constexpr static int XEnd(AreaIndex axi, unsigned int Nx, unsigned int margin){return axi == AreaIndex::NEGATIVE ? 0 : axi == AreaIndex::ZERO ? Nx : Nx + margin;}
		constexpr static int YBegin(AreaIndex ayi, unsigned int Ny, unsigned int margin){return ayi == AreaIndex::NEGATIVE ? -margin : ayi == AreaIndex::ZERO ? 0 : Ny;}
		constexpr static int YEnd(AreaIndex ayi, unsigned int Ny, unsigned int margin){return ayi == AreaIndex::NEGATIVE ? 0 : ayi == AreaIndex::ZERO ? Ny : Ny + margin;}
	};
}

#endif
