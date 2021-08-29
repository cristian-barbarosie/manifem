
#include <map>
#include <iostream>

class Function
{	public:
	static bool less_for_map ( const Function &, const Function & )
	{	return true;  }
};

class A { };  class B { };

B operator< ( const A & a, const B & b  ) { std::cout << "1"; return b; }
B operator< ( const B & b, const B & )    { std::cout << "2"; return b; }
A operator< ( const B & b, const A & a )  { std::cout << "3"; return a; }

	
int main ()

{ A a;  B b;
	b < a < b;  //  equivalent to  ( b < a ) < b
	std::cout << "end of main" << std::endl << std::flush;
}

// output 31
