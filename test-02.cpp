
#include <vector>
#include <iostream>


	
int main ()

{		size_t size_of_round = 0;
		std::vector < std::vector < short int > > directions
			{ { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 } };
		int ii = 0, jj = 0;
		// while ( true )
		for ( size_t g = 0; g < 20; g++ )
		{	size_of_round++;
			for ( size_t d = 0; d < 4; d++ )
			{	if ( d == 2 ) size_of_round++;
				for ( size_t i = 0; i < size_of_round; i++ )
				{	std::cout << ii << " " << jj << std::endl;
					ii += directions[d][0];
					jj += directions[d][1];                      }  }  }
}

