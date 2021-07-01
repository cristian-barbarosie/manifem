
// example presented in paragraph 10.16 of the manual
// http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf
// builds a MetricTree over a cloud of points in RR2

#include <fstream>
#include <random>
std::ofstream file_ps ("metric-tree.eps");

#include "metric-tree-verbose.h"
// metric-tree-verbose.h draws a PS file
// for normal use, #include "metric-tree.h"

// I'm trying to mimic STL containers but I'm not good at it
// someone more skilled should improve this

class SqDistanceOnRn

{	public:
	inline double operator() ( const std::vector<double> &, const std::vector<double> & );
};
	
//-----------------------------------------------------------------------------------------------//

inline double SqDistanceOnRn::operator()
( const std::vector<double> & u, const std::vector<double> & v )

{	double res = 0.;
	size_t n = u.size();
	assert ( n == v.size() );
	for ( size_t i = 0; i < n; i++ )
	{	double tmp = u[i] - v[i];
		res += tmp*tmp;           }
	return res;                        }

//-----------------------------------------------------------------------------------------------//


int main ()

{	SqDistanceOnRn sq_dist_Rn;
	MetricTree < std::vector < double >, SqDistanceOnRn > cloud ( sq_dist_Rn, 1., 6. );

	double scale_factor = 10.;
	double border = 0;
	double xmin = 0, xmax = 60, ymin = 0, ymax = 30;
	double translation_x = -scale_factor*xmin;
	double translation_y = -scale_factor*ymin;
	file_ps << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl;
	file_ps << "%%Title:                     MetricTree" << std::endl;
	file_ps << "%%BoundingBox:  0 0 " << " " << scale_factor*(xmax-xmin+2*border)
	        << "   " << scale_factor*(ymax-ymin+2*border) << std::endl;
	file_ps << "%%EndComments" << std::endl;
	file_ps << "%%BeginSetup" << std::endl;
	file_ps << "<< /PageSize [" << scale_factor*(xmax-xmin+2*border) << " "
	        << scale_factor*(ymax-ymin+2*border) << "] >> setpagedevice" << std::endl;
	file_ps << "%%EndSetup" << std::endl << std::endl;
	file_ps << "gsave" << std::endl;
	file_ps << translation_x + scale_factor*border << " "
	        << translation_y + scale_factor*border << " translate" << std::endl;
	file_ps << scale_factor << " dup scale" << std::endl << std::endl;
	file_ps << "gsave " << 0.9 / scale_factor << " setlinewidth" << std::endl;
	file_ps << "0.5 setgray" << std::endl;
		
	int n = 30;
	std::default_random_engine random_generator;
	std::set<int> set_of_theta;
	std::uniform_real_distribution<double> distr ( 0., n );
	distr(random_generator); distr(random_generator); distr(random_generator);
	distr(random_generator); distr(random_generator); distr(random_generator);
	distr(random_generator); distr(random_generator); distr(random_generator); 
	for ( int i = 0; i < 75; i++ )
	{ distr(random_generator); distr(random_generator); distr(random_generator);
		double x = 2*distr(random_generator), y = distr(random_generator);
		cloud.add ( { x, y } );  }
	cloud.draw_ps ( file_ps );
	
	std::vector < double > P { 10.9, 4.6 };
	cloud.find_close_neighbours_of ( P, 5 );
	file_ps << "1 0 0 setrgbcolor ";
	file_ps << P[0] << " " << P[1] << " moveto ";
	file_ps << P[0] << " " << P[1] << " 0.3 0 360 arc fill" << std::endl;
	
	file_ps << "grestore" << std::endl;
	file_ps << "grestore" << std::endl << std::endl;
	file_ps << "showpage" << std::endl;
	file_ps << "%%Trailer" << std::endl;
	file_ps << "%EOF" << std::endl;

	cloud.remove ( cloud.root );

	std::cout << "arrows not drawn, see comments in main-10.15.cpp" << std::endl;
	// in order to get arrows, go to 'draw_arrows' in metric-tree-verbose.h and
	// replace 'lineto' by 'Lineto^', then 'make run-10.15',
	// then copy the definitions of 'ArrowHead' and 'Lineto^' from 3d.ps into metric-tree.eps

	std::cout << "produced file metric-tree.eps" << std::endl;
}
