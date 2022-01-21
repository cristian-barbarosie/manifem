


#include "maniFEM.h"

using namespace maniFEM;


std::map < std::pair < int, int >, Mesh > import_msh ( const std::string filename );

int main ( )

{	std::map < std::pair < int, int >, Mesh > meshes = import_msh ("queijo.msh");
	for ( std::map < std::pair < int, int >, Mesh > ::const_iterator it = meshes .begin();
				it != meshes .end(); it++                                                       )
		std::cout << it->first .first << " " << it->first .second << " "
							<< it->second .dim() << std::endl;
	return 0;
}

int main_2 ( )

{	Mesh disk ( tag::import, tag::msh, "queijo.msh");
	disk .export_to_file ( tag::msh, "ttor.msh");
	return 0;
}
