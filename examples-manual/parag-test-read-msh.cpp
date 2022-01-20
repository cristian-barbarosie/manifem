


#include "maniFEM.h"

using namespace maniFEM;

int main ( )

{	Mesh disk ( tag::import, tag::msh, "queijo.msh");
	disk .export_to_file ( tag::msh, "ttor.msh");
}
