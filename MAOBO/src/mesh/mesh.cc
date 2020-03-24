#include "mesh.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                                   TETRAHEDRON METHODS                                //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Tetrahedron::Tetrahedron(){
}

/**
 *Constructor receiving 4 nodes
 *@param[in] nod1 Is the first node of tetrahedron
 *@param[in] nod2 Is the second node of tetrahedron
 *@param[in] nod3 Is the third node of tetrahedron
 *@param[in] nod4 Is the fourth node of tetrahedron
 */
Tetrahedron::Tetrahedron( int nod1 , int nod2 , int nod3 , int nod4 ){
	nodes_[ 0 ] = nod1;
	nodes_[ 1 ] = nod2;
	nodes_[ 2 ] = nod3;
	nodes_[ 3 ] = nod4;
	vol_ = 0.0;
}

/**
 *Defalut destructor
 */
Tetrahedron::~Tetrahedron(){

}

//GETS
int Tetrahedron::GetNodeIndex( int node ){
	assert( node >= 0);
	assert( node < 4 );

	return nodes_[ node ];
}

double Tetrahedron::GetVolume(){
	return vol_;
}

//SETS

//UTILITIES
void Tetrahedron::CalculateTetraedronVolume( double* N0 , double* N1 , double* N2 , double* N3 ){
	double A[ 3 ];
	double B[ 3 ];
	double C[ 3 ];
	CalcVector( A , N0 , N1 );
	CalcVector( B , N0 , N2 );
	CalcVector( C , N0 , N3 );
	vol_ = 0.166666666 * Det3( A , B , C );
}


//////////////////////////////////////////////////////////////////////////////////////////
//                                      MESH METHODS                                    //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 */
Mesh::Mesh(){
	scale_ = NULL;
}

/**
 *Constructor receiving a scaler class
 *@param[in] scale Is a scaler class instanced
 */
Mesh::Mesh( Scaler* scale , int n_nodes ){
	scale_ = scale;
	nodes_.assign( n_nodes , NULL );
}

/** 
 *Default destructor
 */
Mesh::~Mesh(){
	if( scale_ ){
		delete scale_;
		scale_ = NULL;
	}
	elems_.clear();
	nodes_.clear();
}



//GETS
/**
 *Getting number of nodes on mesh
 *@return A size_t value indicating the number of nodes on the mesh
 */
size_t Mesh::GetNNodes(){
	return nodes_.size();
}

/**
 *Getting number of elements on mesh
 *@return A size_t value indicating the amount of elements on mesh
 */
size_t Mesh::GetNElements(){
	return elems_.size();
}

/**
 *Getting element from mesh
 *@param[in] i_elem Is the index element requested
 *@return A tetrahedron pointer to the element requested
 */
Tetrahedron* Mesh::GetElement( int i_elem ){
	return elems_[ i_elem ];
}

/**
 *Getting node from mesh
 *@param[in] i_node Is the index of node requested
 *@return A node pointer to the node requested
 */
Node* Mesh::GetNode( int i_node ){
	return nodes_[ i_node ];
}

/**
 *Getting coordinate from node
 *@param[in] node Is the node index 
 *@param[in] dim Is the dimention where the coordinate is requested 0-X, 1-Y and 2-Z
 *@return A double value with the coordinate requested
 */
double Mesh::GetCoord( int node , int dim ){
	return nodes_[ node ]->GetCoord( dim );
}

/**
 *Getting index from node
 *@param[in] i_eleme Is the element index on the list
 *@param[in] i_node Is the node position on element to obtain its index
 *@return A size_t value with the index requested
 */
size_t Mesh::GetIndex( int i_elem , int i_node ){
	return elems_[ i_elem ]->GetNodeIndex( i_node );
}

/**
 *Getting index related to node color
 *Sign -1 means inside node, its index is 3 in order to be printed as green
 *Sign 0  means bounday node, its index is 4 in order to be printed as blue
 *Sign 1 means outside node, its index is 2 in order to be printed as red node
 *Sign 5 means no setted node, its index is 1 in order to be printed as black node
 */
int Mesh::GetIndexRelatedToNodeColor( size_t node_index ){
	assert( node_index >= 0 );	
	assert( node_index < this->GetNNodes() );
	char sign = nodes_[ node_index ]->GetSign();
	int value;
	switch( sign ){
		case -1:
			value = 3;
			break;
		case 0:
			value = 4;
			break;
		case 1:
			value = 2;
			break;
		case -5:
			value = 1;
			break;
		default:
			assert( 0 );
			break;
	}
	return value;
}


//SETS
/**
 *Setting element on mesh
 *@param[in]  elem Is the pointer to the element that will be set
 */
void Mesh::SetElement( Tetrahedron* elem ){
	elems_.push_back( elem );
}

/**
 *Setting node on mesh
 *@param[in] node is the pointer to the node that will be set
 */
void Mesh::SetNode( Node* node ){
	nodes_.push_back( node );
}

/**
 *Setting node on a position of the list
 *@param[in] node Is the pointer to the node that will be setted
 *@param[in] index Is the position where the node will be setted
 */
void Mesh::SetNode( Node* node , int index ){
	nodes_[ index ] = node;
}

//UTILITIES
/**
 *Testing if a mesh node is setted or not
 *@param[in] node Is the node index to be testet to know if has been setted or not
 *@return A bool value indicating if node is setted or not
 */
bool Mesh::IsNodeSetted( int node ){
	return nodes_[ node ];
}

/**
 *Adding a Node on tetrahedral mesh
 *@param[in] index Is the index of the node to be added
 *@param[in] coord Are the coordinates of the node to be added on the mesh
 */
void Mesh::AddNode( int index , double* coord , bool interfaz , char sign ){
	Node* node = new Node( coord , index );
	//node->SetInterface( interfaz );
	node->SetSign( sign );
	this->SetNode( node , index );
}

/**
 *Adding element to the mesh
 @param[in] index0 Is the first node of the tetrahedron to be set
 @param[in] index1 Is the first node of the tetrahedron to be set
 @param[in] index2 Is the first node of the tetrahedron to be set
 @param[in] index3 Is the first node of the tetrahedron to be set
 */
void Mesh::AddElement( int index0 , int index1 , int index2 , int index3 ){
	Tetrahedron* tet = new Tetrahedron( index0 , index1 , index2 , index3 );
	this->SetElement( tet );
}


/**
 *Getting all nodes and elements on a body fitted mesh, in this method also are actualized 
 *the indexes of the elements and nodes belonging to the embeded mesh
 *@param[in] nodes_index Is the array where the node indexes will be actualized
 *@param[in] elems_index Is the array where the elements indexes will be actualized
 */
void Mesh::GetAllNodesAndElemsOnBodyFittedMesh(  int* nodes_index , int* elems_index  ){
	size_t n_nodes = this->GetNNodes();
	size_t n_elems = this->GetNElements();
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		bool flag = false;
		for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){
			int index = elems_[ i_elem ]->GetNodeIndex( i_node );
			int sign = nodes_[ index ]->GetSign();
			if(  sign == -1  ){
				flag = true;
				break;
			}
		}

		if(  flag  ){
			nodes_index[ elems_[ i_elem ]->GetNodeIndex( 0 ) ] = 1;
			nodes_index[ elems_[ i_elem ]->GetNodeIndex( 1 ) ] = 1;
			nodes_index[ elems_[ i_elem ]->GetNodeIndex( 2 ) ] = 1;
			nodes_index[ elems_[ i_elem ]->GetNodeIndex( 3 ) ] = 1;
			elems_index[ i_elem ] = 1;
		}
	}
	int cont = 1;
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if(  nodes_index[ i_node ] == 1  ){
			nodes_index[ i_node ] = cont;
			cont++;
		}
	}
	cont = 1;
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		if(  elems_index[ i_elem ] == 1  ){
			elems_index[ i_elem ] = cont;
			cont++;
		}
	}
}

/**
 *Pproject all exterior and not setted nodes to the triangles defining boundary
 *This function projects all outer and not setted nodes to the triangles defining boundary
 *@param[in] root Is the cell containing the octree
 *@param[in] nodes_index Is the list of nodes belonging to body fitted mesh, the value is bigger 
 *                       than zero is node belongs to boddy fitted mesh and it indicates the real
 *											 nodal index
 *@param[in] r_level Is the refinement level of root
 */
void Mesh::ProjectNodesToboundary( OctreeCell* root , int* nodes_index , int r_level ){
	key_type root_max_key[ 3 ];
	root->GetKey( 6 , root_max_key );
	size_t n_nodes = this->GetNNodes();
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if( nodes_index[ i_node ] > 0 ){//Node belongs to body fitted mesh
			if(  !(nodes_[ i_node ]->GetSign() == -1)  ){//If node is outside or not setted then is projected
				double coord[ 3 ];
				coord[ 0 ] = nodes_[ i_node ]->GetCoord( 0 );
				coord[ 1 ] = nodes_[ i_node ]->GetCoord( 1 );
				coord[ 2 ] = nodes_[ i_node ]->GetCoord( 2 );
				key_type point_key[ 3 ];
        point_key[ 0 ] = static_cast<key_type> ( (1 << ROOT_LEVEL) * coord[ 0 ] );
        point_key[ 1 ] = static_cast<key_type> ( (1 << ROOT_LEVEL) * coord[ 1 ] );
        point_key[ 2 ] = static_cast<key_type> ( (1 << ROOT_LEVEL) * coord[ 2 ] );


				int kind = GetKindOfCaseOnLocalRoot( root_max_key , point_key );
				PerturbateKeysDependingOnKind( point_key , kind );

				OctreeCell* cell = root;
				for(  size_t i = 0  ;  i < ROOT_LEVEL  ;  i++  ){
					if(  cell->IsLeaf()  ){
						break;
					}
					cell = cell->pGetChild( point_key[ 0 ] , point_key[ 1 ] , point_key[ 2 ] );
				}
				double min_corner[ 3 ];
				cell->GetMinPointNormalized( min_corner );

				if( !cell->GetNObjects() ){
					bool local_flag = false;
					key_type neigh_key[ 3 ];
					OctreeCell* neigh_cell;
					for(  int i_direction = 0  ;  i_direction < 18  ;  i_direction++  ){
						if(  cell->GetNeighbourKey( i_direction , neigh_key )  ){
							neigh_cell = root;
							for(  size_t i = 0  ;  i < ROOT_LEVEL  ;  i++  ){
								if(  neigh_cell->IsLeaf()  ){
									break;
								}
								neigh_cell = neigh_cell->pGetChild( neigh_key[ 0 ] , neigh_key[ 1 ] , neigh_key[ 2 ] );
							}
							if( neigh_cell->GetNObjects() ){
								local_flag = neigh_cell->ProjectNodeToCloserTriangle( coord );
							}
							if(  local_flag  ){
								break;
							}
						}
					}
					if(  !local_flag  ){
						double min_dist = 1e20;
						double min_point[ 3 ] = {0.0,0.0,0.0};
						for(  int i_direction = 0  ;  i_direction < 18  ;  i_direction++  ){
							if(  cell->GetNeighbourKey( i_direction , neigh_key )  ){
								neigh_cell = root;
								for(  size_t i = 0  ;  i < ROOT_LEVEL  ;  i++  ){
									if(  neigh_cell->IsLeaf()  ){
										break;
									}
									neigh_cell = neigh_cell->pGetChild( neigh_key[ 0 ] , neigh_key[ 1 ] , neigh_key[ 2 ] );
								}
								if( !( neigh_cell->GetNObjects() == 0 ) ){
									double copy_coord[ 3 ];
									copy_coord[ 0 ] = coord[ 0 ];
									copy_coord[ 1 ] = coord[ 1 ];
									copy_coord[ 2 ] = coord[ 2 ];
									neigh_cell->ProjectNodeToCloserSide( copy_coord );
									double dist = pow( ( copy_coord[ 0 ] - coord[ 0 ] )*( copy_coord[ 0 ] - coord[ 0 ] ) +
								     								 ( copy_coord[ 1 ] - coord[ 1 ] )*( copy_coord[ 1 ] - coord[ 1 ] ) +
								     								 ( copy_coord[ 2 ] - coord[ 2 ] )*( copy_coord[ 2 ] - coord[ 2 ] ) , 0.5 );
									if(  dist < min_dist  ){
										min_dist = dist;
										min_point[ 0 ] = copy_coord[ 0 ];
										min_point[ 1 ] = copy_coord[ 1 ];
										min_point[ 2 ] = copy_coord[ 2 ];
									}
								}
							}
						}
						coord[ 0 ] = min_point[ 0 ];
						coord[ 1 ] = min_point[ 1 ];
						coord[ 2 ] = min_point[ 2 ];
						local_flag = true;
					}
					assert( local_flag );
					nodes_[ i_node ]->SetCoord( 0 , coord[ 0 ] );
					nodes_[ i_node ]->SetCoord( 1 , coord[ 1 ] );
					nodes_[ i_node ]->SetCoord( 2 , coord[ 2 ] );
					nodes_[ i_node ]->SetSign( 0 );
				}else{
					bool local_flag = false;
					local_flag = cell->ProjectNodeToCloserTriangle( coord );
					if(  !local_flag  ){
						local_flag = cell->ProjectNodeToCloserSide( coord );
					}
					assert( local_flag );
					nodes_[ i_node ]->SetCoord( 0 , coord[ 0 ] );
					nodes_[ i_node ]->SetCoord( 1 , coord[ 1 ] );
					nodes_[ i_node ]->SetCoord( 2 , coord[ 2 ] );
					nodes_[ i_node ]->SetSign( 0 );
				}
			}
		}
	}
}

/**
 *This function calculates the quality of tetrahedra on boundary, it means at less one node inside 
 *and other outside or not setted
 */
void Mesh::CalculateBoundaryTetrahedraVolume(  ){
	size_t n_elements = this->GetNElements();
	for(  size_t i_elem = 0  ;  i_elem < n_elements  ;  i_elem++  ){
		int index[ 4 ];
		index[ 0 ] = elems_[ i_elem ]->GetNodeIndex( 0 );
		index[ 1 ] = elems_[ i_elem ]->GetNodeIndex( 1 );
		index[ 2 ] = elems_[ i_elem ]->GetNodeIndex( 2 );
		index[ 3 ] = elems_[ i_elem ]->GetNodeIndex( 3 );
		bool flag = false;
		for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){	
			if(  nodes_[ index[ i_node ] ]->GetSign() == 0  ){
				flag = true;
				break;
			}
		}
		if(  flag  ){
			double N0[ 3 ];
			double N1[ 3 ];
			double N2[ 3 ];
			double N3[ 3 ];
			for(  int i_pos = 0  ;  i_pos < 3  ;  i_pos++  ){
				N0[ i_pos ] = nodes_[ index[ 0 ] ]->GetCoord( i_pos );
				N1[ i_pos ] = nodes_[ index[ 1 ] ]->GetCoord( i_pos );
				N2[ i_pos ] = nodes_[ index[ 2 ] ]->GetCoord( i_pos );
				N3[ i_pos ] = nodes_[ index[ 3 ] ]->GetCoord( i_pos );				
			}
			elems_[ i_elem ]->CalculateTetraedronVolume( N0 , N1 , N2 , N3 ); 
		}
	}
}

/**
 *Obtaining a valid mesh
 *This method erase all elements with negative or almost zero volume
 *@param[in] n_nodes Number of nodes in the whooe octree
 *@param[in] n_elems Number of elements in the whole octree
 *@param[oout] nodes_index Real index of valid nodes
 *@param[oout] elems_index Real index of valid nodes  
 */
void Mesh::EraseNotValidTetrahedra( size_t n_nodes , size_t n_elems , int* nodes_index , int* elems_index ){
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		nodes_index[ i_node ] = 0;
	}
	int acumulado = 0;
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		if(  elems_index[ i_elem ] > 0  ){
			int index[ 4 ];
			index[ 0 ] = elems_[ i_elem ]->GetNodeIndex( 0 );
			index[ 1 ] = elems_[ i_elem ]->GetNodeIndex( 1 );
			index[ 2 ] = elems_[ i_elem ]->GetNodeIndex( 2 );
			index[ 3 ] = elems_[ i_elem ]->GetNodeIndex( 3 );
			bool flag = false;
			for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){	
				if(  nodes_[ index[ i_node ] ]->GetSign() == 0  ){
					flag = true;
					break;
				}
			}
			if(  flag  ){
				if(  elems_[ i_elem ]->GetVolume() <= 1e-11  ){
					elems_index[ i_elem ] = 0;
					acumulado++;
				}else{
					nodes_index[ elems_[ i_elem ]->GetNodeIndex( 0 ) ] = 1;
					nodes_index[ elems_[ i_elem ]->GetNodeIndex( 1 ) ] = 1;
					nodes_index[ elems_[ i_elem ]->GetNodeIndex( 2 ) ] = 1;
					nodes_index[ elems_[ i_elem ]->GetNodeIndex( 3 ) ] = 1;
					elems_index[ i_elem ] -= acumulado;
				}
			}else{
				nodes_index[ elems_[ i_elem ]->GetNodeIndex( 0 ) ] = 1;
				nodes_index[ elems_[ i_elem ]->GetNodeIndex( 1 ) ] = 1;
				nodes_index[ elems_[ i_elem ]->GetNodeIndex( 2 ) ] = 1;
				nodes_index[ elems_[ i_elem ]->GetNodeIndex( 3 ) ] = 1;
				elems_index[ i_elem ] -= acumulado;
			}
		}
	}
	int cont = 1;
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if(  nodes_index[ i_node ] == 1  ){
			nodes_index[ i_node ] = cont;
			cont++;
		}
	}
}




//SAVING MESH ON FILE
/**
 *Saving body fittedmesh on GiD File format
 *@param[in] mpi_rank Is the process that will save the mesh
 */
void Mesh::SaveBodyFittedMeshOnGiDFile(  char* output_ , size_t n_nodes , size_t n_elems , 
																				 int* nodes_index , int* elems_index ){
	FILE* fp = NULL;
	char name[1000];
	strcpy( name , "tet_body_" );
	strcat( name , output_ );

	fp = fopen( name , "w" );	
	if(  !fp  ){
		std::cout << " Error, output file not opened " << std::endl;
		std::cout << name << std::endl;
		assert( fp );
	}	

	fprintf( fp , "MESH \"tet\" dimension 3 ElemType Tetrahedra Nnode 4\n" );
	fprintf( fp , "# color 96 96 96\n" );
	fprintf( fp , "Coordinates\n" );
	fprintf( fp , "# node number coordinate_x coordinate_y coordinate_z  \n" );
	size_t counter = 0;
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if( nodes_index[ i_node ] > 0 ){
			fprintf( fp , "%d ", nodes_index[ i_node ] );
			for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
				double coord = this->GetCoord( (int)i_node , i_dim );
				fprintf( fp , "%lf ", scale_->UnscaleCoord( i_dim , coord ) );
			}
			fprintf( fp, "\n" ); 
			counter++;
		}
	}
	std::cout << "         Nodes: " << counter << std::endl;
	
	counter = 0;
	fprintf( fp , "end coordinates\n" );
	fprintf( fp , "Elements\n" );
	fprintf( fp , "# element node_1 node_2 node_3 node_4 material_number\n" );
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		if(  elems_index[ i_elem ] > 0  ){
			fprintf( fp , "%d ", elems_index[ i_elem ] );
			for(  int i_node = 0 ;  i_node < 4  ;  i_node++  ){
				size_t index = this->GetIndex( i_elem , i_node );
				fprintf( fp , "%d ", nodes_index[ index ] );
			}
			fprintf( fp , "%d\n", 100 );
			counter++;
		}
	}
	fprintf( fp , "end elements\n" );
	fclose( fp );
	std::cout << "         Elements: " << counter << std::endl;
}

/**
 *Saving nodes and its sign on a txt file
 *@paran[in] name Is the name of file to save the nodal information
 *@param[in] n_nodes Number of nodes in the whole octree
 *@param[in] n_elems Number of tetrahedra in whole octree
 *@param[in] nodes_index Number of node in the body fitted mesh
 *@param[in] elems_index Number of tetrahedra in the body fitted mesh 
 */
void Mesh::SaveNodesSignOnTxtFile( char* name , size_t n_nodes , size_t n_elems , int* nodes_index , 
																 int* elems_index ){
	FILE* fp = NULL;
	fp = fopen( name , "w" );	
	if(  !fp  ){
		std::cout << " Error, output file not opened " << std::endl;
		std::cout << name << std::endl;
		assert( fp );
	}	
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if( nodes_index[ i_node ] > 0 ){
			fprintf( fp , "%d ", nodes_index[ i_node ] );
			for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
				double coord = this->GetCoord( (int)i_node , i_dim );
				fprintf( fp , "%lf ", scale_->UnscaleCoord( i_dim , coord ) );
			}
			fprintf( fp, "%d\n" , this->GetIndexRelatedToNodeColor( i_node ) ); 
		}
	}
	

	fclose( fp ); 
}

//DEBUG
void Mesh::PrintInfo(){
	size_t n_nodes = this->GetNNodes();
	size_t n_elems = this->GetNElements();
	std::cout<<" Nodes = "<<n_nodes<<std::endl;
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		if(  nodes_[ i_node ] != NULL  ){
			std::cout<<i_node<<"-"<<this->GetCoord( i_node , 0 )<<" "<<this->GetCoord( i_node , 1 )<<" "<<this->GetCoord( i_node , 2 )<<std::endl;
		}else{
			std::cout<<i_node<<" Nodes not setted"<<std::endl;
			//break;
		}
	}
	std::cout<<" Elements = "<<n_elems<<std::endl;
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		std::cout<<i_elem<<"-"<<this->GetIndex( i_elem , 0 )<<" "<<this->GetIndex( i_elem , 1 )<<" "<<this->GetIndex( i_elem , 2 )<<" "<<this->GetIndex( i_elem , 3 )<<std::endl;
	}
	scale_->PrintScaler();
}





















