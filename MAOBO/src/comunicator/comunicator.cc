#include "comunicator.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                                   COMUNICATOR METHODS                                //
//CONSTRUCTOR AND DESTRUCTOR
/**
 *Default constructor
 *Initializes all variables in 0 or NULL 
 */
Mesher::Mesher(){
	nodes_label_ = 0;
	r_level_ = 0;
	octree_ =  NULL;
	local_root_ = NULL;
	bnd_ = NULL;
	msh_ = NULL;
}

/**
 *Constructor receiving program arguments
 *Read the boundary meshe and creates the octree
 *@param[in] argc Is the number inidicating the program arguments
 *@param[in] argv Is the list of arguments for the program
 */
Mesher::Mesher( int argc , char** argv ){
	input_ = argv[ 1 ];
	output_ = argv[ 2 ];
	//Instantiating a boundary class and reading boundary mesh
	bnd_ = new Boundary( input_ );
	r_level_ = atoi( argv[ 3 ] );
	octree_ = new OctreeDriver();
	local_root_ = octree_->GetRoot();
	srand( time(NULL) );
	nodes_label_ = 0;
}

/**
 *Default destructor 
 */
Mesher::~Mesher(){

}

//GETS
/**
 *Getting all leaves on local root
 *@param[out] leaves Is the vector that will contain all leaves on local root
 */
void Mesher::GetAllLocalRootLeaves( OctreeCell_vector& leaves ){
	OctreeCell_vector cells_stack;
	cells_stack.push_back( local_root_ );
	while( !cells_stack.empty() ){
		OctreeCell* cell = cells_stack.back();
		cells_stack.pop_back();
		if( cell->HasChildren() ){
			for(  size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
				cells_stack.push_back( cell->pGetChild( i_child ) );
			}
		}else{
			leaves.push_back( cell );
		}
	}
}

/**
 *Getting cell on local octree using key
 *This method takes a key and returns a pointer to the cell that contains the key
 *@param[in] keys Are the keys to look for the cell
 *@return A OctreeCell pointer to the cell containning keys
 */
OctreeCell* Mesher::GetCellOnLocalOctree( key_type* keys ){
	OctreeCell* cell = local_root_;
	for(  size_t i = 0  ;  i < ROOT_LEVEL  ;  i++  ){
		if(  cell->IsLeaf()  ){
			return cell;
		}
		cell = cell->pGetChild( keys[ 0 ] , keys[ 1 ] , keys[ 2 ] );
	}
	return cell;
}

/**
 *Look for the minimum refinement on the local octree and returning a integer
 *This method obtains the minimum refinement level used in the local root and substracts
 *1 in order to can go over all octree nodes because there are not linear nodes
 *
 *		X______X______X
 *		|							| 
 *		|							|
 *		|							|
 *		X							X
 *		|							|
 *		|							|
 *		|_____________|
 *    X      X      X
 */
int Mesher::GetMinRefinementLevelUsed(){
	//Obtiene el nivel minimo de refinamiento existente en el octree y le 
	//resta 1, ya que es el valor que se requiere para poder recorrer todos los
	//nodos existentes en el octree
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );
	size_t n_leaves = leaves.size();
	unsigned char min_refinement = 30;
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		if(  leaves[ i_leaf ]->GetLevel() < min_refinement  ){
			min_refinement = leaves[ i_leaf ]->GetLevel();
		}
	}
	return ( min_refinement - 1 );
}

/**
 *Getting min_key from local root
 *@param[out] keys Is the variable where the min key will be saved
 */
void Mesher::GetMinKeyRoot( key_type* keys ){
	local_root_->GetMinKey( keys[ 0 ] , keys[ 1 ] , keys[ 2 ] );
}

/**
 *Getting max key from local root
 *@param[out] keys Is the variable where the max key will be saved
 */
void Mesher::GetMaxKeyRoot( key_type* keys ){
	local_root_->GetKey( 6 , keys );
}

/**
 *Getting number of nodes in the tetrahedral mesh
 *@return A size_t value indicating the number of nodes on tetrahedral mesh
 */
size_t Mesher::GetNumberOfTetrahedralNodes(){
	return msh_->GetNNodes();
}

/**
 *Getting number of elements in the tetrahedral mesh
 *@return A size_t value indicating the number of tetrahedra on mesh
 */
size_t Mesher::GetNumberOfTetrahedraOnMesh(){
	return msh_->GetNElements();
}


//SETS
/**
 *Setting nodes position respect top boundary
 *This method sets the nodes position respect to boundary on cells that intersects the 
 *boundary domain, in this case is used the coloring algorithm using the octree structure
 */
void Mesher::SetNodesPositionRespectToBoundary(){

	//Min refinement on local octree
	int min_refinement = this->GetMinRefinementLevelUsed();

	std::cout<<"         Refinamiento: "<<min_refinement<<std::endl;
	//Number of rays to trace on each dimension
	int n_rays = ( (1<<(ROOT_LEVEL)) / (1<<min_refinement) )+1;

	//Getting min and max keys from root
	key_type root_min_key[ 3 ];
	this->GetMinKeyRoot( root_min_key );
	key_type root_max_key[ 3 ];
	this->GetMaxKeyRoot( root_max_key );

	double t1;
	//performing ray casting over three directions
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		t1 = omp_get_wtime();
		this->RayCastingOverDirection( root_min_key , root_max_key , min_refinement , n_rays , i_dim );
		fprintf(stdout, "                   [%16lf]   **Ray casting algorithm    : \n", omp_get_wtime() - t1 );
		t1 = omp_get_wtime();
	}

	this->CheckNodesSign();
	fprintf(stdout, "                   [%16lf]   **Checking coloring alg    : \n", omp_get_wtime() - t1 );

}

/**
 *Setting node sign on cell and all its local neighbours
 *@param[in] cell Is the cell where the information is set
 *@param[in] node Is the local position where the sign will be set
 *@param[in] value Is the node sign that will be set
 *@param[in] dim Is the axis where the sign is set
 */
void Mesher::SetSignOnNodeAndNeighbours( OctreeCell* cell , int node , int value , int dim ){
	assert( node >= 0 );      
	assert( node <  27);      
  assert( cell->IsLeaf() ); 
	//Checking if node had been added before
	DataInsideOctreeCell *data = cell->pGetData();
	if( data->_IsSetNodeSign( dim , node ) ){
		return ;
	}
	data->_SetNodePositionRespectToBoundary( dim , node , value );
	if(  node == 8  ){//this is the center node, does not have neighbours
		return ;
	}
	//obtainig key of node
	key_type point_key[ 3 ];
	cell->GetKey( (size_t)node , point_key );
	//adding node to neighbours, for each node exist 8 possible cell
	for(  size_t i_direction = 0  ;  i_direction < 8  ;  i_direction++  ){
		key_type neighbour_key[ 3 ];
		if(  cell->GetNeighbourKey( node , i_direction , neighbour_key ) ){
			OctreeCell* neighbour_cell =	this->GetCellOnLocalOctree( neighbour_key );
			if( !neighbour_cell ||  ( neighbour_cell == cell ) ){//if neighbour does not exist or is the actual cell
				continue;                                     //node is not added
			}
			DataInsideOctreeCell *neighbour_data = neighbour_cell->pGetData();
			int position;
			if(  neighbour_cell->ExistNodeOnCell( point_key , &position )  ){
				neighbour_data->_SetNodePositionRespectToBoundary( dim , position , value );
			}
		}
	}
}

/**
 *Setting all nodes on interfaces with other ranks
 *This method set all valid octree nodes, in this case the nodes valid are the ones that 
 *belongs to a cell containing at less one interior or boundary node
 *@param[in] recv Is the list of cells received from other ranks
 */
void Mesher::SetAllOctreeNodes(  ){
	this->SetAllLocalOctreeNodes();
}

void Mesher::SetAllLocalOctreeNodes(){
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );	
	size_t n_leaves = leaves.size();
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		OctreeCell* cell = leaves[ i_leaf ];
		key_type f_neigh_key[ 3 ];
		this->SetCenterNodeOnCell( cell , false );
		//neighbors by face
		for(  int f_direction = 0  ;  f_direction < 6  ;  f_direction++  ){
			if(  cell->GetNeighbourKey( f_direction , 	f_neigh_key )  ){
				OctreeCell* f_neigh_cell = this->GetCellOnLocalOctree( f_neigh_key );
				int current_level = cell->GetLevel();
				int f_neigh_level = f_neigh_cell->GetLevel();
				switch( current_level - f_neigh_level ){
					case -1://neighbour is bigger ading all nodes on that face
						this->SetAllNodesOnCellFace( f_direction , f_neigh_cell , false );
						break;
					case 0://neighbor is same refined checking neighbour by edge
						this->SetAllNodesOnLocalNeighbourSameRefined( f_direction , cell , false  );
						break;
					case 1://neighbour is lower getting all neighbours on that face and adding all linear nodes 
						this->SetAllLinearNodesOnNeighborsCellsByFace( f_direction , cell , false );
						break;
					default:
						assert( 0 );
						break;	
				}
			}else{//Boundary face, add linear and center
				this->SetAllLinearAndCenterNodesOnBoundaryCellFace( f_direction , cell , false );
			}
		}
	}
}

/**
 *Setting center node on cell
 *@param[in] cell Is the cell where center node will be set
 *@param[in] interfaz A bool vale indicating if node is on interfaz between processors
 */
void Mesher::SetCenterNodeOnCell( OctreeCell* cell , bool interfaz ){
	this->SetNodeInCellAndNeighbors( 8 , cell , NULL , interfaz );
}

/**
 *Setting a node on a cell and all its neighbors
 *This method set a node on a cell and all its neighbors, in this case exist at most 8 
 *neighbours
 *@param[in] node Is the node index to be set on the local cell
 *@param[in] cell Is the cell where the node have to be set
 *@param[in] coord Is the coordinate to be set 
 *@param[in] interfaz Is a bool value indicating if a node is interfaz beetwen process
 */
bool Mesher::SetNodeInCellAndNeighbors( int node , OctreeCell* cell , double* coord , bool interfaz ){

	DataInsideOctreeCell* data = cell->pGetData();
	if(  data->_IsSetNode( node )  ){
		return false; 
	}
	this->SetNodeInCell( node , cell , coord , interfaz );
 	if(  node == 8  ){
		nodes_label_++;
		return true;
	}
	key_type point_key[ 3 ];
	cell->GetKey( (size_t)node , point_key );
	for(  size_t i_direction = 0  ;  i_direction < 8  ;  i_direction++  ){
		key_type neighbour_key[ 3 ];
		if(  cell->GetNeighbourKey( (size_t)node , i_direction , neighbour_key )  ){
			OctreeCell* neighbour_cell = this->GetCellOnLocalOctree( neighbour_key );
			if(  !neighbour_cell || ( neighbour_cell == cell )  ){
				continue;
			}
			size_t position = neighbour_cell->GetLocalPosition( point_key );
			this->SetNodeInCell( (int)position , neighbour_cell , coord , interfaz );         
		}
	}
	nodes_label_++;
	return true;
}

/**
 *Setting a node on a cell
 *This method set a node on a cell
 *@param[in] node Is the local position where the node will be set
 *@param[in] cell Is the cell where the node will be set
 *@param[in] coord Is the pointer to the coord that will be set 
 *@param[in] interfaz Is a bool value indicating if a node is interfaz beetwen process
 */
void Mesher::SetNodeInCell( int node , OctreeCell* cell , double* coord , bool interfaz ){
	DataInsideOctreeCell* data = cell->pGetData();
	data->_SetNodeIndex( node , nodes_label_ );
}

/**
 *Setting All nodes on boundary cell face
 *@param[in] direction Is the direction where the face is located
 *@param[in] cell Is the cell where nodes will be set
 *@param[in] interfaz A bool value indicating if node is or not on interfaz between process 
 */
void Mesher::SetAllLinearAndCenterNodesOnBoundaryCellFace(  int direction , OctreeCell* cell , bool interfaz ){
	int index[ 5 ];
	GetLinearAndCenterIndexFromBoundaryFace( direction , index );
	for(  int i_node = 0  ;  i_node < 5  ;  i_node++  ){
		this->SetNodeInCellAndNeighbors( index[ i_node ] , cell , NULL , interfaz );
	}
}

/**
 *Setting all nodes on a cell face
 *@param[in] src_dir Is the source direction from where the nodes need to be set
 *@param[in] cell Is the cell where the nodes need to be add
 *@param[in] interfaz Is a bool value indicating if nodes to be set are interfaz or not
 */
void Mesher::SetAllNodesOnCellFace( int src_dir , OctreeCell* cell , bool interfaz ){
	this->SetCenterNodeOnCellFace( src_dir , cell , interfaz );
	this->SetAllLinearNodesOnCellFace( src_dir , cell , interfaz );
	this->SetAllQuadraticNodesOnCellFace( src_dir , cell , interfaz );
}

/**
 *Setting the center node of a face on a cell and all its neighbors
 *@param[in] src_dir Is the source direction from where the nodes need to be set
 *@param[in] cell Is the cell where the node will be set
 *@param[in] interfaz Is a bool value indicating if a node is interfaz beetwen process
 */
void Mesher::SetCenterNodeOnCellFace( int src_dir , OctreeCell* cell , bool interfaz ){
	int index = GetCenterIndexFromFace( src_dir );
	this->SetNodeInCellAndNeighbors( index , cell , NULL , interfaz );
}

/**
 *Setting all linear nodes on a cell face
 *@param[in] src_dir Is the source direction from where the nodes need to be set
 *@param[in] cell Is the cell where the nodes need to be add
 *@param[in] interfaz Is a bool value indicating if a node is interfaz beetwen process
 */
void Mesher::SetAllLinearNodesOnCellFace( int src_dir , OctreeCell* cell , bool interfaz ){
	int index[ 4 ];
	GetLinearIndexFromFace( src_dir , index );

	for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){
		this->SetNodeInCellAndNeighbors( index[ i_node ] , cell , NULL , interfaz );
	}
}

/**
 *Setting all quadratic nodes on a cell face
 *@param[in] src_dir Is the source direction from where the nodes need to be set
 *@param[in] cell Is the cell where the nodes need to be add
 *@param[in] interfaz Is a bool value indicating if a node is interfaz beetwen process
 */
void Mesher::SetAllQuadraticNodesOnCellFace( int src_dir , OctreeCell* cell , bool interfaz ){
	int index[ 4 ];
	GetQuadraticIndexFromFace( src_dir , index );
	for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){
		this->SetNodeInCellAndNeighbors( index[ i_node ] , cell , NULL , interfaz );
	}
}

/**
 *Setting local nodes when neighbours by face are same refined, then are checked 
 *neighbours by edge and needed extra nodes are added
 *@param[in] direction Is the direction of neighbour
 *@param[in] src_cell Is the cell used to search the neighbours
 *@param[in] interfaz Is a bool value indicating is nodes are interfaz between process
 */
void Mesher::SetAllNodesOnLocalNeighbourSameRefined( int direction , OctreeCell* src_cell , bool interfaz  ){
	int f_index[ 4 ];
	GetLinearIndexOnFace( direction , f_index );
	for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){
		this->SetNodeInCellAndNeighbors( f_index[ i_node ] , src_cell , NULL , interfaz );
	}
	int e_index[ 4 ];

	GetEdgesIndexesFromFace( direction , e_index );
	for(  int i_edge = 0  ;  i_edge < 4  ;  i_edge++  ){
		key_type e_neigh_key[ 3 ];
		if(  src_cell->GetNeighbourKey( (size_t)e_index[ i_edge ] , e_neigh_key )  ){
			OctreeCell* e_cell = this->GetCellOnLocalOctree( e_neigh_key );				
			int current_level = src_cell->GetLevel();
			int e_neigh_level = e_cell->GetLevel();
			if(  e_neigh_level < current_level  ){
				int index = GetCenterIndexOnEdge( e_index[ i_edge ] );
				this->SetNodeInCellAndNeighbors( index , src_cell , NULL , interfaz );
			}
		}
	}
}

/**
 *Setting all linear nodes on 4 neighbours of a cell 
 *@param[in] direction Is the direction where neighbors need to bee searched
 *@param[in] src_cell Is the cell used to search for neighbours
 *@param[in] interfaz Is a bool value indicating if nodes to be set are interfaz or not
 */
void Mesher::SetAllLinearNodesOnNeighborsCellsByFace( int direction , OctreeCell* src_cell , bool interfaz ){
	//Getting 4 neighbor cells
  OctreeCell* cell[ 4 ];
	key_type keys0[ 3 ] , keys1[ 3 ] , keys2[ 3 ] , keys3[ 3 ];
	src_cell->GetNeighbourKey( direction , 	keys0 );
	GetAllNeighbourKeysFromFaceLessRefined( direction , keys0 , keys1 , keys2 , keys3 , src_cell->GetLevel() );	
	cell[ 0 ] = this->GetCellOnLocalOctree( keys0 );
	cell[ 1 ] = this->GetCellOnLocalOctree( keys1 );
	cell[ 2 ] = this->GetCellOnLocalOctree( keys2 );
	cell[ 3 ] = this->GetCellOnLocalOctree( keys3 );
	for(  int i_cell = 0  ;  i_cell < 4  ;  i_cell++  ){
		this->SetAllLinearNodesOnCellFace( direction , cell[ i_cell ] , interfaz );		
	}
}




//UTILITIES
/**
 *Scaling boundary nodes
 *This method calls to the Boundary::ScaleNodes method from Boundary class
 */
void Mesher::ScaleBoundaryNodes(){
	bnd_->ScaleNodes();
}

/**
 *Obtaining AABB for mesh elements
 */
void Mesher::CalculateAABBForBoundaryTriangles(){
	bnd_->CalculateAABBForTriangles();
}

/**
 *Refining local cells interseted by boundary
 *This method subdivide the leaves of the local root until reach the r_level_, in each 
 *cell are stored the triangles intersected
 */
void Mesher::RefineLocalRoot(){


	//Setting triangles intersected on local_root cell
	this->SetTrianglesOnLocalRoot();

	//Refining intersected cells
	double t1  = omp_get_wtime();
	this->SubdivideIntersectedCells();
	fprintf(stdout, "                   [%16lf]   **Intersected cells        : \n", omp_get_wtime() - t1 );
}

/**
 *This method set the triangles that intersects on the local root
 *NOTE: The local root is not subdivided yet
 */
void Mesher::SetTrianglesOnLocalRoot( ){
	double min_corner[ 3 ], max_corner[ 3 ];
	double tolerance = 0.001 * double( 1 << MIN_LEVEL ) / double( 1 << ROOT_LEVEL );//0.1% of minimum size
	local_root_->GetMinPointNormalized( min_corner ); 
	local_root_->GetMaxPointNormalized( max_corner );
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		min_corner[ i_dim ] -= tolerance;
		max_corner[ i_dim ] += tolerance;
	}
	size_t n_elements = bnd_->GetNElements();
	for(  size_t i_tri = 0  ;  i_tri < n_elements  ;  i_tri++  ){
		Triangle* tri = bnd_->GetElement( (int)i_tri );
		if(  tri->IntersectsBox( min_corner , max_corner )  ){
			local_root_->Insert( tri );
		}
	}
}

/**
 *Refining local octree based on r_level and  
 *This method refines the octree and transfers the objects on parent cell to childs
 */
void Mesher::SubdivideIntersectedCells(){
	int i_level = 0;
	while( i_level < r_level_ ){
		i_level++;
		OctreeCell_vector leaves;
		this->GetAllLocalRootLeaves( leaves );
		size_t n_leaves = leaves.size();
		for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
			OctreeCell* leaf = leaves[ i_leaf ];
			if( leaf->GetNObjects() ){//Cell is divided in case of intersect a boundary triangle
				if(  leaf->SubdivideCell(  )  ){//Subdivide cell
					std::cout << " Error while subdividing cell " << std::endl;
					assert( 0 ); 
				}
				if(  leaf->TransferObjectsToSons()  ){//transfer triangles to children
					std::cout << " Error while transfering objects to sons " << std::endl;
					assert( 0 ); 
				}
			}

		}
		leaves.clear();
	}
}

/**
 *Balancing octree
 *Make constrain 2:1 be valid (For each octree leaf, all its neighbours must have at more one 	
 *refinement level of difference)
 */
void Mesher::BalanceOctree(  ){

	//Balancing local octree
	double t1 = omp_get_wtime();
	this->BalanceLocalOctree();

	//Allocating neighbour data on cells
	this->AllocateAllLeavesData();
	fprintf(stdout, "                   [%16lf]   **Balancing octree   : \n", omp_get_wtime() - t1 );
}

/**
 *Balancing local octree
 */
void Mesher::BalanceLocalOctree(){

	//getting leaves on local root
	OctreeCell_vector leaves;
	OctreeCell_vector next_leaves;
	this->GetAllLocalRootLeaves( leaves );
	for(  char i_level = MIN_LEVEL  ;  i_level < ROOT_LEVEL - 1  ;  i_level++  ){
		for(  size_t i_cell = 0  ;  i_cell < leaves.size()  ;  i_cell++  ){
			OctreeCell* current_leaf = leaves[ i_cell ];
			if(  current_leaf->GetLevel() == i_level  ){
				key_type neighbour_key[ 3 ];
				//18 is the number of neighbours counting faces and edges of the cell
				for(  int i_direction = 0  ;  i_direction < 18  ;  i_direction++  ){
					if(  current_leaf->GetNeighbourKey( i_direction , neighbour_key )  ){
						OctreeCell* neighbour_cell = this->GetCellOnLocalOctree( neighbour_key );
						if(  neighbour_cell->GetLevel() > i_level + 1  ){
							OctreeCell* temp_neighbour_cell = neighbour_cell;
							for(  char j_level = neighbour_cell->GetLevel()  ;  j_level > i_level + 1  ;  j_level--  ){
								if(  temp_neighbour_cell->SubdivideCell()  ){
									std::cout << " Error while subdividing cell " << std::endl;
									assert( 0 ); 
								}
								if(  temp_neighbour_cell->GetNObjects()  ){//if cell has objects then are transfered to sons
									if(  temp_neighbour_cell->TransferObjectsToSons()  ){
										std::cout << " Error while transfering objects to sons " << std::endl;
										assert( 0 ); 
									}
								}
								//moving to the new cell subdivided
								size_t child_index = temp_neighbour_cell->GetChildIndex( neighbour_key[ 0 ] , neighbour_key[ 1 ] , neighbour_key[ 2 ] );
								for(  size_t j = 0  ;  j < 8  ;  j++  ){
									if(  j != child_index  ){
										next_leaves.push_back( temp_neighbour_cell->GetChildren() + j );
									}
								}
								temp_neighbour_cell = temp_neighbour_cell->GetChildren() + child_index;
								if(  j_level == neighbour_cell->GetLevel() - 1  ) // the last loop we add all the child as leaf
									next_leaves.push_back( temp_neighbour_cell );
							}
						}
					}
				}
			}else{
				if(  current_leaf->IsLeaf()  ){ // becuase it may be divided
					next_leaves.push_back( current_leaf );
				}
			}
		}
		leaves.swap( next_leaves );
		next_leaves.clear();
	}
	if(  !this->CheckLocalConstrain2to1()  ){
		std::cout << " Error when checking constrain 2 to 1 on local octree " << std::endl; 
		assert( 0 );
	}
}

/**
 *Allocating data on leaves
 */
void Mesher::AllocateAllLeavesData(){
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );
	size_t n_leaves = leaves.size();
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		OctreeCell* cell = leaves[ i_leaf ];
		DataInsideOctreeCell* new_data = new DataInsideOctreeCell;
		DataInsideOctreeCell** data_cell = cell->pGetDataPointer();
		(*data_cell) = new_data;
	}
}

/**
 *Tracing cartesian rays over a direction 
 *@param[in] min_refinement Is the minimum refinement-1 on the octree
 *@param[in] n_rays The amount of rays to trace over each dimension
 *@param[in] dim Is the dimention through which the rays will be testes (0-X, 1-Y and 2-Z)
 */
void Mesher::RayCastingOverDirection( key_type* root_min_key , key_type* root_max_key , 
																					 int min_refinement , int n_rays , int dim ){
	int indexes[ 3 ][ 2 ] = { { 1 , 2 } , { 0 , 2 } , { 0 , 1 } };
	#pragma omp parallel for schedule(dynamic)
	for(  int i_ray = 0  ;  i_ray < n_rays  ;  i_ray++  ){
		key_type keys[ 3 ];
		keys[ 0 ] = 0;
		keys[ 1 ] = 0;
		keys[ 2 ] = 0;
		keys[ indexes[ dim ][ 0 ] ] = root_min_key[ indexes[ dim ][ 0 ] ] +  ( i_ray * ( 1<<min_refinement ) );
		keys[ indexes[ dim ][ 1 ] ] = root_min_key[ indexes[ dim ][ 1 ] ];
		for(  int j_ray = 0  ;  j_ray < n_rays  ;  j_ray++  ){
			if(  this->RayIntersectsModelBBox( keys , dim )  ){
				this->RayCasting( root_min_key , root_max_key , keys , dim );
			}
			keys[ indexes[ dim ][ 1 ] ] += ( 1<<min_refinement );
		}
	}
}

/**
 *Fast test to know if ray intersects model bounding box
 *@param[in] keys Are the ray keys, NOTE: key on dim is not tested
 *@param[in] dim Is the axis parallel to the ray
 *@return A bool value indicating if intersection exist or not
 */
bool Mesher::RayIntersectsModelBBox( key_type* keys , int dim ){
	double scale = 1.00 / ( 1 << ROOT_LEVEL );	
	double coord0, coord1;
	switch( dim ){
		case 0:
			coord0 = static_cast<double>(keys[1]) * scale;
			coord1 = static_cast<double>(keys[2]) * scale;
			if(  !( bnd_->IsCoordInsideDomain( coord0 , 1 ) && bnd_->IsCoordInsideDomain( coord1 , 2 ) )  ){
				return false;
			}
			break;
		case 1:
			coord0 = static_cast<double>(keys[0]) * scale;
			coord1 = static_cast<double>(keys[2]) * scale;
			if(  !( bnd_->IsCoordInsideDomain( coord0 , 0 ) && bnd_->IsCoordInsideDomain( coord1 , 2 ) )  ){
				return false;
			}
			break;
		case 2:
			coord0 = static_cast<double>(keys[0]) * scale;
			coord1 = static_cast<double>(keys[1]) * scale;
			if(  !( bnd_->IsCoordInsideDomain( coord0 , 0 ) && bnd_->IsCoordInsideDomain( coord1 , 1 ) )  ){
				return false;
			}
			break;
	}
	return true;
}

/**
 *Tracing a ray over a dimention and applying the coloring algorithm
 *This method trace a ray on dimention dim using the keys and obtains all intersections on
 *boundary and then set the node position respect to the boundary
 *@param[in] keys Is the binary coding of the ray to be traced
 *@param[in] dim Is the dimention through which the ray is tested
 */
void Mesher::RayCasting( key_type* root_min_key , key_type* root_max_key , key_type* keys , int dim ){
	std::vector<double> points;
	points.push_back(-1.00);
	points.push_back( 2.00);
	int n_triangles = (int)bnd_->GetNElements();
	double l_coord[ 3 ], r_coord[ 3 ];
	for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
		l_coord[ i_dim ] =  (double)keys[ i_dim ] / (double)( 1 << ROOT_LEVEL );
    r_coord[ i_dim ] = l_coord[ i_dim ];
	}
	l_coord[ dim ] = -1.00;
	r_coord[ dim ] =  2.00;
	for(  int i_tri = 0  ;  i_tri < n_triangles  ;  i_tri++  ){
		double X[ 2 ] = { 1e20 , 1e20 };
		Triangle* tri = bnd_->GetElement( i_tri );
		switch( tri->GetAllIntersectionsWithLine( l_coord , r_coord , X , dim ) ){

			case 0:
				break;
			case 1:
				if(  !DoesExistCoordOnVector( points , X[ 0 ] )  ){
					points.push_back( X[ 0 ] );
				}
				break;
			case 2:
				if(  !DoesExistCoordOnVector( points , X[ 0 ] )  ){
					points.push_back( X[ 0 ] );
				}

				if(  !DoesExistCoordOnVector( points , X[ 1 ] )  ){
					points.push_back( X[ 1 ] );
				}
				break;
		}
	}

	//Sorting intersections and obtaining the position respect to boundary of each segment
	sort( points.begin() , points.end() );
	int *signs = new int [ points.size()-1 ];
	signs[ 0 ] = 1;
	signs[ points.size()-2 ] = 1;
	for(  size_t i_pos = 1  ;  i_pos < (points.size() - 2)  ;  i_pos++  ){
		signs[ i_pos ] = -signs[ i_pos-1 ]; 
	}
	//Checking for the cells on ray
	key_type aux_keys[ 3 ];
	aux_keys[ 0 ] = keys[ 0 ];
	aux_keys[ 1 ] = keys[ 1 ];
	aux_keys[ 2 ] = keys[ 2 ];
	aux_keys[ dim ] = root_min_key[ dim ];
	while( aux_keys[ dim ] < root_max_key[ dim ] ){
		int kind = GetKindOfCaseOnLocalRoot( root_max_key , aux_keys );
		PerturbateKeysDependingOnKind( aux_keys , kind );
		OctreeCell* cell;

		cell = this->GetCellOnLocalOctree( aux_keys );
		int level = (int)cell->GetLevel();
		UnPerturbateKeysDependingOnKind( aux_keys , kind );
		int position;
		if(  cell->ExistNodeOnCell( aux_keys , &position )  ){
			int add_nodes[ 3 ];
			GetNodesToAddByRayCasting( position , dim , add_nodes );
			for(  int i_node = 0  ;  i_node < 3  ;i_node++  ){
				//ray_index is -1 if node is an intersection
				int ray_index = GetPositionOverRay( points ,( aux_keys[ dim ] + ( i_node * ( 1<<(level-1) ) ) ) );
				if(  ray_index == -1  ){
					this->SetSignOnNodeAndNeighbours( cell , add_nodes[ i_node ] , 0 , dim );
				}else{
					this->SetSignOnNodeAndNeighbours( cell , add_nodes[ i_node ] , signs[ ray_index ] , dim );
				}
			}
		}
		aux_keys[ dim ] +=  ( 1<<level );
	}
	delete[] signs;
}

/**
 *Creating tetrahedral mesh from octree 
 *@param[in] recv Is the list of cells received from other ranks
 */
void Mesher::CreateTetrahedralMeshFromOctree( ){
	msh_ = new Mesh( bnd_->GetScaler() , nodes_label_ );

	this->CreateTetrahedralMeshOnLocalOctree( );
}

/**
 *Creating the tetrahedra pattern on the local octree, in this case are added tetrahedra 
 *based on the tetraheda patterns presented on Robust Volume Mesh Generation for 
 *non-watertight geometries, Dr. Abel Coll
 */
void Mesher::CreateTetrahedralMeshOnLocalOctree( ){
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );
	size_t n_leaves = leaves.size();
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		OctreeCell* cell = leaves[ i_leaf ];
		key_type f_neigh_key[ 3 ];
		for(  int i_direction = 0  ;  i_direction < 6  ;  i_direction++  ){
			if(  cell->GetNeighbourKey( i_direction , f_neigh_key )  ){
				OctreeCell* f_neigh_cell = this->GetCellOnLocalOctree( f_neigh_key );
				int current_level = cell->GetLevel();
				int f_neigh_level = f_neigh_cell->GetLevel();
				switch( current_level - f_neigh_level ){
					case -1://neighbour is bigger, get position on face and and add two tetrahedra
						this->CreateTetrahedraWhenNeighbourIsLessRefined(  cell , f_neigh_cell , i_direction  );
						break;
					case 0://neighbor is same refined, if actual cell is lower, add 1 or 2 tetrahedra on each edge
						this->CreateTetrahedraWhenNeighbourIsSameRefined(  cell , f_neigh_cell , i_direction  );
						break;
					case 1://neighbour is lower add 8 tetrahedra on that face
						this->CreateTetrahedraWhenNeighbourIsMoreRefined(  cell , i_direction  );
						break;
					default:
						assert( 0 ); 
						break;	
				}
			}else{//Neighbour does not exist, are added the tetrahedra on boundary face

				this->CreateTetrahedraOnBoundaryFace( i_direction , cell );
			}
		}
	}
}

/**
 *Creating tetrahedra on face located on boundary of the octree, in this case are trested 
 *tested the quadratic nodes on the direction by face, if a node of the quadratic nodes is
 *setted, then are created two tetrahedra on that edge, otherwise just one is created.
 *@param[in] direction Is the direction where tetrahedra will be created on cell
 *@param[in] cell Is the cell used to apply the tetrahedra pattern
 */
void Mesher::CreateTetrahedraOnBoundaryFace( int direction , OctreeCell* cell ){
	double scale = 1.00 / ( 1 << ROOT_LEVEL );
	int e_indexes[ 4 ];
	GetQuadraticIndexOnFace( direction , e_indexes );
	DataInsideOctreeCell* data = cell->pGetData(  );
	for(  int i_edge = 0 ;  i_edge < 4  ;  i_edge++  ){
		if(  data->_IsSetNode( e_indexes[ i_edge ] )  ){
			int tet_index[ 8 ];
			int real_index[ 8 ];
			GetLocalTetIndexOnBoundaryCellWhenEdgeIsMoreRefined( direction , i_edge , tet_index );
			for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
				real_index[ i_node ] = data->_GetNodeIndex( tet_index[ i_node ] );
				if(  !msh_->IsNodeSetted( real_index[ i_node ] )  ){
					double* coord = new double[ 3 ];
					key_type keys[ 3 ];
					cell->GetKey( tet_index[ i_node ] , keys );
					for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
						coord[ i_dim ] = ( static_cast<double>( keys[ i_dim ] ) * scale );
					}
					msh_->AddNode( real_index[ i_node ] , coord , false , 
												 data->_GetSign( tet_index[ i_node ] , 0 )  );
				}
			}
			//Adding elements to mesh in this case are added 2 elements
			for(  int i_elem = 0  ;  i_elem < 2  ;  i_elem++  ){
				msh_->AddElement( real_index[ i_elem * 4 + 0 ] , real_index[ i_elem * 4 + 1 ] , 

													real_index[ i_elem * 4 + 2 ] , real_index[ i_elem * 4 + 3 ] );
			}
		}else{
			int tet_index[ 4 ];
			int real_index[ 4 ];
			GetLocalTetIndexOnBoundaryCellWhenEdgeIsSameRefined( direction , i_edge , tet_index );
			for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){
				real_index[ i_node ] = data->_GetNodeIndex( tet_index[ i_node ] );
				if(  !msh_->IsNodeSetted( real_index[ i_node ] )  ){
					double* coord = new double[ 3 ];
					key_type keys[ 3 ];
					cell->GetKey( tet_index[ i_node ] , keys );
					for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
						coord[ i_dim ] = ( static_cast<double>( keys[ i_dim ] ) * scale );
					}
					msh_->AddNode( real_index[ i_node ] , coord , 
												 false , 
												 data->_GetSign( tet_index[ i_node ] , 0 )  );
				}
			}
			//Adding elements to mesh in this case are added 2 elements
			msh_->AddElement( real_index[ 0 ] , real_index[ 1 ] ,	real_index[ 2 ] , real_index[ 3 ] );
		}
	}
}

/**
 *In this case the src_cell is more refined and neet to be searched its position on 
 *neighbour cell face
 *@param[in] src_cell Is the actual cell
 *@param[in] neigh_cell Is the neighbour cell less refined
 *@param[in] direction Is the search direction where neigh_cell was founded
 */
void Mesher::CreateTetrahedraWhenNeighbourIsLessRefined(  OctreeCell* src_cell , OctreeCell* neigh_cell , int direction  ){
	double scale = 1.00 / ( 1 << ROOT_LEVEL );
	int position = GetRelativePositionOverNeighbourLessRefined( src_cell , neigh_cell , direction );
	int tet_index[ 8 ];
	GetLocalTetIndexWhenNeighbourIsLessRefined( position , direction , tet_index );
	//Checking if nodes has been added on mesh
	DataInsideOctreeCell* data =  src_cell->pGetData();
	int real_index[ 8 ];
	for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
		real_index[ i_node ] = data->_GetNodeIndex( tet_index[ i_node ] );
		if(  !msh_->IsNodeSetted( real_index[ i_node ] )  ){
			double* coord = new double[ 3 ];
			key_type keys[ 3 ];
			src_cell->GetKey( tet_index[ i_node ] , keys );
			for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
				coord[ i_dim ] = ( static_cast<double>( keys[ i_dim ] ) * scale );
			}
			msh_->AddNode( real_index[ i_node ] , coord , 
										 false , 
										 data->_GetSign( tet_index[ i_node ] , 0 )  );
		}
	}
	//Adding elements to mesh in this case are added 2 elements
	for(  int i_elem = 0  ;  i_elem < 2  ;  i_elem++  ){
		msh_->AddElement( real_index[ i_elem * 4 + 0 ] , real_index[ i_elem * 4 + 1 ] , 
											real_index[ i_elem * 4 + 2 ] , real_index[ i_elem * 4 + 3 ] );
	}
} 

/**
 *Adding tetrahedra when neighbour by face are same refined.
 *In this case the tetrahedra are added if and only if the src_cell is minimum than 
 *neighbour cell, it means that at less one of its min_keys is lower than the other 
 *If src_cell is lower than neighbour cell then need to be checked the quadratic nodes 
 *over edges in order to determine the amount of tetrahedra created over each edge
 *@param[in] src_cell Is the actual cell
 *@param[in] neigh_cell Is the neighbour cell located on the search direction.
 *@param[in] direction Is the neighbour search direction   
 */
void Mesher::CreateTetrahedraWhenNeighbourIsSameRefined(  OctreeCell* src_cell , OctreeCell* neigh_cell , int direction  ){
	double scale = 1.00 / ( 1 << ROOT_LEVEL );
	if(  !SourceCellIsLower( src_cell , neigh_cell )  ){
		return ;
	}
	int e_indexes[ 4 ];
	GetQuadraticIndexOnFace( direction , e_indexes );
	DataInsideOctreeCell* src_data = src_cell->pGetData();
	DataInsideOctreeCell* neigh_data = neigh_cell->pGetData();
	for(  int i_edge = 0 ;  i_edge < 4  ;  i_edge++  ){
		if(  src_data->_IsSetNode( e_indexes[ i_edge ] )  ){
			int tet_index[ 8 ];
			int real_index[ 8 ];
			GetLocalTetIndexOnCellWhenFaceIsSameRefinedAndEdgeMore( direction , i_edge , tet_index );
			for(  int i_node = 0  ;  i_node < 8  ;  i_node++  ){
				if(  ( i_node == 3 ) || ( i_node == 7 )  ){
					real_index[ i_node ] = neigh_data->_GetNodeIndex( tet_index[ i_node ] );
					if(  !msh_->IsNodeSetted( real_index[ i_node ] )  ){
						double* coord = new double[ 3 ];
						key_type keys[ 3 ];
						neigh_cell->GetKey( tet_index[ i_node ] , keys );
						for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
							coord[ i_dim ] = ( static_cast<double>( keys[ i_dim ] ) * scale );
						}
						msh_->AddNode( real_index[ i_node ] , coord , 
													 false , 
													 neigh_data->_GetSign( tet_index[ i_node ] , 0 )  );

					}
				}else{
					real_index[ i_node ] = src_data->_GetNodeIndex( tet_index[ i_node ] );
					if(  !msh_->IsNodeSetted( real_index[ i_node ] )  ){
						double* coord = new double[ 3 ];
						key_type keys[ 3 ];
						src_cell->GetKey( tet_index[ i_node ] , keys );
						for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
							coord[ i_dim ] = ( static_cast<double>( keys[ i_dim ] ) * scale );
						}
						msh_->AddNode( real_index[ i_node ] , coord , 
													 false , 
													 src_data->_GetSign( tet_index[ i_node ] , 0 )  );
					}
				}
			}
			//Adding elements to mesh in this case are added 2 elements
			for(  int i_elem = 0  ;  i_elem < 2  ;  i_elem++  ){
				msh_->AddElement( real_index[ i_elem * 4 + 0 ] , real_index[ i_elem * 4 + 1 ] , 
													real_index[ i_elem * 4 + 2 ] , real_index[ i_elem * 4 + 3 ] );
			}
		}else{
			int tet_index[ 4 ];
			int real_index[ 4 ];
			GetLocalTetIndexOnCellWhenFaceIsSameRefinedAndEdgeSameOrLess( direction , i_edge , tet_index );
			for(  int i_node = 0  ;  i_node < 4  ;  i_node++  ){
				if(  i_node == 3  ){
					real_index[ i_node ] = neigh_data->_GetNodeIndex( tet_index[ i_node ] );
					if(  !msh_->IsNodeSetted( real_index[ i_node ] )  ){
						double* coord = new double[ 3 ];
						key_type keys[ 3 ];
						neigh_cell->GetKey( tet_index[ i_node ] , keys );
						for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
							coord[ i_dim ] = ( static_cast<double>( keys[ i_dim ] ) * scale );
						}
						msh_->AddNode( real_index[ i_node ] , coord , 
													 false , 
													 neigh_data->_GetSign( tet_index[ i_node ] , 0 )  );
					}
				}else{
					real_index[ i_node ] = src_data->_GetNodeIndex( tet_index[ i_node ] );
					if(  !msh_->IsNodeSetted( real_index[ i_node ] )  ){
						double* coord = new double[ 3 ];
						key_type keys[ 3 ];
						src_cell->GetKey( tet_index[ i_node ] , keys );
						for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
							coord[ i_dim ] = ( static_cast<double>( keys[ i_dim ] ) * scale );
						}
						msh_->AddNode( real_index[ i_node ] , coord , 
													 false , 
													 src_data->_GetSign( tet_index[ i_node ] , 0 )  );
					}
				}

			}
			//Adding elements to mesh in this case are added 2 elements
			msh_->AddElement( real_index[ 0 ] , real_index[ 1 ] , real_index[ 2 ] , real_index[ 3 ] );
		}
	}
}

/**
 *Creating eight tetrahedra on cell face
 *@param[in] cell Is the cell where elements are created
 *@param[in] direction Is the cell face where elements will be created

 */
void Mesher::CreateTetrahedraWhenNeighbourIsMoreRefined(  OctreeCell* cell , int direction  ){
	double scale = 1.00 / ( 1 << ROOT_LEVEL );
	int tet_index[ 32 ];
	GetLocalTetIndexOnFaceWhenNeighbourIsMoreRefined( direction , tet_index );
	//Checking if nodes has been added on mesh
	DataInsideOctreeCell* data =  cell->pGetData();
	int real_index[ 32 ];
	for(  size_t i_node = 0  ;  i_node <= 31  ;  i_node++  ){
		real_index[ i_node ] = data->_GetNodeIndex( tet_index[ i_node ] );
		if(  !msh_->IsNodeSetted( real_index[ i_node ] )  ){
			double* coord = new double[ 3 ];
			key_type keys[ 3 ];
			cell->GetKey( tet_index[ i_node ] , keys );
			for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
				coord[ i_dim ] = ( static_cast<double>( keys[ i_dim ] ) * scale );
			}
			msh_->AddNode( real_index[ i_node ] , coord , 
										 false , 
										 data->_GetSign( tet_index[ i_node ] , 0 )  );
		}
	}

	//Adding elements to mesh in this case are added 2 elements
	for(  int i_elem = 0  ;  i_elem < 8  ;  i_elem++  ){
		msh_->AddElement( real_index[ i_elem * 4 + 0 ] , real_index[ i_elem * 4 + 1 ] , 
											real_index[ i_elem * 4 + 2 ] , real_index[ i_elem * 4 + 3 ] );
	}
}

/**
 *function to adjust mesh created to the original boundary
 *This function just projet exterior(sign =1) and not setted nodes (sign = -5) to the closer triangle
 *@param[out] n_nodes  Is the pointer to set the number of nodes in the tetrahedral mesh
 *@param[out] n_elems  Is the pointer to set the number of elements in the tetrahedral mesh
 *@param[out] nodes_index  Is the pointer to set the real indexes of the nodes to be used
 *@param[out] elems_index  Is the pointer to set the real indexes of the elements to be used
 */
void Mesher::FitMeshToBoundary(   size_t n_nodes , size_t n_elems , int* nodes_index , int* elems_index  ){
	for(  size_t i_node = 0  ;  i_node < n_nodes  ;  i_node++  ){
		nodes_index[ i_node ] = 0;
	}
	for(  size_t i_elem = 0  ;  i_elem < n_elems  ;  i_elem++  ){
		elems_index[ i_elem ] = 0;
	}
	msh_->GetAllNodesAndElemsOnBodyFittedMesh(  nodes_index , elems_index  );
	msh_->ProjectNodesToboundary( local_root_ , nodes_index , r_level_ );
	msh_->CalculateBoundaryTetrahedraVolume(  );
	msh_->EraseNotValidTetrahedra( n_nodes , n_elems , nodes_index , elems_index );
}




//SENDERS AND RECEIVERS

//DEBUG
/**
 *Checking if constrain 2 to 1 is accomplished
 *This method checks on each cell if all its neighbours have at most one level of 
 *difference in the refinement level
 *@return A bool value indicating if the constrain is accomplished or not
 */
bool Mesher::CheckLocalConstrain2to1(){
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );
	for(  char i_level = MIN_LEVEL  ;  i_level < ROOT_LEVEL - 1  ;  i_level++  ){
		for(  size_t i_cell = 0  ;  i_cell < leaves.size()  ;  i_cell++  ){
			if(  leaves[ i_cell ]->GetLevel() == i_level  ){
				key_type neighbour_key[ 3 ];
				for(  int i_direction = 0  ;  i_direction < 18  ;  i_direction++  ){
					if(  leaves[ i_cell ]->GetNeighbourKey( i_direction , neighbour_key )  ){
						OctreeCell* neighbour_cell = this->GetCellOnLocalOctree( neighbour_key );
						if(  neighbour_cell->GetLevel() > i_level + 1  ){
							leaves.clear();
							return false;
						}
					}
				}
			}
		}
	}
	return true;
}

/**
 *Checking if all linear and center nodes have the same sign on the nodes
 */
void Mesher::CheckNodesSign(){
	OctreeCell_vector leaves;
	this->GetAllLocalRootLeaves( leaves );
	size_t n_leaves = leaves.size();
	for(  size_t i_leaf = 0  ;  i_leaf < n_leaves  ;  i_leaf++  ){
		DataInsideOctreeCell *data = leaves[ i_leaf ]->pGetData();
		for(  size_t i_node = 0  ;  i_node < 27  ;  i_node++  ){
			int signs[ 3 ];
			signs[ 0 ] = data->_GetSign( i_node , 0 );
			signs[ 1 ] = data->_GetSign( i_node , 1 );
			signs[ 2 ] = data->_GetSign( i_node , 2 );
			if(  !( (signs[ 0 ] == signs[ 1 ]) && (signs[ 1 ] == signs[ 2 ]) )  ){
				int real_sign;
				if(  NodeSignCanBeSet( signs , &real_sign )  ){
				    data->_SetNodePositionRespectToBoundary( 0 , i_node , real_sign );
				    data->_SetNodePositionRespectToBoundary( 1 , i_node , real_sign );
				    data->_SetNodePositionRespectToBoundary( 2 , i_node , real_sign );					
				}
			}
		}
	}
}



//SAVING ON FILE
/**
 *Saving body fitted tetrahedral mesh on GiD file
 */
void Mesher::PrintLocalBodyFittedTetrahedralMeshOnGiDFile( size_t n_nodes , size_t n_elems , 
																													 int* nodes_index , int* elems_index ){
	msh_->SaveBodyFittedMeshOnGiDFile( output_ , n_nodes , n_elems , nodes_index , elems_index );
}


//TO BE TESTED




/**
 *Testing if a value exists on a vector of doubles
 *@param[in] values Is the vector of elements
 *@param[in] test Is the value to test if exists of not on the vector values
 *@return A bool value indicating if tests exist on vector values or not
 */
bool DoesExistCoordOnVector( std::vector<double> values , double test ){
	for(  size_t i_elem = 0  ;  i_elem < values.size()  ;  i_elem++  ){
		if(  fabs( values[ i_elem ] - test ) < 1e-5  ){
			return true;
		}
	}
	return false;
}

















