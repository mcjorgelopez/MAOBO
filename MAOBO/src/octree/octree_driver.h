#pragma once

//C system files
#include <math.h>
#include <assert.h>
#include <omp.h>

//C++ system files
#include <vector>
#include <algorithm>

//Other libraries
#include "octree_binary.h"
#include "../boundary/boundary.h"



using namespace std;

//Defines
typedef std::size_t key_type;

//relation between lineal position of nodes and position of sons. Lineal positions are numbered 
//following GiD criteria for hexas (right hand rule).
const int sons_positions[ 8 ] = { 0 , 1 , 3 , 2 , 4 , 5 , 7 , 6 };

const int OCTREE_MAX_LEVELS = 30;//max number of levels for octree. This is to avoid infinite recursions.
const int ROOT_LEVEL = OCTREE_MAX_LEVELS - 1;//max number of levels for octree. This is to avoid infinite recursions.
const int MIN_LEVEL = 2;//max number of levels for octree. This is to avoid infinite recursions.

/**
 *Data Inside octree cell structure
 *This structure contains the octree cell information.
 *@param is_outer_cell indicates if cell is inside, outside or on the boundary
 *@param NL_indexes Is the list of neighbours belonging to other rank process
 *@param points Is the list of points setted on the cell
 *@param node_sign Is the node position respect to boundary, equal to is_outer_cell
 *@param interfaz A flag indicating if a node is located on the interface between 2 ranks
 *@param r_index Is a list containing the real node index on the local nodes labelling
 */
struct DataInsideOctreeCell{
	private:
    char node_sign[ 3 ][ 27 ];
		int r_index[ 27 ];

	public:
		//CONSTRUCTOR AND DESTRUCTOR
  	DataInsideOctreeCell();
  	~DataInsideOctreeCell();

		//GETS
		int _GetSign( int node , int dim );
		int _GetNodeIndex( int node );

		//SETS
    void _SetNodePositionRespectToBoundary( int dim , int node , char value );
		void _SetNodeIndex( int node , int index );

		//UTILITIES
		bool _IsSetNodeSign( int dim , int node );
		bool _IsSetNode( int node );

		//DEBUG
		void _PrintInfo();
		bool _PrintNotSetSigns();
};

/**
 *Is the octree configuration class
 *@param CHILDREN_NUMBER Indicates the amount of childrens by cell
 *@param DIMENSION Is the mesh dimension 1-1D, 2-2D and 3-3D
 *@param MAX_LEVEL Is the maximum refinement level allowed by the octree
 *@param MIN_LEVEL Is the minimum refinement level allowed by the octree (to avoid overflow)
 */
class OctreeConfigure {
public:
  enum{
    CHILDREN_NUMBER = 8,
    DIMENSION = 3,
    MAX_LEVEL = OCTREE_MAX_LEVELS,
    MIN_LEVEL = 2 // this cannot be less than 2!!! Pooyan
  };

  typedef DataInsideOctreeCell data_type;
  typedef double* point_coords;
  typedef Triangle* pointer_type;


  static data_type* AllocateData(){
    return new data_type;
  }

  static void CopyData(data_type* source, data_type* destination){
    destination = source;
  }

  static void DeleteData(data_type* data){
    delete data;
  }        

  static int IsIntersected( const pointer_type elem , const double tolerance , const point_coords min_coord , const point_coords max_coord ){
      return elem->Intersects( min_coord , max_coord , tolerance );
  }

  static void GetBoundingBox( pointer_type elem , point_coords min_coord , point_coords max_coord ){
    elem->CalcBoundingBox(min_coord,max_coord);
  }
};


typedef OctreeBinaryCell<OctreeConfigure> OctreeCell;
typedef OctreeBinary<OctreeCell> Octree;
typedef std::vector<OctreeCell*> OctreeCell_vector;

/**
 *Class OctreeDriver
 *This class is an interface between the OctreeDriver and the OctreeBinary class.
 *@param octree_bin Is the pointer to the class OctreeBinary
 */
class OctreeDriver{

	private: 
  	Octree* octree_bin;

	public:
		//CONSTRUCTOR AND DESTRUCTOR
  	OctreeDriver();
  	~OctreeDriver();

		//GETS
  	Octree* GetOctree()const;
  	OctreeCell* GetRoot();
  	OctreeCell* GetOctreeCell( const double *coord )const;
  	OctreeCell* GetOctreeCell( Octree::key_type* keys )const;
  	int CalcAllLeavesVector( OctreeCell_vector*& all_leaves );
  	void GetLeavesInBoundingBox( const double* coord1 , const double* coord2 , OctreeCell_vector& leaves )const;
  	OctreeCell* GetNeighbour( const OctreeCell* cell , const int idir )const;
	  OctreeCell* GetLeftNeighbour( const OctreeCell* cell )const;
  	OctreeCell* GetRightNeighbour( const OctreeCell* cell )const;
  	OctreeCell* GetTopNeighbour( const OctreeCell* cell )const;
  	OctreeCell* GetBottomNeighbour( const OctreeCell* cell )const;
  	OctreeCell* GetFrontNeighbour( const OctreeCell* cell )const;
  	OctreeCell* GetBackNeighbour( const OctreeCell* cell )const;
  	void GetCoordOctreePosition( const OctreeCell* cell , int ipos , double* coord_point )const;
	  double GetCoordinate( Octree::key_type key )const;

		//SETS

		//UTILITIES
  	bool CheckIsBalanced()const;
  	int BalanceOctree();
  	int SubdivideCell( OctreeCell* cell )const; 
  	int RefineOctreeWithSize( const double size );
  	bool IsCellEmpty( const OctreeCell* cell )const;
  	double CalcSize( const OctreeCell* cell )const;

    //DEBUG
		void PrintInfo();
};


int GetKindOfCaseOnLocalRoot( key_type* root_max_key , key_type* keys );
void PerturbateKeysDependingOnKind( key_type* keys , int kind );
void UnPerturbateKeysDependingOnKind( key_type* keys , int kind );
void GetNodesToAddByRayCasting( int position , int dim , int* nodes );
void GetNodesToAddByRayCastingX( int position , int* nodes );
void GetNodesToAddByRayCastingY( int position , int* nodes );
void GetNodesToAddByRayCastingZ( int position , int* nodes );
int GetPositionOverRay( std::vector<double> points , key_type key );
bool NodeSignCanBeSet( int* signs , int* real_sign );
int GetCenterIndexFromFace( int src_dir );
void GetLinearIndexFromFace( int src_dir , int* index );
void GetQuadraticIndexFromFace(int src_dir , int* index );
void GetAllNeighbourKeysFromFaceLessRefined( int direction , key_type* keys0 , key_type* keys1 , 
																						 key_type* keys2 , key_type* keys3 , int level );
void GetLinearIndexFromEdge( int src_dir , int* index );
int GetCenterIndexFromEdge( int src_dir );
void GetAllNeighbourKeysFromEdgeLessRefined( int direction , key_type* keys0 , key_type* keys1 , 
																						 int level );
void GetLinearAndCenterIndexFromBoundaryFace( int direction , int* index );
void GetLinearIndexOnFace( int direction , int* index );
void GetEdgesIndexesFromFace( int direction , int* index );
int GetCenterIndexOnEdge( int direction );
int GetRelativePositionOverNeighbourLessRefined( OctreeCell* src_cell , OctreeCell* cell , int direction );
bool AreEquals(  key_type* key0 , key_type* key1  );
void GetLocalTetIndexWhenNeighbourIsMoreRefined( int position , int src_dir , int* tet_index );
void GetLocalTetIndexOnInterfaceWhenNeighbourIsSameRefinedAndEdgeSame( int src_dir , int i_edge , int* tet_index );
void GetLocalTetIndexOnInterfaceWhenNeighbourIsSameRefinedAndEdgeMore( int src_dir , int i_edge , int* tet_index );
void GetLocalTetIndexOnInterfaceWhenNeighbourIsLessRefined( int i_cell , int src_dir , int* tet_index );
void GetQuadraticIndexOnFace( int direction , int* index );
void GetLocalTetIndexOnBoundaryCellWhenEdgeIsMoreRefined( int direction , int i_edge , int* tet_index );
void GetLocalTetIndexOnBoundaryCellWhenEdgeIsSameRefined( int direction , int i_edge , int* tet_index );
void GetLocalTetIndexWhenNeighbourIsLessRefined( int position , int direction , int* tet_index );
bool SourceCellIsLower( OctreeCell* src_cell , OctreeCell* other_cell );
void GetLocalTetIndexOnCellWhenFaceIsSameRefinedAndEdgeSameOrLess( int direction , int i_edge , int* tet_index );
void GetLocalTetIndexOnCellWhenFaceIsSameRefinedAndEdgeMore( int direction , int i_edge , int* tet_index );
void GetLocalTetIndexOnFaceWhenNeighbourIsMoreRefined( int direction , int* tet_index );




