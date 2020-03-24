#pragma once

//C system files
#include <math.h>
#include <assert.h>
#include <omp.h>

//C++ system files
#include <vector>
#include <iostream>

//Other libraries
#include "../mesh/mesh.h"
#include "../octree/octree_driver.h"

using namespace std;

/**
 *mesher class
 *This class contains the information requiered to create a mesh based on a octree using
 *parallel computing with MPI
 *@param nodes_label Is the value that controls the actual node index setted
 *@param input_ Is the name of the file containing the input boundary mesh
 *@param output_ Is the name of the file where the tetrahedral mesh will be stored
 *@param r_level_ Is the refinement level that the user requires in the mesh
 *@param octree_ Is the octree pointer that the master has created
 *@param local_root_ Is the local cell that each process will use to create its own mesh
 *@param bnd_ Is the pointer to the boundary mesh information
 *@param msh_ Is the pointer to the tetrahedral mesh
 */
class Mesher{ 

	int nodes_label_;
	char* input_;
	char* output_;
	int r_level_;
	OctreeDriver* octree_;
	OctreeCell* local_root_;
	Boundary* bnd_;
	Mesh* msh_;

	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Mesher();
		Mesher( int argc , char** argv );
		~Mesher();

		//GETS
		void GetAllLocalRootLeaves( OctreeCell_vector& leaves );
		OctreeCell* GetCellOnLocalOctree( key_type* keys );
		int  GetMinRefinementLevelUsed();
		void GetMinKeyRoot( key_type* keys );
		void GetMaxKeyRoot( key_type* keys );
		size_t GetNumberOfTetrahedralNodes();
		size_t GetNumberOfTetrahedraOnMesh();


		//SETS
		void SetNodesPositionRespectToBoundary();
		void SetSignOnNodeAndNeighbours( OctreeCell* cell , int node , int value , int dim );
		void SetAllOctreeNodes(  );
		void SetAllLocalOctreeNodes();
		void SetCenterNodeOnCell( OctreeCell* cell , bool interfaz );
		bool SetNodeInCellAndNeighbors( int node , OctreeCell* cell , double* coord , bool interfaz );
		void SetNodeInCell( int node , OctreeCell* cell , double* coord , bool interfaz );
		void SetAllLinearAndCenterNodesOnBoundaryCellFace(  int direction , OctreeCell* cell , bool interfaz );
		void SetAllNodesOnCellFace( int src_dir , OctreeCell* cell , bool interfaz );
		void SetCenterNodeOnCellFace( int src_dir , OctreeCell* cell , bool interfaz );
		void SetAllLinearNodesOnCellFace( int src_dir , OctreeCell* cell , bool interfaz );
		void SetAllQuadraticNodesOnCellFace( int src_dir , OctreeCell* cell , bool interfaz );
		void SetAllNodesOnLocalNeighbourSameRefined( int direction , OctreeCell* src_cell , bool interfaz  );
		void SetAllLinearNodesOnNeighborsCellsByFace( int direction , OctreeCell* src_cell , bool interfaz );

		//UTILITIES
		void ScaleBoundaryNodes();
		void CalculateAABBForBoundaryTriangles();
		void RefineLocalRoot();
		void SetTrianglesOnLocalRoot();
		void SubdivideIntersectedCells();\
		void BalanceOctree(    );
		void BalanceLocalOctree();
		void AllocateAllLeavesData();
		void RayCastingOverDirection( key_type* root_min_key , key_type* root_max_key , 
																	int min_refinement , int n_rays , int dim );
		bool RayIntersectsModelBBox( key_type* keys , int dim );
		void RayCasting( key_type* root_min_key , key_type* root_max_key , key_type* keys , int dim );
		void CreateTetrahedralMeshFromOctree( );
		void CreateTetrahedralMeshOnLocalOctree( );
		void CreateTetrahedraOnBoundaryFace( int direction , OctreeCell* cell );
		void CreateTetrahedraWhenNeighbourIsLessRefined(  OctreeCell* src_cell , OctreeCell* neigh_cell , int direction  );
		void CreateTetrahedraWhenNeighbourIsSameRefined(  OctreeCell* src_cell , OctreeCell* neigh_cell , int direction  );
		void CreateTetrahedraWhenNeighbourIsMoreRefined(  OctreeCell* cell , int direction  );
		void FitMeshToBoundary(  size_t n_nodes , size_t n_elems , int* nodes_index , int* elems_index );


		//SENDERS AND RECEIVERS

		//DEBUG
		bool CheckLocalConstrain2to1();
		void CheckNodesSign();

		//SAVING ON FILE
		void PrintLocalBodyFittedTetrahedralMeshOnGiDFile( size_t n_nodes , size_t n_elems , 
																											 int* nodes_index , int* elems_index );
};

bool DoesExistCoordOnVector( std::vector<double> values , double test );







