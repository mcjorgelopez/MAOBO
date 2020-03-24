#pragma once

//C system files
#include <math.h>
#include <assert.h>
#include <omp.h>

//C++ system files
#include <vector>
#include <iostream>
#include <string.h>

//Project includes 
#include "../octree/octree_driver.h"

/**
 *Tterahedron class
 *This class contains the information from tetrahedron
 *@param nodes_ Is the list of four nodes that form the tetrahedron
 */
class Tetrahedron{

	private:
	  int nodes_[ 4 ];
		double vol_;

	public:
		//CONSTRUCTOR AND DESTRUCTOR
		Tetrahedron();
  	Tetrahedron( int nod1 , int nod2 , int nod3 , int nod4 );
  	~Tetrahedron();
		
		//GETS
		int GetNodeIndex( int node );
		double GetVolume();

		//SETS

		//UTILITIES
		void CalculateTetraedronVolume( double* N0 , double* N1 , double* N2 , double* N3 );

		//DEBUG

};

/**
 *Tetrahedron mesh
 *This class contains a tetrahedral mesh
 *
 *@param nodes_ Is the list of nodes belonging to the boundary mesh
 *@param elems_ Is the list of tetrahedrons creating the mesh
 *@param scale_ Is the class to mapp between the AABB and the original domain
 */
class Mesh{

	private:
		std::vector<Node*> nodes_;
		std::vector<Tetrahedron*> elems_;
		Scaler* scale_;

	public:
		//CONSTRUCTOR, DESTRUCTOR
		Mesh();
		Mesh( Scaler* scale , int n_nodes );
		~Mesh();

		//GETS
		size_t GetNNodes();
		size_t GetNElements();
		Tetrahedron* GetElement( int i_elem );
		Node* GetNode( int i_node );
		double GetCoord( int node , int dim );
		size_t GetIndex( int i_elem , int i_node );
		int GetIndexRelatedToNodeColor( size_t node_index );

		//SETS
		void SetElement( Tetrahedron* elem );
		void SetNode( Node* node );
		void SetNode( Node* node , int index );

		//UTILITIES
		bool IsNodeSetted( int node );
		void AddNode( int index , double* coord , bool interfaz , char sign );
		void AddElement( int index0 , int index1 , int index2 , int index3 );
		void GetAllNodesAndElemsOnBodyFittedMesh(  int* nodes_index , int* elems_index  );
		void ProjectNodesToboundary( OctreeCell* root , int* nodes_index , int r_level );
		void CalculateBoundaryTetrahedraVolume(  );
		void EraseNotValidTetrahedra( size_t n_nodes , size_t n_elems , int* nodes_index , int* elems_index );

		//SAVING MESH ON FILE
		void SaveBodyFittedMeshOnGiDFile(  char* output_ , size_t n_nodes , size_t n_elems , 
																				 int* nodes_index , int* elems_index );
		void SaveNodesSignOnTxtFile( char* name , size_t n_nodes , size_t n_elems , int* nodes_index , 
																 int* elems_index );


		//DEBUG
		void PrintInfo();
};

