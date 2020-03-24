#pragma once


//System includes
#include <string>
#include <iostream> 
#include <assert.h>
#include <omp.h>



//Project includes
#include "octree_binary_cell.h"
template <class TCellType>
class OctreeBinary {

  public:
    typedef TCellType cell_type;
    typedef typename cell_type::key_type key_type;
    typedef typename cell_type::configuration_type configuration_type;
    typedef double coordinate_type;

    //Default constructor.
    OctreeBinary() : root_( new cell_type ) , number_of_cells_( 8 + 1 ) , number_of_leaves_( 1 ) , levels_( 0 ){}

    //Destructor.
   	virtual ~OctreeBinary(){
			delete root_;
		}

		double CalcSizeNormalized( const cell_type* cell )const{
			const double scale = 1.0 / ( 1 << 29 );
			return ( 1 << cell->GetLevel( ) ) * scale;
		}

		double CalcMinCellNormalizedSize( )const{          
			const double scale = 1.0 / ( 1 << 29 );    
			return ( 1 << 2 ) * scale;
		}

		key_type CalcKeyNormalized( coordinate_type coordinate )const{
			assert( coordinate>=0.0 ); 
			assert( coordinate<=1.0 );
			return static_cast<key_type> ( ( 1 << 29 ) * coordinate );
		}

		void InsertNormalized( coordinate_type* point ){
			key_type x_key = CalcKeyNormalized( point[ 0 ] );
			key_type y_key = CalcKeyNormalized( point[ 1 ] );
			key_type z_key = CalcKeyNormalized( point[ 2 ] );
			cell_type* cell = root_;
			for(  std::size_t i = 29  ;  i > 2  ;  i--  ){
				if(  cell->IsLeaf()  ){
					SubdivideCell( cell );
				}
				cell = cell->pGetChild( x_key , y_key , z_key );
			}
		}

		bool BoundingBoxesCollide( const double min_cornerA[ 3 ] , const double max_cornerA[ 3 ] ,
			const double min_cornerB[ 3 ] , const double max_cornerB[ 3 ] , const double tolerance )const{
			for(  int idim = 0  ;  idim < 3  ;  idim++ ){
				if(  min_cornerA[ idim ] > max_cornerB[ idim ] + tolerance ) 
					return false;
				if(  max_cornerA[ idim ] < min_cornerB[ idim ] - tolerance ) 
					return false;
			}
			return true;
		}

		int GetAllLeavesVector( std::vector<cell_type*>& all_leaves )const{
			std::vector<cell_type*> cells_stack;
			cells_stack.push_back( root_ );
			while(  !cells_stack.empty()  ){
				cell_type* cell = cells_stack.back( );
				cells_stack.pop_back( );
				if(  cell->HasChildren( )  ){
					for(  std::size_t i = 0  ;  i < 8  ;  i++  ){
						cells_stack.push_back( cell->pGetChild( i ) );
          }
        }else{
        	all_leaves.push_back( cell );
				}
      }
      return 0;
    }

		cell_type * pGetCellNormalized( const coordinate_type * point )const{
			key_type keys[ 3 ];
			keys[ 0 ] = CalcKeyNormalized( point[ 0 ] );
			keys[ 1 ] = CalcKeyNormalized( point[ 1 ] );
			keys[ 2 ] = CalcKeyNormalized( point[ 2 ] );
			return pGetCell( keys );
		}

		cell_type * pGetCell( key_type * keys )const{
			cell_type* cell = root_;
			for(  std::size_t i = 0  ;  i < 29  ;  i++  ){
				if(  cell->IsLeaf( )  ){
					return cell;
				}
				cell = cell->pGetChild( keys[ 0 ] , keys[ 1 ] , keys[ 2 ] );
			}
			return cell;
		}

		cell_type * pGetCell( key_type* keys , std::size_t level )const{
			cell_type* cell = root_;
			for(  std::size_t i = 29  ;  i > level  ;  i--  ){
				if(  cell->IsLeaf( )  ){
					return cell;
				}
				cell = cell->pGetChild( keys[ 0 ] , keys[ 1 ] , keys[ 2 ] );
			}
			return cell;
		}

		cell_type * pGetLeftCell( const cell_type * p_cell ){
			key_type keys[ 3 ];
			if( p_cell->GetLeftKey( keys ) ){
				return pGetCell( keys );
			}
			return NULL; //No neighbor
		}

		cell_type * pGetLeftCell( cell_type* p_cell , std::size_t level ){
			key_type keys[ 3 ];
			if(  p_cell->GetLeftKey( keys ) ){
				return pGetCell( keys , level );
			}
			return NULL; //No neighbor
		}

		cell_type * pGetRightCell( const cell_type * p_cell ){
			key_type keys[ 3 ];
			if(  p_cell->GetRightKey( keys )  ){
				return pGetCell( keys );
			}
			return NULL; //No neighbor
		}

		cell_type * pGetRightCell( cell_type* p_cell , std::size_t level ){
			key_type keys[ 3 ];
			if(  p_cell->GetRightKey( keys )  ){
				return pGetCell( keys , level );
			}
			return NULL; //No neighbor
		}

		cell_type * pGetBackCell( const cell_type * p_cell ){
			key_type keys[ 3 ];
			if(  p_cell->GetBackKey( keys ) ){
				return pGetCell( keys );
			}
			return NULL; //No neighbor
		}

		cell_type * pGetBackCell( cell_type* p_cell , std::size_t level ){
			key_type keys[ 3 ];
			if(  p_cell->GetBackKey( keys )  ){
				return pGetCell( keys , level );
			}
			return NULL; //No neighbor
		}

		cell_type * pGetFrontCell( const cell_type * p_cell ){
			key_type keys[ 3 ];
			if(  p_cell->GetFrontKey( keys )  ){
				return pGetCell( keys );
			}
			return NULL; //No neighbor
		}

		cell_type * pGetFrontCell( cell_type* p_cell , std::size_t level ){
			key_type keys[ 3 ];
			if(  p_cell->GetFrontKey( keys )  ){
				return pGetCell( keys , level );
			}
			return NULL; //No neighbor
		}

		cell_type * pGetTopCell( const cell_type * p_cell ){
			key_type keys[ 3 ];       
			if(  p_cell->GetTopKey( keys )  ){
				return pGetCell( keys );
			}
			return NULL; //No neighbor
		}

		cell_type * pGetTopCell(  cell_type* p_cell , std::size_t level ){
			key_type keys[ 3 ];
			if(  p_cell->GetTopKey( keys )  ){
				return pGetCell( keys , level );
			}
			return NULL; //No neighbor
		}

		cell_type * pGetBottomCell( const cell_type * p_cell ){
			key_type keys[ 3 ];
			if(  p_cell->GetBottomKey( keys )  ){
				return pGetCell( keys );
			}
			return NULL; //No neighbor
		}

		cell_type * pGetBottomCell( cell_type* p_cell , std::size_t level ){
			key_type keys[ 3 ];
			if(  p_cell->GetBottomKey( keys )  ){
				return pGetCell( keys , level );
			}
			return NULL; //No neighbor	
		}

		cell_type * pGetNeighbourCell( const cell_type* p_cell , std::size_t direction ){
			key_type keys[ 3 ];

			if(  p_cell->GetNeighbourKey(  direction , keys )  ){
				return pGetCell( keys );
      }
      return NULL; //No neighbor
		}

		cell_type * pGetNeighbourCell( cell_type* p_cell , std::size_t position , std::size_t direction ){
			key_type keys[ 3 ];
			if(  p_cell->GetNeighbourKey( position , direction , keys )  ){
				return pGetCell( keys );
			}
			return NULL; //No neighbor
		}

		int SubdivideCell( cell_type* p_cell ){
			number_of_cells_ += 8;
			number_of_leaves_ += 8 - 1;
			return p_cell->SubdivideCell();
		}

		int SubvidiveUntilSizeNormalized( double* coord , const double desired_size ){
			key_type x_key = CalcKeyNormalized( coord[ 0 ] );
			key_type y_key = CalcKeyNormalized( coord[ 1 ] );
			key_type z_key = CalcKeyNormalized( coord[ 2 ] );
			cell_type* cell = root_;
			std::size_t scaled_size = std::size_t( desired_size * ( 1 << 29 ) );
			if(  scaled_size < ( 1 << 2 )  )
				scaled_size = ( 1 << 2 );

			for(  std::size_t i = 29  ;  ( std::size_t( 1 ) << i ) > scaled_size  ; i--  ){
				if(  cell->IsLeaf()  ){
					SubdivideCell( cell );
				}
				cell = cell->pGetChild( x_key , y_key , z_key );
			}
			return 0;
		}

		int RefineWithUniformSizeNormalized( const double uniform_size ){
			const double min_size = double( 1 << 2 ) / double( 1 << 29 );
			double cell_size = uniform_size;
			if(  cell_size < min_size  ){
				cell_size = min_size;
			}
			std::vector<cell_type*> cells_stack;
			cells_stack.push_back( root_ );
			while(  !cells_stack.empty( )  ){
				cell_type* cell = cells_stack.back( );
				cells_stack.pop_back( );
				if(  CalcSizeNormalized( cell ) > cell_size  ){
					if(  cell->IsLeaf( )  ){
						SubdivideCell( cell );
					}
					for(  std::size_t i = 0  ;  i < 8  ;  i++  )
						cells_stack.push_back( cell->pGetChild( i ) );
				}
			}
			return 0;
		}
           

		cell_type* pGetRoot( ){
			return root_;
		}

		double GetCoordinateNormalized( key_type key )const{
			const double scale = 1.00 / ( 1 << 29 );
			return ( static_cast<double>( key ) * scale );
		}

	private:

		cell_type* root_;
		std::size_t number_of_cells_;
		std::size_t number_of_leaves_;
		std::size_t levels_;

		//Assignment operator.
		OctreeBinary & operator=( OctreeBinary const& rOther ){
			return *this;
		}

		//Copy constructor.
		OctreeBinary( OctreeBinary const& rOther ){
		}



};





