#pragma once


// System includes
#include <string>
#include <iostream> 
#include <omp.h>


	//This class represents a cell in an octree to be used with Octree class

	/** This class represents a cell in an octree and holds its level, min_key
		*  ,the children of the cell and pointer to a data class defined by
		*  configuration class.
		*  The level_ start from ROOT_LEVEL and each children has 1 level less than
		*  its parents but more than MIN_LEVEL.
		*  The keys are the binary bisection pattern in each dimension for the
		*  min point of the cell.
		*  The children are stored in in an array of 8 in following order
		*
    *                    Top
    *
    *              ________________
    *             /   6   /   7   /|
    *            /_______/_______/ |
    *           /   4   /   5   /| |
    *          /_______/_______/ | /  Front
    *          |       |       | |/
    *          |       |       | /
    *          |_______|_______|/
    *              ________________
    *    Left     /   2   /   3   /|    Right
    *            /_______/_______/ |
    *           /   0   /   1   /| |
    *          /_______/_______/ | /
    *          |       |       | |/
    *          |       |       | /
    *          |_______|_______|/
    *  z   Back
    *  | y           Bottom
    *  |/__x
    *
   */
template<class TConfiguration>
class OctreeBinaryCell {

	public:
		typedef typename TConfiguration::data_type data_type;
		typedef TConfiguration configuration_type;                
		typedef typename TConfiguration::pointer_type pointer_type;
		typedef std::vector<pointer_type> object_container_type;
		typedef std::size_t key_type;

		//Default constructor.
		OctreeBinaryCell( char Level = 29 ) : level_( Level ) , children_( NULL ), data_( NULL ){
			for(  std::size_t i = 0  ;  i < 3  ;  i++  )
				min_key_[ i ] = 0;
		}

		//Destructor.
		virtual ~OctreeBinaryCell(){
			if( data_ ) configuration_type::DeleteData( data_ );
				delete []children_;
		}

		void DeleteChildren(){
			delete []children_;
			children_ = NULL;
		}

		void DeleteData(){
			if( data_ ){
				configuration_type::DeleteData( data_ );
				data_ = NULL;
			}
		}

		std::size_t GetChildIndex( key_type x_key , key_type y_key , key_type z_key )const{
			char next_level = (char)( level_ - 1);
			key_type level_bit = 1 << next_level;
			return (((x_key & level_bit) >> next_level) + (((y_key & level_bit) >> next_level) << 1) + 
							(((z_key & level_bit) >> next_level) << 2));
		}

		int SubdivideCell( ){
			if(  level_ == 0  )
				return 1;
			if(  children_  )
				return 1;

			children_ = new OctreeBinaryCell[ 8 ];
			char next_level = (char)( level_ - 1 );
			for(  std::size_t i = 0  ;  i < 8  ;  i++  ){
				children_[ i ].SetMinKey( min_key_[ 0 ] | ( ( i & 1 ) << next_level ), 
																	min_key_[ 1 ] | ( ( ( i & 2 ) >> 1 ) ) << next_level, 
																	min_key_[ 2 ] | ( ( ( i & 4 ) >> 2 ) ) << next_level );
				children_[ i ].SetLevel( next_level );
			}
			return 0;
		}

		void GetMinPointNormalized( double* min_point )const{
			for(  std::size_t i = 0  ;  i < 3  ;  i++  ){
				min_point[ i ] = GetCoordinateNormalized( min_key_[ i ] );
			}
		}

		void GetMaxPointNormalized( double* max_point )const{
			double size = CalcSizeNormalized( );
			for(  std::size_t i = 0  ;  i < 3  ;  i++  ){
				max_point[ i ] = GetCoordinateNormalized( min_key_[ i ] ) + size;
			}
		}

		int GetLeftKey( key_type* keys )const{
			if(  min_key_[ 0 ] >= 1  ){
				keys[ 0 ] = min_key_[ 0 ] - 1;
				keys[ 1 ] = min_key_[ 1 ];
				keys[ 2 ] = min_key_[ 2 ];
				return 1; //Neighbor founded
			}
			return 0; //No neighbor
		}

		int GetRightKey( key_type* keys )const{
			if(  min_key_[ 0 ] < static_cast<key_type> ( ( 1 << 29 ) - ( 1 << level_ ) ) ){
				keys[ 0 ] = min_key_[ 0 ] + (static_cast<key_type>(1) << level_);
				keys[ 1 ] = min_key_[ 1 ];
				keys[ 2 ] = min_key_[ 2 ];
				return 1; //Neighbor founded
			}
			return 0; //No neighbor
		}

		int GetBackKey( key_type* keys )const{
			if(  min_key_[ 1 ] >= 1  ){
				keys[ 0 ] = min_key_[ 0 ];
				keys[ 1 ] = min_key_[ 1 ] - 1;
				keys[ 2 ] = min_key_[ 2 ];
				return 1; //Neighbor founded
			}
			return 0; //No neighbor
		}

		int GetFrontKey( key_type* keys )const{
			if( min_key_[ 1 ] < static_cast<key_type> ( ( 1 << 29 ) - ( 1 << level_ ) ) ){
				keys[ 0 ] = min_key_[ 0 ];
				keys[ 1 ] = min_key_[ 1 ] + (static_cast<key_type>(1) << level_);
				keys[ 2 ] = min_key_[ 2 ];
				return 1; //Neighbor founded
			}
			return 0; //No neighbor
		}

		int GetBottomKey( key_type* keys )const{
			if(  min_key_[ 2 ] >= 1  ) {
				keys[ 0 ] = min_key_[ 0 ];
				keys[ 1 ] = min_key_[ 1 ];
				keys[ 2 ] = min_key_[ 2 ] - 1;
				return 1; //Neighbor founded
			}
			return 0; //No neighbor
		}

		int GetTopKey( key_type* keys )const{
			if(  min_key_[ 2 ] < static_cast<key_type> ( ( 1 << 29 ) - ( 1 << level_ ) )  ){
				keys[ 0 ] = min_key_[ 0 ];
				keys[ 1 ] = min_key_[ 1 ];
				keys[ 2 ] = min_key_[ 2 ] + (static_cast<key_type>( 1 ) << level_);
				return 1; //Neighbor founded
			}
			return 0; //No neighbor
		}

		int GetKey( std::size_t position , key_type* keys )const{
			//in total, there are a maximum of 27 nodes (as the Quadratic9 Hexa in GiD). The only difference is that here
			//the central node of the hexa (cell) is just after the lineal nodes.
			//the node index              					:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
			static const std::size_t x_position[ ] = { 0, 2, 2, 0, 0, 2, 2, 0, 1, 1, 2, 1, 0, 0, 2, 2, 0, 1, 2, 1, 0, 1, 1, 2, 1, 0, 1};
			static const std::size_t y_position[ ] = { 0, 0, 2, 2, 0, 0, 2, 2, 1, 0, 1, 2, 1, 0, 0, 2, 2, 0, 1, 2, 1, 1, 0, 1, 2, 1, 1};
			static const std::size_t z_position[ ] = { 0, 0, 0, 0, 2, 2, 2, 2, 1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 0, 1, 1, 1, 1, 2};

			keys[ 0 ] = min_key_[ 0 ] + ( ( x_position[ position ] ) << ( level_ - 1 ) );
			keys[ 1 ] = min_key_[ 1 ] + ( ( y_position[ position ] ) << ( level_ - 1 ) );
			keys[ 2 ] = min_key_[ 2 ] + ( ( z_position[ position ] ) << ( level_ - 1 ) );
			return 1;
		}

		int GetNeighbourKey( std::size_t direction , key_type* keys )const{
			//direction 0 : X < 0
			//direction 1 : X > 0
			//direction 2 : Y < 0
			//direction 3 : Y > 0
			//direction 4 : Z < 0
			//direction 5 : Z > 0
			//from 6 to 18 are the directions coorresponding to edges of the cell

			//the key returned is inside the cell (minkey +1), to ensure that the corresponding
			//cell get in pGetCell is the right one.

			assert( direction < 18 );						  //0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17
			static const std::size_t x_offset[ ] = { 0 , 2 , 1 , 1 , 1 , 1 , 0 , 2 , 0 , 2 , 0 , 2 , 0 , 2 , 1 , 1 , 1 , 1 };
			static const std::size_t y_offset[ ] = { 1 , 1 , 0 , 2 , 1 , 1 , 0 , 0 , 2 , 2 , 1 , 1 , 1 , 1 , 0 , 2 , 0 , 2 };
			static const std::size_t z_offset[ ] = { 1 , 1 , 1 , 1 , 0 , 2 , 1 , 1 , 1 , 1 , 0 , 0 , 2 , 2 , 0 , 0 , 2 , 2 };
			static const std::size_t x_coef[ ]   = { 0 , 1 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 1 , 0 , 0 , 0 , 0 };
			static const std::size_t y_coef[ ]   = { 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 1 };
			static const std::size_t z_coef[ ]   = { 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 0 , 0 , 1 , 1 };

			std::size_t size = ( 1<<level_ );

			// here i'm adding 2 to each dimension and it won't overflow
			keys[ 0 ] = min_key_[ 0 ] + x_offset[ direction ] + x_coef[ direction ] * size;
			keys[ 1 ] = min_key_[ 1 ] + y_offset[ direction ] + y_coef[ direction ] * size;
			keys[ 2 ] = min_key_[ 2 ] + z_offset[ direction ] + z_coef[ direction ] * size;

			for(  int i = 0  ;  i < 3  ;  i++  ){
				if(  keys[ i ] == 0  ){
					return 0; //No neighbor
				}else{
					( keys[ i ] )--;
				}
				if(  keys[ i ] > static_cast<key_type>( 1 << 29 )  )
					return 0; //No neighbor
			}
			return 1; //Neighbor founded
		}

		int GetNeighbourKey( std::size_t position , std::size_t direction , key_type* keys )const{

			GetKey( position , keys );
			keys[ 0 ] += ( direction & 1 ) << 1;
			keys[ 1 ] += ( direction & 2 );
			keys[ 2 ] += ( direction & 4 ) >> 1;

			for(  int i = 0  ;  i < 3  ;  i++  ){
				if(  keys[ i ] == 0  )
					return 0; //No neighbor
				else
					( keys[ i ] )--;

				if(  keys[ i ] > static_cast<key_type>(1 << 29)  )
					return 0; //No neighbor
			}
			return 1; //Neighbor founded
		}

		std::size_t GetLocalPosition( key_type* keys ){
			key_type position[ 3 ];
			//in total, there are a maximum of 27 nodes (as the Quadratic9 Hexa in GiD). The only difference is that here
			//the central node of the hexa (cell) is just after the lineal nodes.
			//the node index                      : 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
			static const std::size_t local_index[]={0, 9, 1,12,21,10, 3,11, 2,13,22,14,25, 8,23,16,24,15, 4,17, 5,20,26,18, 7,19,6};

			for(  std::size_t i = 0  ;  i < 3  ;  i++  ){
				position[ i ] = ( keys[ i ] - min_key_[ i ] ) >> ( level_ - 1 );
			}
			std::size_t index = position[ 0 ] + position[ 1 ] * 3 + position[ 2 ] * 9;
			assert( index <= 26 );
			return local_index[ index ];

		}

		bool ExistNodeOnCell( key_type* keys , int* index ){
			int local_index[ 27 ] = { 0 , 13 ,  4 , 12 , 25 , 20 ,  3 , 16 ,  7 , 
																9 , 22 , 17 , 21 ,  8 , 26 , 11 , 24 , 19 , 
																1 , 14 ,  5 , 10 , 23 , 18 ,  2 , 15 ,  6 };
			for(  int x_pos = 0  ;  x_pos < 3  ;  x_pos++  ){
				for(  int y_pos = 0  ;  y_pos < 3  ;  y_pos++  ){
					for(  int z_pos = 0  ;  z_pos <3  ;  z_pos++  ){
						size_t diff[ 3 ];
						diff[ 0 ] = ( min_key_[ 0 ] + ( x_pos * ( 1<<( level_ - 1 ) ) ) ) - keys[ 0 ]; 
						diff[ 1 ] = ( min_key_[ 1 ] + ( y_pos * ( 1<<( level_ - 1 ) ) ) ) - keys[ 1 ];
						diff[ 2 ] = ( min_key_[ 2 ] + ( z_pos * ( 1<<( level_ - 1 ) ) ) ) - keys[ 2 ];
						if(  ( diff[ 0 ] == 0 ) && ( diff[ 1 ] == 0 ) && ( diff[ 2 ] == 0 )  ){
							(*index) = local_index[ x_pos * 9 + y_pos * 3 + z_pos ];
						  return true;
						}
					}
				}
			}
			return false;
		}


		void Insert( pointer_type object ){
			objects_.push_back( object );
		}

		void TransferObjectsToSonsNormalized(){
          
			if( !objects_.size() ) return;
				assert( this->HasChildren() );

			const double tolerance = 0.001 * double( 1 << 2 ) / double( 1 << 29 ) ; // 0.1% of the min size
			double min_coord[ 3 ] = { 0.00 , 0.00 , 0.00 };
			double max_coord[ 3 ] = { 0.00 , 0.00 , 0.00 };

			for(  std::size_t i = 0  ;  i < 8  ;  i++  ){            
				OctreeBinaryCell* son = pGetChild( i );
				if(  son->HasChildren( )  ){
					son->TransferObjectsToSonsNormalized( );
					continue;
				}
				son->GetMinPointNormalized( min_coord );
				son->GetMaxPointNormalized( max_coord );
				pointer_type object;
				for(  int i = 0  ;  i < (int)objects_.size( )  ;  i++  ){
					object=objects_[ i ];
					const int is_intersected = configuration_type::IsIntersected( object , tolerance , 
																																				min_coord , max_coord );
					if(  is_intersected  )
						son->Insert( object );
				}
			}
			object_container_type temp;
			objects_.swap( temp );                   
		}

		bool TransferObjectsToSons(){
			if(  !objects_.size()  ){
				return true;
			}
			if(  !this->HasChildren()  ){
				return true;
			}
			const double tolerance = 0.001 * double( 1 << 2 ) / double( 1 << 29 ) ; // 0.1% of the min size
			double min_coord[ 3 ] = { 0.00 , 0.00 , 0.00 };
			double max_coord[ 3 ] = { 0.00 , 0.00 , 0.00 };

			for(  std::size_t i_child = 0  ;  i_child < 8  ;  i_child++  ){
				OctreeBinaryCell* son = this->pGetChild( i_child );
				son->GetMinPointNormalized( min_coord );
				son->GetMaxPointNormalized( max_coord );
				for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
					min_coord[ i_dim ] -= tolerance;
					max_coord[ i_dim ] += tolerance;
				}
				pointer_type object;
				for(  int i_obj = 0  ;  i_obj < (int)objects_.size( )  ;  i_obj++  ){
					object = objects_[ i_obj ];
					bool is_intersected = object->IntersectsBox( min_coord , max_coord );
					if(  is_intersected  )
						son->Insert( object );
				}
			}
			object_container_type temp;
			objects_.swap( temp ); 
			return false;
		}
	
		bool ProjectNodeToCloserTriangle( double* coord ){
			size_t n_objects = objects_.size( );
	
			//Projecting
			double new_coord[ 3 ];
			bool flag = false;
			for(  size_t i_obj = 0  ;  i_obj < n_objects  ;  i_obj++  ){
				flag = objects_[ i_obj ]->ProjectPointToTriangle( coord , new_coord );
				if( flag ){
					coord[ 0 ] = new_coord[ 0 ];
					coord[ 1 ] = new_coord[ 1 ];
					coord[ 2 ] = new_coord[ 2 ];
					return flag;
				}
			} 
			return flag;
		}

		bool ProjectNodeToCloserSide( double* coord ){
			size_t n_objects = objects_.size();
			double min_dist = 1e20;
			double possible_point[ 3 ];
			double min_point[ 3 ];
			min_point[ 0 ] = objects_[ 0 ]->GetCoord( 0 , 0 );
			min_point[ 1 ] = objects_[ 0 ]->GetCoord( 0 , 1 );
			min_point[ 2 ] = objects_[ 0 ]->GetCoord( 0 , 2 );
			for(  size_t i_obj = 0  ;  i_obj < n_objects  ;  i_obj++  ){
				for(  int i_side = 0  ;  i_side < 3  ;  i_side++  ){
					double distance = objects_[ i_obj ]->DistanceToSide( i_side , coord , possible_point );
					if(  distance < min_dist  ){
						min_dist = distance;
						min_point[ 0 ] = possible_point[ 0 ];
						min_point[ 1 ] = possible_point[ 1 ];
						min_point[ 2 ] = possible_point[ 2 ];
					}
				}
			}
			coord[ 0 ] = min_point[ 0 ];
			coord[ 1 ] = min_point[ 1 ];
			coord[ 2 ] = min_point[ 2 ];
			return true;
		}

		double GetMinDistanceToNodes( double* coord , size_t* local_obj , int* local_node ){
			double min_dist = 1e20;
			(*local_obj) = 0;
			(*local_node) = 0;
			size_t n_objects = objects_.size();
			for(  size_t i_obj = 0  ;  i_obj < n_objects  ;  i_obj++  ){
				for(  int i_node = 0  ;  i_node < 3  ;  i_node++  ){
					double distance = objects_[ i_obj ]->DistanceToNode( i_node , coord );
					if(  distance < min_dist  ){
						(*local_obj) = i_obj;
						(*local_node) = i_node;
						min_dist = distance;
					}
				}
			}
			return min_dist;
		}

		void CopyNode( size_t min_obj , int min_node , double* coord ){
			objects_[ min_obj ]->Copynode( min_node , coord );
		}


		unsigned char GetLevel( )const{
			return level_;
		}

		size_t GetNNodes( ){
			double min_coord[ 3 ];
			double max_coord[ 3 ];
			double size = this->CalcSizeNormalized( );
			for(  int i_dim = 0  ;  i_dim < 3  ;  i_dim++  ){
				min_coord[ i_dim ] = this->GetCoordinateNormalized( min_key_[ i_dim ] );
				max_coord[ i_dim ] = min_coord[ i_dim ] + size;
			}
			std::vector<size_t> index;
			size_t n_objects = this->GetNObjects( );
			if(  !n_objects  ){
				return 0;
			}
			for(  size_t i_obj = 0  ;  i_obj < n_objects  ;  i_obj++  ){
				for(  int i_node = 0  ;  i_node < 3  ;  i_node++  ){
					//checking if node intersect cell
					if(  objects_[ i_obj ]->NodeIntersectsCell( i_node , min_coord , max_coord )  ){	
						size_t value = objects_[ i_obj ]->GetNodeIndex(  i_node  );
						bool flag = false;
						for(  size_t i_pos = 0  ;  i_pos < index.size( )  ;  i_pos++   ){
							if(  index[ i_pos ] == value  ){
								flag = true;
								break;
							}
						}
						if(  !flag  ){
							index.push_back( value );
							sort( index.begin( ) , index.end( ) , [](const size_t a, const size_t b) {return a > b; });
						}
					}
				}
			}
			return index.size( );
		}

		char SetLevel( char level ){
			level_ = level;
			return level_;
		}

		void GetMinKey( key_type& min_key_x , key_type& min_key_y , key_type& min_key_z )const{
			min_key_x = min_key_[ 0 ];
			min_key_y = min_key_[ 1 ];
			min_key_z = min_key_[ 2 ];
		}

		void SetMinKey( key_type min_key_x , key_type min_key_y , key_type min_key_z ){
			min_key_[ 0 ] = min_key_x;
			min_key_[ 1 ] = min_key_y;
			min_key_[ 2 ] = min_key_z;
		}

		OctreeBinaryCell& rGetChild( std::size_t pos )const{
			return children_[ pos ];
		}

		OctreeBinaryCell* pGetChild( std::size_t pos )const{
			return children_ + pos;
		}

		OctreeBinaryCell* pGetChild( key_type x_key , key_type y_key , key_type z_key )const{
			return pGetChild( GetChildIndex( x_key , y_key , z_key ) );
		}

		OctreeBinaryCell* GetChildren(){
			return children_;
		}

		OctreeBinaryCell const* GetChildren()const{
			return children_;
		}

		data_type* pGetData() const{
			return data_;
		}

		data_type** pGetDataPointer(){
			return &data_;
		}

		const std::vector<pointer_type>* pGetObjects() const{
			return &objects_;
		}

		std::vector<pointer_type>* pGetObjects(){
			return &objects_;
		}

		void EmptyObjects(){
			object_container_type tmp;     
			tmp.swap( objects_ );
		}
	
		size_t GetNObjects(){
			return objects_.size();
		}

		bool IsLeaf()const{
			return ( children_ == NULL );
		}

		bool HasChildren()const{
			return( children_ != NULL );
		}


	private:
		char level_;
		key_type min_key_[ 3 ];
		OctreeBinaryCell* children_;
		data_type* data_;
		object_container_type objects_;
       

		double CalcSizeNormalized()const{
			const double scale = 1.00 / ( 1 << 29 );

			return ( 1 << level_ ) * scale;
		}

		double GetCoordinateNormalized( key_type key )const{
			const double scale = 1.00 / ( 1 << 29 );

			return (static_cast<double>(key) * scale);
		}


		//Assignment operator.
		OctreeBinaryCell & operator=(OctreeBinaryCell const& rOther) {
			return *this;
		}

		//Copy constructor.
		OctreeBinaryCell(OctreeBinaryCell const& rOther) {
		}

};





