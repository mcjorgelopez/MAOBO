#include "./comunicator/comunicator.h"

int main( int argc , char** argv ){

	omp_set_num_threads( atoi(argv[4]) );

	Mesher* com = new Mesher( argc , argv );

	double t0 , t1 , t2; 
	t0 = omp_get_wtime();
	t1 = omp_get_wtime();
	t2 = omp_get_wtime() - t0;
	com->ScaleBoundaryNodes();
	com->CalculateAABBForBoundaryTriangles();
	t2 = omp_get_wtime() - t0;
	fprintf(stdout, "[%16lf]-[%16lf]   **Scaling boundary         : \n", t2 , omp_get_wtime() - t1 );


	t1 = omp_get_wtime();
	com->RefineLocalRoot();
	t2 = omp_get_wtime() - t0;
	fprintf(stdout, "[%16lf]-[%16lf]   **Refining octree          : \n", t2 , omp_get_wtime() - t1 );


	t1 = omp_get_wtime();
	com->BalanceOctree( );
	t2 = omp_get_wtime() - t0;
	fprintf(stdout, "[%16lf]-[%16lf]   **Constrain 2-1            : \n", t2 , omp_get_wtime() - t1 );



	t1 = omp_get_wtime();	
	com->SetNodesPositionRespectToBoundary();
	t2 = omp_get_wtime() - t0;
	fprintf(stdout, "[%16lf]-[%16lf]   **Coloring algorithm       : \n", t2 , omp_get_wtime() - t1 );




	t1 = omp_get_wtime();	
	com->SetAllOctreeNodes(  );
	t2 = omp_get_wtime() - t0;
	fprintf(stdout, "[%16lf]-[%16lf]   **Setting nodes            : \n", t2 , omp_get_wtime() - t1 );



	t1 = omp_get_wtime();
	com->CreateTetrahedralMeshFromOctree( );
	t2 = omp_get_wtime() - t0;
	fprintf(stdout, "[%16lf]-[%16lf]   **Setting tetrahedrons     : \n", t2 , omp_get_wtime() - t1 );




	t1 = omp_get_wtime();
	size_t n_nodes = com->GetNumberOfTetrahedralNodes();
	size_t n_elems = com->GetNumberOfTetrahedraOnMesh();
	int* nodes_index = new int[ n_nodes ];
	int* elems_index = new int[ n_elems ];
	com->FitMeshToBoundary( n_nodes , n_elems , nodes_index , elems_index );
	t2 = omp_get_wtime() - t0;
	fprintf(stdout, "[%16lf]-[%16lf]   **Fitting mesh to boundary : \n", t2 , omp_get_wtime() - t1 );

	t1 = omp_get_wtime();
	com->PrintLocalBodyFittedTetrahedralMeshOnGiDFile( n_nodes , n_elems , nodes_index , elems_index );
	t2 = omp_get_wtime() - t0;
	fprintf(stdout, "[%16lf]-[%16lf]   **Printing file            : \n", t2 , omp_get_wtime() - t1 );


	delete[] nodes_index;
	delete[] elems_index;
	delete com;
  return 0;
}
