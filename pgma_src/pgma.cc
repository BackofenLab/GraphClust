/*
  Cluster according to given distance matrix by PGMA

  Implementation of {U,W}PGMA
  
  reads distance list
  writes NEWICK tree format

  This programm reads a list of triples <i> <j> <d(i,j)>
  (one triple per line)

*/

#include "pgma_aux.hh"

using namespace std;


int
main(int argc, char **argv) {
    
    if (argc!=3 && argc!=2) {
	cerr << "USAGE: "<<argv[0]<<" <names> [<distances>]"<<endl;
	exit(-1);
    }
    
    ifstream names_in(argv[1]);
    
    Descriptions descriptions(names_in);
    
    DistMatrix dist_mat;
    dist_mat.resize(descriptions.size());
    
    if (argc==3) {
	ifstream dists_in(argv[2]);
	dist_mat.read(dists_in);
    } else {
	dist_mat.read(cin);
    }
    
    Tree tree;

    tree.pgma(dist_mat,descriptions,true);
    
    tree.print();
    
    cout <<endl;
    exit(0);
}
