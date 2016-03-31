#include "pgma_aux.hh"
#include "iostream"

std::string itoa(int i) {
    std::stringstream s;
    s<<i;
    return s.str();
}

// actually reads a list of index pairs with similarity/distance
// LINE FORMAT <id1> <id2> <sim>
void DistMatrix::read(std::istream &in) {
    int i; 
    int j;
    float dist;
 
    while(in >> i >> j >> dist) {
	ref(i-1,j-1) = dist;
    }
}

void Descriptions::read(std::istream &in) {
    std::string name;
    while (in >> name) push_back(name);
}

void Tree::Node::print() const {
    if (left==NULL && right==NULL) 
	std::cout<<label;
    else {
	std::cout <<"(";
	if (left!=NULL)
	    left->print();
	else std::cout<<"()";
	std::cout <<",";
	if (right!=NULL)
	    right->print();
	else std::cout<<"()";
	std::cout <<")";
    }
    std::cout <<":"<<length;
}


// ATTENTION: we use distmat destructively!
// apply pgma to distance matrix and return tree
void Tree::pgma(DistMatrix &d, Descriptions &desc, bool wpgma)
{
    
    std::vector<Node*> cluster;
    std::vector<unsigned int> cluster_idxs;

    //init
    unsigned int ncluster = d.size(); // number of clusters
    
    cluster.resize(ncluster);
    cluster_idxs.resize(ncluster);
    
    for (unsigned int i=0; i<ncluster; ++i) {
	cluster[i]=new Node(desc[i],i);
	cluster_idxs[i]=i;
    }

    while (ncluster>1) {
	// search for pair of cluster with minimal distance
	int min_i=-1;
	int min_j=-1;
	int min_xj=-1;
	double min_d=1.0e100;
	
	for (unsigned int xi=0; xi<ncluster; ++xi) {
	  int i=cluster_idxs[xi];
	  for (unsigned int xj=0; xj<ncluster; ++xj) {
	    int j=cluster_idxs[xj];
	    if (j<i && d.ref(i,j)<min_d) {
	      min_d = d.ref(i,j);
	      min_i = i;
	      min_j = j;
	      min_xj = xj;
	    }
	  }
	}
	
	// merge the two clusters min_i and min_j
	
	double height=d.ref(min_i,min_j)/2.0;
	
	int size_i=cluster[min_i]->size;
	int size_j=cluster[min_j]->size;
	
	cluster[min_i]->length = height-cluster[min_i]->height;
	cluster[min_j]->length = height-cluster[min_j]->height;

	cluster[min_i]=new Node(cluster[min_i],
				cluster[min_j],
				height);
	
	
	// remove cluster min_j
	ncluster--;
	cluster_idxs[min_xj]=cluster_idxs[ncluster];
	cluster_idxs.resize(ncluster);
	
	
	// compute new distances
	// row and column for i and j are not needed anymore
        // we reuse the row/column i for distances to the new cluster
	for (unsigned int xk=0; xk < ncluster; ++xk) {
	  unsigned int k=cluster_idxs[xk];
	  if (k!=(unsigned int)min_i)
	    if (wpgma)
	      d.ref_o(min_i,k)=
		(size_i*d.ref_o(min_i,k)
		 + size_j*d.ref_o(min_j,k))
		/ (size_i + size_j);
	    else
	      d.ref_o(min_i,k) = (d.ref_o(min_i,k) + d.ref_o(min_j,k))/2.0;
	}
	d.ref(min_i,min_j)=0;	
    }
    
    root = cluster[cluster_idxs[0]];
    root->length=0;
}
