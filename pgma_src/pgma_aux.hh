#ifndef PGMA_AUX_HH
#define PGMA_AUX_HH

#include <string>
#include <vector>
#include <iostream>
#include <limits.h>
#include <ctype.h>

#include <sstream>
#include <fstream>

#include <cstdlib>
#include <assert.h>

std::string itoa(int i);

class Descriptions : public std::vector<std::string> {
protected:
    void read(std::istream &in);
public:
    Descriptions(std::istream &in) {
	read(in);
    }
    Descriptions(int x=10) : std::vector<std::string>(x) {
    }
};

// store only entries i,j where j<=i
class DistMatrix : private std::vector<std::vector<float> > {
public:
  unsigned int size() {return std::vector<std::vector<float> >::size();}

// actually reads a list of index pairs with similarity/distance
// LINE FORMAT <id1> <id2> <sim>
  void read(std::istream &in);

  void resize(unsigned int d1) {
    std::vector<std::vector<float> >::resize(d1);
    for (unsigned int i=0; i<d1; ++i)
      (*this)[i].resize(i+1);
  }

  float &ref(int i,int j) {
    assert(j<=i);
    return (*this)[i][j];
  }

  const float &ref(int i,int j) const {
    assert(j<=i);
    return (*this)[i][j];
  }

  float &ref_o(int i,int j) {
    if (j<=i)
      return (*this)[i][j];
    else
      return (*this)[j][i];
  }

 const float &ref_o(int i,int j) const {
    if (j<=i)
      return (*this)[i][j];
    else
      return (*this)[j][i];
  }

  void print(const std::vector<std::string> &labels);
};



class Tree {

    class Node {
    public:
      //int label;
      std::string label;
	Node *left;
	Node *right;
	double height;
	double length;
	int size;
	
	Node(const std::string &label_="", double length_=0):
	  label(label_),left(NULL),right(NULL),height(0),length(length_),size(1) {}
	Node(Node *_left, Node *_right, double _height, double length_=0): 
	    left(_left),right(_right),height(_height),length(length_),
	    size(_left->size+right->size) {}

	bool isLeave() const {return left==NULL && right==NULL;}
	
	void print() const;
    };

    Node *root;
    
private:
    void deepDelete(Node *n) {
	if (n!=NULL) {
	    deepDelete(n->left);
	    deepDelete(n->right);
	    delete n;
	}
    }

    double readFloat(std::istream &in) {
	std::string t="";
	char c;
	while (in >> c && (isdigit(c) || c=='.'  || c=='+'  || c=='-'  || c=='e' )) {
	  t+=c;
	}
	in.unget();
	return atof(t.c_str());
    }
    
    std::string readLabel(std::istream &in) {
	std::string t="";
	char c;
	while (in >> c && (isalnum(c) || c=='_' || c=='-')) {t+=c;}
	in.unget();
	return t;
    }
    
    Node *read_tree(std::istream &in) {
      char c;
      in >> c;
      if(c=='(') {
	    Node *left;
	    Node *right;
	    double length;
	    
	    if ((left=read_tree(in))==NULL) return NULL;
	    if (!(in >> c && c==',')) {
	      std::cerr<<"expected ',', read "<<c <<std::endl;
	      delete left;
	      return NULL;
	    } 
	    if ((right=read_tree(in))==NULL) {delete left; return NULL;}
	    if (!(in >> c && c==')')) {
	      std::cerr<<"expected ')', read "<<c <<std::endl;
	      delete left;
	      delete right;
	      return NULL;
	    } 
	    
	    if (in >> c && c==':') {
		length=readFloat(in);
	    } else {
	      in.unget();
	      length=1;
	    }
	    
	    double height=std::max(left->height+left->length,right->height+right->length);
	    
	    return new Node(left,right,height,length);
	} else {
	  in.unget();
	    std::string label=readLabel(in);
	    double length;
	    if (in>>c && c==':') {
	        length=readFloat(in);
	    } else {
	      in.unget();
	      length=1;
	    }
	    return new Node(label,length);
	}
    }
    
    void leaveLabels(Node *n,std::vector<std::string> &s) {
	if (n==NULL) return;
	if (n->isLeave()) s.push_back(n->label);
	leaveLabels(n->left,s);
	leaveLabels(n->right,s);
    }


public:
    Tree():root(NULL) {
    }
    
    ~Tree() {
	deepDelete(root);
    }

    void pgma(DistMatrix &d, Descriptions &desc, bool wpgma=true);
    
    void print() const {
	if (root==NULL) {
	    std::cout << "()";
	} else {
	    root->print();
	}
    };
    double height() const {if(root==NULL) return 0; else return root->height;}
    
    // reads a newick tree 
    void read(std::istream &in) {
        root=read_tree(in);
	if (root==NULL) {
	    std::cerr << "Cannot parse tree!"<<std::endl;
	};
    }

    std::vector<std::string> leaveLabels() {
	std::vector<std::string> labels;
	leaveLabels(root,labels);
	return labels;
    }

};



#endif // PGMA_AUX_HH
