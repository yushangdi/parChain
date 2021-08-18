
class pointset {
  private:
    // we assume points in R^dimension
    int dimension;

    // number of points in Pointset
    int size_;

    // Points[i][j] gives the jth coordinate of pint i 
    float** Points;

  public:

    pointset(int dimension_) {
      dimension=dimension_;
      size_=0;
    }

    pointset(char* name) {
      std::ifstream is(name);
      is>>dimension; is>>size_;
      float h;
      Points=new float*[size_];
      for(int i=0; i<size_; i++) {
        Points[i]=new float[dimension];
        for(int j=0; j<dimension; j++) {
          is>>h; Points[i][j]=h;
        }
      }
    }

    // we create a pointset of the given size of point sampled uniformely at random
    // We assume that this function is called exactly once
    void uniform_random(int s) {
      size_=s;
      Points=new float*[size_];
      for(int i=0; i<size_; i++) {
       Points[i]=new float[dimension];
       for(int j=0; j<dimension; j++) {
         Points[i][j]=(float) random()/RAND_MAX;
       }
      }
    }

    float getx(int i) {
      return Points[i][0];
    }
    
    float gety(int i) {
      return Points[i][1];
    }
    
    float getz(int i) {
      return Points[i][2];
    } 
    // We compute the distance of the points i and j
    float althausDistance(int i, int j) const {
      float d=0;
      for(int k=0; k<dimension; k++) d+=(Points[i][k]-Points[j][k])*(Points[i][k]-Points[j][k]);
      return d; //sqrt(d);
    }

    // Return number of Points in the set
    int size() const {
      return size_;
    }
    
    void print_point(int i) {
      std::cout<<"("<<getx(i)<<","<<gety(i)<<","<<getz(i)<<")";
    };
    
    void print_file(char* s) {
      std::ofstream os(s);
      os<<dimension<<" "<<size_<<"\n";
      for(int i=0; i<size_; i++) {
        for(int j=0; j<dimension; j++) {
          os<<Points[i][j];
          if(j<dimension-1) os<<" "; else os<<std::endl;
        }
      }
    }

    void check_clustering(std::list<std::set<int> >& Clustering, double threshold) {
      std::cout<<"Check Clustering\n";
      for(std::list<std::set<int> >::iterator C1=Clustering.begin(); C1!=Clustering.end(); C1++) {
        //std::cout<<(*C1).size()<<" new cluster\n";
        for(std::list<std::set<int> >::iterator C2=C1; C2!=Clustering.end(); C2++) {
          if(C2==C1) {
            //std::cout<<"same cluster\n";
            for(std::set<int>::iterator P1=(*C1).begin(); P1!=(*C1).end(); P1++) {
              //std::cout<<"second\n";
              for(std::set<int>::iterator P2=P1; P2!=(*C1).end(); P2++) {
                //std::cout<<"hier\n";
                //std::cout<<*P1<<" "<<*P2<<" pair of points\n";
                if(althausDistance(*P1,*P2)>threshold) {
                  std::cout<<*P1<<" "<<*P2<<" ERROR: distance within cluster of "<<sqrt(althausDistance(*P1, *P2))<<"\n";
                  //return;
                }
              }
            }
          } else {
            //std::cout<<"different cluster\n";
            double d1, d=0;
            for(std::set<int>::iterator P1=(*C1).begin(); P1!=(*C1).end(); P1++) {
              for(std::set<int>::iterator P2=(*C2).begin(); P2!=(*C2).end(); P2++) {
                d1=althausDistance(*P1, *P2);
                if(d1>d) d=d1;
               }
            }
            if(d<=threshold) {
              for(std::set<int>::iterator P1=(*C1).begin(); P1!=(*C1).end(); P1++) { print_point(*P1); std::cout<<" "; }
              std::cout<<std::endl;
              for(std::set<int>::iterator P2=(*C2).begin(); P2!=(*C2).end(); P2++) { print_point(*P2); std::cout<<" "; }
              std::cout<<"\n ERROR: distance between clusters of "<<sqrt(d)<<"\n";
              //return;
            }
          }
        }
      }
    }
    
    
    // we delete the arrays in which the points were stored
    // We do not check (but assume), whether there were created  
    ~pointset() {
      for(int i=0; i<size(); i++) delete[] Points[i];
      delete[] Points;
    }
};
