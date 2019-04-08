#ifndef PARTITION_HPP
#define PARTITION_HPP

#include <vector>
#include <unordered_map>
#include <utility> 
#include "matrix.h"

class partition
{
  
  std::function<double()>& rnd;
  
  int length;
  int max_number_of_partitions;
  
  // hash_table which stores the monomials of the bell polynomials and the bell factor; the key is a vector {m,k}
  
  std::vector< std::vector< std::vector< std::pair< std::vector<int>,double > > > > possible_partitions;
  
public:
  
  partition(int size, int n_partitions, std::function<double()>& base ): rnd(base), length(size), max_number_of_partitions(n_partitions)
  {
     // explicied enumeration of possible_partitions if max_number_of_partitions=2; must be imported from e.g. mathematica 
    
    // B_1,1
    int number_recursion=1;
    int number_partitions=1;
    // monomials with their factor
    std::vector<int> partition_sizes={1};
    double bell_factor=1;
    
    std::vector<int> n_rec_n_part={number_recursion,number_partitions};
    std::vector< std::pair< std::vector<int>,double > > aux1={std::make_pair(partition_sizes,bell_factor)};
    std::vector< std::vector< std::pair< std::vector<int>,double > > > aux2(max_number_of_partitions+1);
    possible_partitions.resize(length+1,aux2);
    
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    
    // B_2,1
    number_recursion=2;
    number_partitions=1;
    // monomials with their factor
    partition_sizes={2};
    bell_factor=1;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    
    // B_2,2
    number_recursion=2;
    number_partitions=2;
    // monomials with their factor
    partition_sizes={1,1};
    bell_factor=1;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    
    // B_3,1
    number_recursion=3;
    number_partitions=1;
    // monomials with their factor
    partition_sizes={3};
    bell_factor=1;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    
    // B_3,2
    number_recursion=3;
    number_partitions=2;
    // monomials with their factor
    partition_sizes={1,2};
    bell_factor=3;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    
    // B_4,1
    number_recursion=4;
    number_partitions=1;
    // monomials with their factor
    partition_sizes={4};
    bell_factor=1;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    // B_4,2
    number_recursion=4;
    number_partitions=2;
    // monomials with their factor
    partition_sizes={2,2};
    bell_factor=3;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    // monomials with their factor
    partition_sizes={1,3};
    bell_factor=4;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    // B_5,1
    number_recursion=5;
    number_partitions=1;
    // monomials with their factor
    partition_sizes={5};
    bell_factor=1;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    // B_5,2
    number_recursion=5;
    number_partitions=2;
    // monomials with their factor
    partition_sizes={2,3};
    bell_factor=10;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    // monomials with their factor
    partition_sizes={1,4};
    bell_factor=5;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    // B_6,1
    number_recursion=6;
    number_partitions=1;
    // monomials with their factor
    partition_sizes={6};
    bell_factor=1;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    // B_6,2
    number_recursion=6;
    number_partitions=2;
    // monomials with their factor
    partition_sizes={3,3};
    bell_factor=10;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    // monomials with their factor
    partition_sizes={2,4};
    bell_factor=15;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    // monomials with their factor
    partition_sizes={1,5};
    bell_factor=6;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    // B_7,1
    number_recursion=7;
    number_partitions=1;
    // monomials with their factor
    partition_sizes={6};
    bell_factor=1;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
    // B_7,2
    number_recursion=7;
    number_partitions=2;
    // monomials with their factor
    partition_sizes={3,4};
    bell_factor=35;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    // monomials with their factor
    partition_sizes={2,5};
    bell_factor=21;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    // monomials with their factor
    partition_sizes={1,6};
    bell_factor=7;
    
    n_rec_n_part={number_recursion,number_partitions};   
    possible_partitions[n_rec_n_part[0]][n_rec_n_part[1]].push_back(std::make_pair(partition_sizes,bell_factor));
    
//    test_possible_partitions(7,2);
//    test_possible_partitions(7,1);
//    test_possible_partitions(2,2);
//    test_possible_partitions(3,2);
//    test_possible_partitions(3,1);
//    test_possible_partitions(4,1);
//    
//    exit(1);
  };
  

  std::pair<std::vector<int>,double>& get_random_monomial(int m, int k)
  {
    
#ifdef DEBUG_BELL_POLY   
    std::cout << "getting random monomial of Bell polynomial B_" << m << "," << k << std::endl;
#endif
    
    std::vector<std::pair<std::vector<int>,double>>& possible_monomials=possible_partitions[m][k];
    
    int monomial=static_cast<int>(rnd()*(possible_monomials.size()));
    
#ifdef DEBUG_BELL_POLY   
    std::cout << "choosing random monomial " << monomial << std::endl;
    std::cout << "This monomial is " << possible_monomials[monomial].second << "* " ;
    for(size_t i=0; i<possible_monomials[monomial].first.size(); i++)
    {
      std::cout << possible_monomials[monomial].first[i] << " " ; 
    }
    std::cout << std::endl;
#endif
      
    return (possible_monomials[monomial]);
    
  }
  
  int get_number_possible_mon(int m, int k)
  {
    return possible_partitions[m][k].size();
  }
  
  void init_partition_hierarchy()
  {
//     int current_recursion=static_cast<int>(rnd()*(max_number_of_partitions+1));
//     int current_monomial = rnd()<0.5 ? 1 : 2;
    
    
  }
  
  
  void test_possible_partitions(int m, int k)
  {
    
    std::vector<std::pair<std::vector<int>,double>> bell=possible_partitions[m][k];
    
    std::cout << "For the Bell Polynomial " << m << " " << k << " There is a contribution from" << std::endl ;
    for(size_t i=0;i<bell.size();i++)
    { 
      std::cout << "monomial " << i << ":";
      
      std::vector<int> contri=bell[i].first;
      for(size_t j=0;j<contri.size();j++)
      {
	std::cout << " u " << contri[j];
      }
      
      std::cout << " with a factor " << bell[i].second << std::endl;
    }
   
  }

  

    


    
};



#endif
