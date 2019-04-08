#ifndef TREE_HPP
#define TREE_HPP

#include <unordered_map>
#include <functional>
#include <random>
#include <vector>
#include <list>
#include <math.h>
#include <iostream>
#include "leaf.hpp"
#include "partition.hpp"

#define max_number_nodes 30
#define cutoff 10
#define dim 1


// the key is the numerical name of the leaf
typedef int key;

class tree
{
  
  std::function<double()>& rnd;
  
  partition& possible_branches;
  
  int number_integration_variables;
  
  std::vector<leaf> mContainer_tree;
  std::list<int> active_leaf_list;
  
  int max_element_count;
  int element_count;
  std::vector<int> namelist;
  // at which position are the names in the mContainer_tree // mContainer_tree[container_pos[namelist[i]]]
  std::vector<int> container_pos;
  
public:
  
  tree(std::function<double()>& base, partition& possible_branches_ref): rnd(base), possible_branches(possible_branches_ref), number_integration_variables(0), max_element_count(max_number_nodes),
  element_count(0), namelist(max_element_count,0), container_pos(max_number_nodes,-1)
  {
    for(size_t i=0;i<namelist.size();i++)
    {
      namelist[i]=i;
      container_pos[i]=i;
    }
    leaf aux(1,1);
    mContainer_tree.resize(max_number_nodes,aux);
    
  };
  
  double factorial(int n)
  {
    double val = 1;
    for(int i=n; i>1; --i)
       val*=i; 
     return val; 
  }
   
  int get_number_int_varialbes() { return number_integration_variables;}
  int get_element_count() { return element_count;}
   
  void make_new_random_tree(int given_m, const std::vector<double>& root_variable);
  void grow_branch(key root_key);
  void make_random_subtree(int given_m, const std::vector<double>& root_variable, double& gen_prob);
  void grow_branch_sub(key root_key, double& gen_prob);
   
  double get_weight();
  double get_sub_weight(leaf& from_this_root);
  void get_subtree_facts(leaf& sub_root, std::vector<double>& root_variable, int& number_integration_variables, std::list<int>& sub_namelist, double& sub_gen_prob, double& sub_weight, int& number_elements);
  void DFS_FUNC_WEIGHT(key vertex_key, double& weight);
  void DFS_FUNC_FACTS(leaf& sub_root, int& number_integration_variables, std::list<int>& sub_namelist, double& sub_gen_prob, double& sub_weight, int& number_elements);
   
  void merge_subtree(tree& new_subtree, int int_variables_old_subtree, const std::list<int>& namelist_old_tree);
  void DFS_FUNC_NAMES(leaf& root, std::list<int>& names, int& number_elements, std::vector<int>& translating_new_names_to_old_keys);
   
  void append(std::vector<double>& new_int_var);
   
  void shrink();
   
  double kernel(const std::vector<double>& x, const std::vector<double>& t) {return x[0]-t[0];}
  double u(int m,const std::vector<double>& x)
  {
     if(m==0) 
     {
       return log(1+x[0])+2*log(2)-2*x[0]*(log(2)-1)*(log(2)-1)-5./4.;
     }
     else if(m==1) 
     {
       return 1./288.*(-773.+6*x[0]*(427.+8.*log(2)*(-182.+log(2)*(189.+8.*(-7.+log(2))*log(2))))-16.*log(2)*(-259+log(4)*(193.+(-28.+log(8))*log(8))));
     }
     else 
     {
       std::cout << "Error: you have no access to u_" << m << std::endl;
       exit(1);
     }
  }
   
  std::vector<int> get_random_leaf();
  leaf& get_single_random_root();
  leaf& get_prev_leaf(leaf& this_leaf);
  void get_next_leaves(leaf& this_leaf, std::vector< key >& next_leaves);
  const std::vector<double>& get_next_int_var(key this_key) {return mContainer_tree[container_pos[mContainer_tree[container_pos[this_key]].get_next_names()[0]]].get_integration_variable();}
  leaf& get_root() {return mContainer_tree[container_pos[0]];}
  leaf& get_leaf(key leaf_key) {return mContainer_tree[container_pos[leaf_key]];}
  std::vector<leaf>& get_tree_container() { return mContainer_tree;}
   
   
  void print_tree()
  {
    std::cout << std::endl << "# printing the current tree with external variable ";
    for(auto iter : mContainer_tree[container_pos[0]].get_integration_variable())
      std::cout << iter << " ";
    std::cout << std::endl;
    
    int current_number_int_variables = 0;
    
    DFS_FUNC_PRINT(mContainer_tree[container_pos[0]], current_number_int_variables);
    std::cout << std::endl;
    
    if( current_number_int_variables != number_integration_variables )
    {
      std::cout << "number of integration_variables don't match " << current_number_int_variables << " " << number_integration_variables << std::endl;
    }
    
  }
  void DFS_FUNC_PRINT(leaf& root, int& current_number_int_variables)
  {
     
    if(root.get_m()<=1)
    {
    std::cout << "# From u_" << root.get_m() << " a new branch can't be grown " << std::endl;
    return;
    }

    std::cout << "# growing new branch from u " << root.get_m() << " with name " << root.get_name() << " according to the Bell polynomial B_" << root.get_m()-1 << "," << root.get_k() << std::endl;//     
    std::vector<int> n_n=root.get_next_names();

    root.print_next_names();

    std::cout << "# The branch type is "; 
    for(size_t i=0;i<n_n.size();i++)
      std::cout << "u_" << mContainer_tree[container_pos[n_n[i]]].get_m() << " ";
    std::cout << "and has a int variable ";
    for(auto iter : mContainer_tree[container_pos[n_n[0]]].get_integration_variable())
      std::cout << iter << " ";
    std::cout << std::endl;
    current_number_int_variables++;
     
    for(size_t i=0;i<n_n.size();i++)
    {
      if(n_n[i]==root.get_name())
      {
        std::cout << " The next name can't be the same as the root name "<<std::endl;
        exit(1);    
      }
      for(size_t j=0; j<mContainer_tree[container_pos[n_n[i]]].get_next_names().size(); j++)
      {
        if(mContainer_tree[container_pos[n_n[i]]].get_next_names()[j]==root.get_name())
        {
          std::cout << " The next name of the next name can't be the root name" << std::endl;
          exit(1);
        }
      }

    DFS_FUNC_PRINT(mContainer_tree[container_pos[n_n[i]]], current_number_int_variables);
    }
  }

};



#endif
