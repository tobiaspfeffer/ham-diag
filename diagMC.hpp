#ifndef DIAGMC_HPP
#define DIAGMC_HPP

#include <vector>
#include <list>
#include <functional>
#include <random>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <functional>
#include "matrix.h"
#include "leaf.hpp"
#include "tree.hpp"

#include "partition.hpp"

#define max_recursions 7
#define dim 1

#define n_measure 50000
#define max_number_bins 100000

class diagMC
{
  
  std::mt19937 MyGenerator;
  std::uniform_real_distribution<> uni_dist;
  std::function<double()> rnd;

  partition possible_branches;
    
  enum UPDATE {TO_FAKE,INTEGRATION_VARIABLE,RANDOM_TREE,EXPAND_TREE,SHRINK_TREE,COUNT};
  int update_count_tree;
  
  enum class SECTORS {FAKE,TREE};
  SECTORS sector;
  
  matrix<double> update_statistic;
  
  std::vector<double> x;
  int max_order;
  tree HAM_tree;
  tree sub_tree;
  
  double current_weight;
  
  long long int fake_normalization;
  long long int MC_sum;
  long long int number_updates_tree;
  
  std::vector<long long int> diagram_order;
  std::vector<long long int> MC_sum_in_given_order;
  std::vector<std::vector<double>> result_order_bins;
  int number_measurements;
  
public:
  
  int number_bins;
  
  diagMC():
  MyGenerator(25),
  uni_dist(0,1),
  rnd(bind(uni_dist,MyGenerator)),
  possible_branches(max_recursions,2,std::ref(rnd)),
  update_count_tree(static_cast<int>(UPDATE::COUNT)),
  sector(SECTORS::TREE),
  update_statistic(4+1,update_count_tree+1),
  x(dim,0.5),
  max_order(max_recursions),
  HAM_tree(std::ref(rnd),possible_branches),
  sub_tree(std::ref(rnd),possible_branches),
  current_weight(0),
  fake_normalization(0),
  MC_sum(0),
  number_updates_tree(0),
  diagram_order(max_recursions+1,0),
  MC_sum_in_given_order(max_recursions+1,0),
  result_order_bins(max_number_bins,std::vector<double>(max_recursions+1,0)),
  number_measurements(0),
  number_bins(0)
  {};
  
  void initialize();
  
  void dostep();
  void dostep_FAKE();
  void dostep_TREE();
  
  int CHANGE_TO_TREE();
  
  int CHANGE_TO_FAKE();
  int CHANGE_INTEGRATION_VARIABLE();
  int CHANGE_RANDOM_TREE();
  int SHRINK_TREE_ORDER();
  int EXPAND_TREE_ORDER();
  
  double get_average_sign()
  {
    return static_cast<double>(MC_sum)/static_cast<double>(number_updates_tree);
  };
  
  double get_result()
  {
    return static_cast<double>(MC_sum)/static_cast<double>(fake_normalization)*fake_function(x);
  };
  
  double fake_function(std::vector<double>& arg)
  {
    return 1e-2*log(1+arg[0]);
  };
  
  void estimator();
  
  int sign() { return rnd()<0.5 ? 1 : -1;}
  
  void output_update_statistics() const
  {
    std::cout << "# printing update statistic " << std::endl;
    
    std::ostringstream filename;
    filename << "./out/update_statistic.txt";
    
    std::ofstream os(filename.str().c_str());
    std::vector<std::string> name;
    os << std::setprecision(8) << std::scientific;
    name.push_back("GO_TO_FAKE				");
    name.push_back("CHANGE_INTEGRATION_VARIABLE	");
    name.push_back("CHANGE_RANDOM_TREE			");
    name.push_back("EXPAND_TREE_ORDER			");
    name.push_back("SHRINK_TREE_ORDER			");
    name.push_back("GO_TO_TREE				");

    
    os << "# UPDATE STATISTICS";
    os << "\n# " << "# col 1 : selected"
      << "\n# " << "# col 2 : not possible"
      << "\n# " << "# col 3 : accepted"
      << "\n# " << "# col 4 : rejected"
      << "\n# " << "# col 5 : acceptance factor with respect to Metropolis ratio only" 
      <<  "\n";
    for (int i = 0; i <= update_count_tree; i++) 
    {
      os << "\n# " << name[i] << "\t" << update_statistic(0,i) << "\t" << update_statistic(1,i) << "\t"
	<< update_statistic(2,i) << "\t" << update_statistic(3,i) << "\t" << std::fixed << update_statistic(2,i) / (update_statistic(0,i) - update_statistic(1,i)) << std::scientific ;
    }
    os.close();
  }
  
  void output_diagram_order_statistics() const
  {
    std::cout << "# printing diagram order statistic " << std::endl;
    long long int sum=0;
    for(auto iter : diagram_order)
    {
     sum+=iter; 
    }
    if(sum!=number_updates_tree)
    {
      std::cout << "Error: Where are the missing diagrams " << number_updates_tree << std::endl;
     exit(1);
    }
    std::ostringstream filename;
    filename << "./out/diagram_order_statistic.txt";
    
    std::ofstream os(filename.str().c_str());
    os << std::setprecision(8) << std::scientific;
    for(size_t i=0;i<diagram_order.size();i++)
    {
     os << i << " " << static_cast<double>(diagram_order[i])/static_cast<double>(number_updates_tree) << std::endl; 
    }
    os.close();
  }
  
  void output_result_in_each_order()
  {
     std::cout << "# printing result in each order " << std::endl;

    std::ostringstream filename;
    filename << "./out/result_in_each_order.txt";
    
    std::ofstream os(filename.str().c_str());
    os << std::setprecision(10) ;
    for(size_t i=0;i<diagram_order.size();i++)
    {
      double mean = 0;
      double mean2 = 0;
      for(size_t j=0; j<number_bins; j++)
      {
        mean += result_order_bins[j][i];
        mean2 += result_order_bins[j][i]*result_order_bins[j][i];
      }
      mean /= static_cast<double>(number_bins);
      mean2 /= static_cast<double>(number_bins);
      os << i << " " << mean << " " << sqrt((mean2-mean*mean)/sqrt(number_bins)) << std::endl;
    }
    os.close();
  }
  
};
#endif
