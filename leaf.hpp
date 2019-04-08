#ifndef LEAF_HPP
#define LEAF_HPP

#include <vector>
#include <iostream>

class leaf
{
  int name;
  int prev_name;
  std::vector<int> next_name;
  
  bool is_root;
  bool is_leaf;
 
  int m;
  int k;
  double bell_factor_2;
  
  std::vector<double> integration_variable;
  
  
public:
  
  
  leaf(int order, int n_blocks): m(order), k(n_blocks), bell_factor_2(1.){};
  
  leaf& operator=(const leaf& rhs)
  {
    if(this != &rhs)
    {
      name=rhs.get_name();
      prev_name=rhs.get_prev_name();
      next_name= rhs.get_next_names();
      m=rhs.get_m();
      k=rhs.get_k();
      bell_factor_2=rhs.get_bell_2();
      integration_variable=rhs.get_integration_variable();
    }
    return *this;
  }
  
  void set_leaf() { is_root=false; is_leaf=true;}
  void set_root() { is_root=true; is_leaf=false;}
  void set_name(int n) { name=n;}
  void set_prev_name(int pn) { prev_name=pn;}
  void set_next_name(int nn) { next_name.push_back(nn);}
  void reset_next_names(int nn) {next_name.clear(); next_name.push_back(nn);}
  void set_m(int order) {m=order;}
  void set_k(int n_blocks) { k=n_blocks;}
  void set_bell_2(double factor) {bell_factor_2=factor;}
  void set_integration_variable(const std::vector<double>& new_variable) {integration_variable=new_variable;}
   
  bool is_this_leaf() { return is_leaf;}
  bool is_this_root() { return is_root;}
  
  int get_name() const { return name;}
  int get_prev_name() const { return prev_name;}
  const std::vector<int>& get_next_names() const { return next_name;}
  std::vector<int>& get_next_names() { return next_name;}
  int get_m() const { return m;}
  int get_k() const { return k;}
  double get_bell_2() const { return bell_factor_2;}
  const std::vector<double>& get_integration_variable() const { return integration_variable;}
  
  void print_next_names()
  {
   std::cout << "Printing next names of " << name << "; There are " << next_name.size() << " nodes: ";
   for(auto iter : next_name)
   {
     std::cout << iter << " ";
   }
   std::cout<< std::endl;
  }
  
};



#endif