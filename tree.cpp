#include "tree.hpp"

void tree::make_new_random_tree(int given_m, const std::vector<double>& root_variable)
{
  mContainer_tree.clear();
  number_integration_variables=0;
  element_count=0;

  int current_heigth=given_m;
  number_integration_variables=0;
   
  // initialize the root of the tree; k is still undermined
  leaf root_leaf(current_heigth,-1);
  root_leaf.set_name(namelist[element_count]);
  root_leaf.set_prev_name(-1);
  // set the root to leaf to grow a branch
  root_leaf.set_leaf();
  // seed a external variable
  root_leaf.set_integration_variable(root_variable);
  container_pos[0]=element_count;
  mContainer_tree[element_count]=root_leaf;
  element_count++;
  active_leaf_list.push_back(0);
  
  while(active_leaf_list.size()!=0)
  {
    auto it = active_leaf_list.begin(); 
    grow_branch(*it);
  }
  
  if(active_leaf_list.size()!=0)
  {
    std::cout << "Error: You want to interrupt the growth of the tree; there are stil active leafs" << std::endl;
    exit(1);
  }
}

void tree::grow_branch(key root_key)
{
  leaf& root_leaf=mContainer_tree[container_pos[root_key]];
  
  if(root_leaf.get_m()<=1)
  {
   root_leaf.set_leaf();
   // this leaf is not a active leaf
   active_leaf_list.pop_back();
   
   return;
  }
  
  // grow from the root a branch with; for m_root=2 there is only B_1,1
  
  int k_branch;
  
  if((root_leaf.get_m()-1)<=1)
    k_branch=1;
  else
     k_branch=rnd()<0.5 ? 1 : 2;
  
  root_leaf.set_k(k_branch);
 
  // get the (m_1,...,m_n) of the new leafs
  
  std::pair<std::vector<int>,double> rnd_mon=possible_branches.get_random_monomial(root_leaf.get_m()-1,root_leaf.get_k());
  std::vector<int> branch_type=rnd_mon.first;
  std::vector<int> v={2,2};

  double bell_factor2=sqrt(rnd_mon.second);
  
  // add the new leafs according to the branch type to the tree
  
  std::vector<double> new_integratio_variable(dim,0);
  for(int i=0;i<dim;i++)
    new_integratio_variable[i]=rnd();
  number_integration_variables++;
  
  root_leaf.set_root();
  // this leaf is not a active leaf
  active_leaf_list.pop_front();
  
  for(size_t i=0; i<branch_type.size(); i++)
  {
   leaf new_leaf(branch_type[i],-1);
   new_leaf.set_name(namelist[element_count]);
   new_leaf.set_prev_name(root_leaf.get_name());
   root_leaf.set_next_name(new_leaf.get_name());
   new_leaf.set_leaf();
   new_leaf.set_integration_variable(new_integratio_variable);
   new_leaf.set_bell_2(bell_factor2);
   container_pos[new_leaf.get_name()]=element_count;
   mContainer_tree[element_count]=new_leaf;
   element_count++;
   active_leaf_list.push_front(new_leaf.get_name());
    
  }
}

void tree::make_random_subtree(int given_m, const std::vector<double>& root_variable, double& gen_prob)
{
  int current_heigth=given_m;
  number_integration_variables=0;
  element_count=0;
  
  // initialize the root of the tree; k is still undermined
  leaf root_leaf(current_heigth,-1);
  root_leaf.set_name(0);
  root_leaf.set_prev_name(-1);
  // set the root to leaf to grow a branch
  root_leaf.set_leaf();
  // seed a external variable
  root_leaf.set_integration_variable(root_variable);
  container_pos[0]=element_count;
  mContainer_tree[element_count]=root_leaf;
  element_count++;
  active_leaf_list.push_front(0);
  
  while(active_leaf_list.size()!=0)
  {   
   auto it = active_leaf_list.begin();
   grow_branch_sub(*it,gen_prob);

  }
  
  if(active_leaf_list.size()!=0)
  {
    std::cout << "Error:You want to interupt the growth of the sub tree; there are stil active leafs" << std::endl;
    exit(1);
  }
}
 
void tree::grow_branch_sub(key root_key, double& gen_prob)
{
  leaf& root_leaf=mContainer_tree[container_pos[root_key]];
  
  if(root_leaf.get_m()<=1)
  {
   root_leaf.set_leaf();
   // this leaf is not a active leaf
   active_leaf_list.pop_front();
   
   return;
  }
  
  // grow from the root a branch with; for m_root=2 there is only B_1,1
  
  int k_branch;
  
  if((root_leaf.get_m()-1)<=1)
  {
    k_branch=1;
    gen_prob*=1;
  }
  else
  {
     k_branch=rnd()<0.5 ? 1 : 2;
     gen_prob *= 0.5;
  }
    
  root_leaf.set_k(k_branch);
  
  // get the (m_1,...,m_n) of the new leafs
  
  std::pair<std::vector<int>,double> rnd_mon=possible_branches.get_random_monomial(root_leaf.get_m()-1,root_leaf.get_k());
  std::vector<int> branch_type=rnd_mon.first;

  double bell_factor2=sqrt(rnd_mon.second);   
  gen_prob *= k_branch==1 ? 1 : 1./static_cast<double>(possible_branches.get_number_possible_mon(root_leaf.get_m()-1,root_leaf.get_k()));

   // set the integration variables
  std::vector<double> new_integratio_variable(dim,0);
  for(int i=0;i<dim;i++)
    new_integratio_variable[i]=rnd();
  number_integration_variables++;
  gen_prob *= 1./(1.-0.);    // the new variable is picked with (dt/intervalsize)
 
  root_leaf.set_root();
  // this leaf is not a active leaf
  active_leaf_list.pop_front();
  
  for(size_t i=0; i<branch_type.size(); i++)
  {
   leaf new_leaf(branch_type[i],-1);
   new_leaf.set_name(element_count);
   new_leaf.set_prev_name(root_leaf.get_name());
   root_leaf.set_next_name(new_leaf.get_name());
   new_leaf.set_leaf();
   new_leaf.set_integration_variable(new_integratio_variable);
   new_leaf.set_bell_2(bell_factor2);
    container_pos[new_leaf.get_name()]=element_count;
   mContainer_tree[element_count]=new_leaf;
   element_count++;
   active_leaf_list.push_front(new_leaf.get_name());
  }
}


double tree::get_weight()
{
  double weight=1;
  int root_key=0;
  
  DFS_FUNC_WEIGHT(root_key,weight);

  return weight;
}

double tree::get_sub_weight(leaf& from_this_root)
{
  double weight=1;
  int root_key=from_this_root.get_name();
  
  DFS_FUNC_WEIGHT(root_key,weight);

  return weight;
}
 
void tree::get_subtree_facts(leaf& sub_root, std::vector<double>& root_variable, int& number_integration_variables, std::list<int>& sub_namelist, double& sub_gen_prob, double& sub_weight, int& number_elements)
{
  root_variable=sub_root.get_integration_variable();
  
  DFS_FUNC_FACTS(sub_root, number_integration_variables, sub_namelist, sub_gen_prob, sub_weight, number_elements);
 }
 
void tree::DFS_FUNC_WEIGHT(key root_key, double& weight)
{
  int root_m=mContainer_tree[container_pos[root_key]].get_m();
  std::vector<double> root_integration_variable=mContainer_tree[container_pos[root_key]].get_integration_variable();
  if(root_m<=1) 
  {
    weight*=u(root_m,root_integration_variable);
    return;
  }
  int branch_k=mContainer_tree[container_pos[root_key]].get_k();
  std::vector<int>& n_n=mContainer_tree[container_pos[root_key]].get_next_names();
  weight*=1./factorial(root_m-1);
  std::vector<double> leaf_integration_variable=mContainer_tree[container_pos[n_n[0]]].get_integration_variable();

  weight*=2.;
  weight*=kernel(root_integration_variable,leaf_integration_variable);
  weight*=branch_k==1 ? u(0,leaf_integration_variable) : 1.;
  
  for(size_t i=0; i<n_n.size(); i++)
  {
    weight*=factorial(mContainer_tree[container_pos[n_n[i]]].get_m());
    weight*=mContainer_tree[container_pos[n_n[0]]].get_bell_2();
    DFS_FUNC_WEIGHT(n_n[i],weight);   
  }  
}

void tree::DFS_FUNC_FACTS(leaf& sub_root, int& number_integration_variables, std::list<int>& sub_namelist, double& sub_gen_prob, double& sub_weight, int& number_elements)
{
  int root_m=sub_root.get_m();
  std::vector<double> root_int_var=sub_root.get_integration_variable();
  int branch_k=sub_root.get_k();
  std::vector<int>& n_n=sub_root.get_next_names();
  
  // functor number_elements
  
  number_elements++;
  
  // functor name
  
  sub_namelist.push_back(sub_root.get_name());
  
  // functor weight
  
  if(root_m<=1) 
  {
      sub_weight*=u(root_m,root_int_var);
      // as root m<=1 there is no further branch from this sub_root
      return;
  }
  
  // now there is a branch
  
  sub_weight*=1./factorial(root_m-1);
  std::vector<double> leaf_integration_variable=mContainer_tree[container_pos[n_n[0]]].get_integration_variable();
  sub_weight*=2.;
  sub_weight*=kernel(root_int_var,leaf_integration_variable);
  sub_weight*=branch_k==1 ? u(0,leaf_integration_variable) : 1.;
  
  // func generation prob for branch
  
  // prob to choose k=1,2 - for u_2 there is only B_1_1 to choose
  sub_gen_prob *= root_m<=2 ? 1. : 0.5;
  // prob to choose monomial - B_n_1=u_n
  sub_gen_prob *= branch_k==1 ? 1 : 1./static_cast<double>(possible_branches.get_number_possible_mon(root_m-1,branch_k));
  // prob to pick integration variable
  sub_gen_prob *= 1./(1.-0);
  number_integration_variables++;
  
  for(size_t i=0; i<n_n.size(); i++)
  {
    sub_weight*=factorial(mContainer_tree[container_pos[n_n[i]]].get_m());
    sub_weight*=mContainer_tree[container_pos[n_n[0]]].get_bell_2();
    DFS_FUNC_FACTS(mContainer_tree[container_pos[n_n[i]]], number_integration_variables, sub_namelist, sub_gen_prob, sub_weight, number_elements);
  }  
  
  
}

std::vector<int> tree::get_random_leaf()
{
  // if there is only one element in the tree CHANGE_INTEGRATION_VARIABLE should return update not possible
  int r_key;
  if(element_count<=1)
  {
    std::vector<int> r_val;
    return r_val={0};
  }
  else
  {
    // the root should not be picked
    r_key=namelist[static_cast<int>(rnd()*(element_count-1))+1];
  }
  
  leaf& r_leaf=mContainer_tree[container_pos[r_key]];
  leaf& prev_leaf=mContainer_tree[container_pos[r_leaf.get_prev_name()]];
  
  if(prev_leaf.get_next_names().size()==1)
  {
    std::vector<int> r_val;
    return r_val={r_leaf.get_name()};
  }
  else
  {
   std::vector<int> r_val;
   for(auto iter : prev_leaf.get_next_names())
     r_val.push_back(mContainer_tree[container_pos[iter]].get_name());
   return r_val;
  }
  
  
  
}

leaf& tree::get_single_random_root()
{
  int r_key=namelist[static_cast<int>(rnd()*(element_count))];
  
  return mContainer_tree[container_pos[r_key]];
  
}

leaf& tree::get_prev_leaf(leaf& this_leaf)
{
    return mContainer_tree[container_pos[this_leaf.get_prev_name()]];
}

void tree::get_next_leaves(leaf& this_leaf, std::vector<key>& next_leaves)
{
    for(size_t i=0; i<this_leaf.get_next_names().size(); i++)
      next_leaves.push_back(this_leaf.get_next_names()[i]);
}

void tree::merge_subtree(tree& new_subtree, int int_variables_old_subtree, const std::list<int>& namelist_old_tree)
{
  number_integration_variables-=(int_variables_old_subtree-new_subtree.get_number_int_varialbes());
  
  std::list<int> names=namelist_old_tree;
  // translating_old_names_to_new_names[old_name]=new_name
  std::vector<int> translating_new_names_to_old_names(max_number_nodes,-1);
  // transfering prev name of old_subtree root
  new_subtree.get_root().set_prev_name(mContainer_tree[container_pos[*namelist_old_tree.begin()]].get_prev_name());
  // transfering bell_factor2 of old_subtree root to new_subtree root
  new_subtree.get_root().set_bell_2(mContainer_tree[container_pos[*namelist_old_tree.begin()]].get_bell_2());
    
  new_subtree.DFS_FUNC_NAMES(new_subtree.get_root(),names,element_count,translating_new_names_to_old_names);
  
  // insert the elements of new_subtree into the elements of this
  
  for(auto iter : namelist_old_tree)
    mContainer_tree[container_pos[iter]]=new_subtree.get_leaf(translating_new_names_to_old_names[iter]);
}

// this should be a FUNC on a subtree - number_elements is the element count of the merged tree
void tree::DFS_FUNC_NAMES(leaf& root, std::list<int>& names, int& number_elements, std::vector<int>& translating_new_names_to_old_names)
{
  key old_root_name=root.get_name();
  
  if(names.size()!=0)
  {
    root.set_name(*names.begin());
    translating_new_names_to_old_names[root.get_name()]=old_root_name;
    names.pop_front();
    
  }
  else
  {
    root.set_name(number_elements);
    translating_new_names_to_old_names[root.get_name()]=old_root_name;
    number_elements++;
  }
  // restoring next relation for the prev of the leaf - should not be done if the selected leaf is the root of the new_subtree
  if(old_root_name!=0)
  {
    std::vector<key>& next_names_of_prev=mContainer_tree[container_pos[translating_new_names_to_old_names[root.get_prev_name()]]].get_next_names();

    for(size_t i=0;i<next_names_of_prev.size();i++) 
    {
        if(next_names_of_prev[i]==old_root_name)
        {
            next_names_of_prev[i]=root.get_name();
            break;
        }
    }
  }
    
  if(root.get_m()<=1) // there is no further branch
    return;
    
  // Here comes the branche from the root
    
  std::vector<key> next_names_of_root=root.get_next_names();
  
  for(size_t i=0;i<next_names_of_root.size();i++)
  {
    // restoring the prev relation of the next leaf
    mContainer_tree[container_pos[next_names_of_root[i]]].set_prev_name(root.get_name());
    DFS_FUNC_NAMES(mContainer_tree[container_pos[next_names_of_root[i]]],names,number_elements,translating_new_names_to_old_names);
  }
}

void tree::append(std::vector<double>& new_int_var)
{
 leaf& old_root=get_root();
 
 // updating name of the old_root
 
 old_root.set_name(namelist[element_count]);
 
 // setting the prev of the old root to be the new_root
 
 old_root.set_prev_name(0);
 
 if(old_root.get_m()!=1)
 {
  // updating the prev relation of the next leafs
  std::vector<int>& n_n=old_root.get_next_names();
  
  for(size_t i=0;i<n_n.size();i++)
  {
    mContainer_tree[container_pos[n_n[i]]].set_prev_name(old_root.get_name());
    
  }
 }
 // insert the new root
 
    // getting the external_variable and set the new integration variable on the old root
    std::vector<double> ex_var=old_root.get_integration_variable();
    old_root.set_integration_variable(new_int_var);
    number_integration_variables++;
  
    // changing key of the old_root
    mContainer_tree[container_pos[old_root.get_name()]]=old_root;

  mContainer_tree[container_pos[0]].set_m(old_root.get_m()+1);
  mContainer_tree[container_pos[0]].set_k(1);
  mContainer_tree[container_pos[0]].reset_next_names(old_root.get_name());
  mContainer_tree[container_pos[0]].set_integration_variable(ex_var);
// the order may not be permutated due to the reuse of the reference to old_root
  mContainer_tree[container_pos[0]].set_name(0);
  mContainer_tree[container_pos[0]].set_prev_name(-1);
  element_count++;
}

void tree::shrink()
{
  // getting the leaf after the root
  leaf& new_root=mContainer_tree[container_pos[mContainer_tree[container_pos[0]].get_next_names()[0]]];
  
  // restore the prev relation of the next leafs after the new_root - must only be done if the new_root has m>=2
  if(new_root.get_m()!=1)
  {
    for(size_t i=0;i<new_root.get_next_names().size();i++)
    {
      mContainer_tree[container_pos[new_root.get_next_names()[i]]].set_prev_name(0);
    }
  }
  // params of the new root
  new_root.set_prev_name(-1);
  new_root.set_integration_variable(mContainer_tree[container_pos[0]].get_integration_variable());
  number_integration_variables--;
  //mContainer_tree.erase(0);
  key new_root_name=new_root.get_name();
  
  // dropping the name of the new_root out of the namelist; There must be another data structure !!!!!!!
  int pos_of_new_root_name_in_the_namelist=100;
  for(int i=0; i<element_count; i++)
  {
    if(namelist[i]==new_root_name) 
    {
      pos_of_new_root_name_in_the_namelist=i;
      break;
    }
  }
  for(int i=pos_of_new_root_name_in_the_namelist;i<element_count-1;i++)
  {
   namelist[i]=namelist[i+1]; 
  }
  element_count--;
  namelist[element_count]=new_root_name;
  
  // identify the new_root
  new_root.set_name(0);
  mContainer_tree[container_pos[0]]=new_root;
}
