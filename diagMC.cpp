#include "diagMC.hpp"
#include <functional>

void diagMC::initialize()
{

  HAM_tree.make_new_random_tree(2,x);
  current_weight=HAM_tree.get_weight();
  
  std::cout << "#################################################" << std::endl;
  std::cout << " Diagrammatic Monte Carlo sampling for the local" << std::endl;
  std::cout << " solution of the integral equation: " << std::endl;
  std::cout << " f(x) = c(x) + int_0^1 d^Dy K(x,y) f(y)^2" << std::endl;
  std::cout << "#################################################" << std::endl;
  
}


void diagMC::dostep()
{
  if(sector==SECTORS::FAKE)
    dostep_FAKE();
  else
    dostep_TREE();
}

void diagMC::dostep_FAKE()
{
  int c=CHANGE_TO_TREE();
  
  update_statistic(0,update_count_tree)++;
  update_statistic(c,update_count_tree)++;
}

int diagMC::CHANGE_TO_TREE()
{
  double suggestion_fake_tree=1.;
  double suggestion_tree_fake=1./static_cast<double>(static_cast<double>(UPDATE::COUNT));
  
  double new_weight=current_weight;
  double old_weight=fake_function(x);
  
  double acceptance_ratio=fabs(suggestion_tree_fake/suggestion_fake_tree*new_weight/old_weight);
  
  if(rnd()<acceptance_ratio)
  {
   sector=SECTORS::TREE;
   return 2;
  }
  else
    return 3;
}

void diagMC::dostep_TREE()
{
  double r=rnd();
  int c=0;
   
  if(r<1./static_cast<double>(UPDATE::COUNT))
  {
    c=CHANGE_TO_FAKE();
    update_statistic(0,UPDATE::TO_FAKE)++;
    update_statistic(c,UPDATE::TO_FAKE)++;
  }
  
  else if(r<2./static_cast<double>(UPDATE::COUNT))
  {
    c=CHANGE_INTEGRATION_VARIABLE();
    update_statistic(0,UPDATE::INTEGRATION_VARIABLE)++;
    update_statistic(c,UPDATE::INTEGRATION_VARIABLE)++;
  }
  
  else if(r<3./static_cast<double>(UPDATE::COUNT))
  {
    c=CHANGE_RANDOM_TREE();
    update_statistic(0,UPDATE::RANDOM_TREE)++;
    update_statistic(c,UPDATE::RANDOM_TREE)++;   
  }
  
  else if(r<4./static_cast<double>(UPDATE::COUNT))
  {
    c=EXPAND_TREE_ORDER();
    update_statistic(0,UPDATE::EXPAND_TREE)++;
    update_statistic(c,UPDATE::EXPAND_TREE)++;   
  }

  else if(r<5./static_cast<double>(UPDATE::COUNT))
  {
    c=SHRINK_TREE_ORDER();
    update_statistic(0,UPDATE::SHRINK_TREE)++;
    update_statistic(c,UPDATE::SHRINK_TREE)++;   
  }
}

int diagMC::CHANGE_TO_FAKE()
{
  if(HAM_tree.get_root().get_m()!=2) return 1;
  
  double suggestion_fake_tree=1.;
  double suggestion_tree_fake=1./static_cast<double>(static_cast<double>(UPDATE::COUNT));
  
  double old_weight=current_weight;
  double new_weight=fake_function(x);
  
  double acceptance_ratio=fabs(suggestion_fake_tree/suggestion_tree_fake*new_weight/old_weight);
  
  if(rnd()<acceptance_ratio)
  {
    sector=SECTORS::FAKE; 
    return 2;
  }
  else
    return 3;
}


int diagMC::CHANGE_INTEGRATION_VARIABLE()
{
  std::vector<int> leafs_to_update=HAM_tree.get_random_leaf();
  
  if(leafs_to_update[0]==0) return 1;
  
  std::vector<double> old_int_var=HAM_tree.get_leaf(leafs_to_update[0]).get_integration_variable();
  std::vector<double> new_int_var(old_int_var.size(),0);
  for(size_t i=0;i<new_int_var.size();i++)
    new_int_var[i]=old_int_var[i]+sign()*rnd()*0.05;
  
  for(auto iter: new_int_var)
    if(iter<0 || iter>1) return 1;
  
  for(auto& iter: leafs_to_update)
   HAM_tree.get_leaf(iter).set_integration_variable(new_int_var); 
  
  double new_weight=HAM_tree.get_weight();
  double acceptance_ratio=fabs(new_weight/current_weight);
  
  if(rnd()<acceptance_ratio)
  {    
    current_weight=new_weight;
    return 2;
  }
  else
  {
    for(auto iter: leafs_to_update)
      HAM_tree.get_leaf(iter).set_integration_variable(old_int_var);
    return 3;
  }
}


int diagMC::CHANGE_RANDOM_TREE()
{
  leaf& root_of_subtree=HAM_tree.get_root();
  
  if(root_of_subtree.get_m()<=2)
    return 1;
  if(root_of_subtree.get_m()>max_recursions)
    return 1;
  
  std::vector<double> root_var(dim,0);
  int number_integration_variables_old_subtree = 0;
  std::list<int> namelist_subtree;
  double prob_to_gen_old_subtree=1;
  double weight_of_old_subtree=1;
  int number_elements_old_subtree=0;
  HAM_tree.get_subtree_facts(root_of_subtree,root_var,number_integration_variables_old_subtree,namelist_subtree,prob_to_gen_old_subtree,weight_of_old_subtree,number_elements_old_subtree);
  
  if(fabs(weight_of_old_subtree-current_weight) > 1e-10)
  {
    std::cout << "Error: Weight of old_subtree not correct " << weight_of_old_subtree << " " << current_weight  << std::endl;
    exit(1);
  }
  
  double prob_to_gen_new_subtree=1;
  sub_tree.make_random_subtree(root_of_subtree.get_m(),root_var,prob_to_gen_new_subtree);
  double weight_of_new_subtree=sub_tree.get_weight();
  double number_elements_merged_tree=HAM_tree.get_element_count()-number_elements_old_subtree+sub_tree.get_element_count();
  
  double suggestion_old_new = prob_to_gen_new_subtree;
  double suggestion_new_old = prob_to_gen_old_subtree;
  
  double acceptance_ratio=fabs( suggestion_new_old/suggestion_old_new * weight_of_new_subtree/weight_of_old_subtree );
  
  if(rnd()<acceptance_ratio)
  {
   HAM_tree.merge_subtree(sub_tree, number_integration_variables_old_subtree, namelist_subtree);
   current_weight*=(weight_of_new_subtree/weight_of_old_subtree);
   return 2;
  }
  else
   return 3;
}


int diagMC::SHRINK_TREE_ORDER()
{
  // due to detailed balance only a root in k=1 config can be reduced
  
  if(HAM_tree.get_root().get_k()!=1) return 1;
  
  // a root representing u_2 can't be reduced further
  
  if(HAM_tree.get_root().get_m()<=2) return 1;
  
  std::vector<int> next_leaf;
  HAM_tree.get_next_leaves(HAM_tree.get_root(),next_leaf);
  
  const std::vector<double>& drop_int_var=HAM_tree.get_leaf(next_leaf[0]).get_integration_variable();
  
  // getting int variable of the branch of the next leaf - this is only possible if the tree size>=2
  double new_weight_old_weight;
  if(HAM_tree.get_root().get_m()==2)
  {
    new_weight_old_weight=HAM_tree.u(1,x)/(HAM_tree.kernel(x,drop_int_var)*2*HAM_tree.u(0,drop_int_var)*HAM_tree.u(1,drop_int_var));
    
  }
  else
  { 
    // next to next int variable
    
    const std::vector<double>& int_var=HAM_tree.get_next_int_var(next_leaf[0]);
    
    new_weight_old_weight=1./(HAM_tree.kernel(x,drop_int_var)*2*HAM_tree.u(0,drop_int_var)*HAM_tree.kernel(drop_int_var,int_var)/HAM_tree.kernel(x,int_var));
  }
  double suggestion_new_old=(1.-0);
  double suggestion_old_new=1;
  
  double acceptance_ratio=fabs(suggestion_new_old/suggestion_old_new*new_weight_old_weight);

  if(rnd()<acceptance_ratio)
  {
    HAM_tree.shrink();
    current_weight *= new_weight_old_weight;
    return 2;
  }
  else
    return 3;
}

int diagMC::EXPAND_TREE_ORDER()
{
  // appending to the current tree a k=1 branch of higher order -  can only be done if root_m<max_recursion
  
  if(HAM_tree.get_root().get_m()>=max_order) return 1;
  
  std::vector<double> added_int_var(dim,0);
  added_int_var[0]=rnd();
  double new_weight_old_weight;
  
  if(HAM_tree.get_root().get_m()!=1)
  {
    // update is also changing the int-kernel of the branch growing from the old root

    const std::vector<double>& int_var=HAM_tree.get_next_int_var(0);
    new_weight_old_weight=HAM_tree.kernel(x,added_int_var)*2*HAM_tree.u(0,added_int_var)*HAM_tree.kernel(added_int_var,int_var)/HAM_tree.kernel(x,int_var);
  }
  else
  {
    std::cout << "Error: configuration not possible" << std::endl;
    exit(1);
  }
  double suggestion_old_new=(1.-0);
  double suggestion_new_old=1;
  
  double acceptance_ratio=fabs(suggestion_new_old/suggestion_old_new*new_weight_old_weight);

  if(rnd()<acceptance_ratio)
  {
    
   HAM_tree.append(added_int_var);
   current_weight*=new_weight_old_weight;
   return 2;
    
  }
  else
    return 3;
  
}


void diagMC::estimator()
{

  if(sector==SECTORS::FAKE)
    fake_normalization++;
  if(sector==SECTORS::TREE)
  {
    number_updates_tree++;
    MC_sum+=current_weight/fabs(current_weight);
    diagram_order[HAM_tree.get_root().get_m()]++;
    MC_sum_in_given_order[HAM_tree.get_root().get_m()]+=current_weight/fabs(current_weight);
  }
  number_measurements++;
  if(number_measurements == n_measure)
  {
    for(size_t i=0;i<diagram_order.size();i++)
    {
      result_order_bins[number_bins][i] = static_cast<double>(MC_sum_in_given_order[i])/static_cast<double>(fake_normalization)*fake_function(x);
      MC_sum_in_given_order[i] = 0;
    }
    fake_normalization = 0;
    number_measurements = 0;
    number_bins++;
  }
}
