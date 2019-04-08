#include "diagMC.hpp"

int main()
{ 
  diagMC HAM;
  
  HAM.initialize();
  
  long long int number_updates = static_cast<long long int>(1e10);
  
  for(long long int i=0;i<number_updates;i++)
  {
    HAM.dostep();
    if (i%10 ==0 && i!=0)
      HAM.estimator();
    if(i%static_cast<int>(number_updates/500)==0 && i!=0)
    {
      std::cout << std::endl << "-------------------------------" << std::endl;
      std::cout << "number mc_steps already executed " << i << " ----> " << HAM.number_bins << " collected bins" << std::endl;
      HAM.output_update_statistics();
      HAM.output_diagram_order_statistics();
      HAM.output_result_in_each_order();
    }
  }
  
  HAM.output_update_statistics();
  HAM.output_diagram_order_statistics();
  HAM.output_result_in_each_order();
  std::cout << "The average sign of the integration is " << HAM.get_average_sign() << std::endl;
  std::cout << "The result of the integration is " << HAM.get_result() << std::endl;

  return 0; 
}
