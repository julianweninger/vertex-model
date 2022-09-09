
``ActiveSLS`` - A model for active standard linear solids
=========================================================

This model includes 
   1. Visco-elastic Maxwell element with stiffness k (`maxwell_stiffness`) and viscous relaxation timescale :math:`\tau_v` (`maxwell_tau_viscous`)
   1. An elastic spring with spring constant B (`spring_constant`) and rest lenght :math:`a` (`spring_restlength`) representing the elastic background.
   1. An active element the generates active tension

For more information see Sknepnek, R.; Djafer-Cherif, I.; Chuai, M.; Weijer, C. J.; Henkes, S. Generating Active T1 Transitions through Mechanochemical Feedback. arXiv April 3, 2022.
