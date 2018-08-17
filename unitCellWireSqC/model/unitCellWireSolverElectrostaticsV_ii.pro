/* -------------------------------------------------------------------
   File "mmgas_solver_mmgas_solver_electrostatics_v_ii_ii.pro"

   Electrostatics - Electric scalar potential v formulation
   ------------------------------------------------------------------- 

   I N P U T
   ---------

   Global Groups :  (Extension '_Ele' is for Electric problem)
   -------------
   Domain_Ele               Whole electric domain (not used)
   DomainCC_Ele             Nonconducting regions
   DomainC_Ele              Conducting regions (not used)

   Function :
   --------
   epsr[]                   Relative permittivity

   Constraint :
   ----------
   ElectricScalarPotential  Fixed electric scalar potential
                            (classical boundary condition)

   Physical constants :
   ------------------                                               */

   eps0 = 8.854187818e-12;

Group {
  DefineGroup[ Domain_Ele, DomainCC_Ele, DomainC_Ele ];
}

Function {
  DefineFunction[ epsr ];
}

FunctionSpace {
  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      // v = v  s   ,  for all nodes
      //      n  n
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support DomainCC_Ele; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; 
        NameOfConstraint ElectricScalarPotential; }
    }
  }
}


Formulation {
  { Name mmgas_solver_electrostatics_v_ii; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
    }
    Equation {
      Galerkin { [ epsr[] * Dof{d v} , {d v} ]; In DomainCC_Ele; 
                 Jacobian Vol; Integration GradGrad; }
    }
  }
}


Resolution {
  { Name mmgas_str_wire_v;
    System {
      { Name Sys_Ele; NameOfFormulation mmgas_solver_electrostatics_v_ii; }
    }
    Operation { 
      Generate[Sys_Ele]; Solve[Sys_Ele]; SaveSolution[Sys_Ele];
    }
  }
}


PostProcessing {
  { Name mmgas_str_wire_v; NameOfFormulation mmgas_solver_electrostatics_v_ii;
    Quantity {
      { Name v; 
        Value { 
          Local { [ {v} ]; In DomainCC_Ele; Jacobian Vol; } 
        }
      }
      { Name e; 
        Value { 
          Local { [ -{d v} ]; In DomainCC_Ele; Jacobian Vol; }
        }
      }
      { Name d; 
        Value { 
          Local { [ -eps0*epsr[] * {d v} ]; In DomainCC_Ele; 
                                             Jacobian Vol; }  
        } 
      }
    }
  }
}
