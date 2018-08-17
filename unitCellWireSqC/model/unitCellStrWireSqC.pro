//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/// MMGAS_STR_WIRE SOLVER MODULE

/* -------------------------------------------------------------------
   File "mmgas_str_wire.pro"

   This file defines the problem dependent data structures for the
   microstrip problem.
   
   To compute the solution: 
       getdp mmgas_str_wire -solve mmgas_str_wire_v

   To compute post-results:
       getdp mmgas_str_wire -pos Map
    or getdp mmgas_str_wire -pos Cut
   ------------------------------------------------------------------- */

Group {

SkinDomainC_Ele += Region[850000];					// upper electrode mesh cathode top boundary surface - dirichlet boundary
SkinDomainC_Ele += Region[870000];					// copper plate surface area - bottom electrode - dirichlet boundary

SkinDomainC_Ele += Region[750000];					// steel wire surface area - mesh surface area - dirichlet boundary
SkinDomainC_Ele += Region[760000];					// steel wire surface area - mesh surface area - dirichlet boundary
SkinDomainC_Ele += Region[770000];					// steel wire surface area - mesh surface area - dirichlet boundary
SkinDomainC_Ele += Region[780000];					// steel wire surface area - mesh surface area - dirichlet boundary

DomainCC_Ele += Region[960000]; 					// gas volume - conducting volume
DomainCC_Ele += Region[990000];						// copper plate volume - plate volume 
DomainCC_Ele += Region[900000];						// steel wire volume - mesh volume
DomainCC_Ele += Region[910000];						// steel wire volume - mesh volume 
DomainCC_Ele += Region[920000];						// steel wire volume - mesh volume 
DomainCC_Ele += Region[930000];						// steel wire volume - mesh volume

DomainQ_Ele += Region[960000]; 						// gas volume - conducting volume

Domain_Inf += Region[700000];						// periodic boundary surface - side 1 - continuous boundary
Domain_Inf += Region[710000];						// periodic boundary surface - side 2 - continuous boundary
Domain_Inf += Region[720000];						// periodic boundary surface - side 3 - continuous boundary
Domain_Inf += Region[730000];						// periodic boundary surface - side 4 - continuous boundary

}

// Include "mmgas_solver_options.pro";
Include "mmgas_solver_add_material_database.pro";

Function {

epsr[Region[960000]] = Air_epsilonr;
epsr[Region[990000]] = Copper_epsilonr;
epsr[Region[900000]] = SteelInd_epsilonr;
epsr[Region[910000]] = SteelInd_epsilonr;
epsr[Region[920000]] = SteelInd_epsilonr;
epsr[Region[930000]] = SteelInd_epsilonr;

}

Constraint { 

{ Name ElectricScalarPotential;

	Case { 

		{ Region Region[850000]; Value -2000; } 
		{ Region Region[870000]; Value 0; } 
		{ Region Region[750000]; Value -400; } 
		{ Region Region[760000]; Value -400; } 
		{ Region Region[770000]; Value -400; } 
		{ Region Region[780000]; Value -400; }

		} 
} 
}

Constraint { 

{ Name GlobalElectricPotential;
	
	Case { } 
} 
}

Constraint { 

{ Name GlobalElectricCharge;

	Case { } 
} 
}


Include "mmgas_solver_jacobian_lib.pro";
Include "mmgas_solver_integration_lib.pro";
// Include "mmgas_solver_electrostatics_v.pro";
Include "mmgas_solver_electrostatics_v_ii.pro";

/* Finally, we can define some operations to output results */

e = 1.e-7;

PostOperation {
  { Name Map; NameOfPostProcessing mmgas_str_wire_v;
     Operation {
       Print [ v, OnElementsOf DomainCC_Ele, File "mmgas_str_wire_v.pos" ];
       Print [ e, OnElementsOf DomainCC_Ele, File "mmgas_str_wire_e.pos" ];
     }
  }
  { Name Cut; NameOfPostProcessing mmgas_str_wire_v;
     Operation {
       Print [ e, OnLine {{e,e,0}{10.e-3,e,0}} {500}, File "Cut_e" ];
     }
  }

}