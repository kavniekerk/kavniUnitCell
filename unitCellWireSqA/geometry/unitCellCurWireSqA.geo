Include "pipe_data.pro";

// *********************************************************************
// gemcell.geo
//
// Description:
// Geometry file for a GEM cell.
// This cell can be repeated any number of times within Garfield 
// to construct an arbitrarily large GEM.
//
// References: 
// 1. This specific form of GEM cell was found in 
//    "A How-to Approach for a 3d Simulation of Charge Transfer
//    Characteristics in a Gas Electron Multiplier (GEM)" by A. Sharma.
//    As of 04/08/10, this document can be found at:
//    www.slac.stanford.edu/pubs/icfa/fall99/paper2/paper2.pdf
//
// See also:
// 2. A. Sharma.  Nucl. Instr. Meth. A 454, 267-271 (2000).
// 3. O. Bouianov et al. Nucl. Instr. Meth. A 450, 277-287 (2000).
// 4. V. Tikhonov and R. Veenhof. Nucl. Instr. Meth. A 478, 452-459 (2002).
// 5. C. Shalem et al. Nucl. Instr. Meth. A, 558, 475–489 (2006).
//
// *********************************************************************

// Parameters

// pillar parameters
r0 = 0.01;                    // the pillar radius, in mm
pil_f_x = 0;                  // pillar co-ordinates, multiplication factor in x, 1.25
pil_f_y = 0;	    	          // pillar co-ordinates, multiplication factor in y, 1.25
a = 0.11;                     // the "pitch", or distance between GEM pillars, in mm
pil_c_x = -0.025*0 - 0/4;     // pillar co-ordinates, constant factor in x, -0.025
pil_c_y = -0.025*0 - 0/4;     // pillar co-ordinates, constant factor in y, -0.025

// vertical parameters
r1 = 0.005;                   // the etching amount (etch radius = r0 + r1), in mm
tC = 0.0035;                  // lower copper thickness, in mm
tD = 0.04;                    // dielectric thickness, in mm
tuC = 0.026;                  // higher copper thickness, in mm
lE = 0.5;                     // distance from GEM plates to upper exterior electrode, in mm
lP = 0.1;                     // distance from lower LEM plate to pad (readout) plane, in mm

// mesh window and wire parameteres
mwf = 1;									                      // mesh_window_factor
mm = 1;                                         // geometrical scaling
r_w = 0.0025 * mm;                              // radius of Wiremesh, in microns
p_0 = 0.025;                                    // pitch of the window, in mm
p = 0.025 * mm - 0*r_w/mwf * mm;                // pitch of the window, in microns
R = (p * p + r_w * r_w)/(2 * r_w);              // radius
alpha = Asin((p/R));                            // angle in radians
Total_Grid_size = (a - 0.01)/2;                 // total grid size, in mm, 0.4

Number_Wires = ((Total_Grid_size)/(p_0))/2;     // number of wires
Wire_length = Total_Grid_size / Number_Wires;   // wire length

Number_Units_x = 0;						                  // number of units, 1
Number_Units_y = 0; 					                  // number of units, 1

geo_wc_xr = 2*r_w;                              // y-direction wire in x radial direction
geo_wc_yr = 2*r_w;                              // x-direction wire in y radial direction

geo_wc_xd = 2*r_w;                              // x-direction wire in x-direction
geo_wc_yd = 2*r_w;                              // y-direction wire in y-direction

n_1 = 0;
m_1 = 0;
n_2 = 0;
m_2 = 0;

i_t_x = Number_Wires;
j_t_x = Number_Wires+1;

i_t_y = Number_Wires;
j_t_y = Number_Wires+1;

mesh_level = 0.0375;                            // mesh level, in mm
mesh_window = 0.05;                             // mesh window, in mm

// shell parameters
geo_f_x = 1;                                    // Geometric_factor
geo_f_y = 1;                                    // Geometric_factor

// Characteristic lengths

  lcDielectricpillar = 0.001;				  // characterization of mesh
  lcEtchingpillar = 0.001;
  lcCopperPlateBdry = 0.001;
  lcExtElectrodeBdry = 0.01;
  LcWiremesh = 0.01;                            
                    
// Extrusion Precision

Geometry.ExtrudeSplinePoints = 3;
Geometry.Points = 0;

// *********************************************************************

// Wire Mesh

// First set of wires

  // x-direction

For i In {1:i_t_x}
  For j In {1:j_t_x}

  k=1;
  l=1;

  sf_x1_1a1[] = {}; // surfaces of all Wiremeshs
  llf_x1_1a1_0[] = {}; // line loops of bottom Wiremesh intersects
  llf_x1_1a1_1[] = {}; // line loops of top Wiremesh intersects

      p0_x1_1a1~{j}~{i} = newp; Point(p0_x1_1a1~{j}~{i}) = {2*p+4*(i-1)*p+2*k*p+geo_wc_xd,4*(j-1)*p+geo_wc_yr,-r_w+mesh_level*mm, LcWiremesh * mm};
      p1_x1_1a1~{j}~{i} = newp; Point(p1_x1_1a1~{j}~{i}) = {2*p+4*(i-1)*p+2*k*p+geo_wc_xd,4*(j-1)*p+geo_wc_yr,-2*r_w+mesh_level*mm, LcWiremesh * mm};
      p2_x1_1a1~{j}~{i} = newp; Point(p2_x1_1a1~{j}~{i}) = {2*p+4*(i-1)*p+2*k*p+geo_wc_xd,4*(j-1)*p+r_w+geo_wc_yr,-r_w+mesh_level*mm, LcWiremesh * mm};
      p3_x1_1a1~{j}~{i} = newp; Point(p3_x1_1a1~{j}~{i}) = {2*p+4*(i-1)*p+2*k*p+geo_wc_xd,4*(j-1)*p+geo_wc_yr,0+mesh_level*mm, LcWiremesh * mm};
      p4_x1_1a1~{j}~{i} = newp; Point(p4_x1_1a1~{j}~{i}) = {2*p+4*(i-1)*p+2*k*p+geo_wc_xd,4*(j-1)*p-r_w+geo_wc_yr,-r_w+mesh_level*mm, LcWiremesh * mm};

      l1_x1_1a1~{j}~{i} = newl; Circle(l1_x1_1a1~{j}~{i}) = {p1_x1_1a1~{j}~{i}, p0_x1_1a1~{j}~{i}, p2_x1_1a1~{j}~{i}};
      l2_x1_1a1~{j}~{i} = newl; Circle(l2_x1_1a1~{j}~{i}) = {p2_x1_1a1~{j}~{i}, p0_x1_1a1~{j}~{i}, p3_x1_1a1~{j}~{i}};
      l3_x1_1a1~{j}~{i} = newl; Circle(l3_x1_1a1~{j}~{i}) = {p3_x1_1a1~{j}~{i}, p0_x1_1a1~{j}~{i}, p4_x1_1a1~{j}~{i}};
      l4_x1_1a1~{j}~{i} = newl; Circle(l4_x1_1a1~{j}~{i}) = {p4_x1_1a1~{j}~{i}, p0_x1_1a1~{j}~{i}, p1_x1_1a1~{j}~{i}};

      ll1_x1_1a1~{j}~{i} = newll; Line Loop(ll1_x1_1a1~{j}~{i}) = {l1_x1_1a1~{j}~{i}, l2_x1_1a1~{j}~{i}, l3_x1_1a1~{j}~{i}, l4_x1_1a1~{j}~{i}};
      s1_x1_1a1~{j}~{i} = news; Plane Surface(s1_x1_1a1~{j}~{i}) = {ll1_x1_1a1~{j}~{i}};
      llf_x1_1a1_0[] += ll1_x1_1a1~{j}~{i};

      st_s1_x1[] += s1_x1_1a1~{j}~{i};

        //phys_fil_x1_1a1[] = {};
        //phys_fil_top_x1_1a1 = {};
        //phys_fil_bot_x1_1a1 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_x1_1a1 (%g in layer %g)", j, i),
        //BND_Wiremesh1_x1_1a1 + 1000 * i + j) = {s1_x1_1a1~{j}~{i}}; // bottom
        //phys_fil_bot_x1_1a1 += BND_Wiremesh1_x1_1a1 + 1000 * i + j;

        v_x1_1a1[] = {};
        s_x1_1a1[] = {};
        tmp_x1_1a1[] = {s1_x1_1a1~{j}~{i}};
          tmp_x1_1a1[] = Extrude {{0,0,0}, {0,1,0}, {2*p+4*(i-1)*p+2*k*p,4*(j-1)*p,R-r_w}, alpha} {
            Surface{ tmp_x1_1a1[0] };
          };

          v_x1_1a1[] += tmp_x1_1a1[1];
          vt_x1_1a1[] += v_x1_1a1[];
          s_x1_1a1[] += tmp_x1_1a1[{2:5}];
          st_x1_1a1[] += s_x1_1a1[];

        //Physical Surface(Sprintf("Wiremesh top boundary_x1_1a1 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_x1_1a1 + 1100 * i + j) = tmp_x1_1a1[0]; // top
        //phys_fil_top_x1_1a1 += BND_Wiremesh2_x1_1a1 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_x1_1a1 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_x1_1a1 + 1200 * i + j) = s_x1_1a1[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_x1_1a1 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_x1_1a1 + 1000 * i + j) = v_x1_1a1[];
        //phys_fil_x1_1a1[] += VOL_Wiremesh1_x1_1a1 + 1000 * i + j;

        sf_x1_1a1[] += s_x1_1a1[];
        ll2_x1_1a1~{j}~{i} = newll; Line Loop(ll2_x1_1a1~{j}~{i}) = Boundary{ Surface{tmp_x1_1a1[0]}; };
        llf_x1_1a1_1[] += ll2_x1_1a1~{j}~{i};

  sf_x1_1a2[] = {}; // surfaces of all Wiremeshs
  llf_x1_1a2_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_x1_1a2[] = {};
        //phys_fil_top_x1_1a2 = {};
        //phys_fil_bot_x1_1a2 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_x1_1a2 (%g in layer %g)", j, i),
        //BND_Wiremesh1_x1_1a2 + 1000 * i + j) = {tmp_x1_1a1[0]}; // bottom
        //phys_fil_bot_x1_1a2 += BND_Wiremesh1_x1_1a2 + 1000 * i + j;

        v_x1_1a2[] = {};
        s_x1_1a2[] = {};
        tmp_x1_1a2[] = {tmp_x1_1a1[0]};
          tmp_x1_1a2[] = Extrude {{0,0,0},{0,-1,0},{4*(i-1)*p+2*k*p,4*(j-1)*p,-R+r_w}, alpha} {
            Surface{ tmp_x1_1a2[0] };
          };

          v_x1_1a2[] += tmp_x1_1a2[1];
          vt_x1_1a2[] += v_x1_1a2[];
          s_x1_1a2[] += tmp_x1_1a2[{2:5}];
          st_x1_1a2[] += s_x1_1a2[];

        //Physical Surface(Sprintf("Wiremesh top boundary_x1_1a2 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_x1_1a2 + 1100 * i + j) = tmp_x1_1a2[0]; // top
        //phys_fil_top_x1_1a2 += BND_Wiremesh2_x1_1a2 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_x1_1a2 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_x1_1a2 + 1200 * i + j) = s_x1_1a2[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_x1_1a2 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_x1_1a2 + 1000 * i + j) = v_x1_1a2[];
        //phys_fil_x1_1a2[] += VOL_Wiremesh1_x1_1a2 + 1000 * i + j;

        sf_x1_1a2[] += s_x1_1a2[];
        ll2_x1_1a2~{j}~{i} = newll; Line Loop(ll2_x1_1a2~{j}~{i}) = Boundary{ Surface{tmp_x1_1a2[0]}; };
        llf_x1_1a2_1[] += ll2_x1_1a2~{j}~{i};

  sf_x1_1b1[] = {}; // surfaces of all Wiremeshs
  llf_x1_1b1_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_x1_1b1[] = {};
        //phys_fil_top_x1_1b1 = {};
        //phys_fil_bot_x1_1b1 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_x1_1b1 (%g in layer %g)", j, i),
        //BND_Wiremesh1_x1_1b1 + 1000 * i + j) = {tmp_x1_1a2[0]}; // bottom
        //phys_fil_bot_x1_1b1 += BND_Wiremesh1_x1_1b1 + 1000 * i + j;

        v_x1_1b1[] = {};
        s_x1_1b1[] = {};
        tmp_x1_1b1[] = {tmp_x1_1a2[0]};
          tmp_x1_1b1[] = Extrude {{0,0,0},{0,-1,0},{4*(i-1)*p+2*k*p,2*p+4*(j-1)*p,-R+r_w}, alpha} {
            Surface{ tmp_x1_1b1[0] };
          };

          v_x1_1b1[] += tmp_x1_1b1[1];
          vt_x1_1b1[] += v_x1_1b1[];
          s_x1_1b1[] += tmp_x1_1b1[{2:5}];
          st_x1_1b1[] += s_x1_1b1[];

        //Physical Surface(Sprintf("Wiremesh top boundary_x1_1b1 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_x1_1b1 + 1100 * i + j) = tmp_x1_1b1[0]; // top
        //phys_fil_top_x1_1b1 += BND_Wiremesh2_x1_1b1 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_x1_1b1 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_x1_1b1 + 1200 * i + j) = s_x1_1b1[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_x1_1b1 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_x1_1b1 + 1000 * i + j) = v_x1_1b1[];
        //phys_fil_x1_1b1[] += VOL_Wiremesh1_x1_1b1 + 1000 * i + j;

        sf_x1_1b1[] += s_x1_1b1[];
        ll2_x1_1b1~{j}~{i} = newll; Line Loop(ll2_x1_1b1~{j}~{i}) = Boundary{ Surface{tmp_x1_1b1[0]}; };
        llf_x1_1b1_1[] += ll2_x1_1b1~{j}~{i};

  sf_x1_1b2[] = {}; // surfaces of all Wiremeshs
  llf_x1_1b2_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_x1_1b2[] = {};
        //phys_fil_top_x1_1b2 = {};
        //phys_fil_bot_x1_1b2 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_x1_1b2 (%g in layer %g)", j, i),
        //BND_Wiremesh1_x1_1b2 + 1000 * i + j) = {tmp_x1_1b1[0]}; // bottom
        //phys_fil_bot_x1_1b2 += BND_Wiremesh1_x1_1b2 + 1000 * i + j;

        v_x1_1b2[] = {};
        s_x1_1b2[] = {};
        tmp_x1_1b2[] = {tmp_x1_1b1[0]};
          tmp_x1_1b2[] = Extrude {{0,0,0},{0,1,0},{-2*p+4*(i-1)*p+2*k*p,4*(j-1)*p,R-r_w}, alpha} {
            Surface{ tmp_x1_1b2[0] };
          };

          v_x1_1b2[] += tmp_x1_1b2[1];
          vt_x1_1b2[] += v_x1_1b2[];
          s_x1_1b2[] += tmp_x1_1b2[{2:5}];
          st_x1_1b2[] += s_x1_1b2[];

        //Physical Surface(Sprintf("Wiremesh top boundary_x1_1b2 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_x1_1b2 + 1100 * i + j) = tmp_x1_1b2[0]; // top
        //phys_fil_top_x1_1b2 += BND_Wiremesh2_x1_1b2 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_x1_1b2 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_x1_1b2 + 1200 * i + j) = s_x1_1b2[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_x1_1b2 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_x1_1b2 + 1000 * i + j) = v_x1_1b2[];
        //phys_fil_x1_1b2[] += VOL_Wiremesh1_x1_1b2 + 1000 * i + j;

        sf_x1_1b2[] += s_x1_1b2[];
        ll2_x1_1b2~{j}~{i} = newll; Line Loop(ll2_x1_1b2~{j}~{i}) = Boundary{ Surface{tmp_x1_1b2[0]}; };
        llf_x1_1b2_1[] += ll2_x1_1b2~{j}~{i};

        st_tmp_x1[] += tmp_x1_1b2[0];

  k += 1;
  l += 1;

  EndFor
EndFor

Physical Surface(physsurf_x1_wire) = { st_s1_x1[], st_x1_1a1[], st_x1_1a2[], st_x1_1b1[], st_x1_1b2[], st_tmp_x1[] };
Physical Volume(physvol_x1_wire) = { vt_x1_1a1[], vt_x1_1a2[], vt_x1_1b1[], vt_x1_1b2[] };

  // y-direction

For i In {1:i_t_y}
  For j In {1:j_t_y}

  k=1;
  l=1;

  sf_y1_2a1[] = {}; // surfaces of all Wiremeshs
  llf_y1_2a1_0[] = {}; // line loops of bottom Wiremesh intersects
  llf_y1_2a1_1[] = {}; // line loops of top Wiremesh intersects

      p0_y1_2a1~{j}~{i} = newp; Point(p0_y1_2a1~{j}~{i}) = {4*(j-1)*p+geo_wc_xr,2*p+4*(i-1)*p+2*k*p+geo_wc_yd,r_w+mesh_level*mm, LcWiremesh * mm};
      p1_y1_2a1~{j}~{i} = newp; Point(p1_y1_2a1~{j}~{i}) = {4*(j-1)*p+geo_wc_xr,2*p+4*(i-1)*p+2*k*p+geo_wc_yd,2*r_w+mesh_level*mm, LcWiremesh * mm};
      p2_y1_2a1~{j}~{i} = newp; Point(p2_y1_2a1~{j}~{i}) = {4*(j-1)*p+r_w+geo_wc_xr,2*p+4*(i-1)*p+2*k*p+geo_wc_yd,r_w+mesh_level*mm, LcWiremesh * mm};
      p3_y1_2a1~{j}~{i} = newp; Point(p3_y1_2a1~{j}~{i}) = {4*(j-1)*p+geo_wc_xr,2*p+4*(i-1)*p+2*k*p+geo_wc_yd,0+mesh_level*mm, LcWiremesh * mm};
      p4_y1_2a1~{j}~{i} = newp; Point(p4_y1_2a1~{j}~{i}) = {4*(j-1)*p-r_w+geo_wc_xr,2*p+4*(i-1)*p+2*k*p+geo_wc_yd,r_w+mesh_level*mm, LcWiremesh * mm};

      l1_y1_2a1~{j}~{i} = newl; Circle(l1_y1_2a1~{j}~{i}) = {p1_y1_2a1~{j}~{i}, p0_y1_2a1~{j}~{i}, p2_y1_2a1~{j}~{i}};
      l2_y1_2a1~{j}~{i} = newl; Circle(l2_y1_2a1~{j}~{i}) = {p2_y1_2a1~{j}~{i}, p0_y1_2a1~{j}~{i}, p3_y1_2a1~{j}~{i}};
      l3_y1_2a1~{j}~{i} = newl; Circle(l3_y1_2a1~{j}~{i}) = {p3_y1_2a1~{j}~{i}, p0_y1_2a1~{j}~{i}, p4_y1_2a1~{j}~{i}};
      l4_y1_2a1~{j}~{i} = newl; Circle(l4_y1_2a1~{j}~{i}) = {p4_y1_2a1~{j}~{i}, p0_y1_2a1~{j}~{i}, p1_y1_2a1~{j}~{i}};

      ll1_y1_2a1~{j}~{i} = newll; Line Loop(ll1_y1_2a1~{j}~{i}) = {l1_y1_2a1~{j}~{i}, l2_y1_2a1~{j}~{i}, l3_y1_2a1~{j}~{i}, l4_y1_2a1~{j}~{i}};
      s1_y1_2a1~{j}~{i} = news; Plane Surface(s1_y1_2a1~{j}~{i}) = {ll1_y1_2a1~{j}~{i}};
      llf_y1_2a1_0[] += ll1_y1_2a1~{j}~{i};

      st_s1_y1[] += s1_y1_2a1~{j}~{i};

        //phys_fil_y1_2a1 = {};
        //phys_fil_top_y1_2a1 = {};
        //phys_fil_bot_y1_2a1 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_y1_2a1 (%g in layer %g)", j, i),
        //BND_Wiremesh1_y1_2a1 + 1000 * i + j) = {s1_y1_2a1~{j}~{i}}; // bottom
        //phys_fil_bot_y1_2a1 += BND_Wiremesh1_y1_2a1 + 1000 * i + j;

        v_y1_2a1[] = {};
        s_y1_2a1[] = {};
        tmp_y1_2a1[] = {s1_y1_2a1~{j}~{i}};
          tmp_y1_2a1[] = Extrude {{0,0,0},{1,0,0},{4*(j-1)*p,2*p+4*(i-1)*p+2*k*p,-R+r_w}, alpha} {
            Surface{ tmp_y1_2a1[0] };
          };

          v_y1_2a1[] += tmp_y1_2a1[1];
          vt_y1_2a1[] += v_y1_2a1[];
          s_y1_2a1[] += tmp_y1_2a1[{2:5}];
          st_y1_2a1[] += s_y1_2a1[];

        //Physical Surface(Sprintf("Wiremesh top boundary_y1_2a1 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_y1_2a1 + 1100 * i + j) = tmp_y1_2a1[0]; // top
        //phys_fil_top_y1_2a1 += BND_Wiremesh2_y1_2a1 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_y1_2a1 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_y1_2a1 + 1200 * i + j) = s_y1_2a1[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_y1_2a1 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_y1_2a1 + 1000 * i + j) = v_y1_2a1[];
        //phys_fil_y1_2a1 += VOL_Wiremesh1_y1_2a1 + 1000 * i + j;

        sf_y1_2a1[] += s_y1_2a1[];
        ll2_y1_2a1~{j}~{i} = newll; Line Loop(ll2_y1_2a1~{j}~{i}) = Boundary{ Surface{tmp_y1_2a1[0]}; };
        llf_y1_2a1_1[] += ll2_y1_2a1~{j}~{i};

  sf_y1_2a2[] = {}; // surfaces of all Wiremeshs
  llf_y1_2a2_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_y1_2a2 = {};
        //phys_fil_top_y1_2a2 = {};
        //phys_fil_bot_y1_2a2 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_y1_2a2 (%g in layer %g)", j, i),
        //BND_Wiremesh1_y1_2a2 + 1000 * i + j) = {tmp_y1_2a1[0]}; // bottom
        //phys_fil_bot_y1_2a2 += BND_Wiremesh1_y1_2a2 + 1000 * i + j;

        v_y1_2a2[] = {};
        s_y1_2a2[] = {};
        tmp_y1_2a2[] = {tmp_y1_2a1[0]};
          tmp_y1_2a2[] = Extrude {{0,0,0},{-1,0,0},{4*(j-1)*p,4*(i-1)*p+2*k*p,R-r_w}, alpha} {
            Surface{ tmp_y1_2a2[0] };
          };

          v_y1_2a2[] += tmp_y1_2a2[1];
          vt_y1_2a2[] += v_y1_2a2[];
          s_y1_2a2[] += tmp_y1_2a2[{2:5}];
          st_y1_2a2[] += s_y1_2a2[];

        //Physical Surface(Sprintf("Wiremesh top boundary_y1_2a2 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_y1_2a2 + 1100 * i + j) = tmp_y1_2a2[0]; // top
        //phys_fil_top_y1_2a2 += BND_Wiremesh2_y1_2a2 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_y1_2a2 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_y1_2a2 + 1200 * i + j) = s_y1_2a2[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_y1_2a2 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_y1_2a2 + 1000 * i + j) = v_y1_2a2[];
        //phys_fil_y1_2a2 += VOL_Wiremesh1_y1_2a2 + 1000 * i + j;

        sf_y1_2a2[] += s_y1_2a2[];
        ll2_y1_2a2~{j}~{i} = newll; Line Loop(ll2_y1_2a2~{j}~{i}) = Boundary{ Surface{tmp_y1_2a2[0]}; };
        llf_y1_2a2_1[] += ll2_y1_2a2~{j}~{i};

  sf_y1_2b1[] = {}; // surfaces of all Wiremeshs
  llf_y1_2b1_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_y1_2b1 = {};
        //phys_fil_top_y1_2b1 = {};
        //phys_fil_bot_y1_2b1 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_y1_2b1 (%g in layer %g)", j, i),
        //BND_Wiremesh1_y1_2b1 + 1000 * i + j) = {tmp_y1_2a2[0]}; // bottom
        //phys_fil_bot_y1_2b1 += BND_Wiremesh1_y1_2b1 + 1000 * i + j;

        v_y1_2b1[] = {};
        s_y1_2b1[] = {};
        tmp_y1_2b1[] = {tmp_y1_2a2[0]};
          tmp_y1_2b1[] = Extrude {{0,0,0},{-1,0,0},{4*(j-1)*p,4*(i-1)*p+2*k*p,R-r_w}, alpha} {
            Surface{ tmp_y1_2b1[0] };
          };

          v_y1_2b1[] += tmp_y1_2b1[1];
          vt_y1_2b1[] += v_y1_2b1[];
          s_y1_2b1[] += tmp_y1_2b1[{2:5}];
          st_y1_2b1[] += s_y1_2b1[];

        //Physical Surface(Sprintf("Wiremesh top boundary_y1_2b1 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_y1_2b1 + 1100 * i + j) = tmp_y1_2b1[0]; // top
        //phys_fil_top_y1_2b1 += BND_Wiremesh2_y1_2b1 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_y1_2b1 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_y1_2b1 + 1200 * i + j) = s_y1_2b1[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_y1_2b1 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_y1_2b1 + 1000 * i + j) = v_y1_2b1[];
        //phys_fil_y1_2b1 += VOL_Wiremesh1_y1_2b1 + 1000 * i + j;

        sf_y1_2b1[] += s_y1_2b1[];
        ll2_y1_2b1~{j}~{i} = newll; Line Loop(ll2_y1_2b1~{j}~{i}) = Boundary{ Surface{tmp_y1_2b1[0]}; };
        llf_y1_2b1_1[] += ll2_y1_2b1~{j}~{i};

  sf_y1_2b2[] = {}; // surfaces of all Wiremeshs
  llf_y1_2b2_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_y1_2b2 = {};
        //phys_fil_top_y1_2b2 = {};
        //phys_fil_bot_y1_2b2 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_y1_2b2 (%g in layer %g)", j, i),
        //BND_Wiremesh1_y1_2b2 + 1000 * i + j) = {tmp_y1_2b1[0]}; // bottom
        //phys_fil_bot_y1_2b2 += BND_Wiremesh1_y1_2b2 + 1000 * i + j;

        v_y1_2b2[] = {};
        s_y1_2b2[] = {};
        tmp_y1_2b2[] = {tmp_y1_2b1[0]};
          tmp_y1_2b2[] = Extrude {{0,0,0},{1,0,0},{4*(j-1)*p,-2*p+4*(i-1)*p+2*k*p,-R+r_w}, alpha} {
            Surface{ tmp_y1_2b2[0] };
          };

          v_y1_2b2[] += tmp_y1_2b2[1];
          vt_y1_2b2[] += v_y1_2b2[];
          s_y1_2b2[] += tmp_y1_2b2[{2:5}];
          st_y1_2b2[] += s_y1_2b2[];

        //Physical Surface(Sprintf("Wiremesh top boundary_y1_2b2 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_y1_2b2 + 1100 * i + j) = tmp_y1_2b2[0]; // top
        //phys_fil_top_y1_2b2 += BND_Wiremesh2_y1_2b2 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_y1_2b2 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_y1_2b2 + 1200 * i + j) = s_y1_2b2[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_y1_2b2 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_y1_2b2 + 1000 * i + j) = v_y1_2b2[];
        //phys_fil_y1_2b2 += VOL_Wiremesh1_y1_2b2 + 1000 * i + j;

        sf_y1_2b2[] += s_y1_2b2[];
        ll2_y1_2b2~{j}~{i} = newll; Line Loop(ll2_y1_2b2~{j}~{i}) = Boundary{ Surface{tmp_y1_2b2[0]}; };
        llf_y1_2b2_1[] += ll2_y1_2b2~{j}~{i};

        st_tmp_y1[] += tmp_y1_2b2[0];   

  k += 1;
  l += 1;

  EndFor
EndFor

Physical Surface(physsurf_y1_wire) = { st_s1_y1[], st_y1_2a1[], st_y1_2a2[], st_y1_2b1[], st_y1_2b2[], st_tmp_y1[] };
Physical Volume(physvol_y1_wire) = { vt_y1_2a1[], vt_y1_2a2[], vt_y1_2b1[], vt_y1_2b2[] };

// Second set of wires

  // x-direction

For i In {1:i_t_x}
  For j In {1:j_t_x-1}

  k=1;
  l=1;

  sf_x2_1b1[] = {}; // surfaces of all Wiremeshs
  llf_x2_1b1_0[] = {}; // line loops of bottom Wiremesh intersects
  llf_x2_1b1_1[] = {}; // line loops of top Wiremesh intersects

      p0_x2_1b1~{j}~{i} = newp; Point(p0_x2_1b1~{j}~{i}) = {-2*p+4*(i-1)*p+6*k*p+geo_wc_xd,4*(j-1)*p+2*l*p+geo_wc_yr,r_w+mesh_level*mm, LcWiremesh * mm};
      p1_x2_1b1~{j}~{i} = newp; Point(p1_x2_1b1~{j}~{i}) = {-2*p+4*(i-1)*p+6*k*p+geo_wc_xd,4*(j-1)*p+2*l*p+geo_wc_yr,2*r_w+mesh_level*mm, LcWiremesh * mm};
      p2_x2_1b1~{j}~{i} = newp; Point(p2_x2_1b1~{j}~{i}) = {-2*p+4*(i-1)*p+6*k*p+geo_wc_xd,r_w+4*(j-1)*p+2*l*p+geo_wc_yr,r_w+mesh_level*mm, LcWiremesh * mm};
      p3_x2_1b1~{j}~{i} = newp; Point(p3_x2_1b1~{j}~{i}) = {-2*p+4*(i-1)*p+6*k*p+geo_wc_xd,4*(j-1)*p+2*l*p+geo_wc_yr,0+mesh_level*mm, LcWiremesh * mm};
      p4_x2_1b1~{j}~{i} = newp; Point(p4_x2_1b1~{j}~{i}) = {-2*p+4*(i-1)*p+6*k*p+geo_wc_xd,-r_w+4*(j-1)*p+2*l*p+geo_wc_yr,r_w+mesh_level*mm, LcWiremesh * mm};

      l1_x2_1b1~{j}~{i} = newl; Circle(l1_x2_1b1~{j}~{i}) = {p1_x2_1b1~{j}~{i}, p0_x2_1b1~{j}~{i}, p2_x2_1b1~{j}~{i}};
      l2_x2_1b1~{j}~{i} = newl; Circle(l2_x2_1b1~{j}~{i}) = {p2_x2_1b1~{j}~{i}, p0_x2_1b1~{j}~{i}, p3_x2_1b1~{j}~{i}};
      l3_x2_1b1~{j}~{i} = newl; Circle(l3_x2_1b1~{j}~{i}) = {p3_x2_1b1~{j}~{i}, p0_x2_1b1~{j}~{i}, p4_x2_1b1~{j}~{i}};
      l4_x2_1b1~{j}~{i} = newl; Circle(l4_x2_1b1~{j}~{i}) = {p4_x2_1b1~{j}~{i}, p0_x2_1b1~{j}~{i}, p1_x2_1b1~{j}~{i}};

      ll1_x2_1b1~{j}~{i} = newll; Line Loop(ll1_x2_1b1~{j}~{i}) = {l1_x2_1b1~{j}~{i}, l2_x2_1b1~{j}~{i}, l3_x2_1b1~{j}~{i}, l4_x2_1b1~{j}~{i}};
      s1_x2_1b1~{j}~{i} = news; Plane Surface(s1_x2_1b1~{j}~{i}) = {ll1_x2_1b1~{j}~{i}};
      llf_x2_1b1_0[] += ll1_x2_1b1~{j}~{i};

      st_s1_x2[] += s1_x2_1b1~{j}~{i};

        //phys_fil_x2_1b1 = {};
        //phys_fil_top_x2_1b1 = {};
        //phys_fil_bot_x2_1b1 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_x2_1b1 (%g in layer %g)", j, i),
        //BND_Wiremesh1_x2_1b1 + 1000 * i + j) = {s1_x2_1b1~{j}~{i}}; // bottom
        //phys_fil_bot_x2_1b1 += BND_Wiremesh1_x2_1b1 + 1000 * i + j;

        v_x2_1b1[] = {};
        s_x2_1b1[] = {};
        tmp_x2_1b1[] = {s1_x2_1b1~{j}~{i}};
          tmp_x2_1b1[] = Extrude {{0,0,0},{0,-1,0},{4*(i-1)*p+4*k*p,4*(j-1)*p+2*l*p,-R+r_w}, alpha} {
            Surface{ tmp_x2_1b1[0] };
          };

          v_x2_1b1[] += tmp_x2_1b1[1];
          vt_x2_1b1[] += v_x2_1b1[];
          s_x2_1b1[] += tmp_x2_1b1[{2:5}];
          st_x2_1b1[] += s_x2_1b1[];

        //Physical Surface(Sprintf("Wiremesh top boundary_x2_1b1 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_x2_1b1 + 1100 * i + j) = tmp_x2_1b1[0]; // top
        //phys_fil_top_x2_1b1 += BND_Wiremesh2_x2_1b1 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_x2_1b1 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_x2_1b1 + 1200 * i + j) = s_x2_1b1[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_x2_1b1 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_x2_1b1 + 1000 * i + j) = v_x2_1b1[];
        //phys_fil_x2_1b1 += VOL_Wiremesh1_x2_1b1 + 1000 * i + j;

        sf_x2_1b1[] += s_x2_1b1[];
        ll2_x2_1b1~{j}~{i} = newll; Line Loop(ll2_x2_1b1~{j}~{i}) = Boundary{ Surface{tmp_x2_1b1[0]}; };
        llf_x2_1b1_1[] += ll2_x2_1b1~{j}~{i};

  sf_x2_1b2[] = {}; // surfaces of all Wiremeshs
  llf_x2_1b2_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_x2_1b2 = {};
        //phys_fil_top_x2_1b2 = {};
        //phys_fil_bot_x2_1b2 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_x2_1b2 (%g in layer %g)", j, i),
        //BND_Wiremesh1_x2_1b2 + 1000 * i + j) = {tmp_x2_1b1[0]}; // bottom
        //phys_fil_bot_x2_1b2 += BND_Wiremesh1_x2_1b2 + 1000 * i + j;

        v_x2_1b2[] = {};
        s_x2_1b2[] = {};
        tmp_x2_1b2[] = {tmp_x2_1b1[0]};
          tmp_x2_1b2[] = Extrude {{0,0,0},{0,1,0},{4*(i-1)*p+2*k*p,4*(j-1)*p+l*p,R-r_w}, alpha} {
            Surface{ tmp_x2_1b2[0] };
          };

          v_x2_1b2[] += tmp_x2_1b2[1];
		  vt_x2_1b2[] += v_x2_1b2[];
          s_x2_1b2[] += tmp_x2_1b2[{2:5}];
          st_x2_1b2[] += s_x2_1b2[];

        //Physical Surface(Sprintf("Wiremesh top boundary_x2_1b2 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_x2_1b2 + 1100 * i + j) = tmp_x2_1b2[0]; // top
        //phys_fil_top_x2_1b2 += BND_Wiremesh2_x2_1b2 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_x2_1b2 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_x2_1b2 + 1200 * i + j) = s_x2_1b2[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_x2_1b2 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_x2_1b2 + 1000 * i + j) = v_x2_1b2[];
        //phys_fil_x2_1b2 += VOL_Wiremesh1_x2_1b2 + 1000 * i + j;

        sf_x2_1b2[] += s_x2_1b2[];
        ll2_x2_1b2~{j}~{i} = newll; Line Loop(ll2_x2_1b2~{j}~{i}) = Boundary{ Surface{tmp_x2_1b2[0]}; };
        llf_x2_1b2_1[] += ll2_x2_1b2~{j}~{i};

  sf_x2_1a1[] = {}; // surfaces of all Wiremeshs
  llf_x2_1a1_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_x2_1a1 = {};
        //phys_fil_top_x2_1a1 = {};
        //phys_fil_bot_x2_1a1 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_x2_1a1 (%g in layer %g)", j, i),
        //BND_Wiremesh1_x2_1a1 + 1000 * i + j) = {tmp_x2_1b2[0]}; // bottom
        //phys_fil_bot_x2_1a1 += BND_Wiremesh1_x2_1a1 + 1000 * i + j;

        v_x2_1a1[] = {};
        s_x2_1a1[] = {};
        tmp_x2_1a1[] = {tmp_x2_1b2[0]};
          tmp_x2_1a1[] = Extrude {{0,0,0}, {0,1,0}, {p+4*(i-1)*p+k*p,p+4*(j-1)*p+l*p,R-r_w}, alpha} {
            Surface{ tmp_x2_1a1[0] };
          };

          v_x2_1a1[] += tmp_x2_1a1[1];
          vt_x2_1a1[] += v_x2_1a1[];
          s_x2_1a1[] += tmp_x2_1a1[{2:5}];
          st_x2_1a1[] += s_x2_1a1[];

        //Physical Surface(Sprintf("Wiremesh top boundary_x2_1a1 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_x2_1a1 + 1100 * i + j) = tmp_x2_1a1[0]; // top
        //phys_fil_top_x2_1a1 += BND_Wiremesh2_x2_1a1 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_x2_1a1 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_x2_1a1 + 1200 * i + j) = s_x2_1a1[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_x2_1a1 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_x2_1a1 + 1000 * i + j) = v_x2_1a1[];
        //phys_fil_x2_1a1 += VOL_Wiremesh1_x2_1a1 + 1000 * i + j;

        sf_x2_1a1[] += s_x2_1a1[];
        ll2_x2_1a1~{j}~{i} = newll; Line Loop(ll2_x2_1a1~{j}~{i}) = Boundary{ Surface{tmp_x2_1a1[0]}; };
        llf_x2_1a1_1[] += ll2_x2_1a1~{j}~{i};

  sf_x2_1a2[] = {}; // surfaces of all Wiremeshs
  llf_x2_1a2_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_x2_1a2 = {};
        //phys_fil_top_x2_1a2 = {};
        //phys_fil_bot_x2_1a2 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_x2_1a2 (%g in layer %g)", j, i),
        //BND_Wiremesh1_x2_1a2 + 1000 * i + j) = {tmp_x2_1a1[0]}; // bottom
        //phys_fil_bot_x2_1a2 += BND_Wiremesh1_x2_1a2 + 1000 * i + j;

        v_x2_1a2[] = {};
        s_x2_1a2[] = {};
        tmp_x2_1a2[] = {tmp_x2_1a1[0]};
          tmp_x2_1a2[] = Extrude {{0,0,0},{0,-1,0},{-p+4*(i-1)*p+k*p,-p+4*(j-1)*p+2*l*p,-R+r_w}, alpha} {
            Surface{ tmp_x2_1a2[0] };
          };

          v_x2_1a2[] += tmp_x2_1a2[1];
          vt_x2_1a2[] += v_x2_1a2[];
          s_x2_1a2[] += tmp_x2_1a2[{2:5}];
		  st_x2_1a2[] += s_x2_1a2[];

        //Physical Surface(Sprintf("Wiremesh top boundary_x2_1a2 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_x2_1a2 + 1100 * i + j) = tmp_x2_1a2[0]; // top
        //phys_fil_top_x2_1a2 += BND_Wiremesh2_x2_1a2 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_x2_1a2 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_x2_1a2 + 1200 * i + j) = s_x2_1a2[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_x2_1a2 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_x2_1a2 + 1000 * i + j) = v_x2_1a2[];
        //phys_fil_x2_1a2 += VOL_Wiremesh1_x2_1a2 + 1000 * i + j;

        sf_x2_1a2[] += s_x2_1a2[];
        ll2_x2_1a2~{j}~{i} = newll; Line Loop(ll2_x2_1a2~{j}~{i}) = Boundary{ Surface{tmp_x2_1a2[0]}; };
        llf_x2_1a2_1[] += ll2_x2_1a2~{j}~{i};

        st_tmp_x2[] += tmp_x2_1a2[0];

  k += 1;
  l += 1;

  EndFor
EndFor

Physical Surface(physsurf_x2_wire) = { st_s1_x2[], st_x2_1b1[], st_x2_1b2[], st_x2_1a1[], st_x2_1a2[], st_tmp_x2[] };
Physical Volume(physvol_x2_wire) = { vt_x2_1b1[], vt_x2_1b2[], vt_x2_1a1[], vt_x2_1a2[] };

  // y-direction

For i In {1:i_t_y}
  For j In {1:j_t_y-1}

  k=1;
  l=1;

  sf_y2_2b1[] = {}; // surfaces of all Wiremeshs
  llf_y2_2b1_0[] = {}; // line loops of bottom Wiremesh intersects
  llf_y2_2b1_1[] = {}; // line loops of top Wiremesh intersects

      p0_y2_2b1~{j}~{i} = newp; Point(p0_y2_2b1~{j}~{i}) = {4*(j-1)*p+2*l*p+geo_wc_xr,2*p+4*(i-1)*p+2*k*p+geo_wc_yd,-r_w+mesh_level*mm, LcWiremesh * mm};
      p1_y2_2b1~{j}~{i} = newp; Point(p1_y2_2b1~{j}~{i}) = {4*(j-1)*p+2*l*p+geo_wc_xr,2*p+4*(i-1)*p+2*k*p+geo_wc_yd,-2*r_w+mesh_level*mm, LcWiremesh * mm};
      p2_y2_2b1~{j}~{i} = newp; Point(p2_y2_2b1~{j}~{i}) = {4*(j-1)*p+2*l*p+r_w+geo_wc_xr,2*p+4*(i-1)*p+2*k*p+geo_wc_yd,-r_w+mesh_level*mm, LcWiremesh * mm};
      p3_y2_2b1~{j}~{i} = newp; Point(p3_y2_2b1~{j}~{i}) = {4*(j-1)*p+2*l*p+geo_wc_xr,2*p+4*(i-1)*p+2*k*p+geo_wc_yd,0+mesh_level*mm, LcWiremesh * mm};
      p4_y2_2b1~{j}~{i} = newp; Point(p4_y2_2b1~{j}~{i}) = {4*(j-1)*p+2*l*p-r_w+geo_wc_xr,2*p+4*(i-1)*p+2*k*p+geo_wc_yd,-r_w+mesh_level*mm, LcWiremesh * mm};

      l1_y2_2b1~{j}~{i} = newl; Circle(l1_y2_2b1~{j}~{i}) = {p1_y2_2b1~{j}~{i}, p0_y2_2b1~{j}~{i}, p2_y2_2b1~{j}~{i}};
      l2_y2_2b1~{j}~{i} = newl; Circle(l2_y2_2b1~{j}~{i}) = {p2_y2_2b1~{j}~{i}, p0_y2_2b1~{j}~{i}, p3_y2_2b1~{j}~{i}};
      l3_y2_2b1~{j}~{i} = newl; Circle(l3_y2_2b1~{j}~{i}) = {p3_y2_2b1~{j}~{i}, p0_y2_2b1~{j}~{i}, p4_y2_2b1~{j}~{i}};
      l4_y2_2b1~{j}~{i} = newl; Circle(l4_y2_2b1~{j}~{i}) = {p4_y2_2b1~{j}~{i}, p0_y2_2b1~{j}~{i}, p1_y2_2b1~{j}~{i}};

      ll1_y2_2b1~{j}~{i} = newll; Line Loop(ll1_y2_2b1~{j}~{i}) = {l1_y2_2b1~{j}~{i}, l2_y2_2b1~{j}~{i}, l3_y2_2b1~{j}~{i}, l4_y2_2b1~{j}~{i}};
      s1_y2_2b1~{j}~{i} = news; Plane Surface(s1_y2_2b1~{j}~{i}) = {ll1_y2_2b1~{j}~{i}};
      llf_y2_2b1_0[] += ll1_y2_2b1~{j}~{i};

      st_s1_y2[] += s1_y2_2b1~{j}~{i};

        //phys_fil_y2_2b1 = {};
        //phys_fil_top_y2_2b1 = {};
        //phys_fil_bot_y2_2b1 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_y2_2b1 (%g in layer %g)", j, i),
        //BND_Wiremesh1_y2_2b1 + 1000 * i + j) = {s1_y2_2b1~{j}~{i}}; // bottom
        //phys_fil_bot_y2_2b1 += BND_Wiremesh1_y2_2b1 + 1000 * i + j;

        v_y2_2b1[] = {};
        s_y2_2b1[] = {};
        tmp_y2_2b1[] = {s1_y2_2b1~{j}~{i}};
          tmp_y2_2b1[] = Extrude {{0,0,0},{-1,0,0},{4*(j-1)*p+2*l*p,4*(i-1)*p+4*k*p,R-r_w}, alpha} {
            Surface{ tmp_y2_2b1[0] };
          };

          v_y2_2b1[] += tmp_y2_2b1[1];
          vt_y2_2b1[] += v_y2_2b1[]; 
          s_y2_2b1[] += tmp_y2_2b1[{2:5}];
          st_y2_2b1[] += s_y2_2b1[];

        //Physical Surface(Sprintf("Wiremesh top boundary_y2_2b1 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_y2_2b1 + 1100 * i + j) = tmp_y2_2b1[0]; // top
        //phys_fil_top_y2_2b1 += BND_Wiremesh2_y2_2b1 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_y2_2b1 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_y2_2b1 + 1200 * i + j) = s_y2_2b1[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_y2_2b1 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_y2_2b1 + 1000 * i + j) = v_y2_2b1[];
        //phys_fil_y2_2b1 += VOL_Wiremesh1_y2_2b1 + 1000 * i + j;

        sf_y2_2b1[] += s_y2_2b1[];
        ll2_y2_2b1~{j}~{i} = newll; Line Loop(ll2_y2_2b1~{j}~{i}) = Boundary{ Surface{tmp_y2_2b1[0]}; };
        llf_y2_2b1_1[] += ll2_y2_2b1~{j}~{i};

  sf_y2_2b2[] = {}; // surfaces of all Wiremeshs
  llf_y2_2b2_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_y2_2b2 = {};
        //phys_fil_top_y2_2b2 = {};
        //phys_fil_bot_y2_2b2 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_y2_2b2 (%g in layer %g)", j, i),
        //BND_Wiremesh1_y2_2b2 + 1000 * i + j) = {tmp_y2_2b1[0]}; // bottom
        //phys_fil_bot_y2_2b2 += BND_Wiremesh1_y2_2b2 + 1000 * i + j;

        v_y2_2b2[] = {};
        s_y2_2b2[] = {};
        tmp_y2_2b2[] = {tmp_y2_2b1[0]};
          tmp_y2_2b2[] = Extrude {{0,0,0},{1,0,0},{4*(j-1)*p+2*l*p,2*p-4*p+4*(i-1)*p+4*k*p,-R+r_w}, alpha} {
            Surface{ tmp_y2_2b2[0] };
          };

          v_y2_2b2[] += tmp_y2_2b2[1];
          vt_y2_2b2[] += v_y2_2b2[];
          s_y2_2b2[] += tmp_y2_2b2[{2:5}];
          st_y2_2b2[] += s_y2_2b2[];

        //Physical Surface(Sprintf("Wiremesh top boundary_y2_2b2 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_y2_2b2 + 1100 * i + j) = tmp_y2_2b2[0]; // top
        //phys_fil_top_y2_2b2 += BND_Wiremesh2_y2_2b2 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_y2_2b2 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_y2_2b2 + 1200 * i + j) = s_y2_2b2[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_y2_2b2 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_y2_2b2 + 1000 * i + j) = v_y2_2b2[];
        //phys_fil_y2_2b2 += VOL_Wiremesh1_y2_2b2 + 1000 * i + j;

        sf_y2_2b2[] += s_y2_2b2[];
        ll2_y2_2b2~{j}~{i} = newll; Line Loop(ll2_y2_2b2~{j}~{i}) = Boundary{ Surface{tmp_y2_2b2[0]}; };
        llf_y2_2b2_1[] += ll2_y2_2b2~{j}~{i};

  sf_y2_2a1[] = {}; // surfaces of all Wiremeshs
  llf_y2_2a1_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_y2_2a1 = {};
        //phys_fil_top_y2_2a1 = {};
        //phys_fil_bot_y2_2a1 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_y2_2a1 (%g in layer %g)", j, i),
        //BND_Wiremesh1_y2_2a1 + 1000 * i + j) = {tmp_y2_2b2[0]}; // bottom
        //phys_fil_bot_y2_2a1 += BND_Wiremesh1_y2_2a1 + 1000 * i + j;

        v_y2_2a1[] = {};
        s_y2_2a1[] = {};
        tmp_y2_2a1[] = {tmp_y2_2b2[0]};
          tmp_y2_2a1[] = Extrude {{0,0,0},{1,0,0},{4*(j-1)*p+2*l*p,p+4*(i-1)*p+k*p,-R+r_w}, alpha} {
            Surface{ tmp_y2_2a1[0] };
          };

          v_y2_2a1[] += tmp_y2_2a1[1];
          vt_y2_2a1[] += v_y2_2a1[];
          s_y2_2a1[] += tmp_y2_2a1[{2:5}];
          st_y2_2a1[] += s_y2_2a1[];

        //Physical Surface(Sprintf("Wiremesh top boundary_y2_2a1 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_y2_2a1 + 1100 * i + j) = tmp_y2_2a1[0]; // top
        //phys_fil_top_y2_2a1 += BND_Wiremesh2_y2_2a1 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_y2_2a1 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_y2_2a1 + 1200 * i + j) = s_y2_2a1[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_y2_2a1 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_y2_2a1 + 1000 * i + j) = v_y2_2a1[];
        //phys_fil_y2_2a1 += VOL_Wiremesh1_y2_2a1 + 1000 * i + j;

        sf_y2_2a1[] += s_y2_2a1[];
        ll2_y2_2a1~{j}~{i} = newll; Line Loop(ll2_y2_2a1~{j}~{i}) = Boundary{ Surface{tmp_y2_2a1[0]}; };
        llf_y2_2a1_1[] += ll2_y2_2a1~{j}~{i};

  sf_y2_2a2[] = {}; // surfaces of all Wiremeshs
  llf_y2_2a2_1[] = {}; // line loops of top Wiremesh intersects

        //phys_fil_y2_2a2 = {};
        //phys_fil_top_y2_2a2 = {};
        //phys_fil_bot_y2_2a2 = {};

        //Physical Surface(Sprintf("Wiremesh bottom boundary_y2_2a2 (%g in layer %g)", j, i),
        //BND_Wiremesh1_y2_2a2 + 1000 * i + j) = {tmp_y2_2a1[0]}; // bottom
        //phys_fil_bot_y2_2a2 += BND_Wiremesh1_y2_2a2 + 1000 * i + j;

        v_y2_2a2[] = {};
        s_y2_2a2[] = {};
        tmp_y2_2a2[] = {tmp_y2_2a1[0]};
          tmp_y2_2a2[] = Extrude {{0,0,0},{-1,0,0},{4*(j-1)*p+2*l*p,-p+4*(i-1)*p+k*p,R-r_w}, alpha} {
            Surface{ tmp_y2_2a2[0] };
          };

          v_y2_2a2[] += tmp_y2_2a2[1];
          vt_y2_2a2[] += v_y2_2a2[];
          s_y2_2a2[] += tmp_y2_2a2[{2:5}];
          st_y2_2a2[] += s_y2_2a2[];

        //Physical Surface(Sprintf("Wiremesh top boundary_y2_2a2 (%g in layer %g)", j, i),
        //  BND_Wiremesh2_y2_2a2 + 1100 * i + j) = tmp_y2_2a2[0]; // top
        //phys_fil_top_y2_2a2 += BND_Wiremesh2_y2_2a2 + 1100 * i + j;
        //Physical Surface(Sprintf("Wiremesh lateral boundary_y2_2a2 (%g in layer %g)", j, i),
        //  BND_Wiremesh3_y2_2a2 + 1200 * i + j) = s_y2_2a2[]; // sides
        //Physical Volume(Sprintf("Wiremesh volume_y2_2a2 (%g in layer %g)", j, i),
        //  VOL_Wiremesh1_y2_2a2 + 1000 * i + j) = v_y2_2a2[];
        //phys_fil_y2_2a2 += VOL_Wiremesh1_y2_2a2 + 1000 * i + j;

        sf_y2_2a2[] += s_y2_2a2[];
        ll2_y2_2a2~{j}~{i} = newll; Line Loop(ll2_y2_2a2~{j}~{i}) = Boundary{ Surface{tmp_y2_2a2[0]}; };
        llf_y2_2a2_1[] += ll2_y2_2a2~{j}~{i};

        st_tmp_y2[] += tmp_y2_2a2[0];

  EndFor
EndFor

Physical Surface(physsurf_y2_wire) = { st_s1_y2[], st_y2_2b1[], st_y2_2b2[], st_y2_2a1[], st_y2_2a2[], st_tmp_y2[] };
Physical Volume(physvol_y2_wire) = { vt_y2_2b1[], vt_y2_2b2[], vt_y2_2a1[], vt_y2_2a2[] };


// *********************************************************************

// SHELL

// *******************************
// Corner 1
// *******************************
pc1_1 = newp; Point(pc1_1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -tuC/2,lcCopperPlateBdry};
pc2_1 = newp; Point(pc2_1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*tD/2,lcCopperPlateBdry};
pc3_1 = newp; Point(pc3_1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -(2*tC+tuC)/2,lcCopperPlateBdry};
pc4_1 = newp; Point(pc4_1) = {0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};

// *******************************
// Corner 2
// *******************************
pc1_2 = newp; Point(pc1_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -tuC/2,lcCopperPlateBdry};
pc2_2 = newp; Point(pc2_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*tD/2,lcCopperPlateBdry};
pc3_2 = newp; Point(pc3_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -(2*tC+tuC)/2,lcCopperPlateBdry};
pc4_2 = newp; Point(pc4_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};

// *******************************
// Corner 3
// *******************************
pc1_3 = newp; Point(pc1_3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -tuC/2,lcCopperPlateBdry};
pc2_3 = newp; Point(pc2_3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*tD/2,lcCopperPlateBdry};
pc3_3 = newp; Point(pc3_3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -(2*tC+tuC)/2,lcCopperPlateBdry};
pc4_3 = newp; Point(pc4_3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};

// *******************************
// Corner 4
// *******************************
pc1_4 = newp; Point(pc1_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -tuC/2,lcCopperPlateBdry};
pc2_4 = newp; Point(pc2_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*tD/2,lcCopperPlateBdry};
pc3_4 = newp; Point(pc3_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -(2*tC+tuC)/2,lcCopperPlateBdry};
pc4_4 = newp; Point(pc4_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};

// *******************************************************
// Copper planes
// *******************************************************

// Points between two half pillars on upper LEM
ptmc_1 = newp; Point(ptmc_1) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -(2*tC+tuC)/2, lcCopperPlateBdry};
ptmd_1 = newp; Point(ptmd_1) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -tuC/2, lcCopperPlateBdry};

ptmc_2 = newp; Point(ptmc_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, -(2*tC+tuC)/2, lcCopperPlateBdry};
ptmd_2 = newp; Point(ptmd_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, -tuC/2, lcCopperPlateBdry};

ptmc_3 = newp; Point(ptmc_3) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -(2*tC+tuC)/2, lcCopperPlateBdry};
ptmd_3 = newp; Point(ptmd_3) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -tuC/2, lcCopperPlateBdry};

ptmc_4 = newp; Point(ptmc_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, -(2*tC+tuC)/2, lcCopperPlateBdry};
ptmd_4 = newp; Point(ptmd_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, -tuC/2, lcCopperPlateBdry};

// Top lower boundary
pcptl1 = newp; Point(pcptl1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tD/2,lcCopperPlateBdry};
pcptl2 = newp; Point(pcptl2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, tD/2,lcCopperPlateBdry};
pcptl3 = newp; Point(pcptl3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tD/2,lcCopperPlateBdry};
pcptl4 = newp; Point(pcptl4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, tD/2,lcCopperPlateBdry};

// Top upper boundary
pcptu1 = newp; Point(pcptu1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, (2*tC+tD)/2,lcCopperPlateBdry};
pcptu2 = newp; Point(pcptu2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, (2*tC+tD)/2,lcCopperPlateBdry};
pcptu3 = newp; Point(pcptu3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, (2*tC+tD)/2,lcCopperPlateBdry};
pcptu4 = newp; Point(pcptu4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, (2*tC+tD)/2,lcCopperPlateBdry};

// Border lines
// Upper boundary
lcptub1a = newc; Line(lcptub1a) = {pc3_1,ptmc_1};
lcptub1b = newc; Line(lcptub1b) = {ptmc_1,pc3_2};
lcptub2a = newc; Line(lcptub2a) = {pc3_2,ptmc_2};
lcptub2b = newc; Line(lcptub2b) = {ptmc_2,pc3_3};
lcptub3a = newc; Line(lcptub3a) = {pc3_3,ptmc_3};
lcptub3b = newc; Line(lcptub3b) = {ptmc_3,pc3_4};
lcptub4a = newc; Line(lcptub4a) = {pc3_4,ptmc_4};
lcptub4b = newc; Line(lcptub4b) = {ptmc_4,pc3_1};

// Lower boundary
lcptlb5a = newc; Line(lcptlb5a) = {pc1_1,ptmd_1};
lcptlb5b = newc; Line(lcptlb5b) = {ptmd_1,pc1_2};
lcptlb6a = newc; Line(lcptlb6a) = {pc1_2,ptmd_2};
lcptlb6b = newc; Line(lcptlb6b) = {ptmd_2,pc1_3};
lcptlb7a = newc; Line(lcptlb7a) = {pc1_3,ptmd_3};
lcptlb7b = newc; Line(lcptlb7b) = {ptmd_3,pc1_4};
lcptlb8a = newc; Line(lcptlb8a) = {pc1_4,ptmd_4};
lcptlb8b = newc; Line(lcptlb8b) = {ptmd_4,pc1_1};

// Connect the upper and lower points with lines to form the plate
lcptib9 = newc; Line(lcptib9) = {pc3_1, pc1_1};
lcptib10 = newc; Line(lcptib10) = {pc3_2, pc1_2};
lcptib11 = newc; Line(lcptib11) = {pc3_3, pc1_3};
lcptib12 = newc; Line(lcptib12) = {pc3_4, pc1_4};

// ---------------------------------------------

// Points between two half pillars on lower LEM
pbmd_1 = newp; Point(pbmd_1) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*tD/2, lcCopperPlateBdry};
pbmc_1 = newp; Point(pbmc_1) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*(2*tC+tD)/2, lcCopperPlateBdry};

pbmd_2 = newp; Point(pbmd_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, -1*tD/2, lcCopperPlateBdry};
pbmc_2 = newp; Point(pbmc_2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, -1*(2*tC+tD)/2, lcCopperPlateBdry};

pbmd_3 = newp; Point(pbmd_3) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*tD/2, lcCopperPlateBdry};
pbmc_3 = newp; Point(pbmc_3) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*(2*tC+tD)/2, lcCopperPlateBdry};

pbmd_4 = newp; Point(pbmd_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, -1*tD/2, lcCopperPlateBdry};
pbmc_4 = newp; Point(pbmc_4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a/2+geo_f_y*n_1*a, -1*(2*tC+tD)/2, lcCopperPlateBdry};

// Bottom lower boundary
pcpbl1 = newp; Point(pcpbl1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};
pcpbl2 = newp; Point(pcpbl2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};
pcpbl3 = newp; Point(pcpbl3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};
pcpbl4 = newp; Point(pcpbl4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*(2*tC+tD)/2,lcCopperPlateBdry};

// Bottom upper boundary
pcpbu1 = newp; Point(pcpbu1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*tD/2,lcCopperPlateBdry};
pcpbu2 = newp; Point(pcpbu2) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*tD/2,lcCopperPlateBdry};
pcpbu3 = newp; Point(pcpbu3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*tD/2,lcCopperPlateBdry};
pcpbu4 = newp; Point(pcpbu4) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*tD/2,lcCopperPlateBdry};

// Border lines
// Upper boundary
lcpbub1a = newc; Line(lcpbub1a) = {pc4_1,pbmc_1};
lcpbub1b = newc; Line(lcpbub1b) = {pbmc_1,pc4_2};
lcpbub2a = newc; Line(lcpbub2a) = {pc4_2,pbmc_2};
lcpbub2b = newc; Line(lcpbub2b) = {pbmc_2,pc4_3};
lcpbub3a = newc; Line(lcpbub3a) = {pc4_3,pbmc_3};
lcpbub3b = newc; Line(lcpbub3b) = {pbmc_3,pc4_4};
lcpbub4a = newc; Line(lcpbub4a) = {pc4_4,pbmc_4};
lcpbub4b = newc; Line(lcpbub4b) = {pbmc_4,pc4_1};

// Lower boundary
lcpblb5a = newc; Line(lcpblb5a) = {pc2_1,pbmd_1};
lcpblb5b = newc; Line(lcpblb5b) = {pbmd_1,pc2_2};
lcpblb6a = newc; Line(lcpblb6a) = {pc2_2,pbmd_2};
lcpblb6b = newc; Line(lcpblb6b) = {pbmd_2,pc2_3};
lcpblb7a = newc; Line(lcpblb7a) = {pc2_3,pbmd_3};
lcpblb7b = newc; Line(lcpblb7b) = {pbmd_3,pc2_4};
lcpblb8a = newc; Line(lcpblb8a) = {pc2_4,pbmd_4};
lcpblb8b = newc; Line(lcpblb8b) = {pbmd_4,pc2_1};

// Connect the upper and lower points with lines to form the plate
lcpbib9 = newc; Line(lcpbib9) = {pc4_1, pc2_1};
lcpbib10 = newc; Line(lcpbib10) = {pc4_2, pc2_2};
lcpbib11 = newc; Line(lcpbib11) = {pc4_3, pc2_3};
lcpbib12 = newc; Line(lcpbib12) = {pc4_4, pc2_4};

// Lines connecting the upper and lower level corners
lcorner1 = newc; Line(lcorner1) = {pc1_1, pc2_1};
lcorner2 = newc; Line(lcorner2) = {pc1_2, pc2_2};
lcorner3 = newc; Line(lcorner3) = {pc1_3, pc2_3};
lcorner4 = newc; Line(lcorner4) = {pc1_4, pc2_4};

// Lines splitting the LEM in half
lmid1_1 = newc; Line(lmid1_1) = {ptmc_1, ptmd_1};
lmid1_2 = newc; Line(lmid1_2) = {ptmd_1, pbmd_1};
lmid1_3 = newc; Line(lmid1_3) = {pbmd_1, pbmc_1};

lmid2_1 = newc; Line(lmid2_1) = {ptmc_2, ptmd_2};
lmid2_2 = newc; Line(lmid2_2) = {ptmd_2, pbmd_2};
lmid2_3 = newc; Line(lmid2_3) = {pbmd_2, pbmc_2};

lmid3_1 = newc; Line(lmid3_1) = {ptmc_3, ptmd_3};
lmid3_2 = newc; Line(lmid3_2) = {ptmd_3, pbmd_3};
lmid3_3 = newc; Line(lmid3_3) = {pbmd_3, pbmc_3};

lmid4_1 = newc; Line(lmid4_1) = {ptmc_4, ptmd_4};
lmid4_2 = newc; Line(lmid4_2) = {ptmd_4, pbmd_4};
lmid4_3 = newc; Line(lmid4_3) = {pbmd_4, pbmc_4};

// **********************************************
// External Electrodes
// **********************************************

// Top electrode
pexet1 = newp; Point(pexet1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};
pexet2 = newp; Point(pexet2) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};
pexet3 = newp; Point(pexet3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};
pexet4 = newp; Point(pexet4) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};
pexet5 = newp; Point(pexet5) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};
pexet6 = newp; Point(pexet6) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, (2*tC+tD)/2+lE,lcExtElectrodeBdry};

// Top electrode lines
lexet1 = newc; Line(lexet1) = {pexet1, pexet2};
lexet2 = newc; Line(lexet2) = {pexet2, pexet3};
lexet3 = newc; Line(lexet3) = {pexet3, pexet4};
lexet4 = newc; Line(lexet4) = {pexet4, pexet5};
lexet5 = newc; Line(lexet5) = {pexet5, pexet6};
lexet6 = newc; Line(lexet6) = {pexet6, pexet1};

// Connect the top electrode to the LEM.
lexetc1 = newc; Line(lexetc1) = {pexet1, pc3_1};
lexetc2 = newc; Line(lexetc2) = {pexet2, ptmc_1};
lexetc3 = newc; Line(lexetc3) = {pexet3, pc3_2};
lexetc4 = newc; Line(lexetc4) = {pexet4, pc3_3};
lexetc5 = newc; Line(lexetc5) = {pexet5, ptmc_3};
lexetc6 = newc; Line(lexetc6) = {pexet6, pc3_4};

// Bottom electrode
// pexeb1 = newp; Point(pexeb1) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*(2*tC+tD)/2-lP,lcExtElectrodeBdry};
// pexeb2 = newp; Point(pexeb2) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*(2*tC+tD)/2-lP,lcExtElectrodeBdry};
// pexeb3 = newp; Point(pexeb3) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*0+geo_f_y*n_1*a, -1*(2*tC+tD)/2-lP,lcExtElectrodeBdry};
// pexeb4 = newp; Point(pexeb4) = {geo_f_x*a+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*(2*tC+tD)/2-lP,lcExtElectrodeBdry};
// pexeb5 = newp; Point(pexeb5) = {geo_f_x*a/2+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*(2*tC+tD)/2-lP,lcExtElectrodeBdry};
// pexeb6 = newp; Point(pexeb6) = {geo_f_x*0+geo_f_x*m_1*a, geo_f_y*a+geo_f_y*n_1*a, -1*(2*tC+tD)/2-lP,lcExtElectrodeBdry};

// Copper plate surfaces
llcp_up_border1 = newreg; Line Loop(llcp_up_border1) = {lcptib9, lcptlb5a, lcptlb5b, -lcptib10, -lcptub1a, -lcptub1b};
pscp_up_border1 = newreg; Plane Surface(pscp_up_border1) = {llcp_up_border1};
llcp_up_border2 = newreg; Line Loop(llcp_up_border2) = {lcptib10, lcptlb6a, lcptlb6b, -lcptib11, -lcptub2a, -lcptub2b};
pscp_up_border2 = newreg; Plane Surface(pscp_up_border2) = {llcp_up_border2};
llcp_up_border3 = newreg; Line Loop(llcp_up_border3) = {lcptib11, lcptlb7a, lcptlb7b, -lcptib12, -lcptub3a, -lcptub3b};
pscp_up_border3 = newreg; Plane Surface(pscp_up_border3) = {llcp_up_border3};
llcp_up_border4 = newreg; Line Loop(llcp_up_border4) = {lcptib12, lcptlb8a, lcptlb8b, -lcptib9, -lcptub4a, -lcptub4b};
pscp_up_border4 = newreg; Plane Surface(pscp_up_border4) = {llcp_up_border4};

llcp_low_border1 = newreg; Line Loop(llcp_low_border1) = {lcpbib9, lcpblb5a, lcpblb5b, -lcpbib10, -lcpbub1a, -lcpbub1b};
pscp_low_border1 = newreg; Plane Surface(pscp_low_border1) = {llcp_low_border1}; 
llcp_low_border2 = newreg; Line Loop(llcp_low_border2) = {lcpbib10, lcpblb6a, lcpblb6b, -lcpbib11, -lcpbub2a, -lcpbub2b};
pscp_low_border2 = newreg; Plane Surface(pscp_low_border2) = {llcp_low_border2};
llcp_low_border3 = newreg; Line Loop(llcp_low_border3) = {lcpbib11, lcpblb7a, lcpblb7b, -lcpbib12, -lcpbub3a, -lcpbub3b};
pscp_low_border3 = newreg; Plane Surface(pscp_low_border3) = {llcp_low_border3};
llcp_low_border4 = newreg; Line Loop(llcp_low_border4) = {lcpbib12, lcpblb8a, lcpblb8b, -lcpbib9, -lcpbub4a, -lcpbub4b};
pscp_low_border4 = newreg; Plane Surface(pscp_low_border4) = {llcp_low_border4};

//llcp_face1 = newreg; Line Loop(llcp_face1) = {lcptub1a, lcptub1b, lcptub2a, lcptub2b, lcptub3a, lcptub3b, lcptub4a, lcptub4b};
//llcp_face3 = newreg; Line Loop(llcp_face3) = {lcpbub1a, lcpbub1b, lcpbub2a, lcpbub2b, lcpbub3a, lcpbub3b, lcpbub4a, lcpbub4b};

// Bounding surfaces
ll_bsurf1 = newreg; Line Loop(ll_bsurf1) = {lexet1, lexetc2, -lcptub1a, -lexetc1};
ps_bsurf1 = newreg; Plane Surface(ps_bsurf1) = {ll_bsurf1};

ll_bsurf2 = newreg; Line Loop(ll_bsurf2) = {lexet2, lexetc3, -lcptub1b, -lexetc2};
ps_bsurf2 = newreg; Plane Surface(ps_bsurf2) = {ll_bsurf2};

ll_bsurf3 = newreg; Line Loop(ll_bsurf3) = {-lexet3, lexetc3, lcptub2a, lcptub2b, -lexetc4};
ps_bsurf3 = newreg; Plane Surface(ps_bsurf3) = {ll_bsurf3};

ll_bsurf4 = newreg; Line Loop(ll_bsurf4) = {lexet4, lexetc5, -lcptub3a, -lexetc4};
ps_bsurf4 = newreg; Plane Surface(ps_bsurf4) = {ll_bsurf4};

ll_bsurf5 = newreg; Line Loop(ll_bsurf5) = {lexet5, lexetc6, -lcptub3b, -lexetc5};
ps_bsurf5 = newreg; Plane Surface(ps_bsurf5) = {ll_bsurf5};

ll_bsurf6 = newreg; Line Loop(ll_bsurf6) = {-lexet6, lexetc6, lcptub4a, lcptub4b, -lexetc1};
ps_bsurf6 = newreg; Plane Surface(ps_bsurf6) = {ll_bsurf6};

ll_bsurf7 = newreg; Line Loop(ll_bsurf7) = {lexet1, lexet2, lexet3, lexet4, lexet5, lexet6};
ps_bsurf7 = newreg; Plane Surface(ps_bsurf7) = {ll_bsurf7};

For m In {0:Number_Units_x}
 For n In {0:Number_Units_y}

// --------------------------------------------------------------------------

// ------------------------------------------------------------
// pillar 1f (full pillar)
// ------------------------------------------------------------

// *******************************
// Center
// *******************************
pc1_1f~{m}~{n} = newp; Point(pc1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, (2*tC+tD)/2,lcDielectricpillar};
pc2_1f~{m}~{n} = newp; Point(pc2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, -1*(tD)/2,lcDielectricpillar};
pc3_1f~{m}~{n} = newp; Point(pc3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tD/2,lcDielectricpillar};
pc4_1f~{m}~{n} = newp; Point(pc4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, -1*(2*tC+tD)/2,lcDielectricpillar};

// *******************************
// Dielectric pillar
// *******************************
// Top

pth1_1f~{m}~{n} = newp; Point(pth1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*r0, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, (2*tC+tD)/2,lcDielectricpillar};
pth2_1f~{m}~{n} = newp; Point(pth2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*r0, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, (2*tC+tD)/2,lcDielectricpillar};
pth3_1f~{m}~{n} = newp; Point(pth3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*r0, (2*tC+tD)/2,lcDielectricpillar};
pth4_1f~{m}~{n} = newp; Point(pth4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*r0, (2*tC+tD)/2,lcDielectricpillar};

cth1_1f~{m}~{n} = newc; Circle(cth1_1f~{m}~{n}) = {pth1_1f~{m}~{n}, pc1_1f~{m}~{n}, pth3_1f~{m}~{n}};
cth2_1f~{m}~{n} = newc; Circle(cth2_1f~{m}~{n}) = {pth3_1f~{m}~{n}, pc1_1f~{m}~{n}, pth2_1f~{m}~{n}};
cth3_1f~{m}~{n} = newc; Circle(cth3_1f~{m}~{n}) = {pth2_1f~{m}~{n}, pc1_1f~{m}~{n}, pth4_1f~{m}~{n}};
cth4_1f~{m}~{n} = newc; Circle(cth4_1f~{m}~{n}) = {pth4_1f~{m}~{n}, pc1_1f~{m}~{n}, pth1_1f~{m}~{n}};

// Bottom
pbh1_1f~{m}~{n} = newp; Point(pbh1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*r0, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, -1*(2*tC+tD)/2,lcDielectricpillar};
pbh2_1f~{m}~{n} = newp; Point(pbh2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*r0, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, -1*(2*tC+tD)/2,lcDielectricpillar};
pbh3_1f~{m}~{n} = newp; Point(pbh3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*r0, -1*(2*tC+tD)/2,lcDielectricpillar};
pbh4_1f~{m}~{n} = newp; Point(pbh4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*r0, -1*(2*tC+tD)/2,lcDielectricpillar};

// *******************************
// Upper Etching
// *******************************

// Bottom
ptue1_1f~{m}~{n} = newp; Point(ptue1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tD/2,lcEtchingpillar};
ptue2_1f~{m}~{n} = newp; Point(ptue2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, tD/2,lcEtchingpillar};
ptue3_1f~{m}~{n} = newp; Point(ptue3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*(r0+r1), tD/2,lcEtchingpillar};
ptue4_1f~{m}~{n} = newp; Point(ptue4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*(r0+r1), tD/2,lcEtchingpillar};

// Top
pbue1_1f~{m}~{n} = newp; Point(pbue1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*(r0+0.5*r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, (1.5*tC+tD)/2,lcEtchingpillar};
pbue2_1f~{m}~{n} = newp; Point(pbue2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*(r0+0.5*r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, (1.5*tC+tD)/2,lcEtchingpillar};
pbue3_1f~{m}~{n} = newp; Point(pbue3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a,pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*(r0+0.5*r1), (1.5*tC+tD)/2,lcEtchingpillar};
pbue4_1f~{m}~{n} = newp; Point(pbue4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a,pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*(r0+0.5*r1), (1.5*tC+tD)/2,lcEtchingpillar};

// Circular boundary
ctue1_1f~{m}~{n} = newc; Circle(ctue1_1f~{m}~{n}) = {ptue1_1f~{m}~{n}, pc3_1f~{m}~{n}, ptue3_1f~{m}~{n}};
ctue2_1f~{m}~{n} = newc; Circle(ctue2_1f~{m}~{n}) = {ptue3_1f~{m}~{n}, pc3_1f~{m}~{n}, ptue2_1f~{m}~{n}};
ctue3_1f~{m}~{n} = newc; Circle(ctue3_1f~{m}~{n}) = {ptue2_1f~{m}~{n}, pc3_1f~{m}~{n}, ptue4_1f~{m}~{n}};
ctue4_1f~{m}~{n} = newc; Circle(ctue4_1f~{m}~{n}) = {ptue4_1f~{m}~{n}, pc3_1f~{m}~{n}, ptue1_1f~{m}~{n}};

lue1_1f~{m}~{n} = newc; Line(lue1_1f~{m}~{n}) = {ptue1_1f~{m}~{n}, pth1_1f~{m}~{n}};
lue2_1f~{m}~{n} = newc; Line(lue2_1f~{m}~{n}) = {ptue2_1f~{m}~{n}, pth2_1f~{m}~{n}};
lue3_1f~{m}~{n} = newc; Line(lue3_1f~{m}~{n}) = {ptue3_1f~{m}~{n}, pth3_1f~{m}~{n}};
lue4_1f~{m}~{n} = newc; Line(lue4_1f~{m}~{n}) = {ptue4_1f~{m}~{n}, pth4_1f~{m}~{n}};

// *******************************
// Lower Etching
// *******************************

// Top
ptle1_1f~{m}~{n} = newp; Point(ptle1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, -1*(tD)/2,lcEtchingpillar};
ptle2_1f~{m}~{n} = newp; Point(ptle2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, -1*(tD)/2,lcEtchingpillar};
ptle3_1f~{m}~{n} = newp; Point(ptle3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*(r0+r1), -1*(tD)/2,lcEtchingpillar};
ptle4_1f~{m}~{n} = newp; Point(ptle4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*(r0+r1), -1*(tD)/2,lcEtchingpillar};

// Bottom
pble1_1f~{m}~{n} = newp; Point(pble1_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a-1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, -1*(2*tC+tD)/2,lcEtchingpillar};
pble2_1f~{m}~{n} = newp; Point(pble2_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a+1*(r0+r1), pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a, -1*(2*tC+tD)/2,lcEtchingpillar};
pble3_1f~{m}~{n} = newp; Point(pble3_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a-1*(r0+r1), -1*(2*tC+tD)/2,lcEtchingpillar};
pble4_1f~{m}~{n} = newp; Point(pble4_1f~{m}~{n}) = {pil_c_x+geo_f_x*a/2+pil_f_x*geo_f_x*m*a, pil_c_y+geo_f_y*a/2+pil_f_y*geo_f_y*n*a+1*(r0+r1), -1*(2*tC+tD)/2,lcEtchingpillar};

// Circular boundaries
ctle1_1f~{m}~{n} = newc; Circle(ctle1_1f~{m}~{n}) = {ptle1_1f~{m}~{n}, pc2_1f~{m}~{n}, ptle3_1f~{m}~{n}};
ctle2_1f~{m}~{n} = newc; Circle(ctle2_1f~{m}~{n}) = {ptle3_1f~{m}~{n}, pc2_1f~{m}~{n}, ptle2_1f~{m}~{n}};
ctle3_1f~{m}~{n} = newc; Circle(ctle3_1f~{m}~{n}) = {ptle2_1f~{m}~{n}, pc2_1f~{m}~{n}, ptle4_1f~{m}~{n}};
ctle4_1f~{m}~{n} = newc; Circle(ctle4_1f~{m}~{n}) = {ptle4_1f~{m}~{n}, pc2_1f~{m}~{n}, ptle1_1f~{m}~{n}};

cble1_1f~{m}~{n} = newc; Circle(cble1_1f~{m}~{n}) = {pble1_1f~{m}~{n}, pc4_1f~{m}~{n}, pble3_1f~{m}~{n}};
cble2_1f~{m}~{n} = newc; Circle(cble2_1f~{m}~{n}) = {pble3_1f~{m}~{n}, pc4_1f~{m}~{n}, pble2_1f~{m}~{n}};
cble3_1f~{m}~{n} = newc; Circle(cble3_1f~{m}~{n}) = {pble2_1f~{m}~{n}, pc4_1f~{m}~{n}, pble4_1f~{m}~{n}};
cble4_1f~{m}~{n} = newc; Circle(cble4_1f~{m}~{n}) = {pble4_1f~{m}~{n}, pc4_1f~{m}~{n}, pble1_1f~{m}~{n}};

lle1_1f~{m}~{n} = newc; Line(lle1_1f~{m}~{n}) = {ptle1_1f~{m}~{n}, pble1_1f~{m}~{n}};
lle2_1f~{m}~{n} = newc; Line(lle2_1f~{m}~{n}) = {ptle2_1f~{m}~{n}, pble2_1f~{m}~{n}};
lle3_1f~{m}~{n} = newc; Line(lle3_1f~{m}~{n}) = {ptle3_1f~{m}~{n}, pble3_1f~{m}~{n}};
lle4_1f~{m}~{n} = newc; Line(lle4_1f~{m}~{n}) = {ptle4_1f~{m}~{n}, pble4_1f~{m}~{n}};

// Lines connecting top and bottom
lconn5_1f~{m}~{n} = newc; Line(lconn5_1f~{m}~{n}) = {ptle1_1f~{m}~{n}, ptue1_1f~{m}~{n}};
lconn6_1f~{m}~{n} = newc; Line(lconn6_1f~{m}~{n}) = {ptle2_1f~{m}~{n}, ptue2_1f~{m}~{n}};
lconn7_1f~{m}~{n} = newc; Line(lconn7_1f~{m}~{n}) = {ptle3_1f~{m}~{n}, ptue3_1f~{m}~{n}};
lconn8_1f~{m}~{n} = newc; Line(lconn8_1f~{m}~{n}) = {ptle4_1f~{m}~{n}, ptue4_1f~{m}~{n}};

// --------------------------------------------------------------------------

// *******************************************************
// Copper planes
// *******************************************************

// Connect the upper and lower points with lines to form the plate

// Connect the upper and lower points with lines to form the plate
lcpbib17~{m}~{n} = newc; Line(lcpbib17~{m}~{n}) = {ptle1_1f~{m}~{n}, pble1_1f~{m}~{n}};
lcpbib18~{m}~{n} = newc; Line(lcpbib18~{m}~{n}) = {ptle2_1f~{m}~{n}, pble2_1f~{m}~{n}};
lcpbib19~{m}~{n} = newc; Line(lcpbib19~{m}~{n}) = {ptle3_1f~{m}~{n}, pble3_1f~{m}~{n}};
lcpbib20~{m}~{n} = newc; Line(lcpbib20~{m}~{n}) = {ptle4_1f~{m}~{n}, pble4_1f~{m}~{n}};

// *************************************************
// Define surfaces
// *************************************************

// Copper plate surfaces

ll_side_gas1a = newreg; Line Loop(ll_side_gas1a) = {lcptlb5a, lmid1_2, -lcpblb5a, -lcorner1};
ps_side_gas1a = newreg; Plane Surface(ps_side_gas1a) = {ll_side_gas1a};
ll_side_gas2a = newreg; Line Loop(ll_side_gas2a) = {lcptlb6a, lmid2_2, -lcpblb6a, -lcorner2};
ps_side_gas2a = newreg; Plane Surface(ps_side_gas2a) = {ll_side_gas2a};
ll_side_gas3a = newreg; Line Loop(ll_side_gas3a) = {lcptlb7a, lmid3_2, -lcpblb7a, -lcorner3};
ps_side_gas3a = newreg; Plane Surface(ps_side_gas3a) = {ll_side_gas3a};
ll_side_gas4a = newreg; Line Loop(ll_side_gas4a) = {lcptlb8a, lmid4_2, -lcpblb8a, -lcorner4};
ps_side_gas4a = newreg; Plane Surface(ps_side_gas4a) = {ll_side_gas4a};
ll_side_gas1b = newreg; Line Loop(ll_side_gas1b) = {lcptlb5b, lcorner2, -lcpblb5b, -lmid1_2};
ps_side_gas1b = newreg; Plane Surface(ps_side_gas1b) = {ll_side_gas1b};
ll_side_gas2b = newreg; Line Loop(ll_side_gas2b) = {lcptlb6b, lcorner3, -lcpblb6b, -lmid2_2};
ps_side_gas2b = newreg; Plane Surface(ps_side_gas2b) = {ll_side_gas2b};
ll_side_gas3b = newreg; Line Loop(ll_side_gas3b) = {lcptlb7b, lcorner4, -lcpblb7b, -lmid3_2};
ps_side_gas3b = newreg; Plane Surface(ps_side_gas3b) = {ll_side_gas3b};
ll_side_gas4b = newreg; Line Loop(ll_side_gas4b) = {lcptlb8b, lcorner1, -lcpblb8b, -lmid4_2};
ps_side_gas4b = newreg; Plane Surface(ps_side_gas4b) = {ll_side_gas4b};

// Surfaces to which voltages will be applied

llcp_low_rim_1a~{m}~{n} = newreg; Line Loop(llcp_low_rim_1a~{m}~{n}) = {-lle1_1f~{m}~{n}, -cble1_1f~{m}~{n}, ctle1_1f~{m}~{n}, lle3_1f~{m}~{n}};
llcp_low_rim_1b~{m}~{n} = newreg; Line Loop(llcp_low_rim_1b~{m}~{n}) = {-lle3_1f~{m}~{n}, -cble2_1f~{m}~{n}, ctle2_1f~{m}~{n}, lle2_1f~{m}~{n}};
llcp_low_rim_1c~{m}~{n} = newreg; Line Loop(llcp_low_rim_1c~{m}~{n}) = {-lle2_1f~{m}~{n}, -cble3_1f~{m}~{n}, ctle3_1f~{m}~{n}, lle4_1f~{m}~{n}};
llcp_low_rim_1d~{m}~{n} = newreg; Line Loop(llcp_low_rim_1d~{m}~{n}) = {-lle4_1f~{m}~{n}, -cble4_1f~{m}~{n}, ctle4_1f~{m}~{n}, lle1_1f~{m}~{n}};

// Surfaces to which voltages will be applied

ps_lower_cp1 = newreg; Surface(newreg) = {llcp_low_rim_1a~{m}~{n}};
surf_lower_cp[] += ps_lower_cp1;
ps_lower_cp2 = newreg; Surface(newreg) = {llcp_low_rim_1b~{m}~{n}};
surf_lower_cp[] += ps_lower_cp2;
ps_lower_cp3 = newreg; Surface(newreg) = {llcp_low_rim_1c~{m}~{n}};
surf_lower_cp[] += ps_lower_cp3;
ps_lower_cp4 = newreg; Surface(newreg) = {llcp_low_rim_1d~{m}~{n}};
surf_lower_cp[] += ps_lower_cp4;

// Gas & dielectric surfaces

ll_cyl_dielectric1b~{m}~{n} = newreg; Line Loop(ll_cyl_dielectric1b~{m}~{n}) = {lconn5_1f~{m}~{n}, ctue1_1f~{m}~{n}, -lconn7_1f~{m}~{n}, -ctle1_1f~{m}~{n}};
ps_cyl_dielectric1 = newreg; Surface(newreg) = {ll_cyl_dielectric1b~{m}~{n}};
surf_cyl_dielectric[] += ps_cyl_dielectric1;
ll_cyl_dielectric2b~{m}~{n} = newreg; Line Loop(ll_cyl_dielectric2b~{m}~{n}) = {lconn7_1f~{m}~{n}, ctue2_1f~{m}~{n}, -lconn6_1f~{m}~{n}, -ctle2_1f~{m}~{n}};
ps_cyl_dielectric2 = newreg; Surface(newreg) = {ll_cyl_dielectric2b~{m}~{n}};
surf_cyl_dielectric[] += ps_cyl_dielectric2; 
ll_cyl_dielectric3b~{m}~{n} = newreg; Line Loop(ll_cyl_dielectric3b~{m}~{n}) = {lconn6_1f~{m}~{n}, ctue3_1f~{m}~{n}, -lconn8_1f~{m}~{n}, -ctle3_1f~{m}~{n}};
ps_cyl_dielectric3 = newreg; Surface(newreg) = {ll_cyl_dielectric3b~{m}~{n}}; 
surf_cyl_dielectric[] += ps_cyl_dielectric3;
ll_cyl_dielectric4b~{m}~{n} = newreg; Line Loop(ll_cyl_dielectric4b~{m}~{n}) = {lconn8_1f~{m}~{n}, ctue4_1f~{m}~{n}, -lconn5_1f~{m}~{n}, -ctle4_1f~{m}~{n}};
ps_cyl_dielectric4 = newreg; Surface(newreg) = {ll_cyl_dielectric4b~{m}~{n}};
surf_cyl_dielectric[] += ps_cyl_dielectric4;

ll_top_cp1a = newreg; Line Loop(newreg) = {ctle1_1f~{m}~{n}, ctle2_1f~{m}~{n}, ctle3_1f~{m}~{n}, ctle4_1f~{m}~{n}}; 
ll_top_cp2a[] += {ll_top_cp1a};
ps_top_cp1b = news; Plane Surface(news) = {ll_top_cp1a};
surf_top_cp1b[] += {ps_top_cp1b};

ll_bottom_cp1a = newreg; Line Loop(newreg) = {cble1_1f~{m}~{n}, cble2_1f~{m}~{n}, cble3_1f~{m}~{n}, cble4_1f~{m}~{n}};
ll_bottom_cp2a[] += {ll_bottom_cp1a};
ps_bottom_cp1b = news; Plane Surface(news) = {ll_bottom_cp1a};
surf_bottom_cp1b[] += {ps_bottom_cp1b};

ll_top_gas2~{m}~{n} = newreg; Line Loop(ll_top_gas2~{m}~{n}) = {cth1_1f~{m}~{n}, cth2_1f~{m}~{n}, cth3_1f~{m}~{n}, cth4_1f~{m}~{n}};
ps_top_gas1 = news; Plane Surface(news) = {ll_top_gas2~{m}~{n}};
surf_top_gas1[] += ps_top_gas1;

ll_top_gas3~{m}~{n} = newreg; Line Loop(ll_top_gas3~{m}~{n}) = {lue1_1f~{m}~{n}, cth1_1f~{m}~{n}, -lue3_1f~{m}~{n}, -ctue1_1f~{m}~{n}};
ps_top_gas2 = news; Surface(news) = {ll_top_gas3~{m}~{n}};
surf_top_gas2[] += ps_top_gas2;

ll_top_gas4~{m}~{n} = newreg; Line Loop(ll_top_gas4~{m}~{n}) = {lue2_1f~{m}~{n}, -cth2_1f~{m}~{n}, -lue3_1f~{m}~{n}, ctue2_1f~{m}~{n}};
ps_top_gas3 = news; Surface(news) = {ll_top_gas4~{m}~{n}};
surf_top_gas3[] += ps_top_gas3;

ll_top_gas5~{m}~{n} = newreg; Line Loop(ll_top_gas5~{m}~{n}) = {lue2_1f~{m}~{n}, cth3_1f~{m}~{n}, -lue4_1f~{m}~{n}, -ctue3_1f~{m}~{n}};
ps_top_gas4 = news; Surface(news) = {ll_top_gas5~{m}~{n}};
surf_top_gas4[] += ps_top_gas4;

ll_top_gas6~{m}~{n} = newreg; Line Loop(ll_top_gas6~{m}~{n}) = {lue4_1f~{m}~{n}, cth4_1f~{m}~{n}, -lue1_1f~{m}~{n}, -ctue4_1f~{m}~{n}};
ps_top_gas5 = news; Surface(news) = {ll_top_gas6~{m}~{n}};
surf_top_gas5[] += ps_top_gas5;

 EndFor
EndFor

ll_top_cp3a = newreg; Line Loop(ll_top_cp3a) = {lcpblb5a, lcpblb5b, lcpblb6a, lcpblb6b, lcpblb7a, lcpblb7b, lcpblb8a, lcpblb8b};
ll_top_cp2a[] += ll_top_cp3a;
ps_top_cp1a = news; Plane Surface(news) = {ll_top_cp2a[]};
surf_top_cp[] = {ps_top_cp1a};

ll_bottom_cp3a = newreg; Line Loop(ll_bottom_cp3a) = {lcpbub1a, lcpbub1b, lcpbub2a, lcpbub2b, lcpbub3a, lcpbub3b, lcpbub4a, lcpbub4b};
ll_bottom_cp2a[] += ll_bottom_cp3a;
ps_bottom_cp1a = news; Plane Surface(news) = {ll_bottom_cp2a[]};
surf_bottom_cp[] = {ps_bottom_cp1a};

// Surface loops

sl_dielectric = newreg; Surface Loop(newreg) = { surf_top_gas1[], surf_top_gas2[], surf_top_gas3[], surf_top_gas4[], surf_top_gas5[], surf_cyl_dielectric[], surf_lower_cp[], surf_bottom_cp1b[] };
total_sl_dielectric[] += sl_dielectric;

sl_gas = newreg; Surface Loop(newreg) = { surf_top_cp[], surf_top_gas1[], surf_top_gas2[], surf_top_gas3[], surf_top_gas4[], surf_top_gas5[], surf_cyl_dielectric[], ps_side_gas1a, ps_side_gas2a, ps_side_gas3a, ps_side_gas4a, ps_side_gas1b, ps_side_gas2b, ps_side_gas3b, ps_side_gas4b, ps_bsurf1, ps_bsurf2, ps_bsurf3, ps_bsurf4, ps_bsurf5, ps_bsurf6, ps_bsurf7, pscp_up_border1, pscp_up_border2, pscp_up_border3, pscp_up_border4, st_s1_x1[], st_x1_1a1[], st_x1_1a2[], st_x1_1b1[], st_x1_1b2[], st_tmp_x1[], st_s1_x2[], st_x2_1b1[], st_x2_1b2[], st_x2_1a1[], st_x2_1a2[], st_tmp_x2[], st_s1_y1[], st_y1_2a1[], st_y1_2a2[], st_y1_2b1[], st_y1_2b2[], st_tmp_y1[], st_s1_y2[], st_y2_2b1[], st_y2_2b2[], st_y2_2a1[], st_y2_2a2[], st_tmp_y2[] };
total_sl_gas[] += sl_gas;

sl_lower_cp = newreg; Surface Loop(newreg) = { surf_top_cp[], surf_lower_cp[], surf_bottom_cp[], pscp_low_border1, pscp_low_border2, pscp_low_border3, pscp_low_border4 };
total_sl_lower_cp[] += sl_lower_cp;

// Volumes

vol_dielectric = newreg; Volume(vol_dielectric) = { total_sl_dielectric[] };
vol_gas = newreg; Volume(vol_gas) = { total_sl_gas[] };
vol_lower_cp = newreg; Volume(vol_lower_cp) = { total_sl_lower_cp[] };
//vol_wire = newreg; Volume(vol_wire) = {sl_wire};

// Physical surfaces

// Surfaces for periodic boundary conditions

Physical Surface(physsurf_gas) = { surf_top_cp[], surf_top_gas1[], surf_top_gas2[], surf_top_gas3[], surf_top_gas4[], surf_top_gas5[], surf_cyl_dielectric[], ps_side_gas1a, ps_side_gas2a, ps_side_gas3a, ps_side_gas4a, ps_side_gas1b, ps_side_gas2b, ps_side_gas3b, ps_side_gas4b, ps_bsurf1, ps_bsurf2, ps_bsurf3, ps_bsurf4, ps_bsurf5, ps_bsurf6, ps_bsurf7, pscp_up_border1, pscp_up_border2, pscp_up_border3, pscp_up_border4, st_s1_x1[], st_x1_1a1[], st_x1_1a2[], st_tmp_x1[], st_s1_x2[], st_x2_1b1[], st_x2_1b2[], st_tmp_x2[], st_s1_y1[], st_y1_2a1[], st_y1_2a2[], st_tmp_y1[], st_s1_y2[], st_y2_2b1[], st_y2_2b2[], st_tmp_y2[] };

Physical Surface(physsurf_bd1h1) = {pscp_up_border1, ps_side_gas1a, ps_side_gas1b, ps_bsurf1, ps_bsurf2}; 
Physical Surface(physsurf_bd1h2) = {pscp_up_border2, ps_side_gas2a, ps_side_gas2b, ps_bsurf3};
Physical Surface(physsurf_bd2h1) = {pscp_up_border3, ps_side_gas3a, ps_side_gas3b, ps_bsurf4, ps_bsurf5};
Physical Surface(physsurf_bd2h2) = {pscp_up_border4, ps_side_gas4a, ps_side_gas4b, ps_bsurf6};

// Physical surfaces to elements

Physical Surface(physsurf_upper_el) = {ps_bsurf7};
Physical Surface(physsurf_lower_cp) = { surf_top_cp[], surf_lower_cp[], surf_bottom_cp[], pscp_low_border1, pscp_low_border2, pscp_low_border3, pscp_low_border4 };
Physical Surface(physsurf_dielectric) = { surf_top_gas1[], surf_top_gas2[], surf_top_gas3[], surf_top_gas4[], surf_top_gas5[], surf_cyl_dielectric[], surf_lower_cp[], surf_bottom_cp1b[] };

// Volume

Physical Volume(physvol_dielectric) = { vol_dielectric };
Physical Volume(physvol_gas) = { vol_gas };
Physical Volume(physvol_lower_cp) = { vol_lower_cp };

Coherence;
Geometry.AutoCoherence = 1;

