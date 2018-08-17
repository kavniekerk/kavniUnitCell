//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/// MESHING OPTIONS MODULE

//----------------------------------------------------------
// Mesh parameters - characterization of mesh

// Recombine Surface "*";                        // Recombine all surfaces
// Mesh.RecombineAll = 1;                        // Recombine all surfaces
General.ExpertMode = 1;                       // Don't complain for hybrid structured/unstructured mesh
Mesh.Optimize = 1;                            // Optimize 3D tet mesh 
Mesh.Algorithm = 5;                           // meshing 2D with MeshAdapt - for complex curved surfaces, 2D mesh algorithm (1=MeshAdapt, 5=Delaunay, 6=Frontal), Default value: 1, Saved in: General.OptionsFileName
Mesh.Algorithm3D = 3;                         // meshing 3D with Frontal - robustness, quality, 3D mesh algorithm (1=Delaunay, 4=Frontal), Default value: 1, Saved in: General.OptionsFileName
Mesh.AngleSmoothNormals = 1;                  // Threshold angle below which normals are not smoothed, Default value: 30, Saved in: General.OptionsFileName
Mesh.AllowSwapAngle = 5;                      // Treshold angle (in degrees) between faces normals under which we allow an edge swap, Default value: 10, Saved in: General.OptionsFileName
Mesh.BdfFieldFormat = 0;                      // Field format for Nastran BDF files (0=free, 1=small, 2=large), Default value: 1, Saved in: General.OptionsFileName
Mesh.Binary = 0;                              // Write mesh files in binary format (if possible), Default value: 0, Saved in: General.OptionsFileName

Mesh.CgnsImportOrder = 1;                     // Enable the creation of high-order mesh from CGNS structured meshes.(1, 2, 4, 8, ...), Default value: 1, Saved in: General.OptionsFileName

Mesh.ChacoArchitecture = 1;                   // (Adv. Chaco): Parallel architecture topology (0=hypercube, 1-3=mesh dimensions), Default value: 1, Saved in: General.OptionsFileName
Mesh.ChacoEigensolver = 1;                    // (Adv. Chaco): Type of eigensolver for a spectral algorithm (0=Lanczos, 1=Multilevel RQI/Symmlq), Default value: 1, Saved in: General.OptionsFileName
Mesh.ChacoEigTol = 1;                         // (Adv. Chaco): Tolerance of the eigensolver for spectral or multilevel-KL algorithms, Default value: 0.001, Saved in: General.OptionsFileName
Mesh.ChacoGlobalMethod = 1;                   // Chaco partitioning algorithm (1=Multilevel-KL, 2=Spectral, 4=Linear, 5=Random, 6=Scattered), Default value: 1, Saved in: General.OptionsFileName
Mesh.ChacoHypercubeDim = 0;                   // (Adv. Chaco): Dimensional partitioning for a hypercube topology, Default value: 0, Saved in: General.OptionsFileName
Mesh.ChacoLocalMethod = 1;                    // (Adv. Chaco): Local partitioning algorithm, Default value: 1, Saved in: General.OptionsFileName
Mesh.ChacoMeshDim1 = 1;                       // (Adv. Chaco): Number of partitions in the first dimension of a mesh topology, Default value: 1, Saved in: General.OptionsFileName
Mesh.ChacoMeshDim2 = 1;                       // (Adv. Chaco): Number of partitions in the second dimension of a mesh topology, Default value: 1, Saved in: General.OptionsFileName
Mesh.ChacoMeshDim3 = 1;                       // (Adv. Chaco): Number of partitions in the third dimension of a mesh topology, Default value: 1, Saved in: General.OptionsFileName
Mesh.ChacoPartitionSection = 1;               // (Adv. Chaco): Partition by (1=bisection, 2=quadrisection, 3=octasection, Default value: 1, Saved in: General.OptionsFileName
Mesh.ChacoSeed = 7.65432e+06;                 // (Adv. Chaco): Seed for random number generator, Default value: 7.65432e+06, Saved in: General.OptionsFileName
Mesh.ChacoVMax = 250;                         // (Adv. Chaco): Maximum vertices in a coarse graph (for multilevel-KL algorithm and Multilevel RQI/Symmlq eigensolver), Default value: 250, Saved in: General.OptionsFileName
Mesh.ChacoParamINTERNAL_VERTICES = 0;         // (Adv. Chaco): Parameter INTERNAL_VERTICES, Default value: 0, Saved in: General.OptionsFileName
Mesh.ChacoParamREFINE_MAP = 1;                // (Adv. Chaco): Parameter REFINE_MAP, Default value: 1, Saved in: General.OptionsFileName
Mesh.ChacoParamREFINE_PARTITION = 0;          // (Adv. Chaco): Parameter REFINE_PARTITION, Default value: 0, Saved in: General.OptionsFileName
Mesh.ChacoParamTERMINAL_PROPOGATION = 0;      // (Adv. Chaco): Parameter TERMINAL_PROPOGATION, Default value: 0, Saved in: General.OptionsFileName

Mesh.CharacteristicLengthExtendFromBoundary = 1; // Extend characteristic lengths from the boundaries inside the surface/volume, Default value: 1, Saved in: General.OptionsFileName
Mesh.CharacteristicLengthFactor = 1;          // Factor applied to all characteristic lengths, Default value: 1, Saved in: General.OptionsFileName
Mesh.CharacteristicLengthMin = 0;             // Minimum characteristic length, Default value: 0, Saved in: General.OptionsFileName
Mesh.CharacteristicLengthMax = 1e+22;         // Maximum characteristic length, Default value: 1e+22, Saved in: General.OptionsFileName
Mesh.CharacteristicLengthFromCurvature = 0;   // Compute characteristic lengths from curvature, Default value: 0, Saved in: General.OptionsFileName
Mesh.CharacteristicLengthFromPoints = 1;      // Compute characteristic lengths from values given at geometry points, Default value: 1, Saved in: General.OptionsFileName

Mesh.Clip = 0;                                // Enable clipping planes? (Plane[i]=2^i, i=0,...,5), Default value: 0, Saved in: - 
Mesh.ColorCarousel = 1;                       // Mesh coloring (0=by element type, 1=by elementary entity, 2=by physical entity, 3=by partition), Default value: 1, Saved in: General.OptionsFileName
Mesh.CpuTime = 0;                             // CPU time (in seconds) for the generation of the current mesh (read-only), Default value: 0, Saved in: -
Mesh.DrawSkinOnly = 0;                        // Draw only the skin of 3D meshes?, Default value: 0, Saved in: General.OptionsFileName
Mesh.Dual = 0;                                // Display the dual mesh obtained by barycentric subdivision, Default value: 0, Saved in: General.OptionsFileName
Mesh.ElementOrder = 1;                        // Element order (1=linear elements, N (<6) = elements of higher order), Default value: 1, Saved in: General.OptionsFileName
Mesh.Explode = 1;                             // Element shrinking factor (between 0 and 1), Default value: 1, Saved in: General.OptionsFileName
Mesh.Format = 1;                              // Mesh output format (1=msh, 2=unv, 19=vrml, 27=stl, 30=mesh, 31=bdf, 32=cgns, 33=med), Default value: 1, Saved in: General.OptionsFileName
Mesh.Hexahedra = 1;                           // Display mesh hexahedra?, Default value: 1, Saved in: General.OptionsFileName
// Mesh.LabelsFrequency = 100;                   // Labels display frequency?, Default value: 100, Saved in: General.OptionsFileName
Mesh.LabelType = 0;                           // Type of element label (0=element number, 1=elementary entity number, 2=physical entity number, 3=partition number, 4=coordinates), Default value: 0, Saved in: General.OptionsFileName
Mesh.LcIntegrationPrecision = 1e-05;          // Accuracy of evaluation of the LC field for 1D mesh generation, Default value: 1e-09, Saved in: General.OptionsFileName
Mesh.Light = 1;                               // Enable lighting for the mesh, Default value: 1, Saved in: General.OptionsFileName
Mesh.LightLines = 1;                          // Enable lighting for mesh lines (element edges), Default value: 1, Saved in: General.OptionsFileName
Mesh.LightTwoSide = 1;                        // Light both sides of surfaces (leads to slower rendering), Default value: 1, Saved in: General.OptionsFileName
Mesh.Lines = 0;                               // Display mesh lines (1D elements)?, Default value: 0, Saved in: General.OptionsFileName
Mesh.LineNumbers = 0;                         // Display mesh line numbers?, Default value: 0, Saved in: General.OptionsFileName
Mesh.LineWidth = 1;                           // Display width of mesh lines (in pixels), Default value: 1, Saved in: General.OptionsFileName
Mesh.MeshOnlyVisible = 0;                     // Mesh only visible entities (experimental: use with caution!), Default value: 0, Saved in: General.OptionsFileName

Mesh.MetisAlgorithm = 1;                      // METIS partitioning algorithm (1=Recursive, 2=K-way, 3=Nodal weight), Default value: 1, Saved in: General.OptionsFileName
Mesh.MetisEdgeMatching = 3;                   // (Adv. METIS): Determines the matching type (1=Random, 2=Heavy-Edge, 3=Sorted Heavy-Edge), Default value: 3, Saved in: General.OptionsFileName
Mesh.MetisRefinementAlgorithm = 3;            // (Adv. METIS): Algorithm for k-way refinement (1=Random, 2=Greedy, 3=Random with minimized connectivity), Default value: 3, Saved in: General.OptionsFileName

Mesh.MinimumCirclePoints = 2;                 // Minimum number of points used to mesh a circle, Default value: 7, Saved in: General.OptionsFileName
Mesh.MinimumCurvePoints = 2;                  // Minimum number of points used to mesh a (non-straight) curve, Default value: 3, Saved in: General.OptionsFileName
Mesh.MshFileVersion = 2.2;                    // Version of the MSH file format to use, Default value: 2.2, Saved in: General.OptionsFileName
Mesh.MshFilePartitioned = 0;                  // Split MSH file by mesh partition, Default value: 0, Saved in: General.OptionsFileName
// Mesh.MultiplePassesMeshes = 0;                // Do a first simple mesh and use it for complex background meshes (curvatures...), Default value: 0, Saved in: General.OptionsFileName
Mesh.PartitionHexWeight = 1;                  // Weight of hexahedral element for METIS load balancing, Default value: 1, Saved in: General.OptionsFileName
Mesh.PartitionPrismWeight = 1;                // Weight of prismatic element (wedge) for METIS load balancing, Default value: 1, Saved in: General.OptionsFileName
Mesh.PartitionPyramidWeight = 1;              // Weight of pyramidal element for METIS load balancing, Default value: 1, Saved in: General.OptionsFileName
Mesh.PartitionQuadWeight = 1;                 // Weight of quadrangle for METIS load balancing, Default value: 1, Saved in: General.OptionsFileName
Mesh.PartitionTetWeight = 1;                  // Weight of tetrahedral element for METIS load balancing, Default value: 1, Saved in: General.OptionsFileName
Mesh.PartitionTriWeight = 1;                  // Weight of triangle for METIS load balancing, Default value: 1, Saved in: General.OptionsFileName
Mesh.NbHexahedra = 0;                         // Number of hexahedra in the current mesh (read-only), Default value: 0, Saved in: -
Mesh.NbNodes = 0;                             // Number of nodes in the current mesh (read-only), Default value: 0, Saved in: -
Mesh.NbPartitions = 4;                        // Number of partitions, Default value: 4, Saved in: General.OptionsFileName
Mesh.NbPrisms = 0;                            // Number of prisms in the current mesh (read-only), Default value: 0, Saved in: -
Mesh.NbPyramids = 0;                          // Number of pyramids in the current mesh (read-only), Default value: 0, Saved in: -
Mesh.NbQuadrangles = 0;                       // Number of quadrangles in the current mesh (read-only), Default value: 0, Saved in: -
Mesh.NbTetrahedra = 0;                        // Number of tetrahedra in the current mesh (read-only), Default value: 0, Saved in: -
Mesh.NbTriangles = 0;                         // Number of triangles in the current mesh (read-only), Default value: 0, Saved in: -
Mesh.Normals = 0;                             // Display size of normal vectors (in pixels), Default value: 0, Saved in: General.OptionsFileName
Mesh.NumSubEdges = 2;                         // Number of edge subdivisions when displaying high order elements, Default value: 2, Saved in: General.OptionsFileName
Mesh.Optimize = 0;                            // Optimize the mesh to improve the quality of tetrahedral elements, Default value: 0, Saved in: General.OptionsFileName
Mesh.OptimizeNetgen = 0;                      // Optimize the mesh using Netgen to improve the quality of tetrahedral elements, Default value: 0, Saved in: General.OptionsFileName
Mesh.Partitioner = 2;                         // Partitioner software (1=Chacho, 2=METIS), Default value: 2, Saved in: General.OptionsFileName
Mesh.Points = 0;                              // Display mesh vertices (nodes)?, Default value: 0, Saved in: General.OptionsFileName
Mesh.PointNumbers = 0;                        // Display mesh node numbers?, Default value: 0, Saved in: General.OptionsFileName
Mesh.PointSize = 4;                           // Display size of mesh vertices (in pixels), Default value: 4, Saved in: General.OptionsFileName
Mesh.PointType = 0;                           // Display mesh vertices as solid color dots (0) or 3D spheres (1), Default value: 0, Saved in: General.OptionsFileName
Mesh.Prisms = 1;                              // Display mesh prisms?, Default value: 1, Saved in: General.OptionsFileName
Mesh.Pyramids = 1;                            // Display mesh pyramids?, Default value: 1, Saved in: General.OptionsFileName
Mesh.Quadrangles = 1;                         // Display mesh quadrangles?, Default value: 1, Saved in: General.OptionsFileName
Mesh.QualityInf = 0;                          // Only display elements whose quality measure is greater than QualityInf, Default value: 0, Saved in: General.OptionsFileName
Mesh.QualitySup = 0;                          // Only display elements whose quality measure is smaller than QualitySup, Default value: 0, Saved in: General.OptionsFileName
Mesh.QualityType = 2;                         // Type of quality measure (0=gamma~vol/sum_face/max_edge, 1=eta~vol^(2/3)/sum_edge^2, 2=rho~min_edge/max_edge), Default value: 2, Saved in: General.OptionsFileName
Mesh.RadiusInf = 0;                           // Only display elements whose longest edge is greater than RadiusInf, Default value: 0, Saved in: General.OptionsFileName
Mesh.RadiusSup = 0;                           // Only display elements whose longest edge is smaller than RadiusSup, Default value: 0, Saved in: General.OptionsFileName
Mesh.RandomFactor = 1e-05;                    // Random factor used in the 2D meshing algorithm (should be increased if RandomFactor * size(triangle)/size(model) approaches machine accuracy), Default value: 1e-09, Saved in: General.OptionsFileName

Mesh.RecombinationAlgorithm = 1;              // Mesh recombination algorithm (0=standard, 1=blossom), Default value: 1, Saved in: General.OptionsFileName
Mesh.RecombineAll = 0;                        // Apply recombination algorithm to all surfaces, ignoring per-surface spec, Default value: 0, Saved in: General.OptionsFileName
Mesh.RemeshAlgorithm = 0;                     // Remeshing algorithm (0=no split, 1=automatic, 2=automatic only with metis), Default value: 0, Saved in: General.OptionsFileName
Mesh.RemeshParametrization = 0;               // Remsh Parametrization (0=harmonic, 1=conformal, 2=rbf harmonic), Default value: 0, Saved in: General.OptionsFileName

Mesh.RefineSteps = 10;                        // Number of refinement steps in the MeshAdapt-based 2D algorithms, Default value: 10, Saved in: General.OptionsFileName
// Mesh.ReverseAllNormals = 0;                   // Reverse all the mesh normals (for display), Default value: 0, Saved in: General.OptionsFileName
Mesh.SaveAll = 0;                             // Ignore Physical definitions and save all elements, Default value: 0, Saved in: -
Mesh.SaveElementTagType = 1;                  // Type of the element tag saved in mesh formats that don't support saving physical or partition ids (1=elementary, 2=physical, 3=partition), Default value: 1, Saved in: General.OptionsFileName
Mesh.SaveParametric = 0;                      // Save parametric coordinates of nodes, Default value: 0, Saved in: General.OptionsFileName
Mesh.SaveGroupsOfNodes = 0;                   // Save groups of nodes for each physical line and surface (UNV mesh format only), Default value: 0, Saved in: General.OptionsFileName
Mesh.ScalingFactor = 1;                       // Global scaling factor applied to the saved mesh, Default value: 1, Saved in: General.OptionsFileName
Mesh.SecondOrderExperimental = 0;             // Use experimental code to generate second order mesh, Default value: 0, Saved in: General.OptionsFileName
Mesh.SecondOrderIncomplete = 1;               // Create incomplete second order elements? (8-node quads, 20-node hexas, etc.), Default value: 1, Saved in: General.OptionsFileName
Mesh.SecondOrderLinear = 0;                   // Should second order vertices simply be created by linear interpolation?, Default value: 0, Saved in: General.OptionsFileName
Mesh.Smoothing = 1;                           // Number of smoothing steps applied to the final mesh, Default value: 1, Saved in: General.OptionsFileName
// Mesh.SmoothInternalEdges = 0;                 // Number of smoothing steps of internal edges for high order meshes, Default value: 0, Saved in: General.OptionsFileName
Mesh.SmoothNormals = 0;                       // Smooth the mesh normals?, Default value: 0, Saved in: General.OptionsFileName
Mesh.SubdivisionAlgorithm = 0;                // Mesh subdivision algorithm (0=none, 1=all quadrangles, 2=all hexahedra), Default value: 0, Saved in: General.OptionsFileName
Mesh.SurfaceEdges = 1;                        // Display edges of surface mesh?, Default value: 1, Saved in: General.OptionsFileName
Mesh.SurfaceFaces = 0;                        // Display faces of surface mesh?, Default value: 0, Saved in: General.OptionsFileName
Mesh.SurfaceNumbers = 0;                      // Display surface mesh element numbers?, Default value: 0, Saved in: General.OptionsFileName
Mesh.SwitchElementTags = 0;                   // Invert elementary and physical tags when reading the mesh, Default value: 0, Saved in: General.OptionsFileName
Mesh.Tangents = 0;                            // Display size of tangent vectors (in pixels), Default value: 0, Saved in: General.OptionsFileName
Mesh.Tetrahedra = 0;                          // Display mesh tetrahedra?, Default value: 1, Saved in: General.OptionsFileName
Mesh.ToleranceEdgeLength = 0;                 // Skip a model edge in mesh generation if its length is less than user's defined tolerance, Default value: 0, Saved in: General.OptionsFileName
Mesh.Triangles = 1;                           // Display mesh triangles?, Default value: 1, Saved in: General.OptionsFileName
Mesh.VolumeEdges = 1;                         // Display edges of volume mesh?, Default value: 1, Saved in: General.OptionsFileName
Mesh.VolumeFaces = 0;                         // Display faces of volume mesh?, Default value: 0, Saved in: General.OptionsFileName
Mesh.VolumeNumbers = 0;                       // Display volume mesh element numbers?, Default value: 0, Saved in: General.OptionsFileName
Mesh.Voronoi = 0;                             // Display the voronoi diagram, Default value: 0, Saved in: General.OptionsFileName
Mesh.ZoneDefinition = 0;                      // Method for defining a zone (0=single zone, 1=by partition, 2=by physical), Default value: 0, Saved in: General.OptionsFileName

//----------------------------------------------------------
// Mesh color scheme

Mesh.Color.Points = {0,0,255};                // Mesh node color, Default value: {0,0,255}, Saved in: General.OptionsFileName
Mesh.Color.PointsSup = {255,0,255};           // Second order mesh node color, Default value: {255,0,255}, Saved in: General.OptionsFileName
Mesh.Color.Lines = {0,0,0};                   // Mesh line color, Default value: {0,0,0}, Saved in: General.OptionsFileName
Mesh.Color.Triangles = {160,150,255};         // Mesh triangle color (if Mesh.ColorCarousel=0), Default value: {160,150,255}, Saved in: General.OptionsFileName
Mesh.Color.Quadrangles = {130,120,225};       // Mesh quadrangle color (if Mesh.ColorCarousel=0), Default value: {130,120,225}, Saved in: General.OptionsFileName
Mesh.Color.Tetrahedra = {160,150,255};        // Mesh tetrahedron color (if Mesh.ColorCarousel=0), Default value: {160,150,255}, Saved in: General.OptionsFileName
Mesh.Color.Hexahedra = {130,120,225};         // Mesh hexahedron color (if Mesh.ColorCarousel=0), Default value: {130,120,225}, Saved in: General.OptionsFileName
Mesh.Color.Prisms = {232,210,23};             // Mesh prism color (if Mesh.ColorCarousel=0), Default value: {232,210,23}, Saved in: General.OptionsFileName
Mesh.Color.Pyramids = {217,113,38};           // Mesh pyramid color (if Mesh.ColorCarousel=0), Default value: {217,113,38}, Saved in: General.OptionsFileName
Mesh.Color.Tangents = {255,255,0};            // Tangent mesh vector color, Default value: {255,255,0}, Saved in: General.OptionsFileName
Mesh.Color.Normals = {255,0,0};               // Normal mesh vector color, Default value: {255,0,0}, Saved in: General.OptionsFileName
Mesh.Color.Zero = {255,120,0};                // Color 0 in color carousel, Default value: {255,120,0}, Saved in: General.OptionsFileName
Mesh.Color.One = {255,160,0};                 // Color 1 in color carousel, Default value: {255,160,0}, Saved in: General.OptionsFileName
Mesh.Color.Two = {255,200,0};                 // Color 2 in color carousel, Default value: {255,200,0}, Saved in: General.OptionsFileName
Mesh.Color.Three = {255,240,0};               // Color 3 in color carousel, Default value: {255,240,0}, Saved in: General.OptionsFileName
Mesh.Color.Four = {228,255,0};                // Color 4 in color carousel, Default value: {228,255,0}, Saved in: General.OptionsFileName
Mesh.Color.Five = {188,255,0};                // Color 5 in color carousel, Default value: {188,255,0}, Saved in: General.OptionsFileName
Mesh.Color.Six = {148,255,0};                 // Color 6 in color carousel, Default value: {148,255,0}, Saved in: General.OptionsFileName
Mesh.Color.Seven = {148,255,0};               // Color 7 in color carousel, Default value: {108,255,0}, Saved in: General.OptionsFileName
Mesh.Color.Eight = {68,255,0};                // Color 8 in color carousel, Default value: {68,255,0}, Saved in: General.OptionsFileName
Mesh.Color.Nine = {0,255,52};                 // Color 9 in color carousel, Default value: {0,255,52}, Saved in: General.OptionsFileName
Mesh.Color.Ten = {0,255,132};                 // Color 10 in color carousel, Default value: {0,255,132},  Saved in: General.OptionsFileName
Mesh.Color.Eleven = {0,255,192};              // Color 11 in color carousel, Default value: {0,255,192}, Saved in: General.OptionsFileName
Mesh.Color.Twelve = {0,216,255};              // Color 12 in color carousel, Default value: {0,216,255}, Saved in: General.OptionsFileName
Mesh.Color.Thirteen = {0,176,255};            // Color 13 in color carousel, Default value: {0,176,255}, Saved in: General.OptionsFileName
Mesh.Color.Fourteen = {0,116,255};            // Color 14 in color carousel, Default value: {0,116,255}, Saved in: General.OptionsFileName
Mesh.Color.Fifteen = {0,76,255};              // Color 15 in color carousel, Default value: {0,76,255}, Saved in: General.OptionsFileName
Mesh.Color.Sixteen = {24,0,255};              // Color 16 in color carousel, Default value: {24,0,255}, Saved in: General.OptionsFileName
Mesh.Color.Seventeen = {84,0,255};            // Color 17 in color carousel, Default value: {84,0,255}, Saved in: General.OptionsFileName
Mesh.Color.Eighteen = {104,0,255};            // Color 18 in color carousel, Default value: {104,0,255}, Saved in: General.OptionsFileName  
Mesh.Color.Nineteen = {184,0,255};            // Color 19 in color carousel, Default value: {184,0,255}, Saved in: General.OptionsFileName