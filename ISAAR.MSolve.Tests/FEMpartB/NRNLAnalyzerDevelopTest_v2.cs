//using ISAAR.MSolve.PreProcessor.Elements;
//using ISAAR.MSolve.PreProcessor.Materials;
using System.Collections.Generic;
// compa
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Analyzers;
//using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Discretization.Providers;
using ISAAR.MSolve.FEM;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;
using Xunit;
using ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses;
using ISAAR.MSolve.LinearAlgebra.Commons;

namespace ISAAR.MSolve.Tests.FEMpartB
{
    public static  class NRNLAnalyzerDevelopTest_v2 //Origin: programElegxoiDdm opou eixan ginei comment out kai den htan updated apo ekdosh feat/prosthiki_allagwn 
    {        
        public static (IVector, IVector) SolveDisplLoadsExample()
        {
            #region dhmiourgia montelou
            //VectorExtensions.AssignTotalAffinityCount();
            Model_v2 model = new Model_v2();
            model.SubdomainsDictionary.Add(1, new Subdomain_v2(1));
            // EPILOGH MONTELOU
            int model__builder_choice;
            model__builder_choice = 1;
            if (model__builder_choice == 1) // 
            { HexaCantileverBuilderDispControl(model, 850); }                        
            model.ConnectDataStructures();
            #endregion

            ElementStructuralStiffnessProvider_v2 elementProvider = new ElementStructuralStiffnessProvider_v2();
            Dictionary<int, EquivalentContributionsAssebler_v2> equivalentContributionsAssemblers = new Dictionary<int, EquivalentContributionsAssebler_v2>();//SUNOLIKA STOIXEIA model.SubdomainsDictionary.Count oi oles tis model.subdomains ekei mallon deginontai access me ID.
            equivalentContributionsAssemblers.Add(model.SubdomainsDictionary[1].ID, new EquivalentContributionsAssebler_v2(model.SubdomainsDictionary[1], elementProvider));
            var solverBuilder = new SkylineSolver.Builder();
            solverBuilder.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            var solver = solverBuilder.BuildSolver(model);
            solver.OrderDofs(false);
            foreach (ILinearSystem_v2 linearSystem in solver.LinearSystems.Values)
            {
                linearSystem.Reset(); 
                linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
            }
            //kai dhmiourgia twn mhdenikwn dianusmatwn forces ligo parakatw

            //DdmCalculationsGeneral_v2.BuildModelInterconnectionData(model);
            //var ordering1=solver.DofOrderer.OrderDofs(model);
            //DdmCalculationsGeneral_v2.UndoModelInterconnectionDataBuild(model);

            #region create boundary nodes and create displacements for 1st increment
            Dictionary<int, IVector> uInitialFreeDOFDisplacementsPerSubdomain = new Dictionary<int, IVector>();
            uInitialFreeDOFDisplacementsPerSubdomain.Add(model.SubdomainsDictionary[1].ID, Vector.CreateZero(44));//ordering1.NumGlobalFreeDofs prosoxh sto Id twn subdomain
            Dictionary<int, Node_v2> boundaryNodes = new Dictionary<int, Node_v2>();
            for (int k = 17; k < 21; k++)
            {
                boundaryNodes.Add(model.NodesDictionary[k].ID, model.NodesDictionary[k]);
            }
            Dictionary<int, Dictionary<DOFType, double>> initialConvergedBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            Dictionary<DOFType, double> initialConvergedBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
            initialConvergedBoundaryNodalDisplacements.Add(DOFType.X, 0);
            for (int k = 17; k < 21; k++)
            {
                initialConvergedBoundaryDisplacements.Add(model.NodesDictionary[k].ID, initialConvergedBoundaryNodalDisplacements);
            }
            Dictionary<int, Dictionary<DOFType, double>> totalBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            double[] prescribedDisplacmentXValues = new double[4] { 7.81614E-01, 7.07355E-01, 7.81614E-01, 7.07355E-01 };
            for (int k = 17; k < 21; k++)
            {
                Dictionary<DOFType, double> totalBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
                totalBoundaryNodalDisplacements.Add(DOFType.X, 0.5*prescribedDisplacmentXValues[k - 17]);
                totalBoundaryDisplacements.Add(model.NodesDictionary[k].ID, totalBoundaryNodalDisplacements);
            }
            #endregion

            #region create nesessary structures and analyzers And Solve 1st increment            
            var linearSystems = solver.LinearSystems; // elegxos me model.subdomainsDictionary[1]
            ProblemStructural_v2 provider = new ProblemStructural_v2(model, solver);
            var subdomainUpdaters = new Dictionary<int, NonLinearSubdomainUpdaterWithInitialConditions_v2>(1);
            subdomainUpdaters.Add( 1 ,new NonLinearSubdomainUpdaterWithInitialConditions_v2(model.Subdomains[0]));
            var increments = 1;

            var childAnalyzer = new MicrostructureBvpNRNLAnalyzer(model,solver, subdomainUpdaters,  provider, increments, uInitialFreeDOFDisplacementsPerSubdomain,
                boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, equivalentContributionsAssemblers);            
            childAnalyzer.SetMaxIterations = 100;
            childAnalyzer.SetIterationsForMatrixRebuild = 1;
            
            MSParentAnalyzer parentAnalyzer = new MSParentAnalyzer(model, solver, provider, childAnalyzer);
            //TODO MS
            //foreach (ILinearSystem_v2 linearSystem in solver.LinearSystems.Values)
            //{
            //    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size); // antistoixo tou subdomain.Forces = linearSystem.CreateZeroVector();
            //}
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            IVector uInitialFreeDOFs_state1 = childAnalyzer.GetConvergedSolutionVectorsOfFreeDofs()[1].Copy();
            #endregion

            #region save state and update structures and vectors for second increment
            foreach (var subdomainUpdater in subdomainUpdaters.Values)
            {
                subdomainUpdater.UpdateState();
            }
            // u (or uplusDu) initial 
            uInitialFreeDOFDisplacementsPerSubdomain = childAnalyzer.GetConvergedSolutionVectorsOfFreeDofs();// ousiastika to u pou twra taftizetai me to uPlusuu

            initialConvergedBoundaryDisplacements = totalBoundaryDisplacements;

            totalBoundaryDisplacements = new Dictionary<int, Dictionary<DOFType, double>>();
            for (int k = 17; k < 21; k++)
            {
                Dictionary<DOFType, double> totalBoundaryNodalDisplacements = new Dictionary<DOFType, double>();
                totalBoundaryNodalDisplacements.Add(DOFType.X, 1.0 * prescribedDisplacmentXValues[k - 17]);
                totalBoundaryDisplacements.Add(model.NodesDictionary[k].ID, totalBoundaryNodalDisplacements);
            }
            #endregion

            #region Creation of nessesary analyzers and solution 
            ElementStructuralStiffnessProvider_v2 elementProvider2 = new ElementStructuralStiffnessProvider_v2();
            Dictionary<int, EquivalentContributionsAssebler_v2> equivalentContributionsAssemblers2 = new Dictionary<int, EquivalentContributionsAssebler_v2>();//SUNOLIKA STOIXEIA model.SubdomainsDictionary.Count oi oles tis model.subdomains ekei mallon deginontai access me ID.
            equivalentContributionsAssemblers2.Add(model.SubdomainsDictionary[1].ID, new EquivalentContributionsAssebler_v2(model.SubdomainsDictionary[1], elementProvider2));
            var solverBuilder2 = new SkylineSolver.Builder();
            solverBuilder2.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            ISolver_v2 solver2 = solverBuilder2.BuildSolver(model);
            solver2.OrderDofs(false);
            foreach (ILinearSystem_v2 linearSystem in solver2.LinearSystems.Values) linearSystem.Reset();
            //kalutera apotelesmata otan to parakatw den kratietai, ara pragmati resets ta subd.forces
            //solver.ResetSubdomainForcesVector();
            var linearSystems2 = solver2.LinearSystems; // elegxos me model.subdomainsDictionary[1]
            foreach (ILinearSystem_v2 linearSystem in linearSystems2.Values)
            {                
                linearSystem.RhsVector = linearSystem.Subdomain.Forces; //TODO MS 
            }
            ProblemStructural_v2 provider2 = new ProblemStructural_v2(model, solver2);
            //var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters2 = new Dictionary<int, NonLinearSubdomainUpdaterWithInitialConditions_v2>(1);
            subdomainUpdaters2.Add(1, new NonLinearSubdomainUpdaterWithInitialConditions_v2(model.Subdomains[0]));
            //var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };
            var increments2 = 1;

            var childAnalyzer2 = new MicrostructureBvpNRNLAnalyzer(model, solver2, subdomainUpdaters2, provider2, increments2,  uInitialFreeDOFDisplacementsPerSubdomain,
                boundaryNodes, initialConvergedBoundaryDisplacements, totalBoundaryDisplacements, equivalentContributionsAssemblers2);
            childAnalyzer2.SetMaxIterations = 100;
            childAnalyzer2.SetIterationsForMatrixRebuild = 1;

            MSParentAnalyzer parentAnalyzer2 = new MSParentAnalyzer(model, solver2, provider2, childAnalyzer2);
            //parentAnalyzer2.BuildMatrices();
            //DdmCalculationsGeneral_v2.UndoModelInterconnectionDataBuild(model);
            childAnalyzer2.Initialize(); //parentAnalyzer2.Initialize();
            parentAnalyzer2.Solve();
            IVector uInitialFreeDOFs_state2 = childAnalyzer.GetConvergedSolutionVectorsOfFreeDofs()[1].Copy();
            #endregion

            return (uInitialFreeDOFs_state1, uInitialFreeDOFs_state2);
        }

        public static bool HexaCantileverBuilderDispControl(Model_v2 model, double load_value)
        {
            ElasticMaterial3D_v2 material1 = new ElasticMaterial3D_v2()
            {
                YoungModulus = 1353000,
                PoissonRatio = 0.3,
            };

            double[,] nodeData = new double[,] { {-0.250000,-0.250000,-1.000000},
            {0.250000,-0.250000,-1.000000},
            {-0.250000,0.250000,-1.000000},
            {0.250000,0.250000,-1.000000},
            {-0.250000,-0.250000,-0.500000},
            {0.250000,-0.250000,-0.500000},
            {-0.250000,0.250000,-0.500000},
            {0.250000,0.250000,-0.500000},
            {-0.250000,-0.250000,0.000000},
            {0.250000,-0.250000,0.000000},
            {-0.250000,0.250000,0.000000},
            {0.250000,0.250000,0.000000},
            {-0.250000,-0.250000,0.500000},
            {0.250000,-0.250000,0.500000},
            {-0.250000,0.250000,0.500000},
            {0.250000,0.250000,0.500000},
            {-0.250000,-0.250000,1.000000},
            {0.250000,-0.250000,1.000000},
            {-0.250000,0.250000,1.000000},
            {0.250000,0.250000,1.000000}};

            int[,] elementData = new int[,] {{1,8,7,5,6,4,3,1,2},
            {2,12,11,9,10,8,7,5,6},
            {3,16,15,13,14,12,11,9,10},
            {4,20,19,17,18,16,15,13,14}, };

            // orismos shmeiwn
            for (int nNode = 0; nNode < nodeData.GetLength(0); nNode++)
            {
                model.NodesDictionary.Add(nNode + 1, new Node_v2() { ID = nNode + 1, X = nodeData[nNode, 0], Y = nodeData[nNode, 1], Z = nodeData[nNode, 2] });

            }

            // orismos elements 
            Element_v2 e1;
            int subdomainID = 1;
            for (int nElement = 0; nElement < elementData.GetLength(0); nElement++)
            {
                e1 = new Element_v2()
                {
                    ID = nElement + 1,
                    ElementType = new Hexa8NonLinear_v2(material1, GaussLegendre3D.GetQuadratureWithOrder(3, 3, 3)) // dixws to e. exoume sfalma enw sto beambuilding oxi//edw kaleitai me ena orisma to Hexa8
                };
                for (int j = 0; j < 8; j++)
                {
                    e1.NodesDictionary.Add(elementData[nElement, j + 1], model.NodesDictionary[elementData[nElement, j + 1]]);
                }
                model.ElementsDictionary.Add(e1.ID, e1);
                model.SubdomainsDictionary[subdomainID].Elements.Add(e1);
            }

            // constraint vashh opou z=-1
            for (int k = 1; k < 5; k++)
            {
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.X });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.Y });
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.Z });
            }

            // thetoume constraint tous prescribed
            //Load load1;
            for (int k = 17; k < 21; k++)
            {
                //load1 = new Load()
                //{
                //    Node = model.NodesDictionary[k],
                //    DOF = DOFType.X,
                //    Amount = 1 * load_value
                //};
                //model.Loads.Add(load1);
                model.NodesDictionary[k].Constraints.Add(new Constraint { DOF = DOFType.X });
            }

            return true;
        }

        public static bool AreDisplacementsSame_v2(double[] expectedValues,
            double[] computedValues)
        {
            var comparer = new ValueComparer(1E-14);
            for (int i1 = 0; i1 < expectedValues.GetLength(0); i1++)
            {                
                    if (!comparer.AreEqual(expectedValues[i1], computedValues[i1]))
                    {
                        return false;
                    }               
            }
            return true;
        }
    }
}
