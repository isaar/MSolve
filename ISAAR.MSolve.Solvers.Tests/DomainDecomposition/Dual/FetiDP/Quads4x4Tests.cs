using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP
{
    public static class Quads4x4Tests
    {
        [Fact]
        public static void TestDofSeparation()
        {
            int[][] cornerDofsExpected = new int[4][];
            int[][] remainderDofsExpected = new int[4][];
            int[][] boundaryRemainderDofsExpected = new int[4][];
            int[][] internalRemainderDofsExpected = new int[4][];

            // Subdomain 0
            cornerDofsExpected[0] = new int[] { 2, 3, 10, 11 };
            remainderDofsExpected[0] = new int[] { 0, 1, 4, 5, 6, 7, 8, 9 };
            boundaryRemainderDofsExpected[0] = new int[] { 4, 5, 6, 7 };
            internalRemainderDofsExpected[0] = new int[] { 0, 1, 2, 3 };

            // Subdomain 1
            cornerDofsExpected[1] = new int[] { 0, 1, 12, 13, 16, 17 };
            remainderDofsExpected[1] = new int[] { 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15 };
            boundaryRemainderDofsExpected[1] = new int[] { 4, 5, 10, 11 };
            internalRemainderDofsExpected[1] = new int[] { 0, 1, 2, 3, 6, 7, 8, 9 };

            // Subdomain 2
            cornerDofsExpected[2] = new int[] { 2, 3, 10, 11 };
            remainderDofsExpected[2] = new int[] { 0, 1, 4, 5, 6, 7, 8, 9 };
            boundaryRemainderDofsExpected[2] = new int[] { 0, 1, 4, 5 };
            internalRemainderDofsExpected[2] = new int[] { 2, 3, 6, 7 };

            // Subdomain 3
            cornerDofsExpected[3] = new int[] { 0, 1, 4, 5, 12, 13 };
            remainderDofsExpected[3] = new int[] { 2, 3, 6, 7, 8, 9, 10, 11, 14, 15, 16, 17 };
            boundaryRemainderDofsExpected[3] = new int[] { 0, 1, 2, 3 };
            internalRemainderDofsExpected[3] = new int[] { 4, 5, 6, 7, 8, 9, 10, 11 };

            // Create model
            Model model = CreateModel();
            Dictionary<int, INode[]> cornerNodes = DefineCornerNodes(model);
            model.ConnectDataStructures();

            // Order free dofs.
            var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            IGlobalFreeDofOrdering globalOrdering = dofOrderer.OrderFreeDofs(model);
            model.GlobalDofOrdering = globalOrdering;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                subdomain.FreeDofOrdering = globalOrdering.SubdomainDofOrderings[subdomain];
            }

            //SkylineSolver solver = new SkylineSolver.Builder().BuildSolver(model); //TODO: It should not be necessary to call another solver.
            //solver.OrderDofs(false);

            var dofSeparator = new FetiDPDofSeparator();
            dofSeparator.SeparateDofs(model, cornerNodes);

            // Check
            for (int s = 0; s < 4; ++s)
            {
                Utilities.CheckEqual(cornerDofsExpected[s], dofSeparator.CornerIntoFreeDofIndices[s]);
                Utilities.CheckEqual(remainderDofsExpected[s], dofSeparator.RemainderIntoFreeDofIndices[s]);
                Utilities.CheckEqual(boundaryRemainderDofsExpected[s], dofSeparator.BoundaryIntoRemainderDofIndices[s]);
                Utilities.CheckEqual(internalRemainderDofsExpected[s], dofSeparator.InternalIntoRemainderDofIndices[s]);
            }
        }

        private static Model CreateModel()
        {
            //                                    Λ P
            //                                    | 
            //                                     
            // |> 20 ---- 21 ---- 22 ---- 23 ---- 24
            //    |  (12) |  (13) |  (14) |  (15) |
            //    |       |       |       |       |
            // |> 15 ---- 16 ---- 17 ---- 18 ---- 19
            //    |  (8)  |  (9)  |  (10) |  (11) |
            //    |       |       |       |       |
            // |> 10 ---- 11 ---- 12 ---- 13 ---- 14
            //    |  (4)  |  (5)  |  (6)  |  (7)  |
            //    |       |       |       |       |
            // |> 5 ----- 6 ----- 7 ----- 8 ----- 9
            //    |  (0)  |  (1)  |  (2)  |  (3)  |
            //    |       |       |       |       |
            // |> 0 ----- 1 ----- 2 ----- 3 ----- 4


            var builder = new Uniform2DModelBuilder();
            builder.DomainLengthX = 4.0;
            builder.DomainLengthY = 4.0;
            builder.NumSubdomainsX = 2;
            builder.NumSubdomainsY = 2;
            builder.NumTotalElementsX = 4;
            builder.NumTotalElementsY = 4;
            builder.YoungModulus = 1.0;
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationX, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LeftSide, StructuralDof.TranslationY, 0.0);
            builder.DistributeLoadAtNodes(Uniform2DModelBuilder.BoundaryRegion.UpperRightCorner, StructuralDof.TranslationY, 10.0);

            return builder.BuildModel();
        }

        private static Dictionary<int, INode[]> DefineCornerNodes(Model model)
        {
            // subdomain 2         subdomain 3                      
            // 20 ---- 21 ---- 22  22---- 23 ---- 24
            // |  (12) |  (13) |   | (14) |  (15) |
            // |       |       |   |      |       |
            // 15 ---- 16 ---- 17  17---- 18 ---- 19
            // |  (8)  |  (9)  |   | (10) |  (11) |
            // |       |       |   |      |       |
            // 10 ---- 11 ---- 12  12---- 13 ---- 14

            // subdomain 0         subdomain 1
            // 10 ---- 11 ---- 12  12---- 13 ---- 14
            // |  (4)  |  (5)  |   | (6)  |  (7)  |
            // |       |       |   |      |       |
            // 5 ----- 6 ----- 7   7 ---- 8 ----- 9
            // |  (0)  |  (1)  |   | (2)  |  (3)  |
            // |       |       |   |      |       |
            // 0 ----- 1 ----- 2   2 ---- 3 ----- 4

            var cornerNodes = new Dictionary<int, INode[]>();
            cornerNodes[0] = new INode[] { model.Nodes[2], model.Nodes[12] };
            cornerNodes[1] = new INode[] { model.Nodes[2], model.Nodes[12], model.Nodes[14] };
            cornerNodes[2] = new INode[] { model.Nodes[12], model.Nodes[22] };
            cornerNodes[3] = new INode[] { model.Nodes[12], model.Nodes[14], model.Nodes[22] };
            return cornerNodes;
        }
    }
}
