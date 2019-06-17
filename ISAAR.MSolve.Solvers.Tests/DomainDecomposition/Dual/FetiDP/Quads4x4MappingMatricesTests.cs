using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Operators;
using ISAAR.MSolve.Solvers.Direct;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual.FetiDP
{
    public static class Quads4x4MappingMatricesTests
    {
        [Fact]
        public static void TestDofSeparation()
        {
            //TODO: These should be available from properties
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
            Dictionary<int, HashSet<INode>> cornerNodes = DefineCornerNodes(model);
            model.ConnectDataStructures();

            // Order free dofs.
            var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            IGlobalFreeDofOrdering globalOrdering = dofOrderer.OrderFreeDofs(model);
            model.GlobalDofOrdering = globalOrdering;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                subdomain.FreeDofOrdering = globalOrdering.SubdomainDofOrderings[subdomain];
            }

            // Separate dofs
            var dofSeparator = new FetiDPDofSeparator();
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int s = subdomain.ID;
                IEnumerable<INode> remainderNodes = subdomain.Nodes.Except(cornerNodes[s]);
                dofSeparator.SeparateCornerRemainderDofs(subdomain, cornerNodes[s], remainderNodes);
                dofSeparator.SeparateBoundaryInternalDofs(subdomain, remainderNodes);
            }
            
            // Check
            for (int s = 0; s < 4; ++s)
            {
                Utilities.CheckEqual(cornerDofsExpected[s], dofSeparator.CornerDofIndices[s]);
                Utilities.CheckEqual(remainderDofsExpected[s], dofSeparator.RemainderDofIndices[s]);
                Utilities.CheckEqual(boundaryRemainderDofsExpected[s], dofSeparator.BoundaryDofIndices[s]);
                Utilities.CheckEqual(internalRemainderDofsExpected[s], dofSeparator.InternalDofIndices[s]);
            }
        }

        [Fact]
        public static void TestSignedBooleanMatrices()
        {
            // Expected results
            int expectedNumLagrangeMultipliers = 8;
            var expectedBr = new Dictionary<int, Matrix>();

            expectedBr[0] = Matrix.CreateZero(8, 8);
            expectedBr[0][0, 4] = +1;
            expectedBr[0][1, 5] = +1;
            expectedBr[0][2, 6] = +1;
            expectedBr[0][3, 7] = +1;

            expectedBr[1] = Matrix.CreateZero(8, 12);
            expectedBr[1][0, 4] = -1;
            expectedBr[1][1, 5] = -1;
            expectedBr[1][4, 10] = +1;
            expectedBr[1][5, 11] = +1;

            expectedBr[2] = Matrix.CreateZero(8, 8);
            expectedBr[2][2, 0] = -1;
            expectedBr[2][3, 1] = -1;
            expectedBr[2][6, 4] = +1;
            expectedBr[2][7, 5] = +1;

            expectedBr[3] = Matrix.CreateZero(8, 12);
            expectedBr[3][4, 0] = -1;
            expectedBr[3][5, 1] = -1;
            expectedBr[3][6, 2] = -1;
            expectedBr[3][7, 3] = -1;

            // Create model
            Model model = CreateModel();
            Dictionary<int, HashSet<INode>> cornerNodes = DefineCornerNodes(model);
            model.ConnectDataStructures();

            // Order free dofs.
            var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            IGlobalFreeDofOrdering globalOrdering = dofOrderer.OrderFreeDofs(model);
            model.GlobalDofOrdering = globalOrdering;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                subdomain.FreeDofOrdering = globalOrdering.SubdomainDofOrderings[subdomain];
            }

            // Separate dofs
            var dofSeparator = new FetiDPDofSeparator();
            dofSeparator.DefineGlobalBoundaryDofs(model, cornerNodes);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int s = subdomain.ID;
                IEnumerable<INode> remainderNodes = subdomain.Nodes.Except(cornerNodes[s]);
                dofSeparator.SeparateCornerRemainderDofs(subdomain, cornerNodes[s], remainderNodes);
                dofSeparator.SeparateBoundaryInternalDofs(subdomain, remainderNodes);
            }

            // Enumerate lagranges
            var crosspointStrategy = new FullyRedundantConstraints();
            var lagrangeEnumerator = new FetiDPLagrangeMultipliersEnumerator(crosspointStrategy, dofSeparator);
            lagrangeEnumerator.DefineBooleanMatrices(model);

            // Check
            double tolerance = 1E-13;
            Assert.Equal(expectedNumLagrangeMultipliers, lagrangeEnumerator.NumLagrangeMultipliers);
            for (int id = 0; id < 4; ++id)
            {
                Matrix Br = lagrangeEnumerator.BooleanMatrices[id].CopyToFullMatrix(false);
                Assert.True(expectedBr[id].Equals(Br, tolerance));
            }
        }

        [Fact]
        public static void TestUnsignedBooleanMatrices()
        {
            // Expected results
            int expectedNumCornerDofs = 8;
            var expectedLc = new Dictionary<int, Matrix>();

            expectedLc[0] = Matrix.CreateZero(4, 8);
            expectedLc[0][0, 0] = 1;
            expectedLc[0][1, 1] = 1;
            expectedLc[0][2, 2] = 1;
            expectedLc[0][3, 3] = 1;

            expectedLc[1] = Matrix.CreateZero(6, 8);
            expectedLc[1][0, 0] = 1;
            expectedLc[1][1, 1] = 1;
            expectedLc[1][2, 2] = 1;
            expectedLc[1][3, 3] = 1;
            expectedLc[1][4, 4] = 1;
            expectedLc[1][5, 5] = 1;

            expectedLc[2] = Matrix.CreateZero(4, 8);
            expectedLc[2][0, 2] = 1;
            expectedLc[2][1, 3] = 1;
            expectedLc[2][2, 6] = 1;
            expectedLc[2][3, 7] = 1;

            expectedLc[3] = Matrix.CreateZero(6, 8);
            expectedLc[3][0, 2] = 1;
            expectedLc[3][1, 3] = 1;
            expectedLc[3][2, 4] = 1;
            expectedLc[3][3, 5] = 1;
            expectedLc[3][4, 6] = 1;
            expectedLc[3][5, 7] = 1;

            // Create model
            Model model = CreateModel();
            Dictionary<int, HashSet<INode>> cornerNodes = DefineCornerNodes(model);
            model.ConnectDataStructures();

            // Order free dofs.
            var dofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            IGlobalFreeDofOrdering globalOrdering = dofOrderer.OrderFreeDofs(model);
            model.GlobalDofOrdering = globalOrdering;
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                subdomain.FreeDofOrdering = globalOrdering.SubdomainDofOrderings[subdomain];
            }

            // Separate dofs
            var dofSeparator = new FetiDPDofSeparator();
            dofSeparator.DefineGlobalCornerDofs(model, cornerNodes);
            foreach (ISubdomain subdomain in model.Subdomains)
            {
                int s = subdomain.ID;
                IEnumerable<INode> remainderAndConstrainedNodes = subdomain.Nodes.Where(node => !cornerNodes[s].Contains(node));
                dofSeparator.SeparateCornerRemainderDofs(subdomain, cornerNodes[s], remainderAndConstrainedNodes);
                dofSeparator.SeparateBoundaryInternalDofs(subdomain, remainderAndConstrainedNodes);
            }
            dofSeparator.CalcCornerMappingMatrices(model, cornerNodes);

            // Check
            double tolerance = 1E-13;
            Assert.Equal(expectedNumCornerDofs, dofSeparator.NumGlobalCornerDofs);
            for (int id = 0; id < 4; ++id)
            {
                UnsignedBooleanMatrix Lc = dofSeparator.CornerBooleanMatrices[id];
                Assert.True(expectedLc[id].Equals(Lc, tolerance));
            }
        }

        public static Model CreateModel()
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

        public static Dictionary<int, HashSet<INode>> DefineCornerNodes(Model model)
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

            var cornerNodes = new Dictionary<int, HashSet<INode>>();
            cornerNodes[0] = new HashSet<INode>(new INode[] { model.Nodes[2], model.Nodes[12] });
            cornerNodes[1] = new HashSet<INode>(new INode[] { model.Nodes[2], model.Nodes[12], model.Nodes[14] });
            cornerNodes[2] = new HashSet<INode>(new INode[] { model.Nodes[12], model.Nodes[22] });
            cornerNodes[3] = new HashSet<INode>(new INode[] { model.Nodes[12], model.Nodes[14], model.Nodes[22] });
            return cornerNodes;
        }
    }
}
