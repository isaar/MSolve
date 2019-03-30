using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.Feti1;
using ISAAR.MSolve.Solvers.DomainDecomposition.Dual.LagrangeMultipliers;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Dual
{
    public static class FullyRedundantLagrangeMultipliersTests
    {
        private const int numLagranges = 20;

        [Fact]
        public static void TestLagrangeMultipliersCreation()
        {
            Model_v2 model = CreateModel();
            Dictionary<int, Node_v2> nodes = model.NodesDictionary;
            Dictionary<int, Subdomain_v2> subdomains = model.SubdomainsDictionary;
            LagrangeMultiplier[] lagrangesComputed = CreateBooleanMultipliers(model, true).LagrangeMultipliers;
            var lagrangesExpected = new LagrangeMultiplier[numLagranges];

            // (4, 5)  (6, 7)    (4,5)   (6,7)
            // 6 ----- 7         7 ----- 8 
            // |  s2   |         |   s3  | 
            // |       |         |       | 
            // 3 ----- 4         4 ----- 5
            // (0, 1)  (2, 3)    (0,1)   (2,3)
            //
            //
            // (2,3)   (4,5)     (3,4)   (5,6)
            // 3 ----- 4         4 ----- 5 
            // |       |         |       | 
            // |  s0   |         |   s1  | 
            // 0 ----- 1         1 ----- 2
            //         (0,1)     (0,1)   (2)

            lagrangesExpected[0] = new LagrangeMultiplier(nodes[1], DOFType.X, subdomains[0], subdomains[1]);
            lagrangesExpected[1] = new LagrangeMultiplier(nodes[1], DOFType.Y, subdomains[0], subdomains[1]);

            lagrangesExpected[2] = new LagrangeMultiplier(nodes[3], DOFType.X, subdomains[0], subdomains[2]);
            lagrangesExpected[3] = new LagrangeMultiplier(nodes[3], DOFType.Y, subdomains[0], subdomains[2]);

            lagrangesExpected[4] = new LagrangeMultiplier(nodes[4], DOFType.X, subdomains[0], subdomains[1]);
            lagrangesExpected[5] = new LagrangeMultiplier(nodes[4], DOFType.Y, subdomains[0], subdomains[1]);
            lagrangesExpected[6] = new LagrangeMultiplier(nodes[4], DOFType.X, subdomains[0], subdomains[2]);
            lagrangesExpected[7] = new LagrangeMultiplier(nodes[4], DOFType.Y, subdomains[0], subdomains[2]);
            lagrangesExpected[8] = new LagrangeMultiplier(nodes[4], DOFType.X, subdomains[0], subdomains[3]);
            lagrangesExpected[9] = new LagrangeMultiplier(nodes[4], DOFType.Y, subdomains[0], subdomains[3]);

            lagrangesExpected[10] = new LagrangeMultiplier(nodes[4], DOFType.X, subdomains[1], subdomains[2]);
            lagrangesExpected[11] = new LagrangeMultiplier(nodes[4], DOFType.Y, subdomains[1], subdomains[2]);
            lagrangesExpected[12] = new LagrangeMultiplier(nodes[4], DOFType.X, subdomains[1], subdomains[3]);
            lagrangesExpected[13] = new LagrangeMultiplier(nodes[4], DOFType.Y, subdomains[1], subdomains[3]);

            lagrangesExpected[14] = new LagrangeMultiplier(nodes[4], DOFType.X, subdomains[2], subdomains[3]);
            lagrangesExpected[15] = new LagrangeMultiplier(nodes[4], DOFType.Y, subdomains[2], subdomains[3]);

            lagrangesExpected[16] = new LagrangeMultiplier(nodes[5], DOFType.X, subdomains[1], subdomains[3]);
            lagrangesExpected[17] = new LagrangeMultiplier(nodes[5], DOFType.Y, subdomains[1], subdomains[3]);

            lagrangesExpected[18] = new LagrangeMultiplier(nodes[7], DOFType.X, subdomains[2], subdomains[3]);
            lagrangesExpected[19] = new LagrangeMultiplier(nodes[7], DOFType.Y, subdomains[2], subdomains[3]);

            for (int lag = 0; lag < numLagranges; ++lag)
            {
                Assert.True(lagrangesComputed[lag].Node == lagrangesExpected[lag].Node);
                Assert.True(lagrangesComputed[lag].DofType == lagrangesExpected[lag].DofType);
                Assert.True(lagrangesComputed[lag].SubdomainPlus == lagrangesExpected[lag].SubdomainPlus);
                Assert.True(lagrangesComputed[lag].SubdomainMinus == lagrangesExpected[lag].SubdomainMinus);
            }
        }

        [Fact]
        public static void TestSignedBooleanMatrices()
        {   
            // Node | Lagrange | Sum
            // 1      2          20
            // 3      2           
            // 4      12          
            // 5      2           
            // 7      2           
            // Dof notation: i, b, c = internal, boundary, corner dof
            // Numbering order: node major, subdomain medium, dof minor
            var booleansExpected = new Dictionary<int, Matrix>();

            // Subdomain 0: 6 free dofs
            // (2,3)   (4,5)
            // 3 ----- 4 
            // |       | 
            // |       | 
            // 0 ----- 1
            //         (0,1)
            booleansExpected[0] = Matrix.CreateFromArray(new double[numLagranges, 6]
            {
                // b   b   b   b   c   c
                {  1,  0,  0,  0,  0,  0 }, // Node 1
                {  0,  1,  0,  0,  0,  0 }, // Node 1
                {  0,  0,  1,  0,  0,  0 }, // Node 3
                {  0,  0,  0,  1,  0,  0 }, // Node 3
                {  0,  0,  0,  0,  1,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  1 }, // Node 4
                {  0,  0,  0,  0,  1,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  1 }, // Node 4
                {  0,  0,  0,  0,  1,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  1 }, // Node 4
                {  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0 }, // Node 5
                {  0,  0,  0,  0,  0,  0 }, // Node 5
                {  0,  0,  0,  0,  0,  0 }, // Node 7
                {  0,  0,  0,  0,  0,  0 }  // Node 7
            });

            // Subdomain 1: 7 free dofs
            // (3,4)   (5,6)
            // 4 ----- 5 
            // |       | 
            // |       | 
            // 1 ----- 2
            // (0,1)   (2)
            booleansExpected[1] = Matrix.CreateFromArray(new double[numLagranges, 7]
            {
                // b   b   i   c   c   b   b
                { -1,  0,  0,  0,  0,  0,  0 }, // Node 1
                {  0, -1,  0,  0,  0,  0,  0 }, // Node 1
                {  0,  0,  0,  0,  0,  0,  0 }, // Node 3
                {  0,  0,  0,  0,  0,  0,  0 }, // Node 3
                {  0,  0,  0, -1,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0, -1,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  1,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  1,  0,  0 }, // Node 4
                {  0,  0,  0,  1,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  1,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  1,  0 }, // Node 5
                {  0,  0,  0,  0,  0,  0,  1 }, // Node 5
                {  0,  0,  0,  0,  0,  0,  0 }, // Node 7
                {  0,  0,  0,  0,  0,  0,  0 }  // Node 7
            });

            // Subdomain 2: 8 free dofs
            // (4,5)   (6,7)
            // 6 ----- 7 
            // |       | 
            // |       | 
            // 3 ----- 4
            // (0,1)   (2,3)
            booleansExpected[2] = Matrix.CreateFromArray(new double[numLagranges, 8]
            {
                // b   b   c   c   i   i   b   b
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 1
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 1
                { -1,  0,  0,  0,  0,  0,  0,  0 }, // Node 3
                {  0, -1,  0,  0,  0,  0,  0,  0 }, // Node 3
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0, -1,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0, -1,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0, -1,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0, -1,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  1,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  1,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 5
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 5
                {  0,  0,  0,  0,  0,  0,  1,  0 }, // Node 7
                {  0,  0,  0,  0,  0,  0,  0,  1 }  // Node 7
            });

            // Subdomain 3: 8 free dofs
            // (4,5)   (6,7)
            // 7 ----- 8 
            // |       | 
            // |       | 
            // 4 ----- 5
            // (0,1)   (2,3)
            booleansExpected[3] = Matrix.CreateFromArray(new double[numLagranges, 8]
            {
                // c   c   b   b   b   b   i   i
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 1
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 1
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 3
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 3
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                { -1,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0, -1,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                { -1,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0, -1,  0,  0,  0,  0,  0,  0 }, // Node 4
                { -1,  0,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0, -1,  0,  0,  0,  0,  0,  0 }, // Node 4
                {  0,  0, -1,  0,  0,  0,  0,  0 }, // Node 5
                {  0,  0,  0, -1,  0,  0,  0,  0 }, // Node 5
                {  0,  0,  0,  0, -1,  0,  0,  0 }, // Node 7
                {  0,  0,  0,  0,  0, -1,  0,  0 }  // Node 7
            });

            Model_v2 model = CreateModel();
            LagrangeMultipliersEnumerator lagrangeMultipliers = CreateBooleanMultipliers(model, false);
            foreach (var id in lagrangeMultipliers.BooleanMatrices.Keys)
            {
                Matrix booleanComputed = lagrangeMultipliers.BooleanMatrices[id].CopyToFullMatrix(false);
                Assert.True(booleansExpected[id].Equals(booleanComputed));
            }
        }

        private static LagrangeMultipliersEnumerator CreateBooleanMultipliers(Model_v2 model, bool createLagranges)
        {
            // Initialize model
            model.ConnectDataStructures();

            // Order freedom degrees
            var orderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
            IGlobalFreeDofOrdering globalOrdering = orderer.OrderFreeDofs(model);
            model.GlobalDofOrdering = globalOrdering;
            foreach (Subdomain_v2 subdomain in model.Subdomains)
            {
                subdomain.FreeDofOrdering = globalOrdering.SubdomainDofOrderings[subdomain];
            }

            // Create boolean matrices
            var dofSeparator = new Feti1DofSeparator();
            dofSeparator.SeparateBoundaryInternalDofs(model);
            var lagrangeEnumerator = new LagrangeMultipliersEnumerator(new FullyRedundantConstraints());
            if (createLagranges) lagrangeEnumerator.DefineLagrangesAndBooleanMatrices(model, dofSeparator);
            else lagrangeEnumerator.DefineBooleanMatrices(model, dofSeparator);
            return lagrangeEnumerator;
        }

        private static Model_v2 CreateModel()
        {
            // 6 ----- 7 ----- 8
            // |  (2)  |  (3)  |
            // |       |       |
            // 3 ----- 4 ----- 5
            // |  (0)  |  (1)  |
            // |       |       |
            // 0 ----- 1 ----- 2
            // Δ               Δ    
            // -               o

            var builder = new Uniform2DModelBuilder();
            builder.DomainLengthX = 2.0;
            builder.DomainLengthY = 2.0;
            builder.NumSubdomainsX = 2;
            builder.NumSubdomainsY = 2;
            builder.NumTotalElementsX = 2;
            builder.NumTotalElementsY = 2;
            builder.YoungModulus = 2.1E7;
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerLeftCorner, DOFType.X, 0.0);            
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerLeftCorner, DOFType.Y, 0.0);
            builder.PrescribeDisplacement(Uniform2DModelBuilder.BoundaryRegion.LowerRightCorner, DOFType.Y, 0.0);

            return builder.BuildModel();
        }
    }
}
