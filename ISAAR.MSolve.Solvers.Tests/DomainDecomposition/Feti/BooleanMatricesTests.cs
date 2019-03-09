using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Preprocessor.Meshes;
using ISAAR.MSolve.Preprocessor.Meshes.Custom;
using ISAAR.MSolve.Solvers.DomainDecomposition.Feti;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;
using Xunit;

namespace ISAAR.MSolve.Solvers.Tests.DomainDecomposition.Feti
{
    public static class BooleanMatricesTests
    {
        [Fact]
        public static void TestFullyRedundantConstraints2D()
        {   // Node - Continuity equations - Sum
            // 1      2                      20
            // 3      2           
            // 4      12          
            // 5      2           
            // 7      2           
            // Dof notation: i, b, c = internal, boundary, corner dof
            // Numbering order: node major, subdomain medium, dof minor
            var booleansExpected = new Dictionary<int, Matrix>();

            // Subdomain 0: 6 free dofs
            // (2,3)    (4,5)
            // 3 ----- 4 
            // |       | 
            // |       | 
            // 0 ----- 1
            //         (0,1)
            booleansExpected[0] = Matrix.CreateFromArray(new double[20, 6]
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
            booleansExpected[1] = Matrix.CreateFromArray(new double[20, 7]
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
            booleansExpected[2] = Matrix.CreateFromArray(new double[20, 8]
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
            booleansExpected[3] = Matrix.CreateFromArray(new double[20, 8]
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
            Dictionary<int, SignedBooleanMatrix> booleanMatrices = CreateBooleanMatrices(model);
            foreach (var id in booleanMatrices.Keys)
            {
                Matrix booleanComputed = booleanMatrices[id].CopyToFullMatrix(false);
                Assert.True(booleansExpected[id].Equals(booleanComputed));
            }
        }

        private static Dictionary<int, SignedBooleanMatrix> CreateBooleanMatrices(Model_v2 model)
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
            var continuityEquations = new ContinuityEquationsCalculator(new FullyRedundantConstraints());
            continuityEquations.CreateBooleanMatrices(model);
            return continuityEquations.BooleanMatrices;
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
