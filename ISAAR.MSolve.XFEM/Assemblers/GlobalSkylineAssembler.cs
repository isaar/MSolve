using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Assemblers
{
    class GlobalSkylineAssembler
    {
        public (SkylineMatrix Kuu, CSRMatrix Kuc) BuildGlobalMatrix(Model2D model, IDofOrderer dofOrderer)
        {
            int numDofsConstrained = dofOrderer.ConstrainedDofsCount;
            int numDofsUnconstrained = dofOrderer.FreeDofsCount + dofOrderer.EnrichedDofsCount;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            SkylineBuilder Kuu = InitializeGlobalUncontrainedMatrix(model, dofOrderer);

            // TODO: perhaps I should return a CSC matrix and do the transposed multiplication. This way I will not have to 
            // transpose the element matrix. Another approach is to add an AddTransposed() method to the DOK.
            var Kuc = DOKRowMajor.CreateEmpty(numDofsUnconstrained, numDofsConstrained);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Build standard element matrices and add it contributions to the global matrices
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                dofOrderer.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapFree, out IReadOnlyDictionary<int, int> mapConstrained);
                Matrix kss = element.BuildStandardStiffnessMatrix();
                Kuu.AddSubmatrixSymmetric(kss, mapFree);
                Kuc.AddSubmatrix(kss, mapFree, mapConstrained);

                // Build enriched element matrices and add it contributions to the global matrices
                IReadOnlyDictionary<int, int> mapEnriched = dofOrderer.MatchElementToGlobalEnrichedDofsOf(element);
                if (mapEnriched.Count > 0)
                {
                    element.BuildEnrichedStiffnessMatrices(out Matrix kes, out Matrix kee);

                    // TODO: options: 1) Only work with upper triangle in all symmetric matrices. Same applies to Elements.
                    // 2) The Elements have two versions of BuildStiffness(). 
                    // 3) The Elements return both (redundant; If someone needs it he can make it himself like here) 
                    Matrix kse = kes.Transpose();
                    Kuu.AddSubmatrix(kse, mapFree, mapEnriched);
                    Kuu.AddSubmatrixSymmetric(kee, mapEnriched);
                    Kuc.AddSubmatrix(kes, mapEnriched, mapConstrained);
                }                
            }

            #region DEBUG code
            //(Matrix expectedKuu, Matrix expectedKuc) = DenseGlobalAssembler.BuildGlobalMatrix(model, dofOrderer);
            //Console.WriteLine("Check Kuu:");
            //CheckMatrix(expectedKuu, Kuu);
            //Console.WriteLine("Check Kuc:");
            //CheckMatrix(expectedKuc, Kuc);
            #endregion

            //TODO: perhaps I should filter the matrices in the concrete class before returning (e.g. sok.Build())
            return (Kuu.BuildSkylineMatrix(), Kuc.BuildCSRMatrix(true));
        }


        private static SkylineBuilder InitializeGlobalUncontrainedMatrix(Model2D model, IDofOrderer dofOrderer)
        {
            int numDofsUnconstrained = dofOrderer.FreeDofsCount + dofOrderer.EnrichedDofsCount;
            int[] colHeights = new int[numDofsUnconstrained]; //only entries above the diagonal count towards the column height
            foreach (XContinuumElement2D element in model.Elements)
            {
                // To determine the col height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0 entries in the element stiffness matrix.
                int minDof = Int32.MaxValue;
                foreach (XNode2D node in element.Nodes)
                {
                    foreach (int dof in dofOrderer.GetFreeDofsOf(node)) minDof = Math.Min(dof, minDof);
                    foreach (int dof in dofOrderer.GetEnrichedDofsOf(node)) Math.Min(dof, minDof);
                }

                // The height of each col is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (XNode2D node in element.Nodes)
                {
                    foreach (int dof in dofOrderer.GetFreeDofsOf(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDof);
                    }
                    foreach (int dof in dofOrderer.GetEnrichedDofsOf(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDof);
                    }
                }
            }
            return SkylineBuilder.Create(numDofsUnconstrained, colHeights);
        }
    }
}
