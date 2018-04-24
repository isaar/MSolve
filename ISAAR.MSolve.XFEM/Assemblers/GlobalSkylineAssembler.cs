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
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Assemblers
{
    class GlobalSkylineAssembler
    {
        public (SkylineMatrix Kuu, CSRMatrix Kuc) BuildGlobalMatrix(Model2D model, IDOFEnumerator dofEnumerator)
        {
            int numDofsConstrained = dofEnumerator.ConstrainedDofsCount;
            int numDofsUnconstrained = dofEnumerator.FreeDofsCount + dofEnumerator.EnrichedDofsCount;

            // Rows, columns = standard free dofs + enriched dofs (aka the left hand side sub-matrix)
            SkylineBuilder Kuu = InitializeGlobalUncontrainedMatrix(model, dofEnumerator);

            // TODO: perhaps I should return a CSC matrix and do the transposed multiplication. This way I will not have to 
            // transpose the element matrix. Another approach is to add an AddTransposed() method to the DOK.
            var Kuc = DOKRowMajor.CreateEmpty(numDofsUnconstrained, numDofsConstrained);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Element matrices
                Matrix kss = element.BuildStandardStiffnessMatrix();
                element.BuildEnrichedStiffnessMatrices(out Matrix kes, out Matrix kee);

                // TODO: options: 1) Only work with upper triangle in all symmetric matrices. Same applies to Elements.
                // 2) The Elements have two versions of BuildStiffness(). 
                // 3) The Elements return both (redundant; If someone needs it he can make it himself like here) 
                Matrix kse = kes.Transpose();

                // Element to global dofs mappings
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                dofEnumerator.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapFree, out IReadOnlyDictionary<int, int> mapConstrained);
                IReadOnlyDictionary<int, int> mapEnriched =
                    dofEnumerator.MatchElementToGlobalEnrichedDofsOf(element);

                // Add the element contributions to the global matrices
                Kuu.AddSubmatrixSymmetric(kss, mapFree);
                Kuu.AddSubmatrix(kse, mapFree, mapEnriched);
                Kuu.AddSubmatrixSymmetric(kee, mapEnriched);

                Kuc.AddSubmatrix(kss, mapFree, mapConstrained);
                Kuc.AddSubmatrix(kes, mapEnriched, mapConstrained);
            }

            #region DEBUG code
            //(Matrix expectedKuu, Matrix expectedKuc) = DenseGlobalAssembler.BuildGlobalMatrix(model, dofEnumerator);
            //Console.WriteLine("Check Kuu:");
            //CheckMatrix(expectedKuu, Kuu);
            //Console.WriteLine("Check Kuc:");
            //CheckMatrix(expectedKuc, Kuc);
            #endregion

            //TODO: perhaps I should filter the matrices in the concrete class before returning (e.g. sok.Build())
            return (Kuu.BuildSkylineMatrix(), Kuc.BuildCSRMatrix(true));
        }


        private static SkylineBuilder InitializeGlobalUncontrainedMatrix(Model2D model, IDOFEnumerator dofEnumerator)
        {
            int numDofsUnconstrained = dofEnumerator.FreeDofsCount + dofEnumerator.EnrichedDofsCount;
            int[] colHeights = new int[numDofsUnconstrained]; //only entries above the diagonal count towards the column height
            foreach (XContinuumElement2D element in model.Elements)
            {
                // To determine the col height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0 entries in the element stiffness matrix.
                int minDOF = Int32.MaxValue;
                foreach (XNode2D node in element.Nodes)
                {
                    foreach (int dof in dofEnumerator.GetFreeDofsOf(node)) minDOF = Math.Min(dof, minDOF);
                    foreach (int dof in dofEnumerator.GetEnrichedDofsOf(node)) Math.Min(dof, minDOF);
                }

                // The height of each col is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (XNode2D node in element.Nodes)
                {
                    foreach (int dof in dofEnumerator.GetFreeDofsOf(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDOF);
                    }
                    foreach (int dof in dofEnumerator.GetEnrichedDofsOf(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDOF);
                    }
                }
            }
            return SkylineBuilder.Create(numDofsUnconstrained, colHeights);
        }
    }
}
