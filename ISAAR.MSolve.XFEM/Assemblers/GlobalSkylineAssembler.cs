using System;
using System.Collections.Generic;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Assemblers
{
    class GlobalSkylineAssembler
    {
        public (SkylineMatrix Kff, DokRowMajor Kfc) BuildGlobalMatrix(Model2D model, IDofOrderer dofOrderer)
        {
            int numDofsConstrained = dofOrderer.NumConstrainedDofs;
            int numDofsFree = dofOrderer.NumStandardDofs + dofOrderer.NumEnrichedDofs;
            SkylineBuilder Kff = FindSkylineColumnHeights(model, dofOrderer);
            var Kfc = DokRowMajor.CreateEmpty(numDofsFree, numDofsConstrained);

            foreach (XContinuumElement2D element in model.Elements)
            {
                // Build standard element matrix and add its contributions to the global matrices.
                // TODO: perhaps that could be done and cached during the dof enumeration to avoid iterating over the dofs twice
                dofOrderer.MatchElementToGlobalStandardDofsOf(element,
                    out IReadOnlyDictionary<int, int> mapStandard, out IReadOnlyDictionary<int, int> mapConstrained);
                Matrix kss = element.BuildStandardStiffnessMatrix();
                Kff.AddSubmatrixSymmetric(kss, mapStandard);
                Kfc.AddSubmatrix(kss, mapStandard, mapConstrained);

                // Build enriched element matrices and add their contributions to the global matrices
                IReadOnlyDictionary<int, int> mapEnriched = dofOrderer.MatchElementToGlobalEnrichedDofsOf(element);
                if (mapEnriched.Count > 0)
                {
                    (Matrix kee, Matrix kse) = element.BuildEnrichedStiffnessMatricesUpper(); ;
                    Kff.AddSubmatrix(kse, mapStandard, mapEnriched);
                    Kff.AddSubmatrixSymmetric(kee, mapEnriched);
                    if (mapConstrained.Count > 0) Kfc.AddSubmatrix(kse.Transpose(), mapEnriched, mapConstrained);
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
            return (Kff.BuildSkylineMatrix(), Kfc);
        }


        private static SkylineBuilder FindSkylineColumnHeights(Model2D model, IDofOrderer dofOrderer)
        {
            int numDofsFree = dofOrderer.NumStandardDofs + dofOrderer.NumEnrichedDofs;
            int[] colHeights = new int[numDofsFree]; //only entries above the diagonal count towards the column height
            foreach (XContinuumElement2D element in model.Elements)
            {
                //TODO: perhaps the 2 outer loops could be done at once to avoid a lot of dof indexing. Could I update minDof
                //      and colHeights[] at once?

                // To determine the col height, first find the min of the dofs of this element. All these are 
                // considered to interact with each other, even if there are 0 entries in the element stiffness matrix.
                int minDof = Int32.MaxValue;
                foreach (XNode2D node in element.Nodes)
                {
                    foreach (int dof in dofOrderer.GetStandardDofsOf(node)) minDof = Math.Min(dof, minDof);
                    foreach (int dof in dofOrderer.GetEnrichedDofsOf(node)) minDof = Math.Min(dof, minDof);
                }

                // The height of each col is updated for all elements that engage the corresponding dof. 
                // The max height is stored.
                foreach (XNode2D node in element.Nodes)
                {
                    foreach (int dof in dofOrderer.GetStandardDofsOf(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDof);
                    }
                    foreach (int dof in dofOrderer.GetEnrichedDofsOf(node))
                    {
                        colHeights[dof] = Math.Max(colHeights[dof], dof - minDof);
                    }
                }
            }
            return SkylineBuilder.Create(numDofsFree, colHeights);
        }
    }
}
