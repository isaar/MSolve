using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Testing.Utilities;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Assemblers
{
    class GlobalSkylineAssembler: GlobalAssemblerBase<SkylineBuilder>
    {
        protected override SkylineBuilder InitializeGlobalUncontrainedMatrix(Model2D model, IDOFEnumerator dofEnumerator)
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
