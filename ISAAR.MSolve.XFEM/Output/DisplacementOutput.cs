using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Output
{
    class DisplacementOutput
    {
        private readonly Model2D model;
        private readonly IDofOrderer dofOrderer;

        public DisplacementOutput(Model2D model, IDofOrderer dofOrderer)
        {
            this.model = model;
            this.dofOrderer = dofOrderer;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="model"></param>
        /// <param name="solution"></param>
        /// <returns>A nodesCount x 2 array, where each row stores the x and y displacements of that node</returns>
        public double[,] FindNodalDisplacements(Vector solution)
        {
            return dofOrderer.GatherNodalDisplacements(model, solution);
        }

        public IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<Vector2>> FindElementWiseDisplacements(
            Vector solution)
        {
            Vector constrainedDisplacements = model.CalculateConstrainedDisplacements(dofOrderer);
            var allDisplacements = new Dictionary<XContinuumElement2D, IReadOnlyList<Vector2>>();
            foreach (var element in model.Elements)
            {
                Vector displacementsUnrolled = dofOrderer.ExtractDisplacementVectorOfElementFromGlobal(
                    element, solution, constrainedDisplacements);
                var displacementsAsVectors = new Vector2[element.Nodes.Count];
                for (int i = 0; i < element.Nodes.Count; ++i) // This only works for continuum elements though.
                {
                    displacementsAsVectors[i] = Vector2.Create(displacementsUnrolled[2 * i], displacementsUnrolled[2 * i + 1]); 
                }
                allDisplacements[element] = displacementsAsVectors;
            }
            return allDisplacements;
        }
    }
}
