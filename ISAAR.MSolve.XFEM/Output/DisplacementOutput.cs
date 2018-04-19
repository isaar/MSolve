using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Output
{
    class DisplacementOutput
    {
        private readonly Model2D model;
        //private DOFEnumerator dofEnumerator;

        public DisplacementOutput(Model2D model)
        {
            this.model = model;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="model"></param>
        /// <param name="solution"></param>
        /// <returns>A nodesCount x 2 array, where each row stores the x and y displacements of that node</returns>
        public double[,] FindNodalDisplacements(Vector solution)
        {
            return model.DofEnumerator.GatherNodalDisplacements(model, solution);
        }

        public IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<Vector>> FindElementWiseDisplacements(
            Vector solution)
        {
            Vector constrainedDisplacements = model.CalculateConstrainedDisplacements();
            var allDisplacements = new Dictionary<XContinuumElement2D, IReadOnlyList<Vector>>();
            foreach (var element in model.Elements)
            {
                Vector displacementsUnrolled = model.DofEnumerator.ExtractDisplacementVectorOfElementFromGlobal(
                    element, solution, constrainedDisplacements);
                var displacementsAsVectors = new Vector[element.Nodes.Count];
                for (int i = 0; i < element.Nodes.Count; ++i)
                {
                    displacementsAsVectors[i] = Vector.CreateFromArray(new double[]
                    {
                        displacementsUnrolled[2 * i], displacementsUnrolled[2 * i + 1] // This only works for continuum elements though.
                    });
                }
                allDisplacements[element] = displacementsAsVectors;
            }
            return allDisplacements;
        }
    }
}
