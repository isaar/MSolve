using System.Collections.Generic;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.FEM.Postprocessing
{
    /// <summary>
    /// Recovers the nodal displacements from the solution of an analysis step. For now it only works for linear analysis.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class DisplacementField2D
    {
        private readonly Dictionary<Node, double[]> data;
        private readonly Model model;

        public DisplacementField2D(Model model)
        {
            this.model = model;
            this.data = new Dictionary<Node, double[]>(model.Nodes.Count);
        }

        public void FindNodalDisplacements(IVectorView solution)
        {
            foreach (var idxNodePair in model.NodesDictionary)
            {
                Node node = idxNodePair.Value;
                //if (nodalDofs.Count != 2) throw new Exception("There must be exactly 2 dofs per node, X and Y");
                bool isFree = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(node, StructuralDof.TranslationX, out int dofXIdx);
                double ux = isFree ? solution[dofXIdx] : 0.0;
                isFree = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(node, StructuralDof.TranslationY, out int dofYIdx);
                double uy = isFree ? solution[dofYIdx] : 0.0;
                data.Add(idxNodePair.Value, new double[] { ux, uy });
            }
        }

        public double[] this[Node node] => data[node];
    }
}
