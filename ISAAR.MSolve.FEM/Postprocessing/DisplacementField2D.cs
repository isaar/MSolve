using ISAAR.MSolve.FEM.Entities;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.FEM.Postprocessing
{
    /// <summary>
    /// Recovers the nodal displacements from the solution of an analysis step. For now it only works for linear analysis.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class DisplacementField2D
    {
        private readonly Dictionary<Node_v2, double[]> data;
        private readonly Model_v2 model;

        public DisplacementField2D(Model_v2 model)
        {
            this.model = model;
            this.data = new Dictionary<Node_v2, double[]>(model.Nodes.Count);
        }

        public void FindNodalDisplacements(IVectorView solution)
        {
            foreach (var idxNodePair in model.NodesDictionary)
            {
                Node_v2 node = idxNodePair.Value;
                //if (nodalDofs.Count != 2) throw new Exception("There must be exactly 2 dofs per node, X and Y");
                bool isFree = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(node, DOFType.X, out int dofXIdx);
                double ux = isFree ? solution[dofXIdx] : 0.0;
                isFree = model.GlobalDofOrdering.GlobalFreeDofs.TryGetValue(node, DOFType.Y, out int dofYIdx);
                double uy = isFree ? solution[dofYIdx] : 0.0;
                data.Add(idxNodePair.Value, new double[] { ux, uy });
            }
        }

        public double[] this[Node_v2 node] => data[node];
    }
}
