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
                int dofXIdx = model.GlobalDofOrdering.GlobalFreeDofs[node, DOFType.X];
                double ux = (dofXIdx != Model.constrainedDofIdx) ? solution[dofXIdx] : 0.0;
                int dofYIdx = model.GlobalDofOrdering.GlobalFreeDofs[node, DOFType.Y];
                double uy = (dofYIdx != Model.constrainedDofIdx) ? solution[dofYIdx] : 0.0;
                data.Add(idxNodePair.Value, new double[] { ux, uy });
            }
        }

        public double[] this[Node_v2 node] => data[node];
    }
}
