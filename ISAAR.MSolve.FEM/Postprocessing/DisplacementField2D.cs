using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

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

        public void FindNodalDisplacements(IVector solution)
        {
            foreach (var idxNodePair in model.NodesDictionary)
            {
                Dictionary<DOFType, int> nodalDofs = model.NodalDOFsDictionary[idxNodePair.Key];
                if (nodalDofs.Count != 2) throw new Exception("There must be exactly 2 dofs per node, X and Y");
                int dofXIdx = nodalDofs[DOFType.X];
                double ux = (dofXIdx != Model.constrainedDofIdx) ? solution[dofXIdx] : 0.0;
                int dofYIdx = nodalDofs[DOFType.Y];
                double uy = (dofYIdx != Model.constrainedDofIdx) ? solution[dofYIdx] : 0.0;
                data.Add(idxNodePair.Value, new double[] { ux, uy });
            }
        }

        public double[] this[Node node] => data[node];
    }
}
