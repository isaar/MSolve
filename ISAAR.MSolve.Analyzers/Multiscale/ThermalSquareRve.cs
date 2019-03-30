using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: How can we ensure that the model has the correct shape and discretization?
namespace ISAAR.MSolve.Analyzers.Multiscale
{
    /// <summary>
    /// This works only for square (2D) RVEs, centered at (0,0), with 1 thermal dof per node and linear boundary conditions.
    /// </summary>
    public class ThermalSquareRve : IReferenceVolumeElement
    {
        private const int numDimensions = 2;

        private readonly Vector2 macroscopicTemperatureGradient;
        private readonly IStructuralModel_v2 model;
        private readonly double xMin, yMin, xMax, yMax;
        private readonly double thickness;
        private readonly HashSet<INode> leftNodes, rightNodes, bottomNodes, topNodes;

        /// <summary>
        /// </summary>
        /// <param name="model"></param>
        /// <param name="bottomLeftCoords"></param>
        /// <param name="topRightCoords"></param>
        /// <param name="thickness"></param>
        /// <param name="macroscopicTemperatureGradient"></param>
        /// <param name="meshTolerance">The default is 1E-10 * min(|xMax-xMin|, |yMax-yMin|)</param>
        public ThermalSquareRve(IStructuralModel_v2 model, Vector2 bottomLeftCoords, Vector2 topRightCoords, double thickness,
            Vector2 macroscopicTemperatureGradient, double meshTolerance)
        {
            this.model = model;
            this.xMin = bottomLeftCoords[0];
            this.yMin = bottomLeftCoords[1];
            this.xMax = topRightCoords[0];
            this.yMax = topRightCoords[1];
            this.thickness = thickness;
            this.macroscopicTemperatureGradient = macroscopicTemperatureGradient;

            // Find the nodes of each edge
            leftNodes = new HashSet<INode>();
            rightNodes = new HashSet<INode>();
            bottomNodes = new HashSet<INode>();
            topNodes = new HashSet<INode>();
            foreach (INode node in model.Nodes)
            {
                // Top and right edges are prioritized for corner nodes. //TODO: should the corner nodes be handled differently?
                if (Math.Abs(node.Y - yMax) <= meshTolerance) topNodes.Add(node);
                else if (Math.Abs(node.X - xMax) <= meshTolerance) rightNodes.Add(node);
                else if (Math.Abs(node.Y - yMin) <= meshTolerance) bottomNodes.Add(node);
                else if (Math.Abs(node.X - xMin) <= meshTolerance) leftNodes.Add(node);
            }
        }

        public ThermalSquareRve(IStructuralModel_v2 model, Vector2 bottomLeftCoords, Vector2 topRightCoords, double thickness,
            Vector2 macroscopicTemperatureGradient) : 
            this(model, bottomLeftCoords, topRightCoords, thickness, macroscopicTemperatureGradient, 
                1E-10 * topRightCoords.Subtract(bottomLeftCoords).MinAbsolute())
        {
        }

        public void ApplyBoundaryConditions()
        {
            double dTdx = macroscopicTemperatureGradient[0];
            double dTdy = macroscopicTemperatureGradient[1];
            foreach (var node in leftNodes) node.Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 0.0 });
            foreach (var node in bottomNodes) node.Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = 0.0 });
            foreach (var node in rightNodes) node.Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = dTdx });
            foreach (var node in topNodes) node.Constraints.Add(new Constraint() { DOF = DOFType.Temperature, Amount = dTdy });
        }

        public double CalculateRveVolume() => (xMax - xMin) * (yMax - yMin) * thickness;

        public IMatrixView CalculateKinematicRelationsMatrix(ISubdomain_v2 subdomain)
        {
            ISubdomainConstrainedDofOrdering constrainedDofOrdering = subdomain.ConstrainedDofOrdering;
            var kinematicRelations = Matrix.CreateZero(numDimensions, constrainedDofOrdering.NumConstrainedDofs);
            CalculateKinematicsOfEdge(constrainedDofOrdering, bottomNodes, kinematicRelations);
            CalculateKinematicsOfEdge(constrainedDofOrdering, leftNodes, kinematicRelations);
            CalculateKinematicsOfEdge(constrainedDofOrdering, rightNodes, kinematicRelations);
            CalculateKinematicsOfEdge(constrainedDofOrdering, topNodes, kinematicRelations);
            return kinematicRelations;
        }

        private void CalculateKinematicsOfEdge(ISubdomainConstrainedDofOrdering constrainedDofOrdering,
            IEnumerable<INode> edgeNodes, Matrix kinematicRelations)
        {
            foreach (INode node in edgeNodes)
            {
                int dofIdx = constrainedDofOrdering.ConstrainedDofs[node, DOFType.Temperature];
                kinematicRelations[0, dofIdx] = node.X;
                kinematicRelations[1, dofIdx] = node.Y;
            }
        }
    }
}
