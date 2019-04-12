using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Elements
{
    public class ConcentratedMass3D : IStructuralFiniteElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes };
        private readonly double massCoefficient;
        private IElementDofEnumerator dofEnumerator = new GenericDofEnumerator();

        public int ID => 998;

        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

        public IElementDofEnumerator DofEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element)
        {
            if (element == null) return dofs;

            var d = new List<IList<DOFType>>();
            foreach (var node in element.Nodes)
            {
                var nodeDofs = new List<DOFType>();
                nodeDofs.AddRange(nodalDOFTypes);
                d.Add(nodeDofs);
            }
            return d;
        }

        public bool MaterialModified => false;

        public ConcentratedMass3D(double massCoefficient)
        {
            this.massCoefficient = massCoefficient;
        }

        public ConcentratedMass3D(double massCoefficient, IElementDofEnumerator dofEnumerator)
            : this(massCoefficient)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IMatrix MassMatrix(IElement element)
        {
            var mass = Matrix.CreateZero(3, 3);
            mass[0, 0] = massCoefficient;
            mass[1, 1] = massCoefficient;
            mass[2, 2] = massCoefficient;
            return mass;
        }

        public IMatrix StiffnessMatrix(IElement element) => Matrix.CreateZero(3, 3);

        public IMatrix DampingMatrix(IElement element) => Matrix.CreateZero(3, 3);

        public void ResetMaterialModified() { }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, 
            double[] localdDisplacements)
            => new Tuple<double[], double[]>(new double[6], new double[6]);

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
            => CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
            => new double[6];

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            var accelerations = new double[3];
            IMatrix massMatrix = MassMatrix(element);

            foreach (MassAccelerationLoad load in loads)
            {
                int index = 0;
                foreach (DOFType[] nodalDOFTypes in dofs)
                    foreach (DOFType dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }
            }

            return massMatrix.Multiply(accelerations);
        }

        public void ClearMaterialState() { }
        public void SaveMaterialState() { }
        public void ClearMaterialStresses() { }
    }
}
