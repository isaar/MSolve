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
    public class ConcentratedMass3D_v2 : IStructuralFiniteElement_v2
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes };
        private readonly double massCoefficient;
        private IElementDofEnumerator_v2 dofEnumerator = new GenericDofEnumerator_v2();

        public int ID => 998;

        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

        public IElementDofEnumerator_v2 DofEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement_v2 element)
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

        public ConcentratedMass3D_v2(double massCoefficient)
        {
            this.massCoefficient = massCoefficient;
        }

        public ConcentratedMass3D_v2(double massCoefficient, IElementDofEnumerator_v2 dofEnumerator)
            : this(massCoefficient)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IMatrix MassMatrix(IElement_v2 element)
        {
            var mass = Matrix.CreateZero(3, 3);
            mass[0, 0] = massCoefficient;
            mass[1, 1] = massCoefficient;
            mass[2, 2] = massCoefficient;
            return mass;
        }

        public IMatrix StiffnessMatrix(IElement_v2 element) => Matrix.CreateZero(3, 3);

        public IMatrix DampingMatrix(IElement_v2 element) => Matrix.CreateZero(3, 3);

        public void ResetMaterialModified() { }

        public Tuple<double[], double[]> CalculateStresses(Element_v2 element, double[] localDisplacements, 
            double[] localdDisplacements)
            => new Tuple<double[], double[]>(new double[6], new double[6]);

        public double[] CalculateForcesForLogging(Element_v2 element, double[] localDisplacements)
            => CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

        public double[] CalculateForces(Element_v2 element, double[] localDisplacements, double[] localdDisplacements)
            => new double[6];

        public double[] CalculateAccelerationForces(Element_v2 element, IList<MassAccelerationLoad> loads)
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
