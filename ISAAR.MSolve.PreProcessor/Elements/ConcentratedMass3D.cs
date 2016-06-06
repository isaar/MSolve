using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.Matrices;

namespace ISAAR.MSolve.PreProcessor.Elements
{
    public class ConcentratedMass3D : IStructuralFiniteElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes };
        private readonly double massCoefficient;
        private IFiniteElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public int ID
        {
            get { return 998; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public IFiniteElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public IList<IList<DOFType>> GetElementDOFTypes(Element element)
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

        public bool MaterialModified
        {
            get { return false; }
        }

        public ConcentratedMass3D(double massCoefficient)
        {
            this.massCoefficient = massCoefficient;
        }

        public ConcentratedMass3D(double massCoefficient, IFiniteElementDOFEnumerator dofEnumerator)
            : this(massCoefficient)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IMatrix2D<double> MassMatrix(Element element)
        {
            return new SymmetricMatrix2D<double>(new double[] { massCoefficient, 0, 0,
                massCoefficient, 0,
                massCoefficient
            });
        }

        public IMatrix2D<double> StiffnessMatrix(Element element)
        {
            return new SymmetricMatrix2D<double>(new double[6]);
        }

        public IMatrix2D<double> DampingMatrix(Element element)
        {
            return new SymmetricMatrix2D<double>(new double[6]);
        }

        public void ResetMaterialModified()
        {
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            return new Tuple<double[], double[]>(new double[6], new double[6]);
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            return new double[6];
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            Vector<double> accelerations = new Vector<double>(3);
            IMatrix2D<double> massMatrix = MassMatrix(element);

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
            double[] forces = new double[3];
            massMatrix.Multiply(accelerations, forces);
            return forces;
        }

        public void ClearMaterialState()
        {
        }

        public void SaveMaterialState()
        {
        }

        public void ClearMaterialStresses()
        {
        }
    }
}
