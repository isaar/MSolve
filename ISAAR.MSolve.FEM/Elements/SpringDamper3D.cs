using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Elements
{
    public enum SpringDirections
    {
        X = 0,
        Y,
        Z,
        XY,
        YZ,
        XZ,
        XYZ
    }

    public class SpringDamper3D : IStructuralFiniteElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private readonly double springCoefficient, dampingCoefficient;
        private readonly SpringDirections springDirections, dampingDirections;
        private IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public int ID
        {
            get { return 999; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public IElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element)
        {
            if (element == null) return dofs;

            var d = new List<IList<DOFType>>();
            foreach (var node in element.INodes)
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

        public SpringDamper3D(double springCoefficient, double dampingCoefficient, SpringDirections springDirections, SpringDirections dampingDirections)
        {
            this.springCoefficient = springCoefficient;
            this.dampingCoefficient = dampingCoefficient;
            this.springDirections = springDirections;
            this.dampingDirections = dampingDirections;
        }

        public SpringDamper3D(double springCoefficient, double dampingCoefficient, SpringDirections springDirections, SpringDirections dampingDirections, IElementDOFEnumerator dofEnumerator)
            : this(springCoefficient, dampingCoefficient, springDirections, dampingDirections)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            double x = (springDirections == SpringDirections.X || springDirections == SpringDirections.XY || springDirections == SpringDirections.XZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
            double y = (springDirections == SpringDirections.Y || springDirections == SpringDirections.XY || springDirections == SpringDirections.YZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
            double z = (springDirections == SpringDirections.Z || springDirections == SpringDirections.XZ || springDirections == SpringDirections.YZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
            return new SymmetricMatrix2D(new double[] { x, 0, 0, -x, 0, 0,
                y, 0, 0, -y, 0, 
                z, 0, 0, -z,
                x, 0, 0,
                y, 0,
                z
            });
        }

        public IMatrix2D MassMatrix(IElement element)
        {
            return new SymmetricMatrix2D(new double[] { 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 
                0, 0, 0, 0,
                0, 0, 0,
                0, 0,
                0
            });
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            double x = (dampingDirections == SpringDirections.X || dampingDirections == SpringDirections.XY || dampingDirections == SpringDirections.XZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;
            double y = (dampingDirections == SpringDirections.Y || dampingDirections == SpringDirections.XY || dampingDirections == SpringDirections.YZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;
            double z = (dampingDirections == SpringDirections.Z || dampingDirections == SpringDirections.XZ || dampingDirections == SpringDirections.YZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;
            return new SymmetricMatrix2D(new double[] { x, 0, 0, -x, 0, 0,
                y, 0, 0, -y, 0, 
                z, 0, 0, -z,
                x, 0, 0,
                y, 0,
                z
            });
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
            IMatrix2D stiffnessMatrix = StiffnessMatrix(element);
            var disps = new Vector(localDisplacements);
            double[] forces = new double[localDisplacements.Length];
            stiffnessMatrix.Multiply(disps, forces);
            return forces;
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            return new double[6];
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
