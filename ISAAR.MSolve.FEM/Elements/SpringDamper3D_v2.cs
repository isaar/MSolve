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
    public enum SpringDirections
    {
        X = 0, Y, Z, XY, YZ, XZ, XYZ
    }

    public class SpringDamper3D_v2 : IStructuralFiniteElement_v2
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private readonly double springCoefficient, dampingCoefficient;
        private readonly SpringDirections springDirections, dampingDirections;
        private IElementDofEnumerator_v2 dofEnumerator = new GenericDofEnumerator_v2();

        public int ID => 999;

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

        public SpringDamper3D_v2(double springCoefficient, double dampingCoefficient, SpringDirections springDirections, 
            SpringDirections dampingDirections)
        {
            this.springCoefficient = springCoefficient;
            this.dampingCoefficient = dampingCoefficient;
            this.springDirections = springDirections;
            this.dampingDirections = dampingDirections;
        }

        public SpringDamper3D_v2(double springCoefficient, double dampingCoefficient, SpringDirections springDirections, 
            SpringDirections dampingDirections, IElementDofEnumerator_v2 dofEnumerator)
            : this(springCoefficient, dampingCoefficient, springDirections, dampingDirections)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IMatrix StiffnessMatrix(IElement_v2 element)
        {
            double x = (springDirections == SpringDirections.X || springDirections == SpringDirections.XY || springDirections == SpringDirections.XZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
            double y = (springDirections == SpringDirections.Y || springDirections == SpringDirections.XY || springDirections == SpringDirections.YZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
            double z = (springDirections == SpringDirections.Z || springDirections == SpringDirections.XZ || springDirections == SpringDirections.YZ || springDirections == SpringDirections.XYZ) ? springCoefficient : 0;
            return Matrix.CreateFromArray(new double[,]
                {
                   { x, 0, 0, -x, 0, 0 },
                   { 0, y, 0, 0, -y, 0 },
                   { 0, 0, z, 0, 0, -z },
                   {-x, 0, 0, x, 0, 0 },
                   { 0,-y, 0, 0, y, 0 },
                   { 0, 0,-z, 0, 0, z }
                }
                );

            //return SymmetricMatrix.CreateFromArray(new double[] 
            //{
            //    x, 0, 0, -x, 0, 0,
            //       y, 0, 0, -y, 0, 
            //          z, 0, 0, -z,
            //             x, 0, 0,
            //                y, 0,
            //                   z
            //});
        }

        public IMatrix MassMatrix(IElement_v2 element) => Matrix.CreateZero(6, 6);

        public IMatrix DampingMatrix(IElement_v2 element)
        {
            double x = (dampingDirections == SpringDirections.X || dampingDirections == SpringDirections.XY || dampingDirections == SpringDirections.XZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;
            double y = (dampingDirections == SpringDirections.Y || dampingDirections == SpringDirections.XY || dampingDirections == SpringDirections.YZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;
            double z = (dampingDirections == SpringDirections.Z || dampingDirections == SpringDirections.XZ || dampingDirections == SpringDirections.YZ || dampingDirections == SpringDirections.XYZ) ? dampingCoefficient : 0;

            return Matrix.CreateFromArray(new double[,]
                {
                   { x, 0, 0, -x, 0, 0 },
                   { 0, y, 0, 0, -y, 0 },
                   { 0, 0, z, 0, 0, -z },
                   {-x, 0, 0, x, 0, 0 },
                   { 0,-y, 0, 0, y, 0 },
                   { 0, 0,-z, 0, 0, z }
                }
                );

            //return SymmetricMatrix.CreateFromArray(new double[] 
            //{
            //    x, 0, 0, -x, 0, 0,
            //    y, 0, 0, -y, 0, 
            //    z, 0, 0, -z,
            //    x, 0, 0,
            //    y, 0,
            //    z
            //});
        }

        public void ResetMaterialModified() { }

        public Tuple<double[], double[]> CalculateStresses(Element_v2 element, double[] localDisplacements, 
            double[] localdDisplacements)
            => new Tuple<double[], double[]>(new double[6], new double[6]);

        public double[] CalculateForcesForLogging(Element_v2 element, double[] localDisplacements)
            => CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);

        public double[] CalculateForces(Element_v2 element, double[] localDisplacements, double[] localdDisplacements)
        {
            IMatrix stiffnessMatrix = StiffnessMatrix(element);
            return stiffnessMatrix.Multiply(localDisplacements);
        }

        public double[] CalculateAccelerationForces(Element_v2 element, IList<MassAccelerationLoad> loads) => new double[6];

        public void ClearMaterialState() { }
        public void SaveMaterialState() { }
        public void ClearMaterialStresses() { }
    }
}
