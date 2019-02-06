using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Problems.Structural.Elements
{
    public class Rod2D_v2 : IStructuralFiniteElement_v2
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[2] { DOFType.X, DOFType.Y };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private readonly double youngModulus;
        private IElementDofEnumerator_v2 dofEnumerator = new GenericDofEnumerator_v2();

        public double Density { get; set; }
        public double SectionArea { get; set; }

        public Rod2D_v2(double youngModulus)
        {
            this.youngModulus = youngModulus;
        }

        public Rod2D_v2(double youngModulus, IElementDofEnumerator_v2 dofEnumerator)
            : this(youngModulus)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IElementDofEnumerator_v2 DofEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        //TODO: this should be either cached, or even better the calculations should be incorporated into Stiffness()
        public IMatrix TransformationMatrix(Element_v2 element)
        {
            double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
            double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);
            double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
            double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;

            // T = [ cos sin 0 0; 0 0 cos sin]
            var transformation = Matrix.CreateZero(2, 4);
            transformation[0, 0] = c;
            transformation[0, 1] = s;
            transformation[1, 2] = c;
            transformation[1, 3] = s;
            return transformation;
        }

        /// <summary>
        /// Stress0         Stress1
        /// -> ------------ ->
        /// </summary>
        /// <param name="element"></param>
        /// <param name="localDisplacements"></param>
        /// <param name="local_d_Displacements"></param>
        /// <returns></returns>
        public double CalculateAxialStress(Element_v2 element, double[] localDisplacements, double[] local_d_Displacements)
        {
            double[] globalStresses = CalculateStresses(element, localDisplacements, local_d_Displacements).Item2; // item1 = strains
            IMatrix transformation = TransformationMatrix(element);
            double[] localStresses = transformation.Multiply(globalStresses); // In local natural system there are 2 dofs
            // If Stress1 = localStresses[1] > 0 => tension. Else compression
            return localStresses[1];
        }

        #region IElementType Members

        public int ID => 1;

        public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

        public IList<IList<DOFType>> GetElementDOFTypes(IElement_v2 element) => dofs;

        public IList<Node_v2> GetNodesForMatrixAssembly(Element_v2 element) => element.Nodes;

        public IMatrix StiffnessMatrix(IElement_v2 element)
        {
            double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
            double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);
            double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
            double c2 = c * c;
            double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
            double s2 = s * s;
            double cs = c * s;
            double E = this.youngModulus;
            double A = SectionArea;

            return dofEnumerator.GetTransformedMatrix(
                Matrix.CreateFromArray(new double[,]
                {
                    {A*E*c2/L, A*E*cs/L, -A*E*c2/L, -A*E*cs/L },
                    {A*E*cs/L, A*E*s2/L, -A*E*cs/L, -A*E*s2/L },
                    {-A*E*c2/L, -A*E*cs/L, A*E*c2/L, A*E*cs/L },
                    {-A*E*cs/L, -A*E*s2/L, A*E*cs/L, A*E*s2/L }
                }));
        }

        public IMatrix MassMatrix(IElement_v2 element)
        {
            double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
            double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);

            double totalMassOver2 = Density * SectionArea * L / 2.0;

            // Lumped mass: M = [m/2 0 0 0; 0 m/2 0 0; 0 0 m/2 0; 0 0 0 m/2]
            int order = 4;
            var lumpedMass = Matrix.CreateZero(order, order);
            for (int i = 0; i < order; ++i) lumpedMass[i, i] = totalMassOver2;
            return lumpedMass;
        }

        public IMatrix DampingMatrix(IElement_v2 element) => throw new NotImplementedException();

        public Tuple<double[], double[]> CalculateStresses(Element_v2 element, double[] local_Displacements, 
            double[] local_d_Displacements)
        {
            // WARNING: 1) No strains are computed 2) localdDisplacements are not used.
            double[] strains = null;
            double[] forces = CalculateForces(element, local_Displacements, local_d_Displacements);
            double[] stresses = Array.ConvertAll<double, double>(forces, x => x / SectionArea);
            return new Tuple<double[], double[]>(strains, stresses);
        }

        public double[] CalculateForcesForLogging(Element_v2 element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public double[] CalculateForces(Element_v2 element, double[] localDisplacements, double[] localdDisplacements)
        {
            IMatrix stiffness = StiffnessMatrix(element);
            return stiffness.Multiply(localdDisplacements);
        }

        public double[] CalculateAccelerationForces(Element_v2 element, IList<MassAccelerationLoad> loads)
        {
            var accelerations = new double[4];
            IMatrix massMatrix = MassMatrix(element);

            int index = 0;
            foreach (MassAccelerationLoad load in loads)
                foreach (DOFType[] nodalDOFTypes in dofs)
                    foreach (DOFType dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }

            return massMatrix.Multiply(accelerations);
        }

        public void SaveMaterialState() { }

        #endregion

        #region IFiniteElement Members


        public bool MaterialModified => false;

        public void ResetMaterialModified() {}

        #endregion

        #region IFiniteElement Members

        public void ClearMaterialState() {}

        public void ClearMaterialStresses() => throw new NotImplementedException();

        #endregion
    }
}
