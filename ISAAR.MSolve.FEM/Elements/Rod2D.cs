using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Problems.Structural.Elements
{
    public class Rod2D : IStructuralFiniteElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[2] { DOFType.X, DOFType.Y };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private readonly double youngModulus;
        private IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public double Density { get; set; }
        public double SectionArea { get; set; }

        public Rod2D(double youngModulus)
        {
            this.youngModulus = youngModulus;
        }

        public Rod2D(double youngModulus, IElementDOFEnumerator dofEnumerator)
            : this(youngModulus)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public IMatrix2D TransformationMatrix(Element element)
        {
            double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
            double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);
            double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
            double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;

            double[,] transformation = { {  c,   s, 0.0, 0.0},
                                         {0.0, 0.0,   c,   s}};
            return new Matrix2D(transformation);
        }

        /// <summary>
        /// Stress0         Stress1
        /// -> ------------ ->
        /// </summary>
        /// <param name="element"></param>
        /// <param name="localDisplacements"></param>
        /// <param name="local_d_Displacements"></param>
        /// <returns></returns>
        public double CalculateAxialStress(Element element, double[] localDisplacements, double[] local_d_Displacements)
        {
            double[] globalStresses = CalculateStresses(element, localDisplacements, local_d_Displacements).Item2; // item1 = strains
            IMatrix2D transformation = TransformationMatrix(element);
            double[] localStresses = new double[2]; // In local natural system there are 2 dofs
            transformation.Multiply(new Vector(globalStresses), localStresses);
            // If Stress1 = localStresses[1] > 0 => tension. Else compression
            return localStresses[1];
        }

        #region IElementType Members

        public int ID
        {
            get { return 1; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.TwoD; }
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element)
        {
            return dofs;
        }

        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            return element.Nodes;
        }

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            double x2 = Math.Pow(element.INodes[1].X - element.INodes[0].X, 2);
            double y2 = Math.Pow(element.INodes[1].Y - element.INodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);
            double c = (element.INodes[1].X - element.INodes[0].X) / L;
            double c2 = c * c;
            double s = (element.INodes[1].Y - element.INodes[0].Y) / L;
            double s2 = s * s;
            double cs = c * s;
            double E = this.youngModulus;
            double A = SectionArea;

            return dofEnumerator.GetTransformedMatrix(
                new Matrix2D(new double[,]
                {
                    {A*E*c2/L, A*E*cs/L, -A*E*c2/L, -A*E*cs/L },
                    {A*E*cs/L, A*E*s2/L, -A*E*cs/L, -A*E*s2/L },
                    {-A*E*c2/L, -A*E*cs/L, A*E*c2/L, A*E*cs/L },
                    {-A*E*cs/L, -A*E*s2/L, A*E*cs/L, A*E*s2/L }
                }));
        }

        public IMatrix2D MassMatrix(IElement element)
        {
            //TODO: This is the lumped mass. In continuum elements 2D and 3D we use the consistent mass.
            double x2 = Math.Pow(element.INodes[1].X - element.INodes[0].X, 2);
            double y2 = Math.Pow(element.INodes[1].Y - element.INodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);

            double totalMass = Density * SectionArea * L;

            return new Matrix2D(new double[,]
            {
                { totalMass/2,0,0,0},
                {0,totalMass/2,0,0 },
                {0,0,totalMass/2,0 },
                {0,0,0,totalMass/2 }
            });
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] local_Displacements, double[] local_d_Displacements)
        {
            // WARNING: 1) No strains are computed 2) localdDisplacements are not used.
            double[] strains = null;
            double[] forces = CalculateForces(element, local_Displacements, local_d_Displacements);
            double[] stresses = Array.ConvertAll<double, double>(forces, x => x / SectionArea);
            return new Tuple<double[], double[]>(strains, stresses);
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            IMatrix2D stiffness = StiffnessMatrix(element);
            double[] forces = new double[localDisplacements.Length];
            stiffness.Multiply(new Vector(localDisplacements), forces);
            return forces;
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            Vector accelerations = new Vector(4);
            IMatrix2D massMatrix = MassMatrix(element);

            int index = 0;
            foreach (MassAccelerationLoad load in loads)
                foreach (DOFType[] nodalDOFTypes in dofs)
                    foreach (DOFType dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }

            double[] forces = new double[4];
            massMatrix.Multiply(accelerations, forces);
            return forces;
        }

        public void SaveMaterialState()
        {
            throw new NotImplementedException();
        }

        #endregion

        #region IFiniteElement Members


        public bool MaterialModified
        {
            get { return false; }
        }

        public void ResetMaterialModified()
        {
        }

        #endregion

        #region IFiniteElement Members

        public void ClearMaterialState()
        {
        }

        public void ClearMaterialStresses()
        {
            throw new NotImplementedException();
        }

        #endregion
    }
}
