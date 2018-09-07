using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM.Elements
{
    public class EulerBeam2D : IStructuralFiniteElement
    {
        private static readonly DOFType[] nodalDOFTypes = new DOFType[3] { DOFType.X, DOFType.Y, DOFType.RotZ };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private readonly double youngModulus;
        private IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        public double Density { get; set; }
        public double SectionArea { get; set; }
        public double MomentOfInertia { get; set; }
        public double RayleighAlpha { get; set; }
        public double RayleighBeta { get; set; }

        public EulerBeam2D(double youngModulus)
        {
            this.youngModulus = youngModulus;
        }

        public EulerBeam2D(double youngModulus, IElementDOFEnumerator dofEnumerator)
            : this(youngModulus)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
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

        //[  c^2*E*A/L+12*s^2*E*I/L^3,  s*E*A/L*c-12*c*E*I/L^3*s,              -6*E*I/L^2*s, -c^2*E*A/L-12*s^2*E*I/L^3, -s*E*A/L*c+12*c*E*I/L^3*s,              -6*E*I/L^2*s]
        //[  s*E*A/L*c-12*c*E*I/L^3*s,  s^2*E*A/L+12*c^2*E*I/L^3,               6*E*I/L^2*c, -s*E*A/L*c+12*c*E*I/L^3*s, -s^2*E*A/L-12*c^2*E*I/L^3,               6*E*I/L^2*c]
        //[              -6*E*I/L^2*s,               6*E*I/L^2*c,                   4*E*I/L,               6*E*I/L^2*s,              -6*E*I/L^2*c,                   2*E*I/L]
        //[ -c^2*E*A/L-12*s^2*E*I/L^3, -s*E*A/L*c+12*c*E*I/L^3*s,               6*E*I/L^2*s,  c^2*E*A/L+12*s^2*E*I/L^3,  s*E*A/L*c-12*c*E*I/L^3*s,               6*E*I/L^2*s]
        //[ -s*E*A/L*c+12*c*E*I/L^3*s, -s^2*E*A/L-12*c^2*E*I/L^3,              -6*E*I/L^2*c,  s*E*A/L*c-12*c*E*I/L^3*s,  s^2*E*A/L+12*c^2*E*I/L^3,              -6*E*I/L^2*c]
        //[              -6*E*I/L^2*s,               6*E*I/L^2*c,                   2*E*I/L,               6*E*I/L^2*s,              -6*E*I/L^2*c,                   4*E*I/L]
        public virtual IMatrix2D StiffnessMatrix(IElement element)
        {
            double x2 = Math.Pow(element.INodes[1].X - element.INodes[0].X, 2);
            double y2 = Math.Pow(element.INodes[1].Y - element.INodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);
            double c = (element.INodes[1].X - element.INodes[0].X) / L;
            double c2 = c * c;
            double s = (element.INodes[1].Y - element.INodes[0].Y) / L;
            double s2 = s * s;
            double EL = this.youngModulus / L;
            double EAL = EL * SectionArea;
            double EIL = EL * MomentOfInertia;
            double EIL2 = EIL / L;
            double EIL3 = EIL2 / L;
            return dofEnumerator.GetTransformedMatrix(new SymmetricMatrix2D(new double[] { c2*EAL+12*s2*EIL3, c*s*EAL-12*c*s*EIL3, -6*s*EIL2, -c2*EAL-12*s2*EIL3, -c*s*EAL+12*c*s*EIL3, -6*s*EIL2,
                s2*EAL+12*c2*EIL3, 6*c*EIL2, -s*c*EAL+12*c*s*EIL3, -s2*EAL-12*c2*EIL3, 6*c*EIL2,
                4*EIL, 6*s*EIL2, -6*c*EIL2, 2*EIL,
                c2*EAL+12*s2*EIL3, s*c*EAL-12*c*s*EIL3, 6*s*EIL2,
                s2*EAL+12*c2*EIL3, -6*c*EIL2,
                4*EIL }));
        }

        ////[ 140*c^2+156*s^2,         -16*c*s,         -22*s*L,   70*c^2+54*s^2,          16*c*s,          13*s*L]
        ////[         -16*c*s, 140*s^2+156*c^2,          22*c*L,          16*c*s,   70*s^2+54*c^2,         -13*c*L]
        ////[         -22*s*L,          22*c*L,           4*L^2,         -13*s*L,          13*c*L,          -3*L^2]
        ////[   70*c^2+54*s^2,          16*c*s,         -13*s*L, 140*c^2+156*s^2,         -16*c*s,          22*s*L]
        ////[          16*c*s,   70*s^2+54*c^2,          13*c*L,         -16*c*s, 140*s^2+156*c^2,         -22*c*L]
        ////[          13*s*L,         -13*c*L,          -3*L^2,          22*s*L,         -22*c*L,           4*L^2]
        
        //public IMatrix2D<double> MassMatrix(Element element)
        //{
        //    double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
        //    double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
        //    double L = Math.Sqrt(x2 + y2);
        //    double L2 = L * L;
        //    double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
        //    double c2 = c * c;
        //    double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
        //    double s2 = s * s;
        //    double dAL420 = Density * SectionArea * L / 420;
        //    return new SymmetricMatrix2D<double>(new double[] { dAL420*(140*c2+156*s2), -16*dAL420*c*s, -22*dAL420*s*L, dAL420*(70*c2+54*s2), 16*dAL420*c*s, 13*dAL420*s*L,
        //        dAL420*(140*s2+156*c2), 22*dAL420*c*L, 16*dAL420*c*s, dAL420*(70*s2+54*c2), -13*dAL420*c*L,
        //        4*dAL420*L2, -13*dAL420*s*L, 13*dAL420*c*L, -3*dAL420*L2,
        //        dAL420*(140*c2+156*s2), -16*dAL420*c*s, 22*dAL420*s*L,
        //        dAL420*(140*s2+156*c2), -22*dAL420*c*L,
        //        4*dAL420*L2 });
        //}

        //[ 140*c^2+156*s^2,         -16*c*s,         -22*s*L,   70*c^2+54*s^2,          16*c*s,          13*s*L]
        //[         -16*c*s, 140*s^2+156*c^2,          22*c*L,          16*c*s,   70*s^2+54*c^2,         -13*c*L]
        //[         -22*s*L,          22*c*L,           4*L^2,         -13*s*L,          13*c*L,          -3*L^2]
        //[   70*c^2+54*s^2,          16*c*s,         -13*s*L, 140*c^2+156*s^2,         -16*c*s,          22*s*L]
        //[          16*c*s,   70*s^2+54*c^2,          13*c*L,         -16*c*s, 140*s^2+156*c^2,         -22*c*L]
        //[          13*s*L,         -13*c*L,          -3*L^2,          22*s*L,         -22*c*L,           4*L^2]
        public IMatrix2D MassMatrix(IElement element)
        {
            double x2 = Math.Pow(element.INodes[1].X - element.INodes[0].X, 2);
            double y2 = Math.Pow(element.INodes[1].Y - element.INodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);
            double L2 = L * L;
            double c = (element.INodes[1].X - element.INodes[0].X) / L;
            double c2 = c * c;
            double s = (element.INodes[1].Y - element.INodes[0].Y) / L;
            double s2 = s * s;
            double dAL420 = Density * SectionArea * L / 420;

            double totalMass = Density * SectionArea * L;
            double totalMassOfDiagonalTerms = 2 * dAL420 * (140 * c2 + 156 * s2) + 2 * dAL420 * (140 * s2 + 156 * c2);
            double scale = totalMass / totalMassOfDiagonalTerms;

            return new SymmetricMatrix2D(new double[] { dAL420*(140*c2+156*s2)*scale, 0, 0, 0, 0, 0,
                dAL420*(140*s2+156*c2)*scale, 0, 0, 0, 0,
                0, 0, 0, 0,
                dAL420*(140*c2+156*s2)*scale, 0, 0,
                dAL420*(140*s2+156*c2)*scale, 0,
                0 });
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            //throw new NotImplementedException();
            var m = MassMatrix(element);
            var lc = m as ILinearlyCombinable;
            lc.LinearCombination(new double[] { RayleighAlpha, RayleighBeta }, new IMatrix2D[] { MassMatrix(element), StiffnessMatrix(element) });
            return m;
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            Vector accelerations = new Vector(6);
            IMatrix2D massMatrix = MassMatrix(element);

            int index = 0;
            foreach (MassAccelerationLoad load in loads)
                foreach (DOFType[] nodalDOFTypes in dofs)
                    foreach (DOFType dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }

            double[] forces = new double[6];
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
            // Method intentionally left empty.
        }

        public void ClearMaterialState()
        {
            // Method intentionally left empty.
        }

        public void ClearMaterialStresses()
        {
            throw new NotImplementedException();
        }
        #endregion
    }
}
