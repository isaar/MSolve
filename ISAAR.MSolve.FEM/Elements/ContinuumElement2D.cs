using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Integration.Points;
using ISAAR.MSolve.FEM.Integration.Quadratures;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Elements
{
    public class ContinuumElement2D: IStructuralFiniteElement
    {
        private readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y };
        private readonly DOFType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
        private readonly Dictionary<GaussPoint2D, ElasticMaterial2D> materialsAtGaussPoints;

        public ContinuumElement2D(IReadOnlyList<Node2D> nodes, IIsoparametricInterpolation2D interpolation,
            IQuadrature2D quadrature, IGaussPointExtrapolation2D gaussPointExtrapolation, 
            Dictionary<GaussPoint2D, ElasticMaterial2D> materialsAtGaussPoints)
        {
            this.materialsAtGaussPoints = materialsAtGaussPoints;
            this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.Quadrature = quadrature;

            dofTypes = new DOFType[nodes.Count][];
            for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new DOFType[] { DOFType.X, DOFType.Y };
        }

        public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

        public int ID => throw new NotImplementedException(
            "Element type codes should be in a settings class. Even then it's a bad design choice");

        public IGaussPointExtrapolation2D GaussPointExtrapolation { get; }
        public IIsoparametricInterpolation2D Interpolation { get; }
        public IReadOnlyList<Node2D> Nodes { get; }
        public IQuadrature2D Quadrature { get; }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        //TODO: this method is probably not necessary at all
        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        //TODO: this method must be changed. It should calculates strains, stresses at GPs or nodes.
        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, 
            double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialState()
        {
            foreach (ElasticMaterial2D m in materialsAtGaussPoints.Values) m.ClearState();
        }

        public void ClearMaterialStresses()
        {
            foreach (ElasticMaterial2D m in materialsAtGaussPoints.Values) m.ClearStresses();
        }

        public IMatrix2D DampingMatrix(Element element)
        {
            throw new NotImplementedException();
        }

        public IFiniteElementDOFEnumerator DOFEnumerator { get; set; } = new GenericDOFEnumerator();

        public IList<IList<DOFType>> GetElementDOFTypes(Element element) => dofTypes;

        public IMatrix2D MassMatrix(Element element)
        {
            throw new NotImplementedException();
        }

        public bool MaterialModified
        {
            get
            {
                foreach (ElasticMaterial2D material in materialsAtGaussPoints.Values)
                    if (material.Modified) return true;
                return false;
            }
        }

        public void ResetMaterialModified()
        {
            foreach (ElasticMaterial2D material in materialsAtGaussPoints.Values) material.ResetModified();
        }

        public void SaveMaterialState()
        {
            foreach (ElasticMaterial2D m in materialsAtGaussPoints.Values) m.SaveState();
        }

        //TODO: why do I need the wrapping element?
        public IMatrix2D StiffnessMatrix(Element element)
        {
            int numDofs = 2 * Nodes.Count;
            var stiffness = new Matrix2D(numDofs, numDofs);
            Dictionary<GaussPoint2D, EvalShapeGradients2D> shapeGradients = 
                Interpolation.EvaluateGradientsAtGaussPoints(Nodes, Quadrature);

            foreach (GaussPoint2D gaussPoint in Quadrature.IntegrationPoints)
            {
                // Calculate the necessary quantities for the integration
                double thickness = materialsAtGaussPoints[gaussPoint].Thickness;
                Matrix2D constitutive = (Matrix2D)(materialsAtGaussPoints[gaussPoint].ConstitutiveMatrix); // ugly cast will be removed along with the retarded legacy Matrix classes
                Matrix2D deformation = CalcDeformationMatrix(shapeGradients[gaussPoint]);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix2D partial = deformation.Transpose() * (constitutive * deformation);
                double dV = thickness * shapeGradients[gaussPoint].Jacobian.Determinant * gaussPoint.Weight;
                stiffness.AxpyIntoThis(partial, dV);
            }
            return stiffness;
        }

        private Matrix2D CalcDeformationMatrix(EvalShapeGradients2D shapeGradients)
        {
            var deformation = new Matrix2D(3, 2 * Nodes.Count);
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                int col1 = 2 * nodeIdx;
                int col2 = 2 * nodeIdx + 1;
                IReadOnlyList<double> dNdX = shapeGradients[nodeIdx];

                deformation[0, col1] = dNdX[0];
                deformation[1, col2] = dNdX[1];
                deformation[2, col1] = dNdX[1];
                deformation[2, col2] = dNdX[0];
            }
            return deformation;
        }
    }
}
