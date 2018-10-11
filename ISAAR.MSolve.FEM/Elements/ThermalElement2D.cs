using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

//TODO: Is there any point in having different material properties per Gauss point?
namespace ISAAR.MSolve.FEM.Elements
{
    public class ThermalElement2D : IFiniteElement
    {
        private readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.Temperature };
        private readonly DOFType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
        private readonly ThermalMaterial material;
        //private readonly Dictionary<GaussPoint2D, ThermalMaterial> materialsAtGaussPoints;


        public ThermalElement2D(double thickness, IReadOnlyList<Node2D> nodes, IIsoparametricInterpolation2D interpolation,
            IQuadrature2D quadratureForStiffness, IQuadrature2D quadratureForConsistentMass,
            IGaussPointExtrapolation2D gaussPointExtrapolation,
            ThermalMaterial material)
        {
            this.material = material;
            this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.QuadratureForConsistentMass = quadratureForConsistentMass;
            this.QuadratureForStiffness = quadratureForStiffness;
            this.Thickness = thickness;

            dofTypes = new DOFType[nodes.Count][];
            for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new DOFType[] { DOFType.Temperature };
        }

        public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

        public int ID => throw new NotImplementedException(
            "Element type codes should be in a settings class. Even then it's a bad design choice");

        public IGaussPointExtrapolation2D GaussPointExtrapolation { get; }
        public IIsoparametricInterpolation2D Interpolation { get; }
        public IReadOnlyList<Node2D> Nodes { get; }
        public IQuadrature2D QuadratureForConsistentMass { get; }
        public IQuadrature2D QuadratureForStiffness { get; }
        public double Thickness { get; }

        public bool MaterialModified => throw new NotImplementedException();

        public IElementDOFEnumerator DOFEnumerator { get; set; } = new GenericDOFEnumerator();

        public IMatrix2D MassMatrix(IElement element)
        {
            return BuildCapacityMatrix();
        }

        public Matrix2D BuildCapacityMatrix()
        {
            int numDofs = Nodes.Count;
            var mass = new Matrix2D(numDofs, numDofs);
            Dictionary<GaussPoint2D, EvalInterpolation2D> evalInterpolations =
                Interpolation.EvaluateAllAtGaussPoints(Nodes, QuadratureForConsistentMass);

            foreach (GaussPoint2D gaussPoint in QuadratureForConsistentMass.IntegrationPoints)
            {
                Matrix2D shapeFunctionMatrix = evalInterpolations[gaussPoint].BuildScalarShapeFunctionMatrix();
                Matrix2D partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
                double dA = evalInterpolations[gaussPoint].Jacobian.Determinant * gaussPoint.Weight;
                mass.AxpyIntoThis(partial, dA);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            mass.Scale(Thickness * material.Density * material.SpecialHeatCoeff);
            return mass;
        }

        public Matrix2D BuildConductivityMatrix()
        {
            int numDofs = Nodes.Count;
            var conductivity = new Matrix2D(numDofs, numDofs);
            Dictionary<GaussPoint2D, EvalShapeGradients2D> shapeGradients =
                Interpolation.EvaluateGradientsAtGaussPoints(Nodes, QuadratureForStiffness);

            foreach (GaussPoint2D gaussPoint in QuadratureForStiffness.IntegrationPoints)
            {
                // Calculate the necessary quantities for the integration
                //Matrix2D constitutive = (Matrix2D)(materialsAtGaussPoints[gaussPoint].ConstitutiveMatrix); // ugly cast will be removed along with the retarded legacy Matrix classes
                Matrix2D deformation = BuildDeformationMatrix(shapeGradients[gaussPoint]);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix2D partialK = deformation.Transpose() * deformation;
                //Matrix2D partialΚ = deformation.Transpose() * (constitutive * deformation);
                //partialK.Scale(materialsAtGaussPoints[gaussPoint].ThermalConductivity);

                double dA = shapeGradients[gaussPoint].Jacobian.Determinant * gaussPoint.Weight;
                conductivity.AxpyIntoThis(partialK, dA * material.ThermalConductivity);
            }
            conductivity.Scale(Thickness);
            return conductivity;
        }

        // Provatidis uses two distinct vectors K = N,x^T * k * N,x + N,y^T * k * N,y
        //private (Matrix2D dNdX, Matrix2D dNdY) CalcdNdx(EvalShapeGradients2D shapeGrad)
        //{
        //    int n = Nodes.Count;
        //    var dNdX = new double[n, 1];
        //    var dNdY = new double[n, 1];
        //    for (int i = 0; i < n; ++i)
        //    {
        //        dNdX[i, 0] = shapeGrad[i][0];
        //        dNdY[i, 0] = shapeGrad[i][1];
        //    }
        //    return (new Matrix2D(dNdX), new Matrix2D(dNdY));
        //}

        private Matrix2D BuildDeformationMatrix(EvalShapeGradients2D shapeGradients)
        {
            var deformation = new Matrix2D(2, Nodes.Count);
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                IReadOnlyList<double> dNdX = shapeGradients[nodeIdx];
                deformation[0, nodeIdx] = dNdX[0];
                deformation[1, nodeIdx] = dNdX[1];
            }
            return deformation;
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element) => dofTypes;

        public void ResetMaterialModified()
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            throw new NotImplementedException();
        }

        public void SaveMaterialState()
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialState()
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialStresses()
        {
            throw new NotImplementedException();
        }

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            return BuildConductivityMatrix();
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }
    }
}
