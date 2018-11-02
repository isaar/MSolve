using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Elements
{
    public class ThermalElement3D : IFiniteElement
    {
        private readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.Temperature };
        private readonly DOFType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
        private readonly ThermalMaterial material;

        public ThermalElement3D(IReadOnlyList<Node3D> nodes, IIsoparametricInterpolation3D interpolation,
        IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
        IGaussPointExtrapolation3D gaussPointExtrapolation, ThermalMaterial material)
        {
            this.material = material;
            this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.QuadratureForConsistentMass = quadratureForMass;
            this.QuadratureForStiffness = quadratureForStiffness;

            dofTypes = new DOFType[nodes.Count][];
            for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new DOFType[] { DOFType.Temperature };
        }
        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

        public int ID => throw new NotImplementedException(
            "Element type codes should be in a settings class. Even then it's a bad design choice");

        public IGaussPointExtrapolation3D GaussPointExtrapolation { get; }
        public IIsoparametricInterpolation3D Interpolation { get; }
        public IReadOnlyList<Node3D> Nodes { get; }
        public IQuadrature3D QuadratureForConsistentMass { get; }
        public IQuadrature3D QuadratureForStiffness { get; }

        public bool MaterialModified => throw new NotImplementedException();

        public IElementDOFEnumerator DOFEnumerator { get; set; } = new GenericDOFEnumerator();

        public IMatrix2D MassMatrix(IElement element)
        {
            return BuildCapacityMatrix();
        }

        public Matrix2D BuildCapacityMatrix()
        {
            int numDofs = Nodes.Count;
            var capacity = new Matrix2D(numDofs, numDofs);
            IReadOnlyList<Vector> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
            IReadOnlyList<Matrix2D> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            {
                Matrix2D shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix2D partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
                capacity.AxpyIntoThis(partial, dA);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            capacity.Scale(material.Density * material.SpecialHeatCoeff);
            return capacity;
        }

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            return BuildConductivityMatrix();
        }

        public Matrix2D BuildConductivityMatrix()
        {
            int numDofs = Nodes.Count;
            var conductivity = new Matrix2D(numDofs, numDofs);
            IReadOnlyList<Matrix2D> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                // Calculate the necessary quantities for the integration
                //Matrix2D constitutive = (Matrix2D)(materialsAtGaussPoints[gp].ConstitutiveMatrix); // ugly cast will be removed along with the legacy Matrix classes
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                Matrix2D shapeGradientsCartesian =
                    jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix2D deformation = BuildDeformationMatrix(shapeGradientsCartesian);

                // Contribution of this gauss point to the element stiffness matrix
                Matrix2D partialK = deformation.Transpose() * deformation;
                //Matrix2D partialΚ = deformation.Transpose() * (constitutive * deformation);
                //partialK.Scale(materialsAtGaussPoints[gaussPoint].ThermalConductivity);

                double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight; //TODO: this is used by all methods that integrate. I should cache it.
                conductivity.AxpyIntoThis(partialK, dA * material.ThermalConductivity);
                //conductivity.AxpyIntoThis(partialK, dA * 1);
            }

            conductivity.Scale(1);
            return conductivity;
        }

        private Matrix2D BuildDeformationMatrix(Matrix2D shapeGradientsCartesian)
        {
            //TODO: isn't this just the transpose of [dNi/dxj]?
            var deformation = new Matrix2D(3, Nodes.Count);
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
                deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
                deformation[2, nodeIdx] = shapeGradientsCartesian[nodeIdx, 2];
            }
            return deformation;
        }

        public Matrix2D BuildShapeFunctionMatrix(Vector shapeFunctions) //TODO: reconsider this. As it is, it just returns the shape functions in a Matrix2D
        {
            var array2D = new double[1, shapeFunctions.Length];
            for (int i = 0; i < shapeFunctions.Length; ++i)
            {
                array2D[0, i] = shapeFunctions[i];
            }
            return new Matrix2D(array2D);
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

        public IMatrix2D DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }
    }
}
