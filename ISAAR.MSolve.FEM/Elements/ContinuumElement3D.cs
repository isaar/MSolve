using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials;

namespace ISAAR.MSolve.FEM.Elements
{
    /// <summary>
    /// Represents a continuum finite element for 3D problems. Specific elements (e.g. Hexa8, Hexa20, ...) can be created using
    /// the appropriate <see cref="IIsoparametricInterpolation3D_OLD"/>, <see cref="IQuadrature3D"/> etc. strategies. 
    /// Authors: Dimitris Tsapetis
    /// </summary>
    public class ContinuumElement3D : IStructuralFiniteElement_v2
    {
        private readonly static DOFType[] nodalDOFTypes = new DOFType[] {DOFType.X, DOFType.Y, DOFType.Z};
        private readonly DOFType[][] dofTypes;
        private DynamicMaterial dynamicProperties;
        private readonly IReadOnlyList<ElasticMaterial3D_v2> materialsAtGaussPoints;

        public ContinuumElement3D(IReadOnlyList<Node_v2> nodes, IIsoparametricInterpolation3D interpolation,
            IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
            IGaussPointExtrapolation3D gaussPointExtrapolation,
            IReadOnlyList<ElasticMaterial3D_v2> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
        {
            this.dynamicProperties = dynamicProperties;
            this.materialsAtGaussPoints = materialsAtGaussPoints;
            this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.QuadratureForConsistentMass = quadratureForMass;
            this.QuadratureForStiffness = quadratureForStiffness;

            dofTypes= new DOFType[nodes.Count][];
            for (int i = 0; i < interpolation.NumFunctions; i++)
                dofTypes[i]=new DOFType[]{DOFType.X, DOFType.Y,DOFType.Z};
        }

        public ElementDimensions ElementDimensions => ElementDimensions.ThreeD;

        public int ID=> throw new NotImplementedException(
            "Element type codes should be in a settings class. Even then it's a bad design choice");

        public IGaussPointExtrapolation3D GaussPointExtrapolation { get; }
        public IIsoparametricInterpolation3D Interpolation { get; }
        public IReadOnlyList<Node_v2> Nodes { get; }
        public IQuadrature3D QuadratureForConsistentMass { get; }
        public IQuadrature3D QuadratureForStiffness { get; }


        public Matrix BuildConsistentMassMatrix()
        {
            int numberOfDofs = 3 * Nodes.Count;
            var mass = Matrix.CreateZero(numberOfDofs,numberOfDofs);
            IReadOnlyList<double[]> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            {
                Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix partial = shapeFunctionMatrix.MultiplyRight(shapeFunctionMatrix, true, false);
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
                mass.AxpyIntoThis(partial,dA);
            }
            mass.ScaleIntoThis(dynamicProperties.Density);
            return mass;
        }

        public Matrix BuildLumpedMassMatrix()
        {
            int numberOfDofs = 3 * Nodes.Count;
            var lumpedMass= Matrix.CreateZero(numberOfDofs,numberOfDofs);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            double area = 0;
            for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            {
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                area += jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
            }

            double nodalMass = area * dynamicProperties.Density / Nodes.Count;
            for (int i = 0; i < numberOfDofs; i++) lumpedMass[i, i] = nodalMass;

            return lumpedMass;
        }


        public Matrix BuildStiffnessMatrix()
        {
            int numberOfDofs = 3 * Nodes.Count;
            var stiffness = Matrix.CreateZero(numberOfDofs, numberOfDofs);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                IMatrixView constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                Matrix shapeGradientsCartesian =
                    jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);

                Matrix partial = deformation.ThisTransposeTimesOtherTimesThis(constitutive);
                double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
                stiffness.AxpyIntoThis(partial, dA);
            }

            return stiffness;
        }

        public double[] CalculateAccelerationForces(Element_v2 element, IList<MassAccelerationLoad> loads)
        {
            int numberOfDofs = 3 * Nodes.Count;
            var accelerations = new double[numberOfDofs];
            IMatrix massMatrix = MassMatrix(element);

            foreach (var load in loads)
            {
                int index = 0;
                foreach (var nodalDOFTypes in dofTypes)
                {
                    foreach (var dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }
                }
            }

            return massMatrix.Multiply(accelerations);
        }

        public double[] CalculateForces(Element_v2 element, double[] localTotalDisplacements, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element_v2 element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public Tuple<double[], double[]> CalculateStresses(Element_v2 element, double[] localDisplacements,
            double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialState()
        {
            foreach (var material in materialsAtGaussPoints) material.ClearState();
        }

        public void ClearMaterialStresses()
        {
            foreach (var material in materialsAtGaussPoints) material.ClearStresses();
        }

        public IMatrix DampingMatrix(IElement_v2 element)
        {
            Matrix damping = BuildStiffnessMatrix();
            damping.ScaleIntoThis(dynamicProperties.RayleighCoeffStiffness);
            damping.AxpyIntoThis(MassMatrix(element), dynamicProperties.RayleighCoeffMass);
            return damping;
        }


        public IElementDofEnumerator_v2 DofEnumerator { get; set; } = new GenericDofEnumerator_v2();

        public IList<IList<DOFType>> GetElementDOFTypes(IElement_v2 element) => dofTypes;

        public IMatrix MassMatrix(IElement_v2 element)
        {
            return BuildLumpedMassMatrix();
        }

        public bool MaterialModified
        {
            get
            {
                foreach (ElasticMaterial3D_v2 material in materialsAtGaussPoints)
                    if (material.Modified)
                        return true;
                return false;
            }
        }

        public void ResetMaterialModified()
        {
            foreach (var material in materialsAtGaussPoints) material.ResetModified();
        }

        public void SaveMaterialState()
        {
            foreach (var m in materialsAtGaussPoints) m.SaveState();
        }

        public IMatrix StiffnessMatrix(IElement_v2 element)
        {
            return BuildStiffnessMatrix();
        }


        public (IReadOnlyList<double[]> strains, IReadOnlyList<double[]> stresses) 
            UpdateStrainStressesAtGaussPoints(double[] localDisplacements)
        {
            int numberOfGPs = QuadratureForStiffness.IntegrationPoints.Count;
            var strains = new double[numberOfGPs][];
            var stresses = new double[numberOfGPs][];
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            for (int gp = 0; gp < numberOfGPs; gp++)
            {
                IMatrixView constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                Matrix shapeGrandientsCartesian =
                    jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix deformation = BuildDeformationMatrix(shapeGrandientsCartesian);

                strains[gp] = deformation.Multiply(localDisplacements);
                stresses[gp] = constitutive.Multiply(strains[gp]);
            }

            return (strains, stresses);
        }

        /// <summary>
        /// Assembles the deformation matrix of a solid element.
        /// The calculation are based on <see cref="https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch08.d/AFEM.Ch08.pdf"/>
        /// paragraph 8.4, equation 8.7
        /// </summary>
        /// <param name="shapeGradients"></param>
        /// <returns></returns>
        private Matrix BuildDeformationMatrix(Matrix shapeGradientsCartesian)
        {
            var deformation = Matrix.CreateZero(6, 3 * Nodes.Count);
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; nodeIdx++)
            {
                int col0 = 3 * nodeIdx;
                int col1 = 3 * nodeIdx + 1;
                int col2 = 3 * nodeIdx + 2;

                deformation[0, col0] = shapeGradientsCartesian[nodeIdx, 0];
                deformation[1, col1] = shapeGradientsCartesian[nodeIdx, 1];
                deformation[2, col2] = shapeGradientsCartesian[nodeIdx, 2];

                deformation[3, col0] = shapeGradientsCartesian[nodeIdx, 1];
                deformation[3, col1] = shapeGradientsCartesian[nodeIdx, 0];

                deformation[4, col1] = shapeGradientsCartesian[nodeIdx, 2];
                deformation[4, col2] = shapeGradientsCartesian[nodeIdx, 1];

                deformation[5, col0] = shapeGradientsCartesian[nodeIdx, 2];
                deformation[5, col2] = shapeGradientsCartesian[nodeIdx, 0];
            }

            return deformation;
        }

        /// <summary>
        /// The shape function matrix is 2-by-2n, where n = is the number of shape functions. Row 0 corresponds to dof X, while
        /// row 1 to dof Y, etc.
        /// </summary>
        private Matrix BuildShapeFunctionMatrix(double[] shapeFunctions)
        {
            var shapeFunctionMatrix = Matrix.CreateZero(3, 3 * shapeFunctions.Length);
            for (int i = 0; i < shapeFunctions.Length; i++)
            {
                shapeFunctionMatrix[0, 3 * i] = shapeFunctions[i];
                shapeFunctionMatrix[1, 2 * i + 1] = shapeFunctions[i];
                shapeFunctionMatrix[2, 3 * i + 2] = shapeFunctions[i];
            }
            return shapeFunctionMatrix;
        }

    }
}
