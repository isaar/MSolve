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
    /// <summary>
    /// Represents a continuum finite element for 3D problems. Specific elements (e.g. Hexa8, Hexa20, ...) can be created using
    /// the appropriate <see cref="IIsoparametricInterpolation3D"/>, <see cref="IQuadrature3D"/> etc. strategies. 
    /// Authors: Dimitris Tsapetis
    /// </summary>
    public class ContinuumElement3D:IStructuralFiniteElement
    {
        private readonly static DOFType[] nodalDOFTypes = new DOFType[] {DOFType.X, DOFType.Y, DOFType.Z};
        private readonly DOFType[][] dofTypes;
        private DynamicMaterial dynamicProperties;
        private readonly IReadOnlyList<ElasticMaterial3D> materialsAtGaussPoints;

        public ContinuumElement3D(IReadOnlyList<Node3D> nodes, IIsoparametricInterpolation3D interpolation,
            IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
            IGaussPointExtrapolation3D gaussPointExtrapolation,
            IReadOnlyList<ElasticMaterial3D> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
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
        public IReadOnlyList<Node3D> Nodes { get; }
        public IQuadrature3D QuadratureForConsistentMass { get; }
        public IQuadrature3D QuadratureForStiffness { get; }


        public Matrix2D BuildConsistentMassMatrix()
        {
            int numberOfDofs = 3 * Nodes.Count;
            var mass= new Matrix2D(numberOfDofs,numberOfDofs);
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
                mass.AxpyIntoThis(partial,dA);
            }
            mass.Scale(dynamicProperties.Density);
            return mass;
        }

        public Matrix2D BuildLumpedMassMatrix()
        {
            int numberOfDofs = 3 * Nodes.Count;
            var lumpedMass= new Matrix2D(numberOfDofs,numberOfDofs);
            IReadOnlyList<Matrix2D> shapeGradientsNatural =
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


        public Matrix2D BuildStiffnessMatrix()
        {
            int numberOfDofs = 3 * Nodes.Count;
            var stiffness=new Matrix2D(numberOfDofs, numberOfDofs);
            IReadOnlyList<Matrix2D> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
            {
                Matrix2D constitutive = (Matrix2D) materialsAtGaussPoints[gp].ConstitutiveMatrix;
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                Matrix2D shapeGradientsCartesian =
                    jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix2D deformation = BuildDeformationMatrix(shapeGradientsCartesian);

                Matrix2D partial = deformation.Transpose() * (constitutive * deformation);
                double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight;
                stiffness.AxpyIntoThis(partial,dA);
            }

            return stiffness;
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            int numberOfDofs = 3 * Nodes.Count;
            Vector accelerations= new Vector(numberOfDofs);
            IMatrix2D massMatrix =MassMatrix(element);

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

            double[] forces = new double[numberOfDofs];
            massMatrix.Multiply(accelerations,forces);
            return forces;
        }

        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements,
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

        public IMatrix2D DampingMatrix(IElement element)
        {
            Matrix2D damping = BuildStiffnessMatrix();
            damping.Scale(dynamicProperties.RayleighCoeffStiffness);
            damping.AxpyIntoThis(MassMatrix(element), dynamicProperties.RayleighCoeffMass);
            return damping;
        }


        public IElementDOFEnumerator DOFEnumerator { get; set; }=new GenericDOFEnumerator();

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element) => dofTypes;

        public IMatrix2D MassMatrix(IElement element)
        {
            return BuildLumpedMassMatrix();
        }

        public bool MaterialModified
        {
            get
            {
                foreach (ElasticMaterial3D material in materialsAtGaussPoints)
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

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            return BuildStiffnessMatrix();
        }


        public (IReadOnlyList<double[]> strains, IReadOnlyList<double[]> stresses) UpdateStrainStressesAtGaussPoints(double[] localDisplacements)
        {
            var localDisplacementVector= new Vector(localDisplacements);
            int numberOfGPs = QuadratureForStiffness.IntegrationPoints.Count;
            var strains = new double[numberOfGPs][];
            var stresses = new double[numberOfGPs][];
            IReadOnlyList<Matrix2D> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            for (int gp = 0; gp < numberOfGPs; gp++)
            {
                IMatrix2D constitutive = materialsAtGaussPoints[gp].ConstitutiveMatrix;
                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                Matrix2D shapeGrandientsCartesian =
                    jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
                Matrix2D deformation = BuildDeformationMatrix(shapeGrandientsCartesian);

                strains[gp]= new double[6];
                deformation.Multiply(localDisplacementVector,strains[gp]);
                stresses[gp]=new double[6];
                constitutive.Multiply(new Vector(strains[gp]),stresses[gp]);
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
        private Matrix2D BuildDeformationMatrix(Matrix2D shapeGradientsCartesian)
        {
            var deformation = new Matrix2D(6, 3 * Nodes.Count);
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
        private Matrix2D BuildShapeFunctionMatrix(Vector shapeFunctions)
        {
            var array2D = new double[3, 3 * shapeFunctions.Length];
            for (int i = 0; i < shapeFunctions.Length; i++)
            {
                array2D[0, 3 * i] = shapeFunctions[i];
                array2D[1, 2 * i + 1] = shapeFunctions[i];
                array2D[2, 3 * i + 2] = shapeFunctions[i];
            }
            return new Matrix2D(array2D);
        }

    }
}
