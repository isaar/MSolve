using System;
using System.Collections.Generic;
using System.Net.Http.Headers;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Elements
{
	/// <summary>
	/// Represents a continuum finite element for 3D problems. Specific elements (e.g. Hexa8, Hexa20, ...) can be created using
	/// the appropriate <see cref="IIsoparametricInterpolation3D"/>, <see cref="IQuadrature2D"/> etc. strategies. 
	/// Authors: Dimitris Tsapetis
	/// </summary>
	public class ContinuumElement3D:IStructuralFiniteElement
	{
		private readonly static DOFType[] nodalDOFTypes = new DOFType[] {DOFType.X, DOFType.Y, DOFType.Z};
		private readonly DOFType[][] dofTypes;
		private DynamicMaterial dynamicProperties;
		private readonly Dictionary<GaussPoint3D, ElasticMaterial3D> materialsAtGaussPoints;

		public ContinuumElement3D(IReadOnlyList<Node3D> nodes, IIsoparametricInterpolation3D interpolation,
			IQuadrature3D quadratureForStiffness, IQuadrature3D quadratureForMass,
			IGaussPointExtrapolation3D gaussPointExtrapolation,
			Dictionary<GaussPoint3D, ElasticMaterial3D> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
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
			Dictionary<GaussPoint3D, EvalInterpolation3D> evalInterpolations =
				Interpolation.EvaluateAllAtGaussPoints(Nodes, QuadratureForConsistentMass);

			foreach (var gaussPoint in QuadratureForConsistentMass.IntegrationPoints)
			{
				Matrix2D shapeFunctionMatrix = evalInterpolations[gaussPoint].BuildShapeFunctionMatrix();
				Matrix2D partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
				double dA = evalInterpolations[gaussPoint].Jacobian.Determinant * gaussPoint.Weight;
				mass.AxpyIntoThis(partial,dA);
			}
			mass.Scale(dynamicProperties.Density);
			return mass;
		}

		public Matrix2D BuildLumpedMassMatrix()
		{
			int numberOfDofs = 3 * Nodes.Count;
			var lumpedMass= new Matrix2D(numberOfDofs,numberOfDofs);
			Dictionary<GaussPoint3D, EvalShapeGradients3D> shapeGradients =
				Interpolation.EvaluateGradientsAtGaussPoints(Nodes, QuadratureForConsistentMass);
			double area = 0;
			foreach (var gaussPoint in QuadratureForConsistentMass.IntegrationPoints)
			{
				area += shapeGradients[gaussPoint].Jacobian.Determinant * gaussPoint.Weight;
			}

			double nodalMass = area * dynamicProperties.Density / Nodes.Count;
			for (int i = 0; i < numberOfDofs; i++) lumpedMass[i, i] = nodalMass;

			return lumpedMass;
		}


		public Matrix2D BuildStiffnessMatrix()
		{
			int numberOfDofs = 3 * Nodes.Count;
			var stiffness=new Matrix2D(numberOfDofs, numberOfDofs);
			Dictionary<GaussPoint3D, EvalShapeGradients3D> shapeGradients =
				Interpolation.EvaluateGradientsAtGaussPoints(Nodes, QuadratureForStiffness);

			foreach (var gaussPoint in QuadratureForStiffness.IntegrationPoints)
			{
				Matrix2D constitutive = (Matrix2D) materialsAtGaussPoints[gaussPoint].ConstitutiveMatrix;
				Matrix2D deformation = BuildDeformationMatrix(shapeGradients[gaussPoint]);

				Matrix2D partial = deformation.Transpose() * (constitutive * deformation);
				double dA = shapeGradients[gaussPoint].Jacobian.Determinant * gaussPoint.Weight;
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
			foreach (var material in materialsAtGaussPoints.Values) material.ClearState();
		}

		public void ClearMaterialStresses()
		{
			foreach (var material in materialsAtGaussPoints.Values) material.ClearStresses();
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
				foreach (ElasticMaterial3D material in materialsAtGaussPoints.Values)
					if (material.Modified)
						return true;
				return false;
			}
		}

		public void ResetMaterialModified()
		{
			foreach (var material in materialsAtGaussPoints.Values) material.ResetModified();
		}

		public void SaveMaterialState()
		{
			foreach (var m in materialsAtGaussPoints.Values) m.SaveState();
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
			Dictionary<GaussPoint3D, EvalShapeGradients3D> shapeGradients =
				Interpolation.EvaluateGradientsAtGaussPoints(Nodes, QuadratureForStiffness);
			for (int i = 0; i < numberOfGPs; i++)
			{
				GaussPoint3D gaussPoint = QuadratureForStiffness.IntegrationPoints[i];
				IMatrix2D constitutive = materialsAtGaussPoints[gaussPoint].ConstitutiveMatrix;
				Matrix2D deformation = BuildDeformationMatrix(shapeGradients[gaussPoint]);
				strains[i]= new double[6];
				deformation.Multiply(localDisplacementVector,strains[i]);
				stresses[i]=new double[6];
				constitutive.Multiply(new Vector(strains[i]),stresses[i]);
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
		private Matrix2D BuildDeformationMatrix(EvalShapeGradients3D shapeGradients)
		{
			var deformation= new Matrix2D(6,3*Nodes.Count);
			for (int nodeIdx = 0; nodeIdx < Nodes.Count; nodeIdx++)
			{
				int col1 = 3 * nodeIdx;
				int col2 = 3 * nodeIdx + 1;
				int col3 = 3 * nodeIdx + 2;
				IReadOnlyList<double> dNdx = shapeGradients[nodeIdx];
				deformation[0, col1] = dNdx[0];
				deformation[1, col2] = dNdx[1];
				deformation[2, col3] = dNdx[2];

				deformation[3, col1] = dNdx[1];
				deformation[3, col2] = dNdx[0];

				deformation[4, col2] = dNdx[2];
				deformation[4, col3] = dNdx[1];

				deformation[5, col1] = dNdx[2];
				deformation[5, col3] = dNdx[0];
			}

			return deformation;
		}

	}
}
